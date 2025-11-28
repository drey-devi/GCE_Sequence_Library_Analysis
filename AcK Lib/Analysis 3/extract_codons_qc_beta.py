#!/usr/bin/env python3
"""
Script: extract_codons_qc.py

A tool with two subcommands:

1. extract – Full extraction mode:
   • Optionally orient reads to the reference (auto-detect forward/reverse)
   • Align (or load) MSA, extract codons → detailed CSV
   • Choose aligner: MUSCLE or MAFFT
   • For either aligner, run full MSA or reference-only mode
   • Optional summary CSV & plots
   • Supports multiple FASTQ inputs: outputs are zipped per FASTQ file as <basename>_analysed.zip

2. post    – Post-process-only mode:
   Take an existing detailed CSV → summary CSV + plots + protein-combo graphs

Run `extract_codons_qc.py <command> --help` for subcommand-specific options.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
import shutil
import zipfile
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(prog='extract_codons_qc.py')
    subparsers = parser.add_subparsers(dest='command', required=True)

    ex = subparsers.add_parser('extract', help='Extract codons from FASTQ(s)')
    ex.add_argument('-r', '--ref-fasta',      required=True, help='Reference FASTA (single record)')
    ex.add_argument('-i', '--input-fastq',    nargs='+', required=True, help='One or more FASTQ files')
    ex.add_argument('-s', '--sites',          required=True, help='Comma-separated codon indices (1-based)')
    ex.add_argument('--ref-start',            type=int, default=1, help='1-based codon index where region starts')
    ex.add_argument('--ref-end',              type=int, help='1-based codon index where region ends')
    ex.add_argument('--aligner',              choices=['muscle','mafft'], default='muscle', help='Alignment tool')
    ex.add_argument('-a', '--algorithm',      choices=['align','super5'], default='super5', help='MUSCLE v5 algorithm')
    ex.add_argument('--muscle-bin',           default='muscle', help='Path to MUSCLE v5')
    ex.add_argument('--profile',              action='store_true', help='Reference-only alignment')
    ex.add_argument('--aligned-fasta',        help='Existing MSA FASTA to load')
    ex.add_argument('--aligned-fasta-out',    help='Write generated MSA FASTA to this path')
    ex.add_argument('--summary-csv',          action='store_true', help='Generate summary CSV per input')
    ex.add_argument('--plots',                action='store_true', help='Generate plots per input')

    ps = subparsers.add_parser('post', help='Post-process existing detailed CSV')
    ps.add_argument('--detailed-csv', required=True, help='Detailed CSV with codon data')
    ps.add_argument('--summary-csv',  help='Path for summary CSV')
    ps.add_argument('--plots',        action='store_true', help='Generate plots')
    ps.add_argument('--plot-dir',     help='Directory for plots')

    return parser.parse_args()


def orient_read(read_seq: str, ref_seq: Seq) -> str:
    fwd = pairwise2.align.globalxx(ref_seq, read_seq, score_only=True)
    rev = pairwise2.align.globalxx(ref_seq, str(Seq(read_seq).reverse_complement()), score_only=True)
    return read_seq if fwd >= rev else str(Seq(read_seq).reverse_complement())


def run_muscle(in_fa, out_fa, args):
    if args.profile:
        cmd = [args.muscle_bin, '-profile', '-in1', args._ref, '-in2', args._reads, '-out', out_fa]
    else:
        cmd = [args.muscle_bin, f'-{args.algorithm}', in_fa, '-output', out_fa]
    subprocess.run(cmd, check=True)


def run_mafft(in_fa, out_fa, args):
    if args.profile:
        cmd = ['mafft', '--addfragments', args._reads, '--keeplength', '--anysymbol', args._ref]
        with open(out_fa, 'w') as wf:
            subprocess.run(cmd, check=True, stdout=wf)
    else:
        cmd = ['mafft', '--auto', in_fa]
        with open(out_fa, 'w') as wf:
            subprocess.run(cmd, check=True, stdout=wf)


def map_refpos_to_alncols(aln_seq):
    mapping, idx = {}, 0
    for col_idx, nt in enumerate(str(aln_seq.seq)):
        if nt != '-': mapping[idx] = col_idx; idx += 1
    return mapping


def get_original_index(gapped_seq, col_idx):
    if gapped_seq[col_idx] == '-': return None
    return sum(1 for c in gapped_seq[:col_idx] if c != '-')


def translate_codon(codon):
    if pd.isna(codon) or '-' in codon or len(codon) != 3: return None
    return str(Seq(codon).translate())


def post_process(detailed_csv, summary_csv=None, generate_plots=False, plot_dir=None):
    df = pd.read_csv(detailed_csv, sep=',', engine='python', on_bad_lines='warn')
    sites = sorted(int(c.split('_')[1]) for c in df.columns if c.endswith('_codon'))
    df['combo'] = df.apply(lambda r: ''.join(r[f'site_{x}_aa'] for x in sites if pd.notna(r[f'site_{x}_aa'])), axis=1)

    if summary_csv:
        rows=[]
        for s in sites:
            cnt=Counter(df[f'site_{s}_aa'].dropna()); tot=sum(cnt.values())
            for rnk,(aa,ct) in enumerate(cnt.most_common(),1): rows.append({'site':s,'aa':aa,'count':ct,'rank':rnk,'freq':ct/tot})
        pd.DataFrame(rows).to_csv(summary_csv,index=False)

    if generate_plots:
        outd = plot_dir or os.getcwd(); os.makedirs(outd, exist_ok=True)
        # Combined AA top5 per site
        combined = {}
        for s in sites:
            combined[s] = Counter(df[f'site_{s}_aa'].dropna()).most_common(5)
        plt.figure()
        for s, vals in combined.items():
            aas, cts = zip(*vals) if vals else ([], [])
            plt.bar([f"{s}:{aa}" for aa in aas], cts)
        plt.xticks(rotation=90)
        plt.title('Top5 AAs per site')
        plt.tight_layout(); plt.savefig(os.path.join(outd, 'aa_top5_per_site.png')); plt.close()

        # Combined quality distribution per site
        qdata=[]
        for s in sites:
            arr = df[[f'site_{s}_q{j}' for j in (1,2,3)]].apply(pd.to_numeric, errors='coerce').values.flatten()
            qdata.append([x for x in arr if not pd.isna(x)])
        plt.figure(); plt.boxplot(qdata, tick_labels=sites)
        plt.xticks(rotation=90)
        plt.title('Quality distribution per site')
        plt.ylabel('PHRED')
        plt.tight_layout(); plt.savefig(os.path.join(outd, 'quality_per_site.png')); plt.close()

    print(f'Post-processing done for {detailed_csv}')


def extract_mode(args):
    refs = list(SeqIO.parse(args.ref_fasta, 'fasta'))
    if len(refs) != 1: sys.exit('Reference FASTA must have exactly one record')
    ref = refs[0]; total = len(ref.seq)//3
    start,end = args.ref_start, args.ref_end or total
    if start<1 or end>total or start> end: sys.exit('Invalid ref-start/end')
    region = ref.seq[(start-1)*3:end*3]

    for fq in args.input_fastq:
        base = os.path.splitext(os.path.basename(fq))[0]
        det_csv = f'{base}_detailed.csv'
        sum_flag = args.summary_csv
        plot_flag = args.plots
        plot_dir = f'{base}_plots' if plot_flag else None
        aln_out = args.aligned_fasta_out or f'{base}_aligned.fa'

        reads = list(SeqIO.parse(fq,'fastq'))
        qual_dict = {r.id:r.letter_annotations.get('phred_quality',[]) for r in reads}
        oriented = [SeqRecord(Seq(orient_read(str(r.seq),region)),id=r.id,description='') for r in reads]

        with tempfile.TemporaryDirectory() as td:
            aln_file = os.path.join(td,'aligned.fa')
            if args.profile:
                ref_fa = os.path.join(td,'ref.fa'); SeqIO.write(SeqRecord(region,id=ref.id,description=''),ref_fa,'fasta')
                reads_fa = os.path.join(td,'reads.fa'); SeqIO.write(oriented,reads_fa,'fasta')
                args._ref,args._reads = ref_fa,reads_fa
                if args.aligner=='muscle': run_muscle(None,aln_file,args)
                else: run_mafft(None,aln_file,args)
            else:
                in_fa = os.path.join(td,'in.fa')
                with open(in_fa,'w') as wf:
                    SeqIO.write(SeqRecord(region,id=ref.id,description=''),wf,'fasta')
                    SeqIO.write(oriented,wf,'fasta')
                if args.aligner=='muscle': run_muscle(in_fa,aln_file,args)
                else: run_mafft(in_fa,aln_file,args)
            shutil.copy(aln_file, aln_out)
            aln = AlignIO.read(aln_file,'fasta')

        ref_aln = next(a for a in aln if a.id==ref.id)
        pos2col = map_refpos_to_alncols(ref_aln)
        sites = [int(x) for x in args.sites.split(',')]
        cols = [pos2col[(s-1)*3+off] for s in sites for off in (0,1,2)]

        with open(det_csv,'w',newline='') as cf:
            w = csv.writer(cf)
            hdr = ['seq_id']
            for s in sites:
                hdr += [f'site_{s}_codon', f'site_{s}_aa'] + [f'site_{s}_q{j}' for j in (1,2,3)]
            w.writerow(hdr)
            for rec in aln:
                seqstr = str(rec.seq); ql = qual_dict.get(rec.id, [])
                row=[rec.id]
                for i in range(len(sites)):
                    sub = cols[3*i:3*i+3]
                    cod,qs = [],[]
                    for cidx in sub:
                        nt=seqstr[cidx]
                        if nt=='-': cod.append('-'); qs.append('NaN')
                        else: cod.append(nt); oi=get_original_index(seqstr,cidx); qs.append(str(ql[oi]) if oi is not None and oi<len(ql) else 'NaN')
                    cd=''.join(cod); aa='NaN'
                    if '-' not in cd: aa=str(Seq(cd).translate())
                    row += [cd,aa] + qs
                w.writerow(row)

        # post-process and generate zip
        post_process(det_csv, summary_csv=(f'{base}_summary.csv' if sum_flag else None), generate_plots=plot_flag, plot_dir=plot_dir)
        zipn = f'{base}_analysed.zip'
        with zipfile.ZipFile(zipn, 'w') as zf:
            zf.write(det_csv)
            if sum_flag:
                sc = f'{base}_summary.csv'; zf.write(sc)
            zf.write(aln_out)
            if plot_flag:
                for root, _, files in os.walk(plot_dir):
                    for fn in files:
                        path=os.path.join(root,fn)
                        arc=os.path.join(os.path.basename(plot_dir),fn)
                        zf.write(path, arcname=arc)
        # cleanup files
        os.remove(det_csv)
        if sum_flag: os.remove(sc)
        os.remove(aln_out)
        if plot_flag: shutil.rmtree(plot_dir)
        print(f'Created {zipn}')


def main():
    args = parse_args()
    if args.command == 'extract':
        extract_mode(args)
    else:
        post_process(args.detailed_csv, summary_csv=args.summary_csv, generate_plots=args.plots, plot_dir=args.plot_dir)

if __name__ == '__main__':
    main()
