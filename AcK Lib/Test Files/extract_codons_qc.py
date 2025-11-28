#!/usr/bin/env python3
"""
Script: extract_codons_qc.py

A tool with two subcommands:

1. extract – Full extraction mode:
   • Optionally orient reads to the reference (auto-detect forward/reverse)
   • Align (or load) MSA, extract codons → detailed CSV
   • Choose aligner: MUSCLE or MAFFT
   • For either aligner, run full MSA or reference-only (profile/addfragments) mode
   • Optional: save MSA FASTA (--aligned-fasta-out), summary CSV, and plots
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

    # extract subcommand
    ex = subparsers.add_parser('extract', help='Extract codons from FASTQ(s)')
    ex.add_argument('-r', '--ref-fasta',      required=True, help='Reference FASTA (single record)')
    ex.add_argument('-i', '--input-fastq',    nargs='+', required=True, help='One or more FASTQ files')
    ex.add_argument('-s', '--sites',          required=True, help='Comma-separated codon indices (1-based)')
    ex.add_argument('--ref-start',            type=int, default=1, help='1-based codon index where region starts')
    ex.add_argument('--ref-end',              type=int, help='1-based codon index where region ends')
    ex.add_argument('--aligner',              choices=['muscle','mafft'], default='muscle', help='Alignment tool')
    ex.add_argument('-a', '--algorithm',      choices=['align','super5'], default='super5', help='MUSCLE v5 algorithm')
    ex.add_argument('--muscle-bin',           default='muscle', help='Path to MUSCLE v5')
    ex.add_argument('--profile',              action='store_true', help='Reference-only alignment (profile/addfragments)')
    ex.add_argument('--min-quality',          type=int, default=0, help='PHRED threshold for low-quality')
    ex.add_argument('--aligned-fasta',        help='Existing MSA FASTA to load')
    ex.add_argument('--aligned-fasta-out',    help='Write generated MSA FASTA to this path')
    ex.add_argument('--summary-csv',          action='store_true', help='Generate summary CSV per input')
    ex.add_argument('--plots',                action='store_true', help='Generate plots per input')

    # post subcommand
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
        cmd = [args.muscle_bin, '-profile', '-in1', args._ref_only_fa, '-in2', args._reads_fa, '-out', out_fa]
    else:
        cmd = [args.muscle_bin, f'-{args.algorithm}', in_fa, '-output', out_fa]
    subprocess.run(cmd, check=True)


def run_mafft(in_fa, out_fa, args):
    if args.profile:
        cmd = ['mafft', '--addfragments', args._reads_fa, '--keeplength', '--anysymbol', args._ref_only_fa]
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
            aa_list = list(df[f'site_{s}_aa'].dropna())
            cnt = Counter(aa_list); tot = sum(cnt.values())
            for rnk,(aa_,ct) in enumerate(cnt.most_common(),1): rows.append({'site':s,'aa':aa_,'count':ct,'rank':rnk,'freq':ct/tot})
        pd.DataFrame(rows).to_csv(summary_csv,index=False)
        combos_file = os.path.splitext(summary_csv)[0] + '_combos.csv'
        combos = df[['seq_id','combo']].copy(); combos['count'] = combos['combo'].map(Counter(combos['combo']))
        combos.to_csv(combos_file,index=False)

    if generate_plots:
        outd = plot_dir or os.getcwd(); os.makedirs(outd, exist_ok=True)
        # per-site plots
        for s in sites:
            lst=list(df[f'site_{s}_aa'].dropna()); top5=Counter(lst).most_common(5)
            if top5:
                aas,cts=zip(*top5); plt.figure(); plt.bar(aas,cts)
                plt.title(f'Top5 AAs at site {s}'); plt.tight_layout(); plt.savefig(os.path.join(outd,f'top5_aa_{s}.png')); plt.close()
            stc=df.get(f'site_{s}_status')
            if stc is not None:
                plt.figure(); stc.value_counts().plot(kind='bar')
                plt.title(f'Status dist site {s}'); plt.tight_layout(); plt.savefig(os.path.join(outd,f'status_{s}.png')); plt.close()
            qcs=[f'site_{s}_q{j}' for j in (1,2,3)]; qdf=df[qcs].apply(pd.to_numeric,errors='coerce').dropna()
            plt.figure(); plt.boxplot([qdf[c] for c in qcs], tick_labels=qcs)
            plt.title(f'Quality scores at site {s}'); plt.ylabel('PHRED'); plt.tight_layout()
            plt.savefig(os.path.join(outd,f'quality_site_{s}.png')); plt.close()
        # protein-combo plots
        combos = df['combo'].dropna()
        top5_combos = Counter(combos).most_common(5)
        if top5_combos:
            labs,cts = zip(*top5_combos)
            plt.figure(); plt.bar(labs,cts)
            plt.title('Top5 protein combos'); plt.xticks(rotation=45,ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outd,'top5_protein_combos.png')); plt.close()
        # avg quality per combo
        df['avg_quality'] = df[[f"site_{s}_q{j}" for s in sites for j in (1,2,3)]].apply(pd.to_numeric, errors='coerce').mean(axis=1)
        if top5_combos:
            sel = [c for c,_ in top5_combos]
            cq = df[df['combo'].isin(sel)].groupby('combo')['avg_quality'].mean()
            plt.figure(); plt.bar(cq.index, cq.values)
            plt.title('Avg quality by top5 combos'); plt.xticks(rotation=45,ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outd,'combo_vs_quality_top5.png')); plt.close()

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
        sum_csv = f'{base}_summary.csv' if args.summary_csv else None
        plot_dir = f'{base}_plots' if args.plots else None
        aln_out  = args.aligned_fasta_out or f'{base}_aligned.fa'

        reads = list(SeqIO.parse(fq,'fastq'))
        if not reads: print(f'No reads in {fq}')
        qual_dict = {r.id:r.letter_annotations.get('phred_quality',[]) for r in reads}
        oriented = [SeqRecord(Seq(orient_read(str(r.seq),region)),id=r.id,description='') for r in reads]

        with tempfile.TemporaryDirectory() as td:
            aln_file = os.path.join(td,'aligned.fa')
            if args.profile:
                ref_fa = os.path.join(td,'ref.fa'); SeqIO.write(SeqRecord(region,id=ref.id,description=''),ref_fa,'fasta')
                reads_fa = os.path.join(td,'reads.fa'); SeqIO.write(oriented,reads_fa,'fasta')
                args._ref_only_fa,args._reads_fa = ref_fa,reads_fa
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
        base_idxs = []
        for s in sites: b0=(s-1)*3; base_idxs += [b0,b0+1,b0+2]
        aln_cols = [pos2col[i] for i in base_idxs]

        with open(det_csv,'w',newline='') as cf:
            w=csv.writer(cf); hdr=['seq_id']
            for s in sites: hdr += [f'site_{s}_codon',f'site_{s}_aa']+[f'site_{s}_q{j}' for j in (1,2,3)]
            w.writerow(hdr)
            for rec in aln:
                seqstr=str(rec.seq); ql=qual_dict.get(rec.id,[])
                row=[rec.id]
                for i,s in enumerate(sites):
                    sub = aln_cols[3*i:3*i+3]
                    cod,qs = [],[]
                    for cidx in sub:
                        nt=seqstr[cidx]
                        if nt=='-': cod.append('-'); qs.append('NaN')
                        else: cod.append(nt); oi=get_original_index(seqstr,cidx); qs.append(str(ql[oi]) if oi is not None and oi<len(ql) else 'NaN')
                    cd=''.join(cod); aa='NaN'
                    if '-' not in cd: aa=str(Seq(cd).translate())
                    row += [cd,aa] + qs
                w.writerow(row)

        post_process(det_csv, summary_csv=sum_csv, generate_plots=args.plots, plot_dir=plot_dir)

        zipname = f'{base}_analysed'
        with zipfile.ZipFile(f'{zipname}.zip','w') as zf:
            zf.write(det_csv)
            if sum_csv:
                zf.write(sum_csv)
                combos_file = os.path.splitext(sum_csv)[0] + '_combos.csv'
                if os.path.exists(combos_file): zf.write(combos_file)
            zf.write(aln_out)
            if plot_dir:
                for root,_,files in os.walk(plot_dir):
                    for fname in files:
                        fpath=os.path.join(root,fname)
                        arc=os.path.join(os.path.basename(plot_dir),fname)
                        zf.write(fpath,arcname=arc)
        print(f'Zipped results to {zipname}.zip')


def main():
    args = parse_args()
    if args.command=='extract': extract_mode(args)
    else: post_process(args.detailed_csv, summary_csv=args.summary_csv, generate_plots=args.plots, plot_dir=args.plot_dir)

if __name__=='__main__': main()
