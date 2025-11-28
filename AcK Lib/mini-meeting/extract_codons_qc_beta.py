#!/usr/bin/env python3
"""
Script: extract_codons_qc.py

Two subcommands:

1. extract – Full extraction:
   • Orient reads to reference (auto–detect strand)
   • Align (MUSCLE or MAFFT), extract codons → detailed CSV
   • Optional: summary CSV, plots (site & combo)
   • Compute extra AA mutations outside specified codon sites → extra_mutations CSV & plot
   • Bundle everything into <basename>_analysed.zip and clean up

2. post – Post-process an existing detailed CSV:
   • Summary CSV, plots, extra-AA-mutations report

Usage:
  extract_codons_qc.py --package
  extract_codons_qc.py extract -r REF.fa -i reads.fastq -s 5,10,15 -o out.csv --summary-csv --plots --plot-dir plots
  extract_codons_qc.py post --detailed-csv out.csv --summary-csv sum.csv --plots --plot-dir plots
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
import shutil
import zipfile
import platform
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(prog='extract_codons_qc.py')
    parser.add_argument(
        '--package', action='store_true',
        help='Install Python libs + MAFFT/MUSCLE via brew/apt-get (manual on Windows)'
    )
    sub = parser.add_subparsers(dest='command')

    ex = sub.add_parser('extract', help='Full extraction + zip per FASTQ')
    ex.add_argument('-r','--ref-fasta',    required=True, help='Reference FASTA (single record)')
    ex.add_argument('-i','--input-fastq',  nargs='+', required=True, help='One or more FASTQ files')
    ex.add_argument('-s','--sites',        required=True, help='Comma-separated codon indices (1-based)')
    ex.add_argument('-o','--output-csv',   required=True, help='Detailed CSV output path')
    ex.add_argument('--ref-start',         type=int, default=1, help='Start codon index (1-based)')
    ex.add_argument('--ref-end',           type=int, help='End codon index (1-based)')
    ex.add_argument('--aligner',           choices=['muscle','mafft'], default='muscle', help='Alignment tool')
    ex.add_argument('-a','--algorithm',    choices=['align','super5'], default='super5', help='MUSCLE algorithm')
    ex.add_argument('--muscle-bin',        default='muscle', help='Path to MUSCLE v5 binary')
    ex.add_argument('--profile',           action='store_true', help='Reference-only alignment mode')
    ex.add_argument('--aligned-fasta',     help='Load existing MSA FASTA instead of running alignment')
    ex.add_argument('--aligned-fasta-out', help='Save generated MSA FASTA')
    ex.add_argument('--summary-csv',       action='store_true', help='Produce summary CSV of AA counts')
    ex.add_argument('--plots',             action='store_true', help='Produce plots')
    ex.add_argument('--plot-dir',          help='Directory to save plots')

    ps = sub.add_parser('post', help='Post-process existing detailed CSV')
    ps.add_argument('--detailed-csv', required=True, help='Existing detailed CSV')
    ps.add_argument('--summary-csv',           help='Summary CSV of AA counts per site')
    ps.add_argument('--plots',                 action='store_true', help='Produce plots')
    ps.add_argument('--plot-dir',              help='Directory for plots')

    return parser.parse_args()


def orient_read(read_seq: str, ref_seq: Seq) -> str:
    fwd = pairwise2.align.globalxx(ref_seq, read_seq, score_only=True)
    rev = pairwise2.align.globalxx(ref_seq, str(Seq(read_seq).reverse_complement()), score_only=True)
    return read_seq if fwd >= rev else str(Seq(read_seq).reverse_complement())


def run_muscle(in_fa, out_fa, args):
    if args.profile:
        cmd = [
            args.muscle_bin, '-profile',
            '-in1', args._ref_only_fa,
            '-in2', args._reads_fa,
            '-out', out_fa
        ]
    else:
        cmd = [args.muscle_bin, f'-{args.algorithm}', in_fa, '-output', out_fa]
    subprocess.run(cmd, check=True)


def run_mafft(in_fa, out_fa, args):
    if args.profile:
        cmd = [
            'mafft', '--addfragments', args._reads_fa,
            '--keeplength', '--anysymbol', args._ref_only_fa
        ]
        with open(out_fa, 'w') as wf:
            subprocess.run(cmd, check=True, stdout=wf)
    else:
        cmd = ['mafft', '--auto', in_fa]
        with open(out_fa, 'w') as wf:
            subprocess.run(cmd, check=True, stdout=wf)


def map_refpos_to_alncols(aln_seq):
    mapping, idx = {}, 0
    for col_idx, nt in enumerate(str(aln_seq.seq)):
        if nt != '-':
            mapping[idx] = col_idx
            idx += 1
    return mapping


def get_original_index(gapped_seq, col_idx):
    if gapped_seq[col_idx] == '-': return None
    return sum(1 for c in gapped_seq[:col_idx] if c != '-')


def post_process(detailed_csv, summary_csv=None, generate_plots=False, plot_dir=None):
    df = pd.read_csv(detailed_csv, engine='python', on_bad_lines='warn')
    sites = sorted(int(c.split('_')[1]) for c in df.columns if c.endswith('_codon'))
    df['combo'] = df.apply(
        lambda r: ''.join(r[f'site_{s}_aa'] for s in sites if pd.notna(r[f'site_{s}_aa'])),
        axis=1
    )

    if summary_csv:
        rows = []
        for s in sites:
            cnt = Counter(df[f'site_{s}_aa'].dropna()); tot = sum(cnt.values())
            for rank, (aa, ct) in enumerate(cnt.most_common(), 1):
                rows.append({'site': s, 'aa': aa, 'count': ct, 'rank': rank, 'freq': ct/tot})
        pd.DataFrame(rows).to_csv(summary_csv, index=False)

    if generate_plots:
        outd = plot_dir or os.getcwd()
        os.makedirs(outd, exist_ok=True)

        # Top5 AAs per site combined
        combined = {s: Counter(df[f'site_{s}_aa'].dropna()).most_common(5) for s in sites}
        plt.figure()
        for s, vals in combined.items():
            aas, cts = zip(*vals) if vals else ([],[])
            plt.bar([f"{s}:{aa}" for aa in aas], cts)
        plt.xticks(rotation=90)
        plt.title('Top5 AAs per site')
        plt.tight_layout()
        plt.savefig(os.path.join(outd, 'aa_top5_per_site.png'))
        plt.close()

        # Quality distribution per site
        qdata = []
        for s in sites:
            arr = df[[f'site_{s}_q{j}' for j in (1,2,3)]].apply(pd.to_numeric, errors='coerce').values.flatten()
            qdata.append([x for x in arr if not pd.isna(x)])
        plt.figure()
        plt.boxplot(qdata, tick_labels=sites)
        plt.xticks(rotation=90)
        plt.title('Quality distribution per site')
        plt.ylabel('PHRED')
        plt.tight_layout()
        plt.savefig(os.path.join(outd, 'quality_per_site.png'))
        plt.close()

        # Full-site combo top5 avg quality & counts
        df['avg_quality'] = df[[f'site_{s}_q{j}' for s in sites for j in (1,2,3)]].apply(
            pd.to_numeric, errors='coerce'
        ).mean(axis=1)
        full = df['combo'].dropna()[df['combo'].str.len() == len(sites)]
        top5_full = [c for c,_ in Counter(full).most_common(5)]

        if top5_full:
            cq = df[df['combo'].isin(top5_full)].groupby('combo')['avg_quality'].mean()
            plt.figure()
            plt.bar(cq.index, cq.values)
            plt.title('Avg quality by top5 full-site combos')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(outd, 'combo_vs_quality_top5_full.png'))
            plt.close()

            cnts = [Counter(full)[c] for c in top5_full]
            plt.figure()
            plt.bar(top5_full, cnts)
            plt.title('Counts of top5 full-site combos')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(outd, 'top5_combos_counts_full.png'))
            plt.close()

    print(f'Post-processing done for {detailed_csv}')


def extract_mode(args):
    fq = args.input_fastq[0]
    sites = [int(x) for x in args.sites.split(',')]

    if args.aligned_fasta:
        aln = AlignIO.read(args.aligned_fasta, 'fasta')
    else:
        refs = list(SeqIO.parse(args.ref_fasta, 'fasta'))
        if len(refs)!=1: sys.exit('ERROR: reference FASTA must have one record')
        ref = refs[0]
        total = len(ref.seq)//3
        start, end = args.ref_start, args.ref_end or total
        if start<1 or end>total or start> end: sys.exit('ERROR: invalid ref-start/end')
        region_seq = ref.seq[(start-1)*3:end*3]
        region_rec = SeqRecord(region_seq, id=ref.id, description='')

        reads = list(SeqIO.parse(fq, 'fastq'))
        if not reads: sys.exit(f'ERROR: no reads in {fq}')
        qual_dict = {r.id:r.letter_annotations.get('phred_quality',[]) for r in reads}
        oriented = [
            SeqRecord(Seq(orient_read(str(r.seq), region_seq)), id=r.id, description='')
            for r in reads
        ]

        with tempfile.TemporaryDirectory() as td:
            aln_file = os.path.join(td,'aligned.fa')
            if args.profile:
                ref_fa = os.path.join(td,'ref.fa'); SeqIO.write(region_rec, ref_fa, 'fasta')
                reads_fa = os.path.join(td,'reads.fa'); SeqIO.write(oriented, reads_fa, 'fasta')
                args._ref_only_fa, args._reads_fa = ref_fa, reads_fa
                in_fa = None
            else:
                in_fa = os.path.join(td,'in.fa')
                with open(in_fa,'w') as wf:
                    SeqIO.write(region_rec, wf, 'fasta')
                    SeqIO.write(oriented, wf, 'fasta')

            if args.aligner=='muscle':
                run_muscle(in_fa, aln_file, args)
            else:
                run_mafft(in_fa, aln_file, args)

            if args.aligned_fasta_out:
                shutil.copy(aln_file, args.aligned_fasta_out)
            aln = AlignIO.read(aln_file, 'fasta')

    # Extract codons
    ref_aln = next(a for a in aln if a.id == ref.id)
    pos2col = map_refpos_to_alncols(ref_aln)
    cols = [pos2col[(s-1)*3+o] for s in sites for o in (0,1,2)]

    det_csv = args.output_csv
    with open(det_csv,'w',newline='') as cf:
        w = csv.writer(cf)
        header = ['seq_id']
        for s in sites:
            header += [f'site_{s}_codon', f'site_{s}_aa'] + [f'site_{s}_q{j}' for j in (1,2,3)]
        w.writerow(header)
        for rec in aln:
            seqs = str(rec.seq)
            ql   = qual_dict.get(rec.id, [])
            row  = [rec.id]
            for i in range(len(sites)):
                sub = cols[3*i:3*i+3]
                cod, qs = [], []
                for cidx in sub:
                    nt = seqs[cidx]
                    if nt=='-':
                        cod.append('-'); qs.append('NaN')
                    else:
                        cod.append(nt)
                        oi = get_original_index(seqs, cidx)
                        qs.append(str(ql[oi]) if oi is not None and oi<len(ql) else 'NaN')
                codon = ''.join(cod)
                aa    = 'NaN' if '-' in codon else str(Seq(codon).translate())
                row  += [codon, aa] + qs
            w.writerow(row)
    print(f"Detailed CSV written to {det_csv}")

    # Post-process
    summary_flag = args.summary_csv
    plots_flag   = args.plots
    plot_dir     = args.plot_dir or f"{os.path.splitext(det_csv)[0]}_plots"
    if summary_flag or plots_flag:
        post_process(
            det_csv,
            summary_csv=(f"{os.path.splitext(det_csv)[0]}_summary.csv" if summary_flag else None),
            generate_plots=plots_flag,
            plot_dir=plot_dir
        )

    # Extra AA mutations outside specified sites
    extra_counts = Counter()
    # build reference AA map for codon indices outside `sites`
    ref_seq = str(ref_aln.seq)
    skip = set(s-1 for s in sites)
    for ci in range(len(ref_seq)//3):
        if ci in skip: continue
        cols_c = [pos2col[ci*3+o] for o in (0,1,2)]
        cod   = ''.join(ref_seq[c] for c in cols_c)
        if '-' in cod or len(cod)!=3: continue
        ref_aa = str(Seq(cod).translate())
        ref_codon_cols = cols_c
        # now scan each read
        for rec in aln:
            seqs = str(rec.seq)
            cod2 = ''.join(seqs[c] for c in ref_codon_cols)
            if '-' in cod2 or len(cod2)!=3: continue
            aa2 = str(Seq(cod2).translate())
            if aa2 != ref_aa:
                extra_counts[(ci+1, ref_aa, aa2)] += 1

    extra_csv = f"{os.path.splitext(det_csv)[0]}_extra_mutations.csv"
    with open(extra_csv,'w',newline='') as ef:
        ew = csv.writer(ef)
        ew.writerow(['site','ref_aa','alt_aa','count'])
        for (site, raa, aaa), cnt in extra_counts.most_common():
            ew.writerow([site, raa, aaa, cnt])
    print(f"Extra AA mutations CSV: {extra_csv}")

    if plots_flag:
        os.makedirs(plot_dir, exist_ok=True)
        top5 = extra_counts.most_common(5)
        labels = [f"{site}:{r}>{a}" for (site,r,a),_ in top5]
        values = [cnt for _,cnt in top5]
        plt.figure(figsize=(6,4))
        plt.bar(labels, values)
        plt.title('Top5 extra AA mutations')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'extra_aa_mutations_top5.png'))
        plt.close()

    # Zip results
    base = os.path.splitext(os.path.basename(det_csv))[0]
    zipn = f"{base}_analysed.zip"
    with zipfile.ZipFile(zipn, 'w') as zf:
        zf.write(det_csv)
        if summary_flag: zf.write(f"{base}_summary.csv")
        zf.write(extra_csv)
        if plots_flag:
            for root,_,files in os.walk(plot_dir):
                for fn in files:
                    arc = os.path.join(os.path.basename(plot_dir), fn)
                    zf.write(os.path.join(root, fn), arcname=arc)

    # Clean up
    os.remove(det_csv)
    os.remove(extra_csv)
    if summary_flag: os.remove(f"{base}_summary.csv")
    if plots_flag: shutil.rmtree(plot_dir)

    print(f"Created {zipn}")


def main():
    args = parse_args()

    if getattr(args, 'package', False):
        # install Python packages
        subprocess.run([sys.executable, '-m', 'pip', 'install',
                        'biopython', 'pandas', 'matplotlib'], check=True)
        os_name = platform.system()
        if os_name == 'Darwin':
            subprocess.run(['brew', 'install', 'mafft', 'muscle'], check=True)
        elif os_name == 'Linux':
            subprocess.run(['sudo', 'apt-get', 'update'], check=True)
            subprocess.run(['sudo', 'apt-get', 'install', '-y', 'mafft', 'muscle'], check=True)
        elif os_name == 'Windows':
            print("Please install MAFFT & MUSCLE manually per the Windows instructions.")
        sys.exit(0)

    if args.command == 'extract':
        extract_mode(args)
    else:
        post_process(
            args.detailed_csv,
            summary_csv=args.summary_csv,
            generate_plots=args.plots,
            plot_dir=args.plot_dir
        )


if __name__ == '__main__':
    main()