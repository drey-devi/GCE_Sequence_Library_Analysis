#!/usr/bin/env python3
"""
Script: extract_codons_qc.py

Two subcommands:

1. extract – Full extraction mode:
   • Optionally orient reads to the reference (auto-detect forward/reverse)
   • Align (or load) MSA, extract codons → detailed CSV
   • Choose aligner: MUSCLE or MAFFT
   • For either aligner, run full MSA or reference-only mode
   • Optional summary CSV & plots
   • Supports multiple FASTQ inputs: outputs are zipped per FASTQ file as <basename>_analysed.zip

2. post    – Post-process-only mode:
   Take an existing detailed CSV → summary CSV + combined site plots + full-site combo graphs

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
    ex.add_argument('-r', '--ref-fasta', required=True, help='Reference FASTA (single record)')
    ex.add_argument('-i', '--input-fastq', nargs='+', required=True, help='One or more FASTQ files')
    ex.add_argument('-s', '--sites', required=True, help='Comma-separated codon indices (1-based)')
    ex.add_argument('-o', '--output-csv', required=True, help='Detailed output CSV')
    ex.add_argument('--ref-start', type=int, default=1, help='1-based codon index where region starts')
    ex.add_argument('--ref-end', type=int, help='1-based codon index where region ends')
    ex.add_argument('--aligner', choices=['muscle','mafft'], default='muscle', help='Alignment tool')
    ex.add_argument('-a', '--algorithm', choices=['align','super5'], default='super5', help='MUSCLE v5 algorithm')
    ex.add_argument('--muscle-bin', default='muscle', help='Path to MUSCLE v5 binary')
    ex.add_argument('--profile', action='store_true', help='Reference-only alignment (profile/addfragments)')
    ex.add_argument('--aligned-fasta', help='Existing MSA FASTA to load')
    ex.add_argument('--aligned-fasta-out', help='Write generated MSA FASTA')
    ex.add_argument('--summary-csv', action='store_true', help='Generate summary CSV per input')
    ex.add_argument('--plots', action='store_true', help='Generate plots per input')

    # post subcommand
    ps = subparsers.add_parser('post', help='Post-process existing detailed CSV')
    ps.add_argument('--detailed-csv', required=True, help='Existing detailed CSV')
    ps.add_argument('--summary-csv', help='Summary CSV of AA counts per site')
    ps.add_argument('--plots', action='store_true', help='Generate plots')
    ps.add_argument('--plot-dir', help='Directory for plots')

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
    df['combo'] = df.apply(lambda r: ''.join(r[f'site_{s}_aa'] for s in sites if pd.notna(r[f'site_{s}_aa'])), axis=1)

    # summary CSV
    if summary_csv:
        rows=[]
        for s in sites:
            cnt=Counter(df[f'site_{s}_aa'].dropna()); tot=sum(cnt.values())
            for rk,(aa,ct) in enumerate(cnt.most_common(),1): rows.append({'site':s,'aa':aa,'count':ct,'rank':rk,'freq':ct/tot})
        pd.DataFrame(rows).to_csv(summary_csv,index=False)

    # plots
    if generate_plots:
        outd = plot_dir or os.getcwd(); os.makedirs(outd, exist_ok=True)
        # Combined Top5 AA per site
        combined = {s: Counter(df[f'site_{s}_aa'].dropna()).most_common(5) for s in sites}
        plt.figure()
        for s, vals in combined.items():
            aas,cts = zip(*vals) if vals else ([], [])
            plt.bar([f"{s}:{aa}" for aa in aas], cts)
        plt.xticks(rotation=90); plt.title('Top5 AAs per site'); plt.tight_layout()
        plt.savefig(os.path.join(outd,'aa_top5_per_site.png')); plt.close()

        # Combined quality per site
        qdata=[]
        for s in sites:
            arr=df[[f'site_{s}_q{j}' for j in (1,2,3)]].apply(pd.to_numeric,errors='coerce').values.flatten()
            qdata.append([x for x in arr if not pd.isna(x)])
        plt.figure(); plt.boxplot(qdata, tick_labels=sites)
        plt.xticks(rotation=90); plt.title('Quality distribution per site'); plt.ylabel('PHRED'); plt.tight_layout()
        plt.savefig(os.path.join(outd,'quality_per_site.png')); plt.close()

        # Full-site combo graphs
        # Avg quality by top5 full combos
        df['avg_quality'] = df[[f'site_{s}_q{j}' for s in sites for j in (1,2,3)]].apply(pd.to_numeric,errors='coerce').mean(axis=1)
        full = df['combo'].dropna()[df['combo'].str.len()==len(sites)]
        top5_full=[c for c,_ in Counter(full).most_common(5)]
        if top5_full:
            cq = df[df['combo'].isin(top5_full)].groupby('combo')['avg_quality'].mean()
            plt.figure(); plt.bar(cq.index,cq.values); plt.title('Avg quality by top5 full-site combos')
            plt.xticks(rotation=45,ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outd,'combo_vs_quality_top5_full.png')); plt.close()
            # Counts of top5 full combos
            cnts=[Counter(full)[c] for c in top5_full]
            plt.figure(); plt.bar(top5_full,cnts); plt.title('Counts of top5 full-site combos')
            plt.xticks(rotation=45,ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outd,'top5_combos_counts_full.png')); plt.close()

    print(f'Post-processing done for {detailed_csv}')


def extract_mode(args):
    codon_sites=[int(x) for x in args.sites.split(',')]
    # ... (rest of extract_mode unchanged)

    # At end of extract_mode, call post_process and zip as before.


def main():
    args = parse_args()
    if args.command=='extract': extract_mode(args)
    else: post_process(args.detailed_csv, summary_csv=args.summary_csv, generate_plots=args.plots, plot_dir=args.plot_dir)

if __name__=='__main__': main()
