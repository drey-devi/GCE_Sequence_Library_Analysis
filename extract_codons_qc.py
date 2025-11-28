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

2. post    – Post-process-only mode:
   Take an existing detailed CSV → summary CSV + plots

Run `extract_codons_qc.py <command> --help` for subcommand-specific options.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
import shutil
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        prog='extract_codons_qc.py',
        description='''
A tool to extract codons from an MSA or FASTQ→reference alignment and post-process results.

Subcommands:
  extract : Orient reads, align (MUSCLE/MAFFT), extract codons, and optionally summary and plots.
  post    : Run summary and plots on an existing detailed CSV.
''',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:

# Full extraction with MUSCLE full MSA
./extract_codons_qc.py extract \
  -r ref_ma_gne.fa -i reads.fastq -s 5,10,15 -o detailed.csv \
  --aligner muscle --algorithm super5 --summary-csv summary.csv --plots --plot-dir plots

# Reference-only alignment with MUSCLE profile
./extract_codons_qc.py extract \
  -r ref_ma_gne.fa -i reads.fastq -s 5,10,15 -o detailed.csv \
  --aligner muscle --profile --muscle-bin muscle5 \
  --aligned-fasta-out aln.fa --summary-csv summary.csv --plots --plot-dir plots

# Full MSA with MAFFT
./extract_codons_qc.py extract \
  -r ref_ma_gne.fa -i reads.fastq -s 5,10,15 -o detailed.csv \
  --aligner mafft --summary-csv summary.csv --plots --plot-dir plots

# Reference-only with MAFFT addfragments
./extract_codons_qc.py extract \
  -r ref_ma_gne.fa -i reads.fastq -s 5,10,15 -o detailed.csv \
  --aligner mafft --profile --summary-csv summary.csv --plots --plot-dir plots

# Post-process an existing detailed CSV
./extract_codons_qc.py post \
  --detailed-csv detailed.csv --summary-csv summary.csv --plots --plot-dir plots
'''
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # extract subcommand
    ex = subparsers.add_parser('extract', help='Full extraction: orient + align + extract codons')
    ex.add_argument('-r', '--ref-fasta',     required=True, help='Reference FASTA (single-record)')
    ex.add_argument('-i', '--input-fastq',   required=True, help='Input FASTQ file')
    ex.add_argument('-s', '--sites',         required=True, help='Comma-separated codon indices (1-based)')
    ex.add_argument('-o', '--output-csv',    required=True, help='Detailed output CSV')
    ex.add_argument('--ref-start', type=int, default=1, help='1-based codon where region starts')
    ex.add_argument('--ref-end',   type=int,             help='1-based codon where region ends')
    ex.add_argument('--aligner',   choices=['muscle','mafft'], default='muscle', help='Alignment tool to use')
    ex.add_argument('-a', '--algorithm', choices=['align','super5'], default='super5',
                    help='MUSCLE v5 algorithm (only for MUSCLE)')
    ex.add_argument('--muscle-bin', default='muscle', help='Path to MUSCLE v5 binary')
    ex.add_argument('--profile',   action='store_true', help='Reference-only alignment (profile/addfragments)')
    ex.add_argument('--min-quality', type=int, default=0, help='PHRED threshold for low-quality')
    ex.add_argument('--aligned-fasta',     help='Existing MSA FASTA to load instead of alignment')
    ex.add_argument('--aligned-fasta-out', help='Output path for generated MSA FASTA')
    ex.add_argument('--summary-csv',       help='Optional summary CSV of AA counts per site')
    ex.add_argument('--plots', action='store_true', help='Generate summary plots')
    ex.add_argument('--plot-dir',          help='Directory to save plots')

    # post subcommand
    ps = subparsers.add_parser('post', help='Post-process existing detailed CSV')
    ps.add_argument('--detailed-csv', required=True, help='Existing detailed CSV with codon data')
    ps.add_argument('--summary-csv',           help='Summary CSV of AA counts per site')
    ps.add_argument('--plots', action='store_true', help='Generate summary plots')
    ps.add_argument('--plot-dir',              help='Directory to save plots')

    return parser.parse_args()


def orient_read(read_seq: str, ref_seq: Seq) -> str:
    """
    Align read_seq vs. ref_seq in both orientations and return best orientation.
    """
    fwd = pairwise2.align.globalxx(ref_seq, read_seq, score_only=True)
    rev = pairwise2.align.globalxx(ref_seq, str(Seq(read_seq).reverse_complement()), score_only=True)
    return read_seq if fwd >= rev else str(Seq(read_seq).reverse_complement())


def run_muscle(in_fa, out_fa, args):
    if args.profile:
        cmd = [args.muscle_bin, '-profile',
               '-in1', args._ref_only_fa,
               '-in2', args._reads_fa,
               '-out', out_fa]
    else:
        cmd = [args.muscle_bin, f'-{args.algorithm}', in_fa, '-output', out_fa]
    subprocess.run(cmd, check=True)


def run_mafft(in_fa, out_fa, args):
    if args.profile:
        cmd = ['mafft', '--addfragments', args._reads_fa,
               '--keeplength', '--anysymbol', args._ref_only_fa]
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


def translate_codon(codon):
    if pd.isna(codon) or '-' in codon or len(codon) != 3: return None
    return str(Seq(codon).translate())


def post_process(detailed_csv, summary_csv=None, generate_plots=False, plot_dir=None):
    df = pd.read_csv(
        detailed_csv,
        sep=',',
        engine='python',
        on_bad_lines='warn'
    )
    sites = sorted(int(c.split('_')[1]) for c in df.columns if c.endswith('_codon'))
    df['combo'] = df.apply(lambda r: ''.join(r[f"site_{x}_aa"] for x in sites if pd.notna(r[f"site_{x}_aa"])), axis=1)

    if summary_csv:
        rows=[]
        for s in sites:
            aa_list = df[f'site_{s}_aa'].dropna(); cnts = Counter(aa_list); tot = sum(cnts.values())
            for rank,(aa,cnt) in enumerate(cnts.most_common(),1): rows.append({'site':s,'aa':aa,'count':cnt,'rank':rank,'freq':cnt/tot})
        pd.DataFrame(rows).to_csv(summary_csv,index=False)
        combo_csv = os.path.splitext(summary_csv)[0] + '_combos.csv'
        combos = df[['seq_id','combo']].copy(); combos['count'] = combos['combo'].map(Counter(combos['combo']))
        combos.to_csv(combo_csv,index=False); print(f"Wrote summary: {summary_csv} and combos: {combo_csv}")

    if generate_plots:
        outd = plot_dir or os.getcwd(); os.makedirs(outd, exist_ok=True)
        for s in sites:
            aa_list = df[f'site_{s}_aa'].dropna(); top5 = Counter(aa_list).most_common(5)
            if top5:
                aas,cts=zip(*top5); plt.figure(); plt.bar(aas,cts)
                plt.title(f'Top5 AAs at site {s}'); plt.tight_layout()
                plt.savefig(os.path.join(outd,f'top5_aa_{s}.png')); plt.close()
        for s in sites:
            col=f'site_{s}_status';
            if col in df:
                plt.figure(); df[col].value_counts().plot(kind='bar')
                plt.title(f'Status distribution site {s}'); plt.tight_layout()
                plt.savefig(os.path.join(outd,f'status_{s}.png')); plt.close()
        for s in sites:
            qcols=[f'site_{s}_q{j}' for j in(1,2,3)]; qdf=df[qcols].apply(pd.to_numeric,errors='coerce').dropna()
            plt.figure(); plt.boxplot([qdf[c] for c in qcols],labels=qcols)
            plt.title(f'Quality scores at site {s}'); plt.ylabel('PHRED'); plt.tight_layout()
            plt.savefig(os.path.join(outd,f'quality_box_site_{s}.png')); plt.close()
        # Top5 full-site combos vs quality
        qcols_all=[f"site_{s}_q{j}" for s in sites for j in (1,2,3)]
        df['avg_quality']=df[qcols_all].apply(pd.to_numeric,errors='coerce').mean(axis=1)
        full_combos = df['combo'].dropna()[df['combo'].str.len()==len(sites)]
        top5_full=[c for c,_ in Counter(full_combos).most_common(5)]
        if top5_full:
            cq=df[df['combo'].isin(top5_full)].groupby('combo')['avg_quality'].mean()
            plt.figure(figsize=(8,4)); plt.bar(cq.index,cq.values)
            plt.title('Average quality by top 5 full-site AA combos'); plt.ylabel('Avg PHRED')
            plt.xticks(rotation=45,ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outd,'combo_vs_quality_top5_full.png')); plt.close()
        if top5_full:
            cnts=[Counter(full_combos)[c] for c in top5_full]
            plt.figure(figsize=(6,4)); plt.bar(top5_full,cnts)
            plt.title('Counts of top 5 full-site AA combos'); plt.xticks(rotation=45,ha='right'); plt.tight_layout()
            plt.savefig(os.path.join(outd,'top5_combos_counts_full.png')); plt.close()
        print(f"Plots saved in {outd}")


def extract_mode(args):
    codon_sites=[int(x) for x in args.sites.split(',')]
    if args.aligned_fasta:
        aln=AlignIO.read(args.aligned_fasta,'fasta')
    else:
        refs=list(SeqIO.parse(args.ref_fasta,'fasta'))
        if len(refs)!=1: sys.exit('ERROR: single-record FASTA required')
        ref=refs[0]; total=len(ref.seq)//3
        start, end = args.ref_start, args.ref_end or total
        if start<1 or end>total or start> end: sys.exit('ERROR: invalid ref-start/end')
        region_seq=ref.seq[(start-1)*3:end*3]; region_ref=SeqRecord(region_seq,id=ref.id,description='region')
        reads=list(SeqIO.parse(args.input_fastq,'fastq'))
        if not reads: sys.exit('ERROR: no reads found')
        qual_dict={r.id:r.letter_annotations.get('phred_quality',[]) for r in reads}
        oriented_reads=[]
        for r in reads:
            best=orient_read(str(r.seq), region_seq)
            oriented_reads.append(SeqRecord(Seq(best),id=r.id,description=''))
        with tempfile.TemporaryDirectory() as td:
            in_fa = None
            aln_file=os.path.join(td,'aligned.fa')
            if args.profile:
                ref_fa=os.path.join(td,'ref.fa'); SeqIO.write(region_ref,ref_fa,'fasta')
                reads_fa=os.path.join(td,'reads.fa'); SeqIO.write(oriented_reads,reads_fa,'fasta')
                args._ref_only_fa, args._reads_fa = ref_fa, reads_fa
            else:
                in_fa=os.path.join(td,'in.fa')
                with open(in_fa,'w') as wf:
                    SeqIO.write(region_ref,wf,'fasta'); SeqIO.write(oriented_reads,wf,'fasta')
            if args.aligner=='muscle':
                run_muscle(in_fa, aln_file, args)
            else:
                run_mafft(in_fa, aln_file, args)
            if args.aligned_fasta_out: shutil.copy(aln_file,args.aligned_fasta_out)
            aln=AlignIO.read(aln_file,'fasta')
    if args.aligned_fasta:
        reads=list(SeqIO.parse(args.input_fastq,'fastq'))
        qual_dict={r.id:r.letter_annotations.get('phred_quality',[]) for r in reads}
    ref_aln=next(a for a in aln if a.id==ref.id)
    pos2col=map_refpos_to_alncols(ref_aln)
    base_idxs=[]
    for s in codon_sites:
        b0=(s-1)*3; base_idxs+=[b0,b0+1,b0+2]
    aln_cols=[pos2col[i] for i in base_idxs]
    with open(args.output_csv,'w',newline='') as outf:
        w=csv.writer(outf); hdr=['seq_id']
        for s in codon_sites: hdr+=[f'site_{s}_codon',f'site_{s}_aa']+[f'site_{s}_q{j}' for j in(1,2,3)]
        w.writerow(hdr)
        for rec in aln:
            seq_str=str(rec.seq); qlist=qual_dict.get(rec.id,[])
            row=[rec.id]; aa_list=[]
            for i,s in enumerate(codon_sites):
                cols=aln_cols[3*i:3*i+3]; codon_chars,quals=[],[]
                for c in cols:
                    nt=seq_str[c]
                    if nt=='-': codon_chars.append('-'); quals.append('NaN')
                    else: codon_chars.append(nt); oi=get_original_index(seq_str,c);
                    quals.append(str(qlist[oi]) if oi is not None and oi<len(qlist) else 'NaN')
                codon=''.join(codon_chars); aa='NaN'
                if '-' not in codon: aa=str(Seq(codon).translate())
                row+=[codon,aa]+quals; aa_list.append(aa)
            w.writerow(row)
    print(f"Detailed CSV written to {args.output_csv}")
    if args.summary_csv or args.plots:
        post_process(detailed_csv=args.output_csv,
                     summary_csv=args.summary_csv,
                     generate_plots=args.plots,
                     plot_dir=args.plot_dir)


def main():
    args=parse_args()
    if args.command=='extract': extract_mode(args)
    else: post_process(detailed_csv=args.detailed_csv,
                       summary_csv=args.summary_csv,
                       generate_plots=args.plots,
                       plot_dir=args.plot_dir)

if __name__=='__main__':
    main()
