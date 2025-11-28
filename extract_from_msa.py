#!/usr/bin/env python3
import argparse
import csv
import sys
from Bio import AlignIO

def parse_args():
    p = argparse.ArgumentParser(
        description="Extract codons at given sites from an existing MSA"
    )
    p.add_argument(
        "-a", "--alignment", required=True,
        help="Input alignment file (FASTA/MFA) from MUSCLE v5"
    )
    p.add_argument(
        "-s", "--sites", required=True,
        help="Comma‑separated list of codon indices (1‑based). e.g. '1,5,10'"
    )
    p.add_argument(
        "-r", "--reference", default=None,
        help=(
            "ID of the reference sequence in the MSA to map codon positions against. "
            "If omitted, uses the first record."
        )
    )
    p.add_argument(
        "-o", "--output-csv", required=True,
        help="Path to output CSV"
    )
    return p.parse_args()

def map_refpos_to_alncols(aligned_ref_seq):
    """
    Given the reference as a gapped string, return a dict
    mapping its ungapped positions (0-based) to alignment columns.
    """
    mapping = {}
    ref_pos = 0
    for col_idx, base in enumerate(aligned_ref_seq):
        if base != "-":
            mapping[ref_pos] = col_idx
            ref_pos += 1
    return mapping

def main():
    args = parse_args()

    # parse codon sites → list of 0-based base indices
    codon_sites = [int(x) for x in args.sites.split(",")]
    base_indices = []
    for site in codon_sites:
        start = (site - 1) * 3
        base_indices.extend([start, start+1, start+2])

    # read the alignment
    aln = AlignIO.read(args.alignment, "fasta")
    if not aln:
        sys.exit("ERROR: No sequences found in alignment.")

    # select reference record
    if args.reference:
        ref_records = [rec for rec in aln if rec.id == args.reference]
        if not ref_records:
            sys.exit(f"ERROR: reference ID '{args.reference}' not found")
        ref = ref_records[0]
    else:
        ref = aln[0]

    # build mapping: ungapped ref pos → alignment column
    pos2col = map_refpos_to_alncols(str(ref.seq))

    # translate our desired base indices → alignment columns
    try:
        aln_cols = [pos2col[i] for i in base_indices]
    except KeyError as e:
        sys.exit(f"ERROR: reference too short to map base position {e}")

    # write CSV
    n_sites = len(codon_sites)
    with open(args.output_csv, "w", newline="") as outf:
        writer = csv.writer(outf)
        writer.writerow(["seq_id"] + [f"site_{s}" for s in codon_sites])

        for rec in aln:
            seq = str(rec.seq)
            codons = []
            for i in range(n_sites):
                cols = aln_cols[3*i : 3*i+3]
                codon = "".join(seq[c] for c in cols)
                codons.append(codon)
            writer.writerow([rec.id] + codons)

    print(f"Done. Results written to {args.output_csv}")

if __name__ == "__main__":
    main()
