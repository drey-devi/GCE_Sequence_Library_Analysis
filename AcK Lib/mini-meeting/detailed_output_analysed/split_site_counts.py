#!/usr/bin/env python3
"""
split_site_counts.py

Given a detailed CSV with columns like site_121_aa, site_121_codon, etc.,
produce one CSV per codon site listing each amino acid and its count.

Usage:
    python split_site_counts.py \
      --input detailed.csv \
      [--outdir counts_by_site]

This will create files:
    counts_by_site/site_121_counts.csv
    counts_by_site/site_125_counts.csv
    ...
"""
import argparse
import os
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(
        description="Split detailed.csv into per-site AA count CSVs"
    )
    p.add_argument(
        "--input", "-i", required=True,
        help="Path to detailed CSV (with site_<n>_aa columns)"
    )
    p.add_argument(
        "--outdir", "-o", default="counts_by_site",
        help="Directory to write per-site count CSVs"
    )
    return p.parse_args()

def main():
    args = parse_args()
    df = pd.read_csv(args.input, engine="python")
    # find all *_aa columns
    aa_cols = [c for c in df.columns if c.endswith("_aa")]
    # extract unique site numbers
    sites = sorted({int(c.split("_")[1]) for c in aa_cols})
    os.makedirs(args.outdir, exist_ok=True)

    for site in sites:
        col = f"site_{site}_aa"
        if col not in df:
            continue
        # count non-null AAs
        counts = df[col].dropna().value_counts().reset_index()
        counts.columns = ["aa", "count"]
        out_path = os.path.join(args.outdir, f"site_{site}_counts.csv")
        counts.to_csv(out_path, index=False)
        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()