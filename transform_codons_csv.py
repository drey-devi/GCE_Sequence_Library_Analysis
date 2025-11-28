#!/usr/bin/env python3
"""
Script: transform_codons_csv.py

Reads an existing codons CSV with columns:
  seq_id, site_<n>, ...

Each `site_<n>` entry is "CODON|STATUS" (e.g. "GCA|OK" or "---|MISSING").
This script splits each into:
  site_<n>_codon   (string; missing→"NaN")
  site_<n>_status  (string; missing→"NaN")
Drops the original combined columns.

Usage:
  python3 transform_codons_csv.py \
    --input-csv path/to/codons.csv \
    --output-csv path/to/codons_split.csv
"""
import argparse
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description="Split CODON|STATUS fields into separate columns"
    )
    parser.add_argument(
        "--input-csv", "-i", required=True,
        help="Path to the original codons CSV"
    )
    parser.add_argument(
        "--output-csv", "-o", required=True,
        help="Path for the transformed CSV"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input_csv)

    # Identify all 'site_<number>' columns
    site_cols = [c for c in df.columns if re.match(r'^site_\d+$', c)]

    # Start new DataFrame with seq_id
    new_df = df[['seq_id']].copy()

    for site in site_cols:
        # Split "CODON|STATUS"
        parts = df[site].astype(str).str.split(r"\|", n=1, expand=True)
        codon = parts[0].replace({'---': 'NaN'}).fillna('NaN')
        status = parts[1].fillna('NaN')

        new_df[f"{site}_codon"] = codon
        new_df[f"{site}_status"] = status

    new_df.to_csv(args.output_csv, index=False)
    print(f"Transformed CSV written to {args.output_csv}")

if __name__ == '__main__':
    main()
