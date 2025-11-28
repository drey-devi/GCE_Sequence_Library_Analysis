# GCE_Sequence_Library_Analysis
Library Sequencing Analysis Pipeline

This repository contains a Python-based workflow for processing next-generation sequencing data from Genetic Code Expansion (GCE) variant libraries. The script performs orientation correction, multiple sequence alignment, codon extraction, PHRED scoring, amino acid translation, and summary/visualization of library composition. It is designed for high-throughput analysis of directed evolution libraries and site-specific mutational scans.

Overview

The pipeline takes a reference sequence (FASTA) and one or more sequencing files (FASTQ), aligns all reads to the reference, extracts user-specified codons, evaluates base quality, and generates both residue-level and full-protein-level summaries. Outputs include CSV tables of extracted codons and amino acids, reconstructed protein variants, alignment files, and graphical summaries.

Input Requirements

Required Inputs

Reference file (FASTA) – expected wild-type or designed reference sequence

Library sequencing data (FASTQ) – can accept single or multiple files

Optional Inputs

Precomputed alignment (.fasta)

Previously generated detailed CSV for post-processing

Generated Outputs
Output Type	Description
Codon Isolated CSV	Extracted codons at specified sites, associated PHRED scores, and corresponding amino acids
Amino Acid Summary CSV	Per-site amino-acid counts and frequencies
Protein Summary CSV	Reconstructed full-length protein variants and counts
Alignment File (FASTA)	MSA produced using MAFFT or MUSCLE5
Plots (PNG)	Per-site amino acid distributions, protein variant distributions, and QC visualizations

These outputs support downstream selection analysis, fidelity assessment, or comparison of pre- and post-selection libraries.

Pipeline Structure

Read Orientation (Biopython)
Ensures all reads are converted to forward direction before alignment.

Multiple Sequence Alignment
Users may choose between:

MAFFT (recommended; supports fast reference-profile alignment)

MUSCLE5 (--align or --super5 modes)

Codon Extraction and QC
Codons at user-specified positions are extracted from aligned reads, with PHRED quality filtering if requested.

Amino Acid Translation and Summaries
Outputs per-site amino acid distributions and complete protein variants.

Visualization
Generates publication-quality PNG plots for rapid interpretation.

Example Commands
MAFFT – Full Alignment
python extract_codons_qc.py extract \
  -r ref.fa \
  -i input.fastq \
  -s 121,125,126,129,168,206,223 \
  -o detailed_mafft_full.csv \
  --aligner mafft \
  --aligned-fasta-out aln_mafft_full.fa \
  --summary-csv aa_mafft_full_summary.csv \
  --plots \
  --plot-dir plots_mafft_full

MAFFT Profile Alignment (Fastest)
python extract_codons_qc.py extract \
  -r ref.fa \
  -i input.fastq \
  -s 121,125,126,129,168,206,223 \
  -o detailed_mafft_profile.csv \
  --aligner mafft \
  --profile \
  --aligned-fasta-out aln_mafft_profile.fa \
  --summary-csv aa_mafft_profile_summary.csv \
  --plots \
  --plot-dir plots_mafft_profile

MUSCLE5 (High Accuracy)
python extract_codons_qc.py extract \
  -r ref.fa \
  -i input.fastq \
  -s 121,125,126,129,168,206,223 \
  -o detailed_muscle_full.csv \
  --aligner muscle \
  --algorithm super5 \
  --aligned-fasta-out aln_muscle_full.fa \
  --summary-csv aa_muscle_full_summary.csv \
  --plots \
  --plot-dir plots_muscle_full

Reuse an Existing Alignment
python extract_codons_qc.py extract \
  -r ref.fa \
  -i reads.fastq \
  -s 5,10,15 \
  -o detailed2.csv \
  --aligned-fasta aln.fa \
  --summary-csv summary2.csv \
  --plots \
  --plot-dir plots2

Post-Process an Existing CSV
python extract_codons_qc.py post \
  --detailed-csv detailed.csv \
  --summary-csv summary.csv \
  --plots \
  --plot-dir plots

Command-line Arguments

Core Options

-r : Reference FASTA path

-i : Input FASTQ files

-s : Codon positions (comma-separated)

-o : Output detailed CSV

--aligner {mafft,muscle}

--profile : Reference-only profile alignment

--summary-csv : Output summary CSV

--aligned-fasta-out : Write MSA output

--plots / --plot-dir : Enable and direct plot output

Additional Options

--min-quality : PHRED threshold

--algorithm : MUSCLE v5 mode (align, super5)

--ref-start, --ref-end : Specify codon boundaries in reference

Software Dependencies

Biopython – read orientation, FASTA/FASTQ parsing, alignment parsing

MAFFT – MSA for most use cases; supports fast profile-based alignment

MUSCLE v5 – high-accuracy MSA using clustering and Super5 algorithms
