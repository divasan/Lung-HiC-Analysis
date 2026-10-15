#!/usr/bin/env bash

# hictk installation
HICTK_IMAGE="/path/to/hictk_2.2.0.sif"

# Merge directories
OUTPUT_DIR="/path/to/merged_hic"
TMP_DIR="/path/to/scratch/tmp"

# Input file list
MERGE_SAMPLES="${REPO_DIR}/config/merge_samples.tsv"

# Merge resolution (bp)
MERGE_RESOLUTION=5000

# Number of threads
THREADS=24

# Output file
MERGED_FILE="merged_5kb.hic"
