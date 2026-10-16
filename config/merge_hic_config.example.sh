#!/usr/bin/env bash

# Software
HICTK_IMAGE="/path/to/hictk_2.2.0.sif"
JUICER_TOOLS="/path/to/juicer_tools_2.20.00.jar"

# Directories
OUTPUT_DIR="/path/to/merged_hic"
TMP_DIR="/path/to/scratch/tmp"

# Input .hic file list
MERGE_SAMPLES="${REPO_DIR}/config/merge_samples.tsv"

# Resolutions (bp)
MERGE_RESOLUTION=5000
ZOOMIFY_RESOLUTIONS=(5000 10000 20000 100000 500000 1000000)

# Number of threads
THREADS=24

# Output filenames
MERGED_FILE="merged_5kb.hic"
ZOOMIFIED_FILE="merged_multires.hic"
FINAL_FILE="merged_multires_KR.hic"

