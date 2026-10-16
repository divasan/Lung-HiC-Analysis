#!/usr/bin/env bash

# Juicer 1.6 installation
JUICER_DIR="/path/to/juicer-1.6/SLURM"
JUICER_EXECUTABLE="${JUICER_DIR}/scripts/juicer.sh"

# Sample processing directories
PROJECT_DIR="/path/to/hic_samples"

# GRCh38 reference
GENOME_ID="hg38"
REFERENCE_DIR="/path/to/reference"
REFERENCE_GENOME="${REFERENCE_DIR}/hg38.fa"
CHROM_SIZES="${REFERENCE_DIR}/hg38.chrom.sizes"

# Arima restriction sites
RESTRICTION_ENZYME="Arima"
RESTRICTION_SITES="${REFERENCE_DIR}/hg38_Arima.txt"

# Juicer internal SLURM queues
QUEUE="standard"
LONG_QUEUE="standard"

