
# Lung Hi-C Analysis

### Reproducible analysis of 3D genome organization in healthy and COPD human lung tissue

This repository contains scripts for processing and analyzing Hi-C sequencing data from human lung tissue, including in-house samples and publicly available datasets from the ENCODE consortium.

The workflow begins with downloading public Hi-C contact matrices and processing in-house sequencing data using Juicer.

## 1. Download publicly available ENCODE lung Hi-C datasets

Six human lung Hi-C datasets were selected from [ENCODE](https://www.encodeproject.org/). These experiments represent four individual donors.

| ENCODE accession | GEO accession | Donor | Lung tissue |
|---|---|---|---|
| [ENCSR305GRJ](https://www.encodeproject.org/experiments/ENCSR305GRJ/) | [GSE237368](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237368) | W76, male, 40 | Left lung |
| [ENCSR749PID](https://www.encodeproject.org/experiments/ENCSR749PID/) | [GSE237921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237921) | W78, male, 60 | Right lower lobe |
| [ENCSR673UZF](https://www.encodeproject.org/experiments/ENCSR673UZF/) | [GSE237897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237897) | W78, male, 60 | Right upper lobe |
| [ENCSR446NJD](https://www.encodeproject.org/experiments/ENCSR446NJD/) | [GSE237470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237470) | W73, female, 59 | Left lower lobe |
| [ENCSR068FXF](https://www.encodeproject.org/experiments/ENCSR068FXF/) | [GSE237333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237333) | W78, male, 60 | Left lower lobe |
| [ENCSR892BZA](https://www.encodeproject.org/experiments/ENCSR892BZA/) | [GSE237904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237904) | W62, female, 16 | Left lung |

The GRCh38 mapping-quality-thresholded `.hic` contact matrices are retrieved using the download script:

`scripts/encode/download_encode_hic.sh`

```bash
bash scripts/encode/download_encode_hic.sh
```
## 2. Hi-C processing using Juicer

In-house human lung Hi-C sequencing data were processed using [Juicer v1.6](https://github.com/aidenlab/juicer) on a SLURM-based high-performance computing cluster.

The processing configuration used:

| Parameter | Setting |
|---|---|
| Reference genome | GRCh38 (hg38) |
| Hi-C protocol | Arima |
| Alignment software | BWA |
| Job scheduler | SLURM |
| CPUs requested for submission job | 24 |
| Memory requested for submission job | 400 GB |

Juicer performs read alignment, processing and filtering of Hi-C contact pairs, duplicate removal, and generation of multiresolution `.hic` contact matrices.

The SLURM submission script is provided in:

`scripts/juicer/run_juicer.slurm`

The script requires paths to the Juicer installation, reference genome, chromosome sizes, restriction-site file, and sample input directories to be configured for the user's computing environment.

## 3. Contact-matrix integration

The ENCODE contact matrices and in-house Juicer-generated matrices are intended for integration into a combined human lung Hi-C contact map.

The merging workflow and downstream loop-calling analysis will be documented after validation.

## Software

- [Juicer v1.6](https://github.com/aidenlab/juicer)
- BWA
- SLURM
- `wget`
- Python 3 (for ENCODE metadata retrieval)

## Status

This repository is under active development. Scripts and instructions will be updated as each analysis step is validated.
