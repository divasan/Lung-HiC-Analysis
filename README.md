
# Lung Hi-C Analysis

### 3D genome organization in healthy and COPD human lung tissue

> **Status: Under active development**
>
> This repository is under active development. Scripts and documentation will be updated as the Hi-C analysis progresses and additional steps are validated.

## Overview

## Overview

This project investigates three-dimensional genome organization in healthy and chronic obstructive pulmonary disease (COPD) human lung tissue using Hi-C sequencing data.

Hi-C data from three healthy control and three COPD lung samples are analyzed alongside publicly available human lung Hi-C datasets from ENCODE. Sequencing reads are processed using Juicer to generate multiresolution chromatin contact matrices, providing the basis for downstream analysis of chromatin architecture.

## 1. Human lung Hi-C datasets

The study includes six publicly available human lung Hi-C experiments from [ENCODE](https://www.encodeproject.org/), representing four individual donors, and six in-house Hi-C samples (three healthy controls and three COPD lungs).

Raw sequencing data for the ENCODE experiments are accessible through the linked GEO records and associated SRA accessions.

| ENCODE accession | GEO accession | Donor | Lung tissue |
|---|---|---|---|
| [ENCSR305GRJ](https://www.encodeproject.org/experiments/ENCSR305GRJ/) | [GSE237368](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237368) | W76, male, 40 | Left lung |
| [ENCSR749PID](https://www.encodeproject.org/experiments/ENCSR749PID/) | [GSE237921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237921) | W78, male, 60 | Right lower lobe |
| [ENCSR673UZF](https://www.encodeproject.org/experiments/ENCSR673UZF/) | [GSE237897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237897) | W78, male, 60 | Right upper lobe |
| [ENCSR446NJD](https://www.encodeproject.org/experiments/ENCSR446NJD/) | [GSE237470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237470) | W73, female, 59 | Left lower lobe |
| [ENCSR068FXF](https://www.encodeproject.org/experiments/ENCSR068FXF/) | [GSE237333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237333) | W78, male, 60 | Left lower lobe |
| [ENCSR892BZA](https://www.encodeproject.org/experiments/ENCSR892BZA/) | [GSE237904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237904) | W62, female, 16 | Left lung |

**In-house Hi-C data (3 healthy, 3 COPD):** GEO accession to be added upon release.

## 2. Hi-C processing using Juicer

Hi-C sequencing data were processed using [Juicer v1.6](https://github.com/aidenlab/juicer) with BWA alignment to GRCh38 (hg38) and restriction-site information corresponding to the Arima Hi-C protocol.

### Software versions

| Software | Version |
|---|---|
| Juicer | 1.6 (installation directory; exact release to be confirmed) |
| BWA | 0.7.17-r1188 |
| Java (OpenJDK) | 25.0.3 |
| samtools | 1.17 |
| SLURM | 25.05.6 |

### Juicer installation
### Install Juicer

Juicer v1.6 was used for Hi-C data processing. The software is available from the [Aiden Lab GitHub repository](https://github.com/aidenlab/juicer).

```bash
git clone https://github.com/aidenlab/juicer.git
```

The Juicer installation path and required reference files are defined in `config/juicer_config.sh`.


## 3. Juicer configuration

The Juicer installation paths, reference files, input directories, and processing parameters are specified in `config/juicer_config.sh`.

The GRCh38 reference genome (`hg38.fa`) was indexed using BWA:

```bash
bwa index /path/to/reference/hg38.fa
```

The Arima restriction-site file (`hg38_Arima.txt`) can be generated using the [Juicer tools](https://github.com/aidenlab/juicer). The chromosome sizes (`hg38.chrom.sizes`) and restriction-site file must correspond to the same reference genome assembly.

## 4. Input FASTQ files

The planned sample manifest, `config/samples.tsv`, contains the sample identifier and full paths to paired-end FASTQ files.

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `r1_fastq` | Full path to R1 FASTQ |
| `r2_fastq` | Full path to R2 FASTQ |

Example:

```text
sample_id	r1_fastq	r2_fastq
sample_A	/data/HiC/sample_A_R1.fastq.gz	/data/HiC/sample_A_R2.fastq.gz
sample_B	/data/HiC/sample_B_R1.fastq.gz	/data/HiC/sample_B_R2.fastq.gz
```

The file uses tab-separated columns. Each row represents one paired-end FASTQ pair.

Manifest-based input preparation is under development. The current Juicer script uses sample directories containing a `fastq/` subdirectory.

## 5. Run Juicer

Submit each sample using `scripts/juicer/run_juicer.slurm`:

```bash
sbatch --export=ALL,SAMPLE_ID=sample_A \
    scripts/juicer/run_juicer.slurm
```

Juicer generates `merged_nodups.txt`, `inter.hic`, and `inter_30.hic` in the sample's `aligned/` directory.


## 6. Merging Hi-C contact matrices

Hi-C contact matrices were merged at 5 kb resolution using [hictk v2.2.0](https://github.com/paulsengroup/hictk) (Rossini and Paulsen, [*Bioinformatics*, 2024](https://doi.org/10.1093/bioinformatics/btae408)).

The hictk Apptainer (`.sif`) container used in this analysis is available on [Zenodo](https://doi.org/10.5281/zenodo.23401176).

Users who already have `.hic` files can skip the Juicer processing steps and start directly with merging.

### Merge configuration

The merging workflow uses three files:

- `config/merge_hic_config.example.sh` — configuration template containing the hictk container path, output directories, merge resolution, and number of threads.
- `config/merge_samples.tsv` — list of input `.hic` files.
- `scripts/encode/merge_hic.slurm` — SLURM script for merging the contact matrices.

All input `.hic` files must use the same reference genome assembly, compatible chromosome names, sizes, and order, and contain the selected merge resolution.

### Run Hi-C merging

Submit the merging job using:

```bash
sbatch scripts/encode/merge_hic.slurm
```

The script merges the input contact matrices at 5 kb resolution and validates the resulting `.hic` file using hictk.

