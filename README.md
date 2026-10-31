
# Lung Hi-C Analysis

![Status](https://img.shields.io/badge/Status-Under%20Active%20Development-orange)

### 3D genome organization in healthy and COPD human lung tissue

This repository is under active development. Scripts and documentation will be updated as the Hi-C analysis progresses and additional steps are validated.

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

### Software versions used in this analysis

| Software | Version |
|---|---|
| Juicer | 1.6 |
| Juicer Tools | 2.20.00 |
| hictk | 2.2.0 |
| Apptainer | Version to be confirmed |
| BWA | 0.7.17-r1188 |
| Java (OpenJDK) | 25.0.3 |
| samtools | 1.17 |
| SLURM | 25.05.6 |
| NumPy | 1.26.4 |

### Juicer installation

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

## 6. Hi-C reproducibility analysis

To assess the reproducibility of chromatin contact patterns across healthy control, COPD, and ENCODE lung samples, Hi-C contact matrices from all 12 individual datasets were prepared for HiCRep analysis.

Reproducibility will be evaluated using both the original full-depth contact matrices and independently downsampled matrices to examine the influence of sequencing depth on the observed similarities between samples.

### Contact-depth calculation

Contact depths were calculated from the original `.hic` files using hictk v2.2.0 with the following parameters:

- **Genome assembly:** GRCh38 (hg38)
- **Resolution:** 100 kb
- **Chromosomes:** Autosomes 1–22
- **Contact type:** Cis interactions
- **Genomic distance:** 0–5 Mb, including diagonal contacts
- **Normalization:** NONE (raw contact counts)

The minimum eligible contact count across the 12 datasets was used as the common downsampling target.

HD124 had the lowest eligible contact depth, with **88,356,364 contacts**.

### Contact downsampling

To control for differences in contact depth, a second set of Hi-C matrices was generated by independently downsampling each sample to the minimum eligible contact depth.

For each sample, the sampling probability was calculated as:

\[
p_i = \frac{N_{\mathrm{target}}}{N_i}
\]

where \(N_i\) is the eligible contact count for sample \(i\), and \(N_{\mathrm{target}}\) is the common target depth.

Contact counts at 100-kb resolution were randomly thinned using binomial sampling implemented in NumPy with a fixed random seed (20261031). The same sampling probability was applied to all extracted contact pixels, including cis and trans interactions.

The resulting eligible contact depths are approximately matched because of stochastic sampling.

Each downsampled contact matrix was reconstructed as an individual 100-kb `.hic` file using hictk. Reconstructed files were validated, and retained contact counts were compared with the counts recovered from the output matrices.

The original `.hic` files were preserved for full-depth analysis. This approach downsamples binned contact counts rather than individual sequencing read pairs.

### Scripts and input files

The depth-matching workflow uses:

- `scripts/hicrep/downsample_hic.slurm` — SLURM array script for contact-depth measurement, aggregation, and downsampling.
- `scripts/hicrep/data/hicrep_samples.example.tsv` — example manifest containing sample IDs and input `.hic` file paths.

Users should create a tab-separated `scripts/hicrep/data/hicrep_samples.tsv` file containing their own sample IDs and input file paths.

### Run contact-depth matching

The script supports three modes: `measure`, `aggregate`, and `downsample`.

**Step 1. Calculate contact depths**

```bash
sbatch --export=ALL,MODE=measure scripts/hicrep/downsample_hic.slurm
```

**Step 2. Calculate the common target depth**

After the measurement jobs complete:

```bash
sbatch --array=0 --export=ALL,MODE=aggregate scripts/hicrep/downsample_hic.slurm
```

**Step 3. Generate downsampled contact matrices**

After aggregation completes:

```bash
sbatch --export=ALL,MODE=downsample scripts/hicrep/downsample_hic.slurm
```

The script produces separate downsampled `.hic` files and corresponding QC records for each sample. Original contact matrices are not modified.

### HiCRep comparison

HiCRep will be used to calculate stratum-adjusted correlation coefficients (SCCs) between individual lung Hi-C datasets.

Two sets of comparisons are planned:

1. **Full-depth HiCRep:** Pairwise reproducibility between the original individual Hi-C contact matrices.
2. **Depth-matched HiCRep:** The same pairwise comparisons using the downsampled contact matrices.

Both analyses will use 100-kb contact matrices, allowing direct comparison of SCCs before and after depth matching.

This approach will help assess whether differences in reproducibility among healthy control, COPD, and ENCODE datasets persist after controlling for contact depth.

## 7. Merging Hi-C contact matrices and KR normalization

Hi-C contact matrices from the in-house and ENCODE datasets were merged at 5 kb resolution using [hictk v2.2.0](https://github.com/paulsengroup/hictk) (Rossini and Paulsen, [*Bioinformatics*, 2024](https://doi.org/10.1093/bioinformatics/btae408)).

The merged contact matrix was then processed using hictk to generate contact matrices at 5 kb, 10 kb, 20 kb, 100 kb, 500 kb, and 1 Mb resolutions. KR normalization was performed using Juicer Tools v2.20.00.

The hictk Apptainer (`.sif`) container used in this analysis is available on [Zenodo](https://doi.org/10.5281/zenodo.23401176).

Users who already have `.hic` files can skip the Juicer processing steps and start directly with merging.

### Merge and normalization configuration

The workflow uses three files:

- `config/merge_hic_config.example.sh` — configuration template containing the software paths, output directories, resolutions, and number of threads.
- `config/merge_samples.tsv` — list of input `.hic` files.
- `scripts/encode/merge_zoomify_KR.slurm` — SLURM script for merging, generating additional resolutions, and KR normalization.

All input `.hic` files must use the same reference genome assembly, compatible chromosome names, sizes, and order, and contain the selected merge resolution.

### Run merging and KR normalization

Submit the job using:

```bash
sbatch scripts/encode/merge_zoomify_KR.slurm
```

The script merges the input contact matrices at 5 kb, generates the additional resolutions, and performs KR normalization. The output files are validated using hictk and Juicer Tools.

The final multiresolution `.hic` file is saved in the output directory specified in the configuration file.
