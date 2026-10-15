
#!/usr/bin/env bash
#SBATCH --job-name=juicer
#SBATCH --time=6-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=400G
#SBATCH --partition=standard
#SBATCH --output=juicer_%j.out
#SBATCH --error=juicer_%j.err

set -euo pipefail

# Load software and configuration
module load bwa

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/juicer_config.sh"

# Select sample from samples.tsv
SAMPLE_ID="${SAMPLE_ID:?Set SAMPLE_ID when submitting the job}"
SAMPLES_FILE="${REPO_DIR}/config/samples.tsv"

SAMPLE_DIR="${PROJECT_DIR}/${SAMPLE_ID}"
FASTQ_DIR="${SAMPLE_DIR}/fastq"

mkdir -p "$FASTQ_DIR"

# Link paired FASTQs into Juicer's input directory
pair=0

while IFS=$'\t' read -r sample r1 r2; do
    [[ "$sample" == "sample_id" ]] && continue
    [[ "$sample" == "$SAMPLE_ID" ]] || continue

    [[ -f "$r1" && -f "$r2" ]] || {
        echo "FASTQ files missing for ${SAMPLE_ID}" >&2
        exit 1
    }

    pair=$((pair + 1))
    prefix=$(printf "%s_pair%03d" "$SAMPLE_ID" "$pair")

    ln -sfn "$r1" "${FASTQ_DIR}/${prefix}_R1.fastq.gz"
    ln -sfn "$r2" "${FASTQ_DIR}/${prefix}_R2.fastq.gz"

done < "$SAMPLES_FILE"

if (( pair == 0 )); then
    echo "Sample ${SAMPLE_ID} not found in samples.tsv" >&2
    exit 1
fi

# Run Juicer 1.6
bash "$JUICER_EXECUTABLE" \
    -g "$GENOME_ID" \
    -d "$SAMPLE_DIR" \
    -s "$RESTRICTION_ENZYME" \
    -p "$CHROM_SIZES" \
    -y "$RESTRICTION_SITES" \
    -z "$REFERENCE_GENOME" \
    -D "$JUICER_DIR" \
    -t "${SLURM_CPUS_PER_TASK:-24}" \
    -q "$QUEUE" \
    -l "$LONG_QUEUE"

