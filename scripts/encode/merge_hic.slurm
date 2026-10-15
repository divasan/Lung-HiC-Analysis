
#!/usr/bin/env bash
#SBATCH --job-name=merge_hic
#SBATCH --time=6-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=350G
#SBATCH --partition=standard
#SBATCH --output=merge_hic_%j.out
#SBATCH --error=merge_hic_%j.err

set -euo pipefail

module load apptainer

# Load configuration
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/merge_hic_config.sh"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

export TMPDIR="$TMP_DIR"
export APPTAINER_TMPDIR="$TMP_DIR"
export APPTAINER_CACHEDIR="${TMP_DIR}/apptainer_cache"
mkdir -p "$APPTAINER_CACHEDIR"

# Read input .hic files
HIC_FILES=()

while IFS=$'\t' read -r sample_id hic_file; do
    [[ "$sample_id" == "sample_id" ]] && continue
    [[ -z "$sample_id" ]] && continue

    [[ -f "$hic_file" ]] || {
        echo "Missing .hic file: $hic_file" >&2
        exit 1
    }

    HIC_FILES+=("$hic_file")

done < "$MERGE_SAMPLES"

# Check input files
if (( ${#HIC_FILES[@]} < 2 )); then
    echo "At least two .hic files are required." >&2
    exit 1
fi

MERGED="${OUTPUT_DIR}/${MERGED_FILE}"

[[ ! -e "$MERGED" ]] || {
    echo "Output already exists: $MERGED" >&2
    exit 1
}

echo "Starting merge: $(date)"
echo "Number of files: ${#HIC_FILES[@]}"
echo "Resolution: ${MERGE_RESOLUTION} bp"

# Merge .hic files at 5 kb
apptainer exec "$HICTK_IMAGE" hictk merge \
    --resolution "$MERGE_RESOLUTION" \
    --threads "$THREADS" \
    --chunk-size 1000000 \
    --tmpdir "$TMP_DIR" \
    --output-file "$MERGED" \
    "${HIC_FILES[@]}"

echo "Merge completed: $(date)"

# Validate merged .hic file
apptainer exec "$HICTK_IMAGE" hictk validate "$MERGED"

echo "Finished: $(date)"
ls -lh "$MERGED"

