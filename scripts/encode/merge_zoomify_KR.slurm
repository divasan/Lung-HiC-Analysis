
#!/usr/bin/env bash
#SBATCH --job-name=hic_processing
#SBATCH --partition=standard
#SBATCH --time=6-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=350G
#SBATCH --output=hic_processing_%j.out
#SBATCH --error=hic_processing_%j.err

set -euo pipefail

module load apptainer
module load java

# Load the paths and processing settings
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/merge_hic_config.sh"

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

TMP="${TMP_DIR}/hic_processing_${SLURM_JOB_ID}"
mkdir -p "$TMP"

export TMPDIR="$TMP"
export APPTAINER_TMPDIR="$TMP"
export APPTAINER_CACHEDIR="${TMP_DIR}/apptainer_cache"
mkdir -p "$APPTAINER_CACHEDIR"

MERGED="${OUTPUT_DIR}/${MERGED_FILE}"
ZOOMIFIED="${OUTPUT_DIR}/${ZOOMIFIED_FILE}"
FINAL="${OUTPUT_DIR}/${FINAL_FILE}"

# Check that the required files are available
for f in "$HICTK_IMAGE" "$JUICER_TOOLS" "$MERGE_SAMPLES"; do
    [[ -f "$f" ]] || { echo "Missing file: $f" >&2; exit 1; }
done

# Read the .hic files listed in the sample manifest
HIC_FILES=()

while IFS=$'\t' read -r sample_id hic_file; do
    [[ -z "$sample_id" || "$sample_id" == "sample_id" ]] && continue

    [[ -f "$hic_file" ]] || {
        echo "Missing .hic file: $hic_file" >&2
        exit 1
    }

    HIC_FILES+=("$hic_file")
done < "$MERGE_SAMPLES"

if (( ${#HIC_FILES[@]} < 2 )); then
    echo "At least two .hic files are required." >&2
    exit 1
fi

# Do not overwrite previous results
for f in "$MERGED" "$ZOOMIFIED" "$FINAL"; do
    [[ ! -e "$f" ]] || { echo "Output already exists: $f" >&2; exit 1; }
done

echo "Hi-C processing started: $(date)"
echo "Number of input files: ${#HIC_FILES[@]}"
echo "Resolutions: ${ZOOMIFY_RESOLUTIONS[*]}"

# Step 1: Merge the input contact matrices at 5 kb
echo
echo "STEP 1: Merging Hi-C files"

apptainer exec "$HICTK_IMAGE" hictk merge \
    --resolution "$MERGE_RESOLUTION" \
    --threads "$THREADS" \
    --chunk-size 1000000 \
    --tmpdir "$TMP" \
    --output-file "$MERGED" \
    "${HIC_FILES[@]}"

apptainer exec "$HICTK_IMAGE" hictk validate "$MERGED"

echo "Merge completed: $(date)"

# Step 2: Generate additional resolutions from the merged matrix
echo
echo "STEP 2: Generating multiresolution Hi-C matrix"

apptainer exec "$HICTK_IMAGE" hictk zoomify \
    "$MERGED" "$ZOOMIFIED" \
    --resolutions "${ZOOMIFY_RESOLUTIONS[@]}" \
    --copy-base-resolution \
    --threads "$THREADS" \
    --chunk-size 1000000 \
    --tmpdir "$TMP"

# Confirm that all requested resolutions are present
apptainer exec "$HICTK_IMAGE" hictk dump \
    "$ZOOMIFIED" --table resolutions > "$TMP/resolutions.txt"

echo "Resolutions in zoomified file:"
cat "$TMP/resolutions.txt"

for res in "${ZOOMIFY_RESOLUTIONS[@]}"; do
    grep -qx "$res" "$TMP/resolutions.txt" || {
        echo "Missing resolution: $res bp" >&2
        exit 1
    }
done

apptainer exec "$HICTK_IMAGE" hictk validate "$ZOOMIFIED"
java -Xmx16g -jar "$JUICER_TOOLS" validate "$ZOOMIFIED"

echo "Zoomify completed: $(date)"

# Step 3: Calculate KR normalization at all available resolutions
echo
echo "STEP 3: KR normalization"

cp "$ZOOMIFIED" "$FINAL"

java -Xmx280g -jar "$JUICER_TOOLS" addNorm \
    -k KR \
    -j "$THREADS" \
    -F \
    "$FINAL"

echo "KR normalization completed: $(date)"

# Check the final file before using it for downstream analysis
echo
echo "STEP 4: Final validation"

apptainer exec "$HICTK_IMAGE" hictk validate "$FINAL"
java -Xmx16g -jar "$JUICER_TOOLS" validate "$FINAL"

# Make sure none of the requested resolutions were lost
apptainer exec "$HICTK_IMAGE" hictk dump \
    "$FINAL" --table resolutions > "$TMP/final_resolutions.txt"

for res in "${ZOOMIFY_RESOLUTIONS[@]}"; do
    grep -qx "$res" "$TMP/final_resolutions.txt" || {
        echo "Missing final resolution: $res bp" >&2
        exit 1
    }
done

# Confirm that a KR vector exists at 5 kb for chromosome 1
java -Xmx16g -jar "$JUICER_TOOLS" dump norm KR \
    "$FINAL" chr1 BP 5000 "$TMP/chr1_KR_5kb.txt"

[[ -s "$TMP/chr1_KR_5kb.txt" ]] || {
    echo "5 kb KR vector missing for chr1" >&2
    exit 1
}

echo
echo "Hi-C processing completed successfully: $(date)"
echo "Final file: $FINAL"

ls -lh "$FINAL"
