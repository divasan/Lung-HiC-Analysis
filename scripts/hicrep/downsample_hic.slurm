#!/usr/bin/env bash
#SBATCH --job-name=hic_downsample
#SBATCH --partition=standard
#SBATCH --array=0-11%4
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=hic_downsample_%A_%a.out
#SBATCH --error=hic_downsample_%A_%a.err

set -euo pipefail

# ============================================================
# Hi-C contact depth measurement and downsampling
#
# MODE=measure:
#   Calculate autosomal cis contact depths at 100 kb
#   for genomic distances of 0-5 Mb.
#
# MODE=aggregate:
#   Combine sample depths and calculate the minimum depth.
#
# MODE=downsample:
#   Randomly downsample individual .hic matrices to
#   the common contact depth using binomial sampling.
#
# Each array task processes one sample.
# Original .hic files are not modified.
# ============================================================

module load apptainer

# Directories and software

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

HICTK_IMAGE="${HICTK_IMAGE:-/path/to/hictk_2.2.0.sif}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

OUTPUT_DIR="${OUTPUT_DIR:-/path/to/downsampled_hic}"
WORK_DIR="${WORK_DIR:-/path/to/hicrep_work}"

SAMPLE_MANIFEST="${SAMPLE_MANIFEST:-${REPO_DIR}/scripts/hicrep/data/hicrep_samples.tsv}"

# Analysis settings

MODE="${MODE:-measure}"

RES=100000
MAX_DIST=5000000
SEED=20261031

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
JOB_ID="${SLURM_JOB_ID:-manual}"

TMP="${WORK_DIR}/tmp_${JOB_ID}_${TASK_ID}"

# Check required files

[[ -f "$HICTK_IMAGE" ]] || {
    echo "Missing hictk image: $HICTK_IMAGE" >&2
    exit 1
}

[[ -f "$SAMPLE_MANIFEST" ]] || {
    echo "Missing sample manifest: $SAMPLE_MANIFEST" >&2
    exit 1
}

command -v "$PYTHON_BIN" >/dev/null 2>&1 || {
    echo "Python not found: $PYTHON_BIN" >&2
    exit 1
}

mkdir -p "$OUTPUT_DIR" "$WORK_DIR" "$TMP"

export TMPDIR="$TMP"
export APPTAINER_TMPDIR="$TMP"
export APPTAINER_CACHEDIR="${WORK_DIR}/apptainer_cache"

mkdir -p "$APPTAINER_CACHEDIR"

echo "=============================================="
echo "Hi-C contact depth matching"
echo "=============================================="
echo "Mode: $MODE"
echo "Task ID: $TASK_ID"
echo "Resolution: $RES"
echo "Maximum distance: $MAX_DIST"
echo "Threads: $THREADS"
echo "Start: $(date)"
echo "=============================================="

# Run contact depth measurement or downsampling

"$PYTHON_BIN" - \
    "$MODE" \
    "$TASK_ID" \
    "$SAMPLE_MANIFEST" \
    "$OUTPUT_DIR" \
    "$WORK_DIR" \
    "$TMP" \
    "$HICTK_IMAGE" \
    "$RES" \
    "$MAX_DIST" \
    "$SEED" \
    "$THREADS" <<'PY'

import sys
import csv
import gzip
import json
import os
import subprocess
from pathlib import Path

import numpy as np


(
    mode,
    task_id,
    manifest,
    output_dir,
    work_dir,
    tmp_dir,
    hictk_image,
    resolution,
    max_distance,
    seed,
    threads,
) = sys.argv[1:]

task_id = int(task_id)
resolution = int(resolution)
max_distance = int(max_distance)
seed = int(seed)
threads = int(threads)

output_dir = Path(output_dir)
work_dir = Path(work_dir)
tmp_dir = Path(tmp_dir)

depth_dir = work_dir / "depths"
depth_table = work_dir / "hicrep_depths_100kb.tsv"

depth_dir.mkdir(parents=True, exist_ok=True)

AUTOSOMES = {str(i) for i in range(1, 23)}


# ------------------------------------------------------------
# Read the sample manifest
# ------------------------------------------------------------

samples = []

with open(manifest, newline="") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if not row or not row[0].strip():
            continue

        if row[0].strip().startswith("#"):
            continue

        if row[0].strip().lower() == "sample_id":
            continue

        if len(row) != 2:
            raise ValueError(
                "Manifest must have two columns: "
                "sample_id and hic_file"
            )

        sample_id = row[0].strip()
        hic_file = row[1].strip()

        if not sample_id or not hic_file:
            raise ValueError("Empty sample ID or file path")

        samples.append((sample_id, hic_file))

if not samples:
    raise RuntimeError("No samples found in manifest")

sample_names = [name for name, _ in samples]

if len(sample_names) != len(set(sample_names)):
    raise RuntimeError("Duplicate sample IDs in manifest")

if mode in ("measure", "downsample"):
    if not 0 <= task_id < len(samples):
        raise RuntimeError(
            f"Invalid task ID {task_id}; "
            f"manifest has {len(samples)} samples"
        )

    sample, input_file = samples[task_id]

    if not Path(input_file).is_file():
        raise FileNotFoundError(input_file)


# ------------------------------------------------------------
# hictk commands
# ------------------------------------------------------------

def hictk_cmd(*args):
    return [
        "apptainer",
        "exec",
        "--bind", "/standard,/scratch",
        hictk_image,
        "hictk",
        *map(str, args),
    ]


def stream_pixels(path):
    """
    Stream unnormalized contact pixels in BG2 format.
    """

    command = hictk_cmd(
        "dump",
        path,
        "--resolution", resolution,
        "--matrix-type", "observed",
        "--balance", "NONE",
        "--join",
    )

    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )

    try:
        for line in proc.stdout:
            fields = line.strip().split()

            if not fields:
                continue

            if len(fields) != 7:
                raise ValueError(
                    f"Unexpected contact record: {line[:200]}"
                )

            yield fields

    finally:
        proc.stdout.close()
        returncode = proc.wait()

        if returncode != 0:
            raise RuntimeError(
                f"hictk dump failed for {path}: {returncode}"
            )


def normalize_chrom(chrom):
    if chrom.lower().startswith("chr"):
        return chrom[3:]

    return chrom


def is_eligible(fields):
    """
    Autosomal cis contacts within 0-5 Mb.
    Includes diagonal contacts.
    """

    chrom1 = normalize_chrom(fields[0])
    chrom2 = normalize_chrom(fields[3])

    if chrom1 != chrom2:
        return False

    if chrom1 not in AUTOSOMES:
        return False

    start1 = int(fields[1])
    start2 = int(fields[4])

    return abs(start2 - start1) <= max_distance


def get_count(fields):
    value = float(fields[6])

    if not np.isfinite(value):
        raise ValueError("Nonfinite contact count")

    if value < 0 or not value.is_integer():
        raise ValueError(
            f"Invalid contact count: {value}"
        )

    return int(value)


def write_chrom_sizes(hic_path, destination):
    """
    Extract chromosome sizes from the input .hic.
    """

    result = subprocess.run(
        hictk_cmd(
            "dump",
            hic_path,
            "--table", "chroms",
        ),
        capture_output=True,
        text=True,
        check=True,
    )

    n_chroms = 0

    with open(destination, "w") as f:
        for line in result.stdout.splitlines():
            fields = line.split()

            if len(fields) != 2:
                continue

            chrom, length = fields

            if chrom.lower() in ("all", "chrall"):
                continue

            try:
                length = int(length)
            except ValueError:
                continue

            if length <= 0:
                continue

            f.write(f"{chrom}\t{length}\n")
            n_chroms += 1

    if n_chroms < 22:
        raise RuntimeError(
            f"Only {n_chroms} chromosomes found"
        )


def write_json(destination, data):
    temporary = destination.with_suffix(".tmp")

    with open(temporary, "w") as f:
        json.dump(data, f, indent=2)

    os.replace(temporary, destination)


# ============================================================
# Step 1: Measure contact depth
# ============================================================

if mode == "measure":

    total = 0
    pixels = 0

    print(f"\nMeasuring {sample}...", flush=True)

    for fields in stream_pixels(input_file):
        if is_eligible(fields):
            total += get_count(fields)
            pixels += 1

    if total <= 0:
        raise RuntimeError(
            f"No eligible contacts for {sample}"
        )

    result = {
        "sample": sample,
        "resolution": resolution,
        "max_distance": max_distance,
        "eligible_contacts": total,
        "eligible_pixels": pixels,
    }

    destination = depth_dir / f"{sample}.json"

    write_json(destination, result)

    print(f"Sample: {sample}")
    print(f"Eligible contacts: {total:,}")
    print(f"Eligible pixels: {pixels:,}")
    print(f"Depth file: {destination}")


# ============================================================
# Step 2: Calculate common contact depth
# ============================================================

elif mode == "aggregate":

    records = []

    for name in sample_names:
        path = depth_dir / f"{name}.json"

        if not path.exists():
            raise FileNotFoundError(
                f"Missing depth file: {path}"
            )

        with open(path) as f:
            record = json.load(f)

        if record["sample"] != name:
            raise RuntimeError(
                f"Sample mismatch: {path}"
            )

        if (
            record["resolution"] != resolution
            or record["max_distance"] != max_distance
        ):
            raise RuntimeError(
                f"Parameter mismatch: {path}"
            )

        if record["eligible_contacts"] <= 0:
            raise RuntimeError(
                f"Invalid contact depth: {name}"
            )

        records.append(record)

    target = min(
        r["eligible_contacts"] for r in records
    )

    temporary = depth_table.with_suffix(".tmp")

    with open(temporary, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        writer.writerow([
            "sample",
            "eligible_contacts",
            "eligible_pixels",
            "sampling_probability",
        ])

        for record in records:
            probability = (
                target / record["eligible_contacts"]
            )

            writer.writerow([
                record["sample"],
                record["eligible_contacts"],
                record["eligible_pixels"],
                f"{probability:.12f}",
            ])

    os.replace(temporary, depth_table)

    print("\nContact depths")

    for record in records:
        print(
            f"{record['sample']:18s} "
            f"{record['eligible_contacts']:,}"
        )

    print(f"\nCommon target: {target:,}")
    print(f"Depth table: {depth_table}")


# ============================================================
# Step 3: Downsample individual .hic files
# ============================================================

elif mode == "downsample":

    if not depth_table.exists():
        raise FileNotFoundError(
            "Depth table missing. Run aggregate first."
        )

    with open(depth_table) as f:
        records = list(
            csv.DictReader(f, delimiter="\t")
        )

    depths = {
        r["sample"]: int(r["eligible_contacts"])
        for r in records
    }

    if (
        len(records) != len(samples)
        or set(depths) != set(sample_names)
    ):
        raise RuntimeError(
            "Depth table does not match sample manifest"
        )

    target = min(depths.values())
    probability = target / depths[sample]

    if not 0 < probability <= 1:
        raise RuntimeError(
            f"Invalid sampling probability: {probability}"
        )

    output = (
        output_dir /
        f"{sample}_downsampled_100kb.hic"
    )

    if output.exists():
        raise FileExistsError(
            f"Output already exists: {output}"
        )

    rng = np.random.default_rng(seed + task_id)

    bg2 = tmp_dir / f"{sample}.bg2.gz"
    chrom_sizes = tmp_dir / f"{sample}.chrom.sizes"
    temporary_hic = tmp_dir / f"{sample}.hic"

    write_chrom_sizes(input_file, chrom_sizes)

    print(f"\nDownsampling {sample}", flush=True)
    print(f"Original eligible contacts: {depths[sample]:,}")
    print(f"Common target: {target:,}")
    print(f"Sampling probability: {probability:.6f}")
    print(f"Random seed: {seed + task_id}")

    original_total = 0
    retained_total = 0
    original_eligible = 0
    retained_eligible = 0
    pixels_kept = 0

    # Randomly retain contacts using binomial sampling

    with gzip.open(bg2, "wt", compresslevel=3) as fout:

        for fields in stream_pixels(input_file):

            count = get_count(fields)
            eligible = is_eligible(fields)

            original_total += count

            if eligible:
                original_eligible += count

            kept = int(
                rng.binomial(count, probability)
            )

            if kept == 0:
                continue

            retained_total += kept

            if eligible:
                retained_eligible += kept

            fields[6] = str(kept)

            fout.write("\t".join(fields) + "\n")

            pixels_kept += 1

    # Confirm that the original depth matches measurement

    if original_eligible != depths[sample]:
        raise RuntimeError(
            f"Depth mismatch for {sample}: "
            f"expected {depths[sample]}, "
            f"observed {original_eligible}"
        )

    if retained_eligible <= 0:
        raise RuntimeError(
            f"No contacts retained for {sample}"
        )

    print("\nSampling complete")
    print(f"Original total: {original_total:,}")
    print(f"Retained total: {retained_total:,}")
    print(f"Retained eligible: {retained_eligible:,}")
    print(f"Pixels retained: {pixels_kept:,}")

    # Reconstruct the downsampled .hic file

    print("\nReconstructing .hic...", flush=True)

    subprocess.run(
        hictk_cmd(
            "load",
            bg2,
            temporary_hic,
            "--chrom-sizes", chrom_sizes,
            "--format", "bg2",
            "--bin-size", resolution,
            "--assembly", "hg38",
            "--chunk-size", 1000000,
            "--threads", threads,
            "--tmpdir", tmp_dir,
        ),
        check=True,
    )

    # Validate the reconstructed file

    print("\nValidating .hic...", flush=True)

    subprocess.run(
        hictk_cmd(
            "validate",
            temporary_hic,
        ),
        check=True,
    )

    # Verify contact counts after reconstruction

    verified_eligible = 0

    for fields in stream_pixels(temporary_hic):
        if is_eligible(fields):
            verified_eligible += get_count(fields)

    if verified_eligible != retained_eligible:
        raise RuntimeError(
            f"Reconstructed contact count mismatch: "
            f"{retained_eligible} vs {verified_eligible}"
        )

    # Move the validated .hic file to the output directory

    os.replace(temporary_hic, output)

    # Save downsampling QC results

    qc = {
        "sample": sample,
        "resolution": resolution,
        "max_distance": max_distance,
        "seed": seed + task_id,
        "original_eligible": original_eligible,
        "target": target,
        "sampling_probability": probability,
        "retained_eligible": retained_eligible,
        "verified_eligible": verified_eligible,
        "original_total": original_total,
        "retained_total": retained_total,
        "pixels_retained": pixels_kept,
        "output_hic": str(output),
    }

    write_json(
        output_dir / f"{sample}_downsampling_qc.json",
        qc,
    )

    bg2.unlink()

    print("\nDownsampling completed successfully")
    print(f"Output: {output}")


else:
    raise ValueError(f"Unknown MODE: {mode}")

PY

echo
echo "=============================================="
echo "Job completed: $(date)"
echo "=============================================="

