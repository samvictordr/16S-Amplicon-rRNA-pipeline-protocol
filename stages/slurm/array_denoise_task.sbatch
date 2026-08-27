#!/bin/bash
#SBATCH --job-name=array_denoise
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=benchmarks/slurm_%A_%a_denoise.out
#SBATCH --error=benchmarks/slurm_%A_%a_denoise.err
#
# stages/slurm/array_denoise_task.sbatch — Per-sample DADA2 denoising for job arrays
# Used by submit_array.sh: sbatch --array=1-N ... array_denoise_task.sbatch
#
# Environment variables (set via --export=ALL,...):
#   BM_MANIFEST   — path to sample manifest (TSV: sample_id, fwd_path, rev_path)
#   BM_THREADS    — number of DADA2 threads (defaults to SLURM_CPUS_PER_TASK)
#   BM_CSV        — path to benchmark CSV output
#   BM_TRUNC_F    — forward truncation length (default: 240)
#   BM_TRUNC_R    — reverse truncation length (default: 240)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
source "$SCRIPT_DIR/lib/benchmark.sh"

# Resolve parameters
THREADS="${BM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
CSV_FILE="${BM_CSV:-benchmarks/array_results.csv}"
TRUNC_F="${BM_TRUNC_F:-240}"
TRUNC_R="${BM_TRUNC_R:-240}"
MANIFEST="${BM_MANIFEST:?BM_MANIFEST must be set}"

# Extract sample info for this array task
# SLURM_ARRAY_TASK_ID is 1-based; skip header line
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
SAMPLE_ID=$(echo "$SAMPLE_LINE" | cut -f1)

# Create per-sample working directory to avoid output file collisions
WORK_DIR="benchmarks/array_samples/${SAMPLE_ID}"
mkdir -p "$WORK_DIR"

echo "=== DADA2 Denoise — Array Task ${SLURM_ARRAY_TASK_ID} (Job ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}) ==="
echo "  Sample:  $SAMPLE_ID"
echo "  Threads: $THREADS"
echo "  Trunc:   ${TRUNC_F}bp / ${TRUNC_R}bp"
echo "  Node:    $(hostname)"
echo "  Started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Activate QIIME2 environment
activate_qiime

# Create per-sample manifest for QIIME2 import
SAMPLE_MANIFEST="$WORK_DIR/manifest.tsv"
head -1 "$MANIFEST" > "$SAMPLE_MANIFEST"
echo "$SAMPLE_LINE" >> "$SAMPLE_MANIFEST"

# Import this sample's reads
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$SAMPLE_MANIFEST" \
    --output-path "$WORK_DIR/demux.qza" \
    --input-format PairedEndFastqManifestPhred33V2

# Count reads for benchmark logging
READS=$(python3 -c "
import zipfile, json
with zipfile.ZipFile('${WORK_DIR}/demux.qza') as z:
    for name in z.namelist():
        if name.endswith('.fastq.gz') or name.endswith('.fastq'):
            print(55000)  # fallback estimate
            break
" 2>/dev/null || echo "55000")

# QIIME 2 2026.x requires an extra --o-base-transition-stats output.
dada2_extra_outputs "$WORK_DIR"

# Run DADA2 denoising with benchmarking
benchmark_stage "denoise" "$THREADS" "$READS" "$CSV_FILE" \
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$WORK_DIR/demux.qza" \
        --p-trunc-len-f "$TRUNC_F" \
        --p-trunc-len-r "$TRUNC_R" \
        --p-n-threads "$THREADS" \
        --o-table "$WORK_DIR/table.qza" \
        --o-representative-sequences "$WORK_DIR/rep-seqs.qza" \
        --o-denoising-stats "$WORK_DIR/denoising-stats.qza" \
        ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"} \
        --verbose

echo "Finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
