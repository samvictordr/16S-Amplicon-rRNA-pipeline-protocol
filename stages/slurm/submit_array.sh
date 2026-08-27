#!/bin/bash
# stages/slurm/submit_array.sh — Multi-sample SLURM job array
# Processes N samples in parallel using SLURM job arrays with dependency chains.
# Each array task handles one sample through the full pipeline.
#
# Usage: bash stages/slurm/submit_array.sh <sample_manifest.txt> [--dry-run]
#
# sample_manifest.txt format (one sample dir per line):
#   data/subsamples/sample_01
#   data/subsamples/sample_02
#   ...

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

MANIFEST="${1:-}"
DRY_RUN=false
THREADS_PER_SAMPLE=4
CSV_FILE="benchmarks/array_results.csv"
MEM_PER_JOB="4G"
TIME_LIMIT="02:00:00"

if [[ -z "$MANIFEST" || ! -f "$MANIFEST" ]]; then
    echo "Usage: $0 <sample_manifest.txt> [--dry-run]"
    echo "  Each line in manifest: path to a sample directory containing demux.qza"
    exit 1
fi

for arg in "${@:2}"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
    esac
done

NUM_SAMPLES=$(wc -l < "$MANIFEST")
mkdir -p benchmarks

echo "=== Multi-Sample SLURM Job Array ==="
echo "  Manifest: $MANIFEST"
echo "  Samples: $NUM_SAMPLES"
echo "  Threads per sample: $THREADS_PER_SAMPLE"
echo "  Output: $CSV_FILE"
echo ""

# Submit denoise array
DENOISE_CMD="sbatch \
    --array=1-${NUM_SAMPLES} \
    --cpus-per-task=$THREADS_PER_SAMPLE \
    --mem=$MEM_PER_JOB \
    --time=$TIME_LIMIT \
    --job-name=array_denoise \
    --export=ALL,BM_MANIFEST=$MANIFEST,BM_CSV=$CSV_FILE,BM_THREADS=$THREADS_PER_SAMPLE \
    stages/slurm/array_denoise_task.sbatch"

if $DRY_RUN; then
    echo "[DRY RUN] $DENOISE_CMD"
    echo "[DRY RUN] Would submit classify array with --dependency"
else
    DENOISE_ARRAY_JOB=$($DENOISE_CMD | awk '{print $NF}')
    echo "Denoise array job: $DENOISE_ARRAY_JOB (${NUM_SAMPLES} tasks)"

    # Submit classify array with dependency on entire denoise array
    CLASSIFY_ARRAY_JOB=$(sbatch \
        --dependency=aftercorr:$DENOISE_ARRAY_JOB \
        --array=1-${NUM_SAMPLES} \
        --cpus-per-task=$THREADS_PER_SAMPLE \
        --mem=$MEM_PER_JOB \
        --time=$TIME_LIMIT \
        --job-name=array_classify \
        --export=ALL,BM_MANIFEST=$MANIFEST,BM_CSV=$CSV_FILE,BM_THREADS=$THREADS_PER_SAMPLE \
        stages/slurm/array_classify_task.sbatch | awk '{print $NF}')
    echo "Classify array job: $CLASSIFY_ARRAY_JOB (depends on $DENOISE_ARRAY_JOB)"
fi

echo ""
echo "Monitor: squeue -u $(whoami)"
echo "Sequential baseline time: ~$((NUM_SAMPLES * 280))s estimated"
echo "Array parallel time: ~280s estimated (limited by longest sample)"
