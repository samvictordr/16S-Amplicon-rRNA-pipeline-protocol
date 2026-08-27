#!/bin/bash
# stages/slurm/submit_strong_scaling.sh — Strong scaling experiment
# Submits DADA2 denoise and classifier jobs across thread counts 1,2,4,8,16
# with fixed dataset size (55K reads). Measures speedup and efficiency.
#
# Usage: bash stages/slurm/submit_strong_scaling.sh [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

THREAD_COUNTS=(1 2 4 8 16)
READS=55000
CSV_FILE="benchmarks/strong_scaling.csv"
MEM_PER_JOB="4G"
TIME_LIMIT="02:00:00"
DRY_RUN=false
REPETITIONS=3

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
    esac
done

mkdir -p benchmarks

echo "=== Strong Scaling Experiment ==="
echo "  Thread counts: ${THREAD_COUNTS[*]}"
echo "  Reads: $READS"
echo "  Repetitions: $REPETITIONS"
echo "  Output: $CSV_FILE"
echo ""

for rep in $(seq 1 $REPETITIONS); do
    echo "--- Repetition $rep/$REPETITIONS ---"
    for threads in "${THREAD_COUNTS[@]}"; do
        echo "Submitting denoise (threads=$threads, rep=$rep)..."

        DENOISE_CMD="sbatch \
            --cpus-per-task=$threads \
            --mem=$MEM_PER_JOB \
            --time=$TIME_LIMIT \
            --job-name=denoise_t${threads}_r${rep} \
            --export=ALL,BM_THREADS=$threads,BM_READS=$READS,BM_CSV=$CSV_FILE \
            stages/slurm/stage_03_denoise.sbatch"

        if $DRY_RUN; then
            echo "  [DRY RUN] $DENOISE_CMD"
        else
            DENOISE_JOB=$($DENOISE_CMD | awk '{print $NF}')
            echo "  Denoise job: $DENOISE_JOB"

            # Submit classifier with dependency on denoise completion
            echo "Submitting classify (threads=$threads, rep=$rep)..."
            CLASSIFY_JOB=$(sbatch \
                --dependency=afterok:$DENOISE_JOB \
                --cpus-per-task=$threads \
                --mem=$MEM_PER_JOB \
                --time=$TIME_LIMIT \
                --job-name=classify_t${threads}_r${rep} \
                --export=ALL,BM_THREADS=$threads,BM_READS=$READS,BM_CSV=$CSV_FILE \
                stages/slurm/stage_04_classify.sbatch | awk '{print $NF}')
            echo "  Classify job: $CLASSIFY_JOB (depends on $DENOISE_JOB)"
        fi
    done
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u $(whoami)"
echo "Results will be written to: $CSV_FILE"
