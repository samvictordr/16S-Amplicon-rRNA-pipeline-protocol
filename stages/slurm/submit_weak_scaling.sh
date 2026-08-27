#!/bin/bash
# stages/slurm/submit_weak_scaling.sh — Weak scaling experiment
# Submits jobs with proportionally increasing dataset size and core count.
# Tests whether time remains constant as both scale together.
#
# Usage: bash stages/slurm/submit_weak_scaling.sh [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

# Weak scaling pairs: (reads, cores) — ~7K reads per core
SCALING_PAIRS=(
    "7000:1"
    "14000:2"
    "28000:4"
    "55000:8"
)

CSV_FILE="benchmarks/weak_scaling.csv"
MEM_PER_JOB="4G"
TIME_LIMIT="02:00:00"
DRY_RUN=false
REPETITIONS=3

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
    esac
done

mkdir -p benchmarks

echo "=== Weak Scaling Experiment ==="
echo "  Pairs (reads:cores): ${SCALING_PAIRS[*]}"
echo "  Repetitions: $REPETITIONS"
echo "  Output: $CSV_FILE"
echo ""

for rep in $(seq 1 $REPETITIONS); do
    echo "--- Repetition $rep/$REPETITIONS ---"
    for pair in "${SCALING_PAIRS[@]}"; do
        reads="${pair%%:*}"
        cores="${pair##*:}"

        echo "Submitting denoise (reads=$reads, cores=$cores, rep=$rep)..."

        DENOISE_CMD="sbatch \
            --cpus-per-task=$cores \
            --mem=$MEM_PER_JOB \
            --time=$TIME_LIMIT \
            --job-name=weak_denoise_r${reads}_c${cores}_rep${rep} \
            --export=ALL,BM_THREADS=$cores,BM_READS=$reads,BM_CSV=$CSV_FILE,BM_OUTPUT_DIR=data/subsamples/${reads} \
            stages/slurm/stage_03_denoise.sbatch"

        if $DRY_RUN; then
            echo "  [DRY RUN] $DENOISE_CMD"
        else
            DENOISE_JOB=$($DENOISE_CMD | awk '{print $NF}')
            echo "  Denoise job: $DENOISE_JOB"

            CLASSIFY_JOB=$(sbatch \
                --dependency=afterok:$DENOISE_JOB \
                --cpus-per-task=$cores \
                --mem=$MEM_PER_JOB \
                --time=$TIME_LIMIT \
                --job-name=weak_classify_r${reads}_c${cores}_rep${rep} \
                --export=ALL,BM_THREADS=$cores,BM_READS=$reads,BM_CSV=$CSV_FILE,BM_OUTPUT_DIR=data/subsamples/${reads} \
                stages/slurm/stage_04_classify.sbatch | awk '{print $NF}')
            echo "  Classify job: $CLASSIFY_JOB (depends on $DENOISE_JOB)"
        fi
    done
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u $(whoami)"
echo "Results will be written to: $CSV_FILE"
