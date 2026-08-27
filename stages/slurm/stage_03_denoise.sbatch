#!/bin/bash
#SBATCH --job-name=denoise_scaling
#SBATCH --partition=main
#SBATCH --output=benchmarks/slurm_%j_denoise.out
#SBATCH --error=benchmarks/slurm_%j_denoise.err
#
# stages/slurm/stage_03_denoise.sbatch — DADA2 denoising under SLURM
# Submit with: sbatch --cpus-per-task=N --mem=4G stages/slurm/stage_03_denoise.sbatch
#
# Environment variables (set via submit scripts or --export):
#   BM_THREADS    — number of DADA2 threads (defaults to SLURM_CPUS_PER_TASK)
#   BM_READS      — read count for benchmark logging
#   BM_CSV        — path to benchmark CSV output
#   BM_TRUNC_F    — forward truncation length (default: 240)
#   BM_TRUNC_R    — reverse truncation length (default: 240)
#   BM_OUTPUT_DIR — directory containing demux.qza and for outputs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
source "$SCRIPT_DIR/lib/benchmark.sh"

# Resolve parameters
THREADS="${BM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
READS="${BM_READS:-55000}"
CSV_FILE="${BM_CSV:-benchmarks/strong_scaling.csv}"
TRUNC_F="${BM_TRUNC_F:-240}"
TRUNC_R="${BM_TRUNC_R:-240}"
WORK_DIR="${BM_OUTPUT_DIR:-.}"

echo "=== DADA2 Denoise (SLURM Job ${SLURM_JOB_ID:-local}) ==="
echo "  Threads: $THREADS"
echo "  Reads:   $READS"
echo "  Trunc:   ${TRUNC_F}bp / ${TRUNC_R}bp"
echo "  Node:    $(hostname)"
echo "  Started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Activate QIIME2 environment
activate_qiime

cd "$WORK_DIR"

# QIIME 2 2026.x requires an extra --o-base-transition-stats output.
dada2_extra_outputs "$WORK_DIR"

# Run DADA2 denoising with benchmarking
benchmark_stage "denoise" "$THREADS" "$READS" "$CSV_FILE" \
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux.qza \
        --p-trunc-len-f "$TRUNC_F" \
        --p-trunc-len-r "$TRUNC_R" \
        --p-n-threads "$THREADS" \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"} \
        --verbose

echo "Finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
