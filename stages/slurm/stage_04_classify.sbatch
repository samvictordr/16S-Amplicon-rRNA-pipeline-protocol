#!/bin/bash
#SBATCH --job-name=classify_scaling
#SBATCH --partition=main
#SBATCH --output=benchmarks/slurm_%j_classify.out
#SBATCH --error=benchmarks/slurm_%j_classify.err
#
# stages/slurm/stage_04_classify.sbatch — Taxonomic classification under SLURM
# Submit with: sbatch --cpus-per-task=N --mem=4G stages/slurm/stage_04_classify.sbatch
#
# Environment variables:
#   BM_THREADS       — number of sklearn jobs (defaults to SLURM_CPUS_PER_TASK)
#   BM_READS         — read count for benchmark logging
#   BM_CSV           — path to benchmark CSV output
#   BM_CLASSIFIER    — path to classifier .qza file
#   BM_OUTPUT_DIR    — directory containing rep-seqs.qza and for outputs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
source "$SCRIPT_DIR/lib/benchmark.sh"

# Resolve parameters
THREADS="${BM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
READS="${BM_READS:-55000}"
CSV_FILE="${BM_CSV:-benchmarks/strong_scaling.csv}"
CLASSIFIER="${BM_CLASSIFIER:-silva-138-99-515-806-nb-classifier.qza}"
WORK_DIR="${BM_OUTPUT_DIR:-.}"

echo "=== Taxonomic Classification (SLURM Job ${SLURM_JOB_ID:-local}) ==="
echo "  Threads: $THREADS"
echo "  Reads:   $READS"
echo "  Classifier: $CLASSIFIER"
echo "  Node:    $(hostname)"
echo "  Started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Activate QIIME2 environment
activate_qiime

cd "$WORK_DIR"

# Run classification with benchmarking
benchmark_stage "classify" "$THREADS" "$READS" "$CSV_FILE" \
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" \
        --i-reads rep-seqs.qza \
        --p-n-jobs "$THREADS" \
        --o-classification taxonomy.qza \
        --verbose

echo "Finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
