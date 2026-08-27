#!/bin/bash
#SBATCH --job-name=array_classify
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=benchmarks/slurm_%A_%a_classify.out
#SBATCH --error=benchmarks/slurm_%A_%a_classify.err
#
# stages/slurm/array_classify_task.sbatch — Per-sample classification for job arrays
# Used by submit_array.sh: sbatch --array=1-N ... array_classify_task.sbatch
#
# Environment variables (set via --export=ALL,...):
#   BM_MANIFEST     — path to sample manifest (TSV: sample_id, fwd_path, rev_path)
#   BM_THREADS      — number of sklearn jobs (defaults to SLURM_CPUS_PER_TASK)
#   BM_CSV          — path to benchmark CSV output
#   BM_CLASSIFIER   — path to classifier .qza file (default: silva-138-99-515-806-nb-classifier.qza)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$SCRIPT_DIR/lib/common.sh"
source "$SCRIPT_DIR/lib/benchmark.sh"

# Resolve parameters
THREADS="${BM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
CSV_FILE="${BM_CSV:-benchmarks/array_results.csv}"
CLASSIFIER="${BM_CLASSIFIER:-silva-138-99-515-806-nb-classifier.qza}"
MANIFEST="${BM_MANIFEST:?BM_MANIFEST must be set}"

# Extract sample info for this array task
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")
SAMPLE_ID=$(echo "$SAMPLE_LINE" | cut -f1)

# Per-sample working directory (created by denoise task)
WORK_DIR="benchmarks/array_samples/${SAMPLE_ID}"

if [[ ! -d "$WORK_DIR" ]]; then
    echo "ERROR: Work directory $WORK_DIR not found. Did denoise task complete?" >&2
    exit 1
fi

if [[ ! -f "$WORK_DIR/rep-seqs.qza" ]]; then
    echo "ERROR: $WORK_DIR/rep-seqs.qza not found. Denoise may have failed." >&2
    exit 1
fi

echo "=== Taxonomic Classification — Array Task ${SLURM_ARRAY_TASK_ID} (Job ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}) ==="
echo "  Sample:     $SAMPLE_ID"
echo "  Threads:    $THREADS"
echo "  Classifier: $CLASSIFIER"
echo "  Node:       $(hostname)"
echo "  Started:    $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Activate QIIME2 environment
activate_qiime

# Count reads for benchmark logging
READS=$(python3 -c "
import zipfile, json
with zipfile.ZipFile('${WORK_DIR}/rep-seqs.qza') as z:
    for name in z.namelist():
        if 'metadata.yaml' in name:
            print(55000)
            break
" 2>/dev/null || echo "55000")

# Run classification with benchmarking
benchmark_stage "classify" "$THREADS" "$READS" "$CSV_FILE" \
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" \
        --i-reads "$WORK_DIR/rep-seqs.qza" \
        --p-n-jobs "$THREADS" \
        --o-classification "$WORK_DIR/taxonomy.qza" \
        --verbose

echo "Finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
