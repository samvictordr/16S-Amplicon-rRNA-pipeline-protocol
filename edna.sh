#!/bin/bash
# ==============================================================================
# QIIME 2 eDNA ANALYSIS PIPELINE — Entry Point
# ==============================================================================
# Runs all 6 pipeline stages in sequence. All arguments are forwarded to every
# stage, allowing config files and parameter overrides to apply globally.
#
# Usage: ./edna.sh [options]
#
# Options:
#   --config FILE         Path to config file (default: pipeline.conf)
#   --sample-id ID        EBI accession ID (default: ERR3444605)
#   --sample-name NAME    Sample label in the QIIME 2 manifest
#   --qiime-env ENV       Conda environment name (default: qiime2-amplicon-2025.7)
#   --trunc-f N           DADA2 forward read truncation length in bp (default: 240)
#   --trunc-r N           DADA2 reverse read truncation length in bp (default: 240)
#   --threads N           Number of threads for DADA2 (default: 1, 0 = all cores)
#   --benchmark MODE      Benchmarking mode: off, instrument, slurm (default: off)
#   --benchmark-dir DIR   Directory for benchmark CSV outputs (default: benchmarks)
#   --classifier-url URL  Download URL for the SILVA classifier artifact
#   --classifier-file F   Local filename for the classifier artifact
#   --output-dir DIR      Directory for all outputs (default: current directory)
#   --url-r1 URL          Manual override for forward read URL
#   --url-r2 URL          Manual override for reverse read URL
#   --start-at N          Resume pipeline from stage N (1-6)
#   --help                Show this message
#
# Examples:
#   ./edna.sh
#   ./edna.sh --sample-id ERR9876543 --sample-name my-sample --trunc-f 220
#   ./edna.sh --config my_sample.conf
#   ./edna.sh --start-at 3 --trunc-f 200 --trunc-r 200   # resume from DADA2
#   ./edna.sh --threads 4                                  # use 4 CPU cores
#   bash stages/03_denoise.sh --trunc-f 200 --trunc-r 200
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

# Handle --help anywhere in args
show_help "$0" "$@"

# Load config and CLI overrides
[ -f "$SCRIPT_DIR/pipeline.conf" ] && load_config "$SCRIPT_DIR/pipeline.conf"
parse_args "$@"
validate_params

# Pre-flight dependency check — install missing tools automatically
install_dependencies wget python

echo "========================================"
echo "   eDNA Metabarcoding Pipeline"
echo "========================================"
echo "  Sample ID:    $SAMPLE_ID"
echo "  Sample name:  $SAMPLE_NAME"
echo "  Output dir:   $OUTPUT_DIR"
echo "  Trunc F/R:    ${TRUNC_LEN_F}bp / ${TRUNC_LEN_R}bp"
echo "  Threads:      $THREADS"
[[ "${BENCHMARK_MODE:-off}" != "off" ]] && echo "  Benchmark:    $BENCHMARK_MODE → $BENCHMARK_DIR"
[[ "$START_AT" -gt 1 ]] && echo "  Resuming from: stage $START_AT"
echo "========================================"

STAGES=(
    "$SCRIPT_DIR/stages/01_download.sh"
    "$SCRIPT_DIR/stages/02_import.sh"
    "$SCRIPT_DIR/stages/03_denoise.sh"
    "$SCRIPT_DIR/stages/04_classify.sh"
    "$SCRIPT_DIR/stages/05_export.sh"
    "$SCRIPT_DIR/stages/06_visualize.sh"
)

for i in "${!STAGES[@]}"; do
    stage_num=$((i + 1))
    if [[ "$stage_num" -lt "$START_AT" ]]; then
        echo "Skipping stage $stage_num (--start-at $START_AT)"
        continue
    fi
    bash "${STAGES[$i]}" "$@"
done
