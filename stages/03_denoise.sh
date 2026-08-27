#!/bin/bash
# stages/03_denoise.sh — Stage 3: Quality control and ASV inference via DADA2
#
# Usage: bash stages/03_denoise.sh [options]
#   --trunc-f N      Forward read truncation length in bp (default: 240)
#   --trunc-r N      Reverse read truncation length in bp (default: 240)
#   --threads N      Number of threads for DADA2 (default: 1, 0 = all cores)
#   --output-dir DIR Directory containing demux.qza and for output artifacts
#   --qiime-env ENV  Conda environment name for QIIME 2
#   --config FILE    Load configuration from file
#   --help           Show this message
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"

show_help "$0" "$@"

[ -f "$SCRIPT_DIR/../pipeline.conf" ] && load_config "$SCRIPT_DIR/../pipeline.conf"
parse_args "$@"

# Verify input
if [[ ! -f "${OUTPUT_DIR}/demux.qza" ]]; then
    echo "Error: ${OUTPUT_DIR}/demux.qza not found. Run stage 02_import.sh first." >&2
    exit 1
fi

echo "STEP 3: Running DADA2 for QC and ASV generation..."
echo "  Truncation: forward=${TRUNC_LEN_F}bp, reverse=${TRUNC_LEN_R}bp"
echo "  Threads:    $THREADS"

activate_qiime

# QIIME 2 2026.x requires an extra --o-base-transition-stats output here.
dada2_extra_outputs "$OUTPUT_DIR"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${OUTPUT_DIR}/demux.qza" \
  --p-trunc-len-f "$TRUNC_LEN_F" \
  --p-trunc-len-r "$TRUNC_LEN_R" \
  --p-n-threads "$THREADS" \
  --o-table "${OUTPUT_DIR}/table.qza" \
  --o-representative-sequences "${OUTPUT_DIR}/rep-seqs.qza" \
  --o-denoising-stats "${OUTPUT_DIR}/denoising-stats.qza" \
  ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"}

echo "DADA2 complete. Outputs: table.qza, rep-seqs.qza, denoising-stats.qza"
echo "-------------------------"
