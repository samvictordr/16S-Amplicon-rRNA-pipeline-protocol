#!/bin/bash
# stages/05_export.sh — Stage 5: Export QIIME 2 artifacts to TSV for downstream use
#
# Usage: bash stages/05_export.sh [options]
#   --output-dir DIR  Directory containing table.qza and taxonomy.qza
#   --qiime-env ENV   Conda environment name for QIIME 2
#   --config FILE     Load configuration from file
#   --help            Show this message
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"

show_help "$0" "$@"

[ -f "$SCRIPT_DIR/../pipeline.conf" ] && load_config "$SCRIPT_DIR/../pipeline.conf"
parse_args "$@"

# Verify inputs
for artifact in table.qza taxonomy.qza; do
    if [[ ! -f "${OUTPUT_DIR}/${artifact}" ]]; then
        echo "Error: ${OUTPUT_DIR}/${artifact} not found. Run earlier stages first." >&2
        exit 1
    fi
done

echo "STEP 5: Exporting data for custom plotting..."

activate_qiime

EXPORT_DIR="${OUTPUT_DIR}/exported-data"
mkdir -p "$EXPORT_DIR"

# Export feature table (BIOM) and convert to TSV
qiime tools export \
  --input-path "${OUTPUT_DIR}/table.qza" \
  --output-path "$EXPORT_DIR"

biom convert \
  -i "${EXPORT_DIR}/feature-table.biom" \
  -o "${EXPORT_DIR}/feature-table.tsv" \
  --to-tsv

# Export taxonomy assignments
qiime tools export \
  --input-path "${OUTPUT_DIR}/taxonomy.qza" \
  --output-path "$EXPORT_DIR"

echo "Data export complete. Files in: $EXPORT_DIR"
echo "-------------------------"
