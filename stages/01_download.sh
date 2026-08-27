#!/bin/bash
# stages/01_download.sh — Stage 1: Download raw paired-end FASTQ reads from EBI
#
# Usage: bash stages/01_download.sh [options]
#   --sample-id ID    EBI accession ID (auto-derives download URLs)
#   --url-r1 URL      Override auto-derived URL for forward reads
#   --url-r2 URL      Override auto-derived URL for reverse reads
#   --output-dir DIR  Directory to write downloaded files
#   --config FILE     Load configuration from file
#   --help            Show this message
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"

show_help "$0" "$@"

[ -f "$SCRIPT_DIR/../pipeline.conf" ] && load_config "$SCRIPT_DIR/../pipeline.conf"
parse_args "$@"

check_dependencies wget

mkdir -p "$OUTPUT_DIR"

# Derive EBI URLs (use manual overrides if provided)
URL_R1="${MANUAL_URL_R1:-$(derive_ebi_url "$SAMPLE_ID" 1)}"
URL_R2="${MANUAL_URL_R2:-$(derive_ebi_url "$SAMPLE_ID" 2)}"

echo "STEP 1: Downloading raw sequence data..."
echo "  Sample:  $SAMPLE_ID"
echo "  Read 1:  $URL_R1"
echo "  Read 2:  $URL_R2"

R1_GZ="${OUTPUT_DIR}/${SAMPLE_ID}_1.fastq.gz"
R2_GZ="${OUTPUT_DIR}/${SAMPLE_ID}_2.fastq.gz"
R1="${OUTPUT_DIR}/${SAMPLE_ID}_1.fastq"
R2="${OUTPUT_DIR}/${SAMPLE_ID}_2.fastq"

# Skip download if uncompressed files already exist
if [[ -f "$R1" && -f "$R2" ]]; then
    echo "  FASTQ files already present, skipping download."
else
    # wget -c enables resume of partial downloads; --tries for retry on transient failures
    wget -c --tries=3 --timeout=60 -O "$R1_GZ" "$URL_R1"
    wget -c --tries=3 --timeout=60 -O "$R2_GZ" "$URL_R2"

    # Decompress, keeping originals (-k) so re-runs don't need re-download
    gunzip -kf "$R1_GZ" "$R2_GZ"
fi

echo "Data download complete."
echo "-------------------------"
