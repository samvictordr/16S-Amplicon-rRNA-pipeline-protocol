#!/bin/bash
# stages/02_import.sh — Stage 2: Import FASTQ reads into QIIME 2
#
# Usage: bash stages/02_import.sh [options]
#   --sample-id ID      EBI accession ID (used to locate FASTQ files)
#   --sample-name NAME  Sample label written into the QIIME 2 manifest
#   --output-dir DIR    Directory containing FASTQ files and for output artifacts
#   --qiime-env ENV     Conda environment name for QIIME 2
#   --config FILE       Load configuration from file
#   --help              Show this message
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"

show_help "$0" "$@"

[ -f "$SCRIPT_DIR/../pipeline.conf" ] && load_config "$SCRIPT_DIR/../pipeline.conf"
parse_args "$@"

echo "STEP 2: Importing data into QIIME 2..."
echo "  Sample name: $SAMPLE_NAME"
echo "  QIIME env:   $QIIME_ENV"

# Verify input files exist before activating conda (faster feedback).
#
# The manifest must reference GZIPPED FASTQs. q2-types opens manifest-referenced
# files with gzip unconditionally (see validate_paired_ends_equal_record_count),
# so a plain .fastq fails import with "BadGzipFile: Not a gzipped file (b'@E')".
# Stage 01 keeps the .gz alongside the decompressed copy; if only the plain file
# is present we recompress it here rather than making the user re-download.
ABS_OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

for n in 1 2; do
    plain="${ABS_OUTPUT_DIR}/${SAMPLE_ID}_${n}.fastq"
    gz="${plain}.gz"
    if [[ -f "$gz" ]]; then
        continue
    elif [[ -f "$plain" ]]; then
        echo "  Compressing ${SAMPLE_ID}_${n}.fastq for import..."
        gzip -kf "$plain"
    else
        echo "Error: neither $gz nor $plain was found." >&2
        echo "Run stage 01_download.sh first, or check --output-dir and --sample-id." >&2
        exit 1
    fi
done

activate_qiime

# Write the manifest file
MANIFEST="${OUTPUT_DIR}/manifest.tsv"
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANIFEST"
printf "%s\t%s/%s_1.fastq.gz\t%s/%s_2.fastq.gz\n" \
    "$SAMPLE_NAME" \
    "$ABS_OUTPUT_DIR" "$SAMPLE_ID" \
    "$ABS_OUTPUT_DIR" "$SAMPLE_ID" >> "$MANIFEST"

# Import into QIIME 2 artifact
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path "${OUTPUT_DIR}/demux.qza" \
  --input-format PairedEndFastqManifestPhred33V2

echo "QIIME 2 import complete. Artifact: ${OUTPUT_DIR}/demux.qza"
echo "-------------------------"
