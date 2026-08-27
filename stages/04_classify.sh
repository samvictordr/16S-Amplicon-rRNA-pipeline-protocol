#!/bin/bash
# stages/04_classify.sh — Stage 4: Assign taxonomy using SILVA Naive Bayes classifier
#
# Usage: bash stages/04_classify.sh [options]
#   --classifier-file FILE  Local filename for the classifier artifact
#   --classifier-url URL    Download URL if classifier file is absent
#   --download-only         Fetch the classifier and exit; skips classification
#                           and does not require rep-seqs.qza to exist yet
#   --threads N             Number of threads for classification (default: 1)
#   --output-dir DIR        Directory containing rep-seqs.qza and for output
#   --qiime-env ENV         Conda environment name for QIIME 2
#   --config FILE           Load configuration from file
#   --help                  Show this message
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../lib/common.sh"

show_help "$0" "$@"

[ -f "$SCRIPT_DIR/../pipeline.conf" ] && load_config "$SCRIPT_DIR/../pipeline.conf"

# Pull --download-only out before parse_args, which rejects unknown flags.
DOWNLOAD_ONLY=0
FILTERED_ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "--download-only" ]]; then
        DOWNLOAD_ONLY=1
    else
        FILTERED_ARGS+=("$arg")
    fi
done
parse_args ${FILTERED_ARGS+"${FILTERED_ARGS[@]}"}

# Verify input (not needed when we are only fetching the reference)
if [[ $DOWNLOAD_ONLY -eq 0 && ! -f "${OUTPUT_DIR}/rep-seqs.qza" ]]; then
    echo "Error: ${OUTPUT_DIR}/rep-seqs.qza not found." >&2
    echo "" >&2
    echo "Run the earlier stages first:" >&2
    echo "  bash stages/01_download.sh" >&2
    echo "  bash stages/02_import.sh" >&2
    echo "  bash stages/03_denoise.sh --threads 40" >&2
    echo "" >&2
    echo "Or, to only fetch the SILVA classifier without classifying anything:" >&2
    echo "  bash stages/04_classify.sh --download-only" >&2
    exit 1
fi

echo "STEP 4: Assigning taxonomy..."

CLASSIFIER_PATH="${OUTPUT_DIR}/${CLASSIFIER_FILE}"

# Download the classifier only if it isn't already cached
if [ ! -f "$CLASSIFIER_PATH" ]; then
    echo "  Downloading pre-trained SILVA classifier..."
    echo "  URL: $CLASSIFIER_URL"
    check_dependencies wget
    wget -c --tries=3 --timeout=120 -O "$CLASSIFIER_PATH" "$CLASSIFIER_URL"
else
    echo "  Using cached classifier: $CLASSIFIER_PATH"
fi

if [[ $DOWNLOAD_ONLY -eq 1 ]]; then
    echo "Classifier ready: $CLASSIFIER_PATH"
    echo "Verify its primer region before using it:"
    echo "  unzip -p \"$CLASSIFIER_PATH\" '*/provenance/artifacts/*/action/action.yaml' | grep -A1 primer"
    echo ""
    echo "Either of these is correct — both target V4 (E. coli 515-806, ~253 bp):"
    echo "  GTGCCAGCMGCCGCGGTAA / GGACTACHVGGGTWTCTAAT   (Caporaso 515F/806R)"
    echo "  GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT   (Parada/Apprill, EMP)"
    echo ""
    echo "WRONG: CCTACGGGNGGCWGCAG / GGACTACNVGGGTWTCTAAT is 341F/805R (V3-V4)."
    exit 0
fi

activate_qiime

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_PATH" \
  --i-reads "${OUTPUT_DIR}/rep-seqs.qza" \
  --p-n-jobs "$THREADS" \
  --o-classification "${OUTPUT_DIR}/taxonomy.qza"

echo "Taxonomy assignment complete. Output: ${OUTPUT_DIR}/taxonomy.qza"
echo "-------------------------"
