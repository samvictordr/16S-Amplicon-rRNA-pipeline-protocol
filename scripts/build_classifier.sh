#!/bin/bash
# ==============================================================================
# scripts/build_classifier.sh — Train a 515F-806R SILVA classifier with RESCRIPt
# ==============================================================================
# Use this when the pre-built artifact from data.qiime2.org is rejected by your
# QIIME 2 install with a scikit-learn version mismatch, or when you want the
# classifier's provenance to be reproducible from within this repository.
#
# Both primer sets below target the SAME V4 region (E. coli 515-806, ~253 bp).
# They differ only in degeneracy, which changes how much of the SILVA reference
# is recovered during extract-reads:
#
#   classic  515F  GTGCCAGCMGCCGCGGTAA    806R   GGACTACHVGGGTWTCTAAT
#            Caporaso et al. 2011/2012. What the official QIIME 2 artifact
#            (silva-138-99-515-806-nb-classifier.qza) is built with.
#
#   emp      515F-Y GTGYCAGCMGCCGCGGTAA   806RB  GGACTACNVGGGTWTCTAAT
#            Parada et al. 2016 / Apprill et al. 2015. Broader coverage of
#            SAR11 and Thaumarchaeota; preferred for marine and soil work.
#
# Default is `classic`, so a locally rebuilt classifier matches the downloaded
# one. Switch to `emp` if your gradient includes marine or soil datasets, but
# then use it consistently across every dataset in the study.
#
# The classifier previously shipped with this repo was built with 341F / 805R
# (V3-V4) despite its filename, which is the defect this script exists to fix.
#
# Runtime: ~30-60 min on 40 cores. Peak memory ~30 GB during fit.
#
# Usage:
#   bash scripts/build_classifier.sh
#   bash scripts/build_classifier.sh --threads 40 --primer-set emp
#
# Options:
#   --threads N       Cores for extraction and fitting (default: physical cores)
#   --primer-set SET  classic (default) or emp
#   --output FILE     Output classifier artifact (default: from pipeline.conf)
#   --keep-intermediates  Retain the SILVA reference artifacts after training
#   --help            Show this message
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

show_help "$0" "$@"
[ -f "$REPO_ROOT/pipeline.conf" ] && load_config "$REPO_ROOT/pipeline.conf"

PRIMER_SET="classic"
NTHREADS=""
OUTPUT="$REPO_ROOT/${CLASSIFIER_FILE}"
KEEP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)     NTHREADS="$2";   shift 2 ;;
        --primer-set)  PRIMER_SET="$2"; shift 2 ;;
        --output)      OUTPUT="$2";     shift 2 ;;
        --keep-intermediates) KEEP=1; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

case "$PRIMER_SET" in
    classic)
        F_PRIMER="GTGCCAGCMGCCGCGGTAA"    # 515F,  Caporaso et al. 2011
        R_PRIMER="GGACTACHVGGGTWTCTAAT"   # 806R,  Caporaso et al. 2011
        ;;
    emp)
        F_PRIMER="GTGYCAGCMGCCGCGGTAA"    # 515F-Y, Parada et al. 2016
        R_PRIMER="GGACTACNVGGGTWTCTAAT"   # 806RB,  Apprill et al. 2015
        ;;
    *)
        echo "Error: --primer-set must be 'classic' or 'emp', got '$PRIMER_SET'" >&2
        exit 1
        ;;
esac

if [[ -z "$NTHREADS" ]]; then
    sockets=$(lscpu | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2}')
    cps=$(lscpu | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2); print $2}')
    NTHREADS=$(( sockets * cps ))
fi

WORK="$REPO_ROOT/.classifier-build"
mkdir -p "$WORK"

activate_qiime

if ! qiime rescript --help &>/dev/null; then
    echo "Error: the RESCRIPt plugin is not installed in '$QIIME_ENV'." >&2
    echo "Install it with:  conda install -c bioconda -c conda-forge -c qiime2 q2-rescript" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Pre-flight: report each action's output options for THIS QIIME 2 version.
#
# 2026.x promoted several outputs to required, and each one only surfaces when
# its step runs -- potentially an hour into a two-hour build. Listing them up
# front costs ~90 s and turns "it died at step 4" into a visible diff.
# ------------------------------------------------------------------------------
echo "Pre-flight: output options for this QIIME 2 build..."
for spec in \
    "rescript:get-silva-data" \
    "rescript:reverse-transcribe" \
    "rescript:cull-seqs" \
    "feature-classifier:extract-reads" \
    "rescript:dereplicate" \
    "feature-classifier:fit-classifier-naive-bayes"
do
    plugin="${spec%%:*}"; action="${spec##*:}"
    outs=$(qiime "$plugin" "$action" --help 2>/dev/null \
             | grep -oE '\-\-o-[a-z-]+' | sort -u | tr '\n' ' ')
    printf "  %-45s %s\n" "${plugin} ${action}" "${outs:-<unavailable>}"
done
echo "  (any --o- flag listed above must be passed; this script supplies the"
echo "   ones it knows about and probes for the version-dependent extras)"
echo

echo "=================================================="
echo " Building 515F-806R (V4) SILVA classifier"
echo "   primer set:     $PRIMER_SET"
echo "   forward primer: $F_PRIMER"
echo "   reverse primer: $R_PRIMER"
echo "   threads:        $NTHREADS"
echo "   output:         $OUTPUT"
echo "=================================================="
echo " Record the primer set in your Methods section. '515F-806R' alone is"
echo " ambiguous between the Caporaso and Parada/Apprill variants."
echo "=================================================="

# --- 1. Fetch SILVA 138.2 SSU NR99 -------------------------------------------
if [[ ! -f "$WORK/silva-seqs.qza" ]]; then
    echo "[1/6] Downloading SILVA 138.2 SSU NR99 (this takes a while)..."
    qiime rescript get-silva-data \
        --p-version '138.2' \
        --p-target 'SSURef_NR99' \
        --o-silva-sequences "$WORK/silva-seqs-rna.qza" \
        --o-silva-taxonomy "$WORK/silva-tax.qza"

    echo "[2/6] Reverse-transcribing RNA to DNA..."
    qiime rescript reverse-transcribe \
        --i-rna-sequences "$WORK/silva-seqs-rna.qza" \
        --o-dna-sequences "$WORK/silva-seqs.qza"
else
    echo "[1-2/6] Using cached SILVA reference."
fi

# --- 2. Quality-filter the reference -----------------------------------------
if [[ ! -f "$WORK/silva-seqs-filt.qza" ]]; then
    echo "[3/6] Culling low-quality reference sequences..."
    qiime rescript cull-seqs \
        --i-sequences "$WORK/silva-seqs.qza" \
        --p-n-jobs "$NTHREADS" \
        --o-clean-sequences "$WORK/silva-seqs-filt.qza"
fi

# --- 3. Extract the V4 amplicon region ---------------------------------------
# This is the step that decides the classifier's region. Getting the primers
# wrong here is invisible downstream except as degraded confidence scores.
if [[ ! -f "$WORK/silva-v4.qza" ]]; then
    echo "[4/6] Extracting 515F-806R region..."
    # QIIME 2 2026.x requires --o-read-extraction-stats here; older ones reject it.
    q2_extra_output feature-classifier extract-reads \
        --o-read-extraction-stats "$WORK/extraction-stats.qza"
    qiime feature-classifier extract-reads \
        --i-sequences "$WORK/silva-seqs-filt.qza" \
        --p-f-primer "$F_PRIMER" \
        --p-r-primer "$R_PRIMER" \
        --p-n-jobs "$NTHREADS" \
        --o-reads "$WORK/silva-v4.qza" \
        ${Q2_EXTRA_OUTPUTS[@]+"${Q2_EXTRA_OUTPUTS[@]}"}
fi

# --- 4. Dereplicate ----------------------------------------------------------
if [[ ! -f "$WORK/silva-v4-derep.qza" ]]; then
    echo "[5/6] Dereplicating extracted region..."
    qiime rescript dereplicate \
        --i-sequences "$WORK/silva-v4.qza" \
        --i-taxa "$WORK/silva-tax.qza" \
        --p-mode 'uniq' \
        --p-threads "$NTHREADS" \
        --o-dereplicated-sequences "$WORK/silva-v4-derep.qza" \
        --o-dereplicated-taxa "$WORK/silva-v4-tax.qza"
fi

# --- 5. Fit the Naive Bayes classifier ---------------------------------------
echo "[6/6] Fitting Naive Bayes classifier..."
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$WORK/silva-v4-derep.qza" \
    --i-reference-taxonomy "$WORK/silva-v4-tax.qza" \
    --o-classifier "$OUTPUT"

echo
echo "Classifier written to: $OUTPUT"
echo
echo "Verify the region recorded in its provenance:"
echo "  unzip -p \"$OUTPUT\" '*/provenance/artifacts/*/action/action.yaml' | grep -A1 primer"
echo

if [[ $KEEP -eq 0 ]]; then
    echo "Removing intermediates (pass --keep-intermediates to retain them)..."
    rm -rf "$WORK"
fi
