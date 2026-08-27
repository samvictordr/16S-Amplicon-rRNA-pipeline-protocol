#!/bin/bash
# ==============================================================================
# scripts/prepare_datasets.sh — Stage multiple public 16S datasets for benchmarking
# ==============================================================================
# Downloads, imports and denoises every accession listed in datasets.tsv, then
# reports the ASV richness of each. Richness is the classifier's real input
# dimension, so a spread of communities is what makes the scaling study
# generalizable beyond a single sample.
#
# Usage:
#   bash scripts/prepare_datasets.sh                      # all datasets
#   bash scripts/prepare_datasets.sh --datasets my.tsv
#   bash scripts/prepare_datasets.sh --threads 40
#   bash scripts/prepare_datasets.sh --list-only          # show plan, download nothing
#
# Options:
#   --datasets FILE   TSV of accessions (default: datasets.tsv)
#   --scratch DIR     Where artifacts are staged (default: .bench-scratch)
#   --threads N       Cores for the one-off denoise pass (default: physical cores)
#   --qiime-env ENV   Conda environment ("current" = use the active one)
#   --list-only       Print the dataset table and exit
#   --help            Show this message
#
# datasets.tsv columns (tab-separated, '#' comments allowed):
#   accession   label   biome   trunc_f   trunc_r
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

show_help "$0" "$@"
[ -f "$REPO_ROOT/pipeline.conf" ] && load_config "$REPO_ROOT/pipeline.conf"

DATASETS="$REPO_ROOT/datasets.tsv"
SCRATCH="$REPO_ROOT/.bench-scratch"
NTHREADS=""
LIST_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --datasets)   DATASETS="$2";  shift 2 ;;
        --scratch)    SCRATCH="$2";   shift 2 ;;
        --threads)    NTHREADS="$2";  shift 2 ;;
        --qiime-env)  QIIME_ENV="$2"; shift 2 ;;
        --list-only)  LIST_ONLY=1;    shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -f "$DATASETS" ]]; then
    echo "Error: dataset table not found: $DATASETS" >&2
    exit 1
fi

if [[ -z "$NTHREADS" ]]; then
    sockets=$(lscpu | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2}')
    cps=$(lscpu | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2); print $2}')
    NTHREADS=$(( sockets * cps ))
fi

echo "Dataset table: $DATASETS"
printf "%-14s %-22s %-16s %s\n" ACCESSION LABEL BIOME TRUNC
grep -vE '^\s*(#|$)' "$DATASETS" | while IFS=$'\t' read -r acc label biome tf tr; do
    printf "%-14s %-22s %-16s %s/%s\n" "$acc" "$label" "$biome" "$tf" "$tr"
done
echo

[[ $LIST_ONLY -eq 1 ]] && exit 0

mkdir -p "$SCRATCH"
activate_qiime

SUMMARY="${SCRATCH}/dataset_richness.tsv"
printf "accession\tlabel\tbiome\treads\tasvs\tmedian_asv_bp\tregion\n" > "$SUMMARY"

# Read the table on FD 3, not stdin. qiime and wget both read stdin, and if the
# loop shares it they swallow the remaining rows — the sweep then silently runs
# on one dataset instead of all of them.
while IFS=$'\t' read -r acc label biome tf tr <&3; do
    [[ -z "$acc" || "$acc" =~ ^[[:space:]]*# ]] && continue

    echo "=================================================="
    echo " $acc — $label ($biome)"
    echo "=================================================="

    # Keep the FASTQs gzipped. q2-types opens manifest-referenced files with
    # gzip unconditionally, so a plain .fastq fails import with BadGzipFile.
    r1="${SCRATCH}/${acc}_1.fastq.gz"
    r2="${SCRATCH}/${acc}_2.fastq.gz"

    # --- download ---
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        url1="$(derive_ebi_url "$acc" 1)"
        url2="$(derive_ebi_url "$acc" 2)"
        echo "  Downloading..."
        wget -c -q --show-progress --tries=3 --timeout=120 -O "$r1" "$url1"
        wget -c -q --show-progress --tries=3 --timeout=120 -O "$r2" "$url2"
    else
        echo "  Using cached FASTQ."
    fi

    reads=$(( $(gunzip -c "$r1" | wc -l) / 4 ))
    echo "  Reads: $reads"

    # DADA2 aborts if a truncation length exceeds the actual read length, and
    # read length varies between studies. Cap instead of failing partway through
    # a multi-dataset run. head closes the pipe early, so contain pipefail.
    len_f=$( ( set +o pipefail; gunzip -c "$r1" | head -2 | tail -1 | tr -d '\n' | wc -c ) )
    len_r=$( ( set +o pipefail; gunzip -c "$r2" | head -2 | tail -1 | tr -d '\n' | wc -c ) )
    len_f=$(( len_f )); len_r=$(( len_r ))
    echo "  Read length: ${len_f} / ${len_r} bp"
    if [[ "$tf" -gt "$len_f" ]]; then
        echo "  Capping trunc_f ${tf} -> ${len_f}"
        tf="$len_f"
    fi
    if [[ "$tr" -gt "$len_r" ]]; then
        echo "  Capping trunc_r ${tr} -> ${len_r}"
        tr="$len_r"
    fi

    # The V4 amplicon is ~253 bp; DADA2 needs >=12 bp of overlap to merge, and
    # 20 bp is the practical floor. Without it, merging silently returns almost
    # nothing and the sample looks empty rather than misconfigured.
    if (( tf + tr < 253 + 20 )); then
        echo "  WARNING: trunc_f + trunc_r = $(( tf + tr )) bp leaves under 20 bp" >&2
        echo "           of overlap for a ~253 bp V4 amplicon. Merging will mostly" >&2
        echo "           fail and this sample will yield few or no ASVs." >&2
    fi

    # --- import ---
    demux="${SCRATCH}/sample_${acc}.qza"
    if [[ ! -f "$demux" ]]; then
        manifest="${SCRATCH}/manifest_${acc}.tsv"
        printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$manifest"
        printf "%s\t%s\t%s\n" "$label" "$r1" "$r2" >> "$manifest"
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "$manifest" \
            --output-path "$demux" \
            --input-format PairedEndFastqManifestPhred33V2 >/dev/null
    fi

    # --- denoise once, to learn richness ---
    repseqs="${SCRATCH}/repseqs_${acc}.qza"
    if [[ ! -f "$repseqs" ]]; then
        echo "  Denoising (trunc ${tf}/${tr}, ${NTHREADS} threads)..."
        tmp="${SCRATCH}/prep_${acc}"; rm -rf "$tmp"; mkdir -p "$tmp"
        dada2_extra_outputs "$tmp"
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$demux" \
            --p-trunc-len-f "$tf" --p-trunc-len-r "$tr" \
            --p-n-threads "$NTHREADS" \
            --o-table "${tmp}/table.qza" \
            --o-representative-sequences "$repseqs" \
            --o-denoising-stats "${tmp}/stats.qza" \
            ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"} >/dev/null
    fi

    # ASV count plus median amplicon length. The length is a primer-region check:
    # 515F-806R (V4) lands near 253 bp, 341F-805R (V3-V4) near 425-465 bp. A
    # dataset whose region disagrees with the classifier must not be pooled in.
    read -r asvs medlen region < <(python3 - "$repseqs" <<'PY'
import sys, zipfile, re, statistics
z = zipfile.ZipFile(sys.argv[1])
fa = [n for n in z.namelist() if n.endswith("dna-sequences.fasta")][0]
txt = z.read(fa).decode()
lens = [len(s.replace("\n", "")) for s in re.split(r">[^\n]*\n", txt)[1:]]
if not lens:
    print("0 0 EMPTY"); raise SystemExit
m = int(statistics.median(lens))
region = "V4" if 230 <= m <= 275 else ("V3-V4" if 400 <= m <= 480 else "OTHER")
print(len(lens), m, region)
PY
)
    echo "  ASVs: $asvs   median length: ${medlen} bp   region: $region"
    if [[ "$region" != "V4" ]]; then
        echo "  WARNING: amplicon looks like '$region', not V4. The 515F-806R"
        echo "           classifier will misclassify it. Exclude or reclassify." >&2
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$acc" "$label" "$biome" "$reads" "$asvs" "$medlen" "$region" >> "$SUMMARY"

done 3< <(grep -vE '^\s*(#|$)' "$DATASETS")

# --- report the richest set, but do NOT stage it as the sweep input ---
#
# An earlier revision copied the richest single sample to rich_rep-seqs.qza here.
# run_real_benchmarks.sh builds the pooled ASV set only when that file is absent,
# so this copy silently pre-empted the pooling: following the documented sequence
# ran the ASV sweep on one sample's ASVs instead of the cohort union, and the
# published sweep could not be reproduced from the released code.
richest=$(tail -n +2 "$SUMMARY" | sort -k5 -n -r | head -1)
rich_acc=$(echo "$richest" | cut -f1)
rich_asvs=$(echo "$richest" | cut -f5)

echo
echo "=================================================="
echo " Richness summary"
echo "=================================================="
column -t -s$'\t' "$SUMMARY"
echo
echo "Richest single sample: $rich_acc ($rich_asvs ASVs)"
echo "The classifier sweep uses the POOLED union across all samples, which"
echo "run_real_benchmarks.sh builds; it is larger than any single sample."
echo
echo "Next: bash scripts/run_real_benchmarks.sh --all"
