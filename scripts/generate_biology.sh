#!/bin/bash
# ==============================================================================
# scripts/generate_biology.sh — Taxonomic and diversity results for the cohort
# ==============================================================================
# The benchmarking arms treat the nine samples as a workload. This script
# produces the biological result they actually contain: taxonomic composition
# and alpha diversity per sample and per body site.
#
# It also runs the check that ties the two halves of the paper together —
# whether the taxonomy assigned to a sample changes with the thread count used
# to produce it. The determinism arm in run_real_benchmarks.sh compares
# representative SEQUENCES; this compares the resulting taxonomic ASSIGNMENTS,
# which is what a practitioner acting on the resource recommendations cares
# about.
#
# Usage:
#   bash scripts/generate_biology.sh                 # all samples in datasets.tsv
#   bash scripts/generate_biology.sh --threads 40
#   bash scripts/generate_biology.sh --skip-determinism
#
# Options:
#   --threads N          Cores for classification (default: physical cores)
#   --scratch DIR        Staging directory (default: .bench-scratch)
#   --outdir DIR         Where results are written (default: biology)
#   --skip-determinism   Skip the taxonomy-vs-thread-count check
#   --help               Show this message
#
# Expects scripts/prepare_datasets.sh to have run first.
# Runtime: roughly 10-15 min for nine samples, plus ~3 min for the
# determinism check.
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$REPO_ROOT/lib/common.sh"

show_help "$0" "$@"
[ -f "$REPO_ROOT/pipeline.conf" ] && load_config "$REPO_ROOT/pipeline.conf"

SCRATCH="$REPO_ROOT/.bench-scratch"
OUTDIR="$REPO_ROOT/biology"
NTHREADS=""
DO_DETERMINISM=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)          NTHREADS="$2"; shift 2 ;;
        --scratch)          SCRATCH="$2";  shift 2 ;;
        --outdir)           OUTDIR="$2";   shift 2 ;;
        --skip-determinism) DO_DETERMINISM=0; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$NTHREADS" ]]; then
    sockets=$(lscpu | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2}')
    cps=$(lscpu | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2); print $2}')
    NTHREADS=$(( sockets * cps ))
fi

mkdir -p "$OUTDIR" "$SCRATCH"
activate_qiime

CLASSIFIER_PATH="${REPO_ROOT}/${CLASSIFIER_FILE}"
[[ -f "$CLASSIFIER_PATH" ]] || {
    echo "Error: classifier not found at $CLASSIFIER_PATH" >&2
    echo "Run: bash stages/04_classify.sh --download-only" >&2
    exit 1
}

COMP="${OUTDIR}/composition.tsv"
ALPHA="${OUTDIR}/alpha_diversity.tsv"
printf "accession\tlabel\tbody_site\tlevel\ttaxon\treads\trel_abundance\n" > "$COMP"
printf "accession\tlabel\tbody_site\treads\tobserved_asvs\tshannon\tsimpson_evenness\tdominant_taxon\tdominant_frac\n" > "$ALPHA"

# ------------------------------------------------------------------------------
# Per-sample taxonomy and diversity
# ------------------------------------------------------------------------------
while IFS=$'\t' read -r acc label biome tf tr <&3; do
    [[ -z "$acc" || "$acc" =~ ^[[:space:]]*# ]] && continue

    repseqs="${SCRATCH}/repseqs_${acc}.qza"
    table="${SCRATCH}/prep_${acc}/table.qza"
    if [[ ! -f "$repseqs" ]]; then
        echo "  skip ${acc}: ${repseqs} missing (run prepare_datasets.sh)" >&2
        continue
    fi

    # prepare_datasets.sh writes the table to a temp dir it may have cleared.
    # Regenerate only if it is genuinely absent.
    if [[ ! -f "$table" ]]; then
        echo "  ${acc}: feature table absent, re-denoising..."
        tmp="${SCRATCH}/prep_${acc}"; mkdir -p "$tmp"
        dada2_extra_outputs "$tmp"
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "${SCRATCH}/sample_${acc}.qza" \
            --p-trunc-len-f "$tf" --p-trunc-len-r "$tr" \
            --p-n-threads "$NTHREADS" \
            --o-table "$table" \
            --o-representative-sequences "${tmp}/rep-seqs.qza" \
            --o-denoising-stats "${tmp}/stats.qza" \
            ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"} >/dev/null
    fi

    echo "=== ${acc} (${label}, ${biome}) ==="
    tax="${SCRATCH}/tax_${acc}.qza"
    if [[ ! -f "$tax" ]]; then
        qiime feature-classifier classify-sklearn \
            --i-classifier "$CLASSIFIER_PATH" \
            --i-reads "$repseqs" \
            --p-n-jobs "$NTHREADS" \
            --o-classification "$tax" >/dev/null
    fi

    exp="${SCRATCH}/exp_${acc}"; rm -rf "$exp"; mkdir -p "$exp"
    qiime tools export --input-path "$tax"   --output-path "$exp" >/dev/null
    qiime tools export --input-path "$table" --output-path "$exp" >/dev/null
    biom convert -i "$exp/feature-table.biom" -o "$exp/feature-table.tsv" --to-tsv 2>/dev/null

    python3 - "$exp" "$acc" "$label" "$biome" "$COMP" "$ALPHA" <<'PY'
import sys, csv, math, os
exp, acc, label, biome, comp_out, alpha_out = sys.argv[1:7]

# taxonomy: FeatureID -> lineage
tax = {}
with open(os.path.join(exp, "taxonomy.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        tax[row["Feature ID"]] = row["Taxon"]

# abundance: FeatureID -> count (single-sample table)
counts = {}
with open(os.path.join(exp, "feature-table.tsv")) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            try:
                counts[parts[0]] = float(parts[1])
            except ValueError:
                pass

total = sum(counts.values()) or 1.0

def rank(lineage, idx):
    """Return the taxon at a rank, falling back to the deepest named level."""
    parts = [p.strip() for p in lineage.split(";")]
    if idx < len(parts):
        p = parts[idx]
        if p and not p.endswith("__"):
            return p
    for p in reversed(parts):
        if p and not p.endswith("__"):
            return p + " (unresolved)"
    return "Unassigned"

with open(comp_out, "a") as fh:
    for level_name, idx in (("phylum", 1), ("genus", 5)):
        agg = {}
        for fid, n in counts.items():
            agg[rank(tax.get(fid, "Unassigned"), idx)] = agg.get(
                rank(tax.get(fid, "Unassigned"), idx), 0.0) + n
        for taxon, n in sorted(agg.items(), key=lambda kv: -kv[1]):
            fh.write(f"{acc}\t{label}\t{biome}\t{level_name}\t{taxon}\t{n:.0f}\t{n/total:.6f}\n")

# alpha diversity from the observed counts
obs = sum(1 for n in counts.values() if n > 0)
props = [n / total for n in counts.values() if n > 0]
shannon = -sum(p * math.log(p) for p in props) if props else 0.0
simpson = 1.0 / sum(p * p for p in props) / obs if obs and props else 0.0
gen = {}
for fid, n in counts.items():
    gen[rank(tax.get(fid, "Unassigned"), 5)] = gen.get(rank(tax.get(fid, "Unassigned"), 5), 0.0) + n
dom, domn = max(gen.items(), key=lambda kv: kv[1]) if gen else ("NA", 0.0)

with open(alpha_out, "a") as fh:
    fh.write(f"{acc}\t{label}\t{biome}\t{total:.0f}\t{obs}\t{shannon:.4f}\t"
             f"{simpson:.4f}\t{dom}\t{domn/total:.4f}\n")
print(f"  {obs} ASVs, Shannon {shannon:.2f}, dominant {dom} ({100*domn/total:.1f}%)")
PY
done 3< <(grep -vE '^\s*(#|$)' "$REPO_ROOT/datasets.tsv")

# ------------------------------------------------------------------------------
# Does the taxonomy change with thread count?
#
# The determinism arm compares representative sequences. This compares the
# assignments themselves, which is the claim a practitioner needs: that adopting
# the resource recommendations does not change the community you report.
# ------------------------------------------------------------------------------
if [[ $DO_DETERMINISM -eq 1 ]]; then
    echo
    echo "=== Does taxonomic assignment depend on thread count? ==="
    DET="${OUTDIR}/taxonomy_determinism.tsv"
    printf "threads\tn_assignments\tassignments_sha256\n" > "$DET"
    ref="${SCRATCH}/repseqs_ERR3444605.qza"
    [[ -f "$ref" ]] || ref="${REPO_ROOT}/rep-seqs.qza"

    for p in 1 4 "$NTHREADS"; do
        out="${SCRATCH}/taxdet_${p}"; rm -rf "$out"; mkdir -p "$out"
        qiime feature-classifier classify-sklearn \
            --i-classifier "$CLASSIFIER_PATH" --i-reads "$ref" \
            --p-n-jobs "$p" --o-classification "${out}/tax.qza" >/dev/null
        qiime tools export --input-path "${out}/tax.qza" --output-path "$out" >/dev/null
        read -r n sha < <(python3 - "${out}/taxonomy.tsv" <<'PY2'
import sys, csv, hashlib
rows = []
with open(sys.argv[1]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.append(f"{r['Feature ID']}\t{r['Taxon']}")
rows.sort()
print(len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest())
PY2
)
        printf "%s\t%s\t%s\n" "$p" "$n" "$sha" >> "$DET"
        echo "  p=${p}: ${n} assignments, sha256 ${sha:0:16}..."
    done

    if [[ $(cut -f3 "$DET" | tail -n +2 | sort -u | wc -l) -eq 1 ]]; then
        echo "  RESULT: taxonomic assignments identical at every thread count."
    else
        echo "  RESULT: assignments DIFFER between thread counts." >&2
        echo "  The resource recommendations change the reported community and" >&2
        echo "  must be presented as a speed/reproducibility trade, not a saving." >&2
    fi
fi

echo
echo "Written to ${OUTDIR}/:"
echo "  composition.tsv          per-sample phylum and genus composition"
echo "  alpha_diversity.tsv      observed ASVs, Shannon, evenness, dominant taxon"
[[ $DO_DETERMINISM -eq 1 ]] && echo "  taxonomy_determinism.tsv assignments vs thread count"
echo
echo "Next: python3 paper/analysis/generate_biology_figures.py"
