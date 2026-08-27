#!/bin/bash
# ==============================================================================
# scripts/run_real_benchmarks.sh — Measured (not simulated) scaling experiments
# ==============================================================================
# Replaces scripts/simulate_benchmarks.py with real wall-clock measurements.
# Emits CSVs in the identical schema, so paper/analysis/generate_figures.py and
# lib/predictor.py consume the output unchanged.
#
# Experiments:
#   1. strong-denoise  DADA2, fixed input, cores 1..NPROC
#   2. weak-denoise    DADA2, reads scaled with cores (constant reads/core)
#   3. strong-classify classify-sklearn, fixed ASV set, cores 1..NPROC
#   4. asv-sweep       classify-sklearn vs. ASV count at fixed cores
#   4b. query-sweep    same, extended to 10,000 real amplicon queries
#   5. depth-pair      shallow vs deep sample, same body site  <-- key experiment
#   6. array           multi-sample: sequential loop vs. SLURM job array
#   7. numa            same core count, one socket vs. spanning sockets (optional)
#
# Usage:
#   bash scripts/run_real_benchmarks.sh --all
#   bash scripts/run_real_benchmarks.sh --experiment asv-sweep --reps 5
#   bash scripts/run_real_benchmarks.sh --experiment query-sweep --reps 3
#   bash scripts/run_real_benchmarks.sh --all --max-cores 40 --reps 3
#
# Options:
#   --experiment NAME  One of: strong-denoise, weak-denoise, strong-classify,
#                      asv-sweep, query-sweep, depth-pair, depth-series,
#                      per-sample, determinism, numa, array  (repeatable)
#   --all              Run every experiment in order
#   --reps N           Repetitions per measurement point (default: 3)
#   --max-cores N      Highest core count to test (default: physical cores)
#   --outdir DIR       Where CSVs are written (default: benchmarks)
#   --scratch DIR      Working directory for intermediates (default: .bench-scratch)
#   --qiime-env ENV    Conda environment (default: from pipeline.conf)
#   --dry-run          Print the measurement plan and exit
#   --help             Show this message
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$REPO_ROOT/lib/common.sh"
source "$REPO_ROOT/lib/benchmark.sh"

show_help "$0" "$@"
[ -f "$REPO_ROOT/pipeline.conf" ] && load_config "$REPO_ROOT/pipeline.conf"

# ------------------------------------------------------------------------------
# Argument parsing (this script has its own flags, not lib/common.sh's)
# ------------------------------------------------------------------------------
EXPERIMENTS=()
REPS=3
MAX_CORES=""
OUTDIR="$REPO_ROOT/benchmarks"
SCRATCH="$REPO_ROOT/.bench-scratch"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --experiment) EXPERIMENTS+=("$2"); shift 2 ;;
        --all)        EXPERIMENTS=(strong-denoise weak-denoise strong-classify asv-sweep depth-pair depth-series per-sample determinism array); shift ;;
        --reps)       REPS="$2";       shift 2 ;;
        --max-cores)  MAX_CORES="$2";  shift 2 ;;
        --outdir)     OUTDIR="$2";     shift 2 ;;
        --scratch)    SCRATCH="$2";    shift 2 ;;
        --qiime-env)  QIIME_ENV="$2";  shift 2 ;;
        --dry-run)    DRY_RUN=1;       shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ ${#EXPERIMENTS[@]} -eq 0 ]]; then
    echo "Error: specify --all or at least one --experiment NAME" >&2
    echo "Run with --help for the list." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Topology detection
#
# Scaling studies must use PHYSICAL cores. Hyperthread siblings share execution
# units, so counting them inflates p and makes efficiency look artificially bad.
# ------------------------------------------------------------------------------
detect_topology() {
    SOCKETS=$(lscpu | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2}')
    CORES_PER_SOCKET=$(lscpu | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2); print $2}')
    THREADS_PER_CORE=$(lscpu | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2}')
    PHYSICAL_CORES=$(( SOCKETS * CORES_PER_SOCKET ))
    LOGICAL_CPUS=$(( PHYSICAL_CORES * THREADS_PER_CORE ))
    [[ -n "$MAX_CORES" ]] || MAX_CORES=$PHYSICAL_CORES

    if [[ "$MAX_CORES" -gt "$PHYSICAL_CORES" ]]; then
        echo "Warning: --max-cores $MAX_CORES exceeds $PHYSICAL_CORES physical cores." >&2
        echo "         Points above $PHYSICAL_CORES measure hyperthreads, not cores." >&2
    fi
}

# Build the core ladder: powers of two, then socket-aligned steps.
build_core_ladder() {
    CORE_LADDER=()
    local c=1
    while [[ $c -le $MAX_CORES ]]; do
        CORE_LADDER+=("$c")
        c=$(( c * 2 ))
    done
    # Add socket-aligned points so NUMA boundaries are sampled explicitly
    for extra in "$CORES_PER_SOCKET" "$MAX_CORES"; do
        if [[ $extra -le $MAX_CORES ]] && [[ ! " ${CORE_LADDER[*]} " =~ " ${extra} " ]]; then
            CORE_LADDER+=("$extra")
        fi
    done
    IFS=$'\n' CORE_LADDER=($(sort -n -u <<<"${CORE_LADDER[*]}")); unset IFS
}

# ------------------------------------------------------------------------------
# Environment hygiene
#
# QIIME 2 sits on numpy/scipy, which spawn their own BLAS thread pools. If those
# are left unbounded they oversubscribe the box and contaminate every timing.
# Pin them to 1 so --p-n-threads / --p-n-jobs is the only parallelism knob.
# ------------------------------------------------------------------------------
pin_blas_threads() {
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
}

# Drop page cache between reps so filesystem I/O is measured cold rather than
# served from RAM. Needs root; degrades gracefully to a warning.
drop_caches() {
    if sudo -n true 2>/dev/null; then
        sync && sudo -n sh -c 'echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null || true
    else
        CACHE_WARNED=1
    fi
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Read a FASTQ that may or may not be gzipped.
cat_fastq() {
    local base="$1"   # path without the .gz
    if [[ -f "${base}.gz" ]]; then gunzip -c "${base}.gz"; else cat "$base"; fi
}

# Deterministically take the first N read pairs. Both mates stay in sync because
# ENA writes _1 and _2 in identical order.
#
# Output is gzipped: q2-types opens manifest-referenced FASTQs with gzip
# unconditionally, so a plain .fastq fails import with BadGzipFile.
subsample_fastq() {
    local n_reads="$1" src_prefix="$2" dst_prefix="$3"
    local lines=$(( n_reads * 4 ))
    local m
    for m in 1 2; do
        # head exits early and SIGPIPEs gunzip, which trips pipefail. The
        # subshell contains the disabled pipefail so callers keep theirs.
        ( set +o pipefail
          cat_fastq "${src_prefix}_${m}.fastq" | head -n "$lines" | gzip -c \
              > "${dst_prefix}_${m}.fastq.gz" )
    done
}

count_fastq_reads() {
    local base="${1%.gz}"
    echo $(( $(cat_fastq "$base" | wc -l) / 4 ))
}

# Import a FASTQ pair into a demux artifact. Accepts a prefix whose mates are
# either <prefix>_N.fastq.gz or <prefix>_N.fastq, and guarantees the manifest
# points at gzipped files.
import_demux() {
    local prefix="$1" sample_name="$2" out_qza="$3"
    local manifest="${SCRATCH}/manifest_$(basename "$prefix").tsv"
    local m
    for m in 1 2; do
        if [[ ! -f "${prefix}_${m}.fastq.gz" ]]; then
            if [[ -f "${prefix}_${m}.fastq" ]]; then
                gzip -kf "${prefix}_${m}.fastq"
            else
                echo "Error: no FASTQ found for ${prefix}_${m}" >&2
                return 1
            fi
        fi
    done
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$manifest"
    printf "%s\t%s_1.fastq.gz\t%s_2.fastq.gz\n" "$sample_name" "$prefix" "$prefix" >> "$manifest"
    rm -f "$out_qza"
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$manifest" \
        --output-path "$out_qza" \
        --input-format PairedEndFastqManifestPhred33V2 >/dev/null
}

# Number of sequences in a FeatureData[Sequence] artifact.
count_asvs() {
    local qza="$1"
    python3 - "$qza" <<'PY'
import sys, zipfile
z = zipfile.ZipFile(sys.argv[1])
fa = [n for n in z.namelist() if n.endswith("dna-sequences.fasta")][0]
print(z.read(fa).decode().count(">"))
PY
}

# Write the first K sequences of a rep-seqs artifact into a new artifact.
subset_rep_seqs() {
    local src_qza="$1" k="$2" out_qza="$3"
    local tmp_fa="${SCRATCH}/subset_${k}.fasta"
    python3 - "$src_qza" "$k" "$tmp_fa" <<'PY'
import sys, zipfile, re
src, k, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
z = zipfile.ZipFile(src)
fa = [n for n in z.namelist() if n.endswith("dna-sequences.fasta")][0]
records = re.findall(r">[^\n]*\n(?:[^>]*)", z.read(fa).decode())
with open(out, "w") as fh:
    fh.write("".join(records[:k]))
PY
    rm -f "$out_qza"
    qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path "$tmp_fa" \
        --output-path "$out_qza" >/dev/null
}

# ------------------------------------------------------------------------------
# Measured stage invocations
# ------------------------------------------------------------------------------

run_denoise() {
    local threads="$1" reads="$2" demux="$3" csv="$4" tag="${5:-denoise}"
    local out="${SCRATCH}/out_${tag}_${threads}_${reads}"
    rm -rf "$out"; mkdir -p "$out"
    dada2_extra_outputs "$out"
    benchmark_stage "$tag" "$threads" "$reads" "$csv" \
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$demux" \
            --p-trunc-len-f "$TRUNC_LEN_F" \
            --p-trunc-len-r "$TRUNC_LEN_R" \
            --p-n-threads "$threads" \
            --o-table "${out}/table.qza" \
            --o-representative-sequences "${out}/rep-seqs.qza" \
            --o-denoising-stats "${out}/stats.qza" \
            ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"}
}

# NOTE: the `reads` column carries the ASV COUNT for classify rows. That is the
# classifier's real input size: classify-sklearn receives one representative
# sequence per ASV, not the read set, so ASV count is what varies its cost.
run_classify() {
    local jobs="$1" n_asvs="$2" repseqs="$3" csv="$4" tag="${5:-classify}"
    local out="${SCRATCH}/out_${tag}_${jobs}_${n_asvs}"
    rm -rf "$out"; mkdir -p "$out"
    benchmark_stage "$tag" "$jobs" "$n_asvs" "$csv" \
        qiime feature-classifier classify-sklearn \
            --i-classifier "$CLASSIFIER_PATH" \
            --i-reads "$repseqs" \
            --p-n-jobs "$jobs" \
            --o-classification "${out}/taxonomy.qza"
}

# ------------------------------------------------------------------------------
# Experiment 1: strong scaling, DADA2
# ------------------------------------------------------------------------------
exp_strong_denoise() {
    echo "=== Experiment: strong scaling (DADA2 denoise) ==="
    local csv="${OUTDIR}/strong_scaling.csv"
    local reads; reads=$(count_fastq_reads "${BASE_PREFIX}_1.fastq")
    for rep in $(seq 1 "$REPS"); do
        for p in "${CORE_LADDER[@]}"; do
            echo "  [rep $rep] denoise p=$p reads=$reads"
            drop_caches
            run_denoise "$p" "$reads" "$BASE_DEMUX" "$csv"
        done
    done
}

# ------------------------------------------------------------------------------
# Experiment 2: weak scaling, DADA2 (constant reads per core)
# ------------------------------------------------------------------------------
exp_weak_denoise() {
    echo "=== Experiment: weak scaling (DADA2 denoise) ==="
    local csv="${OUTDIR}/weak_scaling.csv"
    local total; total=$(count_fastq_reads "${BASE_PREFIX}_1.fastq")
    # Hold reads/core fixed at total/max_cores
    local per_core=$(( total / MAX_CORES ))
    echo "  Holding ${per_core} reads/core constant"
    for rep in $(seq 1 "$REPS"); do
        for p in "${CORE_LADDER[@]}"; do
            local n=$(( per_core * p ))
            [[ $n -gt $total ]] && n=$total
            local prefix="${SCRATCH}/weak_${n}"
            [[ -f "${prefix}_1.fastq.gz" ]] || subsample_fastq "$n" "$BASE_PREFIX" "$prefix"
            local demux="${SCRATCH}/weak_${n}.qza"
            [[ -f "$demux" ]] || import_demux "$prefix" "weak-${n}" "$demux"
            echo "  [rep $rep] denoise p=$p reads=$n"
            drop_caches
            run_denoise "$p" "$n" "$demux" "$csv"
        done
    done
}

# ------------------------------------------------------------------------------
# Experiment 3: strong scaling, classifier
# ------------------------------------------------------------------------------
exp_strong_classify() {
    echo "=== Experiment: strong scaling (classify-sklearn) ==="
    local csv="${OUTDIR}/strong_scaling.csv"
    local n_asvs; n_asvs=$(count_asvs "$BASE_REPSEQS")
    echo "  ASV set size: ${n_asvs}"
    if [[ $n_asvs -lt 500 ]]; then
        echo "  NOTE: ${n_asvs} ASVs is small. Expect the serial model load to"
        echo "        dominate and speedup to be near-flat. That is the finding,"
        echo "        not a bug. Run --experiment asv-sweep to quantify it."
    fi
    for rep in $(seq 1 "$REPS"); do
        for p in "${CORE_LADDER[@]}"; do
            echo "  [rep $rep] classify p=$p asvs=$n_asvs"
            drop_caches
            run_classify "$p" "$n_asvs" "$BASE_REPSEQS" "$csv"
        done
    done
}

# ------------------------------------------------------------------------------
# Experiment 4: classifier cost vs. ASV count   <-- the key experiment
#
# Isolates the classifier's true input dimension. k=1 measures the serial floor
# (model deserialization); the slope above it is the parallelizable work.
# ------------------------------------------------------------------------------
exp_asv_sweep() {
    echo "=== Experiment: classifier cost vs. ASV count ==="
    local csv="${OUTDIR}/asv_sweep.csv"
    local available; available=$(count_asvs "$RICH_REPSEQS")
    echo "  Richest available ASV set: ${available}"

    local ks=()
    for k in 1 10 25 50 100 250 500 1000 2500 5000 10000; do
        [[ $k -le $available ]] && ks+=("$k")
    done
    [[ " ${ks[*]} " =~ " ${available} " ]] || ks+=("$available")

    # Two core counts: serial floor, and a mid-range parallel setting.
    local mid=$(( CORES_PER_SOCKET < 8 ? CORES_PER_SOCKET : 8 ))
    for rep in $(seq 1 "$REPS"); do
        for p in 1 "$mid"; do
            for k in "${ks[@]}"; do
                local sub="${SCRATCH}/repseqs_${k}.qza"
                [[ -f "$sub" ]] || subset_rep_seqs "$RICH_REPSEQS" "$k" "$sub"
                echo "  [rep $rep] classify p=$p asvs=$k"
                drop_caches
                run_classify "$p" "$k" "$sub" "$csv" "classify_asvsweep"
            done
        done
    done
}

# ------------------------------------------------------------------------------
# Experiment 4b: where does per-query work overtake the fixed cost?
#
# asv-sweep is bounded by the ASV count the nine samples actually pool to (218),
# which sits well inside the flat region. Bokulich et al. (2018) fitted the same
# multinomial Naive Bayes implementation out to 10,000 queries at 23 ms each;
# between 218 and 10,000 the linear term must overtake the fixed cost, and at
# some query count parallelism must begin to pay. Neither paper locates that
# transition. This arm measures it directly.
#
# Denoising collapses the nine samples to 218 ASVs, so the query pool is built
# instead from the deep sample's reads before denoising. See build_query_pool
# below for how, and for why the queries are fixed-length.
#
# Sequences are taken at an even stride rather than as a prefix. Reads come off
# the flowcell in tile order, so a prefix would draw one physical region of the
# flowcell; striding spreads the draw across the whole run.
# ------------------------------------------------------------------------------
subset_rep_seqs_strided() {
    local src_qza="$1" k="$2" out_qza="$3"
    local tmp_fa="${SCRATCH}/qsubset_${k}.fasta"
    python3 - "$src_qza" "$k" "$tmp_fa" <<'PYSTRIDE'
import sys, zipfile, re
src, k, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
z = zipfile.ZipFile(src)
fa = [n for n in z.namelist() if n.endswith("dna-sequences.fasta")][0]
records = re.findall(r">[^\n]*\n(?:[^>]*)", z.read(fa).decode())
if k >= len(records):
    picked = records
else:
    stride = len(records) / k          # even stride over the whole artifact
    picked = [records[int(i * stride)] for i in range(k)]
with open(out, "w") as fh:
    fh.write("".join(picked))
print("  drew %d of %d pool sequences" % (len(picked), len(records)),
      file=sys.stderr)
PYSTRIDE
    rm -f "$out_qza"
    qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path "$tmp_fa" \
        --output-path "$out_qza" >/dev/null
}

# Build the query pool from real amplicon reads.
#
# Deliberately dependency-free: it reads the FASTQ out of the demux artifact and
# does the work in Python. The first version of this chained join-pairs,
# quality-filter and vsearch dereplicate-sequences, which failed on this host for
# reasons the redirected stderr hid. Those plugins are not needed -- a query set
# is a FASTA of distinct sequences, and the reads already are that.
#
# Queries are R1 reads trimmed to a fixed 250 bp. V4 with the Caporaso primers is
# ~253 bp, so an R1 read is very nearly the whole amplicon, and a fixed trim
# makes every query exactly the same length. That matters for this measurement:
# Naive Bayes scores k-mer profiles, so per-query cost tracks sequence length,
# and holding length constant means the sweep varies query COUNT and nothing
# else. Reads containing N are dropped.
#
# These are real amplicons from the study carrying real sequencing error, which
# is what a deeply sequenced community actually hands the classifier. Assignment
# quality is not a result here and is not recorded; only wall time and peak RSS.
# Whether a query is novel to the model does not change what it costs to score.
build_query_pool() {
    local acc="$1" out_qza="$2"
    local demux="${SCRATCH}/sample_${acc}.qza"
    local work="${SCRATCH}/querypool_${acc}"
    local fa="${work}/pool.fasta"

    if [[ ! -f "$demux" ]]; then
        echo "  Cannot build a query pool: ${demux} is missing." >&2
        echo "  Run scripts/prepare_datasets.sh first." >&2
        return 1
    fi

    mkdir -p "$work"
    echo "  Building query pool from ${acc}..."

    # NOTE: no >/dev/null here. If this breaks, the reason must be visible.
    if ! python3 - "$demux" "$fa" <<'PYPOOL'
import gzip, io, re, sys, zipfile

src, out = sys.argv[1], sys.argv[2]
WANT, MAX_READS, TRIM_CAP = 60000, 600000, 250

z = zipfile.ZipFile(src)
fq = [n for n in z.namelist()
      if n.endswith((".fastq.gz", ".fastq")) and "_R1_" in n]
if not fq:
    fq = [n for n in z.namelist() if n.endswith((".fastq.gz", ".fastq"))]
if not fq:
    sys.exit("  no FASTQ found inside %s" % src)
member = sorted(fq)[0]
print("    reading %s" % member.split("/")[-1])

def reader():
    raw = z.open(member)
    return io.TextIOWrapper(gzip.open(raw) if member.endswith(".gz") else raw)

# Trim length is measured, not assumed. This study is 2x251 bp so the cap
# applies, but a hardcoded length that silently discards every read is exactly
# the failure this function was rewritten to stop hiding.
probe = []
with reader() as fh:
    for i, line in enumerate(fh):
        if i % 4 == 1:
            probe.append(len(line.strip()))
        if len(probe) >= 1000:
            break
if not probe:
    sys.exit("    no reads found in %s" % member)
modal = max(set(probe), key=probe.count)
TRIM = min(TRIM_CAP, modal)
print("    modal read length %d bp -> trimming queries to %d bp" % (modal, TRIM))

fh = reader()
seen, order, n, short, ambig = set(), [], 0, 0, 0
for i, line in enumerate(fh):
    if i % 4 != 1:
        continue
    n += 1
    if n > MAX_READS or len(seen) >= WANT:
        break
    sq = line.strip()
    if len(sq) < TRIM:
        short += 1
        continue
    sq = sq[:TRIM]
    if "N" in sq:
        ambig += 1
        continue
    if sq not in seen:
        seen.add(sq)
        order.append(sq)

if len(order) < 10000:
    print("    WARNING: only %d distinct reads from %d scanned "
          "(%d too short, %d with N)" % (len(order), n, short, ambig),
          file=sys.stderr)

with open(out, "w") as f:
    for i, sq in enumerate(order):
        f.write(">q%06d\n%s\n" % (i, sq))
fh.close()
print("    %d distinct %d bp queries from %d reads "
      "(%d short, %d with N)" % (len(order), TRIM, n, short, ambig))
PYPOOL
    then
        echo "  Query pool extraction failed." >&2
        return 1
    fi

    [[ -s "$fa" ]] || { echo "  Query pool FASTA is empty: ${fa}" >&2; return 1; }

    rm -f "$out_qza"
    if ! qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path "$fa" \
        --output-path "$out_qza"; then
        echo "  Importing the query pool failed." >&2
        return 1
    fi

    [[ -f "$out_qza" ]] || { echo "  Import produced no artifact." >&2; return 1; }
    echo "  Query pool: $(count_asvs "$out_qza") distinct queries -> ${out_qza}"
}

exp_query_sweep() {
    echo "=== Experiment: classifier cost vs. query count, 1 to 10,000 ==="
    local csv="${OUTDIR}/query_sweep.csv"
    local pool="${QUERY_POOL:-${SCRATCH}/query_pool.qza}"

    # A SILVA extract left behind by build_classifier.sh --keep-intermediates is
    # a usable pool too, but real amplicons are the better query set, so it is
    # only picked up when the caller asks for it via QUERY_POOL.
    if [[ ! -f "$pool" ]]; then
        if ! build_query_pool "${QUERY_POOL_SAMPLE:-ERR3444628}" "$pool"; then
            echo "  Could not build a query pool; skipping this experiment." >&2
            return 1
        fi
    fi
    # Belt and braces: count_asvs on a missing artifact throws an unreadable
    # zipfile traceback, which is how the previous failure presented.
    if [[ ! -f "$pool" ]]; then
        echo "  Query pool artifact absent after build: ${pool}" >&2
        return 1
    fi

    local available; available=$(count_asvs "$pool")
    echo "  Query pool: ${pool}"
    echo "  Available:  ${available} sequences"

    local ks=() k
    for k in 1 10 50 100 218 500 1000 2000 5000 10000; do
        [[ $k -le $available ]] && ks+=("$k")
    done
    if [[ $available -lt 10000 ]]; then
        echo "  NOTE: pool holds ${available}, short of 10,000. The ladder stops" >&2
        echo "  at ${ks[-1]}; the crossover may lie beyond it." >&2
    fi

    # p=1 is the serial floor; p=8 is the mid-range setting used by asv-sweep, so
    # the two sweeps are directly comparable at the query counts they share.
    local mid=$(( CORES_PER_SOCKET < 8 ? CORES_PER_SOCKET : 8 ))
    local rep p
    for rep in $(seq 1 "$REPS"); do
        for p in 1 "$mid"; do
            for k in "${ks[@]}"; do
                local sub="${SCRATCH}/qseqs_${k}.qza"
                [[ -f "$sub" ]] || subset_rep_seqs_strided "$pool" "$k" "$sub"
                echo "  [rep $rep] classify p=$p queries=$k"
                drop_caches
                run_classify "$p" "$k" "$sub" "$csv" "classify_querysweep"
            done
        done
    done

    echo
    echo "  Wrote ${csv}."
    echo "  The crossover is the query count at which p=${mid} first beats p=1;"
    echo "  the linear term overtakes the fixed cost where the sweep stops being"
    echo "  flat. Both are read off this one CSV."
}

# ------------------------------------------------------------------------------
# Experiment: depth invariance  <-- the controlled comparison the paper turns on
#
# Two nasopharyngeal samples from the same study and body site, differing ~27x in
# sequencing depth but only ~7% in ASV count:
#
#   ERR3444605   66,993 reads -> 54 ASVs
#   ERR3444628  1,804,054 reads -> 58 ASVs
#
# Prediction: denoise time scales with reads, classify time does not. Both stages
# are measured at a single fixed core count so depth is the only variable.
# ------------------------------------------------------------------------------
exp_depth_pair() {
    echo "=== Experiment: depth invariance (shallow vs deep, same body site) ==="
    local csv="${OUTDIR}/depth_pair.csv"
    write_csv_header "$csv"
    local p=4   # fixed; this experiment varies depth, not parallelism
    local richness="${SCRATCH}/dataset_richness.tsv"

    local acc
    for acc in ERR3444605 ERR3444628; do
        local demux="${SCRATCH}/sample_${acc}.qza"
        local repseqs="${SCRATCH}/repseqs_${acc}.qza"
        if [[ ! -f "$demux" || ! -f "$repseqs" ]]; then
            echo "  Skipping ${acc}: not staged. Run scripts/prepare_datasets.sh first." >&2
            continue
        fi

        local reads=0 asvs
        [[ -f "$richness" ]] && reads=$(awk -F'\t' -v a="$acc" '$1==a{print $4}' "$richness")
        [[ -n "$reads" ]] || reads=0
        asvs=$(count_asvs "$repseqs")
        echo "  ${acc}: ${reads} reads, ${asvs} ASVs"

        local rep
        for rep in $(seq 1 "$REPS"); do
            # denoise: input size is the read count
            echo "  [rep $rep] denoise ${acc} (${reads} reads, p=$p)"
            drop_caches
            run_denoise "$p" "$reads" "$demux" "$csv" "denoise_depth_${acc}"

            # classify: input size is the ASV count, NOT the read count
            echo "  [rep $rep] classify ${acc} (${asvs} ASVs, p=$p)"
            drop_caches
            run_classify "$p" "$asvs" "$repseqs" "$csv" "classify_depth_${acc}"
        done
    done

    echo
    echo "  Compare denoise_depth_* against classify_depth_* in ${csv}."
    echo "  Denoise should differ ~27x between the two accessions; classify"
    echo "  should not differ appreciably."
}


# ------------------------------------------------------------------------------
# Experiment: within-sample depth series  <-- removes the patient confound
#
# The depth-pair experiment compares ERR3444605 (D-003-ANS) against ERR3444628
# (D-006-ANS). Those are DIFFERENT PATIENTS, so community composition differs
# alongside depth and the comparison is confounded. This experiment subsamples a
# single library to a range of depths, so the community is held exactly fixed and
# sequencing depth is genuinely the only variable.
#
# Uses ERR3444628 because at 1.8M read pairs it can be subsampled across more
# than an order of magnitude without leaving the range of real samples.
# ------------------------------------------------------------------------------
exp_depth_series() {
    echo "=== Experiment: within-sample depth series (community held fixed) ==="
    local csv="${OUTDIR}/depth_series.csv"
    write_csv_header "$csv"
    local acc="${DEPTH_SERIES_ACC:-ERR3444628}"
    local p=4
    local src="${SCRATCH}/${acc}"

    if [[ ! -f "${src}_1.fastq.gz" ]]; then
        echo "  Skipping: ${src}_1.fastq.gz not staged. Run prepare_datasets.sh." >&2
        return 0
    fi

    local total; total=$(count_fastq_reads "${src}_1.fastq")
    echo "  Source ${acc}: ${total} read pairs"

    local depths=()
    for frac in 27 14 7 4 2 1; do
        local n=$(( total / frac ))
        [[ $n -ge 5000 ]] && depths+=("$n")
    done

    local rep n
    for rep in $(seq 1 "$REPS"); do
        for n in "${depths[@]}"; do
            local prefix="${SCRATCH}/ds_${acc}_${n}"
            [[ -f "${prefix}_1.fastq.gz" ]] || subsample_fastq "$n" "$src" "$prefix"
            local demux="${SCRATCH}/ds_${acc}_${n}.qza"
            [[ -f "$demux" ]] || import_demux "$prefix" "ds-${n}" "$demux"

            echo "  [rep $rep] denoise ${acc} @ ${n} reads (p=$p)"
            drop_caches
            run_denoise "$p" "$n" "$demux" "$csv" "denoise_series"

            # classify the ASVs this depth produced -- the point is that the ASV
            # count, not the read count, is what the classifier sees.
            local rs="${SCRATCH}/out_denoise_series_${p}_${n}/rep-seqs.qza"
            if [[ -f "$rs" ]]; then
                local asvs; asvs=$(count_asvs "$rs")
                echo "  [rep $rep] classify ${acc} @ ${n} reads -> ${asvs} ASVs"
                drop_caches
                run_classify "$p" "$asvs" "$rs" "$csv" "classify_series_${n}reads"
            fi
        done
    done
    echo "  Depth and ASV count are both recorded, so cost can be regressed on each."
}

# ------------------------------------------------------------------------------
# Experiment: per-sample cost against measured richness
#
# Nine samples at comparable sequencing depth but spanning 5-70 ASVs. Regressing
# each stage's cost on that sample's own ASV count tests the richness claim
# across real communities rather than by subsetting one pooled set.
# ------------------------------------------------------------------------------
exp_per_sample() {
    echo "=== Experiment: per-sample cost vs richness ==="
    local csv="${OUTDIR}/per_sample.csv"
    write_csv_header "$csv"
    local richness="${SCRATCH}/dataset_richness.tsv"
    local p=4 rep acc

    for rep in $(seq 1 "$REPS"); do
        while read -r acc; do
            local demux="${SCRATCH}/sample_${acc}.qza"
            local rs="${SCRATCH}/repseqs_${acc}.qza"
            [[ -f "$demux" && -f "$rs" ]] || { echo "  skip ${acc}: not staged" >&2; continue; }
            local reads=0 asvs
            [[ -f "$richness" ]] && reads=$(awk -F'\t' -v a="$acc" '$1==a{print $4}' "$richness")
            [[ -n "$reads" ]] || reads=0
            asvs=$(count_asvs "$rs")
            echo "  [rep $rep] ${acc}: ${reads} reads, ${asvs} ASVs"
            drop_caches
            run_denoise  "$p" "$reads" "$demux" "$csv" "denoise_ps_${acc}"
            drop_caches
            run_classify "$p" "$asvs"  "$rs"    "$csv" "classify_ps_${acc}"
        done < <(grep -vE '^\s*(#|$)' "$REPO_ROOT/datasets.tsv" | cut -f1)
    done
}

# ------------------------------------------------------------------------------
# Experiment: output determinism across core counts
#
# The paper recommends changing --p-n-threads. That recommendation is only safe
# if the biological output does not change with it. DADA2 learns its error model
# from a subsample of reads, so identical output across thread counts is a
# property to be demonstrated, not assumed.
# ------------------------------------------------------------------------------
exp_determinism() {
    echo "=== Experiment: is DADA2 output identical across thread counts? ==="
    local report="${OUTDIR}/determinism.tsv"
    printf "threads\tasvs\ttotal_reads_in_table\trep_seqs_sha256\n" > "$report"

    local p out
    for p in 1 4 "$MAX_CORES"; do
        out="${SCRATCH}/det_${p}"; rm -rf "$out"; mkdir -p "$out"
        dada2_extra_outputs "$out"
        echo "  denoising at p=${p}..."
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$BASE_DEMUX" \
            --p-trunc-len-f "$TRUNC_LEN_F" --p-trunc-len-r "$TRUNC_LEN_R" \
            --p-n-threads "$p" \
            --o-table "${out}/table.qza" \
            --o-representative-sequences "${out}/rep-seqs.qza" \
            --o-denoising-stats "${out}/stats.qza" \
            ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"} >/dev/null

        # Hash the SEQUENCES, not the .qza: artifacts embed a fresh UUID and
        # timestamp on every run and would never compare equal.
        local info
        info=$(python3 - "${out}/rep-seqs.qza" "${out}/table.qza" <<'PY2'
import sys, zipfile, re, hashlib
z = zipfile.ZipFile(sys.argv[1])
fa = [n for n in z.namelist() if n.endswith("dna-sequences.fasta")][0]
seqs = sorted(s.replace("\n", "").strip().upper()
              for s in re.split(r">[^\n]*\n", z.read(fa).decode())[1:] if s.strip())
print(len(seqs), hashlib.sha256("".join(seqs).encode()).hexdigest()[:16])
PY2
)
        local asvs sha
        asvs=$(echo "$info" | cut -d' ' -f1); sha=$(echo "$info" | cut -d' ' -f2)
        printf "%s\t%s\t%s\t%s\n" "$p" "$asvs" "-" "$sha" >> "$report"
        echo "    p=${p}: ${asvs} ASVs, sequence-set sha256 ${sha}"
    done

    echo
    if [[ $(cut -f4 "$report" | tail -n +2 | sort -u | wc -l) -eq 1 ]]; then
        echo "  RESULT: identical ASV sequence sets at every thread count."
    else
        echo "  RESULT: ASV sequence sets DIFFER between thread counts." >&2
        echo "  This must be reported: the core-count recommendation changes the" >&2
        echo "  biological output, not just the runtime." >&2
    fi
    echo "  Written to ${report}"
}

# ------------------------------------------------------------------------------
# Experiment 5: NUMA locality
#
# Same core count, two placements: confined to socket 0 vs. spread across both.
# Amdahl cannot express this, which is precisely why it is worth measuring.
# ------------------------------------------------------------------------------
exp_numa() {
    echo "=== Experiment: NUMA locality ==="
    if ! command -v numactl &>/dev/null; then
        echo "  Skipping: numactl not installed (sudo apt-get install numactl)" >&2
        return 0
    fi
    if [[ ${SOCKETS:-1} -lt 2 ]]; then
        echo "  Skipping: single-socket system, nothing to compare" >&2
        return 0
    fi

    local csv="${OUTDIR}/numa_comparison.csv"
    write_csv_header "$csv"
    local p=$CORES_PER_SOCKET
    local reads; reads=$(count_fastq_reads "${BASE_PREFIX}_1.fastq")
    local out

    for rep in $(seq 1 "$REPS"); do
        # (a) pinned to socket 0, memory local
        out="${SCRATCH}/out_numa_local"; rm -rf "$out"; mkdir -p "$out"
        dada2_extra_outputs "$out"
        echo "  [rep $rep] denoise p=$p pinned to socket 0"
        drop_caches
        benchmark_stage "denoise_numa_local" "$p" "$reads" "$csv" \
            numactl --cpunodebind=0 --membind=0 \
            qiime dada2 denoise-paired \
                --i-demultiplexed-seqs "$BASE_DEMUX" \
                --p-trunc-len-f "$TRUNC_LEN_F" --p-trunc-len-r "$TRUNC_LEN_R" \
                --p-n-threads "$p" \
                --o-table "${out}/table.qza" \
                --o-representative-sequences "${out}/rep-seqs.qza" \
                --o-denoising-stats "${out}/stats.qza" \
                ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"}

        # (b) same core count, interleaved across both sockets
        out="${SCRATCH}/out_numa_split"; rm -rf "$out"; mkdir -p "$out"
        dada2_extra_outputs "$out"
        echo "  [rep $rep] denoise p=$p spread across sockets"
        drop_caches
        benchmark_stage "denoise_numa_split" "$p" "$reads" "$csv" \
            numactl --interleave=all \
            qiime dada2 denoise-paired \
                --i-demultiplexed-seqs "$BASE_DEMUX" \
                --p-trunc-len-f "$TRUNC_LEN_F" --p-trunc-len-r "$TRUNC_LEN_R" \
                --p-n-threads "$p" \
                --o-table "${out}/table.qza" \
                --o-representative-sequences "${out}/rep-seqs.qza" \
                --o-denoising-stats "${out}/stats.qza" \
                ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"}
    done
}

# ------------------------------------------------------------------------------
# Experiment 6: multi-sample throughput (sequential loop vs. SLURM array)
# ------------------------------------------------------------------------------
exp_array() {
    echo "=== Experiment: multi-sample throughput ==="
    local csv="${OUTDIR}/array_results.csv"
    write_csv_header "$csv"

    local samples=("${SCRATCH}"/sample_*.qza)
    if [[ ! -e "${samples[0]}" ]]; then
        echo "  Skipping: no per-sample demux artifacts found in ${SCRATCH}." >&2
        echo "  Run scripts/prepare_datasets.sh first to stage multiple samples." >&2
        return 0
    fi

    local n=${#samples[@]}
    local per_sample_cores=4
    echo "  ${n} samples, ${per_sample_cores} cores each"

    # (a) sequential loop — the baseline everyone actually uses
    local richness="${SCRATCH}/dataset_richness.tsv"
    local loop_start=$SECONDS
    for s in "${samples[@]}"; do
        local acc reads
        acc="$(basename "$s" .qza)"; acc="${acc#sample_}"
        # Recover the true read count so the CSV's `reads` column stays meaningful.
        reads=0
        [[ -f "$richness" ]] && reads=$(awk -F'\t' -v a="$acc" '$1==a{print $4}' "$richness")
        [[ -n "$reads" ]] || reads=0
        run_denoise "$per_sample_cores" "$reads" "$s" "$csv" "denoise_loop_${acc}"
    done
    local loop_elapsed=$(( SECONDS - loop_start ))
    echo "  Sequential loop: ${loop_elapsed}s"

    # (b) SLURM array, if a scheduler is present
    if command -v sbatch &>/dev/null; then
        echo "  Submitting SLURM array (see stages/slurm/submit_array.sh)"
        bash "$REPO_ROOT/stages/slurm/submit_array.sh" --samples "$n" --wait || true
    else
        echo "  sbatch not found — install SLURM via stages/slurm/setup_slurm.sh"
        echo "  to measure the array arm. Loop timing recorded regardless."
    fi

    echo "loop_total_sec,${loop_elapsed},${n}" >> "${OUTDIR}/array_loop_total.csv"
}

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
detect_topology
build_core_ladder
pin_blas_threads

BASE_PREFIX="${REPO_ROOT}/${SAMPLE_ID}"
BASE_DEMUX="${SCRATCH}/base_demux.qza"
BASE_REPSEQS="${SCRATCH}/base_rep-seqs.qza"
# Pooled ASV set for the classifier sweep. Rebuilt every run unless the caller
# points RICH_REPSEQS somewhere else: a cached file here previously masked the
# pooling step and made the sweep un-reproducible.
POOLED_REPSEQS_DEFAULT="${SCRATCH}/pooled_rep-seqs.qza"
RICH_REPSEQS="${RICH_REPSEQS:-$POOLED_REPSEQS_DEFAULT}"
CLASSIFIER_PATH="${REPO_ROOT}/${CLASSIFIER_FILE}"
CACHE_WARNED=0

cat <<EOF

========================================================
  Measured scaling benchmarks
========================================================
  Host topology:   ${SOCKETS} socket(s) x ${CORES_PER_SOCKET} cores
                   = ${PHYSICAL_CORES} physical, ${LOGICAL_CPUS} logical
  Core ladder:     ${CORE_LADDER[*]}
  Repetitions:     ${REPS}
  Experiments:     ${EXPERIMENTS[*]}
  Output CSVs:     ${OUTDIR}
  Scratch:         ${SCRATCH}
  QIIME 2 env:     ${QIIME_ENV}
  Classifier:      ${CLASSIFIER_FILE}
========================================================

EOF

if [[ $DRY_RUN -eq 1 ]]; then
    echo "Dry run — no measurements taken."
    exit 0
fi

mkdir -p "$OUTDIR" "$SCRATCH"
activate_qiime

if [[ ! -f "$CLASSIFIER_PATH" ]]; then
    echo "Error: classifier not found at $CLASSIFIER_PATH" >&2
    echo "Run: bash stages/04_classify.sh --download-only, or set CLASSIFIER_FILE." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Pre-flight: prove the classifier actually runs before committing to the sweep
#
# A scikit-learn version mismatch only surfaces when classify-sklearn is invoked,
# which for --all is well over an hour in. Worse, the failure still writes a CSV
# row carrying the seconds spent loading the model before raising. Classify one
# sequence up front and stop here if it does not work.
# ------------------------------------------------------------------------------
preflight_classifier() {
    local probe="${SCRATCH}/preflight_1asv.qza"
    local out="${SCRATCH}/preflight_out"
    [[ -f "$probe" ]] || subset_rep_seqs "$BASE_REPSEQS" 1 "$probe"
    rm -rf "$out"; mkdir -p "$out"

    echo "Pre-flight: verifying the classifier runs in this environment..."
    if qiime feature-classifier classify-sklearn \
            --i-classifier "$CLASSIFIER_PATH" \
            --i-reads "$probe" \
            --p-n-jobs 1 \
            --o-classification "${out}/tax.qza" >"${out}/log" 2>&1; then
        echo "  OK"
        return 0
    fi

    echo >&2
    echo "Error: the classifier failed on a single sequence. Aborting before" >&2
    echo "the sweep, since every classify measurement would be invalid." >&2
    echo >&2
    sed -n '1,12p' "${out}/log" >&2
    echo >&2
    if grep -qi 'scikit-learn version' "${out}/log"; then
        echo "This is a scikit-learn version mismatch. QIIME 2 stopped shipping" >&2
        echo "region-specific classifiers after 2024.2, so there is no prebuilt" >&2
        echo "515F-806R artifact for current scikit-learn. Rebuild it locally:" >&2
        echo >&2
        echo "  bash scripts/build_classifier.sh --threads ${MAX_CORES}" >&2
        echo >&2
        echo "That trains against your own environment, using the same Caporaso" >&2
        echo "515F/806R primers, so results stay comparable." >&2
    fi
    exit 1
}

# Only relevant when a classify-based experiment was requested.
NEEDS_CLASSIFIER=0
for _exp in "${EXPERIMENTS[@]}"; do
    case "$_exp" in
        strong-classify|asv-sweep|query-sweep|depth-pair|depth-series|per-sample) NEEDS_CLASSIFIER=1 ;;
    esac
done

# Stage the base artifacts once; every experiment reuses them.
if [[ ! -f "$BASE_DEMUX" ]]; then
    echo "Staging base demux artifact..."
    import_demux "$BASE_PREFIX" "$SAMPLE_NAME" "$BASE_DEMUX"
fi
if [[ ! -f "$BASE_REPSEQS" ]]; then
    if [[ -f "${REPO_ROOT}/rep-seqs.qza" ]]; then
        cp "${REPO_ROOT}/rep-seqs.qza" "$BASE_REPSEQS"
    else
        echo "Staging base rep-seqs (one denoise pass)..."
        tmp="${SCRATCH}/stage_repseqs"; rm -rf "$tmp"; mkdir -p "$tmp"
        dada2_extra_outputs "$tmp"
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$BASE_DEMUX" \
            --p-trunc-len-f "$TRUNC_LEN_F" --p-trunc-len-r "$TRUNC_LEN_R" \
            --p-n-threads "$MAX_CORES" \
            --o-table "${tmp}/table.qza" \
            --o-representative-sequences "$BASE_REPSEQS" \
            --o-denoising-stats "${tmp}/stats.qza" \
            ${DADA2_EXTRA_OUTPUTS[@]+"${DADA2_EXTRA_OUTPUTS[@]}"} >/dev/null
    fi
fi
[[ "${NEEDS_CLASSIFIER:-0}" -eq 1 ]] && preflight_classifier

# ------------------------------------------------------------------------------
# Pooled ASV set for the classifier sweep
#
# No single sample here carries enough ASVs to trace a cost curve (5-70). The
# union across the cohort does, and it is also what a real multi-sample study
# actually classifies: denoise every sample, then classify the set of distinct
# representative sequences once. Sequences are deduplicated exactly, since the
# same ASV recurs across samples from the same body site.
# ------------------------------------------------------------------------------
build_pooled_repseqs() {
    local out_qza="$1"
    local fasta="${SCRATCH}/pooled_rep_seqs.fasta"
    local accs=()
    while read -r acc; do accs+=("${SCRATCH}/repseqs_${acc}.qza"); done \
        < <(grep -vE '^\s*(#|$)' "$REPO_ROOT/datasets.tsv" | cut -f1)

    python3 - "$fasta" "${accs[@]}" <<'PY'
import sys, zipfile, re, os, hashlib
out, srcs = sys.argv[1], sys.argv[2:]
seen, kept, missing = {}, 0, []
for src in srcs:
    if not os.path.isfile(src):
        missing.append(os.path.basename(src)); continue
    z = zipfile.ZipFile(src)
    fa = [n for n in z.namelist() if n.endswith("dna-sequences.fasta")][0]
    for rec in re.split(r">[^\n]*\n", z.read(fa).decode())[1:]:
        s = rec.replace("\n", "").strip().upper()
        if s and s not in seen:
            seen[s] = hashlib.md5(s.encode()).hexdigest()
            kept += 1
with open(out, "w") as fh:
    for s, sid in seen.items():
        fh.write(f">{sid}\n{s}\n")
if missing:
    print(f"  (skipped {len(missing)} uncached: {', '.join(missing[:3])}"
          f"{'...' if len(missing) > 3 else ''})", file=sys.stderr)
print(kept)
PY
}

if [[ "$RICH_REPSEQS" == "$POOLED_REPSEQS_DEFAULT" || ! -f "$RICH_REPSEQS" ]]; then
    rm -f "$RICH_REPSEQS"
    echo "Pooling ASVs across the cohort for the classifier sweep..."
    n_pooled=$(build_pooled_repseqs "$RICH_REPSEQS")
    rm -f "$RICH_REPSEQS"
    qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path "${SCRATCH}/pooled_rep_seqs.fasta" \
        --output-path "$RICH_REPSEQS" >/dev/null
    echo "  Pooled ASV set: ${n_pooled} distinct sequences"
fi

# ------------------------------------------------------------------------------
# Dispatch
# ------------------------------------------------------------------------------
START_TS=$SECONDS
for exp in "${EXPERIMENTS[@]}"; do
    case "$exp" in
        strong-denoise)   exp_strong_denoise ;;
        weak-denoise)     exp_weak_denoise ;;
        strong-classify)  exp_strong_classify ;;
        asv-sweep)        exp_asv_sweep ;;
        query-sweep)      exp_query_sweep ;;
        depth-pair)       exp_depth_pair ;;
        depth-series)     exp_depth_series ;;
        per-sample)       exp_per_sample ;;
        determinism)      exp_determinism ;;
        numa)             exp_numa ;;
        array)            exp_array ;;
        *) echo "Unknown experiment: $exp" >&2; exit 1 ;;
    esac
done
TOTAL=$(( SECONDS - START_TS ))

echo
echo "========================================================"
printf "  Done in %02d:%02d:%02d\n" $((TOTAL/3600)) $((TOTAL%3600/60)) $((TOTAL%60))
echo "  CSVs written to ${OUTDIR}"
if [[ ${CACHE_WARNED} -eq 1 ]]; then
    echo
    echo "  NOTE: page cache was NOT dropped between repetitions (no passwordless"
    echo "        sudo). I/O figures reflect warm-cache behaviour. Re-run under"
    echo "        sudo for cold-cache numbers, and say which you used in Methods."
fi
echo "========================================================"
echo
echo "Next:"
echo "  python3 lib/predictor.py ${OUTDIR}/strong_scaling.csv --save-models ${OUTDIR}/fitted_models.json"
echo "  python3 paper/analysis/generate_figures.py"
echo "  python3 paper/analysis/analyse_query_sweep.py   # if query-sweep was run"
