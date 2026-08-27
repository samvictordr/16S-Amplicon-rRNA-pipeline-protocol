#!/bin/bash
# lib/benchmark.sh — Stage-level instrumentation for HPC scaling experiments
# Wraps pipeline stages with /usr/bin/time -v to capture wall-clock, CPU time,
# peak RSS, context switches, and I/O metrics. Outputs structured CSV rows.
#
# Usage (source from any stage or experiment script):
#   source lib/benchmark.sh
#   benchmark_stage "denoise" 4 55000 "benchmarks/results.csv" \
#       qiime dada2 denoise-paired ...

set -euo pipefail

BENCHMARK_DIR="${BENCHMARK_DIR:-benchmarks}"
BENCHMARK_TIME_BIN="${BENCHMARK_TIME_BIN:-/usr/bin/time}"

# Ensure benchmark output directory exists
ensure_benchmark_dir() {
    mkdir -p "$BENCHMARK_DIR"
}

# Write CSV header if the file doesn't exist or is empty
write_csv_header() {
    local csv_file="$1"
    if [[ ! -f "$csv_file" ]] || [[ ! -s "$csv_file" ]]; then
        echo "stage,threads,reads,wall_sec,user_sec,sys_sec,peak_rss_kb,ctx_voluntary,ctx_involuntary,fs_inputs,fs_outputs,timestamp,exit_code" \
            > "$csv_file"
    fi
}

# Parse GNU time -v output into CSV fields
# Arguments: $1 = path to time output file
# Sets global variables: _BM_WALL, _BM_USER, _BM_SYS, _BM_RSS, _BM_CTX_VOL,
#                        _BM_CTX_INV, _BM_FS_IN, _BM_FS_OUT, _BM_EXIT
parse_time_output() {
    local time_file="$1"

    # Wall clock (h:mm:ss or m:ss.ss format)
    local wall_raw
    wall_raw=$(grep "Elapsed (wall clock)" "$time_file" | sed 's/.*: //')
    _BM_WALL=$(parse_wall_time "$wall_raw")

    # User time (seconds)
    _BM_USER=$(grep "User time" "$time_file" | sed 's/.*: //')

    # System time (seconds)
    _BM_SYS=$(grep "System time" "$time_file" | sed 's/.*: //')

    # Peak resident set size (KB)
    _BM_RSS=$(grep "Maximum resident" "$time_file" | sed 's/.*: //')

    # Context switches
    _BM_CTX_VOL=$(grep "Voluntary context" "$time_file" | sed 's/.*: //')
    _BM_CTX_INV=$(grep "Involuntary context" "$time_file" | sed 's/.*: //')

    # Filesystem I/O
    _BM_FS_IN=$(grep "File system inputs" "$time_file" | sed 's/.*: //')
    _BM_FS_OUT=$(grep "File system outputs" "$time_file" | sed 's/.*: //')

    # Exit status
    _BM_EXIT=$(grep "Exit status" "$time_file" | sed 's/.*: //')
}

# Convert wall clock time string to seconds
# Handles: "h:mm:ss", "m:ss.ss", "ss.ss"
parse_wall_time() {
    local raw="$1"
    local hours=0 mins=0 secs=0

    # Count colons to determine format
    local colons
    colons=$(echo "$raw" | tr -cd ':' | wc -c)

    if [[ "$colons" -eq 2 ]]; then
        # h:mm:ss or h:mm:ss.ss
        IFS=: read -r hours mins secs <<< "$raw"
    elif [[ "$colons" -eq 1 ]]; then
        # m:ss.ss
        IFS=: read -r mins secs <<< "$raw"
    else
        secs="$raw"
    fi

    echo "$(echo "$hours * 3600 + $mins * 60 + $secs" | bc -l)"
}

# Main benchmarking wrapper
# Arguments:
#   $1 = stage name (e.g., "denoise", "classify")
#   $2 = thread count
#   $3 = read count
#   $4 = CSV output file path
#   $5... = command to benchmark
benchmark_stage() {
    local stage_name="$1"
    local threads="$2"
    local reads="$3"
    local csv_file="$4"
    shift 4

    ensure_benchmark_dir
    write_csv_header "$csv_file"

    local time_output
    time_output=$(mktemp "${BENCHMARK_DIR}/time_${stage_name}_XXXXXX.txt")

    local exit_code=0
    # Run the command under /usr/bin/time -v, capturing timing to file
    "$BENCHMARK_TIME_BIN" -v -o "$time_output" "$@" || exit_code=$?

    # Parse the timing output
    parse_time_output "$time_output"

    # Append CSV row (with flock to avoid race conditions in concurrent SLURM jobs)
    local timestamp
    timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    local csv_row="${stage_name},${threads},${reads},${_BM_WALL},${_BM_USER},${_BM_SYS},${_BM_RSS},${_BM_CTX_VOL},${_BM_CTX_INV},${_BM_FS_IN},${_BM_FS_OUT},${timestamp},${_BM_EXIT}"
    (
        flock -x 200
        echo "$csv_row" >> "$csv_file"
    ) 200>"${csv_file}.lock"

    # Clean up temp file
    rm -f "$time_output"

    echo "[benchmark] ${stage_name} (${threads} threads, ${reads} reads): ${_BM_WALL}s wall, ${_BM_RSS}KB RSS"

    return "$exit_code"
}

# Lightweight timer for stages where /usr/bin/time is unavailable
# Uses bash SECONDS variable — less precise but always available
benchmark_stage_lite() {
    local stage_name="$1"
    local threads="$2"
    local reads="$3"
    local csv_file="$4"
    shift 4

    ensure_benchmark_dir
    write_csv_header "$csv_file"

    local start_time=$SECONDS
    local exit_code=0
    "$@" || exit_code=$?
    local elapsed=$(( SECONDS - start_time ))

    local timestamp
    timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    # Lite mode: user/sys/rss/ctx/io fields set to 0
    echo "${stage_name},${threads},${reads},${elapsed},0,0,0,0,0,0,0,${timestamp},${exit_code}" \
        >> "$csv_file"

    echo "[benchmark-lite] ${stage_name}: ${elapsed}s wall"
    return "$exit_code"
}
