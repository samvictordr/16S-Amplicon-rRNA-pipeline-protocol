#!/bin/bash
# lib/common.sh — Sourced by all stages and edna.sh
# Provides: defaults, load_config(), parse_args(), derive_ebi_url(), activate_qiime(),
#           check_dependencies(), validate_params()

set -euo pipefail

# ------------------------------------------------------------------------------
# DEFAULTS (match original hardcoded values exactly — backward compatible)
# ------------------------------------------------------------------------------
SAMPLE_ID="${SAMPLE_ID:-ERR3444605}"
SAMPLE_NAME="${SAMPLE_NAME:-nasopharynx-D003}"
QIIME_ENV="${QIIME_ENV:-qiime2-amplicon-2025.7}"
TRUNC_LEN_F="${TRUNC_LEN_F:-240}"
TRUNC_LEN_R="${TRUNC_LEN_R:-240}"
# 515F-806R (V4) classifier — must match the amplicon region. See pipeline.conf.
CLASSIFIER_URL="${CLASSIFIER_URL:-https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza}"
CLASSIFIER_FILE="${CLASSIFIER_FILE:-silva-138-99-515-806-nb-classifier.qza}"
OUTPUT_DIR="${OUTPUT_DIR:-.}"
MANUAL_URL_R1="${MANUAL_URL_R1:-}"
MANUAL_URL_R2="${MANUAL_URL_R2:-}"
THREADS="${THREADS:-1}"
START_AT="${START_AT:-1}"
BENCHMARK_MODE="${BENCHMARK_MODE:-off}"
BENCHMARK_DIR="${BENCHMARK_DIR:-benchmarks}"
SLURM_PARTITION="${SLURM_PARTITION:-main}"

# ------------------------------------------------------------------------------
# load_config(file) — Parse KEY="value" pairs from a config file
# Uses eval-free parsing. Handles values containing '=' (e.g. URLs).
# ------------------------------------------------------------------------------
load_config() {
    local cfg_file="$1"
    if [ ! -f "$cfg_file" ]; then
        echo "Warning: config file not found: $cfg_file" >&2
        return 1
    fi
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Skip blank lines and comments
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        # Split on first '=' only (preserves '=' in values like URLs)
        local key="${line%%=*}"
        local val="${line#*=}"
        # Trim leading/trailing whitespace from key
        key="$(echo "$key" | xargs)"
        # Strip surrounding quotes from value
        val="${val%\"}"
        val="${val#\"}"
        val="${val%\'}"
        val="${val#\'}"
        # Only set recognized config keys (security: no arbitrary variable injection)
        case "$key" in
            SAMPLE_ID|SAMPLE_NAME|QIIME_ENV|TRUNC_LEN_F|TRUNC_LEN_R|\
            CLASSIFIER_URL|CLASSIFIER_FILE|OUTPUT_DIR|MANUAL_URL_R1|\
            MANUAL_URL_R2|THREADS|START_AT|\
            BENCHMARK_MODE|BENCHMARK_DIR|SLURM_PARTITION)
                printf -v "$key" '%s' "$val"
                ;;
            *)
                echo "Warning: ignoring unknown config key: $key" >&2
                ;;
        esac
    done < "$cfg_file"
}

# ------------------------------------------------------------------------------
# parse_args("$@") — Override config values via CLI flags
# ------------------------------------------------------------------------------
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --config)           load_config "$2";          shift 2 ;;
            --sample-id)        SAMPLE_ID="$2";            shift 2 ;;
            --sample-name)      SAMPLE_NAME="$2";          shift 2 ;;
            --qiime-env)        QIIME_ENV="$2";            shift 2 ;;
            --trunc-f)          TRUNC_LEN_F="$2";          shift 2 ;;
            --trunc-r)          TRUNC_LEN_R="$2";          shift 2 ;;
            --classifier-url)   CLASSIFIER_URL="$2";       shift 2 ;;
            --classifier-file)  CLASSIFIER_FILE="$2";      shift 2 ;;
            --output-dir)       OUTPUT_DIR="$2";           shift 2 ;;
            --url-r1)           MANUAL_URL_R1="$2";        shift 2 ;;
            --url-r2)           MANUAL_URL_R2="$2";        shift 2 ;;
            --threads)          THREADS="$2";              shift 2 ;;
            --start-at)         START_AT="$2";             shift 2 ;;
            --benchmark)        BENCHMARK_MODE="$2";       shift 2 ;;
            --benchmark-dir)    BENCHMARK_DIR="$2";        shift 2 ;;
            --help|-h)          return 0 ;;
            *) echo "Unknown argument: $1" >&2; exit 1 ;;
        esac
    done
}

# ------------------------------------------------------------------------------
# show_help(script_path) — Print help text from a script's header comments
# Works regardless of where --help appears in the argument list.
# ------------------------------------------------------------------------------
show_help() {
    local script="$1"
    shift
    local arg
    for arg in "$@"; do
        if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then
            sed -n '/^# Usage:/,/^[^#]/{ /^#/s/^# \?//p }' "$script"
            exit 0
        fi
    done
}

# ------------------------------------------------------------------------------
# derive_ebi_url(sample_id, read_num) — Auto-construct EBI HTTPS URL
#
# Uses HTTPS (not FTP) to avoid firewall issues common in institutional networks.
# EBI URL pattern:
#   https://ftp.sra.ebi.ac.uk/vol1/fastq/{first6}/{subdir}/{id}/{id}_{N}.fastq.gz
# Rules:
#   first6  = first 6 chars of accession
#   subdir  = last digit zero-padded to 3 chars, only if len(accession) >= 10
# ------------------------------------------------------------------------------
derive_ebi_url() {
    local sample="$1"
    local read_num="$2"
    local prefix="${sample:0:6}"
    local len="${#sample}"
    local base="https://ftp.sra.ebi.ac.uk/vol1/fastq"

    if [ "$len" -le 9 ]; then
        echo "${base}/${prefix}/${sample}/${sample}_${read_num}.fastq.gz"
    else
        local subdir
        subdir=$(printf '%03d' "${sample: -1}")
        echo "${base}/${prefix}/${subdir}/${sample}/${sample}_${read_num}.fastq.gz"
    fi
}

# ------------------------------------------------------------------------------
# activate_qiime() — Source conda and activate the QIIME 2 environment
# Handles both conda and mamba, and checks if already active.
# Suppresses the pkg_resources deprecation warning from QIIME 2 internals
# (https://setuptools.pypa.io/en/latest/pkg_resources.html).
# ------------------------------------------------------------------------------
activate_qiime() {
    # Suppress pkg_resources deprecation warning from QIIME 2's q2_demux
    if [[ "${PYTHONWARNINGS:-}" != *"pkg_resources is deprecated"* ]]; then
        export PYTHONWARNINGS="${PYTHONWARNINGS:+${PYTHONWARNINGS},}ignore:pkg_resources is deprecated:UserWarning"
    fi

    # Skip if already in the target environment
    if [[ "${CONDA_DEFAULT_ENV:-}" == "$QIIME_ENV" ]]; then
        return 0
    fi

    # QIIME_ENV="current" (or empty) means: use whatever environment is already
    # active. Useful when the local env name differs from the one in
    # pipeline.conf, which is the common case across machines.
    if [[ -z "$QIIME_ENV" || "$QIIME_ENV" == "current" ]]; then
        if command -v qiime &>/dev/null; then
            echo "  Using active environment: ${CONDA_DEFAULT_ENV:-<none>}"
            return 0
        fi
        echo "Error: QIIME_ENV is 'current' but no 'qiime' command is on PATH." >&2
        echo "Activate your QIIME 2 environment first, or set QIIME_ENV in pipeline.conf." >&2
        exit 1
    fi

    local conda_base
    if command -v conda &>/dev/null; then
        conda_base="$(conda info --base 2>/dev/null)"
    elif command -v mamba &>/dev/null; then
        conda_base="$(mamba info --base 2>/dev/null)"
    else
        echo "Error: neither conda nor mamba found. Install conda first:" >&2
        echo "  https://docs.conda.io/en/latest/miniconda.html" >&2
        exit 1
    fi

    # shellcheck disable=SC1091
    source "${conda_base}/etc/profile.d/conda.sh"
    conda activate "$QIIME_ENV" 2>/dev/null || {
        echo "Error: could not activate environment '$QIIME_ENV'" >&2
        if [[ -n "${CONDA_DEFAULT_ENV:-}" ]] && command -v qiime &>/dev/null; then
            echo "" >&2
            echo "You are currently in '${CONDA_DEFAULT_ENV}', which does have qiime." >&2
            echo "Either set QIIME_ENV=\"${CONDA_DEFAULT_ENV}\" in pipeline.conf," >&2
            echo "or set QIIME_ENV=\"current\" to always use the active environment." >&2
        fi
        echo "" >&2
        echo "Available environments:" >&2
        conda env list >&2
        exit 1
    }
}

# ------------------------------------------------------------------------------
# q2_extra_output(plugin, action, flag, path) — Version-dependent QIIME 2 outputs
#
# QIIME 2 2026.x turned several previously non-existent outputs into REQUIRED
# ones (dada2 denoise-paired gained --o-base-transition-stats, feature-classifier
# extract-reads gained --o-read-extraction-stats). Older releases reject those
# same flags as unknown, so neither hardcoding nor omitting them works across
# versions. Probe the action's help text and emit the flag only where supported.
#
# Results are cached per plugin/action/flag: `qiime ... --help` costs ~15 s and
# the benchmark sweep invokes these actions on the order of a hundred times.
#
# Usage:
#   q2_extra_output dada2 denoise-paired --o-base-transition-stats "$dir/bts.qza"
#   qiime dada2 denoise-paired ... ${Q2_EXTRA_OUTPUTS[@]+"${Q2_EXTRA_OUTPUTS[@]}"}
# ------------------------------------------------------------------------------
q2_extra_output() {
    local plugin="$1" action="$2" flag="$3" path="$4"

    local key="_Q2_FLAG_${plugin}_${action}_${flag}"
    key="${key//[^a-zA-Z0-9_]/_}"

    if [[ -z "${!key:-}" ]]; then
        if qiime "$plugin" "$action" --help 2>/dev/null | grep -q -- "$flag"; then
            printf -v "$key" '%s' yes
        else
            printf -v "$key" '%s' no
        fi
    fi

    Q2_EXTRA_OUTPUTS=()
    [[ "${!key}" == "yes" ]] && Q2_EXTRA_OUTPUTS=("$flag" "$path")
    return 0
}

# ------------------------------------------------------------------------------
# dada2_extra_outputs(out_dir) — Convenience wrapper for denoise-paired
# Sets DADA2_EXTRA_OUTPUTS (empty on older QIIME 2).
# ------------------------------------------------------------------------------
dada2_extra_outputs() {
    local out_dir="${1:-.}"
    q2_extra_output dada2 denoise-paired \
        --o-base-transition-stats "${out_dir}/base-transition-stats.qza"
    DADA2_EXTRA_OUTPUTS=(${Q2_EXTRA_OUTPUTS[@]+"${Q2_EXTRA_OUTPUTS[@]}"})
}

# ------------------------------------------------------------------------------
# check_dependencies(tool...) — Verify required commands are available
# Usage: check_dependencies wget python qiime biom
# ------------------------------------------------------------------------------
check_dependencies() {
    local missing=()
    for cmd in "$@"; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "Error: missing required tools: ${missing[*]}" >&2
        echo "Install them before running the pipeline." >&2
        exit 1
    fi
}

# ------------------------------------------------------------------------------
# detect_pkg_manager() — Detect the system package manager
# Sets PKG_MANAGER and PKG_INSTALL globals. Returns 1 if none found.
# ------------------------------------------------------------------------------
detect_pkg_manager() {
    if command -v apt-get &>/dev/null; then
        PKG_MANAGER="apt-get"
        PKG_INSTALL="sudo apt-get install -y"
    elif command -v dnf &>/dev/null; then
        PKG_MANAGER="dnf"
        PKG_INSTALL="sudo dnf install -y"
    elif command -v yum &>/dev/null; then
        PKG_MANAGER="yum"
        PKG_INSTALL="sudo yum install -y"
    elif command -v pacman &>/dev/null; then
        PKG_MANAGER="pacman"
        PKG_INSTALL="sudo pacman -S --noconfirm"
    elif command -v brew &>/dev/null; then
        PKG_MANAGER="brew"
        PKG_INSTALL="brew install"
    else
        return 1
    fi
}

# ------------------------------------------------------------------------------
# resolve_pkg_name(tool) — Map a command name to its package name
# Some tools have different package names depending on the manager.
# ------------------------------------------------------------------------------
resolve_pkg_name() {
    local tool="$1"
    case "$tool" in
        python)
            case "$PKG_MANAGER" in
                apt-get)  echo "python3" ;;
                pacman)   echo "python" ;;
                brew)     echo "python3" ;;
                *)        echo "python3" ;;
            esac
            ;;
        *)
            echo "$tool"
            ;;
    esac
}

# ------------------------------------------------------------------------------
# install_dependencies(tool...) — Check for missing tools and install them
#
# Detects the system package manager, maps tool names to package names,
# and installs anything that is missing. Falls back to an error if no
# package manager is found or if installation fails.
# Usage: install_dependencies wget python
# ------------------------------------------------------------------------------
install_dependencies() {
    local missing=()
    for cmd in "$@"; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done

    [[ ${#missing[@]} -eq 0 ]] && return 0

    echo "Missing system tools: ${missing[*]}"

    if ! detect_pkg_manager; then
        echo "Error: no supported package manager found (tried apt-get, dnf, yum, pacman, brew)." >&2
        echo "Please install manually: ${missing[*]}" >&2
        exit 1
    fi

    echo "Detected package manager: $PKG_MANAGER"

    # Update package index for apt-get (others don't require this)
    if [[ "$PKG_MANAGER" == "apt-get" ]]; then
        echo "Updating package index..."
        sudo apt-get update -qq
    fi

    for cmd in "${missing[@]}"; do
        local pkg
        pkg="$(resolve_pkg_name "$cmd")"
        echo "Installing $pkg..."
        if ! $PKG_INSTALL "$pkg"; then
            echo "Error: failed to install '$pkg'. Please install it manually." >&2
            exit 1
        fi
    done

    # Verify installation succeeded
    for cmd in "${missing[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            echo "Error: '$cmd' still not found after installation. Please install it manually." >&2
            exit 1
        fi
    done

    echo "All system dependencies installed successfully."
}

# ------------------------------------------------------------------------------
# install_python_packages(package...) — Check for missing Python packages
# and install them via pip.
# Usage: install_python_packages pandas plotly
# ------------------------------------------------------------------------------
install_python_packages() {
    local missing=()
    for pkg in "$@"; do
        if ! python -c "import $pkg" 2>/dev/null; then
            missing+=("$pkg")
        fi
    done

    [[ ${#missing[@]} -eq 0 ]] && return 0

    echo "Missing Python packages: ${missing[*]}"

    local pip_cmd=""
    if command -v pip3 &>/dev/null; then
        pip_cmd="pip3"
    elif command -v pip &>/dev/null; then
        pip_cmd="pip"
    elif python -m pip --version &>/dev/null 2>&1; then
        pip_cmd="python -m pip"
    else
        echo "Error: pip not found. Install pip first, then run:" >&2
        echo "  pip install ${missing[*]}" >&2
        exit 1
    fi

    echo "Installing Python packages via $pip_cmd..."
    if ! $pip_cmd install "${missing[@]}"; then
        echo "Error: failed to install Python packages. Try manually:" >&2
        echo "  $pip_cmd install ${missing[*]}" >&2
        exit 1
    fi

    # Verify installation
    for pkg in "${missing[@]}"; do
        if ! python -c "import $pkg" 2>/dev/null; then
            echo "Error: Python package '$pkg' still not importable after installation." >&2
            exit 1
        fi
    done

    echo "All Python packages installed successfully."
}

# ------------------------------------------------------------------------------
# validate_params() — Validate key parameters before pipeline execution
# ------------------------------------------------------------------------------
validate_params() {
    # Validate SAMPLE_ID format (ENA/SRA accessions: letters + digits)
    if [[ ! "$SAMPLE_ID" =~ ^[A-Z]{3}[0-9]+$ ]]; then
        echo "Error: SAMPLE_ID '$SAMPLE_ID' doesn't look like a valid ENA accession (e.g., ERR3444605)" >&2
        exit 1
    fi

    # Validate truncation lengths are positive integers
    if [[ ! "$TRUNC_LEN_F" =~ ^[0-9]+$ ]] || [[ "$TRUNC_LEN_F" -lt 1 ]]; then
        echo "Error: --trunc-f must be a positive integer, got '$TRUNC_LEN_F'" >&2
        exit 1
    fi
    if [[ ! "$TRUNC_LEN_R" =~ ^[0-9]+$ ]] || [[ "$TRUNC_LEN_R" -lt 1 ]]; then
        echo "Error: --trunc-r must be a positive integer, got '$TRUNC_LEN_R'" >&2
        exit 1
    fi

    # Validate thread count (0 = all available cores in DADA2)
    if [[ ! "$THREADS" =~ ^[0-9]+$ ]]; then
        echo "Error: --threads must be a non-negative integer, got '$THREADS'" >&2
        exit 1
    fi

    # Validate START_AT
    if [[ ! "$START_AT" =~ ^[1-6]$ ]]; then
        echo "Error: --start-at must be between 1 and 6, got '$START_AT'" >&2
        exit 1
    fi
}
