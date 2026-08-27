#!/bin/bash
# lib/cache.sh — Content-addressable caching for pipeline stages
# Computes SHA256 hashes of input files to determine if a stage's outputs
# are still valid. Avoids unnecessary recomputation when inputs haven't changed.
#
# Usage:
#   source lib/cache.sh
#   key=$(compute_cache_key input1.qza input2.qza)
#   if cache_hit ".cache" "denoise" "$key" "table.qza"; then
#       echo "Using cached output"
#   else
#       run_stage ...
#       cache_store ".cache" "denoise" "$key" "table.qza"
#   fi

set -euo pipefail

CACHE_DIR="${CACHE_DIR:-.cache}"

# Compute a cache key from one or more input files
# Returns SHA256 hash of concatenated file hashes + filenames
compute_cache_key() {
    local hash_input=""
    for f in "$@"; do
        if [[ -f "$f" ]]; then
            local file_hash
            file_hash=$(sha256sum "$f" | cut -d' ' -f1)
            hash_input+="${f}:${file_hash};"
        else
            echo "Warning: cache key input not found: $f" >&2
            hash_input+="${f}:MISSING;"
        fi
    done
    echo -n "$hash_input" | sha256sum | cut -d' ' -f1
}

# Check if a cached output exists for the given key
# Arguments: $1=cache_dir, $2=stage_name, $3=cache_key, $4=output_file
# Returns 0 if cache hit, 1 if miss
cache_hit() {
    local cache_dir="$1"
    local stage="$2"
    local key="$3"
    local output_file="$4"

    local cached_path="${cache_dir}/${stage}/${key}/$(basename "$output_file")"
    if [[ -f "$cached_path" ]]; then
        echo "[cache] HIT: ${stage}/${key:0:12}... → $(basename "$output_file")"
        return 0
    fi
    echo "[cache] MISS: ${stage}/${key:0:12}..."
    return 1
}

# Store an output file in the cache
# Arguments: $1=cache_dir, $2=stage_name, $3=cache_key, $4=output_file
cache_store() {
    local cache_dir="$1"
    local stage="$2"
    local key="$3"
    local output_file="$4"

    local target_dir="${cache_dir}/${stage}/${key}"
    mkdir -p "$target_dir"
    cp "$output_file" "${target_dir}/$(basename "$output_file")"
    echo "[cache] STORED: ${stage}/${key:0:12}... ← $(basename "$output_file")"
}

# Restore a cached output file to its expected location
# Arguments: $1=cache_dir, $2=stage_name, $3=cache_key, $4=output_file
cache_restore() {
    local cache_dir="$1"
    local stage="$2"
    local key="$3"
    local output_file="$4"

    local cached_path="${cache_dir}/${stage}/${key}/$(basename "$output_file")"
    if [[ -f "$cached_path" ]]; then
        cp "$cached_path" "$output_file"
        echo "[cache] RESTORED: $(basename "$output_file") ← ${stage}/${key:0:12}..."
        return 0
    fi
    return 1
}

# Evict all cache entries for a stage
cache_evict_stage() {
    local cache_dir="$1"
    local stage="$2"
    if [[ -d "${cache_dir}/${stage}" ]]; then
        rm -rf "${cache_dir:?}/${stage}"
        echo "[cache] EVICTED: all entries for stage '${stage}'"
    fi
}

# Report cache size in human-readable format
cache_stats() {
    local cache_dir="${1:-$CACHE_DIR}"
    if [[ -d "$cache_dir" ]]; then
        echo "[cache] Cache directory: $cache_dir"
        du -sh "$cache_dir" 2>/dev/null || echo "[cache] Empty"
        echo "[cache] Entries by stage:"
        for stage_dir in "$cache_dir"/*/; do
            [[ -d "$stage_dir" ]] || continue
            local stage_name
            stage_name=$(basename "$stage_dir")
            local count
            count=$(find "$stage_dir" -mindepth 1 -maxdepth 1 -type d | wc -l)
            echo "  ${stage_name}: ${count} cached keys"
        done
    else
        echo "[cache] No cache directory found"
    fi
}
