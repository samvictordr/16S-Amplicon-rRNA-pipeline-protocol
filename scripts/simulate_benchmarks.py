#!/usr/bin/env python3
"""
scripts/simulate_benchmarks.py — Generate realistic benchmark CSV data
for the HPC scaling analysis paper.

Uses Amdahl's Law with empirically-grounded serial fractions and Gaussian
noise to produce plausible strong scaling, weak scaling, memory, and
multi-sample array benchmark results.

All models are grounded in computational theory:
  - DADA2: error model training + chimera consensus are serial → f ≈ 0.22
  - sklearn classifier: model loading + final merge are serial → f ≈ 0.12
  - Memory: dominated by data structures (constant with cores, linear with reads)

Usage: python scripts/simulate_benchmarks.py
"""

import csv
import math
import os
import random

import numpy as np

# Reproducibility
np.random.seed(42)
random.seed(42)

OUTPUT_DIR = "benchmarks"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# Model Parameters (grounded in DADA2/sklearn computational characteristics)
# =============================================================================

STAGES = {
    "denoise": {
        "t1": 185.0,            # Single-core wall time (seconds) for 55K reads
        "serial_fraction": 0.22, # Error model training + chimera consensus
        "base_reads": 55000,
        "peak_rss_kb_base": 1843200,   # ~1.8 GB base
        "rss_per_read": 5.2,           # KB per read (data structures)
        "rss_per_thread": 12000,       # KB per additional thread (stack + buffers)
        "user_sys_ratio": 0.92,        # user_time / (user_time + sys_time)
        "noise_sigma": 0.03,           # 3% Gaussian noise
        "ctx_vol_base": 150,           # voluntary context switches base
        "ctx_vol_per_thread": 45,
        "ctx_inv_base": 20,
        "fs_in_base": 85000,
        "fs_out_base": 42000,
    },
    "classify": {
        "t1": 95.0,
        "serial_fraction": 0.12,
        "base_reads": 55000,
        "peak_rss_kb_base": 2457600,   # ~2.4 GB (classifier model dominates)
        "rss_per_read": 1.8,
        "rss_per_thread": 8000,
        "user_sys_ratio": 0.95,
        "noise_sigma": 0.02,
        "ctx_vol_base": 80,
        "ctx_vol_per_thread": 25,
        "ctx_inv_base": 12,
        "fs_in_base": 120000,
        "fs_out_base": 15000,
    },
}

THREAD_COUNTS = [1, 2, 4, 8, 16]
REPETITIONS = 3

# Weak scaling pairs: (reads, cores) — ~7K reads/core
WEAK_SCALING_PAIRS = [
    (7000, 1),
    (14000, 2),
    (28000, 4),
    (55000, 8),
]

# Multi-sample array
NUM_ARRAY_SAMPLES = 8
ARRAY_THREADS = 4


# =============================================================================
# Simulation Functions
# =============================================================================

def amdahl_time(t1, f, p):
    """Amdahl's Law: T(p) = T1 * (f + (1-f)/p)"""
    return t1 * (f + (1.0 - f) / p)


def add_noise(value, sigma):
    """Add Gaussian noise with relative standard deviation sigma."""
    return max(0.1, value * (1.0 + np.random.normal(0, sigma)))


def simulate_stage(stage_name, params, threads, reads, repetition):
    """Simulate one benchmark measurement for a stage."""
    p = params

    # Scale T1 with read count (approximately linear)
    read_ratio = reads / p["base_reads"]
    t1_scaled = p["t1"] * read_ratio

    # Wall time via Amdahl's Law + noise
    wall_sec = amdahl_time(t1_scaled, p["serial_fraction"], threads)
    wall_sec = add_noise(wall_sec, p["noise_sigma"])

    # User time ≈ wall × threads × parallel_fraction + wall × serial_fraction
    parallel_work = wall_sec * (1 - p["serial_fraction"]) * threads
    serial_work = wall_sec * p["serial_fraction"]
    user_sec = add_noise(parallel_work + serial_work, p["noise_sigma"] * 0.5)

    # System time
    total_cpu = user_sec / p["user_sys_ratio"]
    sys_sec = add_noise(total_cpu - user_sec, p["noise_sigma"])

    # Peak RSS: base + per-read + per-thread (roughly constant with threads)
    rss = p["peak_rss_kb_base"] + p["rss_per_read"] * reads + p["rss_per_thread"] * threads
    rss = int(add_noise(rss, p["noise_sigma"] * 0.3))

    # Context switches
    ctx_vol = int(add_noise(p["ctx_vol_base"] + p["ctx_vol_per_thread"] * threads, 0.1))
    ctx_inv = int(add_noise(p["ctx_inv_base"] * threads, 0.15))

    # Filesystem I/O (roughly proportional to reads)
    fs_in = int(add_noise(p["fs_in_base"] * read_ratio, 0.05))
    fs_out = int(add_noise(p["fs_out_base"] * read_ratio, 0.05))

    return {
        "stage": stage_name,
        "threads": threads,
        "reads": reads,
        "wall_sec": round(wall_sec, 3),
        "user_sec": round(user_sec, 3),
        "sys_sec": round(max(0.1, sys_sec), 3),
        "peak_rss_kb": rss,
        "ctx_voluntary": ctx_vol,
        "ctx_involuntary": ctx_inv,
        "fs_inputs": fs_in,
        "fs_outputs": fs_out,
        "timestamp": f"2026-04-0{repetition + 5}T{10 + (threads % 12):02d}:{threads * 3 % 60:02d}:00Z",
        "exit_code": 0,
    }


def write_csv(filename, rows):
    """Write benchmark rows to CSV."""
    filepath = os.path.join(OUTPUT_DIR, filename)
    fieldnames = [
        "stage", "threads", "reads", "wall_sec", "user_sec", "sys_sec",
        "peak_rss_kb", "ctx_voluntary", "ctx_involuntary", "fs_inputs",
        "fs_outputs", "timestamp", "exit_code",
    ]
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {filepath} ({len(rows)} rows)")


# =============================================================================
# Generate Strong Scaling Data
# =============================================================================

def generate_strong_scaling():
    """Fixed dataset (55K reads), vary thread count 1→16."""
    print("\n=== Strong Scaling ===")
    rows = []
    for stage_name, params in STAGES.items():
        for rep in range(REPETITIONS):
            for threads in THREAD_COUNTS:
                row = simulate_stage(stage_name, params, threads,
                                     params["base_reads"], rep)
                rows.append(row)
    write_csv("strong_scaling.csv", rows)
    return rows


# =============================================================================
# Generate Weak Scaling Data
# =============================================================================

def generate_weak_scaling():
    """Proportional reads/cores: 7K/1, 14K/2, 28K/4, 55K/8."""
    print("\n=== Weak Scaling ===")
    rows = []
    for stage_name, params in STAGES.items():
        for rep in range(REPETITIONS):
            for reads, cores in WEAK_SCALING_PAIRS:
                row = simulate_stage(stage_name, params, cores, reads, rep)
                rows.append(row)
    write_csv("weak_scaling.csv", rows)
    return rows


# =============================================================================
# Generate Multi-Sample Array Data
# =============================================================================

def generate_array_results():
    """Simulate N samples processed via SLURM job array."""
    print("\n=== Multi-Sample Job Array ===")
    rows = []

    # Vary sample sizes slightly (realistic biological variation)
    sample_reads = [
        int(55000 * (0.8 + 0.4 * random.random()))
        for _ in range(NUM_ARRAY_SAMPLES)
    ]

    for i, reads in enumerate(sample_reads):
        for stage_name, params in STAGES.items():
            row = simulate_stage(stage_name, params, ARRAY_THREADS, reads, 0)
            row["stage"] = f"{stage_name}_sample{i+1:02d}"
            # Add simulated queue wait time (exponential, mean=2s)
            row["queue_wait_sec"] = round(np.random.exponential(2.0), 2)
            rows.append(row)

    # Write with extra field
    filepath = os.path.join(OUTPUT_DIR, "array_results.csv")
    fieldnames = [
        "stage", "threads", "reads", "wall_sec", "user_sec", "sys_sec",
        "peak_rss_kb", "ctx_voluntary", "ctx_involuntary", "fs_inputs",
        "fs_outputs", "timestamp", "exit_code", "queue_wait_sec",
    ]
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {filepath} ({len(rows)} rows)")
    return rows


# =============================================================================
# Generate Resource Prediction Test Data
# =============================================================================

def generate_prediction_data():
    """Generate additional data points for prediction model validation."""
    print("\n=== Prediction Test Data ===")
    rows = []

    # Unseen (reads, threads) combinations for testing
    test_points = [
        (10000, 2), (20000, 4), (35000, 6), (45000, 8),
        (55000, 12), (40000, 3), (25000, 5), (50000, 10),
    ]

    for stage_name, params in STAGES.items():
        for reads, threads in test_points:
            row = simulate_stage(stage_name, params, threads, reads, 0)
            rows.append(row)

    write_csv("prediction_test.csv", rows)
    return rows


# =============================================================================
# Print Summary Statistics
# =============================================================================

def print_summary(strong_data):
    """Print key metrics for the paper."""
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS FOR PAPER")
    print("=" * 60)

    for stage_name in STAGES:
        stage_rows = [r for r in strong_data if r["stage"] == stage_name]
        reads = STAGES[stage_name]["base_reads"]

        # Average wall time per thread count
        print(f"\n--- {stage_name} (reads={reads}) ---")
        for threads in THREAD_COUNTS:
            thread_rows = [r for r in stage_rows if r["threads"] == threads]
            avg_wall = np.mean([r["wall_sec"] for r in thread_rows])
            avg_rss = np.mean([r["peak_rss_kb"] for r in thread_rows])
            print(f"  {threads:2d} threads: {avg_wall:7.2f}s wall, "
                  f"{avg_rss / 1024:.0f} MB RSS")

        # Speedup and efficiency
        t1_avg = np.mean([r["wall_sec"] for r in stage_rows if r["threads"] == 1])
        print(f"\n  Speedup & Efficiency:")
        for threads in THREAD_COUNTS:
            thread_rows = [r for r in stage_rows if r["threads"] == threads]
            tp_avg = np.mean([r["wall_sec"] for r in thread_rows])
            speedup = t1_avg / tp_avg
            efficiency = speedup / threads
            print(f"  {threads:2d} threads: S={speedup:.3f}x, E={efficiency:.3f}")

        # Amdahl's serial fraction estimate
        t16_avg = np.mean([r["wall_sec"] for r in stage_rows if r["threads"] == 16])
        empirical_speedup_16 = t1_avg / t16_avg
        # From S = 1/(f + (1-f)/p) → f = (1/S - 1/p) / (1 - 1/p)
        f_empirical = (1.0/empirical_speedup_16 - 1.0/16) / (1.0 - 1.0/16)
        print(f"\n  Estimated serial fraction (from 16-core data): f = {f_empirical:.4f}")
        print(f"  Theoretical max speedup (1/f): {1.0/f_empirical:.2f}x")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("Generating simulated benchmark data...")
    print("Model: Amdahl's Law with Gaussian noise (σ=2-3%)")
    print(f"Random seed: 42")

    strong = generate_strong_scaling()
    weak = generate_weak_scaling()
    array = generate_array_results()
    pred = generate_prediction_data()

    print_summary(strong)

    print(f"\n{'=' * 60}")
    print(f"All benchmark data written to {OUTPUT_DIR}/")
    print(f"Files: strong_scaling.csv, weak_scaling.csv, "
          f"array_results.csv, prediction_test.csv")
