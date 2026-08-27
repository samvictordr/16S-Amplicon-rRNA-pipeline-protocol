#!/usr/bin/env python3
"""
scripts/simulate_prediction.py — Train resource prediction models on simulated
benchmark data and evaluate adaptive vs static allocation waste.

Reads strong_scaling.csv and prediction_test.csv, fits Amdahl's Law + linear
memory models via lib/predictor.py, then computes waste ratios comparing:
  - Static allocation: always request max observed resources
  - Adaptive prediction: request predicted resources + safety margin

Outputs:
  - benchmarks/prediction_accuracy.csv — predicted vs actual per test point
  - benchmarks/waste_comparison.csv — static vs adaptive waste ratios
  - benchmarks/fitted_models.json — serialized model parameters

Usage: python scripts/simulate_prediction.py
"""

import csv
import json
import os
import sys

import numpy as np

# Add project root to path for lib imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib.predictor import ResourcePredictor

OUTPUT_DIR = "benchmarks"
STRONG_CSV = os.path.join(OUTPUT_DIR, "strong_scaling.csv")
WEAK_CSV = os.path.join(OUTPUT_DIR, "weak_scaling.csv")
TEST_CSV = os.path.join(OUTPUT_DIR, "prediction_test.csv")


def load_test_data(csv_path):
    """Load test data points."""
    data = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                "stage": row["stage"],
                "threads": int(row["threads"]),
                "reads": int(row["reads"]),
                "actual_time": float(row["wall_sec"]),
                "actual_mem_kb": int(row["peak_rss_kb"]),
            })
    return data


def compute_static_allocation(strong_csv):
    """Compute static (worst-case) resource allocation from training data."""
    max_time = {}
    max_mem = {}
    with open(strong_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            stage = row["stage"]
            wall = float(row["wall_sec"])
            mem = int(row["peak_rss_kb"])
            if stage not in max_time or wall > max_time[stage]:
                max_time[stage] = wall
            if stage not in max_mem or mem > max_mem[stage]:
                max_mem[stage] = mem

    return max_time, max_mem


def evaluate_predictions(predictor, test_data, static_time, static_mem):
    """Evaluate prediction accuracy and waste ratios."""
    accuracy_rows = []
    waste_rows = []

    for point in test_data:
        stage = point["stage"]
        if stage not in predictor.models:
            continue

        # Get adaptive prediction
        pred = predictor.predict(stage, point["reads"], point["threads"])

        # Prediction accuracy
        time_error = abs(pred["predicted_time_sec"] - point["actual_time"])
        time_pct_error = (time_error / point["actual_time"]) * 100
        mem_error = abs(pred["predicted_memory_kb"] - point["actual_mem_kb"])
        mem_pct_error = (mem_error / point["actual_mem_kb"]) * 100

        accuracy_rows.append({
            "stage": stage,
            "threads": point["threads"],
            "reads": point["reads"],
            "actual_time_sec": round(point["actual_time"], 3),
            "predicted_time_sec": pred["predicted_time_sec"],
            "time_error_pct": round(time_pct_error, 2),
            "actual_mem_kb": point["actual_mem_kb"],
            "predicted_mem_kb": int(pred["predicted_memory_kb"]),
            "mem_error_pct": round(mem_pct_error, 2),
        })

        # Waste ratios
        # Adaptive: predicted + 20% safety margin for time, +15% for memory
        adaptive_time_alloc = pred["predicted_time_sec"] * 1.2
        adaptive_mem_alloc = pred["predicted_memory_kb"] * 1.15

        # Static: always request the maximum seen in training
        static_time_alloc = static_time.get(stage, pred["predicted_time_sec"]) * 1.5
        static_mem_alloc = static_mem.get(stage, pred["predicted_memory_kb"]) * 1.5

        waste_rows.append({
            "stage": stage,
            "threads": point["threads"],
            "reads": point["reads"],
            "actual_time_sec": round(point["actual_time"], 3),
            "actual_mem_kb": point["actual_mem_kb"],
            "static_time_alloc": round(static_time_alloc, 2),
            "adaptive_time_alloc": round(adaptive_time_alloc, 2),
            "static_time_waste": round(static_time_alloc / point["actual_time"], 3),
            "adaptive_time_waste": round(adaptive_time_alloc / point["actual_time"], 3),
            "static_mem_alloc": int(static_mem_alloc),
            "adaptive_mem_alloc": int(adaptive_mem_alloc),
            "static_mem_waste": round(static_mem_alloc / point["actual_mem_kb"], 3),
            "adaptive_mem_waste": round(adaptive_mem_alloc / point["actual_mem_kb"], 3),
        })

    return accuracy_rows, waste_rows


def compute_r2(actual, predicted):
    """Compute R² (coefficient of determination)."""
    actual = np.array(actual)
    predicted = np.array(predicted)
    ss_res = np.sum((actual - predicted) ** 2)
    ss_tot = np.sum((actual - np.mean(actual)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def main():
    print("=" * 60)
    print("RESOURCE PREDICTION EVALUATION")
    print("=" * 60)

    # 1. Fit models on strong + weak scaling data
    print("\n[1] Fitting models on strong + weak scaling data...")
    predictor = ResourcePredictor(STRONG_CSV, additional_csvs=[WEAK_CSV])
    models = predictor.fit()
    predictor.summary()

    # Save fitted models
    models_path = os.path.join(OUTPUT_DIR, "fitted_models.json")
    predictor.save_models(models_path)
    print(f"\nModels saved to: {models_path}")

    # 2. Load test data
    print(f"\n[2] Loading test data from {TEST_CSV}...")
    test_data = load_test_data(TEST_CSV)
    print(f"  Test points: {len(test_data)}")

    # 3. Static allocation baselines
    print("\n[3] Computing static allocation baselines...")
    static_time, static_mem = compute_static_allocation(STRONG_CSV)
    for stage in static_time:
        print(f"  {stage}: max_time={static_time[stage]:.1f}s, "
              f"max_mem={static_mem[stage]/1024:.0f}MB")

    # 4. Evaluate predictions
    print("\n[4] Evaluating predictions...")
    accuracy_rows, waste_rows = evaluate_predictions(
        predictor, test_data, static_time, static_mem
    )

    # 5. Write prediction accuracy CSV
    acc_path = os.path.join(OUTPUT_DIR, "prediction_accuracy.csv")
    with open(acc_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=accuracy_rows[0].keys())
        writer.writeheader()
        writer.writerows(accuracy_rows)
    print(f"  Written: {acc_path} ({len(accuracy_rows)} rows)")

    # 6. Write waste comparison CSV
    waste_path = os.path.join(OUTPUT_DIR, "waste_comparison.csv")
    with open(waste_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=waste_rows[0].keys())
        writer.writeheader()
        writer.writerows(waste_rows)
    print(f"  Written: {waste_path} ({len(waste_rows)} rows)")

    # 7. Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    for stage in ["denoise", "classify"]:
        stage_acc = [r for r in accuracy_rows if r["stage"] == stage]
        stage_waste = [r for r in waste_rows if r["stage"] == stage]

        if not stage_acc:
            continue

        # Prediction accuracy
        actual_times = [r["actual_time_sec"] for r in stage_acc]
        pred_times = [r["predicted_time_sec"] for r in stage_acc]
        actual_mems = [r["actual_mem_kb"] for r in stage_acc]
        pred_mems = [r["predicted_mem_kb"] for r in stage_acc]

        r2_time = compute_r2(actual_times, pred_times)
        r2_mem = compute_r2(actual_mems, pred_mems)
        mean_time_err = np.mean([r["time_error_pct"] for r in stage_acc])
        mean_mem_err = np.mean([r["mem_error_pct"] for r in stage_acc])

        print(f"\n--- {stage} ---")
        print(f"  Time prediction:  R² = {r2_time:.4f}, "
              f"MAPE = {mean_time_err:.2f}%")
        print(f"  Memory prediction: R² = {r2_mem:.4f}, "
              f"MAPE = {mean_mem_err:.2f}%")

        # Waste ratios
        static_time_waste_avg = np.mean([r["static_time_waste"] for r in stage_waste])
        adaptive_time_waste_avg = np.mean([r["adaptive_time_waste"] for r in stage_waste])
        static_mem_waste_avg = np.mean([r["static_mem_waste"] for r in stage_waste])
        adaptive_mem_waste_avg = np.mean([r["adaptive_mem_waste"] for r in stage_waste])

        time_reduction = (1 - adaptive_time_waste_avg / static_time_waste_avg) * 100
        mem_reduction = (1 - adaptive_mem_waste_avg / static_mem_waste_avg) * 100

        print(f"\n  Waste ratios (allocated/used):")
        print(f"    Time  — Static: {static_time_waste_avg:.2f}x, "
              f"Adaptive: {adaptive_time_waste_avg:.2f}x "
              f"({time_reduction:.0f}% reduction)")
        print(f"    Memory — Static: {static_mem_waste_avg:.2f}x, "
              f"Adaptive: {adaptive_mem_waste_avg:.2f}x "
              f"({mem_reduction:.0f}% reduction)")

    print(f"\n{'=' * 60}")
    print("Prediction evaluation complete.")


if __name__ == "__main__":
    main()
