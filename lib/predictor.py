#!/usr/bin/env python3
"""
lib/predictor.py — Adaptive resource prediction for SLURM job scheduling.

Reads benchmark CSV data from pipeline scaling experiments, fits Amdahl's Law
and linear memory models, and generates right-sized SLURM job scripts.

Usage:
    from predictor import ResourcePredictor
    pred = ResourcePredictor("benchmarks/strong_scaling.csv")
    pred.fit()
    resources = pred.predict(n_reads=100000, n_cores=8)
    pred.generate_sbatch(resources, "denoise", "jobs/denoise.sbatch")
"""

import csv
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit


class ResourcePredictor:
    """Fits Amdahl's Law and memory models from benchmark data."""

    def __init__(self, csv_path, additional_csvs=None):
        self.csv_path = csv_path
        self.data = defaultdict(list)
        self.models = {}
        self._load_data()
        if additional_csvs:
            for path in additional_csvs:
                self._load_csv(path)

    def _load_csv(self, csv_path):
        """Load a benchmark CSV and append to per-stage data.

        Rows whose exit_code is non-zero are discarded. A failed stage still
        produces a plausible-looking wall time -- a classifier that dies on a
        scikit-learn version mismatch, for instance, records the seconds it
        spent loading the model before raising -- and those numbers would
        otherwise be fitted as if they were successful runs.
        """
        skipped = 0
        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    if int(row.get("exit_code", 0) or 0) != 0:
                        skipped += 1
                        continue
                except ValueError:
                    skipped += 1
                    continue
                stage = row["stage"]
                self.data[stage].append({
                    "threads": int(row["threads"]),
                    "reads": int(row["reads"]),
                    "wall_sec": float(row["wall_sec"]),
                    "peak_rss_kb": int(row["peak_rss_kb"]),
                    "user_sec": float(row["user_sec"]),
                    "sys_sec": float(row["sys_sec"]),
                })
        if skipped:
            print(f"Note: skipped {skipped} failed run(s) (exit_code != 0) "
                  f"in {os.path.basename(csv_path)}", file=sys.stderr)

    def _load_data(self):
        """Load benchmark CSV into per-stage data dictionaries."""
        self._load_csv(self.csv_path)

    @staticmethod
    def amdahl_model(p, t1, f):
        """Amdahl's Law: T(p) = T1 * (f + (1-f)/p)

        Parameters:
            p: number of processors
            t1: single-core execution time
            f: serial fraction (0 <= f <= 1)

        Monotonic by construction: predicted time can only fall as p rises. Real
        stages measured here slow down past a peak core count, which this form
        cannot represent — see usl_model.
        """
        return t1 * (f + (1.0 - f) / p)

    @staticmethod
    def usl_model(p, t1, alpha, beta):
        """Universal Scalability Law (Gunther): T(p) = T1 * (1 + a(p-1) + b*p(p-1)) / p

        Speedup is S(p) = p / (1 + alpha(p-1) + beta*p(p-1)).

        Parameters:
            p:     number of processors
            t1:    single-core execution time
            alpha: contention / serialization coefficient (Amdahl-like term)
            beta:  coherency coefficient (cross-talk between workers)

        beta > 0 produces a peak followed by decline, which is what DADA2
        denoising actually does: measured speedup tops out near 8-16 threads and
        falls off by 40. Amdahl fits that curve at R^2 ~ 0.95; this fits at
        ~0.997 and recovers the optimum core count instead of implying that more
        cores are always at least as good.
        """
        return t1 * (1.0 + alpha * (p - 1) + beta * p * (p - 1)) / p

    @staticmethod
    def usl_optimum(alpha, beta, max_p=1024):
        """Core count maximizing USL speedup: p* = sqrt((1-alpha)/beta).

        Returns max_p when beta is ~0 (no coherency penalty, speedup saturates
        rather than declining).
        """
        if beta <= 0:
            return max_p
        p_star = math.sqrt(max(0.0, (1.0 - alpha)) / beta)
        return max(1, min(max_p, int(round(p_star))))

    @staticmethod
    def memory_model(n_reads, a, b):
        """Linear memory model: M(n) = a + b * n_reads"""
        return a + b * n_reads

    @staticmethod
    def _input_scaling(records):
        """Fit T ~ a * size^b over distinct input sizes AT A FIXED THREAD COUNT.

        Returns (exponent, r2, n_distinct_sizes).

        The fixed-thread restriction is the whole point. A weak-scaling series
        varies input size and thread count together by construction, so a fit
        pooled across it cannot separate the two: on this repository's own data
        that pooling yields an exponent of 0.17 for denoising and under-predicts
        a 1.8M-read run sevenfold. Restricted to the thread count with the most
        distinct input sizes, the same data give 0.82, matching the measured
        14.3-fold response over a 27-fold depth range.

        Under-prediction is the dangerous direction -- it produces a --time that
        kills the job -- so this refuses to fit rather than fit something
        confounded. A depth-series style calibration (input size varied at fixed
        threads) is what makes this identifiable.
        """
        by_threads = defaultdict(list)
        for r in records:
            by_threads[r["threads"]].append(r)

        best = max(by_threads.values(),
                   key=lambda rs: len({r["reads"] for r in rs}), default=[])
        uniq = sorted({r["reads"] for r in best})
        if len(uniq) < 3:
            return None, 0.0, len(uniq)

        means = np.array([np.mean([r["wall_sec"] for r in best if r["reads"] == u])
                          for u in uniq], dtype=float)
        x = np.log(np.array(uniq, dtype=float))
        y = np.log(means)
        if np.ptp(x) == 0 or np.ptp(y) == 0:
            return 0.0, 0.0, len(uniq)
        b, a = np.polyfit(x, y, 1)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2 = 1.0 - np.sum((y - (a + b * x)) ** 2) / ss_tot if ss_tot > 0 else 0.0
        return float(b), float(r2), len(uniq)

    def fit(self):
        """Fit Amdahl's Law (time) and linear (memory) models per stage."""
        for stage, records in self.data.items():
            # --- Time model (Amdahl's Law) ---
            # Use records with the same read count (strong scaling data)
            read_counts = set(r["reads"] for r in records)
            # Pick the most common read count for strong scaling fit
            main_reads = max(read_counts, key=lambda rc: sum(
                1 for r in records if r["reads"] == rc
            ))
            strong_data = [r for r in records if r["reads"] == main_reads]

            # Defaults for the USL half, overwritten below when the fit succeeds.
            alpha_fit, beta_fit, r2_usl, p_star = 0.0, 0.0, 0.0, 1
            best_model = "amdahl"

            if len(strong_data) >= 3:
                cores = np.array([r["threads"] for r in strong_data], dtype=float)
                times = np.array([r["wall_sec"] for r in strong_data], dtype=float)
                t1_init = times[np.argmin(cores)]

                def _r2(predicted):
                    ss_res = np.sum((times - predicted) ** 2)
                    ss_tot = np.sum((times - np.mean(times)) ** 2)
                    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

                try:
                    popt, pcov = curve_fit(
                        self.amdahl_model, cores, times,
                        p0=[t1_init, 0.2],
                        bounds=([0, 0], [np.inf, 1.0]),
                        maxfev=10000
                    )
                    t1_fit, f_fit = popt
                    r2_time = _r2(self.amdahl_model(cores, t1_fit, f_fit))
                except RuntimeError:
                    t1_fit = t1_init
                    f_fit = 0.5
                    r2_time = 0.0

                # USL needs at least four distinct core counts to separate the
                # contention term from the coherency term.
                if len(set(cores)) >= 4:
                    try:
                        popt_u, _ = curve_fit(
                            self.usl_model, cores, times,
                            p0=[t1_init, max(f_fit, 1e-3), 1e-4],
                            bounds=([0, 0, 0], [np.inf, 1.0, 1.0]),
                            maxfev=20000
                        )
                        t1_u, alpha_fit, beta_fit = popt_u
                        r2_usl = _r2(self.usl_model(cores, t1_u, alpha_fit, beta_fit))
                        p_star = self.usl_optimum(alpha_fit, beta_fit)
                        if r2_usl > r2_time:
                            best_model = "usl"
                            t1_fit = t1_u
                    except RuntimeError:
                        pass
            else:
                t1_fit = strong_data[0]["wall_sec"] if strong_data else 100.0
                f_fit = 0.5
                r2_time = 0.0

            # --- Memory model (linear in reads) ---
            all_reads = np.array([r["reads"] for r in records], dtype=float)
            all_mem = np.array([r["peak_rss_kb"] for r in records], dtype=float)

            # A linear fit needs at least two distinct read counts. A strong-
            # scaling CSV alone holds exactly one, so the regression is
            # unidentifiable and the slope collapses to zero. Flag that rather
            # than emitting a flat memory model that looks like a valid fit.
            n_distinct_reads = len(set(all_reads))
            mem_valid = False

            if n_distinct_reads >= 2:
                try:
                    popt_mem, _ = curve_fit(
                        self.memory_model, all_reads, all_mem,
                        p0=[500000, 10.0],
                        bounds=([0, 0], [np.inf, np.inf]),
                    )
                    a_fit, b_fit = popt_mem
                    predicted_mem = self.memory_model(all_reads, a_fit, b_fit)
                    ss_res_m = np.sum((all_mem - predicted_mem) ** 2)
                    ss_tot_m = np.sum((all_mem - np.mean(all_mem)) ** 2)
                    r2_mem = 1.0 - ss_res_m / ss_tot_m if ss_tot_m > 0 else 0.0
                    mem_valid = True
                except RuntimeError:
                    a_fit = np.mean(all_mem)
                    b_fit = 0.0
                    r2_mem = 0.0
            else:
                a_fit = np.mean(all_mem)
                b_fit = 0.0
                r2_mem = 0.0

            # How does this stage's cost respond to its own input size? The
            # `reads` column holds whatever that stage consumes -- read count
            # for denoising, ASV count for classification -- so this is fitted
            # per stage rather than assumed to be reads.
            exp_fit, exp_r2, n_sizes = self._input_scaling(records)

            self.models[stage] = {
                "t1": float(t1_fit),
                "input_exponent": exp_fit,
                "input_exponent_r2": float(exp_r2),
                "input_sizes_seen": int(n_sizes),
                "input_scaling_valid": exp_fit is not None,
                "serial_fraction": float(f_fit),
                "r2_time": float(r2_time),
                "usl_alpha": float(alpha_fit),
                "usl_beta": float(beta_fit),
                "r2_usl": float(r2_usl),
                "usl_optimum_cores": int(p_star),
                "best_model": best_model,
                "mem_intercept": float(a_fit),
                "mem_slope": float(b_fit),
                "r2_memory": float(r2_mem),
                "memory_model_valid": bool(mem_valid),
                "mem_distinct_read_counts": int(n_distinct_reads),
                "main_reads": int(main_reads),
            }

        return self.models

    def predict(self, stage, n_reads, n_cores):
        """Predict wall time and peak memory for given parameters.

        Returns dict with predicted_time, predicted_memory_kb,
        recommended SLURM parameters.
        """
        if stage not in self.models:
            raise ValueError(f"No model fitted for stage '{stage}'")

        model = self.models[stage]

        # Project T1 to the requested input size using the FITTED exponent.
        # If the calibration data contained too few distinct input sizes to fit
        # one, do not extrapolate at all -- returning the calibration time is
        # wrong by however much the stage really scales, but inventing an
        # exponent is wrong by orders of magnitude. See the warning below.
        ratio = n_reads / model["main_reads"] if model["main_reads"] else 1.0
        if model.get("input_scaling_valid") and model.get("input_exponent") is not None:
            t1_scaled = model["t1"] * (ratio ** model["input_exponent"])
        else:
            t1_scaled = model["t1"]
        # Use whichever form fit the strong-scaling sweep better. USL matters
        # when the stage slows down past a peak: Amdahl would under-predict the
        # wall time there and the job would hit its --time limit.
        if model.get("best_model") == "usl":
            predicted_time = self.usl_model(
                n_cores, t1_scaled, model["usl_alpha"], model["usl_beta"])
        else:
            predicted_time = self.amdahl_model(
                n_cores, t1_scaled, model["serial_fraction"])

        # Predict memory
        predicted_mem_kb = self.memory_model(n_reads, model["mem_intercept"], model["mem_slope"])

        # Add safety margins for SLURM: 20% time, 15% memory
        safe_time_sec = predicted_time * 1.2
        safe_mem_mb = (predicted_mem_kb * 1.15) / 1024.0

        # Convert to SLURM time format (HH:MM:SS)
        hours = int(safe_time_sec // 3600)
        mins = int((safe_time_sec % 3600) // 60)
        secs = int(safe_time_sec % 60)
        slurm_time = f"{hours:02d}:{mins:02d}:{secs:02d}"

        result = {
            "stage": stage,
            "n_reads": n_reads,
            "n_cores": n_cores,
            "predicted_time_sec": round(predicted_time, 2),
            "predicted_memory_kb": round(predicted_mem_kb, 0),
            "slurm_time": slurm_time,
            "slurm_mem_mb": int(math.ceil(safe_mem_mb)),
            "slurm_cpus": n_cores,
            "time_model": model.get("best_model", "amdahl"),
            "input_exponent": model.get("input_exponent"),
            "input_scaling_valid": bool(model.get("input_scaling_valid")),
            "memory_model_valid": model.get("memory_model_valid", True),
        }

        if not result["input_scaling_valid"]:
            msg = (f"WARNING [{stage}]: input scaling was not fitted -- the "
                   f"calibration data held only {model.get('input_sizes_seen', 1)} "
                   f"distinct input size(s). predicted_time_sec is the "
                   f"calibration time at {model['main_reads']} input units, NOT a "
                   f"projection to {n_reads}. Supply calibration data spanning a "
                   f"range of input sizes for this stage.")
            print(msg, file=sys.stderr)
            result.setdefault("warnings", []).append(msg)

        # An unfitted memory model returns the training mean at every input
        # size. Harmless at the calibration size, badly wrong away from it, so
        # say so loudly rather than letting it reach an sbatch --mem line.
        if not result["memory_model_valid"]:
            msg = (f"WARNING [{stage}]: memory model was not fitted "
                   f"(needs >= 2 distinct read counts). "
                   f"slurm_mem_mb={result['slurm_mem_mb']} is the calibration "
                   f"mean, not a prediction for {n_reads} reads.")
            print(msg, file=sys.stderr)
            result.setdefault("warnings", []).append(msg)

        return result

    def compute_waste_ratio(self, predicted, actual_time, actual_mem_kb):
        """Compute resource waste ratio: allocated / used.

        Values close to 1.0 indicate efficient allocation.
        Values >> 1.0 indicate over-provisioning (waste).
        """
        time_waste = (predicted["predicted_time_sec"] * 1.2) / actual_time if actual_time > 0 else float("inf")
        mem_waste = (predicted["predicted_memory_kb"] * 1.15) / actual_mem_kb if actual_mem_kb > 0 else float("inf")
        return {
            "time_waste_ratio": round(time_waste, 3),
            "memory_waste_ratio": round(mem_waste, 3),
            "combined_waste": round(time_waste * mem_waste, 3),
        }

    def generate_sbatch(self, prediction, job_command, output_path):
        """Generate a right-sized SLURM sbatch script."""
        content = f"""#!/bin/bash
#SBATCH --job-name={prediction['stage']}_adaptive
#SBATCH --cpus-per-task={prediction['slurm_cpus']}
#SBATCH --mem={prediction['slurm_mem_mb']}M
#SBATCH --time={prediction['slurm_time']}
#SBATCH --output=benchmarks/slurm_%j_{prediction['stage']}.out
#SBATCH --error=benchmarks/slurm_%j_{prediction['stage']}.err
#SBATCH --partition=main

# Auto-generated by predictor.py — adaptive resource allocation
# Predicted: {prediction['predicted_time_sec']:.1f}s wall, {prediction['predicted_memory_kb']:.0f}KB RSS
# Safety margins: +20% time, +15% memory

echo "=== Adaptive SLURM Job: {prediction['stage']} ==="
echo "Cores: {prediction['slurm_cpus']}, Memory: {prediction['slurm_mem_mb']}MB, Time limit: {prediction['slurm_time']}"
echo "Started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

{job_command}

echo "Finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
"""
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w") as f:
            f.write(content)
        os.chmod(output_path, 0o755)

    def save_models(self, output_path):
        """Save fitted models to JSON for reproducibility."""
        with open(output_path, "w") as f:
            json.dump(self.models, f, indent=2)

    def summary(self):
        """Print a human-readable summary of fitted models."""
        for stage, m in self.models.items():
            print(f"\n{'='*50}")
            print(f"Stage: {stage}")
            print(f"  Amdahl's Law:  T1={m['t1']:.2f}s, f={m['serial_fraction']:.4f} (R\u00b2={m['r2_time']:.4f})")
            print(f"  Max speedup:   {1.0/m['serial_fraction']:.2f}x (1/f limit)")

            if m.get("r2_usl"):
                print(f"  USL:           \u03b1={m['usl_alpha']:.4f}, \u03b2={m['usl_beta']:.6f} (R\u00b2={m['r2_usl']:.4f})")
                if m["usl_beta"] > 0:
                    print(f"  USL optimum:   {m['usl_optimum_cores']} cores (speedup declines beyond this)")
                print(f"  Better fit:    {m['best_model'].upper()}")

            # How cost responds to this stage's own input size. Fitted, never
            # assumed -- assuming 1.0 here is the error the paper is about.
            if m.get("input_scaling_valid"):
                print(f"  Input scaling: T ~ size^{m['input_exponent']:.3f} "
                      f"(R\u00b2={m['input_exponent_r2']:.4f}, over {m['input_sizes_seen']} distinct sizes)")
                if abs(m["input_exponent"]) < 0.05:
                    print(f"                 exponent \u2248 0: cost does NOT track input size")
                    print(f"                 over the range measured.")
            else:
                print(f"  Input scaling: NOT FITTED \u2014 only {m.get('input_sizes_seen', 1)} "
                      f"distinct input size(s)")
                print(f"                 predict() will NOT extrapolate on input size;")
                print(f"                 it returns the calibration time and warns.")

            if m.get("memory_model_valid", True):
                print(f"  Memory model:  {m['mem_intercept']:.0f} + {m['mem_slope']:.2f}\u00b7n_reads KB "
                      f"(R\u00b2={m['r2_memory']:.4f})")
            else:
                print(f"  Memory model:  NOT FITTED \u2014 only "
                      f"{m.get('mem_distinct_read_counts', 1)} distinct input size(s)")
                print(f"                 The reported intercept ({m['mem_intercept']:.0f} KB) is the")
                print(f"                 mean RSS and the slope is zero; it will not generalize.")
                print(f"                 Add a weak-scaling or depth-series CSV.")

            print(f"  Calibrated at: {m['main_reads']} input units")

def main():
    """CLI entry point for resource prediction."""
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        print("Usage: predictor.py <benchmark.csv> [more.csv ...] "
              "[--predict stage reads cores] [--save-models out.json]")
        print()
        print("Pass BOTH a strong-scaling and a weak-scaling CSV. The first")
        print("pins down the time model; only the second varies input size, so")
        print("without it the memory regression cannot be fitted at all.")
        sys.exit(1)

    # Everything before the first flag is an input CSV.
    csv_paths = []
    for a in args:
        if a.startswith("--"):
            break
        csv_paths.append(a)

    if not csv_paths:
        print("Error: no input CSV given.", file=sys.stderr)
        sys.exit(1)

    missing = [p for p in csv_paths if not os.path.isfile(p)]
    if missing:
        print(f"Error: input file(s) not found: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    pred = ResourcePredictor(csv_paths[0], additional_csvs=csv_paths[1:])
    models = pred.fit()
    pred.summary()

    if "--predict" in sys.argv:
        idx = sys.argv.index("--predict")
        stage = sys.argv[idx + 1]
        reads = int(sys.argv[idx + 2])
        cores = int(sys.argv[idx + 3])
        result = pred.predict(stage, reads, cores)
        print(f"\nPrediction for {stage} ({reads} reads, {cores} cores):")
        print(json.dumps(result, indent=2))

    if "--save-models" in sys.argv:
        idx = sys.argv.index("--save-models")
        out_path = sys.argv[idx + 1]
        pred.save_models(out_path)
        print(f"\nModels saved to {out_path}")


if __name__ == "__main__":
    main()
