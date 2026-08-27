#!/usr/bin/env python3
"""
scripts/test_predictor.py — Regression tests for lib/predictor.py.

Run against the deposited benchmark CSVs:

    python3 scripts/test_predictor.py

These exist because the predictor previously shipped with a defect that made it
recommend `--time=449:47:28` for a job measured at 40.6 s: it scaled the classify
stage's runtime linearly with read count, which is exactly the heuristic this
project exists to show is wrong. The stage's `reads` column holds ASV count, not
reads, so the scaling factor was 1,804,054/54.

The tests below fix the contract that prevents that class of error returning:

  1. An exponent relating cost to input size is FITTED, never assumed.
  2. It is fitted only where thread count is held constant, because a
     weak-scaling design confounds input size with thread count.
  3. Where it cannot be identified, predict() does not extrapolate at all, and
     says so. Silence would be worse than a wrong number.
"""

import io
import os
import sys
import contextlib

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from predictor import ResourcePredictor  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BENCH = os.path.join(ROOT, "benchmarks")

# (stage, input size, measured wall seconds at 4 cores) taken from the CSVs.
MEASURED = [
    ("denoise_series",    1804054, 623.12),
    ("denoise_series",     902027, 308.58),
    ("denoise_series",      66816,  43.67),
    ("classify",          1804054,  40.64),
    ("classify",            66993,  36.68),
    ("classify_asvsweep",     218,  37.27),
    ("classify_asvsweep",       1,  36.20),
]

TOLERANCE = 0.35   # predictions must land within 35% of measurement


def build():
    csvs = ["strong_scaling.csv", "weak_scaling.csv",
            "asv_sweep.csv", "depth_series.csv"]
    paths = [os.path.join(BENCH, c) for c in csvs]
    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        print(f"SKIP: benchmark data absent: {', '.join(map(os.path.basename, missing))}")
        sys.exit(0)
    pred = ResourcePredictor(paths[0], additional_csvs=paths[1:])
    pred.fit()
    return pred


def main():
    pred = build()
    failures = []

    print("=" * 68)
    print("1. Predictions land near measurement")
    print("=" * 68)
    for stage, n, measured in MEASURED:
        if stage not in pred.models:
            continue
        with contextlib.redirect_stderr(io.StringIO()):
            r = pred.predict(stage, n, 4)
        ratio = r["predicted_time_sec"] / measured
        ok = abs(ratio - 1.0) <= TOLERANCE
        print(f"  {'PASS' if ok else 'FAIL'}  {stage:<20}{n:>9} -> "
              f"{r['predicted_time_sec']:>8.1f}s  measured {measured:>7.2f}s  "
              f"ratio {ratio:.2f}x")
        if not ok:
            failures.append(f"{stage}@{n}: predicted {r['predicted_time_sec']:.1f}s "
                            f"vs measured {measured:.2f}s ({ratio:.1f}x)")

    print()
    print("=" * 68)
    print("2. The original defect specifically: classification must not scale")
    print("   with read count")
    print("=" * 68)
    with contextlib.redirect_stderr(io.StringIO()):
        shallow = pred.predict("classify", 66993, 4)["predicted_time_sec"]
        deep = pred.predict("classify", 1804054, 4)["predicted_time_sec"]
    spread = max(shallow, deep) / min(shallow, deep)
    ok = spread < 1.5
    print(f"  {'PASS' if ok else 'FAIL'}  27x more reads changes the prediction "
          f"{spread:.2f}x (must stay below 1.5x)")
    if not ok:
        failures.append(f"classify prediction scales {spread:.0f}x with read count")

    print()
    print("=" * 68)
    print("3. An unidentifiable exponent is refused, not invented")
    print("=" * 68)
    for stage in ("classify", "denoise"):
        if stage not in pred.models:
            continue
        m = pred.models[stage]
        ok = not m["input_scaling_valid"]
        print(f"  {'PASS' if ok else 'FAIL'}  {stage:<10} "
              f"{m['input_sizes_seen']} distinct size(s) at fixed threads -> "
              f"{'refused' if ok else 'FITTED ANYWAY'}")
        if not ok:
            failures.append(f"{stage}: fitted an exponent from "
                            f"{m['input_sizes_seen']} input size(s)")

    print()
    print("=" * 68)
    print("4. Refusal is reported to the caller, not swallowed")
    print("=" * 68)
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        r = pred.predict("classify", 1804054, 4)
    warned = "warnings" in r and any("input scaling" in w for w in r["warnings"])
    on_stderr = "WARNING" in buf.getvalue()
    print(f"  {'PASS' if warned else 'FAIL'}  warning attached to the result dict")
    print(f"  {'PASS' if on_stderr else 'FAIL'}  warning written to stderr")
    if not (warned and on_stderr):
        failures.append("silent refusal: caller is not told extrapolation was skipped")

    print()
    print("=" * 68)
    if failures:
        print(f"FAILED ({len(failures)})")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("All predictor regression tests passed.")


if __name__ == "__main__":
    main()
