# Measured resource scaling of a QIIME 2 amplicon pipeline

Code, data and instructions accompanying the manuscript *"Taxonomic
classification cost is invariant to sequencing depth: measured resource scaling
of a QIIME 2 amplicon pipeline"* (submitted to *Frontiers in Microbiology*,
Systems Microbiology).

Measured benchmark data can be found here: [`benchmarks/`](benchmarks/) 

---

## 0. What you need

| | |
|---|---|
| OS | Linux (measurements were made on Ubuntu 22.04.5, kernel 5.15.0) |
| QIIME 2 | Amplicon distribution. Measured under **2026.7** |
| CPU | Multi-core. Reported runs used 2 × Xeon Gold 6230 = 40 physical cores |
| RAM | ≥ 8 GB to run the pipeline; ~35 GB if you rebuild the classifier. We used a box with 256GB. |
| Disk | ~30 GB for the full benchmark set |
| Network | Downloads from ENA and SILVA |

### Install QIIME 2 — not `environment.yml`

```bash
# Use QIIME 2's own published environment file for the release you want (Recommended for furture protocol):
#   https://docs.qiime2.org  ->  Installing QIIME 2
conda activate <your-qiime2-env>
```

`environment.yml` in this repository is a `conda env export` of the exact
environment the benchmarks ran in. It exists so the reported numbers can be
**audited** against exact versions — not so the environment can be rebuilt from
it. All 684 of its dependencies carry linux-64 build hashes, so
`conda env create -f environment.yml` will fail on macOS or ARM, and will fail on
linux-64 once those builds rotate out of the channels.

Versions used: `qiime2` 2026.7.0, `rachis` 2026.7.0, `q2cli` 2026.7.0,
`q2-dada2` 2026.7.0 over `bioconductor-dada2` 1.38.0, scikit-learn 1.7.1,
Python 3.12.

`pipeline.conf` sets `QIIME_ENV="current"`, so the scripts use whatever
environment is active. No environment name has to match.

---

## 1. Get the reference classifier

QIIME 2 no longer distributes region-specific classifiers of its own — its
position is that full-length classifiers perform comparably — and the last SILVA
515F-806R artifact it hosted was serialised under scikit-learn 0.24.1, which
current releases refuse to load. From release 2026.4 QIIME 2 refers users to
external providers (SILVA, GTDB) instead. Try the download first:

```bash
bash stages/04_classify.sh --download-only
```

**Verify the primer region before going further.** A region mismatch is invisible
downstream except as quietly degraded classification:

```bash
unzip -p silva-138-99-515-806-nb-classifier.qza \
  '*/provenance/artifacts/*/action/action.yaml' | grep -A1 primer
```

A **full-length** classifier is also fine; QIIME 2 reports little difference
between full-length and region-specific training. What matters is that the
classifier is not extracted for a *different* region. If it is region-extracted,
either of these is correct — both are V4 (*E. coli* 515–806, ~253 bp):

```
GTGCCAGCMGCCGCGGTAA / GGACTACHVGGGTWTCTAAT   Caporaso 515F/806R   <- used here
GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT   Parada/Apprill (EMP)
```

`CCTACGGGNGGCWGCAG` is **341F/805R (V3–V4)** — a different region from the
amplicon, which is the failure this check exists for. Stop.

If loading fails on a scikit-learn version mismatch, train one locally
(~1.5–2.5 h, ~30 GB peak RAM). This is what was done for the manuscript:

```bash
tmux new -s classifier
bash scripts/build_classifier.sh --threads 40     # defaults to Caporaso primers
```

---

## 2. Stage the datasets

```bash
bash scripts/prepare_datasets.sh --list-only      # show the plan, download nothing
bash scripts/prepare_datasets.sh --threads 40     # ~60-90 min
```

Downloads, imports and denoises the nine runs in
[`datasets.tsv`](datasets.tsv), then writes
`.bench-scratch/dataset_richness.tsv`. **Check it before continuing** — the
`region` column must read `V4` for every row.

Expected (from ENA study PRJEB33591):

| Accession | Source | Read pairs | ASVs |
|---|---|---:|---:|
| ERR3444623 | middle ear | 62,444 | 5 |
| ERR3444633 | middle ear | 54,938 | 11 |
| ERR3444597 | nasopharynx | 70,630 | 16 |
| ERR3444606 | ear canal | 44,604 | 38 |
| ERR3444605 | nasopharynx | 66,993 | 54 |
| ERR3444642 | nasopharynx | 66,212 | 55 |
| ERR3444628 | nasopharynx | 1,804,054 | 58 |
| ERR3444680 | adenoid | 62,686 | 67 |
| ERR3444641 | adenoid | 82,257 | 70 |

---

## 3. Run the benchmarks

```bash
tmux new -s bench
sudo -v && bash scripts/run_real_benchmarks.sh --all --reps 3
```

`sudo -v` lets the harness drop the page cache between repetitions. Without it
you get warm-cache I/O and the script says so at the end — either is fine, but
Methods must state which. Use `tmux`; a dropped SSH session kills the run.

Individual arms via `--experiment NAME` (repeatable):

| Name | Measures | Manuscript |
|---|---|---|
| `depth-pair` | Two samples, same body site, ~27× depth difference | Table 2, Fig. 1 |
| `depth-series` | **One library** subsampled across depths, community fixed | Table 5, Fig. 5 |
| `per-sample` | Both stages per dataset, cost vs. its own richness | Fig. 6 |
| `asv-sweep` | Classifier cost vs. ASV count | Table 3, Fig. 2 |
| `strong-denoise` | DADA2 speedup vs. core count | Table 4, Fig. 3 |
| `strong-classify` | Classifier speedup vs. core count | Table 4, Fig. 3 |
| `determinism` | Are ASV sequences identical at 1, 4 and 40 threads? | "Output is unchanged by thread count" |
| `weak-denoise` | Wall time at constant reads/core; memory model | Fig. 4 |
| `array` | Multi-sample: sequential loop vs. SLURM job array | "Resource recommendations" |
| `numa` | One socket vs. spanning sockets (not in `--all`) | — |

Two controls that a naive timing loop omits, and without which the numbers are
wrong rather than merely noisy:

- **BLAS thread pools are pinned to 1.** numpy/scipy otherwise spawn their own
  pool regardless of `--p-n-threads`. Unpinned, a nominally single-threaded
  denoise consumed 105 s of CPU in 58 s wall — concurrency 1.8, which would
  understate every speedup reported.
- **The core ladder is built from physical cores** via `lscpu`, not logical
  CPUs. Hyperthread siblings share execution units.

`depth-series` exists because `depth-pair` compares two *different participants*
(D-003 and D-006), so community composition varies alongside depth. Subsampling
one library removes that confound — and it also shows the effect is not simply
richness saturating: ASV count rose 2.8× over the range while classification time
did not move.

---

## 3b. Generate the biological results

```bash
bash scripts/generate_biology.sh --threads 40      # ~15 min
```

Classifies every sample, writes composition and diversity per body site, and
checks whether taxonomic *assignments* change with thread count. Expected:

| Site | ASVs | Shannon | Dominant genus |
|---|---|---|---|
| middle ear (n=3) | 5–38 | 0.01–1.13 | *Alloiococcus* / *Haemophilus* |
| nasopharynx (n=4) | 16–58 | 1.01–1.49 | *Moraxella* |
| adenoid (n=2) | 67–70 | 2.58–2.62 | *Haemophilus* |

`biology/taxonomy_determinism.tsv` should show one identical SHA-256 across all
three thread counts. If it does not, the resource recommendations change the
reported community and must be reframed as a trade rather than a saving.

---

## 4. Fit models and regenerate figures

```bash
python3 lib/predictor.py benchmarks/strong_scaling.csv benchmarks/weak_scaling.csv \
        benchmarks/asv_sweep.csv benchmarks/depth_series.csv \
        --save-models benchmarks/fitted_models.json
python3 scripts/test_predictor.py          # regression tests
python3 paper/analysis/generate_figures.py
```

Pass **all four** CSVs. A strong-scaling sweep holds one read count by
construction, so the memory regression is unidentifiable from it alone and the
predictor reports `NOT FITTED` instead of emitting a flat model that looks valid.

The weak-scaling series is *not* the right source for the memory model either:
it varies read count in proportion to thread count, so the two cannot be
separated (`corr = 1.000`). The manuscript fits memory on `depth_series.csv`,
where threads are held at four, giving 1.006 GiB + 2.254 MiB per 1k reads at
R² = 0.995 against R² = 0.873 for the confounded fit.

---

## 5. Check you got the same answer

| Quantity | Expected | Source |
|---|---|---|
| Denoise wall, 66,993 reads @ 4 cores | 48.84 s | `depth_pair.csv` |
| Denoise wall, 1,804,054 reads @ 4 cores | 626.14 s | `depth_pair.csv` |
| Classify wall, both of the above | 41.02 / 40.64 s | `depth_pair.csv` |
| Denoise `T₁` @ 66,993 reads | 56.61 s | `strong_scaling.csv` |
| Denoise peak speedup | 1.184× at 8 threads | `strong_scaling.csv` |
| Classify speedup at 40 jobs | **0.782×** (slower than serial) | `strong_scaling.csv` |
| Classify wall, 1 ASV → 218 ASVs | 36.20 → 37.27 s | `asv_sweep.csv` |
| **Depth series**, denoise 66,816 → 1,804,054 reads | 43.67 → 623.12 s (**14.3×**) | `depth_series.csv` |
| **Depth series**, classify over the same range | 40.63 → 41.18 s (**1.01×**) | `depth_series.csv` |
| **Depth series**, ASVs recovered | 21 → 58 (2.8×) | `depth_series.csv` |
| Per-sample denoise vs reads | exponent 0.75 (8 shallow alone: R²=0.67) | `per_sample.csv` |
| Per-sample classify vs ASVs | R² = 0.003, 0.46 s spread over 5–70 ASVs | `per_sample.csv` |
| ASV sequences at 1 / 4 / 40 threads | identical, sha256 `c759115c59be27c5` | `determinism.tsv` |
| Taxonomic assignments at 1 / 4 / 40 threads | identical, sha256 `11b2812c…` | `biology/taxonomy_determinism.tsv` |
| Classify cost, 5-ASV effusion vs 70-ASV adenoid | 40.81 vs 40.79 s (**1.000×**) | `per_sample.csv` + `biology/` |
| Classifier memory, any condition | 3.11 GB | `asv_sweep.csv` |
| Memory model (denoise) | 1.006 GiB + 2.254 MiB/1k reads, R²=0.995 | `depth_series.csv` |
| USL fit (denoise) | α = 0.818, β = 0.000916, p\* = 14 | `fitted_models.json` |

Absolute times are hardware-specific. The **ratios** are the claims: denoise and
classify response to depth (12.8× vs 0.99×), and classify speedup below 1.0 at
every core count above one.

Failed runs are excluded from all fits, but check for them explicitly — a stage
that dies still records a plausible-looking wall time:

```bash
awk -F, 'FNR>1 && $NF!=0 {print FILENAME": "$1", exit="$NF}' benchmarks/*.csv
```

---

## 6. Repository map

```
edna.sh                       pipeline driver (stages 1-6)
stages/                       one script per pipeline stage
  slurm/                      sbatch templates and submission scripts
lib/
  common.sh                   config, arg parsing, QIIME 2 activation
  benchmark.sh                /usr/bin/time -v wrapper -> CSV
  predictor.py                Amdahl + USL time models, linear memory model
  cache.sh                    keyed cache over stage outputs
scripts/
  prepare_datasets.sh         download, import, denoise, report richness
  run_real_benchmarks.sh      the measurement harness
  build_classifier.sh         train a 515F-806R classifier with RESCRIPt
  simulate_benchmarks.py      SUPERSEDED - simulation only, do not use
benchmarks/                   measured CSVs + fitted_models.json
datasets.tsv                  the nine runs, with exclusion reasoning
environment.yml               provenance record - see section 0
```

---

## 7. Troubleshooting

**`BadGzipFile: Not a gzipped file (b'@E')`** — `q2-types` opens
manifest-referenced FASTQs with gzip unconditionally. The manifests here point at
`.fastq.gz`; `stages/02_import.sh` recompresses if only the plain file is present.

**`Missing option '--o-base-transition-stats'`** or
**`'--o-read-extraction-stats'`** — QIIME 2 2026.x promoted these outputs to
required; older releases reject the flags as unknown. `q2_extra_output()` in
`lib/common.sh` probes for them, so the same code runs on both.

**`The scikit-learn version (0.24.1) ... does not match`** — see section 1;
rebuild with `scripts/build_classifier.sh`.

**`QIIME_ENV is 'current' but no 'qiime' command is on PATH`** — activate your
QIIME 2 environment first.

**Denoise slower with more threads** — expected, and one of the paper's results.

---

## License

See [LICENSE](LICENSE). Sequence data belong to the original depositors and are
subject to ENA terms.
