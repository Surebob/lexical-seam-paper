# S2 v3 Windows Port

This bundle packages the decoupled S2 v3 five-gate sweep for a Windows Codex instance. It includes:

- the 9-parameter decoupled smoother with `w_gate` and `w_tail`
- the five gate families: `logistic`, `tanh`, `erf`, `algebraic`, `arctan`
- the 25 canonical English corpora
- Mac Step A reference outputs for reproducibility checking
- a smoke-test script and a checkpointed parallel full-sweep runner

## Dependencies

Create a fresh virtual environment and install the exact versions used in the Mac reference environment:

- `numpy==2.4.2`
- `scipy==1.17.1`
- `psutil==7.1.3`

`pandas` is intentionally not listed because it was not installed in the Mac reference environment and is not required by this bundle.

## Scripts

- `s2_v3_smoke_test.py`
  - runs Shakespeare only
  - seed `20260415`
  - `100` starts
  - checks all five BICs against `reference/expected_shakespeare_bic.json`
- `run_s2_v3_windows.py`
  - runs a pre-flight check
  - runs the Shakespeare smoke test unless `--skip-smoke` is supplied
  - stops after the smoke test unless `--full-sweep-confirmed` is passed
  - with `--full-sweep-confirmed`, launches the full `25`-corpus parallel sweep

## Parallelism

The full runner uses:

- `multiprocessing.Pool(processes=20, initializer=worker_init)`
- each worker fits one corpus and all five gates sequentially
- the worker initializer caps BLAS/OpenMP threads to `2` before importing NumPy/SciPy

## Outputs

Smoke test outputs:

- `outputs/smoke/s2_v3_smoke_per_fit_results.csv`
- `outputs/smoke/s2_v3_smoke_per_start_dispersion.csv`
- `outputs/smoke/s2_v3_smoke_comparison.csv`
- `outputs/smoke/s2_v3_smoke_summary.json`

Full sweep outputs:

- `outputs/full/s2_v3_per_fit_results.csv`
- `outputs/full/s2_v3_per_start_dispersion.csv`
- `outputs/full/s2_v3_per_corpus_results.csv`
- `outputs/full/s2_v3_tanh_calibration.csv`
- `outputs/full/s2_v3_aggregate_statistics.csv`
- `outputs/full/s2_v3_runtime_log.txt`
- `outputs/full/progress.log`

## Resume behavior

The full sweep is checkpointed per corpus. If the process crashes or the machine reboots, re-running the same command resumes from the already completed corpus rows in `outputs/full/`.

