# Launch Prompt For Windows Codex

Use the provided bundle mechanically. Do not redesign the experiment.

## Goal

Run the decoupled S2 v3 five-gate sweep on Windows using the packaged code and corpora.

## Working directory

Change into the bundle root:

```powershell
cd s2_v3_windows_port
```

## Step 1: create a fresh virtual environment

```powershell
py -3 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Step 2: verify dependency versions

Run:

```powershell
python -c "import numpy, scipy, psutil; print(numpy.__version__, scipy.__version__, psutil.__version__)"
```

Expected:

- `2.4.2`
- `1.17.1`
- `7.1.3`

If versions do not match, stop and report.

## Step 3: run the smoke test first

Run:

```powershell
python s2_v3_smoke_test.py
```

This writes:

- `outputs/smoke/s2_v3_smoke_per_fit_results.csv`
- `outputs/smoke/s2_v3_smoke_per_start_dispersion.csv`
- `outputs/smoke/s2_v3_smoke_comparison.csv`
- `outputs/smoke/s2_v3_smoke_summary.json`

## Step 4: smoke-test pass criteria

The Shakespeare smoke test must match `reference/expected_shakespeare_bic.json` to `1e-6` absolute BIC tolerance for all five gates:

- logistic: `-107202.19167713833`
- tanh: `-107202.19168431366`
- erf: `-108101.40890745526`
- algebraic: `-106146.18281255414`
- arctan: `-105034.13460379837`

Also confirm the comparison CSV marks all rows as passing.

If any gate differs by more than `1e-6`, stop and report. Do not run the full sweep.

For extended reproducibility checks, use `reference/mac_step_a_results.csv` and compare the seed `20260415` rows for:

- `shakespeare`
- `aesops_fables`
- `moby_dick`

## Step 5: stop and ask for go signal

Do not launch the full sweep automatically after the smoke test.

Wait for explicit confirmation from the human, or use the explicit CLI flag below only if the human has already approved:

```powershell
python run_s2_v3_windows.py --full-sweep-confirmed
```

## Step 6: full sweep invocation

After approval, run:

```powershell
python run_s2_v3_windows.py --full-sweep-confirmed
```

Defaults:

- base seed `20260415`
- `100` starts
- `12000` max function evaluations
- `20` workers
- `2` BLAS threads per worker

## Step 7: outputs to inspect

Full sweep writes:

- `outputs/full/s2_v3_per_fit_results.csv`
- `outputs/full/s2_v3_per_start_dispersion.csv`
- `outputs/full/s2_v3_per_corpus_results.csv`
- `outputs/full/s2_v3_tanh_calibration.csv`
- `outputs/full/s2_v3_aggregate_statistics.csv`
- `outputs/full/s2_v3_runtime_log.txt`
- `outputs/full/progress.log`

## Step 8: stop conditions

Stop and report immediately if any of the following occur:

1. dependency versions do not match the pinned versions
2. smoke test fails the `1e-6` BIC tolerance on any gate
3. the full run writes non-finite or missing BIC values
4. `outputs/full/s2_v3_tanh_calibration.csv` shows any corpus with `tanh_calibration_pass = False`
5. any gate hits bounds suspiciously often and the aggregate file suggests a broad clipping issue
6. the process crashes repeatedly on the same corpus

## Step 9: resume behavior

The full sweep is checkpointed per corpus.

If the machine crashes or the process is interrupted, resume with the same command:

```powershell
python run_s2_v3_windows.py --full-sweep-confirmed
```

The runner reads existing CSVs in `outputs/full/` and skips completed corpora.

