# Experiment 2c README

## Experiment id

`2c_mechanism_absorption_residual_space`

## Research question

Does the smooth two-regime model absorb the Bregman seam where MOEZipf leaves symbolic residual structure?

## Why this is one experiment

The core historical bundle compares symbolic step-2 search on residuals after MOEZipf and after the smooth model, and also measures the low-dimensional head basis. The v4 claim map additionally assigns several Euclidean-gap comparison aggregates to this experiment to avoid builder-time joins; those are imported from the historical simulation-recovery bundle as saved values.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_breakthrough_probe/`
- `/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery/` for five v4-mapped Euclidean-gap comparison aggregates only.

## Outputs produced

### `outputs/residual_search_per_corpus.csv`

One row per English corpus comparing MOEZipf residual search and smooth-model residual search.

### `outputs/head_formula_competition.csv`

One row per English corpus for full/head-window formula competition against Euclidean.

### `outputs/head_basis_per_corpus.csv`

One row per English corpus for the top-200 cubic head-basis fit.

### `outputs/aggregate_statistics.csv`

Rows:
- `moe_residual_step2_help_count`
- `moe_residual_step2_gain_median`
- `smooth_residual_step2_gain_median`
- `smooth_residual_step2_hurt_count`
- `smooth_residual_step2_help_count`
- `smooth_residual_step2_nonhelp_count`
- `head_basis_r2_median_top200`
- `slant_vs_c_correlation`
- `median_winner_minus_euclidean_gap_full_rmse`
- `median_winner_minus_euclidean_gap_top100_rmse`
- `empirical_euclidean_gap_full_rmse_median`
- `empirical_euclidean_gap_top100_rmse_median`
- `smooth_synthetic_euclidean_gap_full_rmse_median`
- `smooth_synthetic_euclidean_gap_top100_rmse_median`
- `zm_control_euclidean_gap_top100_rmse_median`

## Manuscript lines fed

- line `13`: MOE residual search persists on `19/25`; smooth residual search hurts on `25/25`.
- line `31`: head-basis median `R^2` and slant/c correlation.
- lines `273`, `279`, `339`, `494`, `544`: mechanism absorption, Euclidean-gap, and residual-search summary claims.

## Methods

This migration reads the saved `zipf_breakthrough_probe/summary.json` for residual-search, formula-competition, head-basis, and phase-summary values. It imports the v4-mapped synthetic/empirical Euclidean-gap comparison rows from `zipf_simulation_recovery/summary.json` so the LaTeX builder will not need to join across experiments.

## AUDIT

The import from `zipf_simulation_recovery` is a claim-map accommodation, not a new experimental run. It should be revisited if the Phase 2 contract is further tightened so that each aggregate lives only in its natural experiment.

## Protocol constants

- `helpful_step2_gain_threshold = 0.001`
- `head_windows = [50, 100, 200]`

## Rerun command

```bash
python3 experiments/2c_mechanism_absorption_residual_space/scripts/run_2c.py \
  --output-dir experiments/2c_mechanism_absorption_residual_space/outputs
```

## Verification mapping

- `residual_search_per_corpus.csv` must have `25` rows.
- `head_formula_competition.csv` must have `25` rows.
- `head_basis_per_corpus.csv` must have `25` rows.
- `moe_residual_step2_help_count` must be `19`.
- `smooth_residual_step2_hurt_count` must be `25`.
