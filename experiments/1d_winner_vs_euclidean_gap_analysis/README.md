# Experiment 1d README

## Experiment id

`1d_winner_vs_euclidean_gap_analysis`

## Research question

How large is the full-RMSE separation between the selected step-2 Bregman-type winner and the Euclidean Bregman candidate across the 25 English corpora?

## Why this is one experiment

The historical bundle computes the same winner-vs-core-generator comparison for all 25 English corpora: winner, Euclidean candidate, Itakura-Saito candidate, and exponential Bregman candidate.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_english_gap_verify/`

## Inputs

- `results/zipf_english_gap_verify/summary.json`
- `results/zipf_english_gap_verify/english_gap_table.csv`

## Outputs produced

### `outputs/gap_analysis_per_corpus.csv`

One row per English corpus.

Columns:
- `slug`
- `corpus`
- `winner_expr`
- `winner_rmse`
- `euclidean_rmse`
- `is_bregman_rmse`
- `exp_bregman_rmse`
- `gap_vs_euclidean`
- `gap_vs_is_bregman`
- `gap_vs_exp_bregman`
- `winner_is_euclidean`
- `max_gap_to_core_bregman_family`
- `core_family_gap_below_0p01`

### `outputs/aggregate_statistics.csv`

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Rows:
- `median_winner_minus_euclidean_gap_full_rmse`
- `euclidean_step2_winner_count`
- `max_core_generator_gap_full_rmse_typical_band`
- `winner_family_gap_below_0p01_count`

### `outputs/manifest.json`

Machine-readable output inventory and source-bundle provenance.

## Manuscript lines fed

- v4 line `206` and v5.1 line `220`: median gap `0.006`, Euclidean never winning, and core candidate family typically within `0.01`.
- v4 line `534`: full-RMSE gaps typically below `0.01`.

## Methods

This migration reads the saved 25-corpus gap verification bundle and normalizes it to canonical CSVs. No symbolic search is rerun.

## AUDIT

No data-claim mismatch was detected for the full-RMSE Euclidean-gap values. Head-only Euclidean-gap claims are not sourced from this experiment; those are imported into `2c` from simulation-recovery/mechanism outputs.

## Canonical claim mapping

- `median_winner_minus_euclidean_gap_full_rmse` maps to `outputs/aggregate_statistics.csv`.
- `euclidean_step2_winner_count` maps to `outputs/aggregate_statistics.csv`.
- `max_core_generator_gap_full_rmse_typical_band` maps to `outputs/aggregate_statistics.csv`.
- `winner_family_gap_below_0p01_count` maps to `outputs/aggregate_statistics.csv`.

## Rerun command

```bash
python3 experiments/1d_winner_vs_euclidean_gap_analysis/scripts/run_1d.py \
  --output-dir experiments/1d_winner_vs_euclidean_gap_analysis/outputs
```

## Verification mapping

- `gap_analysis_per_corpus.csv` must have `25` rows.
- `aggregate_statistics.csv` must have exactly the four claim-map-facing metric rows listed above.
- `euclidean_step2_winner_count` must be `0`.
- `median_winner_minus_euclidean_gap_full_rmse` must be `0.005951682350291393`.
