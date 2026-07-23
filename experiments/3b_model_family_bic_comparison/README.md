# Experiment 3b README

## Experiment id

`3b_model_family_bic_comparison`

## Research question

Which rank-curve model family is preferred by RMSE-based BIC when the baseline is the paper's actual 3-parameter single Zipf-Mandelbrot model?

## Why this is one experiment

The historical model-family bundles compare the same English corpus set under related rank-curve model families: single ZM, MOEZipf, hard and continuous piecewise ZM, fixed-`k` smooth ZM, and free-`k` smooth ZM. The migration consolidates those saved results into one per-corpus model-family table.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_bic_comparison/`
- `/Volumes/External2TB/emlexperiment/results/zipf_moezipf_comparison/`
- `/Volumes/External2TB/emlexperiment/results/zipf_continuous_piecewise/`
- `/Volumes/External2TB/emlexperiment/results/zipf_sqrt_v_all_corpora/`
- `/Volumes/External2TB/emlexperiment/results/zipf_bic_landscape/`

## Inputs

- `results/zipf_bic_comparison/summary.json`
- `results/zipf_moezipf_comparison/summary.json`
- `results/zipf_continuous_piecewise/summary.json`
- `results/zipf_sqrt_v_all_corpora/summary.json`

## Outputs produced

### `outputs/table2_model_family.csv`

One row per English corpus.

Columns:
- `slug`
- `corpus`
- `vocabulary_size`
- `single_zm_bic`
- `moe_bic`
- `hard_piecewise_bic`
- `continuous_piecewise_bic`
- `smooth_ksqrtv_bic`
- `smooth_freek_bic`
- `winner_family`
- `legacy_v4_zipf_rank_bic`
- `legacy_v4_zipf_rank_bic_delta_vs_single_zm_bic`

### `outputs/model_family_rmse_per_corpus.csv`

One row per English corpus.

Columns:
- `slug`
- `corpus`
- `vocabulary_size`
- `zipf_rank_rmse_legacy`
- `moezipf_rank_rmse`
- `single_zm_rank_rmse`
- `hard_piecewise_rank_rmse`
- `continuous_piecewise_rank_rmse`
- `smooth_ksqrtv_rank_rmse`
- `smooth_free_k_rank_rmse`
- `sqrtv_minus_free_k_rmse`

### `outputs/aggregate_statistics.csv`

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Rows:
- `piecewise_loses_to_smooth_on_all_25`
- `winner_count_smooth_family_total`
- `moezipf_bic_winner_count`
- `winner_count_moezipf`
- `median_rmse_cost_sqrtv_vs_free_k`

### `outputs/manifest.json`

Machine-readable output inventory and audit flags.

## Manuscript lines fed

- v4 lines `232-243`: historical Table 2 model-family table.
- v4 line `243`: smooth-family BIC winner count.
- v4 line `247`: Moby Dick MOEZipf and smooth rank-RMSE comparison.
- v4 line `257`: median RMSE cost of fixed `k = sqrt(V)` relative to free `k`.
- v4 line `275`: MOEZipf BIC winner count.

In v5.1, Section 3.5 Table 2 has already moved to S2 v3 decoupled gate-family outputs, so this legacy `3b` table is no longer the manuscript's canonical Table 2 source.

## Methods

This migration reads saved historical summaries only. It does not refit any model and does not rerun likelihood or rank-curve optimization.

The canonical BIC formula is:

```text
BIC = p * ln(n) + n * ln(MSE)
MSE = RMSE^2
n = vocabulary size
```

The canonical baseline column is `single_zm_bic`, the true 3-parameter single-ZM BIC from `zipf_bic_comparison`. MOEZipf BIC remains the rank-curve RMSE-BIC from `zipf_moezipf_comparison`, as imported through `zipf_continuous_piecewise`.

## AUDIT

The v4 Table 2 pathway mislabeled `zipf_continuous_piecewise.zipf_bic` as a ZM baseline. That field is actually the 1-parameter pure-Zipf rank-BIC from `zipf_moezipf_comparison.zipf_rank_bic`, not true single-ZM BIC.

Canonical `3b` corrects the structure by emitting true 3-parameter `single_zm_bic` from `zipf_bic_comparison`. The legacy v4 pure-Zipf value is retained only in `legacy_v4_zipf_rank_bic` for drift auditing.

Manuscript revision flag: every v5.1 citation of "single-ZM BIC" needs re-verification against `outputs/table2_model_family.csv`. If a cited value was carried over from v4's mislabeled pure-Zipf table, it must be updated to the true ZM-BIC convention or removed.

`results/zipf_bic_landscape/` contains only a PNG and no machine-readable table. It is preserved as provenance but supplies no canonical numeric value.

## Canonical claim mapping

- `single_zm_bic` maps to `outputs/table2_model_family.csv`.
- `moe_bic`, `hard_piecewise_bic`, `continuous_piecewise_bic`, `smooth_ksqrtv_bic`, and `smooth_freek_bic` map to `outputs/table2_model_family.csv`.
- Winner counts and fixed-`k` RMSE cost map to `outputs/aggregate_statistics.csv`.
- Moby Dick rank-RMSE claims map to `outputs/model_family_rmse_per_corpus.csv`.

## Rerun command

```bash
python3 experiments/3b_model_family_bic_comparison/scripts/run_3b.py \
  --output-dir experiments/3b_model_family_bic_comparison/outputs
```

## Verification mapping

- `table2_model_family.csv` must have `25` rows.
- `model_family_rmse_per_corpus.csv` must have `25` rows.
- `aggregate_statistics.csv` must contain the claim-map-facing metric rows listed above.
- `winner_count_moezipf` should be `13`, and `winner_count_smooth_family_total` should be `12`, even after correcting the first column to true single-ZM BIC.
