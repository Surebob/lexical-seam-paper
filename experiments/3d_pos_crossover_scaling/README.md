# Experiment 3d README

## Experiment id

`3d_pos_crossover_scaling`

## Research question

How does the POS closed-class crossover rank scale with vocabulary size across the 25 English corpora?

## Why this is one experiment

The historical POS runs all apply the same crossover operation to the same English corpus set. The automatic spaCy pass supplies the 25-corpus scaling law, while the manual POS bundle supplies validation points for a small subset.

## Upstream dependencies

- Historical automatic POS bundle: `/Volumes/External2TB/emlexperiment/results/zipf_pos_all_corpora/`
- Historical manual POS bundle: `/Volumes/External2TB/emlexperiment/results/zipf_pos_manual_v2/`

## Inputs

- `results/zipf_pos_all_corpora/pos_all_corpora_points.csv`
- `results/zipf_pos_all_corpora/summary.json`
- `results/zipf_pos_manual_v2/manual_alpha_points.csv`
- `results/zipf_pos_manual_v2/summary.json`

## Outputs produced

### `outputs/pos_scaling_points.csv`

One row per English corpus for the automatic POS crossover dataset.

Columns:
- `slug`
- `corpus`
- `vocabulary_size`
- `sqrt_v`
- `k_crossover`
- `alpha_per_corpus`
- `log_vocabulary_size`
- `log_k_crossover`
- `crossover_deviation_from_sqrt_v`
- `censored_at_top_n`

### `outputs/manual_validation.csv`

One row per manual validation alpha point.

Columns:
- `slug`
- `corpus`
- `source`
- `vocabulary_size`
- `manual_k_crossover`
- `manual_alpha_per_corpus`
- `spacy_k_crossover`
- `abs_manual_minus_spacy_rank_gap`
- `has_top300_manual_csv`

### `outputs/aggregate_statistics.csv`

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Rows:
- `top_pos_window`
- `crossover_fraction`
- `manual_validation_window`
- `manual_validation_corpus_count`
- `max_manual_vs_spacy_rank_gap`
- `forced_alpha`
- `forced_alpha_ci_low`
- `forced_alpha_ci_high`
- `forced_alpha_vs_half_pvalue`

### `outputs/manifest.json`

Machine-readable output inventory and audit flags.

## Manuscript lines fed

- line `13`: POS exponent and CI.
- line `119`: top-500 and top-300 POS protocol and manual validation summary.
- lines `261-263`: POS exponent, CI, and p-value.
- line `516`: POS exponent restatement.
- line `534`: manual cross-check count.
- line `660`: appendix verification checklist for POS alpha.

## Methods

This migration reads the historical automatic POS scaling output and manual validation output without rerunning POS tagging. The reported 25-corpus POS exponent, confidence interval, and p-value are the saved values from `zipf_pos_all_corpora/summary.json`.

CI method: Student-`t` 95% confidence interval on the forced alpha using the historical covariance standard error and `df = n - 1`.

P-value method: two-sided one-sample Student-`t` test on per-corpus POS alpha values against `0.5`, as saved in the historical bundle.

## AUDIT

The manuscript and claim map refer to a "manual cross-check on four of the 25 English corpora" and a rank-gap bound of `15` ranks. The historical data has four alpha points in `manual_alpha_points.csv` because Shakespeare is included as `manual POS v2`, but only three top-300 manual POS CSV files are present: Moby Dick, Federalist Papers, and Wealth of Nations.

The available data gives `max_manual_vs_spacy_rank_gap = 10`, not `15`. This migration preserves the data value and flags the mismatch rather than silently adjusting either the data or the manuscript.

## Canonical claim mapping

- `forced_alpha`, `forced_alpha_ci_low`, `forced_alpha_ci_high`, and `forced_alpha_vs_half_pvalue` map to `outputs/aggregate_statistics.csv`.
- `top_pos_window`, `crossover_fraction`, `manual_validation_window`, `manual_validation_corpus_count`, and `max_manual_vs_spacy_rank_gap` map to `outputs/aggregate_statistics.csv`.
- Per-corpus POS rows map to `outputs/pos_scaling_points.csv`.

## Rerun command

```bash
python3 experiments/3d_pos_crossover_scaling/scripts/run_3d.py \
  --output-dir experiments/3d_pos_crossover_scaling/outputs
```

## Verification mapping

- `pos_scaling_points.csv` must have `25` rows.
- `manual_validation.csv` must have `4` rows.
- `aggregate_statistics.csv` must contain exactly the claim-map-facing metric rows listed above.
- Audit flags in `manifest.json` must remain until manuscript text and claim-map expectations are reconciled.
