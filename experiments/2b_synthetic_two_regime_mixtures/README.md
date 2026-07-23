# Experiment 2b README

## Experiment id

`2b_synthetic_two_regime_mixtures`

## Research question

Can simple two-regime synthetic mixtures reproduce the Bregman residual pattern, and does the planted mixture continue to expose residual structure after a MOEZipf fit?

## Why this is one experiment

The saved synthetic mixture bundle contains the primary two-regime mixture tests. The planted-mixture follow-up was historically saved in `zipf_breakthrough_probe`, but the manuscript narrative attributes it to the synthetic-mixture story, so Phase 2 migrates it here as a named subexperiment.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_synthetic_mixture/`
- `/Volumes/External2TB/emlexperiment/results/zipf_breakthrough_probe/`

## Outputs produced

### `outputs/synthetic_mixture_runs.csv`

Rows from the saved exact step-2 synthetic-mixture bundle.

### `outputs/planted_mixture_runs.csv`

The planted-mixture row historically saved in `zipf_breakthrough_probe`.

### `outputs/aggregate_statistics.csv`

Rows:
- `synthetic_mixture_configuration_count`

### `outputs/manifest.json`

Machine-readable output inventory, source-bundle provenance, audit flags, and blocked items.

## Manuscript lines fed

- lines `115`, `222`, `224`, and `226`: synthetic mixture design parameters.
- line `277`: planted-mixture RMSE values and `log(x)^2` winner.
- line `534`: synthetic mixture configuration count.

## Methods

This migration reads historical JSON summaries and normalizes them into canonical CSVs. It does not regenerate synthetic data and does not rerun symbolic search.

## AUDIT

The v4 manuscript/map claims a same-exponent control with `(alpha_1 = alpha_2 = 1.5)`. The historical source bundle does not contain that row. It contains `control_two_copies_alpha0.8`, migrated as `historical_same_exponent_control_alpha0p8`.

The planted-mixture line `277` values live in `zipf_breakthrough_probe`, not in the original synthetic-mixture bundle. This migration intentionally places those values in `2b/outputs/planted_mixture_runs.csv` to make the Section 3.7 narrative attribution correct in the Phase 2 structure.

## BLOCKED

See `BLOCKED.md` for the missing `alpha_1 = alpha_2 = 1.5` same-exponent control.

## Rerun command

```bash
python3 experiments/2b_synthetic_two_regime_mixtures/scripts/run_2b.py \
  --output-dir experiments/2b_synthetic_two_regime_mixtures/outputs
```

## Verification mapping

- `synthetic_mixture_runs.csv` must contain the saved `small_gap`, `medium_gap`, `large_gap`, and historical `alpha0.8/alpha0.8` control rows.
- `planted_mixture_runs.csv` must contain one row with ZM RMSE `0.1810649142555675`, ZM+step2 RMSE `0.17405431027139098`, MOE RMSE `0.3257257480868418`, MOE+step2 RMSE `0.262697670409533`, and winner display `log(x)^2`.
