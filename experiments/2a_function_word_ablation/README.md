# Experiment 2a README

## Experiment id

`2a_function_word_ablation`

## Research question

Does removing or isolating high-frequency function-word-like tokens destroy the Bregman residual pattern in Shakespeare?

## Why this is one experiment

The historical bundle applies the same ablation logic to three Shakespeare variants: remove top 50, remove top 100, and keep only top 100.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_function_word_test/`

## Inputs

- `results/zipf_function_word_test/summary.json`
- `results/zipf_function_word_test/report.md`
- Historical subrun summaries under `remove_top50`, `remove_top100`, and `only_top100`

## Outputs produced

### `outputs/function_word_ablation.csv`

One row per ablation.

Columns:
- `ablation`
- `historical_case`
- `token_count`
- `unique_words`
- `zm_a`
- `zm_b`
- `zm_c`
- `zm_rmse`
- `step2_winner`
- `step2_winner_display`
- `step2_rmse`
- `step2_helpful`
- `rmse_delta_step2_minus_zm`

### `outputs/aggregate_statistics.csv`

Header-only file. The v4 claim map assigns no aggregate rows to this experiment; manuscript claims map directly to `function_word_ablation.csv`.

### `outputs/manifest.json`

Machine-readable output inventory and source-bundle provenance.

## Manuscript lines fed

- line `216`: `c = 607.1`, `c = 454.1`, and `c = 5.3 with RMSE 0.044` for the three Shakespeare ablations.

## Methods

This migration reads the saved historical `zipf_function_word_test/summary.json` and normalizes the three historical cases into a single canonical CSV. No corpus text is reread and no symbolic search is rerun.

## AUDIT

No data-claim mismatch was detected during migration. The three manuscript values are present in the historical summary.

## Canonical claim mapping

- `remove_top_100` maps to `outputs/function_word_ablation.csv`, row `ablation=remove_top_100`.
- `remove_top_50` maps to `outputs/function_word_ablation.csv`, row `ablation=remove_top_50`.
- `top_100_only` maps to `outputs/function_word_ablation.csv`, row `ablation=top_100_only`.

## Rerun command

```bash
python3 experiments/2a_function_word_ablation/scripts/run_2a.py \
  --output-dir experiments/2a_function_word_ablation/outputs
```

## Verification mapping

- `function_word_ablation.csv` must have `3` rows.
- The `top_100_only` row must have `zm_c = 5.320352416124492` and `zm_rmse = 0.04393090866192364`.
- The `remove_top_100` row must have `zm_c = 607.061047283561`.
- The `remove_top_50` row must have `zm_c = 454.069761577837`.
