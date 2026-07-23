# Experiment 1b README

## Experiment id

`1b_non_linguistic_control`

## Research question

Does the Bregman residual pattern found in language appear in a non-linguistic Zipfian dataset, specifically world city populations?

## Why this is one experiment

The historical city-population bundles apply the English enriched-search pipeline to a single non-language control dataset. The step-2-only bundle supplies the canonical Section 3.2 diagnostic, while the full-search bundle is retained as terminal-search provenance.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_enriched_city_populations/`
- `/Volumes/External2TB/emlexperiment/results/zipf_enriched_city_populations_step2_only/`

## Inputs

- `results/zipf_enriched_city_populations/summary.json`
- `results/zipf_enriched_city_populations_step2_only/summary.json`

## Outputs produced

### `outputs/city_population_control.csv`

One row for the world city population control.

Columns:
- `dataset`
- `city_count`
- `vocab_size`
- `population_sum`
- `zm_a`
- `zm_b`
- `zm_c`
- `zm_rmse`
- `step2_winner_expr`
- `step2_winner_math`
- `step2_composite_rmse`
- `step2_helpful`
- `step2_delta_vs_zm`
- `is_bregman_rank_in_step2_beam`
- `is_bregman_rmse`
- `exp_bregman_rank_in_step2_beam`
- `exp_bregman_rmse`
- `terminal_winner_expr`
- `terminal_winner_math`
- `terminal_composite_rmse`
- `terminal_helpful`

### `outputs/aggregate_statistics.csv`

Header-only file. The v4 claim map assigns no aggregate rows to this experiment; Section 3.2 claims map directly to `city_population_control.csv`.

### `outputs/manifest.json`

Machine-readable output inventory and source-bundle provenance.

## Manuscript lines fed

- v4 line `210` and v5.1 line `224`: `33,535` cities, `V = 33,535`, single-ZM `c = 100.4`, ZM RMSE `0.1206`, and step-2 composite RMSE `0.1225`.

## Methods

This migration reads saved city-population summaries only. It does not download GeoNames data and does not rerun symbolic search.

## AUDIT

No data-claim mismatch was detected for the Section 3.2 line-item values. The step-2 winner is present in the saved beam and is not helpful relative to the single-ZM baseline.

## Canonical claim mapping

- Section 3.2 line-item claims map to `outputs/city_population_control.csv`, row `dataset=world_city_populations`.

## Rerun command

```bash
python3 experiments/1b_non_linguistic_control/scripts/run_1b.py \
  --output-dir experiments/1b_non_linguistic_control/outputs
```

## Verification mapping

- `city_population_control.csv` must have `1` row.
- `city_count` and `vocab_size` must both equal `33535`.
- `zm_c` must equal `100.42856857594688`.
- `zm_rmse` must equal `0.12063481439379493`.
- `step2_composite_rmse` must equal `0.12245636082251322`.
