# Experiment 1c README

## Experiment id

`1c_search_robustness_ablations`

## Research question

Is the step-2 Bregman winner stable under search/scoring knobs: normalization boundary, numeric guards, and weighted least-squares variants of the ZM baseline?

## Why this is one experiment

The historical boundary, guard, and WLS bundles all perturb the symbolic-search or baseline-fitting protocol while asking whether the core step-2 winner remains stable. They are different knobs on the same robustness question.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_boundary_ablation/`
- `/Volumes/External2TB/emlexperiment/results/zipf_guard_ablation/`
- `/Volumes/External2TB/emlexperiment/results/zipf_wls_test/`
- Canonical 1a table for comparison against baseline step-2 winners.

## Inputs

- `results/zipf_boundary_ablation/summary.json`
- `results/zipf_guard_ablation/summary.json`
- `results/zipf_wls_test/summary.json`
- `experiments/1a_per_corpus_enriched_search/outputs/table1_per_corpus.csv`

## Outputs produced

### `outputs/robustness_ablation_by_corpus.csv`

One row per saved ablation setting.

Columns:
- `ablation_family`
- `slug`
- `corpus`
- `setting`
- `x_low`
- `exp_clamp`
- `value_abs_limit`
- `weighting`
- `token_count`
- `vocabulary_size`
- `zm_c`
- `zm_rmse`
- `step2_winner_expr`
- `step2_winner_math`
- `step2_rmse`
- `step2_delta_vs_zm`
- `step2_helpful`
- `canonical_step2_winner_expr`
- `matches_canonical_step2_winner`
- `notes`

### `outputs/aggregate_statistics.csv`

Header-only file. The v4 claim map assigns no aggregate rows to this experiment; manuscript robustness prose maps to the per-setting CSV rows.

### `outputs/manifest.json`

Machine-readable output inventory and audit flags.

## Manuscript lines fed

- v4 line `206` and v5.1 line `220`: robustness claims involving numeric guards, `x_low`, and weighted least-squares variants.

## Methods

This migration reads saved historical ablation summaries only. It does not rerun symbolic search or weighted fitting.

## AUDIT

The boundary and guard bundles save fresh step-2 winner expressions for each setting. Their scope is diagnostic subsets, not all 25 corpora.

The WLS bundle did not run and save a fresh symbolic-search winner for each weighted fit. The script applied the canonical IS correction `(x-1)-log(x)` after fitting each weighted ZM baseline and saved the resulting unweighted RMSE delta. Therefore, canonical migration can support the magnitude/effect of the fixed correction under WLS, but the manuscript phrase "did not alter winner identity on any corpus" should be read as a fixed-correction diagnostic unless a fresh WLS search is rerun.

## Canonical claim mapping

- Boundary rows map to `outputs/robustness_ablation_by_corpus.csv`, `ablation_family=boundary_x_low`.
- Guard rows map to `outputs/robustness_ablation_by_corpus.csv`, `ablation_family=numeric_guard`.
- WLS rows map to `outputs/robustness_ablation_by_corpus.csv`, `ablation_family=weighted_least_squares`.

## Rerun command

```bash
python3 experiments/1c_search_robustness_ablations/scripts/run_1c.py \
  --output-dir experiments/1c_search_robustness_ablations/outputs
```

## Verification mapping

- `robustness_ablation_by_corpus.csv` must have `30` rows: 9 boundary, 9 guard, and 12 WLS rows.
- `aggregate_statistics.csv` must be header-only unless the claim map gains explicit aggregate metrics for this experiment.
