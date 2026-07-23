# Experiment 3a README

## Experiment id

`3a_smooth_two_regime_fits`

## Research question

Across the 25 English corpora, what are the canonical bounded 8-parameter smooth two-regime fits, and how do those fits change when the parameter box is relaxed?

## Why this is one experiment

The constrained and relaxed runs are two reporting views on the same fitted model family applied to the same 25-corpus English sample. The constrained fit is the canonical smooth-model result; the relaxed rerun is its bounds-robustness check. They belong together because the manuscript compares them directly.

## Upstream dependencies

- Shared English corpus registry in [`/Volumes/External2TB/emlexperiment/zipf_analysis_common.py`](/Volumes/External2TB/emlexperiment/zipf_analysis_common.py)
- Smooth two-regime base implementation in [`/Volumes/External2TB/emlexperiment/zipf_correct_model_reranked.py`](/Volumes/External2TB/emlexperiment/zipf_correct_model_reranked.py)
- Piecewise comparison model in [`/Volumes/External2TB/emlexperiment/zipf_correct_model.py`](/Volumes/External2TB/emlexperiment/zipf_correct_model.py)
- Historical comparison bundles:
  - [`/Volumes/External2TB/emlexperiment/results/zipf_correct_model_reranked_all/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_correct_model_reranked_all/summary.json)
  - [`/Volumes/External2TB/emlexperiment/results/zipf_reranked_model_all_corpora_relaxed/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_reranked_model_all_corpora_relaxed/summary.json)

## Inputs

- The 25 canonical English corpora
- Canonical constrained smooth-model protocol:
  - `8` parameters
  - `100` random starts
  - `max_nfev = 12000`
  - bounds:
    - `a1, a2 in [5, 50]`
    - `b1, b2 in [0.5, 3.0]`
    - `c1, c2 in [0, 1000]`
    - `k in [20, 2000]`
    - `w in [0.1, 3.0]`
- Reranked tail coordinate shift:
  - `r - floor(0.8 k)` approximation documented via `RERANK_SHIFT_FRACTION = 0.8`
- Relaxed-bounds robustness protocol:
  - same optimizer settings
  - relaxed bounds:
    - `a1, a2 in [-100, 100]`
    - `b1, b2 in [0.5, 5.0]`
    - `c1, c2 in [0, 5000]`
    - `k in [20, 2000]`
    - `w in [0.1, 5.0]`

## Outputs produced

### `outputs/smooth_fit_per_corpus.csv`
One row per corpus for the canonical constrained 8-parameter smooth fit.

Columns:
- `slug`
- `corpus`
- `token_count`
- `vocab_size`
- `single_zm_rmse`
- `piecewise_rmse`
- `step2_rmse`
- `step2_expression`
- `reranked_rmse`
- `beats_single_zm`
- `beats_piecewise`
- `beats_step2`
- `best_start_index`
- `best_nfev`
- `a1`
- `b1`
- `c1`
- `a2`
- `b2`
- `c2`
- `k`
- `w`
- `transition_fraction`
- `tried_count`

### `outputs/bounds_robustness_per_corpus.csv`
One row per corpus comparing the canonical constrained fit to the relaxed-bounds rerun.

Columns:
- `slug`
- `corpus`
- `original_rmse`
- `relaxed_rmse`
- `rmse_delta_vs_original`
- `original_transition_fraction`
- `relaxed_transition_fraction`
- `transition_fraction_delta`
- `original_k`
- `relaxed_k`
- `original_w`
- `relaxed_w`
- `best_start_index`
- `best_nfev`
- `tried_count`
- `relaxed_a1`
- `relaxed_b1`
- `relaxed_c1`
- `relaxed_a2`
- `relaxed_b2`
- `relaxed_c2`

### `outputs/aggregate_statistics.csv`
Single-cell aggregates used by manuscript prose.

Required rows:
- `english_corpus_count`
- `smooth_beats_single_zm_count`
- `smooth_beats_piecewise_count`
- `smooth_beats_step2_count`
- `mean_k`
- `sd_k`
- `min_k`
- `max_k`
- `mean_w`
- `sd_w`
- `min_w`
- `max_w`
- `mean_transition_fraction`
- `sd_transition_fraction`
- `min_transition_fraction`
- `max_transition_fraction`
- `w_lt_0_2_count`
- `w_lt_0_2_corpora`
- `relaxed_bounds_improves_rmse_count`
- `relaxed_bounds_unchanged_count`
- `relaxed_bounds_degraded_count`
- `relaxed_rmse_delta_mean`
- `relaxed_rmse_delta_median`
- `relaxed_rmse_delta_min`
- `relaxed_rmse_delta_max`
- `relaxed_rmse_delta_sd`
- `relaxed_transition_fraction_delta_mean`
- `relaxed_transition_fraction_delta_median`

### Historical diff outputs
- `outputs/historical_smooth_diff.csv`
- `outputs/historical_relaxed_diff.csv`
- `outputs/historical_diff_summary.csv`
- `outputs/manifest.json`

## Manuscript lines fed

Primary:
- lines `230–257` (Sections 3.5 and 3.6 smooth-fit side)
- line `245` (`w < 0.2` corpora and mean/SD of `w`)
- line `538` (`25/25` relaxed-bounds improvement count)

Secondary:
- Supplementary parameter table for the 25 constrained smooth fits
- English-side totals later joined by `J1`

## Methods

This experiment reruns the canonical bounded 8-parameter smooth two-regime model on the full 25-corpus English sample using the original reranked implementation, then reruns the same family under relaxed bounds as a robustness check. The constrained fit is the canonical scientific result; the relaxed fit is a sensitivity analysis. Both are compared back to their historical saved summaries to verify reproducibility. No statistical inference is performed in this experiment.

## Canonical claim mapping

Examples:
- line `245` mean width and sharp-transition corpus list -> `outputs/aggregate_statistics.csv`, rows `mean_w`, `sd_w`, `w_lt_0_2_corpora`
- line `538` relaxed-bounds `25/25` claim -> `outputs/aggregate_statistics.csv`, row `relaxed_bounds_improves_rmse_count`
- per-corpus smooth parameters for supplement -> `outputs/smooth_fit_per_corpus.csv`

## Rerun command

```bash
python3 experiments/3a_smooth_two_regime_fits/scripts/run_3a.py \
  --output-dir experiments/3a_smooth_two_regime_fits/outputs
```

## Verification mapping

Post-run verification checks:
- `smooth_fit_per_corpus.csv` row count must be `25`
- `bounds_robustness_per_corpus.csv` row count must be `25`
- `smooth_beats_single_zm_count`, `smooth_beats_piecewise_count`, and `smooth_beats_step2_count` must each equal `25`
- `w_lt_0_2_count` must equal the number of names listed in `w_lt_0_2_corpora`
- `historical_diff_summary.csv` must report max numeric drift `<= 1e-4` versus the two subsumed historical bundles

## Upstream / downstream provenance

Upstream:
- English corpus preprocessing
- piecewise ZM fit helper
- reranked smooth-model implementation

Downstream:
- `3b` comparison context
- `J1` English corpus rollups
- manuscript Sections 3.5 and 3.6

## Notes on non-goals

This experiment does **not** include:
- BIC family competition against MOE and piecewise variants (`3b`)
- statistical scaling inference for `k_stat` (`3c`)
- POS crossover scaling (`3d`)

Those are separate experiments because they answer different research questions.
