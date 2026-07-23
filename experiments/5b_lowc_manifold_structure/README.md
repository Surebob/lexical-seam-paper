# Experiment 5b README

## Experiment id

`5b_lowc_manifold_structure`

## Research question

What is the geometry of the low-c manifold, and how does winner identity change when the scoring metric shifts from full-curve RMSE to head-focused RMSE?

## Why this is one experiment

The low-c manifold analysis and the phase-coordinate sweep are two views of the same question: whether the empirical low-c exp winner is a stable discrete phase or a near-degenerate head manifold whose named winner changes under metric emphasis.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_lowc_manifold_analysis`](/Volumes/External2TB/emlexperiment/results/zipf_lowc_manifold_analysis)
- [`/Volumes/External2TB/emlexperiment/results/zipf_phase_coordinate`](/Volumes/External2TB/emlexperiment/results/zipf_phase_coordinate)
- [`/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery`](/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery) for low-c smooth-vs-ZM top-100 modal-match rates

## Outputs produced

### `outputs/lowc_manifold_per_corpus.csv`

One row per low-c English corpus (`14` rows).

Key columns:
- `slug`
- `corpus`
- `single_zm_c`
- `full_winner`
- `top100_winner`
- `top200_winner`
- `cos_exp_xpow_top200`
- `r2_span_exp_xpow_top200`
- `r2_span_exp_is_top200`
- `xpow_minus_exp_full_rmse`
- `xpow_minus_exp_top100_rmse`

### `outputs/phase_coordinate_per_corpus.csv`

One row per English corpus (`25` rows), with phase-coordinate and lambda-sweep summaries.

### `outputs/aggregate_statistics.csv`

Claim-facing rows include:
- `lowc_median_cosine_exp_vs_xpow_top200`
- `lowc_median_span_r2_exp_xpow_top200`
- `lowc_median_span_r2_exp_is_top200`
- `lowc_full_rmse_exp_winner_count`
- `lowc_top100_xpow_winner_count`
- `lowc_top200_xpow_winner_count`
- `lowc_median_delta_xpow_minus_exp_full_rmse`
- `lowc_median_delta_xpow_minus_exp_top100_rmse`
- `lowc_smooth_modal_vs_empirical_top100_match_rate`
- `lowc_zm_modal_vs_empirical_top100_match_rate`

## Manuscript claims fed

Mapped claims include:
- median cosine `0.965` between `exp(x-1)-x` and `x^x-sqrt(x)` over the low-c head
- full-RMSE exp winner count `14/14`
- top-100 xpow winner count `11/14`
- smooth synthetic modal top-100 match rate `10/14 = 71.4%`
- single-ZM control modal top-100 match rate `8/14 = 57.1%`

## Methods

This migration reads saved historical summaries and recomputes only deterministic aggregations from saved rows. It does not rerun symbolic search.

Inference documentation: no confidence intervals or p-values are computed. The outputs are descriptive counts, medians, cosine similarities, R² values, and Pearson correlations saved by the historical bundles.

## Rerun command

```bash
python3 experiments/5b_lowc_manifold_structure/scripts/run_5b.py
```

## Verification mapping

- `lowc_manifold_per_corpus.csv` must contain `14` rows.
- `phase_coordinate_per_corpus.csv` must contain `25` rows.
- `lowc_smooth_modal_vs_empirical_top100_match_rate` must equal `10/14`.
- `lowc_zm_modal_vs_empirical_top100_match_rate` must equal `8/14`.

## AUDIT

This experiment is primarily empirical residual geometry and metric-sensitivity analysis. It depends on historical simulation-recovery output only for the modal-match comparison. If the simulation-recovery source is later rerun under a strict decoupled-erf canon, the modal-match aggregate rows should be revisited.

