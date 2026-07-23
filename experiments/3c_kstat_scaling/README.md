# Experiment 3c README

## Experiment id

`3c_kstat_scaling`

## Research question

How does the transition centre inferred from the canonical smooth two-regime fit scale with vocabulary size `V` across the 25 English corpora?

## Why this is one experiment

This is the same statistical-scaling operation applied to the same 25-corpus English sample after the canonical smooth-fit experiment has already emitted one `k_stat` value per corpus. Historically the saved bundle also carried a POS comparison, but under the consolidated design the `k_stat` scaling result stays here and the cross-experiment `k_stat` vs `k_POS` comparison moves to `J2`.

## Upstream dependencies

- Canonical smooth-fit output from [`/Volumes/External2TB/emlexperiment/experiments/3a_smooth_two_regime_fits/outputs/smooth_fit_per_corpus.csv`](/Volumes/External2TB/emlexperiment/experiments/3a_smooth_two_regime_fits/outputs/smooth_fit_per_corpus.csv)
- Historical comparison bundle:
  - [`/Volumes/External2TB/emlexperiment/results/zipf_kstat_scaling/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_kstat_scaling/summary.json)
  - [`/Volumes/External2TB/emlexperiment/results/zipf_kstat_scaling/report.md`](/Volumes/External2TB/emlexperiment/results/zipf_kstat_scaling/report.md)

## Inputs

- The `25` canonical English corpora as represented by `3a`'s per-corpus smooth-fit output
- For each corpus:
  - `V` = vocabulary size
  - `k_stat` = smooth-fit transition centre from experiment `3a`
- Forced scaling model:
  - `log(k_stat) = alpha * log(V)`
- Reference null value for inferential checks:
  - `alpha = 0.5`

## Outputs produced

### `outputs/kstat_scaling_points.csv`
One row per English corpus for the statistical scaling dataset.

Columns:
- `slug`
- `corpus`
- `vocabulary_size`
- `k_stat`
- `log_vocabulary_size`
- `log_k_stat`
- `alpha_stat_per_corpus`

### `outputs/aggregate_statistics.csv`
Single-cell aggregates used by manuscript prose.

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Required rows:
- `english_corpus_count`
- `forced_alpha`
- `forced_alpha_se`
- `forced_alpha_ci_low`
- `forced_alpha_ci_high`
- `forced_alpha_df`
- `mean_alpha_per_corpus`
- `median_alpha_per_corpus`
- `mean_alpha_ci_low`
- `mean_alpha_ci_high`
- `alpha_stat_vs_half_t_statistic`
- `alpha_stat_vs_half_pvalue`
- `forced_alpha_vs_half_pvalue`
- `alpha_stat_vs_half_df`

### Historical diff outputs
- `outputs/historical_point_diff.csv`
- `outputs/historical_aggregate_diff.csv`
- `outputs/historical_diff_summary.csv`
- `outputs/manifest.json`

## Manuscript lines fed

Primary:
- line `13` (`V^0.521` statistical-scaling sentence, via direct `3c` outputs)
- lines `255–257` (Section 3.6 `k_stat` scaling claim)
- line `516` (conclusion restatement of the free-fit statistical exponent)
- line `661` (appendix verification checklist)

Downstream only:
- `J2` uses `3c` together with `3d` to emit the explicit `k_stat` vs `k_POS` comparison claims.

## Methods

This experiment constructs the statistical-scaling dataset by reading the canonical `k_stat` values already emitted by experiment `3a`, then fitting the forced model `log(k_stat) = alpha * log(V)` with `scipy.optimize.curve_fit`.

The reported **95% confidence interval for the forced-fit alpha** is computed from the fitted parameter covariance returned by `curve_fit`: the standard error is `sqrt(pcov[0,0])`, and the interval is `alpha ± t_(0.975, df) * SE` with `df = n - 1 = 24`.

The reported **p-value against 0.5** is **not** a test on the forced-fit alpha parameter itself. It is the two-sided one-sample Student-`t` test used in the historical bundle: `scipy.stats.ttest_1samp` applied to the per-corpus values `alpha_stat_per_corpus = log(k_stat) / log(V)` with null mean `0.5`. The corresponding `t` statistic and degrees of freedom are emitted alongside the p-value.

## Canonical claim mapping

Examples:
- line `255` / `661` alpha `0.5214` -> `outputs/aggregate_statistics.csv`, row `forced_alpha`
- line `255` / `661` CI `[0.4915, 0.5513]` -> `outputs/aggregate_statistics.csv`, rows `forced_alpha_ci_low`, `forced_alpha_ci_high`
- line `255` / `661` `p = 0.142` -> `outputs/aggregate_statistics.csv`, row `forced_alpha_vs_half_pvalue`

## Rerun command

```bash
python3 experiments/3c_kstat_scaling/scripts/run_3c.py \
  --output-dir experiments/3c_kstat_scaling/outputs
```

## Verification mapping

Post-run verification checks:
- `kstat_scaling_points.csv` row count must be `25`
- every corpus present in `3a` must appear exactly once here
- `forced_alpha`, `forced_alpha_ci_low`, `forced_alpha_ci_high`, and `forced_alpha_vs_half_pvalue` must match the historical `zipf_kstat_scaling` bundle to within `1e-4`
- `historical_diff_summary.csv` must report max numeric drift `<= 1e-4` on the overlapping `k_stat`-scaling outputs
- any old POS-comparison content from the historical bundle must be documented as deferred to `J2`, not silently recomputed here

## Upstream / downstream provenance

Upstream:
- `3a_smooth_two_regime_fits`

Downstream:
- `J2_scaling_comparison_summary`
- manuscript Section 3.6 `k_stat` claims

## Notes on non-goals

This experiment does **not** include:
- POS crossover scaling (`3d`)
- the explicit `k_stat` vs `k_POS` comparison (`J2`)

Those are separate outputs under the consolidated design because the comparison is a genuine cross-experiment join.
