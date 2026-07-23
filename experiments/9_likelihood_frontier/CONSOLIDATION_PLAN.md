# Experiment 9 Consolidation Plan

## Experiment id

`9_likelihood_frontier`

## Research question

How far can likelihood be pushed beyond the mechanism-first soft-k Seam-Mandelbrot PMF, and what mechanism cost appears when the architecture is relaxed?

## Source bundles inspected

### `results/zipf_hybrid_headtail_splitfit`

- `hybrid_headtail_table.csv`: 25 rows.
- Columns: `slug`, `name`, `lambda`, `k`, `w`, `top100_hybrid`, `top100_moe`, `top100_softk`, `full_hybrid`, `full_moe`, `full_softk`, `full_step2_help`, `full_step2_gain`.
- `summary.json`: per-corpus held-out windows and aggregate cutoffs for top50/top100/top200/top500/top1000/full.
- `report.md`: states the main hybrid result.
- Full-cutoff headline values from this source:
- hybrid beats MOE: `21/25`
- hybrid beats soft-k: `22/25`
- median hybrid minus soft-k: `-0.001677944396615949`
- hybrid step-2 help count: `9/25`

### `results/zipf_hybrid_mechpenalty`

- `hybrid_mechpenalty_table.csv`: 125 rows = 25 corpora x 5 penalty strengths.
- Columns: `slug`, `name`, `lambda_mech`, `test_avg_nll`, `train_step2_gain`, `heldout_step2_gain`, `heldout_step2_help`, `k`, `w`, `alpha_tail`, `beta_tail`.
- `summary.json`: aggregate by `lambda_mech`.
- `report.md`: states no penalty strength satisfied the intended decision rule.
- Penalty grid: `{0.01, 0.1, 1, 10, 100}`.
- Held-out step-2 help count is `25/25` for every penalty strength.

### `results/zipf_hybrid_vs_softk_analysis`

- `hybrid_vs_softk_table.csv`: 25 rows.
- Columns: `slug`, `name`, `V`, `token_count`, `zm_c`, `hybrid_full_nll`, `softk_full_nll`, `delta_hybrid_minus_softk`, `winner`, `head_bregman_winner`, `anthology_like`, `softk_lambda`, `softk_transition_fraction`, `hybrid_transition_fraction`, `hybrid_step2_help`, `softk_step2_help`, `moe_full_nll`, `zm_full_nll`.
- `summary.json`: correlations and by-winner summaries.
- `report.md`: states hybrid wins on `16/25`, soft-k wins on `9/25`.

### `results/zipf_fully_ours_variants`

- `fully_ours_variants_table.csv`: 50 rows = 25 corpora x 2 variants.
- Columns: `slug`, `name`, `variant`, `test_avg_nll`, `bic`, `rmse`, `step2_gain`, `step2_help`.
- `summary.json`: aggregate rows for `nested` and `three_regime`.
- `report.md`: states nested and three-regime headline results.
- Nested headline values:
- beats MOE: `21/25`
- beats soft-k: `18/25`
- beats ZM: `16/25`
- four-way winner count: `14/25`
- step-2 help count: `11/25`
- Three-regime headline values:
- beats soft-k: `12/25`
- median minus soft-k: `0.00010526860518478287`

## Source discrepancy: hybrid-vs-soft-k convention

`zipf_hybrid_headtail_splitfit` and `zipf_hybrid_vs_softk_analysis` agree on the hybrid full NLL for almost all corpora, but they disagree on the soft-k comparator used for the hybrid-vs-soft-k decision. The discrepancy changes the headline:

- `zipf_hybrid_headtail_splitfit`: hybrid beats soft-k on `22/25`, median hybrid minus soft-k `-0.001677944396615949`.
- `zipf_hybrid_vs_softk_analysis`: hybrid wins on `16/25`, soft-k wins on `9/25`.

Spot check:

- Aesop's Fables: headtail soft-k `6.22074737680541`, analysis soft-k `6.2105248148329`; same hybrid value, different winner.
- Moby Dick: headtail soft-k `6.948129203324875`, analysis soft-k `6.944063683425718`; same hybrid value, different winner.
- Ulysses: headtail soft-k `7.353701494174882`, analysis soft-k `7.335024268866188`; same hybrid value, different margin.

Mechanism confirmed by row-level source tracing: this is a comparator-convention divergence, not a different hybrid model. The analysis bundle uses the legacy pre-splitfit derivative branch soft-k comparator. Its `softk_full_nll` values match `results/zipf_hybrid_headtail`, `results/zipf_angle1_head_windows`, and `results/zipf_angle3_asymmetric_gate` exactly on all `25/25` rows, and match `0/25` rows for v4 verification, `softk_splitfit`, `softkw` repeated soft-k, or the patched splitfit derivative bundles.

The producer script `zipf_hybrid_vs_softk_analysis.py` loads `results/zipf_hybrid_headtail/summary.json`; its `softk_full_nll` column is taken from `h["heldout"]["full"]["avg_nll"]["softk"]` in that legacy hybrid-headtail row. The separate `zipf_seam_mandelbrot_softk` summary is loaded for metadata fields such as lambda and transition fraction, not for the soft-k NLL comparator column.

See `SOFTK_COMPARATOR_DIAGNOSTIC.md` and `softk_comparator_source_diagnostic.csv` in this experiment directory.

## Canonical-source decision proposed

Use `zipf_hybrid_headtail_splitfit` as the canonical source for manuscript-facing hybrid-vs-baseline metrics because:

- It is the primary hybrid fit bundle.
- Its metadata explicitly says the soft-k comparator uses split-fit parameters.
- Its headline values match manuscript v5.1 and the claim map for the hybrid Pareto-frontier paragraph.
- It aligns with the saved report for the hybrid model itself.

Preserve `zipf_hybrid_vs_softk_analysis` as a diagnostic/provenance sidecar, not as the canonical source for hybrid-vs-soft-k counts, because it uses the legacy pre-splitfit comparator convention. This is legitimate historical research output, but it is not the current post-cleanup canonical comparator convention.

This follows the general principle approved for methodological-path differences: preserve both sources with clear diagnostics. Here, however, the manuscript-facing primary hybrid metric should not silently switch to the analysis bundle's alternate comparator convention.

## Planned canonical outputs

### `outputs/hybrid_headtail_per_corpus.csv`

Source: `zipf_hybrid_headtail_splitfit/hybrid_headtail_table.csv` and `summary.json`.

Rows: 25.

Canonical columns:

- corpus identifiers and counts
- selected `lambda`, `k`, `w`
- `top100_hybrid`, `top100_moe`, `top100_softk`
- `full_hybrid`, `full_moe`, `full_softk`
- `hybrid_minus_moe`
- `hybrid_minus_softk`
- `full_step2_help`
- `full_step2_gain`

Feeds:

- `hybrid_beats_moe_count`
- `hybrid_beats_softk_count`
- `median_hybrid_minus_softk`
- `hybrid_step2_help_count`

### `outputs/fully_ours_variants.csv`

Source: `zipf_fully_ours_variants/fully_ours_variants_table.csv`.

Rows: 50.

Canonical columns:

- corpus identifiers
- `variant` (`nested`, `three_regime`)
- held-out NLL
- BIC
- RMSE
- step-2 gain/help

Feeds:

- `nested_beats_moe_count`
- `nested_beats_softk_count`
- `nested_beats_zm_count`
- `nested_fourway_winner_count`
- `nested_step2_help_count`
- three-regime beats-softk and median-minus-softk context

### `outputs/mechanism_penalty_sweep.csv`

Source: `zipf_hybrid_mechpenalty/hybrid_mechpenalty_table.csv`.

Rows: 125.

Canonical columns:

- corpus identifiers
- `lambda_mech`
- held-out NLL
- train step-2 gain
- held-out step-2 gain/help
- fitted `k`, `w`, `alpha_tail`, `beta_tail`

Feeds:

- `mechanism_penalty_grid`
- `mechanism_penalty_step2_help_count_all_lambdas`
- `mechanism_penalty_frontier_result`

### `outputs/hybrid_vs_softk_diagnostic.csv`

Source: `zipf_hybrid_vs_softk_analysis/hybrid_vs_softk_table.csv` plus comparison against `hybrid_headtail_per_corpus.csv`.

Rows: 25.

Purpose: preserve the alternate comparator analysis and explicitly document the convention divergence.

Columns:

- source hybrid NLL
- analysis soft-k NLL
- headtail/canonical soft-k NLL
- analysis delta
- canonical delta
- analysis winner
- canonical winner
- soft-k comparator delta

Feeds: no Paper 1 claim directly. This is provenance/audit support.

### `outputs/hybrid_structure_summary.csv`

Source: aggregate dictionaries from all four bundles plus protocol constants.

Rows: one row per named model or diagnostic:

- `hybrid_headtail`
- `nested_seam`
- `three_regime`
- `mechanism_penalty`
- `nested_protocol_constants`

This file is the primary claim-map target for the Pareto-frontier paragraph.

### `outputs/aggregate_statistics.csv`

Claim-facing key/value rows including:

- `hybrid_beats_moe_count = 21`
- `hybrid_beats_softk_count = 22`
- `median_hybrid_minus_softk = -0.001677944396615949`
- `hybrid_step2_help_count = 9`
- `nested_beats_moe_count = 21`
- `nested_beats_softk_count = 18`
- `nested_beats_zm_count = 16`
- `nested_fourway_winner_count = 14`
- `nested_step2_help_count = 11`
- `three_regime_beats_softk_count = 12`
- `three_regime_median_minus_softk = 0.00010526860518478287`
- `mechanism_penalty_grid = [0.01, 0.1, 1, 10, 100]`
- `mechanism_penalty_step2_help_count_all_lambdas = 25/25 for every lambda`
- `mechanism_penalty_frontier_result = no lambda satisfies help <= 4/25 and beats MOE >= 18/25`
- `nested_seam_free_parameter_count = 8`
- `nested_mini_seam_floor = 5`
- `nested_mini_seam_logV_multiplier = 0.5`
- `nested_mini_seam_width = 0.22`

## Precedence rules

1. Hybrid model fit values come from `zipf_hybrid_headtail_splitfit`.
2. Hybrid-vs-soft-k manuscript-facing counts use the `full_softk` comparator stored in `zipf_hybrid_headtail_splitfit`, not the alternate `zipf_hybrid_vs_softk_analysis` comparator.
3. `zipf_hybrid_vs_softk_analysis` is preserved as diagnostic only and must not overwrite the headline hybrid counts.
4. Nested and three-regime metrics come from `zipf_fully_ours_variants`.
5. Mechanism-penalty metrics come from `zipf_hybrid_mechpenalty`.
6. Protocol constants for the nested seam architecture come from manuscript Section 2.19 / claim map protocol rows and will be emitted as constants in `aggregate_statistics.csv`.

## Union vs select rules

- Union all source rows into separate canonical sidecar outputs.
- Select only one primary source for each manuscript-facing metric.
- If future Paper 2 revisits the PMF arc, it should rerun or re-document the hybrid-vs-soft-k comparator convention rather than inheriting this selection without review.

## README AUDIT entries to include after execution

- The hybrid-vs-soft-k analysis bundle uses a different soft-k comparator convention and yields `16/25` hybrid wins rather than the canonical `22/25`; preserved as diagnostic, not used for Paper 1 headline counts.
- The PMF arc is queued for Paper 2; this experiment's current canonical role is research-record preservation and v5.1 provenance support.
- Mechanism-penalty results are negative: held-out step-2 help remains `25/25` for all penalty strengths.
