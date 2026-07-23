# Experiment 3b Consolidation Plan

## Stop Point

This document is a consolidation plan only. Experiment `3b` has not been migrated or executed. Approval is required before creating canonical outputs.

## Research Question

Which model family wins the BIC comparison across the 25 English corpora, and how large is the fixed-`k = sqrt(V)` smooth-model cost relative to the free-`k` smooth model?

## Source Bundles

### `results/zipf_continuous_piecewise`

Files:
- `continuous_piecewise_table.csv`
- `summary.json`
- `report.md`

Saved content:
- 25 rows.
- Continuous-piecewise RMSE and chosen break `k`.
- Six-family rank-BIC columns:
  - `zipf_bic`
  - `moezipf_bic`
  - `piecewise_k500_bic`
  - `continuous_piecewise_bic`
  - `reranked_7param_sqrtv_bic`
  - `reranked_8param_bic`
- `best_model`.
- `winner_counts` for the six-family comparison.

Canonical use:
- Primary source for `outputs/table2_model_family.csv`.
- Primary source for `winner_family` and BIC winner counts.
- Primary source for hard/continuous piecewise comparison rows.

### `results/zipf_moezipf_comparison`

Files:
- `moezipf_table.csv`
- `summary.json`
- `report.md`

Saved content:
- 25 rows.
- Token and vocabulary counts.
- Zipf and MOEZipf likelihood/AIC/BIC fields.
- Rank-BIC fields:
  - `zipf_rank_bic`
  - `moe_rank_bic`
  - `sqrt_v_rank_bic`
  - `reranked_rank_bic`
- Rank-RMSE fields:
  - `zipf_rmse`
  - `moe_rmse`
  - `sqrt_v_rmse`
  - `reranked_rmse`
- `rank_bic_winner`.
- `moby_paper_check`.

Canonical use:
- Source for `outputs/model_family_rmse_per_corpus.csv`.
- Source for the Moby Dick RMSE claim: MOEZipf RMSE `0.050` versus smooth free-`k` RMSE `0.119`.
- Cross-check source for `zipf_bic`, `moezipf_bic`, `reranked_7param_sqrtv_bic`, and `reranked_8param_bic` in `zipf_continuous_piecewise`.

Precedence:
- For BIC table values shared with `zipf_continuous_piecewise`, use `zipf_continuous_piecewise` as the primary table source and use this bundle only as a drift check.

### `results/zipf_bic_comparison`

Files:
- `bic_table.csv`
- `summary.json`
- `report.md`

Saved content:
- 25 rows.
- BIC columns:
  - `single_zm_bic`
  - `piecewise_k500_bic`
  - `reranked_8param_bic`
  - `reranked_7param_sqrtv_bic`
  - `reranked_6param_sqrtv_w12_bic`
- `best_model`.
- `winner_counts`.
- `param_counts`.

Canonical use:
- Cross-check for `piecewise_k500_bic`, `reranked_8param_bic`, and `reranked_7param_sqrtv_bic`.
- Archive/provenance source for the exploratory 6-parameter fixed-width variant.

Important conflict:
- `bic_table.csv`'s `single_zm_bic` does **not** match `zipf_continuous_piecewise`'s `zipf_bic` or `zipf_moezipf_comparison`'s `zipf_rank_bic`; all 25 rows differ by more than `1e-4`, with max absolute difference `38590.100157910994`.
- The shared piecewise and smooth BIC fields match exactly (`max drift = 0.0`).
- Supplementary diagnosis: see `BIC_CONVENTION_DIAGNOSTIC.md`. The divergence is a baseline-model mismatch, not a BIC-formula mismatch: `zipf_bic_comparison/single_zm_bic` is a 3-parameter single-ZM rank-BIC, while `zipf_continuous_piecewise/zipf_bic` is the 1-parameter pure-Zipf rank-BIC imported from `zipf_moezipf_comparison`.

Precedence:
- User decision after `BIC_CONVENTION_DIAGNOSTIC.md`: use true 3-parameter `single_zm_bic` from `zipf_bic_comparison` as the canonical baseline.
- Drop `zipf_rank_bic` from canonical outputs because it is a 1-parameter pure-Zipf value and is not the paper's baseline.
- Preserve the v4 pure-Zipf value only as a legacy audit field, not as a canonical model-family column.

### `results/zipf_bic_landscape`

Files:
- `bic_winner_counts.png`

Saved content:
- No CSV or JSON data was found.
- The bundle has no machine-readable source table.

Canonical use:
- None for canonical numeric outputs.

Blocked/missing:
- Any claim-map metric requiring a machine-readable BIC-landscape table is not reconstructible from this bundle without rerunning or manual image extraction.
- No such specific metric was identified in the v4 map beyond the model-family BIC comparison, which can be migrated from other bundles.

### `results/zipf_sqrt_v_all_corpora`

Files:
- `summary.json`
- `report.md`

Saved content:
- 25 rows in JSON.
- Per-corpus `sqrt_v_rmse`, free-`k` smooth RMSE, delta relative to free-`k`, fixed `k`, and parameters.
- Aggregate `rmse_gap_to_8param` with `mean`, `median`, `min`, `max`, and `std`.

Canonical use:
- Source for `median_rmse_cost_sqrtv_vs_free_k`.
- Optional per-corpus fixed-`k` RMSE cross-check for `model_family_rmse_per_corpus.csv`.

## Planned Canonical Outputs

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

Source mapping:
- `single_zm_bic` <- `zipf_bic_comparison.single_zm_bic`
- `moe_bic` <- `zipf_continuous_piecewise.moezipf_bic`
- `hard_piecewise_bic` <- `zipf_continuous_piecewise.piecewise_k500_bic`
- `continuous_piecewise_bic` <- `zipf_continuous_piecewise.continuous_piecewise_bic`
- `smooth_ksqrtv_bic` <- `zipf_continuous_piecewise.reranked_7param_sqrtv_bic`
- `smooth_freek_bic` <- `zipf_continuous_piecewise.reranked_8param_bic`
- `winner_family` <- normalized `zipf_continuous_piecewise.best_model`

### `outputs/model_family_rmse_per_corpus.csv`

One row per English corpus.

Columns:
- `slug`
- `corpus`
- `vocabulary_size`
- `zipf_rank_rmse`
- `moezipf_rank_rmse`
- `smooth_ksqrtv_rank_rmse`
- `smooth_free_k_rank_rmse`

Source mapping:
- Primary source: `zipf_moezipf_comparison/moezipf_table.csv`.
- `smooth_free_k_rank_rmse` <- `reranked_rmse`.
- `smooth_ksqrtv_rank_rmse` <- `sqrt_v_rmse`.

### `outputs/bic_source_crosscheck.csv`

One row per corpus and shared field.

Columns:
- `slug`
- `field`
- `primary_source`
- `secondary_source`
- `primary_value`
- `secondary_value`
- `abs_diff`
- `status`

Purpose:
- Preserve the cross-source BIC convention audit.
- Mark `single_zm_bic` as a convention conflict.
- Confirm shared piecewise/smooth BIC values match to `<= 1e-4`.

### `outputs/aggregate_statistics.csv`

Rows required by v4 claim map:
- `piecewise_loses_to_smooth_on_all_25`
- `moezipf_bic_winner_count`
- `winner_count_moezipf`
- `winner_count_smooth_family_total`
- `median_rmse_cost_sqrtv_vs_free_k`

Potential protocol constants in `source_config.py`:
- `piecewise_breakpoint_rank = 500`
- `fixed_k_law = sqrt(V)`
- `smooth_model_parameter_count = 7`

### `outputs/manifest.json`

Machine-readable inventory:
- source bundles
- output schemas
- row counts
- precedence rules
- conflict flags
- blocked items, if any

## Precedence Rules

1. For the canonical six-family BIC table, use true 3-parameter `single_zm_bic` from `zipf_bic_comparison`.
2. Use `zipf_moezipf_comparison` to cross-check shared rank-BIC values and to provide rank-RMSE columns.
3. Use `zipf_bic_comparison` for the canonical single-ZM baseline and shared smooth/piecewise BIC cross-checks.
4. Use `zipf_sqrt_v_all_corpora` for the fixed-`k` RMSE-cost aggregate.
5. Do not use `zipf_bic_landscape` for any numeric canonical value because it has only a PNG.

## Union vs Select Rules

Select:
- The canonical table keeps exactly the six model-family columns named in the v4 claim map.
- The canonical RMSE table keeps only columns needed for manuscript claims and BIC interpretation.

Archive / do not merge into canonical table:
- `reranked_6param_sqrtv_w12_bic` from `zipf_bic_comparison`.
- `moby_paper_check` from `zipf_moezipf_comparison`, unless a future manuscript claim maps to it.
- `bic_winner_counts.png` from `zipf_bic_landscape`.

## Planned Audit Notes

The README should include an AUDIT section with:

- `zipf_bic_comparison/single_zm_bic` conflicts with the rank-BIC convention used by `zipf_continuous_piecewise` and `zipf_moezipf_comparison`.
- The historical v4 Table 2 builder used `zipf_continuous_piecewise.zipf_bic`, which is pure-Zipf rank-BIC, while labeling the rendered column `ZM`. This is a naming/provenance mismatch documented in `BIC_CONVENTION_DIAGNOSTIC.md`; canonical 3b corrects it.
- `zipf_bic_landscape` is not empty, but contains only a PNG and no machine-readable table.
- Model-canon note: Adjustment 1 updates `3a`, `3c`, `4`, and `5c` to decoupled-erf; this `3b` migration remains based on saved historical BIC comparison bundles unless the user explicitly requests decoupled-erf BIC recomputation or remapping.

## Approval Needed

The user approved the true 3-parameter single-ZM convention. Migration should proceed with `single_zm_bic` as canonical and flag any manuscript values carried over from the v4 pure-Zipf table for revision.
