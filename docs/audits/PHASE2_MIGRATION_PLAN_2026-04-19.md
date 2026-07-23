# Phase 2 Migration Plan — 2026-04-19

## Scope and Stop Point

This is a plan document only. No experiments have been migrated, no historical `results/zipf_*` bundles have been modified, and no new experiment has been run.

Packaging is paused. The prior dry-run manifest remains as a record, but this plan supersedes it as the next operational document.

Planning inputs:
- `results/CONSOLIDATION_PROPOSAL_v3.md`
- `results/MANUSCRIPT_CLAIM_TO_CSV_MAP_v4.md`
- Current canonical completed experiments under `experiments/`: `1a`, `3a`, `3c`
- Historical source bundles under `results/zipf_*`

Target after approved migration:
- `25/25` v3 experiments complete in canonical `experiments/` structure.
- `2/2` joined-output producers complete in canonical `joins/` structure.
- Additional `1a` threshold-characterization output complete.
- New `11_aligned_coordinate_robustness` experiment complete.

Important model-canon warning:
- This plan migrates the Phase 2 v3/v4 claim-map contract. It does **not** automatically resolve the later manuscript v5.1 migration to the decoupled 9-parameter erf-gate model. The current canonical `3a` and `3c` outputs still represent the older coupled 8-parameter logistic model; S2 v3 decoupled-erf outputs currently live outside the 25-experiment contract.

## Global Migration Method

For each experiment:
1. Create canonical directory `experiments/<id>_<slug>/`.
2. Read source bundle(s) from `results/zipf_*` read-only.
3. Emit exactly the v4 claim-map-facing canonical outputs for that experiment.
4. Emit `aggregate_statistics.csv` with the `metric_name` rows cited by `MANUSCRIPT_CLAIM_TO_CSV_MAP_v4.md`, no speculative extra rows.
5. Emit a standardized `README.md` following the existing `1a`, `3a`, `3c` schema.
6. Emit `source_config.py`.
7. Emit `scripts/run_<id>.py` or `run_experiment.py` as a thin, rerunnable migration script reading historical bundles and regenerating canonical outputs.
8. Emit `outputs/manifest.json` describing files, schemas, row counts, source bundles, and any drift checks.
9. If any required value is missing from historical outputs and cannot be reconstructed without rerunning expensive computation, emit `BLOCKED.md` with exact missing fields and stop for review.

Numerical verification rule:
- Any value directly cited in the manuscript or claim map should match the historical saved value or manuscript cited value to within `1e-4`; larger differences become a drift investigation, not a silent migration.

## Existing Complete Experiments To Preserve

### 1a. Per-corpus English symbolic search

Current status: complete, but will be extended.

Current canonical outputs:
- `outputs/table1_per_corpus.csv`
- `outputs/table1_step2_beam_top10.csv`
- `outputs/aggregate_statistics.csv`
- historical diff CSVs
- `outputs/manifest.json`

Additional migration task A:
- Add `outputs/threshold_characterization.csv`.
- Add claim-map aggregate rows for rigorous threshold characterization.

Data reconstruction:
- No historical rerun needed. `outputs/table1_step2_beam_top10.csv` already contains both IS and exp Bregman expressions for all 25 corpora.

Planned new `threshold_characterization.csv` columns:
- `corpus`
- `fitted_c`
- `is_rmse`
- `exp_rmse`
- `signed_gap_is_minus_exp`
- `winner`
- `is_rank`
- `exp_rank`

Planned aggregate rows:
- `is_exp_gap_crossing_c_linear_interp`
- `is_exp_gap_abs_lt_0p01_c_min`
- `is_exp_gap_abs_lt_0p01_c_max`
- `is_exp_gap_abs_lt_0p01_count`
- `is_exp_gap_negative_count`
- `is_exp_gap_positive_count`

Risk:
- This extends the v4 map; the new rows are not currently in v4 but are needed to make Section 3.1 threshold language rigorous after T3v2.

### 3a. Smooth two-regime fits

Current status: complete by v3/v4.

No migration planned unless the user chooses to revise Phase 2 around the decoupled-erf model.

### 3c. Statistical k scaling

Current status: complete by v3/v4.

No migration planned unless the user chooses to revise Phase 2 around the decoupled-erf model.

## Foundation Migrations

### 3d. POS crossover scaling

Priority: foundation; blocks `J2`.

Source bundles:
- `results/zipf_pos_all_corpora/`
- `results/zipf_pos_manual_v2/`

Historical files present:
- `results/zipf_pos_all_corpora/pos_all_corpora_points.csv`
- `results/zipf_pos_all_corpora/summary.json`
- `results/zipf_pos_all_corpora/report.md`
- 25 per-corpus `*_top500_pos.csv` files
- `results/zipf_pos_manual_v2/manual_alpha_points.csv`
- `results/zipf_pos_manual_v2/summary.json`
- `results/zipf_pos_manual_v2/report.md`
- 3 manual `*_top300_manual_pos.csv` files

Expected canonical outputs:
- `outputs/pos_scaling_points.csv`
- `outputs/manual_validation.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `forced_alpha`
- `forced_alpha_ci_low`
- `forced_alpha_ci_high`
- `forced_alpha_vs_half_pvalue`
- `manual_validation_corpus_count`

Data reconstruction plan:
- Copy/normalize `pos_all_corpora_points.csv` into `pos_scaling_points.csv`.
- Copy/normalize manual rows from `manual_alpha_points.csv` and manual top-300 files into `manual_validation.csv`.
- Read inferential values from `summary.json`; if summary lacks CI/p-value components, compute from saved point rows using the same method described in the historical report and document the method in README.

Potential blockers:
- Manual validation appears to have 3 saved manual top-300 files, while the map references a four-corpus manual cross-check. Need inspect `summary.json` to confirm whether the fourth corpus is represented only in aggregate form. If absent, emit `BLOCKED.md` for `manual_validation_corpus_count`.

### 6. Multilingual extension

Priority: foundation; blocks `J1`.

Source bundles:
- `results/zipf_multilang_romance/`
- `results/zipf_multilang_verify/`
- Legacy `results/zipf_multilang/` may be used only as superseded provenance if a row is missing from the canonical two bundles.

Historical files present:
- `results/zipf_multilang_romance/multilang_table.csv`
- `results/zipf_multilang_romance/summary.json`
- `results/zipf_multilang_romance/report.md`
- `results/zipf_multilang_verify/summary.json`
- `results/zipf_multilang_verify/report.md`

Expected canonical outputs:
- `outputs/table3_multilang.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `non_english_corpus_count`
- `smooth_beats_single_zm_count`
- `multilang_bregman_winner_count`
- `transition_fraction_min`
- `transition_fraction_max`
- `transition_fraction_near_half_count`
- `latin_french_spanish_dutch_c_lt_5_count`
- `very_low_c_multilang_count`

Data reconstruction plan:
- Join the 7-row romance summary/table with `zipf_multilang_verify` c/winner verification fields.
- Emit one row per non-English corpus with Table 3 columns: `language`, `vocabulary_size`, `zm_c`, `single_zm_rmse`, `smooth_rmse`, `step2_winner_display`, `k_over_sqrt_v`.
- Compute aggregate counts directly from the 7 canonical rows.

Potential blockers:
- Need verify that `zm_c` exists in `zipf_multilang_verify/summary.json` for all 7 rows. If not, recomputing c from cached texts is possible but would be new computation, not pure migration; decide whether allowed.

## Core Mechanism Migrations

### 2a. Function-word ablation

Priority: core mechanism.

Source bundle:
- `results/zipf_function_word_test/`

Historical files present:
- top-level `summary.json` and `report.md`
- subruns `remove_top50/`, `remove_top100/`, `only_top100/` each with `summary.json` and `report.md`

Expected canonical outputs:
- `outputs/function_word_ablation.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- No `metric_name` rows parsed from v4 for this experiment; values are table/row-cell claims from `function_word_ablation.csv`.

Data reconstruction plan:
- One row per ablation condition, carrying `condition`, `token_count`, `vocabulary_size`, `zm_a`, `zm_b`, `zm_c`, `zm_rmse`, `step2_winner`, and any top-5 beam fields saved.
- Aggregate file may contain no map-cited metric rows unless v4 inspection finds direct `metric_name` row keys during implementation.

Potential blockers:
- Low risk; needed c/RMSE values appear in saved summaries.

### 2b. Synthetic two-regime mixtures

Priority: core mechanism.

Source bundles:
- `results/zipf_synthetic_mixture/`
- `results/zipf_breakthrough_probe/` for planted-mixture follow-up

Historical files present:
- `zipf_synthetic_mixture/summary.json`
- `zipf_synthetic_mixture/report.md`
- `zipf_breakthrough_probe/summary.json`
- `zipf_breakthrough_probe/report.md`

Expected canonical outputs:
- `outputs/synthetic_mixture_runs.csv`
- `outputs/planted_mixture_runs.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `synthetic_mixture_configuration_count`

Direct row claims to satisfy:
- Four synthetic configurations from line `115`.
- Same-exponent control from line `222`.
- Small-gap mixture from line `224`.
- Larger exponent-gap rows from line `226`.
- Planted-mixture row from line `277`: single-ZM RMSE `0.181`, step-2 residual RMSE `0.174`, MOEZipf RMSE `0.326`, MOE + step-2 RMSE `0.263`, winner `log(x)^2`.

Data reconstruction plan:
- Extract synthetic-mixture rows from `zipf_synthetic_mixture/summary.json`.
- Extract planted-mixture row from `zipf_breakthrough_probe/summary.json`.
- If planted fields are nested differently, search exact values in the breakthrough summary/report and normalize into one `mixture_id=planted_small_gap` row.

Potential blockers:
- Stop if any of the exact line-277 numbers cannot be located in saved outputs. This was a known historical provenance mismatch risk.

### 2c. Mechanism absorption in residual space

Priority: core mechanism.

Source bundle:
- `results/zipf_breakthrough_probe/`

Historical files present:
- `summary.json`
- `report.md`

Expected canonical outputs:
- `outputs/residual_absorption_by_corpus.csv`
- `outputs/head_basis_by_corpus.csv`
- `outputs/generator_rescore_by_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `smooth_residual_step2_gain_median`
- `smooth_residual_step2_help_count`
- `smooth_residual_step2_hurt_count`
- `smooth_residual_step2_nonhelp_count`
- `moe_residual_step2_help_count`
- `moe_residual_step2_gain_median`
- `head_basis_r2_median_top200`
- `slant_vs_c_correlation`
- `median_winner_minus_euclidean_gap_full_rmse`
- `median_winner_minus_euclidean_gap_top100_rmse`
- `empirical_euclidean_gap_full_rmse_median`
- `empirical_euclidean_gap_top100_rmse_median`
- `smooth_synthetic_euclidean_gap_full_rmse_median`
- `smooth_synthetic_euclidean_gap_top100_rmse_median`
- `zm_control_euclidean_gap_top100_rmse_median`

Data reconstruction plan:
- Split `breakthrough_probe` summary into three canonical tables: residual absorption, head basis, generator rescore.
- Compute cited aggregate rows from the saved per-corpus values if present; otherwise use saved aggregate fields.

Potential blockers:
- This is a dense, multi-purpose bundle. If per-corpus rows are not saved for all 25, output aggregate rows can still be migrated from saved aggregates, but per-corpus CSVs may need `BLOCKED.md` detail.

### 3b. Model-family BIC comparison

Priority: core mechanism.

Source bundles:
- `results/zipf_continuous_piecewise/`
- `results/zipf_moezipf_comparison/`
- `results/zipf_bic_comparison/`
- `results/zipf_bic_landscape/`
- `results/zipf_sqrt_v_all_corpora/`

Historical files present:
- `zipf_continuous_piecewise/continuous_piecewise_table.csv`, `summary.json`
- `zipf_moezipf_comparison/moezipf_table.csv`, `summary.json`
- `zipf_bic_comparison/bic_table.csv`, `summary.json`
- `zipf_sqrt_v_all_corpora/summary.json`
- `zipf_bic_landscape/` currently has no saved files

Expected canonical outputs:
- `outputs/table2_model_family.csv`
- `outputs/model_family_rmse_per_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `piecewise_loses_to_smooth_on_all_25`
- `winner_count_smooth_family_total`
- `moezipf_bic_winner_count`
- `winner_count_moezipf`
- `median_rmse_cost_sqrtv_vs_free_k`

Data reconstruction plan:
- Use `zipf_bic_comparison/bic_table.csv` as the primary Table 2 source.
- Use `moezipf_table.csv`, `continuous_piecewise_table.csv`, and `zipf_sqrt_v_all_corpora/summary.json` to fill RMSE columns or verify table values.
- Emit Table 2 exactly as v4 map requires: `single_zm_bic`, `moe_bic`, `hard_piecewise_bic`, `continuous_piecewise_bic`, `smooth_ksqrtv_bic`, `smooth_freek_bic`, `winner_family`.

Potential blockers:
- `results/zipf_bic_landscape/` is empty. If v4 contains a landscape-specific claim not already represented in `bic_table.csv`, this experiment must emit `BLOCKED.md` for that claim. Current parsed v4 metrics do not show a landscape-specific metric.

### 1b. Non-linguistic Zipfian control

Priority: controls/robustness.

Source bundles:
- `results/zipf_enriched_city_populations/`
- `results/zipf_enriched_city_populations_step2_only/` as verification/source fallback

Historical files present:
- each has `summary.json`, `report.md`

Expected canonical outputs:
- `outputs/city_population_control.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- No parsed `metric_name` rows; line `210` is satisfied by row cells in `city_population_control.csv`.

Data reconstruction plan:
- Emit one row for `world_city_populations`, with city count, V, ZM c, ZM RMSE, step-2/composite RMSE, winner expression, and whether the correction helps.

Potential blockers:
- Low risk; exact values should be in `summary.json`.

### 1c. Search-robustness ablations

Priority: controls/robustness.

Source bundles:
- `results/zipf_boundary_ablation/`
- `results/zipf_guard_ablation/`
- `results/zipf_wls_test/`

Historical files present:
- boundary per-condition summaries plus top-level `summary.json`
- guard per-condition summaries plus top-level `summary.json`
- WLS `summary.csv`, `summary.json`, `report.md`

Expected canonical outputs:
- `outputs/robustness_ablation_by_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- No parsed `metric_name` rows; claims are prose robustness statements.

Data reconstruction plan:
- Normalize all ablation rows into one schema: `ablation_family`, `corpus`, `setting`, `x_low`/`guard_mode`/`weighting`, `winner_expression`, `winner_family`, `rmse`, `matches_canonical`.
- Emit aggregate rows only if v4 implementation inspection finds explicit robustness metrics.

Potential blockers:
- Historical ablations cover selected diagnostic corpora, not all 25. README must state scope exactly.

### 1d. Winner-vs-Euclidean gap analysis

Priority: controls/robustness.

Source bundle:
- `results/zipf_english_gap_verify/`

Historical files present:
- `english_gap_table.csv`
- `summary.json`
- `report.md`

Expected canonical outputs:
- `outputs/gap_analysis_per_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `euclidean_step2_winner_count`
- `median_winner_minus_euclidean_gap_full_rmse`
- `max_core_generator_gap_full_rmse_typical_band`
- `winner_family_gap_below_0p01_count`

Data reconstruction plan:
- Copy/normalize `english_gap_table.csv`.
- Compute aggregate rows from per-corpus columns.
- Coordinate with `1a` threshold-characterization so IS-vs-exp and Euclidean-gap claims are not duplicated inconsistently.

Potential blockers:
- Low risk if `english_gap_table.csv` has all required columns.

## Structural and Scope Migrations

### 4. Analytic seam theory

Source bundles:
- `results/zipf_seam_sign_theory/`
- `results/zipf_seam_projection_theory/`
- `results/zipf_seam_second_order_theory/`

Expected canonical outputs:
- `outputs/seam_sign_checks.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `raw_sign_prediction_match_count_internal`
- `internal_predicted_sign_match_count`
- `raw_seam_full_sign_match_head10`
- `raw_seam_full_sign_match_head50`
- `raw_seam_term_match_head10_linear`
- `raw_seam_term_match_head10_quadratic`
- `raw_seam_term_match_head10_cubic`
- `raw_seam_term_match_head50_linear`
- `raw_seam_term_match_head50_quadratic`
- `raw_seam_term_match_head50_cubic`
- `projected_vs_empirical_full_match_count`
- `projected_vs_empirical_full_match_head50`
- `projected_vs_empirical_match_head50_linear`
- `projected_vs_empirical_match_head50_quadratic`
- `projected_vs_empirical_match_head50_cubic`
- `projected_vs_numerical_refit_full_match_count`
- `projected_vs_numerical_full_match_head50`
- `projected_vs_numerical_match_head50_linear`
- `projected_vs_numerical_match_head50_quadratic`
- `projected_vs_numerical_match_head50_cubic`
- `numerical_refit_full_sign_match_head50`
- `numerical_refit_match_head50_linear`
- `numerical_refit_match_head50_quadratic`
- `numerical_refit_match_head50_cubic`
- `second_order_full_sign_match_head50`

Data reconstruction plan:
- Extract aggregate counts from the three summaries.
- If per-corpus sign-check rows are saved, emit them in `seam_sign_checks.csv`; otherwise emit the aggregate rows and mark per-corpus detail as unavailable.

Potential blockers:
- If v5 decoupled-erf model remains canonical, old coupled-logistic seam theory is scientifically superseded. Under v4 migration, values are reconstructible.

### 5a. Simulation recovery from fitted smooth models

Source bundle:
- `results/zipf_simulation_recovery/`

Expected canonical outputs:
- `outputs/simulation_recovery_per_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `smooth_exact_winner_match_rate`
- `single_zm_control_exact_winner_match_rate`
- `smooth_step2_help_rate`
- `single_zm_control_step2_help_rate`
- `single_zm_control_step2_help_count`
- `high_c_block_size`
- `high_c_total_count`
- `high_c_exact_match_count`
- `high_c_exact_match_rate_smooth`
- `high_c_help_rate_smooth`
- `high_c_help_rate_zm_control`
- `basis_corr_linear_smooth`
- `basis_corr_quadratic_smooth`
- `basis_corr_cubic_smooth`
- `basis_corr_linear_zm`
- `basis_corr_quadratic_zm`
- `basis_corr_cubic_zm`
- `smooth_basis_corr_linear`
- `smooth_basis_corr_quadratic`
- `smooth_basis_corr_cubic`
- `zm_control_basis_corr_linear`
- `zm_control_basis_corr_quadratic`
- `zm_control_basis_corr_cubic`

Data reconstruction plan:
- Copy/normalize `simulation_recovery_table.csv`.
- Extract aggregate rates/correlations from `summary.json`.

Potential blockers:
- Low risk; source has both table and summary.

### 5b. Low-c manifold structure

Source bundles:
- `results/zipf_lowc_manifold_analysis/`
- `results/zipf_phase_coordinate/`

Expected canonical outputs:
- `outputs/lowc_manifold_per_corpus.csv`
- `outputs/phase_coordinate_per_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `lowc_median_cosine_exp_vs_xpow_top200`
- `lowc_full_rmse_exp_winner_count`
- `lowc_top100_xpow_winner_count`
- `lowc_top200_xpow_winner_count`
- `lowc_median_span_r2_exp_xpow_top200`
- `lowc_median_span_r2_exp_is_top200`
- `lowc_smooth_modal_vs_empirical_top100_match_rate`
- `lowc_zm_modal_vs_empirical_top100_match_rate`
- `lowc_smooth_modal_vs_empirical_full_match_rate`
- `lowc_median_delta_xpow_minus_exp_full_rmse`
- `lowc_median_delta_xpow_minus_exp_top100_rmse`

Data reconstruction plan:
- Extract manifold rows from `zipf_lowc_manifold_analysis/summary.json`.
- Extract phase-coordinate rows from `zipf_phase_coordinate/summary.json`.
- Compute aggregate rows from saved per-corpus data or migrate saved aggregate values.

Potential blockers:
- If per-corpus phase rows are not saved, top-100/top-200 flip counts may still be in aggregate; per-corpus CSV would need a `BLOCKED.md` note.

### 5c. Smooth-parameter control of low-c manifold

Source bundle:
- `results/zipf_smooth_parameter_sweep/`

Expected canonical outputs:
- `outputs/parameter_sweep_rows.csv`
- `outputs/parameter_sweep_correlations.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `synthetic_full_winner_count_exp`
- `synthetic_full_winner_count_is`
- `synthetic_full_winner_count_xpow`
- `synthetic_top100_winner_count_exp`
- `synthetic_top100_winner_count_is`
- `synthetic_top100_winner_count_xpow`

Parameter correlation file rows:
- `parameter=transition_fraction`, fields for `r = -0.540`
- `parameter=sigmoid_width`, fields for `r = +0.554`

Data reconstruction plan:
- Extract sweep rows and correlation table from summary/report.
- If no explicit CSV is saved, reconstruct rows from JSON arrays in `summary.json`.

Potential blockers:
- Current v5 drift tracker states these are coupled-family diagnostics, not decoupled-erf canonical claims.

## PMF and Discrete-Likelihood Migrations

### 7a. Canonical PMF family

Source bundles:
- `results/zipf_seam_mandelbrot_pmf/`
- `results/zipf_seam_mandelbrot_regularized/`
- `results/zipf_seam_mandelbrot_softk_splitfit/`
- `results/zipf_seam_mandelbrot_softkw/`
- Verification-only comparison: `results/zipf_v4_verification/table_a_fourway_pmf.csv`

Expected canonical outputs:
- `outputs/splitfit/table4_fourway.csv`
- `outputs/splitfit/pmf_variant_per_corpus.csv`
- `outputs/splitfit/aggregate_statistics.csv`
- `outputs/fullrefit/fourway_per_corpus.csv`
- `outputs/fullrefit/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `softk_beats_moe_count`
- `softk_beats_zm_count`
- `softk_beats_zipf_count`
- `softk_beats_free_count`
- `softk_beats_fixedk_count`
- `winner_count_zm`
- `winner_count_softk`
- `winner_count_moe`
- `winner_count_zipf`
- `median_softk_minus_moe`
- `median_softk_minus_zm`
- `median_softk_minus_zipf`
- `median_softk_minus_free`
- `median_softk_minus_fixedk`
- `softk_step2_help_count`
- `softk_step2_nonhelp_count`
- `softk_step2_gain_median`
- `free_pmf_step2_help_count`
- `free_pmf_step2_gain_median`
- `free_pmf_beats_moe_count`
- `free_pmf_beats_zm_count`
- `free_pmf_train_bic_win_count`
- `free_pmf_heldout_winner_count_zm`
- `free_pmf_heldout_winner_count_moe`
- `free_pmf_heldout_winner_count_seam`
- `free_pmf_heldout_winner_count_zipf`
- `fixedk_beats_free_count`
- `fixedkw_beats_free_count`
- `median_fixedk_minus_free`
- `median_fixedkw_minus_free`
- `fixedk_step2_help_count`
- `fixedkw_step2_help_count`
- `fixedk_step2_gain_median`
- `fixedkw_step2_gain_median`
- `lambda_k_min`
- `lambda_k_max`
- `lambda_count_1e`
- `lambda_count_3e`
- `softk_parameter_count`

Data reconstruction plan:
- Use splitfit soft-k bundle as canonical held-out source.
- Use PMF/free and regularized bundle CSVs to populate variant comparison table.
- Use `zipf_v4_verification/table_a_fourway_pmf.csv` only as a historical verification target for Table 4, not as primary data if upstream bundle values are available.
- README must state splitfit metrics are canonical and fullrefit diagnostics are non-canonical.

Potential blockers:
- If table A values are hand-assembled and not reproducible from upstream CSVs, stop with drift/gap report rather than copying hand assembly silently.

### 7b. PMF head-window evaluation

Source bundle:
- `results/zipf_angle1_head_windows_splitfit/`

Expected canonical outputs:
- `outputs/head_window_per_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `top50_modal_winner`
- `top100_modal_winner`
- `top200_modal_winner`
- `top500_modal_winner`
- `top1000_modal_winner`
- `top100_winner_count_zm`
- `top100_softk_beats_moe_count`
- `top100_softk_step2_help_count`
- `top200_winner_count_zipf`
- `top200_winner_count_zm`
- `top200_winner_count_moe`
- `top200_winner_count_softk`
- `top200_winner_count_summary`
- `top200_softk_step2_help_count`
- `top1000_winner_count_zipf`
- `top1000_winner_count_zm`
- `top1000_winner_count_moe`
- `top1000_winner_count_softk`
- `top1000_winner_count_summary`
- `top1000_softk_step2_help_count`
- `top50_softk_step2_help_count`
- `top500_softk_step2_help_count`
- `full_softk_step2_help_count`
- `max_head_window_rank`

Data reconstruction plan:
- Normalize `head_window_table.csv`.
- Compute modal winners and counts by `K`.

Potential blockers:
- Low risk if the splitfit head-window table contains all K levels.

### 7c. PMF lambda metadata

Source bundle:
- `results/zipf_angle2_lambda_metadata/`

Expected canonical outputs:
- `outputs/lambda_metadata_summary.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `max_abs_metadata_correlation`

Data reconstruction plan:
- Copy/normalize `lambda_metadata_table.csv`.
- Compute max absolute correlation from table rows.

Potential blockers:
- Low risk.

### 7d. PMF asymmetric gate

Source bundle:
- `results/zipf_angle3_asymmetric_gate_splitfit/`

Expected canonical outputs:
- `outputs/asymmetric_gate_per_corpus.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `median_width_ratio`
- `material_asymmetry_count`
- `asymmetric_beats_softk_count`
- `median_asymmetric_minus_softk`
- `asymmetric_step2_help_count`

Data reconstruction plan:
- Copy/normalize `asymmetric_gate_table.csv`.
- Compute aggregates from per-corpus rows.

Potential blockers:
- Low risk.

### 7e. PMF hierarchical pooling on k

Source bundle:
- `results/zipf_angle4_hierk/`

Expected canonical outputs:
- `outputs/hierk_summary.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `hierk_vs_softk_summary`

Direct row claims:
- line `411`: `alpha = 0.505`, `sigma = 0.05` lower-bound behavior.

Data reconstruction plan:
- Normalize `hierk_head_window_table.csv` and `summary.json` into `hierk_summary.csv`.
- Emit textual summary row in `aggregate_statistics.csv` if the map expects a compact statement.

Potential blockers:
- Low risk.

### 8. Bible per-book decomposition

Source bundle:
- `results/zipf_angle6_bible_books/`

Expected canonical outputs:
- `outputs/bible_per_book.csv`
- `outputs/table5_bible_summary.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `books_analyzed`
- `per_book_step2_help_count`
- `per_book_step2_nonhelp_count`
- `per_book_softk_beats_moe_count`
- `per_book_softk_beats_zm_count`
- `whole_bible_single_fit_softk_heldout_nll`
- `aggregate_per_book_softk_heldout_nll`
- `whole_bible_singlefit_softk_nll`
- `aggregate_per_book_softk_nll`
- `improvement_from_decomposition`
- `median_per_book_step2_gain`
- `median_per_book_softk_minus_moe`
- `median_per_book_softk_minus_zm`

Data reconstruction plan:
- Normalize `bible_books_table.csv` into per-book rows.
- Build `table5_bible_summary.csv` directly from saved summary fields.
- Ensure Table 5 summary rows are not builder-hardcoded.

Potential blockers:
- Low risk if summary contains whole-Bible baseline; otherwise need retrieve from table/report and flag if absent.

### 9. Likelihood frontier beyond soft-k

Source bundles:
- `results/zipf_hybrid_headtail_splitfit/`
- `results/zipf_hybrid_mechpenalty/`
- `results/zipf_hybrid_vs_softk_analysis/`
- `results/zipf_fully_ours_variants/`

Expected canonical outputs:
- `outputs/likelihood_frontier_per_corpus.csv`
- `outputs/hybrid_structure_summary.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `hybrid_beats_moe_count`
- `hybrid_beats_softk_count`
- `hybrid_step2_help_count`
- `median_hybrid_minus_softk`
- `nested_beats_moe_count`
- `nested_beats_softk_count`
- `nested_beats_zm_count`
- `nested_step2_help_count`
- `nested_fourway_winner_count`
- `three_regime_beats_softk_count`
- `median_three_regime_minus_softk`
- `mechanism_penalty_grid`
- `mechanism_penalty_frontier_result`
- `mechanism_penalty_step2_help_count_all_lambdas`

Data reconstruction plan:
- Join the four historical CSVs on corpus where needed.
- Use `hybrid_vs_softk_table.csv` for pairwise hybrid-vs-softk claims.
- Use `fully_ours_variants_table.csv` for nested/three-regime claims.
- Use `hybrid_mechpenalty_table.csv` for penalty-grid/frontier claims.

Potential blockers:
- If variant names or row keys differ between historical tables, create explicit mapping in README and manifest.

## Robustness Migrations

### 10a. Search-depth robustness

Source bundles:
- `results/zipf_step10_ablation/`
- `results/zipf_head_poly_decomposition/`
- `results/zipf_head_poly_transfer/`

Expected canonical outputs:
- `outputs/step10_poly_per_corpus.csv`
- `outputs/poly_transfer.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- Parser did not extract named `metric_name` rows for this experiment; claims are row/value claims around degree-5 R² and transfer RMSE values.

Direct row claims:
- line `448`: degree-5 R² `0.977`, `0.993`; RMSE comparisons `0.1506 vs 0.1599`, `0.1458 vs 0.1453`.
- lines `506 / 538`: Shakespeare→War and Peace and War and Peace→Bible transfer restatements.

Data reconstruction plan:
- Extract per-corpus step-10 outcomes from `zipf_step10_ablation/summary.json`.
- Extract polynomial decompositions from `zipf_head_poly_decomposition/summary.json`.
- Extract transfer rows from `zipf_head_poly_transfer/summary.json`.

Potential blockers:
- If historical summaries only save aggregate values and not per-corpus rows, `step10_poly_per_corpus.csv` may be partially blocked.

### 10b. Grammar-width robustness

Source bundles:
- `results/zipf_widened_grammar_diagnostic/`
- `results/zipf_widened_grammar_extended/`

Historical files present:
- `zipf_widened_grammar_extended/widened_lowc_summary.csv`
- both summaries and reports

Expected canonical outputs:
- `outputs/widened_beam_top10.csv`
- `outputs/table6_widened_manifold.csv`
- `outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `english_lowc_flip_count`
- `english_lowc_total_count`
- `lowc_corpus_count_extended`
- `lowc_widened_english_flip_count`
- `lowc_widened_nonenglish_xpow_count`
- `lowc_widened_retains_exp_count`
- `lowc_widened_manifold_yes_count`
- `lowc_widened_total_count`
- `lowc_widened_span_r2_threshold`
- `lowc_manifold_yes_count`
- `lowc_manifold_total_count`
- `lowc_manifold_min_span_r2`
- `flip_case_rmse_gap_min`
- `flip_case_rmse_gap_max`
- `flip_case_rmse_gap_median`
- `flipped_exp_gap_min`
- `flipped_exp_gap_max`
- `flipped_exp_gap_median`
- `flipped_exp_rank2_count`
- `flipped_exp_rank3_count`
- `flip_cases_old_exp_rank2_count`
- `flip_cases_old_exp_rank3_count`

Data reconstruction plan:
- Use extended `widened_lowc_summary.csv` as primary Table 6 source.
- Use diagnostic summary for Shakespeare, War and Peace, Pride and Prejudice high-c/diagnostic rows.
- Build Table 6 body from extended rows.

Potential blockers:
- Need verify Table 6 expected column names (`cos_vs_exp`, `cos_vs_xpow`, `span_r2`) are present or reconstructible from widened summary.

## New Experiment 11

### 11. Aligned-coordinate robustness

Status: new experiment; not historical migration.

Source data:
- Existing 25 English corpus registry and raw corpora under `data/zipf/`
- Existing single-ZM fit code
- Existing step-2 SR machinery

Expected canonical outputs:
- `experiments/11_aligned_coordinate_robustness/outputs/aligned_coordinate_comparison.csv`
- `experiments/11_aligned_coordinate_robustness/outputs/aggregate_statistics.csv`
- `outputs/manifest.json`

Planned per-corpus columns:
- `corpus`
- `vocabulary_size`
- `fitted_c`
- `logrank_winner`
- `aligned_winner`
- `logrank_top5_json`
- `aligned_top5_json`
- `logrank_cos_is`
- `logrank_cos_exp`
- `logrank_cos_euclidean`
- `aligned_cos_is`
- `aligned_cos_exp`
- `aligned_cos_euclidean`
- `winner_same`
- `new_attractor_flag`

Planned aggregate rows:
- `aligned_coordinate_same_winner_count`
- `aligned_coordinate_same_winner_fraction`
- `logrank_is_winner_count`
- `logrank_exp_winner_count`
- `aligned_is_winner_count`
- `aligned_exp_winner_count`
- `aligned_other_winner_count`
- `aligned_coordinate_preserves_is_exp_split`
- `aligned_new_attractor_count`

Data reconstruction need:
- Cannot be reconstructed from historical bundles because aligned coordinate `x' = 0.05 + 0.95 log(r+c_fitted)/log(V+c_fitted)` was not saved. Requires a new light step-2 run after approval.

Stop risk:
- If aligned coordinate produces materially different winners, this may require manuscript interpretation rather than simple robustness support.

## Joined Producers

### J1. Paper-scope rollup

Priority: after `1a`, `3a`, and `6`.

Dependencies:
- `1a` complete
- `3a` complete by v3
- `6` pending

Expected canonical outputs:
- `joins/J1_paper_scope_rollup/outputs/paper_scope_summary.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `english_corpus_count`
- `non_english_corpus_count`
- `total_corpus_count`
- `language_family_count`
- `smooth_beats_single_zm_total_count`
- `bregman_indoeuropean_alphabetic_count`
- `indoeuropean_alphabetic_total_count`

Data reconstruction plan:
- Read only canonical outputs from `1a`, `3a`, and `6`.
- Emit one row per metric.

Potential blockers:
- Blocked until `6` completes.
- If using v5 decoupled-erf manuscript, `smooth_beats_single_zm_total_count` should specify which smooth model family is being counted.

### J2. Scaling comparison

Priority: after `3c` and `3d`.

Dependencies:
- `3c` complete
- `3d` pending

Expected canonical outputs:
- `joins/J2_scaling_comparison/outputs/scaling_comparison_summary.csv`
- `outputs/manifest.json`

v4 aggregate metrics to emit:
- `alpha_kstat`
- `alpha_pos`
- `alpha_difference_pos_minus_kstat`
- `alpha_difference_ci_low`
- `alpha_difference_ci_high`
- `alpha_difference_pvalue`

Data reconstruction plan:
- Read canonical aggregate rows from `3c` and `3d`.
- Compute explicit paired comparison only in producer script; LaTeX builder computes nothing.

Potential blockers:
- Blocked until `3d` completes.
- If v5 decoupled-erf k-scaling is canonical, J2 must compare POS alpha to the decoupled-erf scaling output instead of v3 `3c`.

## Priority Execution Order After Approval

1. Extend `1a` threshold-characterization output.
2. Migrate `3d`.
3. Migrate `6`.
4. Migrate `2a`, `2b`, `2c`.
5. Migrate `3b`.
6. Migrate `1b`, `1c`, `1d`.
7. Migrate `4`, `5a`, `5b`, `5c`.
8. Migrate `7a`, `7b`, `7c`, `7d`, `7e`.
9. Migrate `8`.
10. Migrate `9`.
11. Migrate `10a`, `10b`.
12. Run new `11_aligned_coordinate_robustness`.
13. Build `J1`.
14. Build `J2`.
15. Rerun Phase 2 audit.
16. Rerun audit-package dry-run manifest.

## Progress Reporting Plan

After approval and during execution:
- Report after each of the first two foundation migrations (`3d`, `6`) because they unblock joins.
- Then report after every 3 experiments migrated.
- For each report include: files written, row counts, metrics emitted, drift checks, and any `BLOCKED.md` files.

## Known Planning Risks

1. `3b`: `results/zipf_bic_landscape/` is empty. Likely not fatal because v4 appears to cite table/BIC counts rather than a landscape-only metric, but this must be verified during migration.
2. `3d`: v4 claims manual cross-check on four corpora, while `zipf_pos_manual_v2` visibly contains three manual top-300 CSVs. Summary may contain the fourth; otherwise this is a canonical gap.
3. `2b`: line-277 planted-mixture numbers have known provenance sensitivity. Stop if exact values cannot be located.
4. `7a`: Table 4 historically involved a hand-assembled verification CSV. Use it as a verification target, not as a silent source, unless upstream rows cannot reconstruct it.
5. `4`, `5c`, `3a`, `3c`: if manuscript v5.1 remains canonical, these older coupled-logistic outputs may be scientifically superseded by decoupled-erf work outside the v3/v4 Phase 2 structure.
6. Experiment `11` is a new run, not migration. It should begin only after explicit approval.

