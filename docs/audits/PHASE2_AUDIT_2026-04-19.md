# Phase 2 Consolidation Audit — 2026-04-19

## Executive Summary

Audit source documents:
- `results/CONSOLIDATION_PROPOSAL_v3.md`, modified `2026-04-17 11:07:02`.
- `results/MANUSCRIPT_CLAIM_TO_CSV_MAP_v3.md`, modified `2026-04-17 11:20:39`.

Important versioning note: the repo also contains `results/MANUSCRIPT_CLAIM_TO_CSV_MAP_v4.md`, modified `2026-04-17 11:31:59`, which declares itself authoritative. This audit follows the user's requested v3 map, reading Part II as superseding Part I where they overlap. Before Phase 2 resumes, document authority should be re-confirmed because v3 is not a clean single-pass map, while v4 claims to be one.

Current Phase 2 status by the finalized v3 structure:
- 25 consolidated experiments total.
- Complete by v3 canonical-output criteria: `3` of `25` (`1a`, `3a`, `3c`).
- Partial by v3 canonical-output criteria: `0` of `25`.
- Pending by v3 canonical-output criteria: `22` of `25`.
- Blocked by missing source data: `0` of `25` identified from this audit, although several pending experiments have incomplete upstream dependencies.
- Joined-output producers complete: `0` of `2`.
- Joined-output producers blocked/pending: `2` of `2` (`J1`, `J2`).

Canonical experiment directories currently present:
- `experiments/1a_per_corpus_enriched_search/`
- `experiments/3a_smooth_two_regime_fits/`
- `experiments/3c_kstat_scaling/`

Canonical `joins/` directory status:
- `joins/` does not exist.

Non-canonical add-on outputs present:
- S2 v3 decoupled five-gate outputs are present under `results/s2_v3_windows_full_outputs_2026-04-18/` and `results/windows_s2_decoupled_v3_results_2026-04-18/`.
- S2 synthetic recovery outputs are present under `results/s2_synthetic_recovery_outputs_workers32_2026-04-18/`.
- Decoupled-erf k-scaling recheck is present at `phase2_addon/k_scaling_decoupled_erf_recheck.csv`.
- T1/T2/T3/T3v2 theorem-track outputs are present under `phase2_addon/`.

These add-on outputs are not yet integrated into the 25-experiment/J1/J2 Phase 2 contract. In particular, `MANUSCRIPT_DRAFT_v5.md` migrated the smooth model toward the decoupled 9-parameter erf-gate result, but `CONSOLIDATION_PROPOSAL_v3.md` and the canonical `experiments/3a_smooth_two_regime_fits/` output still describe the older bounded 8-parameter coupled logistic model.

## Section 3.1 Threshold Flag

Experiment `1a` currently emits enough raw canonical data to reconstruct an IS-vs-exp RMSE gap because `outputs/table1_step2_beam_top10.csv` contains both `sub[sub[x,1],log[x]]` and `eml[sub[x,1],eml[x,1]]` for all 25 corpora. However, `1a` does **not** emit a dedicated per-corpus `IS_rmse`, `exp_rmse`, `is_minus_exp_rmse`, or threshold/gap-crossing CSV.

Given T3v2's finding that symbolic winner transitions arise from smooth RMSE-gap crossings, the Section 3.1 threshold language (`c` approximately `79` and `66`) should not close on winner identity alone. Before Phase 2 is considered complete for Section 3.1, add or extend a canonical output, preferably in `1a` or a narrowly scoped `1a` threshold-diagnostic sub-output, to emit:
- one row per English corpus with `c`, `is_rmse`, `exp_rmse`, `is_minus_exp_rmse`, winner, and rank of each generator in the step-2 beam;
- aggregate threshold rows for the zero-crossing/gap band;
- notes connecting the empirical threshold to T3v2's smooth gap-crossing result.

Experiments needing T3v2 integration or threshold awareness:
- `1a`: direct owner of Section 3.1 c-band and IS/exp winner split.
- `1d`: related owner of winner-vs-alternate Bregman / Euclidean gap analysis.
- `5b` and `10b`: low-c manifold and grammar-width claims may cite the same c-band framing, but they are not direct owners of the Section 3.1 threshold.

## Per-Experiment Audit

### 1a. Per-corpus enriched search on the 25 English corpora

Status: **complete by v3**, with a Section 3.1 threshold-extension recommendation.

Canonical directory: `experiments/1a_per_corpus_enriched_search/`

Outputs present:
- `outputs/table1_per_corpus.csv` — 25 rows.
- `outputs/table1_step2_beam_top10.csv` — 250 rows.
- `outputs/aggregate_statistics.csv`.
- `outputs/historical_per_corpus_diff.csv`.
- `outputs/historical_top10_diff.csv`.
- `outputs/historical_diff_summary.csv`.
- `outputs/manifest.json`.

Outputs missing:
- None required by `CONSOLIDATION_PROPOSAL_v3.md` or the `1a` README.
- Missing for the new threshold concern: no dedicated IS-vs-exp RMSE-gap per-corpus/aggregate output.

Last modified:
- Most recent canonical output: `2026-04-17 14:17:37`.

README status:
- **exists / follows schema**.
- Includes research question, upstream dependencies, outputs and schemas, manuscript line ranges, rerun command, verification mapping, and no-inference methods note.

Manuscript claims currently satisfied:
- English corpus count: line `13` English side and line `29` (`25 English corpora`).
- Table 1 body: line `202 / Table 1 body`, from `outputs/table1_per_corpus.csv`.
- Step-2 winner family split: line `13`, line `199`, line `200`, and related Section 3.1 rows from `outputs/aggregate_statistics.csv`.
- Both Bregman forms in beam: line `202`, row `both_bregman_forms_present_in_step2_beam_count = 25`.
- Terminal enriched-search correction for line `197`: rows `terminal_all_improve_count = 25`, `terminal_improvement_pct_min = 3.6978883162152028`, `terminal_improvement_pct_max = 25.294536761006274`, `terminal_improvement_pct_mean = 11.023397106886518`, `terminal_improvement_pct_median = 9.798300633636057`.
- Step-2 diagnostic correction for line `202`: row `mean_step2_winner_runner_gap = 0.002009831063671038`, not the prose-only `0.003`.
- Protocol constants for lines `57` and `69` are represented in `source_config.py` / README.

Manuscript claims still unsatisfied:
- The paper-scope rollup claims using `1a + 6 + 3a` remain unsatisfied until `J1` is produced.
- Rigorous Section 3.1 c-threshold/gap-crossing characterization remains unsatisfied until an IS-vs-exp RMSE-gap output is emitted.

### 1b. Non-linguistic Zipfian control: city populations

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/city_population_control.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `210`: city-population control values (`33,535` cities, `V = 33,535`, `c = 100.4`, RMSE `0.1206`, composite RMSE `0.1225`).

### 1c. Search-robustness ablations

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/robustness_ablation_by_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Section 3.1 robustness statements about boundary, guard, and WLS stability remain unsatisfied as consolidated Phase 2 claims.

### 1d. Winner-vs-Euclidean gap analysis

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/gap_analysis_per_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `206`: Euclidean never wins step-2.
- Line `206`: median winner-vs-Euclidean full-RMSE gap.
- Line `206` and line `534`: full-RMSE tight-family / gap-below-`0.01` claims.
- Threshold-related integration note: this experiment should coordinate with the `1a` IS-vs-exp RMSE-gap extension so Section 3.1 has one coherent gap story.

### 2a. Function-word ablation

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/function_word_ablation.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `216`: function-word ablation c/RMSE claims (`c = 607.1`, `454.1`, `5.3`, RMSE `0.044`).

### 2b. Synthetic two-regime mixtures

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/synthetic_mixture_runs.csv`.
- `outputs/planted_mixture_runs.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `115`: four synthetic configurations.
- Lines `222`, `224`, `226`: same-exponent and exponent-gap mixture results.
- Line `277`: planted-mixture values (`0.181`, `0.174`, `0.326`, `0.263`, winner `log(x)^2`).
- Line `534`: synthetic mixture configuration count.

### 2c. Mechanism absorption in residual space

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/residual_absorption_by_corpus.csv`.
- `outputs/head_basis_by_corpus.csv`.
- `outputs/generator_rescore_by_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `13`: smooth residual no-help / MOE residual help claims.
- Line `31`: head-basis median top-200 R² and slant correlation.
- Lines `123` and `127`: helpful threshold and top-window protocol constants.
- Lines `273`, `279`, `339`, `492`, `494`, `544`: mechanism absorption, MOE residual, Euclidean-gap, head-basis, and generator-rescoring claims.

### 3a. Canonical smooth two-regime model fits

Status: **complete by v3**, with a model-canon integration warning.

Canonical directory: `experiments/3a_smooth_two_regime_fits/`

Outputs present:
- `outputs/smooth_fit_per_corpus.csv` — 25 rows.
- `outputs/bounds_robustness_per_corpus.csv` — 25 rows.
- `outputs/aggregate_statistics.csv`.
- `outputs/historical_smooth_diff.csv`.
- `outputs/historical_relaxed_diff.csv`.
- `outputs/historical_diff_summary.csv`.
- `outputs/manifest.json`.

Outputs missing:
- None required by `CONSOLIDATION_PROPOSAL_v3.md` or the `3a` README.
- Missing for current manuscript v5 model canon: no canonical `experiments/3a` output for the decoupled 9-parameter erf-gate model. S2 v3 outputs exist outside this experiment tree.

Last modified:
- Most recent canonical output: `2026-04-17 14:42:40`.

README status:
- **exists / follows schema** for the v3 8-parameter coupled-logistic model.
- No statistical inference method needed.

Manuscript claims currently satisfied by v3 map:
- Lines `95`, `97`, `99`, `103`: original smooth-model protocol constants.
- Lines `230–257`: smooth-fit side of Sections 3.5/3.6.
- Line `245`: mean/SD of `w`, `w < 0.2` list. Current canonical list is `Grimm's Fairy Tales; Aesop's Fables; Dubliners; Critique of Pure Reason`.
- Line `538`: relaxed bounds improve RMSE on `25/25`.
- English-side smooth coverage contribution to `J1`, pending the join.

Manuscript claims still unsatisfied:
- Any decoupled-erf canonical claims now in v5 are not represented in this v3 experiment directory.
- Any cross-experiment rollup using `3a` remains unsatisfied until `J1` exists.

S2 v3 add-on note:
- Decoupled five-gate real-data results are present at `results/s2_v3_windows_full_outputs_2026-04-18/`.
- Decoupled logistic-vs-coupled comparison and parameter-shift files are present at `results/windows_s2_decoupled_v3_results_2026-04-18/outputs/`.
- Synthetic logistic-gate recovery is present at `results/s2_synthetic_recovery_outputs_workers32_2026-04-18/outputs/`.
- These files are not in `experiments/3a_smooth_two_regime_fits/outputs/` and are not described by the `3a` README.

### 3b. Model-family BIC comparison

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/table2_model_family.csv`.
- `outputs/model_family_rmse_per_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `13`: smooth family dominates hard/continuous piecewise alternatives.
- Lines `85`, `101`: model-family protocol constants.
- Lines `232–243 / Table 2 body`.
- Line `243`: smooth family winner count.
- Line `247`: Moby Dick MOE/RMSE comparison.
- Line `257`: median RMSE cost of `k = sqrt(V)`.
- Lines `273`, `275`, `492`: MOE BIC winner counts and family comparison claims.

### 3c. Statistical scaling of transition centre `k_stat`

Status: **complete by v3**, with a decoupled-erf scaling warning.

Canonical directory: `experiments/3c_kstat_scaling/`

Outputs present:
- `outputs/kstat_scaling_points.csv` — 25 rows.
- `outputs/aggregate_statistics.csv`.
- `outputs/historical_point_diff.csv`.
- `outputs/historical_aggregate_diff.csv`.
- `outputs/historical_diff_summary.csv`.
- `outputs/manifest.json`.

Outputs missing:
- None required by `CONSOLIDATION_PROPOSAL_v3.md` or the `3c` README.
- Missing for current manuscript v5 model canon: no canonical `3c` rerun using decoupled erf `k`. The recheck exists only at `phase2_addon/k_scaling_decoupled_erf_recheck.csv`.

Last modified:
- Most recent canonical output: `2026-04-17 15:12:37`.

README status:
- **exists / follows schema**.
- Explicitly documents CI construction and p-value computation.

Manuscript claims currently satisfied by v3 map:
- Line `13`: `V^0.521` statistical scaling claim.
- Lines `255–257`: `k_stat` scaling alpha, CI, p-value.
- Line `516`: conclusion restatement of free-fit statistical exponent.
- Line `661`: appendix checklist alpha/CI/p-value.

Manuscript claims still unsatisfied:
- Line `265`: POS-vs-statistical paired comparison remains unsatisfied until `J2`.
- Any decoupled-erf k-scaling claim in v5 is not represented in canonical `3c`; it is present only in `phase2_addon/k_scaling_decoupled_erf_recheck.csv`.

### 3d. POS crossover scaling `k_POS`

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/pos_scaling_points.csv`.
- `outputs/manual_validation.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `119`: POS top-500/top-300 and manual-validation claims.
- Lines `261–263`: POS exponent, CI, p-value.
- Line `516`: POS exponent conclusion restatement.
- Line `534`: manual cross-check on four corpora.
- Line `660`: appendix checklist POS alpha/CI/p-value.
- `J2` is blocked until this experiment exists.

### 4. Analytic seam theory

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/seam_sign_checks.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Lines `299–313`: sign-pattern, tangent-projection, numerical-refit, and second-order seam claims.
- Line `500`: discussion restatement.

Model-canon warning:
- `DRIFT_TRACKER.md` notes the old coupled-logistic seam derivation does not transfer directly to the decoupled erf model. If v5 remains canonical, this experiment should be redesigned before closure.

### 5a. Simulation recovery from fitted smooth models

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/simulation_recovery_per_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Lines `13`, `31`, `143`, `319–321`, `335`, `502`, `662`, `663`: smooth-vs-ZM simulation recovery, high-c block exact match, help rates, Poisson count sampling, and basis-correlation claims.

### 5b. Low-c manifold structure

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/lowc_manifold_per_corpus.csv`.
- `outputs/phase_coordinate_per_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Lines `13`, `147`, `149`, `323`, `335`, `502`, `664`: low-c cosine, full-vs-top100 flips, top-100 match rates, phase-coordinate, and low-c manifold claims.

### 5c. Smooth-parameter control of the low-c manifold

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/parameter_sweep_rows.csv`.
- `outputs/parameter_sweep_correlations.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Lines `13`, `31`, `327–333`, `502`: parameter correlation and winner-count sweep claims.

Model-canon warning:
- `DRIFT_TRACKER.md` notes these width/transition-correlation claims are coupled-family diagnostics, not final decoupled-erf parameter-correlation claims.

### 6. Multilingual extension

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/table3_multilang.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `13`: non-English corpus/language-family side of scope via `J1`.
- Line `153`: multilingual token-count protocol claims.
- Lines `347–365 / Table 3 body`.
- Lines `359`, `361`, `363`, `365`, `498`, `518`: multilingual extension counts, transition fractions, low-c/c≈0 claims, and `29/29` contribution to `J1`.

### 7a. Canonical PMF family

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/splitfit/table4_fourway.csv`.
- `outputs/splitfit/pmf_variant_per_corpus.csv`.
- `outputs/splitfit/aggregate_statistics.csv`.
- `outputs/fullrefit/*` diagnostics.
- README note separating split-fit canonical held-out metrics from full-refit diagnostics.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `163`: 80/20 train/test split protocol.
- Lines `371–407 / Table 4 summary and body`.
- Lines `514`, `536`, `544–546`.
- Appendix lines `644–652`.

### 7b. PMF head-window evaluation

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/head_window_per_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `409`: top-100/top-200/top-1000 head-window winner counts.
- Line `504`: head-restricted held-out NLL discussion.
- Line `536`: top-K cutoff and PMF protocol claims.

### 7c. PMF regularization diagnostic: metadata predictability of `lambda_k`

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/lambda_metadata_summary.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `415`: metadata correlations.
- Line `504`: `|r| < 0.27` discussion restatement.

### 7d. PMF regularization diagnostic: asymmetric gate

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/asymmetric_gate_per_corpus.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `413 / 504`: asymmetric-gate width ratio, counts, median NLL delta, and step-2 help claims.

### 7e. PMF regularization diagnostic: hierarchical pooling on `k`

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/hierk_summary.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `411`: hierarchical `α = 0.505`, `σ = 0.05` lower-bound claim.
- Line `504`: hierarchical pooling does not improve over per-corpus soft-k.

### 8. Bible per-book decomposition

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/bible_per_book.csv`.
- `outputs/table5_bible_summary.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `173`: `66` books.
- Lines `417–432 / Table 5 body`.
- Appendix lines `653–656`.
- Discussion lines `504`, `544–546` involving Bible per-book soft-k and NLL reductions.

### 9. Likelihood frontier beyond soft-k

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/likelihood_frontier_per_corpus.csv`.
- `outputs/hybrid_structure_summary.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Lines `175`, `179`: hybrid/mechanism protocol constants.
- Lines `436–438`, `504`, `657–658`: hybrid, nested, three-regime, and mechanism-penalty claims.

### 10a. Search-depth robustness

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/step10_poly_per_corpus.csv`.
- `outputs/poly_transfer.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `183`: 1000-point normalized grid protocol.
- Line `448`: degree-5 R², step-10/poly decomposition, transfer RMSE values.
- Lines `506 / 538`: discussion restatements of search-depth robustness and transfer claims.

### 10b. Grammar-width robustness

Status: **pending**.

Canonical directory: not present.

Outputs present:
- None in canonical `experiments/`.

Outputs missing:
- `outputs/widened_beam_top10.csv`.
- `outputs/table6_widened_manifold.csv`.
- `outputs/aggregate_statistics.csv`.

Last modified:
- No canonical output.

README status:
- **missing**.

Manuscript claims currently satisfied:
- None from the canonical Phase 2 output structure.

Manuscript claims still unsatisfied:
- Line `337`: 11/14 and 17/17 widened-manifold claims.
- Line `450`: widened beam rows for Shakespeare, War and Peace, Pride and Prejudice.
- Lines `452`, `456–476 / Table 6 body`, `476`, `506`, `538`, `544`, `546`, `659`: widened grammar corpus counts, span R² threshold, rank/gap values, and manifold-membership checklist claim.

## Joined-Output Producers

### J1. Paper-scope rollup

Status: **blocked / pending**.

Canonical directory:
- `joins/J1_paper_scope_rollup/` is missing.

Upstream dependencies:
- `1a`: complete by v3.
- `3a`: complete by v3, but decoupled-erf model-canon integration unresolved if v5 is the target.
- `6`: pending.

Upstream deps complete:
- **No**. `6` is missing.

Producer script:
- Missing: `joins/J1_paper_scope_rollup/scripts/build_paper_scope_rollup.py`.

Joined CSVs present:
- Missing: `joins/J1_paper_scope_rollup/outputs/paper_scope_summary.csv`.

README status:
- **missing**.

Claims currently satisfied:
- None as joined-output claims. `1a` and `3a` upstream cells exist, but no single J1 CSV cell exists.

Claims still unsatisfied:
- Line `13`: full `25 English + 7 non-English + four language families` scope.
- Lines `13 / 31 / 518`: `29/29` Indo-European alphabetic corpus claim.
- Lines `31 / 544`: `32/32` smooth-vs-single-ZM / all-corpus generalization claim.
- Line `534`: full `25 + 7` corpus-scope statement.

### J2. Scaling comparison

Status: **blocked / pending**.

Canonical directory:
- `joins/J2_scaling_comparison/` is missing.

Upstream dependencies:
- `3c`: complete by v3.
- `3d`: pending.

Upstream deps complete:
- **No**. `3d` is missing.

Producer script:
- Missing: `joins/J2_scaling_comparison/scripts/build_scaling_comparison.py`.

Joined CSVs present:
- Missing: `joins/J2_scaling_comparison/outputs/scaling_comparison_summary.csv`.

README status:
- **missing**.

Claims currently satisfied:
- None as joined-output claims. `3c` upstream cells exist, but no paired comparison CSV exists.

Claims still unsatisfied:
- Line `265`: exponent difference, CI, and p-value.
- Any conclusion/discussion sentence that directly compares `k_stat` and `k_POS`.

## Satisfied Manuscript Claims Ready for Builder Refactor

These claims have canonical Phase 2 outputs in `experiments/*/outputs/` and can be wired into a future LaTeX builder after document-version authority is resolved.

From `1a`:
- English corpus count `25`.
- Table 1 body from `outputs/table1_per_corpus.csv`.
- Step-2 top-10 beam from `outputs/table1_step2_beam_top10.csv`.
- Winner family counts: IS `11`, exp `14`, other `0`.
- Current v3 c-band rows: exp-dominant below `65.66394188576379`, IS-dominant above `78.7645572814913`.
- Both Bregman forms present in the step-2 beam on `25/25`.
- Verified replacement for line `197`: terminal enriched-search improves `25/25`, range `3.6978883162152028%` to `25.294536761006274%`, mean `11.023397106886518%`, median `9.798300633636057%`.
- Verified replacement for line `202`: step-2 winner-runner mean gap `0.002009831063671038`; terminal final-step mean gap `0.0007512966946449018`.

From `3a`:
- Smooth beats single ZM on `25/25`.
- Smooth beats hard piecewise on `25/25`.
- Smooth beats step-2 on `25/25`.
- Mean `w = 1.218292636877446`, SD `0.5075054512710272`.
- `w < 0.2` count `4`; corpora: Grimm's Fairy Tales, Aesop's Fables, Dubliners, Critique of Pure Reason.
- Relaxed bounds improves RMSE on `25/25`.
- Smooth-fit per-corpus table and bounds-robustness table are present.

From `3c`:
- `k_stat` forced alpha `0.5214222874234132`.
- CI `[0.4914985326409799, 0.5513460422058466]`.
- p-value against `0.5`: `0.14215982943526212`.
- 25-row `kstat_scaling_points.csv` is present.

## Unsatisfied Manuscript Claims

All claims mapped to missing canonical experiments remain unsatisfied by Phase 2 outputs. Grouped by missing owner:

- `1b`: line `210` city-population control.
- `1c`: search/scoring robustness paragraphs for boundary, guard, and WLS stability.
- `1d`: line `206`, line `534` Euclidean and alternate-generator gap claims.
- `2a`: line `216` function-word ablation c/RMSE claims.
- `2b`: lines `115`, `222`, `224`, `226`, `277`, `534` synthetic and planted-mixture claims.
- `2c`: lines `13`, `31`, `123`, `127`, `273`, `279`, `339`, `492`, `494`, `544` mechanism-absorption claims.
- `3b`: lines `13`, `85`, `101`, `232–243`, `243`, `247`, `257`, `273`, `275`, `492` BIC/model-family claims.
- `3d`: lines `119`, `261–263`, `516`, `534`, `660` POS scaling and manual-validation claims.
- `4`: lines `299–313`, `500` analytic seam-theory claims.
- `5a`: lines `13`, `31`, `143`, `319–321`, `335`, `502`, `662`, `663` simulation-recovery claims.
- `5b`: lines `13`, `147`, `149`, `323`, `335`, `502`, `664` low-c manifold and phase-coordinate claims.
- `5c`: lines `13`, `31`, `327–333`, `502` smooth-parameter control claims.
- `6`: lines `13`, `153`, `347–365`, `359`, `361`, `363`, `365`, `498`, `518` multilingual claims.
- `7a`: lines `163`, `371–407`, `514`, `536`, `544–546`, `644–652` canonical PMF claims.
- `7b`: lines `409`, `504`, `536` head-window PMF claims.
- `7c`: lines `415`, `504` lambda-metadata claims.
- `7d`: lines `413`, `504` asymmetric-gate claims.
- `7e`: lines `411`, `504` hierarchical-pooling claims.
- `8`: lines `173`, `417–432`, `504`, `544–546`, `653–656` Bible per-book claims.
- `9`: lines `175`, `179`, `436–438`, `504`, `657–658` likelihood-frontier claims.
- `10a`: lines `183`, `448`, `506`, `538` search-depth robustness claims.
- `10b`: lines `337`, `450`, `452`, `456–476`, `476`, `506`, `538`, `544`, `546`, `659` grammar-width robustness claims.
- `J1`: lines `13`, `31`, `365`, `518`, `534`, `544` cross-experiment paper-scope rollups.
- `J2`: line `265` direct exponent-comparison claims.

## Add-On / Theorem-Track Outputs Not Yet in Phase 2 Contract

The following work is valuable but currently outside the v3 25-experiment/J1/J2 canonical structure:

- T1 branch test: `phase2_addon/t1_branch_test/`.
- T2 cross-family and coordinate-alignment tests: `phase2_addon/t2_*`.
- T3/T3v2 bifurcation sweeps: `phase2_addon/t3_bifurcation_sweep/` and `phase2_addon/t3_bifurcation_sweep_v2/`.
- Subclaim 3 λ-shift fit: `phase2_addon/t_subclaim3_shift_parameter/`.
- S1/S2 alternative-gate add-ons: `phase2_addon/s1_erf_selection_verification/`, `phase2_addon/s2_*`, plus Windows result folders under `results/`.

If any of these become manuscript-facing, they need either:
- a new consolidated experiment id and README under `experiments/`, or
- an explicit supplement/theorem-track output spec separate from the original Phase 2 plan.

## Recommended Next Planning Decisions

1. Resolve document authority: `MANUSCRIPT_CLAIM_TO_CSV_MAP_v3.md` was requested here, but `v4` exists and declares itself authoritative.
2. Decide whether Phase 2 continues under the original v3 8-parameter coupled-logistic smooth-model contract, or under the newer decoupled 9-parameter erf model used by `MANUSCRIPT_DRAFT_v5.md`.
3. Before closing Section 3.1, add an IS-vs-exp RMSE-gap output to `1a` or a narrowly scoped threshold diagnostic, because T3v2 shows winner flips should be characterized as score crossings rather than just thresholded winner identities.
4. Resume canonical Phase 2 in dependency order only after the model-canon decision:
   - If v3 coupled-logistic remains the target: proceed to `3d`, then `J2`, and `6`, then `J1`.
   - If v5 decoupled-erf is the target: revise `3a`, `3c`, `3b`, `4`, and `5c` specs before running more downstream experiments.

