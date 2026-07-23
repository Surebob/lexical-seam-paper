# Drift Tracker

## 2026-04-17

- **Line 197 range `0.15% to 2.0%`**
  - Status: `PROSE-ONLY`
  - Finding: no computational source was found via exhaustive repo search or prior session-history reconstruction. The phrase appears in manuscript prose and later provenance/mapping documents, but not in a saved experiment aggregate.
  - Verified replacement options:
    - terminal semantics: `3.7%` to `25.3%`, mean `11.0%`
    - step-2 semantics: step-2 helps on `12/25` corpora with gains `0.3%` to `2.5%`; on the remaining `13/25` corpora step-2 slightly hurts fit
  - Note: choosing between these changes the framing of Section 3.1 and must be resolved in Phase 3 manuscript edits.

- **Line 197 `25 of 25 corpora`**
  - Status: `VALID under terminal semantics`
  - Finding: the full enriched symbolic search beats the single-ZM baseline on `25/25` corpora when evaluated at the terminal search step using `zm_search.best.composite_rmse`.
  - Clarification for Phase 3 prose: “the full enriched symbolic search … beats the single-ZM baseline on 25 of 25 corpora at the terminal search step.”

- **Line 202 `0.003 mean RMSE gap`**
  - Status: `PROSE-ONLY`
  - Finding: no computational source was found for `0.003` as written.
  - Verified replacement options:
    - step-2 winner-runner mean gap: `0.00201`
    - terminal final-step winner-runner mean gap: `0.00075`
  - Note: the replacement should follow whichever semantics Section 3.1 ultimately adopts in Phase 3.

- **Line 245 `w < 0.2` corpus list**
  - Status: `DRIFT`
  - Finding: experiment `3a` confirms the manuscript currently names the wrong fourth corpus in the sharp-transition list.
  - Verified canonical list: `Grimm's Fairy Tales; Aesop's Fables; Dubliners; Critique of Pure Reason`
  - Correction needed in Phase 3: replace `Principia Ethica` with `Critique of Pure Reason`.

## Validated Claims

- **Experiment 3a clean validations**
  - `smooth_beats_single_zm_count = 25`
  - `smooth_beats_piecewise_count = 25`
  - `smooth_beats_step2_count = 25`
  - `relaxed_bounds_improves_rmse_count = 25`
  - Note: these Section 3.5/3.6 smooth-model dominance claims hold exactly under the canonical Phase 2 rerun.

## 2026-04-18 Phase 3 v5 Manuscript Update

- **Section 3.1 terminal enriched-search improvement**
  - Updated claim: terminal enriched symbolic search improves over single ZM on `25/25` English corpora with percentage RMSE improvement range `3.70%` to `25.29%` and mean `11.02%`.
  - Source: `experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv`
  - Source rows: `terminal_all_improve_count`, `terminal_improvement_pct_min`, `terminal_improvement_pct_max`, `terminal_improvement_pct_mean`.
  - Justification: replaces prose-only `0.15% to 2.0%` range with verified terminal-search metric.

- **Section 3.1 step-2 diagnostic improvement**
  - Updated claim: step-2 winner improves over single ZM on `12/25` corpora; positive-only gains range `0.33%` to `2.50%`; all-corpus median is `-0.23%`.
  - Source: `experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv`
  - Source rows: `step2_all_improve_count`, `step2_positive_improvement_pct_min`, `step2_positive_improvement_pct_max`, `step2_improvement_pct_median`.
  - Justification: clarifies that step 2 is an interpretable diagnostic layer, not the terminal fit-improvement claim.

- **Section 3.1 winner-runner gaps**
  - Updated claim: step-2 mean winner-runner gap is `0.00201`; terminal final-step mean gap is `0.00075`.
  - Source: `experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv`
  - Source rows: `step2_winner_runner_gap_mean`, `terminal_winner_runner_gap_mean`.
  - Justification: replaces prose-only `0.003` gap.

- **Section 3.5 smooth dominance counts and corrected sharp-transition corpus list**
  - Updated claim: smooth two-regime fit beats single ZM, hard/piecewise alternatives, and ZM-plus-step-2 diagnostics on `25/25`; relaxed bounds improve RMSE on `25/25`; the `w < 0.2` coupled-width list is `Grimm's Fairy Tales; Aesop's Fables; Dubliners; Critique of Pure Reason`.
  - Source: `experiments/3a_smooth_two_regime_fits/outputs/aggregate_statistics.csv`
  - Source rows: `smooth_beats_single_zm_count`, `smooth_beats_piecewise_count`, `smooth_beats_step2_count`, `relaxed_bounds_improves_rmse_count`, `w_lt_0_2_corpora`.
  - Justification: preserves validated 3a dominance claims and corrects the prior Principia Ethica drift.

- **Section 3.5 gate-family specificity**
  - New claim: under decoupled 9-parameter fitting, `erf` wins `24/25`, `arctan` wins `1/25`, `logistic` and `algebraic` win `0/25`; tanh/logistic calibration passes `25/25`; median independent-gate BIC spread is `677.37`, mean spread is `864.66`.
  - Source: `results/s2_v3_windows_full_outputs_2026-04-18/s2_v3_aggregate_statistics.csv`
  - Source rows: `erf_bic_wins`, `arctan_bic_wins`, `logistic_bic_wins`, `algebraic_bic_wins`, `tanh_calibration_pass_count`, `median_bic_spread`, `mean_bic_spread`.
  - Per-corpus table source: `results/s2_v3_windows_full_outputs_2026-04-18/s2_v3_per_corpus_results.csv`.

- **Section 3.5 synthetic gate-recovery validation**
  - New claim: on logistic-generated synthetic data, logistic is recovered on `25/25`; erf wins `0/25`; erf never beats the true logistic generator.
  - Source: `results/s2_synthetic_recovery_outputs_workers32_2026-04-18/outputs/synthetic_gate_recovery_aggregate_statistics.csv`
  - Source rows: `logistic_bic_wins`, `erf_bic_wins`, `logistic_recovered_count`, `erf_beats_logistic_count`.
  - Justification: validates that empirical erf preference is not a generic five-gate fitting artifact.

- **Section 3.6 k-scaling weakened**
  - Updated claim: decoupled erf k-scaling OLS slope is `β = 0.7584`, 95% CI `[0.3457, 1.1711]`, `R² = 0.3858`.
  - Source: `phase2_addon/k_scaling_decoupled_erf_recheck.csv`
  - Source rows: `ols_beta`, `ols_r2`.
  - Justification: replaces the prior tight `0.5214 [0.4915, 0.5513]` coupled-logistic statistical scaling claim with the decoupled-erf recheck.

- **Section 3.8 erf seam derivation**
  - Updated claim: the product-rule expansion remains exact, but the coupled-logistic `(−, −, −)` sign-count validation is not transferred to the decoupled erf model.
  - Source: analytic derivation from `σ_erf(z) = [1 - erf(z / w_gate)] / 2`; implementation reference `s2_v3_windows_port/shared/gate_functions.py`.
  - Justification: for `τ = 1 - σ_erf`, `τ(0)=1/2`, `τ′(0)=1/(√π w_gate)`, `τ″(0)=0`, `τ‴(0)=-2/(√π w_gate^3)`, changing the third-order coefficient structure.

- **Section 3.9 parameter-correlation caveat**
  - Updated claim: low-c parameter correlations involving transition width are described as coupled-family diagnostics, not final decoupled-erf parameter-correlation claims.
  - Source: existing coupled-family smooth parameter sweep; no decoupled-erf parameter-correlation rerun was found among provided canonical outputs.
  - Justification: prevents overclaiming w_gate/w_tail interpretation from the older single-width sweep.
