# Experiment 3e README

## Experiment id

`3e_gate_family_bic_sweep`

## Research question

Under the decoupled 9-parameter two-regime ZM model, is the empirical transition specifically erf-shaped, or would any smooth bounded monotone gate fit equally well?

## Why this is one experiment

The real-data five-gate sweep and the logistic-generated synthetic recovery check answer one gate-family specificity question. The real-data sweep tests which gate is empirically preferred. The synthetic recovery sidecar tests whether the fitting protocol would recover a known non-erf generator rather than spuriously preferring erf.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/s2_v3_windows_full_outputs_2026-04-18/`
- `/Volumes/External2TB/emlexperiment/results/s2_synthetic_recovery_outputs_workers32_2026-04-18/outputs/`
- `/Volumes/External2TB/emlexperiment/results/windows_s2_decoupled_v3_results_2026-04-18/outputs/s2_v3_logistic_local_params.csv`

## Inputs

- `s2_v3_per_fit_results.csv`
- `s2_v3_per_corpus_results.csv`
- `s2_v3_per_start_dispersion.csv`
- `s2_v3_tanh_calibration.csv`
- `s2_v3_aggregate_statistics.csv`
- `synthetic_gate_recovery_per_corpus.csv`
- `synthetic_gate_recovery_per_fit.csv`
- `synthetic_gate_recovery_aggregate_statistics.csv`
- `synthetic_gate_recovery_status.json`

## Outputs produced

### `outputs/per_fit_results.csv`

One row per empirical `(corpus, gate_family)` fit.

Columns include:
- `slug`
- `corpus`
- `gate_family`
- `seed`
- `a1`, `b1`, `c1`, `a2`, `b2`, `c2`
- `k`, `w_gate`, `w_tail`
- `BIC`
- `final_rmse`
- bound-hit flags
- optimizer diagnostics
- per-start dispersion summary
- `parameter_source`

### `outputs/per_corpus_results.csv`

One row per English corpus.

Columns include:
- `slug`
- `corpus`
- `winning_gate`
- `winner_bic`
- `independent_gate_spread`
- `tanh_calibration_pass`
- per-independent-gate BICs

### `outputs/aggregate_statistics.csv`

Claim-facing empirical gate-family aggregates:
- `erf_wins_count`
- `arctan_wins_count`
- `logistic_wins_count`
- `algebraic_wins_count`
- `tanh_calibration_pass_count`
- `median_independent_gate_bic_spread`
- `mean_independent_gate_bic_spread`
- `gates_indistinguishable_strict_count`
- `gates_indistinguishable_positive_count`
- `gates_indistinguishable_strong_count`

### `outputs/synthetic_recovery.csv`

One row per logistic-generated synthetic corpus.

Columns include:
- `synthetic_corpus_id`
- `true_gate`
- `recovered_gate_by_fitter`
- `recovery_match`
- `erf_beats_true`
- `logistic_bic`
- `erf_bic`

### `outputs/synthetic_recovery_aggregate.csv`

Claim-facing synthetic-recovery aggregates:
- `logistic_recovery_count`
- `erf_beats_true_count`
- `synthetic_completed_corpus_count`
- `synthetic_completed_fit_count`

### `outputs/manifest.json`

Machine-readable output inventory, source provenance, and audit limitations.

## Manuscript claims fed

- Abstract gate-family sentence: erf wins `24/25`; arctan wins only Dubliners; median BIC spread `677`; mean spread `865`; logistic synthetic recovery `25/25`; erf beats true generator `0/25`.
- Introduction pillar `c″`: gate-family specificity and synthetic recovery validation.
- Section 3.5 Table 2: decoupled gate-family BIC winner counts.
- Section 3.5 prose: tanh/logistic calibration and independent-gate spread.
- Conclusion paragraph 2: gate-family specificity restatement.
- Appendix verification checklist gate-family BIC wins.

## Methods

No fitting is rerun during migration. The canonical outputs are regenerated from the saved Windows S2 v3 full sweep and the saved logistic synthetic recovery bundle.

The model is the decoupled 9-parameter two-regime ZM:

`a1, b1, c1, a2, b2, c2, k, w_gate, w_tail`

The five gates are logistic, tanh, erf, algebraic, and arctan. Tanh is a calibration control because it is mathematically equivalent to logistic under `w_tanh = 2 * w_logistic`. Independent-gate winner counts and BIC spreads exclude tanh.

## AUDIT

The Mac prototype writer (`run_s2_v3_full.py`) would have saved all-gate coefficients: line 754 defines `gate_param_fieldnames` with fields such as `logistic_a1`, `erf_a1`, `algebraic_a1`, and `arctan_a1`; line 791 writes `s2_v3_gate_params_per_corpus.csv`; and line 853 populates those fields from `fits[gate]`. The corresponding Mac output directory, `results/s2_v3_decoupled_five_gate_2026-04-18/`, is absent, consistent with an interrupted run.

The completed Windows port (`run_s2_v3_windows.py`) had coefficients in memory via `decoupled_smooth_model.py`, but wrote only BIC/RMSE/k/w diagnostics to `s2_v3_per_fit_results.csv`: line 117 defines the saved fields, line 216 appends rows with the reduced schema, and line 408 writes the final outputs. Therefore non-logistic `a1,b1,c1,a2,b2,c2` are unavailable from saved outputs. Recovering them requires either a fresh rerun of the Mac script or a patched Windows script and rerun of the full five-gate sweep.

The canonical `per_fit_results.csv` includes coefficient columns to preserve the planned schema, but non-logistic cells are blank and `parameter_source` states that the coefficients were not saved. Logistic coefficients are filled from `s2_v3_logistic_local_params.csv`, which is a companion decoupled-logistic output and agrees with the full sweep to numerical precision on BIC/RMSE. Manuscript-cited claims (BIC values, winner counts, spreads, and synthetic recovery) do not depend on these coefficients.

This limitation does not affect manuscript-cited gate-family aggregates, which depend on BIC, RMSE, winner identity, calibration, bounds, and synthetic recovery rows that are present in the source bundles.

The synthetic recovery bundle marks tanh/logistic calibration as failed because logistic-generated synthetic data are fit almost exactly and tiny numerical RMSE differences produce large BIC swings. The validation claim in the manuscript does not use synthetic tanh calibration; it uses independent-gate recovery: logistic wins `25/25`, erf wins `0/25`, and erf never beats the true logistic generator.

## Canonical claim mapping

- Gate winner counts map to `outputs/aggregate_statistics.csv`.
- Independent-gate BIC spread maps to `outputs/aggregate_statistics.csv`.
- Per-corpus Table 2 body maps to `outputs/per_corpus_results.csv`.
- Supplementary full gate-fit diagnostics map to `outputs/per_fit_results.csv`.
- Logistic synthetic recovery maps to `outputs/synthetic_recovery.csv` and `outputs/synthetic_recovery_aggregate.csv`.

## Rerun command

```bash
python3 experiments/3e_gate_family_bic_sweep/scripts/run_3e.py \
  --output-dir experiments/3e_gate_family_bic_sweep/outputs
```

This command consolidates saved outputs only. It does not rerun the expensive five-gate fits.

## Verification mapping

- `per_fit_results.csv` must have `125` rows.
- `per_corpus_results.csv` must have `25` rows.
- `synthetic_recovery.csv` must have `25` rows.
- `erf_wins_count = 24`.
- `arctan_wins_count = 1`.
- `logistic_wins_count = 0`.
- `algebraic_wins_count = 0`.
- `tanh_calibration_pass_count = 25`.
- `median_independent_gate_bic_spread ≈ 677.37`.
- `mean_independent_gate_bic_spread ≈ 864.66`.
- `logistic_recovery_count = 25`.
- `erf_beats_true_count = 0`.
