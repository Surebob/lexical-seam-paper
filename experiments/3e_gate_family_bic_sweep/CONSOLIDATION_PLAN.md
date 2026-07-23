# Experiment 3e Consolidation Plan

## Purpose

Experiment `3e_gate_family_bic_sweep` migrates the S2 v3 decoupled five-gate sweep into the canonical Phase 2 experiment structure. This closes the gate-family data-access gap noted in the GPT-5 Pro audit finding G.4 / proposed experiment E1.

## Source bundles

### Real-data gate-family sweep

- Source directory: `results/s2_v3_windows_full_outputs_2026-04-18/`
- Key files:
  - `s2_v3_per_fit_results.csv`
  - `s2_v3_per_corpus_results.csv`
  - `s2_v3_per_start_dispersion.csv`
  - `s2_v3_tanh_calibration.csv`
  - `s2_v3_aggregate_statistics.csv`
  - `s2_v3_runtime_log.txt`

This bundle contains the 25-corpus, five-gate decoupled sweep with logistic, tanh, erf, algebraic, and arctan gates. Tanh is a calibration control and is excluded from independent-gate winner counts and spreads.

### Synthetic logistic-gate recovery check

- Source directory: `results/s2_synthetic_recovery_outputs_workers32_2026-04-18/outputs/`
- Key files:
  - `synthetic_gate_recovery_per_corpus.csv`
  - `synthetic_gate_recovery_per_fit.csv`
  - `synthetic_gate_recovery_per_start_dispersion.csv`
  - `synthetic_gate_recovery_aggregate_statistics.csv`
  - `synthetic_gate_recovery_status.json`

This bundle tests whether the same fitting protocol recovers the true gate on logistic-generated synthetic corpora. It is a validation sidecar for the empirical erf preference.

### Producer-script provenance

The result bundles do not contain their producer scripts directly. The source scripts are preserved from the working tree:

- `phase2_addon/s2_decoupled_v3/run_s2_v3_full.py`
- `s2_v3_windows_port/run_s2_v3_windows.py`
- `phase2_addon/s2_logistic_synthetic_gate_recovery/run_synthetic_gate_recovery.py`

The canonical `run_3e.py` does not rerun fitting; it regenerates canonical CSVs from saved outputs.

## Canonical outputs

- `outputs/per_fit_results.csv`: normalized one-row-per `(corpus, gate_family)` real-data fit table.
- `outputs/per_corpus_results.csv`: one-row-per-corpus winner/spread/calibration summary.
- `outputs/aggregate_statistics.csv`: claim-facing gate-family aggregate rows.
- `outputs/synthetic_recovery.csv`: one-row-per synthetic corpus logistic-recovery summary.
- `outputs/synthetic_recovery_aggregate.csv`: claim-facing synthetic-recovery aggregates.
- `outputs/manifest.json`: output schemas and source provenance.

## Precedence and conflict rules

- Real-data BIC, RMSE, `k`, `w_gate`, `w_tail`, bounds, runtime, and optimizer diagnostics come from `s2_v3_windows_full_outputs_2026-04-18`.
- Synthetic recovery results come from `s2_synthetic_recovery_outputs_workers32_2026-04-18`.
- There is no row-level conflict between the real-data and synthetic bundles because they feed disjoint canonical outputs.
- Logistic regime coefficients `a1,b1,c1,a2,b2,c2` are recoverable from `windows_s2_decoupled_v3_results_2026-04-18/outputs/s2_v3_logistic_local_params.csv`; non-logistic regime coefficients were not saved in the Windows full-sweep CSVs. The canonical schema includes the columns, leaves unavailable cells blank, and flags this as an audit limitation.

## Known limitations

- Non-logistic fitted regime coefficients are not reconstructible without rerunning the expensive fits. The BIC/RMSE/winner claims do not depend on these coefficients.
- Synthetic recovery has tanh/logistic BIC calibration marked false because the synthetic data are nearly exactly logistic and tiny absolute RMSE differences produce huge BIC differences. This does not affect the requested validation rows: logistic is recovered on 25/25 independent-gate comparisons and erf never beats the true generator.
