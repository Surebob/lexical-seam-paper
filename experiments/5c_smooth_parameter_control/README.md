# Experiment 5c README

## Experiment id

`5c_smooth_parameter_control`

## Research question

Which smooth-model parameters control the orientation of the low-c manifold and the threshold at which winner identity flips under head-weighted scoring?

## Current canonical status

**BLOCKED for canonical decoupled-erf output.**

The current model canon is decoupled-erf, with separate `w_gate` and `w_tail`. The available historical source bundle, `zipf_smooth_parameter_sweep`, varies the old coupled-logistic/sigmoid-width parameterization. That sweep is real prior research and remains scientifically useful as a prior-canon comparison, but it cannot be presented as the current decoupled-erf canonical parameter-control experiment.

## Outputs produced

Canonical:

- `BLOCKED.md`
- `outputs/manifest.json`

Legacy archive:

- `archive/legacy_coupled_logistic/README.md`
- `archive/legacy_coupled_logistic/scripts/run_legacy_coupled_logistic.py`
- `archive/legacy_coupled_logistic/outputs/parameter_sweep_rows.csv`
- `archive/legacy_coupled_logistic/outputs/parameter_sweep_correlations.csv`
- `archive/legacy_coupled_logistic/outputs/aggregate_statistics.csv`
- `archive/legacy_coupled_logistic/outputs/manifest.json`

## Manuscript claims fed

The old claim-map rows for `r = -0.540`, `r = +0.554`, and the synthetic full/top-100 winner counts are satisfied only by the legacy coupled-logistic archive, not by canonical decoupled-erf outputs.

## Methods

No fresh computation is run. The canonical experiment records missing decoupled-erf outputs. The archive producer reads the historical `zipf_smooth_parameter_sweep/summary.json` and emits normalized CSVs for audit/reproduction of the prior-canon result.

Inference documentation: no confidence intervals or p-values are computed in the historical sweep. It reports descriptive winner counts and Pearson/standardized-regression correlations saved by the historical bundle.

## Rerun command

```bash
python3 experiments/5c_smooth_parameter_control/run_experiment.py
```

## Verification mapping

- Canonical `outputs/manifest.json` must report `blocked_missing_decoupled_erf_parameter_sweep`.
- Legacy `parameter_sweep_rows.csv` must contain `216` rows.
- Legacy `parameter_sweep_correlations.csv` must contain the parameter rows `transition_fraction`, `sigmoid_width`, and `tail_head_slope_contrast`.

## AUDIT

Archive means "superseded by later canonical choice," not "bad data." The historical sweep is preserved as a first-class model-canon-comparison record. Current manuscript uses of 5c-style parameter-control claims need revision or a fresh decoupled-erf parameter sweep.

