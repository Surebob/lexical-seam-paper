# Legacy Coupled-Logistic Smooth-Parameter Sweep Archive

## Status

This archive preserves the historical coupled-logistic smooth-parameter sweep from `zipf_smooth_parameter_sweep`. It was the prior canonical parameter-control analysis before the project moved to decoupled-erf.

## Outputs

- `outputs/parameter_sweep_rows.csv`: all 216 synthetic sweep configurations
- `outputs/parameter_sweep_correlations.csv`: claim-facing parameter correlations and standardized regression weights
- `outputs/aggregate_statistics.csv`: synthetic winner counts and other aggregate rows
- `outputs/manifest.json`: provenance and row counts

## Rerun command

```bash
python3 experiments/5c_smooth_parameter_control/archive/legacy_coupled_logistic/scripts/run_legacy_coupled_logistic.py
```

