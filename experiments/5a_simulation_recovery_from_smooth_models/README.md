# Experiment 5a README

## Experiment id

`5a_simulation_recovery_from_smooth_models`

## Research question

Does synthetic data sampled from fitted smooth rank-frequency models reproduce the empirical symbolic residual signal?

## Why this is one experiment

This experiment asks one causal-sufficiency question: if the fitted smooth model generated the corpus, would the symbolic residual search recover the same Bregman-family correction seen empirically? The single-ZM synthetic control is included because it is the negative control for the same question.

## Upstream dependencies

- Historical source bundle: [`/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery`](/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery)
- Main source table: [`/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery/simulation_recovery_table.csv`](/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery/simulation_recovery_table.csv)
- Summary source: [`/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_simulation_recovery/summary.json)

## Outputs produced

### `outputs/simulation_recovery_per_corpus.csv`

One row per English corpus.

Key columns:
- `slug`
- `corpus`
- `empirical_winner`
- `empirical_zm_c`
- `smooth_modal_winner`
- `smooth_exact_match_rate`
- `smooth_help_rate`
- `single_zm_control_modal_winner`
- `single_zm_control_exact_match_rate`
- `single_zm_control_help_rate`

### `outputs/aggregate_statistics.csv`

Single-cell aggregates used by manuscript prose.

Required rows include:
- `single_zm_control_step2_help_count`
- `single_zm_control_step2_help_rate`
- `smooth_step2_help_rate`
- `smooth_exact_winner_match_rate`
- `single_zm_control_exact_winner_match_rate`
- `basis_corr_linear_smooth`
- `basis_corr_quadratic_smooth`
- `basis_corr_cubic_smooth`
- `basis_corr_linear_zm`
- `basis_corr_quadratic_zm`
- `basis_corr_cubic_zm`
- `high_c_total_count`
- `high_c_block_size`
- `high_c_exact_match_count`
- `high_c_exact_match_rate_smooth`
- `high_c_help_rate_smooth`
- `high_c_help_rate_zm_control`

## Manuscript claims fed

Mapped claims include:
- single-ZM synthetic controls never produce helpful step-2 correction (`0/25`)
- smooth-generated synthetic corpora reproduce the empirical IS winner on the high-c block (`11/11`)
- smooth-generated high-c block step-2 help rate is `1.00`
- replicate-level exact-match/help-rate and head-basis correlation claims in Section 3.9 and conclusion restatements

## Methods

This migration reads the saved historical simulation-recovery table and summary. It does not regenerate synthetic corpora or rerun symbolic search.

Inference documentation: no confidence intervals or p-values are computed in this migration. All reported quantities are counts, proportions, medians, and correlations saved by the historical bundle.

## Rerun command

```bash
python3 experiments/5a_simulation_recovery_from_smooth_models/scripts/run_5a.py
```

## Verification mapping

- `simulation_recovery_per_corpus.csv` must contain `25` rows.
- `aggregate_statistics.csv` must contain the high-c block count `11`, high-c smooth exact-match rate `1.0`, and single-ZM-control help count `0`.
- Source files are copied under `archive/provenance/`.

## AUDIT

This is a historical smooth-generation experiment. The source bundle predates the decoupled-erf model canon, so if the paper requires a strict decoupled-erf simulation-recovery version, that should be run as a new experiment rather than inferred during migration. The historical output remains legitimate prior research and is preserved as such.

