# J2 README

## Join id

`J2_scaling_comparison`

## Research question

How does the statistical transition-centre scaling exponent compare with the POS crossover scaling exponent under the current decoupled-erf model canon?

## Why this is a joined output

The manuscript cites a numerical comparison between a smooth-model `k` scaling exponent and a POS crossover exponent. Those estimates come from separate producers. J2 emits the comparison as a single CSV so the LaTeX builder does not compute the join.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/experiments/3e_gate_family_bic_sweep/outputs/per_fit_results.csv`
- `/Volumes/External2TB/emlexperiment/experiments/3e_gate_family_bic_sweep/outputs/per_corpus_results.csv`
- `/Volumes/External2TB/emlexperiment/experiments/3d_pos_crossover_scaling/outputs/aggregate_statistics.csv`

## Inputs

J2 reads canonical decoupled-erf `k` values from experiment `3e` and joins them to canonical experiment `3e` vocabulary sizes. It then computes the OLS slope in `log(k_erf) = alpha + beta log(V)` and compares that statistical transition exponent with the POS crossover exponent from experiment `3d`.

## Outputs produced

### `outputs/scaling_comparison_summary.csv`

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Rows:
- `alpha_kstat`
- `alpha_pos`
- `alpha_difference_pos_minus_kstat`
- `alpha_difference_ci_low`
- `alpha_difference_ci_high`
- `alpha_difference_pvalue`

### `outputs/manifest.json`

Machine-readable output inventory, upstream provenance, and inference-method metadata.

## Manuscript lines fed

- line `265`: exponent-difference comparison and uncertainty.
- Any v5.1 restatement comparing decoupled-erf statistical k scaling against POS crossover scaling.

## Methods

Canonical statistical exponent chain: experiment `3e` erf `k` values from `outputs/per_fit_results.csv` plus experiment `3e` vocabulary sizes from `outputs/per_corpus_results.csv` → J2 OLS fit `log(k_erf) = alpha + beta log(V)` → joined with experiment `3d` POS crossover exponent → `outputs/scaling_comparison_summary.csv`.

POS exponent: forced-alpha estimate from experiment `3d`.

CI method: independent/Welch propagation of uncertainty. The decoupled-erf beta uses its OLS standard error and `df = 23`. The POS alpha standard error is recovered from the saved 95% Student-`t` interval with `df = 24`. The difference interval uses a Welch-Satterthwaite degrees-of-freedom approximation.

P-value method: two-sided Welch `t` test for `alpha_pos - alpha_kstat = 0` using the same propagated standard error and Welch degrees of freedom.

## AUDIT

The claim map v4 was written before the decoupled-erf model became canonical and describes a paired comparison between the old coupled-logistic `3c` exponent and POS exponent. This J2 output follows the user's 2026-04-19 model-canon adjustment by using canonical experiment `3e` as the authoritative erf-k source. Manuscript line `265` and any old `0.024`, CI, or p-value text from the coupled-logistic comparison must be revised against this CSV.

## Canonical claim mapping

- `alpha_kstat` maps to the decoupled-erf statistical scaling exponent.
- `alpha_pos` maps to the POS crossover exponent.
- `alpha_difference_pos_minus_kstat`, `alpha_difference_ci_low`, `alpha_difference_ci_high`, and `alpha_difference_pvalue` map to the joined comparison claim.

## Rerun command

```bash
python3 joins/J2_scaling_comparison/scripts/build_scaling_comparison.py \
  --output-dir joins/J2_scaling_comparison/outputs
```

## Verification mapping

- `scaling_comparison_summary.csv` must contain exactly the six metric rows listed above.
- `alpha_kstat` must match the OLS slope recomputed from experiment `3e` rows with `gate_family=erf`.
- `alpha_pos` must match `forced_alpha` in experiment `3d`.
