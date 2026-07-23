# Legacy Coupled-Logistic Seam Theory Archive

## Status

This archive preserves the original coupled-logistic seam-theory experiments. These outputs are legitimate prior research from when the coupled-logistic smooth model was the canonical model. They are not drift, not errors, and not quarantined data.

The parent experiment's canonical decoupled-erf status is blocked only because the equivalent decoupled-erf seam-theory computation has not been run.

## Research question

Under the coupled-logistic smooth model, what sign pattern does the smooth seam predict analytically, and how far do tangent-space projection and second-order correction move that prediction toward the empirical residual?

## Source bundles

- `results/zipf_seam_sign_theory`
- `results/zipf_seam_projection_theory`
- `results/zipf_seam_second_order_theory`

## Outputs

### `outputs/seam_sign_checks.csv`

One row per English corpus. Columns include corpus metadata, local sign vectors, empirical sign vectors, smooth-model sign vectors, projected sign vectors, second-order sign vectors, and match booleans.

### `outputs/aggregate_statistics.csv`

Historical coupled-logistic aggregate metrics using the old claim-map row names.

### `outputs/source_consistency_checks.csv`

Agreement checks for metrics repeated across source bundles.

### `outputs/manifest.json`

Output inventory and provenance.

## Rerun command

```bash
python3 experiments/4_analytic_seam_theory/archive/legacy_coupled_logistic/scripts/run_legacy_coupled_logistic.py
```

This script reads historical saved outputs only. It does not rerun seam derivations or modify historical bundles.

