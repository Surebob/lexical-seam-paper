# Experiment 7d README

## Experiment id

`7d_asymmetric_gate_diagnostic`

## Research question

Does allowing the Seam-Mandelbrot PMF seam to have asymmetric left/right widths materially improve held-out likelihood or mechanism capture over the symmetric soft-k gate?

## Current role

Per the 2026-04-19 scoping decision, the PMF arc (`7a` through `9` and related diagnostics) is queued for a separate Paper 2. In the current Phase 2 structure, 7d preserves the asymmetric-gate diagnostic used by manuscript v5.1.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_angle3_asymmetric_gate_splitfit`](/Volumes/External2TB/emlexperiment/results/zipf_angle3_asymmetric_gate_splitfit)

## Outputs produced

### `outputs/asymmetric_gate_per_corpus.csv`

One row per English corpus (`25` rows), including selected `lambda_k`, fitted `w_left`, `w_right`, width ratio, width-imbalance flag, full held-out average NLL values, asymmetric-vs-baseline deltas, and full residual step-2 gain/help.

### `outputs/aggregate_statistics.csv`

Claim-facing aggregate rows including `median_width_ratio`, `material_asymmetry_count`, `asymmetric_beats_softk_count`, `median_asymmetric_minus_softk`, and `asymmetric_step2_help_count`.

### `outputs/manifest.json`

Machine-readable inventory of source bundles and output schemas.

## Manuscript claims fed

- Section 3.11: median `w_right / w_left = 1.00`.
- Section 3.11: `5/25` corpora have material asymmetry above `1.5x`.
- Section 3.11: asymmetric gate beats symmetric soft-k on full held-out NLL on `22/25` corpora.
- Section 3.11: median asymmetric minus soft-k held-out delta is `-0.000272496315`.
- Section 3.11: asymmetric-gate residual step-2 helps on `6/25` corpora.

## Methods

No asymmetric-gate fitting is rerun. This migration reads the saved split-fit asymmetric-gate bundle, normalizes the full-cutoff per-corpus table, recomputes only simple display fields such as `width_imbalance_ratio`, and archives the historical source files.

## Rerun command

```bash
python3 experiments/7d_asymmetric_gate_diagnostic/scripts/run_7d.py
```

This command regenerates the consolidated CSVs from saved historical bundles.

## Verification mapping

- `asymmetric_gate_per_corpus.csv` must contain `25` rows.
- `median_width_ratio` must be `1.0`.
- `material_asymmetry_count` must be `5`.
- `asymmetric_beats_softk_count` must be `22`.
- `median_asymmetric_minus_softk` must be `-0.00027249631523407203`.
- `asymmetric_step2_help_count` must be `6`.

## AUDIT

This is a single-source migration, so no multi-bundle conflict resolution was needed. The asymmetric model is a valid historical diagnostic for the PMF family. It is not a Paper 1 canonical mechanism result, and Paper 2 should decide whether to rerun it under whatever PMF source convention is chosen for that manuscript.
