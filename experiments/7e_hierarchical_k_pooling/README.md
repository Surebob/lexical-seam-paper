# Experiment 7e README

## Experiment id

`7e_hierarchical_k_pooling`

## Research question

Does empirical-Bayes pooling of the PMF seam centre `k` across corpora improve the Seam-Mandelbrot likelihood/mechanism tradeoff relative to per-corpus soft-k regularization?

## Current role

Per the 2026-04-19 scoping decision, the PMF arc (`7a` through `9` and related diagnostics) is queued for a separate Paper 2. In the current Phase 2 structure, 7e preserves the hierarchical-`k` pooling diagnostic and flags a manuscript wording mismatch discovered during migration.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_angle4_hierk`](/Volumes/External2TB/emlexperiment/results/zipf_angle4_hierk)
- [`/Volumes/External2TB/emlexperiment/experiments/7a_canonical_pmf_family/outputs/splitfit/table4_fourway.csv`](/Volumes/External2TB/emlexperiment/experiments/7a_canonical_pmf_family/outputs/splitfit/table4_fourway.csv)

The historical 7e bundle contains Zipf/ZM/MOE/hierarchical-k comparisons but does not contain a soft-k column. This migration therefore adds a cross-experiment comparison against the current 7a canonical v4 soft-k Table 4 values.

## Outputs produced

### `outputs/hierk_summary.csv`

One-row model summary with `model=hierarchical_k`, fitted `alpha`, fitted `sigma`, lower-bound status, and iteration count.

### `outputs/hierk_iteration_history.csv`

Four-row iteration trace for the empirical-Bayes alpha/sigma fit.

### `outputs/hierk_head_window_per_corpus.csv`

One row per corpus per cutoff (`150` rows), preserving the saved head-window held-out NLL and step-2 diagnostics.

### `outputs/hierk_vs_softk_comparison.csv`

One row per corpus comparing full held-out NLL from hierarchical-k against the current 7a canonical v4 soft-k value.

### `outputs/aggregate_statistics.csv`

Claim-facing rows including `alpha`, `sigma`, `sigma_at_lower_bound`, `hierk_vs_softk_summary`, and the underlying soft-k comparison counts/deltas.

## Manuscript claims fed

- Section 3.11: `log k_i ~ Normal(alpha log V_i, sigma^2)`.
- Section 3.11: `alpha = 0.5049174340956519`, rendered as `0.505`.
- Section 3.11: `sigma = 0.05`, pinned at the lower bound.
- Discussion restatement: hierarchical pooling compared to per-corpus soft-k.

## Methods

No hierarchical-k fitting is rerun. This migration reads the saved `zipf_angle4_hierk` outputs and normalizes them into canonical CSVs. The comparison to soft-k is a direct join against 7a's canonical Table 4 soft-k column; no model fitting or optimization is performed.

## Rerun command

```bash
python3 experiments/7e_hierarchical_k_pooling/scripts/run_7e.py
```

This command regenerates the consolidated CSVs from saved historical bundles and the already-migrated 7a canonical Table 4.

## Verification mapping

- `hierk_summary.csv` must contain one row with `model=hierarchical_k`.
- `alpha` must be `0.5049174340956519`.
- `sigma` must be `0.05`.
- `sigma_at_lower_bound` must be `True`.
- `hierk_vs_softk_comparison.csv` must contain `25` rows.

## AUDIT

Manuscript v5.1 says: "Empirical-Bayes hierarchical pooling of k across corpora does not improve over per-corpus soft-k regularization." The saved 7e source does not directly contain soft-k values. When joined to the current 7a canonical v4 soft-k Table 4, hierarchical-k beats soft-k on `20/25` corpora with median `hierk_minus_softk = -0.0009123035534583934`.

This is a data-claim mismatch. The canonical migration preserves the data-level comparison and flags the manuscript wording for revision rather than forcing the output to match the prose. A more conservative Paper 2 wording would distinguish statistical held-out improvement from the broader likelihood/mechanism tradeoff and from source-convention sensitivity in 7a.

The historical hierarchical-k work remains legitimate prior PMF research. If Paper 2 uses this result, it should rerun or restate the comparison under the final Paper 2 PMF canonical-source convention.
