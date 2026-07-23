# Experiment 7a README

## Experiment id

`7a_canonical_pmf_family`

## Research question

Within the Seam-Mandelbrot PMF family, what is the canonical held-out tradeoff among the free, fixed-`k`, fixed-`k,w`, soft-`k`, and soft-`k,w` variants under the train/test protocol used in the paper?

## Current role

Per the 2026-04-19 scoping decision, the PMF arc (`7a` through `9` and related diagnostics) is queued for a separate Paper 2. In the current Phase 2 structure, 7a's role is research-record preservation and provenance cleanup, not Paper 1 support. Paper 2 should revisit canonical-source decisions when that manuscript is drafted.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_pmf`](/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_pmf)
- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_regularized`](/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_regularized)
- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_softk`](/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_softk)
- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_softk_splitfit`](/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_softk_splitfit)
- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_softkw`](/Volumes/External2TB/emlexperiment/results/zipf_seam_mandelbrot_softkw)
- [`/Volumes/External2TB/emlexperiment/results/zipf_v4_verification/table_a_fourway_pmf.csv`](/Volumes/External2TB/emlexperiment/results/zipf_v4_verification/table_a_fourway_pmf.csv)

## Outputs produced

### `outputs/splitfit/table4_fourway.csv`

Canonical manuscript-facing Table 4 body. This is copied from `zipf_v4_verification/table_a_fourway_pmf.csv` with normalized column names and recomputed `softk_minus_zipf`.

### `outputs/splitfit/table4_provenance.csv`

One row per corpus documenting which source each Table 4 column traces to and whether it matches upstream PMF/soft-k sources.

### `outputs/splitfit/softk_source_diagnostic.csv`

One row per corpus comparing:

- v4 verification soft-`k`
- legacy soft-`k`
- soft-k,w repeated soft-`k`
- patched `softk_splitfit`

Three rows differ between v4/legacy and patched splitfit because the optimizer lands in different local optima under different warm-start behavior.

### `outputs/splitfit/pmf_variant_per_corpus.csv`

Union table for free, fixed-`k`, fixed-`k,w`, soft-`k`, and soft-`k,w` variants.

### `outputs/splitfit/aggregate_statistics.csv`

Claim-facing aggregates including Table 4 winner counts, soft-`k` median deltas, regularization variant summaries, lambda distribution, and step-2 diagnostic counts.

### `outputs/fullrefit/*`

Full-refit diagnostics only. These are not canonical held-out metrics.

## Manuscript claims fed

This experiment preserves the PMF Table 4 and associated verification-checklist rows:

- soft-k beats MOE on held-out: `17/25`
- soft-k beats ZM on held-out: `13/25`
- soft-k beats Zipf on held-out: `24/25`
- four-way winner counts: Zipf `0`, ZM `11`, MOE `4`, soft-k `10`
- median soft-k minus MOE, ZM, and Zipf held-out deltas
- regularization study counts for fixed-`k`, fixed-`k,w`, and soft-`k`
- residual step-2 help counts for free PMF and soft-`k`

## Methods

No PMF optimization is rerun. This migration consolidates saved historical outputs.

Canonical Table 4 uses the v4 verification CSV because it is what the manuscript and verification docs cite. It is **canonical-by-documentation, not canonical-by-reproducibility**: `zipf_v4_verification` has no dedicated reproducible producer script. Fully reproducing Table 4 from scratch would require reconstructing the exact soft-k optimization path and initialization/warm-start behavior.

## Optimization-path sensitivity

Three corpora sit on different soft-k local optima depending on optimization path:

- Critique of Pure Reason
- Jane Eyre
- King James Bible

The λ_k choice is the same (`0.0003`) across sources. The difference is the fitted local optimum, especially `k,w`, under different warm-start behavior. This is a characteristic of the model family, not a provenance bug. Future reruns should expect similar-but-not-identical local optima on sensitivity-prone corpora unless the initialization strategy is made part of the canonical protocol.

The `softkw` bundle's repeated soft-k column inherits from legacy soft-k rather than independently recomputing soft-k, so it does **not** serve as independent confirmation of v4's values.

## Rerun command

```bash
python3 experiments/7a_canonical_pmf_family/scripts/run_7a.py
```

This command regenerates the consolidated CSVs from saved historical bundles. It does not rerun model fitting.

## Verification mapping

- `table4_fourway.csv` must contain `25` rows.
- Winner counts must be Zipf `0`, ZM `11`, MOE `4`, soft-k `10`.
- `softk_source_diagnostic.csv` must flag exactly `3` local-optimum discrepancies against `softk_splitfit`.
- Full-refit outputs must remain under `outputs/fullrefit/` and must not be used for held-out Table 4 claims.

## AUDIT

When source bundles differ on methodological-path grounds, preserve both with explicit documentation rather than declaring one "right" by fiat. This differs from the 3b BIC issue, where one source had false provenance and a clearer correction existed.

When Paper 2 is drafted, revisit the canonical Table 4 source decision. Options:

- rerun the PMF optimization with a documented initialization and warm-start strategy to produce canonical-by-reproducibility outputs;
- or retain the documented-canonical-by-citation approach, but make initialization/path sensitivity explicit in the manuscript.

