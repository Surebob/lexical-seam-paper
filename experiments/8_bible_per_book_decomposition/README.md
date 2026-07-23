# Experiment 8 README

## Experiment id

`8_bible_per_book_decomposition`

## Research question

Does decomposing the King James Bible into its 66 canonical books resolve the aggregate-level Seam-Mandelbrot PMF failure signature?

## Current role

This experiment feeds the manuscript's anthology-scope claim and Table 5. It remains relevant to Paper 1 as a scope validation for the discrete PMF family, even though the broader PMF arc is queued for Paper 2.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_angle6_bible_books`](/Volumes/External2TB/emlexperiment/results/zipf_angle6_bible_books)
- [`/Volumes/External2TB/emlexperiment/experiments/7a_canonical_pmf_family/outputs/splitfit/pmf_variant_per_corpus.csv`](/Volumes/External2TB/emlexperiment/experiments/7a_canonical_pmf_family/outputs/splitfit/pmf_variant_per_corpus.csv)

The 7a dependency supplies the whole-Bible soft-k residual step-2 gain and the median soft-k step-2 gain across the other 24 English corpora, because those context values are not stored in the Bible per-book bundle itself.

## Outputs produced

### `outputs/bible_per_book.csv`

One row per canonical Bible book (`66` rows), with token/vocabulary counts, held-out average NLL values, soft-k deltas versus ZM/MOE, and residual step-2 gain/help/expression.

### `outputs/table5_bible_summary.csv`

Two-column table-body source for manuscript Table 5. Each row has `metric_label`, rendered value, raw numeric value, and display format.

### `outputs/aggregate_statistics.csv`

Claim-facing aggregate rows reused by manuscript prose, Table 5, and the appendix verification checklist.

### `outputs/manifest.json`

Machine-readable inventory of source bundles and output schemas.

## Manuscript claims fed

- Bible books analyzed: `66`.
- Per-book step-2 helps: `6/66`.
- Per-book soft-k beats ZM: `33/66`.
- Per-book soft-k beats MOE: `45/66`.
- Median per-book soft-k minus ZM: `-0.00013339837693404633`.
- Median per-book soft-k minus MOE: `-0.0028037828229647843`.
- Median per-book step-2 gain: `-0.0036131801096206256`.
- Whole-Bible single-fit soft-k held-out NLL: `6.016015736044842`.
- Aggregate per-book soft-k held-out NLL: `5.604219509298761`.
- Per-book step-2 non-help count: `60/66`.

## Methods

No per-book model fitting is rerun. This migration reads the saved `zipf_angle6_bible_books` summary and table, builds the manuscript Table 5 source CSV, and archives the historical source files. The aggregate per-book held-out NLL is the test-token-weighted average of the per-book soft-k held-out likelihoods as stored in the historical bundle.

## Rerun command

```bash
python3 experiments/8_bible_per_book_decomposition/scripts/run_8.py
```

This command regenerates the consolidated CSVs from saved historical bundles and the already-migrated 7a canonical PMF variant table.

## Verification mapping

- `bible_per_book.csv` must contain `66` rows.
- `per_book_step2_help_count` must be `6`.
- `per_book_softk_beats_zm_count` must be `33`.
- `per_book_softk_beats_moe_count` must be `45`.
- `aggregate_per_book_softk_nll` must be `5.604219509298761`.
- `whole_bible_singlefit_softk_nll` must be `6.016015736044842`.

## AUDIT

This is a single-source per-book migration with one explicit cross-experiment context dependency on 7a for whole-Bible residual-gain comparison. The historical per-book output itself supports the Table 5 counts and aggregate held-out NLL values directly.
