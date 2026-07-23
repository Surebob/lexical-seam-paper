# Experiment 7c README

## Experiment id

`7c_lambda_metadata_predictability`

## Research question

Can the per-corpus soft-k regularization strength `lambda_k` be predicted from coarse corpus-heterogeneity metadata?

## Current role

Per the 2026-04-19 scoping decision, the PMF arc (`7a` through `9` and related diagnostics) is queued for a separate Paper 2. In the current Phase 2 structure, 7c preserves the metadata-predictability diagnostic used by manuscript v5.1 as a limitation on the parsimony of soft-k.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_angle2_lambda_metadata`](/Volumes/External2TB/emlexperiment/results/zipf_angle2_lambda_metadata)

The source bundle records that best `lambda_k` values were taken from the existing soft-k sweep, and that structural metadata are coarse hand-coded estimates intended to test whether lambda varies with corpus heterogeneity at all.

## Outputs produced

### `outputs/lambda_metadata_per_corpus.csv`

One row per English corpus (`25` rows), including `lambda_k`, `log10_lambda_k`, token/vocabulary counts, hand-coded structure fields, log-transformed metadata fields, and the composite heterogeneity score.

### `outputs/lambda_metadata_summary.csv`

One row per predictor (`7` rows), including Pearson and Spearman correlations with `log10(lambda_k)`.

### `outputs/lambda_by_structure.csv`

Median `lambda_k` by coarse corpus-structure category.

### `outputs/aggregate_statistics.csv`

Claim-facing aggregate rows including the four manuscript-listed Pearson correlations and `max_abs_metadata_correlation`.

### `outputs/manifest.json`

Machine-readable inventory of source bundles and output schemas.

## Manuscript claims fed

- Section 3.11: all manuscript-listed Pearson correlations with `log10(lambda_k)` are weak:
- log unit count: `r = -0.0875`
- log author count: `r = -0.2673`
- log era span: `r = -0.1007`
- heterogeneity score: `r = -0.1579`
- Discussion restatement: all manuscript-listed correlations have `|r| < 0.27`.

## Methods

No metadata coding or correlation computation is rerun beyond reconstructing saved table outputs from `summary.json`. This migration reads the saved correlation dictionary, normalizes it into CSV form, and archives the historical source files.

## Rerun command

```bash
python3 experiments/7c_lambda_metadata_predictability/scripts/run_7c.py
```

This command regenerates the consolidated CSVs from saved historical bundles.

## Verification mapping

- `lambda_metadata_per_corpus.csv` must contain `25` rows.
- `lambda_metadata_summary.csv` must contain `7` predictor rows.
- `max_abs_metadata_correlation` must equal `0.26729532763453223` for the manuscript-listed metadata predictors.
- `correlation_pearson_log_author_count` must be the largest absolute manuscript-listed correlation.

## AUDIT

The metadata fields are coarse hand-coded estimates, not independently measured corpus annotations. This limits the conclusion to: these particular coarse metadata summaries do not predict the selected soft-k `lambda_k`. It does not rule out more detailed predictors.

The source uses historical soft-k sweep lambdas. 7a remains the canonical Table 4 source for Paper 1/Paper 2 PMF table provenance, while this experiment preserves the metadata-predictability diagnostic tied to those lambda choices.
