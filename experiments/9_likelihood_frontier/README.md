# Experiment 9 README

## Experiment id

`9_likelihood_frontier`

## Research question

How far can likelihood be pushed beyond the mechanism-first soft-k Seam-Mandelbrot PMF, and what mechanism cost appears when the architecture is relaxed?

## Current role

Per the 2026-04-19 scoping decision, the PMF arc (`7a` through `9` and related diagnostics) is queued for a separate Paper 2. Paper 1 will remove the Pareto-frontier discussion. In the current Phase 2 structure, 9's role is research-record preservation and provenance cleanup, not Paper 1 support.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_hybrid_headtail_splitfit`](/Volumes/External2TB/emlexperiment/results/zipf_hybrid_headtail_splitfit)
- [`/Volumes/External2TB/emlexperiment/results/zipf_hybrid_mechpenalty`](/Volumes/External2TB/emlexperiment/results/zipf_hybrid_mechpenalty)
- [`/Volumes/External2TB/emlexperiment/results/zipf_hybrid_vs_softk_analysis`](/Volumes/External2TB/emlexperiment/results/zipf_hybrid_vs_softk_analysis)
- [`/Volumes/External2TB/emlexperiment/results/zipf_fully_ours_variants`](/Volumes/External2TB/emlexperiment/results/zipf_fully_ours_variants)
- [`/Volumes/External2TB/emlexperiment/experiments/9_likelihood_frontier/softk_comparator_source_diagnostic.csv`](/Volumes/External2TB/emlexperiment/experiments/9_likelihood_frontier/softk_comparator_source_diagnostic.csv)

## Outputs produced

### `outputs/hybrid_headtail_per_corpus.csv`

Canonical manuscript-facing hybrid head-tail source. One row per English corpus (`25` rows) with splitfit held-out Zipf, ZM, MOE, soft-k, and hybrid average NLL values; hybrid deltas; selected `lambda`, `k`, and `w`; and residual step-2 help metrics.

### `outputs/fully_ours_variants.csv`

One row per corpus and architecture (`50` rows = `25` corpora x `nested` and `three_regime`), preserving held-out NLL, BIC, RMSE, residual step-2 gain/help, and fitted parameters as JSON.

### `outputs/mechanism_penalty_sweep.csv`

One row per corpus and penalty strength (`125` rows = `25` corpora x `5` lambda values), preserving the mechanism-penalty hybrid sweep.

### `outputs/hybrid_vs_softk_diagnostic.csv`

One row per corpus documenting the comparator divergence between the current splitfit soft-k comparator and the legacy pre-splitfit comparator used by `zipf_hybrid_vs_softk_analysis`.

### `outputs/hybrid_structure_summary.csv`

Long-form summary table for the likelihood-frontier components: canonical hybrid, nested seam, three-regime, mechanism-penalty sweep, legacy comparator diagnostic, and nested protocol constants.

### `outputs/aggregate_statistics.csv`

Claim-facing aggregate rows for hybrid, nested, three-regime, and mechanism-penalty metrics.

### `outputs/manifest.json`

Machine-readable inventory of sources, canonical decisions, output schemas, row counts, and claim-map rows satisfied.

## Manuscript claims fed

This experiment preserves the Paper 2 PMF/likelihood-frontier research record:

- canonical splitfit hybrid beats MOE on `21/25` corpora;
- canonical splitfit hybrid beats soft-k on `22/25` corpora;
- median canonical hybrid minus soft-k held-out NLL is `-0.001677944396615949`;
- hybrid residual step-2 search helps on `9/25` corpora;
- nested seam beats MOE on `21/25`, soft-k on `18/25`, and ZM on `16/25`;
- nested seam is the four-way winner on `14/25` corpora and has residual step-2 help on `11/25`;
- three-regime beats soft-k on `12/25`, with median three-regime minus soft-k `0.00010526860518478287`;
- mechanism penalty preserves held-out step-2 help on `25/25` corpora at every penalty strength and does not satisfy the intended Pareto decision rule.

## Methods

No model fitting is rerun. This migration reads the saved historical outputs, selects canonical manuscript-facing rows according to the approved consolidation plan, and preserves alternate branches as diagnostics.

`zipf_hybrid_headtail_splitfit` is canonical for hybrid-vs-soft-k counts because it uses the current splitfit soft-k comparator. `zipf_hybrid_vs_softk_analysis` is preserved as a legitimate legacy diagnostic because it uses the pre-splitfit derivative-branch soft-k comparator.

The `22/25` canonical hybrid-win count and the `16/25` legacy hybrid-win count reflect different soft-k comparator generations: splitfit/current versus pre-splitfit/legacy. This is not a methodological dispute. Both are valid within their own historical time frame of canonical soft-k.

The mechanism-penalty result is a primary negative finding: held-out residual step-2 help remains `25/25` for all penalty strengths, and no penalty satisfies the intended decision rule of keeping residual help low while retaining strong MOE improvement. This result is worth preserving for Paper 2's eventual likelihood-mechanism Pareto-frontier discussion.

## Rerun command

```bash
python3 experiments/9_likelihood_frontier/scripts/run_9.py
```

This command regenerates the consolidated CSVs from saved historical bundles. It does not rerun model fitting.

## Verification mapping

- `hybrid_headtail_per_corpus.csv` must contain `25` rows.
- `fully_ours_variants.csv` must contain `50` rows.
- `mechanism_penalty_sweep.csv` must contain `125` rows.
- `hybrid_vs_softk_diagnostic.csv` must contain `25` rows and label the analysis comparator as `legacy_pre_splitfit_derivative_branch`.
- `hybrid_beats_softk_count` must be `22` in `aggregate_statistics.csv`.
- The legacy diagnostic must preserve the `16/25` hybrid-vs-soft-k result as a sidecar rather than replacing the canonical count.
- `mechanism_penalty_step2_help_count_all_lambdas` must report `25/25` help for every penalty strength.

## AUDIT

When Paper 2 is drafted, if the comparator convention is revisited, the `16/25` value should be re-examined. The legacy comparator used a different soft-k variant; whether that is the right choice for Paper 2's narrative depends on Paper 2's chosen canonical.

The PMF arc including this experiment is not currently Paper 1 support. It is preserved as a first-class research record so Paper 2 can make an explicit canonical-source decision rather than inheriting an undocumented historical branch.
