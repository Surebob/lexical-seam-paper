# J1 README

## Join id

`J1_paper_scope_rollup`

## Research question

Which manuscript-wide corpus-scope counts combine results from the English, smooth-fit, and multilingual experiments?

## Why this is a joined output

The mapped manuscript claims combine outputs from experiments `1a`, `3a`, and `6`. The LaTeX builder must not compute cross-experiment joins, so J1 emits a single CSV cell for each cross-experiment count.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv`
- `/Volumes/External2TB/emlexperiment/experiments/3a_smooth_two_regime_fits/outputs/aggregate_statistics.csv`
- `/Volumes/External2TB/emlexperiment/experiments/6_multilingual_extension/outputs/aggregate_statistics.csv`
- `/Volumes/External2TB/emlexperiment/experiments/6_multilingual_extension/outputs/table3_multilang.csv`

## Inputs

J1 reads only canonical Phase 2 aggregate outputs. It does not read historical bundles directly and does not rerun any experiment.

## Outputs produced

### `outputs/paper_scope_summary.csv`

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Rows:
- `english_corpus_count`
- `non_english_corpus_count`
- `total_corpus_count`
- `language_family_count`
- `smooth_beats_single_zm_total_count`
- `bregman_indoeuropean_alphabetic_count`
- `indoeuropean_alphabetic_total_count`

### `outputs/manifest.json`

Machine-readable output inventory and upstream provenance.

## Manuscript lines fed

- line `13`: corpus-scope counts and language-family count.
- lines `13`, `31`, `365`, `498`, `518`: `29/29` Indo-European alphabetic Bregman-scope restatements.
- line `534`: corpus-scope restatement.

## Methods

J1 performs integer rollups from upstream aggregate CSVs. The `29/29` Indo-European alphabetic count follows the manuscript's analysis-scope convention: 25 English corpora plus Latin, French, Spanish, and Dutch. Russian is preserved in the multilingual table but excluded from this alphabetic/analysis-scope denominator.

## AUDIT

Experiment `3a` is still represented in the canonical experiment directory by the earlier coupled-logistic migration, while v5.1 manuscript prose has moved toward the decoupled-erf model canon. J1 uses only the count that is stable across this model-canon change: the smooth model beats single-ZM on all 25 English corpora.

## Canonical claim mapping

- `english_corpus_count`, `non_english_corpus_count`, and `language_family_count` map to line-13 and line-534 corpus-scope claims.
- `bregman_indoeuropean_alphabetic_count` and `indoeuropean_alphabetic_total_count` map to `29/29` claims.
- `smooth_beats_single_zm_total_count` maps to joined smooth-model scope claims.

## Rerun command

```bash
python3 joins/J1_paper_scope_rollup/scripts/build_paper_scope_rollup.py \
  --output-dir joins/J1_paper_scope_rollup/outputs
```

## Verification mapping

- `paper_scope_summary.csv` must contain exactly the seven metric rows listed above.
- `english_corpus_count + non_english_corpus_count` must equal `total_corpus_count`.
- `bregman_indoeuropean_alphabetic_count` must equal `indoeuropean_alphabetic_total_count`.
