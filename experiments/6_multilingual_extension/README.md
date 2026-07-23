# Experiment 6 README

## Experiment id

`6_multilingual_extension`

## Research question

Does the residual/Bregman pattern and smooth-model improvement extend to seven non-English corpora?

## Why this is one experiment

The historical multilingual runs all apply the same operation to additional non-English corpora. `zipf_multilang_romance` contains the canonical 7-corpus table; `zipf_multilang_verify` supplies fitted ZM parameters and formula checks for the same rows.

## Upstream dependencies

- `/Volumes/External2TB/emlexperiment/results/zipf_multilang_romance/`
- `/Volumes/External2TB/emlexperiment/results/zipf_multilang_verify/`

## Inputs

- `results/zipf_multilang_romance/summary.json`
- `results/zipf_multilang_romance/multilang_table.csv`
- `results/zipf_multilang_verify/summary.json`

## Outputs produced

### `outputs/multilang_per_corpus.csv`

One row per non-English corpus, merging the romance and verification bundles.

### `outputs/table3_multilang.csv`

Table-body data for manuscript Table 3.

Columns:
- `language`
- `corpus`
- `vocabulary_size`
- `zm_c`
- `single_zm_rmse`
- `smooth_rmse`
- `step2_winner_display`
- `k_over_sqrt_v`

### `outputs/aggregate_statistics.csv`

Rows:
- `non_english_corpus_count`
- `smooth_beats_single_zm_count`
- `transition_fraction_min`
- `transition_fraction_max`
- `transition_fraction_near_half_count`
- `multilang_bregman_winner_count`
- `latin_french_spanish_dutch_c_lt_5_count`
- `very_low_c_multilang_count`

### `outputs/manifest.json`

Machine-readable inventory and source-bundle provenance.

## Manuscript lines fed

- line `153`: Russian/Mandarin/Arabic token counts.
- lines `345-355`: Table 3 body.
- line `359`: `7/7`, transition-fraction interval, and near-half count.
- line `361`: four of seven Bregman-style winners and low-`c` Romance/Germanic rows.
- line `363`: very-low-`c` multilingual rows and formula metadata.
- line `498`: non-English `7/7` restatement.

## Methods

This is a read-only migration of historical multilingual outputs. It joins rows by slug between `zipf_multilang_romance/summary.json` and `zipf_multilang_verify/summary.json`, then emits canonical table and aggregate CSVs.

The smooth-model values here are the historical multilingual `smooth_8param` values. The decoupled-erf model-canon adjustment in this migration applies to `3a`, `3c`, `4`, and `5c`, not to this already-saved multilingual bundle.

## AUDIT

The current manuscript v5.1 incorporates T1's later result that `x^x - sqrt(x)` can be recentered into Bregman form. That theorem-track result is outside this Phase 2 migration and is not mixed into the canonical multilingual CSVs here. This experiment preserves the saved Phase 2 multilingual outputs.

## Canonical claim mapping

- Table 3 body maps to `outputs/table3_multilang.csv`.
- Per-corpus token, vocabulary, RMSE, winner, and transition rows map to `outputs/multilang_per_corpus.csv`.
- Summary counts and interval bounds map to `outputs/aggregate_statistics.csv`.

## Rerun command

```bash
python3 experiments/6_multilingual_extension/scripts/run_6.py \
  --output-dir experiments/6_multilingual_extension/outputs
```

## Verification mapping

- `multilang_per_corpus.csv` must have `7` rows.
- `table3_multilang.csv` must have `7` rows.
- `smooth_beats_single_zm_count` must be `7`.
- `transition_fraction_near_half_count` must be `6`.
- `very_low_c_multilang_count` must be `3`.
