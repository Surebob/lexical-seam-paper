# Experiment 1a README

## Experiment id

`1a_per_corpus_enriched_search`

## Research question

For each of the 25 English corpora, what does deterministic step-2 symbolic search recover on the residual of a single-ZM fit, and how do the resulting winner identities, RMSE improvements, and step-2 beams vary corpus by corpus?

## Why this is one experiment

This is the same symbolic-search operation applied to more data points. Historically the 25 English corpora were saved as separate per-corpus bundles, but under the consolidation rule they belong to one experiment with one script and one tabular output containing one row per corpus.

## Upstream dependencies

- Shared English corpus registry in [`/Volumes/External2TB/emlexperiment/zipf_analysis_common.py`](/Volumes/External2TB/emlexperiment/zipf_analysis_common.py)
- Deterministic enriched-search implementation in [`/Volumes/External2TB/emlexperiment/eml_zipf_enriched_search.py`](/Volumes/External2TB/emlexperiment/eml_zipf_enriched_search.py)
- Single-ZM fitting and corpus preprocessing in [`/Volumes/External2TB/emlexperiment/eml_zipf_experiment.py`](/Volumes/External2TB/emlexperiment/eml_zipf_experiment.py)

## Inputs

- Tokenized rank-frequency data for the 25 English corpora
- Normalized rank coordinate `x = 0.05 + 0.95 * (log r - min log r) / (max log r - min log r)`
- Canonical search config:
  - beam width `50`
  - historical full search depth `10`
  - replay depth for this experiment `2` (exactly sufficient for step-2 outputs)
  - `keep_all_until_step_2 = true`
  - seeds `{x, 1}`
  - no continuous coefficient fitting
  - structural-diversity weight `0.35`
  - `EXP_CLAMP = 30`
  - `VALUE_ABS_LIMIT = 1e6`

## Outputs produced

### `outputs/table1_per_corpus.csv`
One row per English corpus. Canonical source for the rendered body of manuscript Table 1.

Columns:
- `corpus`
- `tokens`
- `vocabulary_size`
- `zm_a`
- `zm_b`
- `zm_c`
- `zm_rmse`
- `step2_winner_expression`
- `step2_winner_display`
- `step2_rmse`
- `rmse_improvement_absolute`
- `rmse_improvement_pct`
- `runner_up_expression`
- `runner_up_rmse`
- `winner_minus_runner_rmse`
- `winner_family`
- `contains_is_in_step2_beam`
- `contains_exp_in_step2_beam`

### `outputs/table1_step2_beam_top10.csv`
One row per `(corpus, beam_rank)` pair for the canonical step-2 top-10 beam.

Columns:
- `corpus`
- `beam_rank`
- `expression`
- `expression_display`
- `rmse`
- `winner_family`
- `uses_new_operator`

### `outputs/aggregate_statistics.csv`
Single-cell aggregates used by manuscript prose.

Columns:
- `metric_name`
- `value`
- `display_format`
- `notes`

Required rows:
- `english_corpus_count`
- `winner_family_is_count`
- `winner_family_exp_count`
- `winner_family_other_count`
- `winner_runner_gap_mean`
- `winner_runner_gap_median`
- `winner_family_gap_below_0p01_count`
- `high_c_is_winner_count`
- `low_c_exp_winner_count`
- `c_transition_band_summary`
- `step2_improvement_pct_min`
- `step2_improvement_pct_max`
- `high_c_is_lower_band`
- `low_c_exp_upper_band`
- `both_bregman_forms_present_in_step2_beam_count`
- `mean_step2_winner_runner_gap`
- `step2_all_improve_count`
- `step2_improvement_pct_mean`
- `step2_improvement_pct_median`
- `step2_positive_improve_count`
- `step2_positive_improvement_pct_min`
- `step2_positive_improvement_pct_max`
- `step2_positive_improvement_pct_mean`
- `step2_positive_improvement_pct_median`
- `step2_winner_runner_gap_mean`
- `terminal_all_improve_count`
- `terminal_improvement_pct_min`
- `terminal_improvement_pct_max`
- `terminal_improvement_pct_mean`
- `terminal_improvement_pct_median`
- `terminal_winner_runner_gap_mean`

### `outputs/historical_per_corpus_diff.csv`
Per-corpus comparison against the 25 historical enriched-search bundles subsumed by this experiment.

Columns:
- `corpus`
- `field`
- `new_value`
- `historical_value`
- `abs_diff`
- `matches`

### `outputs/historical_top10_diff.csv`
Per-rank comparison of the canonical step-2 top-10 beam against each historical bundle.

Columns:
- `corpus`
- `beam_rank`
- `field`
- `new_value`
- `historical_value`
- `abs_diff`
- `matches`

### `outputs/historical_diff_summary.csv`
Aggregate reproducibility summary over the subsumed bundles.

Columns:
- `metric_name`
- `value`
- `notes`

## Manuscript lines fed

Primary:
- line `13` (English-corpus side of the abstract scope statement, via `J1`)
- lines `197–206` (Section 3.1)
- line `202` / Table 1 body
- line `365` / `518` corpus-family totals indirectly via `J1`

Secondary:
- any later restatement of English corpus count or step-2 winner-family totals

## Methods

This experiment reruns the canonical deterministic enriched search across the full 25-corpus English sample using the same search machinery and hyperparameters as the historical per-corpus bundles. The emitted scientific object is the step-2 beam. Because candidate generation through step 2 is independent of later search steps, the consolidated replay stops at step 2 rather than replaying historical steps 3-10; this is a runtime optimization, not a different scientific protocol, and the historical diff outputs explicitly verify that the emitted step-2 rows match the historical 10-step bundles. No statistical inference is performed in this experiment.

## Canonical claim mapping

Examples:
- `line 202 / Table 1 body` -> `outputs/table1_per_corpus.csv`
- `line 202 mean RMSE gap 0.003` -> `outputs/aggregate_statistics.csv`, row `mean_step2_winner_runner_gap`
- `line 13 English corpus count 25` -> `outputs/aggregate_statistics.csv`, row `english_corpus_count`, then joined by `J1`
- `line 197 terminal 25/25 claim` -> `outputs/aggregate_statistics.csv`, row `terminal_all_improve_count`
- `line 197 replacement range candidates` -> `outputs/aggregate_statistics.csv`, rows `terminal_improvement_pct_*` and `step2_positive_improvement_pct_*`

## Rerun command

```bash
python3 experiments/1a_per_corpus_enriched_search/scripts/run_1a.py \
  --output-dir experiments/1a_per_corpus_enriched_search/outputs
```

## Verification mapping

Post-run verification checks:
- row count in `table1_per_corpus.csv` must be `25`
- row count in `table1_step2_beam_top10.csv` must be `250`
- every corpus in the canonical English registry must appear exactly once in `table1_per_corpus.csv`
- `winner_family_is_count + winner_family_exp_count + winner_family_other_count = 25`
- `historical_diff_summary.csv` must report zero expression mismatches and max numeric drift `<= 1e-4` for the subsumed historical bundles
- manuscript Table 1 rendering must be a direct read of `table1_per_corpus.csv` with no builder-time recomputation

## Upstream / downstream provenance

Upstream:
- raw corpora preprocessing
- single-ZM fitting module
- deterministic enriched-search module

Downstream:
- `J1_paper_scope_rollup`
- manuscript Table 1
- Section 3.1 prose checks

## Notes on non-goals

This experiment does **not** include:
- the city-population control (`1b`)
- robustness ablations over `x_low`, guards, or WLS (`1c`)
- winner-vs-Euclidean gap analysis (`1d`)

Those are separate experiments because they answer distinct research questions even though they use related symbolic-search machinery.
