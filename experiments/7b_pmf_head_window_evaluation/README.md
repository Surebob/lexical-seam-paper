# Experiment 7b README

## Experiment id

`7b_pmf_head_window_evaluation`

## Research question

Does restricting held-out PMF evaluation to progressively larger head windows change which discrete family wins, and does head-window scoring suppress or preserve the symbolic seam in soft-k residuals?

## Current role

Per the 2026-04-19 scoping decision, the PMF arc (`7a` through `9` and related diagnostics) is queued for a separate Paper 2. In the current Phase 2 structure, 7b's role is research-record preservation and provenance cleanup for the head-window diagnostic used by manuscript v5.1.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_angle1_head_windows_splitfit`](/Volumes/External2TB/emlexperiment/results/zipf_angle1_head_windows_splitfit)
- This source bundle itself records the upstream PMF sources:
- `results/zipf_seam_mandelbrot_softk_splitfit/summary.json`
- `results/zipf_seam_mandelbrot_pmf/summary.json`

## Outputs produced

### `outputs/head_window_per_corpus.csv`

One row per corpus per cutoff (`top50`, `top100`, `top200`, `top500`, `top1000`, `full`), for `150` rows total. Columns include the four held-out average NLL values (`zipf`, `zm`, `moe`, `softk`), the held-out winner, soft-k deltas versus MOE/ZM/Zipf, and the soft-k residual step-2 winner/gain/help flag under that cutoff.

### `outputs/aggregate_statistics.csv`

Claim-facing aggregate rows. The manuscript currently cites:

- `top100_winner_count_zm`
- `top100_softk_beats_moe_count`
- `top200_winner_count_zipf`, `top200_winner_count_zm`, `top200_winner_count_moe`, `top200_winner_count_softk`
- `top1000_winner_count_zipf`, `top1000_winner_count_zm`, `top1000_winner_count_moe`, `top1000_winner_count_softk`
- `top50_softk_step2_help_count`, `top100_softk_step2_help_count`, `top200_softk_step2_help_count`, `top500_softk_step2_help_count`, `top1000_softk_step2_help_count`, `full_softk_step2_help_count`
- `top50_modal_winner`, `top100_modal_winner`, `top200_modal_winner`, `top500_modal_winner`, `top1000_modal_winner`
- `max_head_window_rank`

### `outputs/manifest.json`

Machine-readable inventory of source bundles and output schemas.

## Manuscript claims fed

- Section 3.11 head-window evaluation: top-100 ZM wins `25/25`; soft-k beats MOE on `11/25`.
- Section 3.11 head-window evaluation: top-200 winner counts Zipf `0`, ZM `15`, MOE `8`, soft-k `2`.
- Section 3.11 head-window evaluation: top-1000 winner counts Zipf `0`, ZM `1`, MOE `14`, soft-k `10`.
- Section 3.11 head-window symbolic seam diagnostic: step-2 help counts `{top-50: 2/25, top-100: 6/25, top-200: 0/25, top-500: 1/25, top-1000: 4/25, full: 4/25}`.
- Discussion restatement: ZM is modal on top-50/top-100/top-200; MOE is modal on top-500/top-1000.
- Protocol restatement: largest finite head-window cutoff is `1000`.

## Methods

No model fitting or fresh head-window evaluation is rerun. This migration reads the saved `zipf_angle1_head_windows_splitfit` summary and table, normalizes column names, recomputes only simple aggregate counts, and archives the historical source files.

The source bundle's metadata states that soft-k held-out windows use split-fit parameters only, while full-refit parameters are used only for full-corpus residual step-2 diagnostics.

## Rerun command

```bash
python3 experiments/7b_pmf_head_window_evaluation/scripts/run_7b.py
```

This command regenerates the consolidated CSVs from saved historical bundles. It does not rerun PMF fitting or held-out scoring.

## Verification mapping

- `head_window_per_corpus.csv` must contain `150` rows.
- `top100_winner_count_zm` must be `25`.
- `top100_softk_beats_moe_count` must be `11`.
- Top-200 winner counts must be Zipf `0`, ZM `15`, MOE `8`, soft-k `2`.
- Top-1000 winner counts must be Zipf `0`, ZM `1`, MOE `14`, soft-k `10`.
- Step-2 help counts must be top-50 `2`, top-100 `6`, top-200 `0`, top-500 `1`, top-1000 `4`, full `4`.

## AUDIT

This is a single-source migration, so no multi-bundle conflict resolution was needed. The PMF arc is queued for Paper 2; if Paper 2 uses this diagnostic, it should decide whether to preserve the saved split-fit head-window source or rerun the full PMF protocol with a documented optimizer initialization strategy shared with 7a.
