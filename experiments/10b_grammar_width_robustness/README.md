# Experiment 10b README

## Experiment id

`10b_grammar_width_robustness`

## Research question

Does widening the symbolic-regression operator dictionary change the step-2 residual finding, or does it preserve the high-c Bregman result and the low-c manifold structure?

## Current role

This experiment feeds the grammar-width robustness discussion. It is a historical migration of saved widened-grammar outputs; no symbolic search is rerun.

## Upstream dependencies

- [`/Volumes/External2TB/emlexperiment/results/zipf_widened_grammar_diagnostic`](/Volumes/External2TB/emlexperiment/results/zipf_widened_grammar_diagnostic)
- [`/Volumes/External2TB/emlexperiment/results/zipf_widened_grammar_extended`](/Volumes/External2TB/emlexperiment/results/zipf_widened_grammar_extended)

## Outputs produced

### `outputs/widened_diagnostic_per_mode.csv`

Six rows from the original deeper widened-grammar diagnostic: three corpora (`Shakespeare`, `War and Peace`, `Pride and Prejudice`) x two residual modes (`zm_residual`, `post_is_residual`).

### `outputs/widened_diagnostic_step_winners.csv`

Step-level widened-grammar winners for the original diagnostic (`54` rows = `6` diagnostic modes x `9` steps).

### `outputs/widened_diagnostic_top10.csv`

Top-10 widened-grammar candidates for every saved diagnostic step (`540` rows).

### `outputs/widened_lowc_manifold_summary.csv`

Seventeen-row extended low-c/multilingual step-2 summary with widened winner identity, cosines against exp / `x^x - sqrt(x)` / IS, and span-R² in the low-c manifold basis.

### `outputs/aggregate_statistics.csv`

Claim-facing aggregates for high-c survival, low-c erf/sin flip counts, and manifold-membership counts.

### `outputs/manifest.json`

Machine-readable inventory of source bundles, output schemas, and claim rows satisfied.

## Manuscript claims fed

- High-c IS-Bregman winner survives widened grammar on Shakespeare and War and Peace: `2/2`.
- English low-c widened winners flip to erf/sin-family forms on `11/14` corpora.
- Low-c/multilingual manifold membership holds on `17/17` rows.
- The minimum weighted centered span R² is `0.9759519786969072`, so all `17/17` rows exceed the `0.975` threshold cited in manuscript prose.

## Methods

No widened-grammar search is rerun. The original diagnostic added nine unary operators: `besselj`, `cos`, `cosh`, `erf`, `gamma`, `sin`, `sinh`, `tan`, and `tanh`. The deterministic beam width is `50`; the deeper diagnostic ran to step `10` and kept all expressions until step `2`. The extended low-c check materialized step-2 outputs only.

## Rerun command

```bash
python3 experiments/10b_grammar_width_robustness/scripts/run_10b.py
```

This command regenerates the consolidated CSVs from saved historical bundles. It does not rerun symbolic search.

## Verification mapping

- `widened_diagnostic_per_mode.csv` must contain `6` rows.
- `widened_diagnostic_step_winners.csv` must contain `54` rows.
- `widened_diagnostic_top10.csv` must contain `540` rows.
- `widened_lowc_manifold_summary.csv` must contain `17` rows.
- `high_c_zm_step2_survives_count` must be `2`.
- `lowc_english_erf_sin_flip_count` must be `11`.
- `lowc_all_manifold_yes_count` must be `17`.

## AUDIT

The high-c widened-grammar survival claim is supported only for Shakespeare and War and Peace. The low-c manifold extension covers 14 English low-c corpora plus three multilingual c≈0 corpora, selected according to the historical low-c family rather than a literal c-threshold.
