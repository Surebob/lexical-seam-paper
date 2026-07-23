# Experiment 10a Consolidation Plan

## Experiment id

`10a_search_depth_robustness`

## Research question

Does deeper symbolic search reveal a new interpretable residual structure, or does it produce syntactically complex formulas that approximate the same smooth/polynomial residual curve more verbosely?

## Source bundles inspected

### `results/zipf_step10_ablation`

- `summary.json`: 2 rows, Shakespeare and War and Peace.
- `report.md`: prose summary.
- `shakespeare_step10_ablation.svg`, `war_and_peace_step10_ablation.svg`: diagnostic plots.
- Stored content:
- step-2 expression for both corpora: `sub[sub[x,1],log[x]]`.
- step-10 dominant expression trees.
- grid comparison of `(step10 - step2)` on 1000 points in normalized `x in [0.05, 0.95]`.
- polynomial fits of degrees 3, 4, and 5 to `(step10 - step2)`.
- degree-5 R² values: Shakespeare `0.9772794594503794`, War and Peace `0.9927384752120223`.

### `results/zipf_head_poly_decomposition`

- `summary.json`: 2 rows, Shakespeare and War and Peace.
- `report.md`: prose summary.
- Stored content:
- ZM RMSE.
- step-2 RMSE.
- step-10 monster RMSE.
- direct residual decomposition fits of `ZM + [(x - 1) - log(x)] + poly_d(x)` for degrees 3, 4, and 5.
- degree-5 direct polynomial RMSE: Shakespeare `0.1505496191235793`, War and Peace `0.14581271318709715`.
- step-10 monster RMSE: Shakespeare `0.1598718529967204`, War and Peace `0.14534740465186036`.

### `results/zipf_head_poly_transfer`

- `summary.json`: `records` for Shakespeare and War and Peace, plus `2` transfer rows.
- `report.md`: prose summary.
- Stored content:
- in-domain degree-5 polynomial coefficients for Shakespeare and War and Peace.
- zero-shot transfer Shakespeare -> War and Peace: RMSE `0.15127747094311222`.
- zero-shot transfer War and Peace -> Shakespeare: RMSE `0.15565236538090854`.

## Cross-source diagnostics

The overlapping RMSE values in `zipf_head_poly_decomposition` and `zipf_head_poly_transfer` match exactly for both corpora:

- `zm_rmse`
- `step2_rmse`
- `step10_rmse`
- in-domain `poly5_rmse`
- degree-5 coefficients

`zipf_step10_ablation` is not a conflicting source for these RMSEs. It fits polynomials to the **incremental function** `(step10 - step2)` on a fixed grid. `zipf_head_poly_decomposition` fits `ZM + fixed step-2 + poly_d` directly to full-corpus residuals. These are related diagnostics for the same question but not duplicate computations of the same metric.

No false-provenance issue was found. The source divergence mechanism is methodological complementarity:

- `step10_ablation`: what shape is the incremental step-10 correction?
- `head_poly_decomposition`: can a simple polynomial reproduce the step-10 improvement in RMSE?
- `head_poly_transfer`: do fitted degree-5 coefficients transfer between the two exemplar high-c corpora?

## Manuscript/support mismatch

Current `MANUSCRIPT_DRAFT_v5_1.md` contains two transfer statements that are not supported by the saved `zipf_head_poly_transfer` bundle:

- Section 2.20 says coefficients were evaluated for `Shakespeare -> War and Peace`, `War and Peace -> Shakespeare and King James Bible`, and `Moby Dick -> Bible`.
- Section 4 limitations says `War and Peace -> Bible RMSE 0.150 vs in-domain 0.149`.

The saved transfer bundle contains only:

- Shakespeare -> War and Peace.
- War and Peace -> Shakespeare.

There are no saved King James Bible or Moby Dick transfer rows in the inspected source bundles. Under the migration rules, these missing transfer claims should be flagged in `BLOCKED.md` / README `AUDIT` after execution rather than reconstructed by running fresh computation.

## Canonical-source decision proposed

Use all three source bundles as separate first-class components of one consolidated experiment:

1. `zipf_step10_ablation` is canonical for incremental step-10-minus-step-2 polynomial R² values and saved step-10 expression trees.
2. `zipf_head_poly_decomposition` is canonical for direct polynomial-decomposition RMSE values.
3. `zipf_head_poly_transfer` is canonical for the two saved zero-shot transfer rows only.

Do not infer or fabricate missing Bible/Moby transfer rows. If those rows are required for manuscript support, this experiment should be marked partially blocked until a fresh transfer run is explicitly authorized.

## Planned canonical outputs

### `outputs/step10_ablation_per_corpus.csv`

Source: `zipf_step10_ablation/summary.json`.

Rows: 2.

Columns:

- `slug`
- `corpus`
- `step2_expression`
- `step10_expression`
- `grid_min`
- `grid_max`
- `difference_min`
- `difference_max`
- `difference_mean_abs`
- `difference_std`
- `best_poly_degree`
- `best_poly_r2`
- `polyfit_degree3_r2`
- `polyfit_degree4_r2`
- `polyfit_degree5_r2`
- `polyfit_degree3_coefficients_json`
- `polyfit_degree4_coefficients_json`
- `polyfit_degree5_coefficients_json`
- `plot_path`
- `canonical_source`

Feeds Section 3.12 and conclusion/verification claims about degree-5 R² `0.977` and `0.993`.

### `outputs/polynomial_decomposition_per_corpus.csv`

Source: `zipf_head_poly_decomposition/summary.json`.

Rows: 2.

Columns:

- `slug`
- `corpus`
- `zm_rmse`
- `step2_rmse`
- `step10_rmse`
- degree-3/4/5 polynomial RMSE and coefficients
- `poly5_minus_step10_rmse`
- `canonical_source`

Feeds Section 3.12 and conclusion claims about direct polynomial decomposition RMSE `0.1506` vs `0.1599` on Shakespeare and `0.1458` vs `0.1453` on War and Peace.

### `outputs/poly_transfer.csv`

Source: `zipf_head_poly_transfer/summary.json`.

Rows: 2.

Columns:

- `source_corpus`
- `target_corpus`
- `source_poly5_coefficients_json`
- `target_zm_rmse`
- `target_step2_rmse`
- `target_in_domain_poly5_rmse`
- `target_step10_rmse`
- `zero_shot_rmse`
- `zero_shot_minus_in_domain_poly5`
- `zero_shot_minus_step10`
- `canonical_source`

Feeds only the two saved transfer claims. Missing Bible/Moby transfer claims must be flagged.

### `outputs/aggregate_statistics.csv`

Claim-facing key/value rows:

- `shakespeare_step10_minus_step2_poly5_r2`
- `war_and_peace_step10_minus_step2_poly5_r2`
- `shakespeare_poly5_direct_rmse`
- `shakespeare_step10_rmse`
- `war_and_peace_poly5_direct_rmse`
- `war_and_peace_step10_rmse`
- `shakespeare_to_war_and_peace_zero_shot_rmse`
- `war_and_peace_to_shakespeare_zero_shot_rmse`
- `saved_transfer_count`
- `missing_transfer_claims_count`

### `outputs/manifest.json`

Machine-readable inventory of source bundles and output schemas.

### `BLOCKED.md`

Emit if migration is approved, documenting that the Bible/Moby transfer rows cited in v5.1 are not reconstructible from current saved bundles without a fresh run.

## Precedence rules

1. Use `zipf_step10_ablation` for step-10 expression trees and incremental polynomial R² only.
2. Use `zipf_head_poly_decomposition` for direct polynomial-decomposition RMSE only.
3. Use `zipf_head_poly_transfer` for transfer rows only.
4. Where `zipf_head_poly_decomposition` and `zipf_head_poly_transfer` duplicate in-domain values, require exact agreement; if execution finds drift above `1e-12`, stop and flag.
5. Do not use manuscript prose to fill missing transfer rows.

## Union vs select rules

- Union the three diagnostics into separate canonical outputs.
- Select one source per metric according to the rules above.
- Preserve all source files in `archive/provenance/`.
- Preserve the SVG plots from `zipf_step10_ablation` as source artifacts.

## README AUDIT entries to include after execution

- The formal polynomial approximability check covers Shakespeare and War and Peace only.
- The saved transfer bundle contains only Shakespeare -> War and Peace and War and Peace -> Shakespeare.
- Manuscript v5.1 statements about War and Peace -> Bible and Moby Dick -> Bible transfer are unsupported by current saved outputs and require either text revision or an authorized fresh run.
- This experiment supports the interpretation that depth-10 search produces syntactic complexity without new interpretable mathematical structure, but it does not establish the polynomial decomposition across all 25 corpora.
