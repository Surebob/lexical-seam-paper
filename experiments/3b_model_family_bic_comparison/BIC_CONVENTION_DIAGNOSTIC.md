# 3b BIC Convention Diagnostic

This note supplements `CONSOLIDATION_PLAN.md`. It is read-only provenance work; no canonical `3b` outputs have been migrated or regenerated.

## Executive Finding

The divergence between `results/zipf_bic_comparison` and `results/zipf_continuous_piecewise` / `results/zipf_moezipf_comparison` is not caused by a different BIC formula, log base, sample-size convention, or likelihood convention for the shared rank-curve table. It is caused by a different baseline model being inserted into the first column:

- `zipf_bic_comparison/single_zm_bic` is a 3-parameter Zipf-Mandelbrot baseline computed from `single_zm_rmse`.
- `zipf_moezipf_comparison/zipf_rank_bic` is a 1-parameter pure Zipf baseline computed from the rank curve implied by a Zipf MLE fit.
- `zipf_continuous_piecewise/zipf_bic` imports `zipf_moezipf_comparison/zipf_rank_bic` exactly, then combines it with MOEZipf and the piecewise/smooth BIC fields.

The old six-family table therefore mixed a pure-Zipf baseline with ZM-family piecewise and smooth alternatives, while the v4 claim map names that column `single_zm_bic`. That name is not correct for the `zipf_continuous_piecewise` convention.

## Mechanism of the Difference

All three rank-curve conventions use the same BIC functional form:

```text
BIC = p * ln(n) + n * ln(MSE)
MSE = RMSE^2
n = vocabulary size
```

The source-code evidence:

- `zipf_bic_comparison.py:13-28` sets `single_zm` parameter count to `3` and defines `bic_from_rmse(rmse, p, n) = p * log(n) + n * log(rmse^2)`.
- `zipf_bic_comparison.py:39-51` reads `single_zm_rmse` from `results/zipf_correct_model_reranked_all/summary.json` and computes `single_zm_bic`.
- `zipf_moezipf_comparison.py:223-240` fits a pure Zipf MLE, converts that distribution to a rank curve, computes `zipf_rmse`, and then computes `zipf_rank_bic` with `p=1`.
- `zipf_continuous_piecewise.py:249-267` builds its table by importing `ref_moe["zipf_rank_bic"]` as the `zipf` column and importing `piecewise_k500_bic`, `reranked_7param_sqrtv_bic`, and `reranked_8param_bic` from `zipf_bic_comparison`.

The shared fields verify the mechanism:

- `zipf_continuous_piecewise.zipf_bic == zipf_moezipf_comparison.zipf_rank_bic` for all 25 corpora.
- `zipf_continuous_piecewise.piecewise_k500_bic`, `reranked_7param_sqrtv_bic`, and `reranked_8param_bic` match `zipf_bic_comparison` exactly for all 25 corpora.
- Only the first baseline column differs, because it is pure Zipf in one table and single ZM in the other.

Representative rows:

| Corpus | n | pure-Zipf RMSE | single-ZM RMSE | continuous/pure-Zipf BIC | bic_comparison/single-ZM BIC | Difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Aesop's Fables | 4,394 | 0.2362856650 | 0.1673532078 | -12,670.1805 | -15,684.6926 | 3,014.5121 |
| Moby Dick | 16,956 | 0.2835704123 | 0.1615487122 | -42,729.3795 | -61,790.6164 | 19,061.2369 |
| Shakespeare | 24,458 | 0.4060965387 | 0.1844304250 | -44,071.2515 | -82,661.3516 | 38,590.1002 |
| War and Peace | 17,445 | 0.4773826321 | 0.1829823012 | -25,789.1882 | -59,226.6839 | 33,437.4957 |

Across all 25 corpora, `continuous_piecewise.zipf_bic - bic_comparison.single_zm_bic` has min `-27344.8912`, median `15522.0087`, mean `13600.1232`, and max `38590.1002`. The negative minimum occurs because the one-parameter pure-Zipf penalty can sometimes offset fit degradation on a small corpus, but most corpora strongly favor the 3-parameter ZM baseline by rank RMSE.

## Which Convention the Manuscript Cites

### Current manuscript v5.1

`MANUSCRIPT_DRAFT_v5_1.md` no longer cites the historical six-family Table 2 as its Section 3.5 table. Section 3.5 now uses the decoupled gate-family comparison:

- `MANUSCRIPT_DRAFT_v5_1.md:248-260` reports logistic/tanh/erf/algebraic/arctan results.
- The table builder `build_manuscript_v5_latex.py:371-394` reads `results/s2_v3_windows_full_outputs_2026-04-18/s2_v3_per_corpus_results.csv`.
- The cited values `erf wins 24/25`, `arctan wins 1/25`, `median BIC spread 677.37`, and `mean BIC spread 864.66` come from the S2 v3 gate-family outputs, not from either historical 3b convention.

Section 3.7 in v5.1 mentions the older coupled-logistic model-family comparison narratively:

- `MANUSCRIPT_DRAFT_v5_1.md:286` says the coupled-logistic BIC comparison showed MOEZipf can be competitive through parsimony.
- It does not cite specific historical BIC values or BIC spreads.

Thus, under v5.1, the historical 3b six-family table is legacy provenance rather than the current Table 2 source.

### Manuscript v4 / legacy claim map

The v4 LaTeX builder did use the `zipf_continuous_piecewise` convention:

- `build_manuscript_v4_latex.py:371-394` reads `results/zipf_continuous_piecewise/continuous_piecewise_table.csv`.
- It renders `row["zipf_bic"]`, `row["moezipf_bic"]`, `row["piecewise_k500_bic"]`, `row["continuous_piecewise_bic"]`, `row["reranked_7param_sqrtv_bic"]`, and `row["reranked_8param_bic"]`.
- The rendered header labels the first column as `ZM`, although the source field is pure-Zipf rank-BIC, not single-ZM BIC.

`results/MANUSCRIPT_CLAIM_TO_CSV_MAP_v4.md` carries this legacy naming mismatch forward:

- It maps Table 2 to `3b/outputs/table2_model_family.csv`.
- It expects a column named `single_zm_bic`.
- The planned source mapping currently pointed that column at `zipf_continuous_piecewise.zipf_bic`, which is pure Zipf.

The v4 prose claims that MOEZipf wins `13/25` and the smooth family wins `12/25` match the `zipf_continuous_piecewise` / `zipf_moezipf_comparison` rank-BIC convention, not the `zipf_bic_comparison` single-ZM-only table. However, the column label "ZM" / `single_zm_bic` does not match that convention.

## Internal Documentation and Git History

Internal documentation found:

- `results/zipf_bic_comparison/report.md` documents the RMSE-BIC formula and explicitly labels the baseline as "Single ZM".
- `results/zipf_moezipf_comparison/report.md` documents two conventions: likelihood AIC/BIC for Zipf and MOEZipf, and rank-curve RMSE-BIC for cross-model winner comparison.
- `results/zipf_continuous_piecewise/report.md` says it uses the same RMSE-based BIC formula as earlier model-comparison runs and labels the first column "Single Zipf".
- `zipf_continuous_piecewise.py` explicitly imports both `MOE_SUMMARY_PATH` and `BIC_SUMMARY_PATH`, then selects pure-Zipf and MOEZipf rank-BIC values from the MOE bundle while selecting piecewise/smooth BIC values from the BIC bundle.

Git history found:

- The relevant scripts appear only in the initial commit (`150cd42`, 2026-04-16) and were not modified by later commits.
- No commit message or diff hunk explains the baseline substitution as a deliberate correction, a bug fix, or a methodological change.

Interpretation:

- The coexistence of `results/zipf_bic_comparison/report.md` ("Single ZM") and `results/zipf_continuous_piecewise/report.md` ("Single Zipf") indicates the source bundles did distinguish the baselines locally.
- The mismatch appears when the continuous-piecewise table is consumed by the manuscript builder and v4 claim map under the label `ZM` / `single_zm_bic`.
- I found no internal note saying this relabeling was deliberate. It should be treated as a naming/provenance divergence until the user chooses the canonical convention.

## Decision Options

### Option A: Use the continuous_piecewise / moezipf convention as canonical for legacy 3b

Use this if the goal is to reproduce the v4 six-family BIC claims exactly.

- Pros: Matches the old Table 2 builder and the v4 prose winner counts (`MOEZipf 13/25`, smooth family `12/25`).
- Cons: The first column must be renamed to pure Zipf, not single ZM. Calling it `single_zm_bic` would be false provenance.

Required schema correction:

- Replace `single_zm_bic` with `zipf_rank_bic` or `single_zipf_rank_bic`.
- If compatibility with the old claim map requires a `single_zm_bic` column, include it only as an explicit legacy alias with an AUDIT warning, not as a scientifically correct label.

### Option B: Use the bic_comparison convention as canonical

Use this if the scientific question is "how do ZM-family rank-curve models compare against a single-ZM baseline?"

- Pros: The baseline is actually single ZM and aligns with the paper's main ZM residual framing.
- Cons: It does not include MOEZipf in the same table unless MOEZipf is joined in from a different convention. It will not reproduce v4's MOEZipf-vs-smooth winner counts without recomputing/renormalizing the table.

### Option C: Emit both conventions

Use this if preserving both historical claims and corrected provenance is the priority.

- Primary table: `table2_model_family_rank_bic.csv` reproduces the legacy `zipf_continuous_piecewise` / `zipf_moezipf_comparison` convention and labels the first column `zipf_rank_bic`.
- Alternative table: `table2_zm_family_rank_bic.csv` preserves `zipf_bic_comparison/single_zm_bic` and the ZM-family piecewise/smooth comparison.
- README AUDIT explicitly states that v4's old Table 2 label conflated pure Zipf with single ZM.

## Recommendation

For v5.1 packaging, do not treat historical 3b as the current Table 2 source. Current Table 2 is S2 v3 gate-family output.

For legacy 3b migration, I recommend Option C unless the user wants minimum output count. It preserves the historical v4 table exactly while preventing the false `single_zm_bic` label from becoming canonical again.

If only one convention is allowed, choose based on target:

- Reproduce v4 manuscript claims: choose Option A, but rename the first column to `zipf_rank_bic`.
- Compare true single-ZM against smooth/piecewise ZM-family models: choose Option B, and do not claim it reproduces v4 MOEZipf-vs-smooth BIC counts.
