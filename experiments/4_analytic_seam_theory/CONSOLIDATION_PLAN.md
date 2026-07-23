# Experiment 4 Consolidation Plan: Analytic Seam Theory

## Research Question

What sign pattern does the smooth seam predict analytically, and how far can tangent-space projection and second-order correction carry that prediction toward the empirical residual?

## Source Bundles

This experiment has three historical source bundles, corresponding to three stages of the same derivation:

| Source bundle | Role | Key files inspected | Canonical-use status |
|---|---|---|---|
| `results/zipf_seam_sign_theory` | Zeroth/local sign-law derivation for the logistic smooth seam | `summary.json`, `report.md` | Legacy coupled-logistic source |
| `results/zipf_seam_projection_theory` | First-order tangent-space projection of the seam against the single-ZM tangent space | `summary.json`, `report.md` | Legacy coupled-logistic source |
| `results/zipf_seam_second_order_theory` | Second-order correction and quadratic surrogate checks | `summary.json`, `report.md` | Legacy coupled-logistic source |

The three bundles are not duplicate runs. They are progressive layers of the same analytic question: raw seam sign law, first-order projection, and second-order correction.

## Model-Canon Conflict

The available source bundles are tied to the historical coupled-logistic smooth model. The reports explicitly derive the seam signs from a logistic gate; for example, `zipf_seam_sign_theory/report.md` states that the smooth seam uses a logistic tail activation and derives the Taylor sign law from that gate.

The current Phase 2 migration instruction says that experiments `3a`, `3c`, `4`, and `5c` should use the decoupled-erf model as canonical where applicable, preserving coupled-logistic results only as legacy archive material.

I do not currently find a decoupled-erf seam-theory source bundle equivalent to `zipf_seam_sign_theory`, `zipf_seam_projection_theory`, or `zipf_seam_second_order_theory`. Under the no-fresh-computation migration rule, the decoupled-erf version of experiment 4 is therefore not reconstructible from existing historical outputs.

Planned handling unless you override:

| Output area | Plan |
|---|---|
| Canonical decoupled-erf outputs | Emit `BLOCKED.md` explaining that decoupled-erf seam theory is missing and cannot be migrated without a new derivation/run. |
| Legacy coupled-logistic outputs | Preserve the reconstructed historical seam-theory tables under `archive/legacy_coupled_logistic/` with manifest and provenance, but do not present them as current canonical decoupled-erf claims. |
| README audit section | Flag manuscript/claim-map rows that still reference coupled-logistic seam theory and need Phase 3 revision or a new decoupled-erf seam derivation. |

If you instead want the historical coupled-logistic seam theory migrated as canonical despite the updated model canon, approve that explicitly and I will emit canonical `outputs/seam_sign_checks.csv` and `outputs/aggregate_statistics.csv` from the three historical bundles with a strong AUDIT note.

## Historical Output Inventory

### `zipf_seam_sign_theory`

Saved aggregate values in `summary.json`:

| Historical key | Value | Proposed legacy row |
|---|---:|---|
| `local_match_counts_head10[0]` | 17 | `raw_seam_term_match_head10_linear` |
| `local_match_counts_head10[1]` | 17 | `raw_seam_term_match_head10_quadratic` |
| `local_match_counts_head10[2]` | 9 | `raw_seam_term_match_head10_cubic` |
| `local_match_counts_head50[0]` | 21 | `raw_seam_term_match_head50_linear` |
| `local_match_counts_head50[1]` | 13 | `raw_seam_term_match_head50_quadratic` |
| `local_match_counts_head50[2]` | 15 | `raw_seam_term_match_head50_cubic` |
| `local_full_match_head10` | 2 | `raw_seam_full_sign_match_head10` |
| `local_full_match_head50` | 2 | `raw_seam_full_sign_match_head50` |
| `smooth_match_counts_head50[0]` | 22 | `numerical_refit_match_head50_linear` |
| `smooth_match_counts_head50[1]` | 16 | `numerical_refit_match_head50_quadratic` |
| `smooth_match_counts_head50[2]` | 15 | `numerical_refit_match_head50_cubic` |
| `smooth_full_match_head50` | 9 | `numerical_refit_full_sign_match_head50` |

The file also contains per-corpus coefficient/sign arrays that can feed a row-level legacy `seam_sign_checks.csv`.

### `zipf_seam_projection_theory`

Saved aggregate values in `summary.json`:

| Historical key | Value | Proposed legacy row |
|---|---:|---|
| `local_baseline_head50_term_matches[0]` | 21 | duplicate check for `raw_seam_term_match_head50_linear` |
| `local_baseline_head50_term_matches[1]` | 13 | duplicate check for `raw_seam_term_match_head50_quadratic` |
| `local_baseline_head50_term_matches[2]` | 15 | duplicate check for `raw_seam_term_match_head50_cubic` |
| `local_baseline_head50_full_matches` | 2 | duplicate check for `raw_seam_full_sign_match_head50` |
| `smooth_exact_head50_term_matches[0]` | 22 | duplicate check for `numerical_refit_match_head50_linear` |
| `smooth_exact_head50_term_matches[1]` | 16 | duplicate check for `numerical_refit_match_head50_quadratic` |
| `smooth_exact_head50_term_matches[2]` | 15 | duplicate check for `numerical_refit_match_head50_cubic` |
| `smooth_exact_head50_full_matches` | 9 | duplicate check for `numerical_refit_full_sign_match_head50` |
| `projected_head50_term_matches[0]` | 20 | `projected_vs_empirical_match_head50_linear` |
| `projected_head50_term_matches[1]` | 15 | `projected_vs_empirical_match_head50_quadratic` |
| `projected_head50_term_matches[2]` | 18 | `projected_vs_empirical_match_head50_cubic` |
| `projected_head50_full_matches` | 9 | `projected_vs_empirical_full_match_head50` |
| `projected_vs_exact_smooth_term_matches[0]` | 23 | `projected_vs_numerical_match_head50_linear` |
| `projected_vs_exact_smooth_term_matches[1]` | 18 | `projected_vs_numerical_match_head50_quadratic` |
| `projected_vs_exact_smooth_term_matches[2]` | 16 | `projected_vs_numerical_match_head50_cubic` |
| `projected_vs_exact_smooth_full_matches` | 14 | `projected_vs_numerical_full_match_head50` |

The file also contains per-corpus local, projected, smooth, and empirical sign vectors.

### `zipf_seam_second_order_theory`

Saved aggregate values in `summary.json`:

| Historical key | Value | Proposed legacy row |
|---|---:|---|
| `first_order_head50_term_matches[0]` | 20 | duplicate check for `projected_vs_empirical_match_head50_linear` |
| `first_order_head50_term_matches[1]` | 15 | duplicate check for `projected_vs_empirical_match_head50_quadratic` |
| `first_order_head50_term_matches[2]` | 18 | duplicate check for `projected_vs_empirical_match_head50_cubic` |
| `first_order_head50_full_matches` | 9 | duplicate check for `projected_vs_empirical_full_match_head50` |
| `exact_smooth_head50_term_matches[0]` | 22 | duplicate check for `numerical_refit_match_head50_linear` |
| `exact_smooth_head50_term_matches[1]` | 16 | duplicate check for `numerical_refit_match_head50_quadratic` |
| `exact_smooth_head50_term_matches[2]` | 15 | duplicate check for `numerical_refit_match_head50_cubic` |
| `exact_smooth_head50_full_matches` | 9 | duplicate check for `numerical_refit_full_sign_match_head50` |
| `second_order_head50_term_matches[0]` | 17 | `second_order_match_head50_linear` |
| `second_order_head50_term_matches[1]` | 15 | `second_order_match_head50_quadratic` |
| `second_order_head50_term_matches[2]` | 16 | `second_order_match_head50_cubic` |
| `second_order_head50_full_matches` | 5 | `second_order_full_sign_match_head50` |
| `second_order_vs_smooth_term_matches[0]` | 20 | optional diagnostic row, not currently claim-map canonical |
| `second_order_vs_smooth_term_matches[1]` | 20 | optional diagnostic row, not currently claim-map canonical |
| `second_order_vs_smooth_term_matches[2]` | 18 | optional diagnostic row, not currently claim-map canonical |
| `second_order_vs_smooth_full_matches` | 13 | optional diagnostic row, not currently claim-map canonical |

The file also contains per-corpus second-order and quadratic-surrogate sign vectors.

## Proposed Legacy Consolidation Rules

If coupled-logistic legacy migration is approved:

1. Join per-corpus rows across all three `summary.json` files by normalized corpus name.
2. Preserve one row per corpus in `archive/legacy_coupled_logistic/seam_sign_checks.csv`.
3. Emit one aggregate row per claim-map metric in `archive/legacy_coupled_logistic/aggregate_statistics.csv`.
4. Use `zipf_seam_projection_theory` as the precedence source for first-order projection metrics, because it is the dedicated projection bundle.
5. Use `zipf_seam_second_order_theory` as the precedence source for second-order metrics, because it is the dedicated second-order bundle.
6. Treat repeated local/smooth counts across bundles as consistency checks. If duplicates disagree by any amount, emit an `AUDIT` warning rather than selecting silently.
7. Do not add non-claim optional diagnostics to canonical `aggregate_statistics.csv`; place them in a separate legacy diagnostics CSV if needed.

## Planned Canonical Output Under Current Decoupled-Erf Contract

Because the decoupled-erf seam-theory derivation is missing, the experiment directory should contain:

| Path | Status | Contents |
|---|---|---|
| `README.md` | to create | Standard schema with an AUDIT/BLOCKED section explaining model-canon mismatch |
| `BLOCKED.md` | to create | Exact missing source: decoupled-erf analytic seam sign/projection/second-order outputs |
| `manifest.json` | to create | Marks canonical outputs as blocked and records legacy source bundles |
| `source_config.py` | to create | Declares legacy bundle sources and canonical blocked status |
| `run_experiment.py` | to create | Thin wrapper that refuses canonical migration unless a decoupled-erf source is supplied |
| `archive/legacy_coupled_logistic/*` | optional on approval | Reconstructed historical coupled-logistic tables |

## Claim-Map Impact

`results/MANUSCRIPT_CLAIM_TO_CSV_MAP_v4.md` maps Section 3.8 and related discussion claims to experiment 4 rows in `4/outputs/aggregate_statistics.csv` and `4/outputs/seam_sign_checks.csv`. Those rows refer to the coupled-logistic historical seam theory. Under the decoupled-erf canon, these claims are not currently satisfied by canonical outputs.

Affected claim families include:

| Claim family | Claim-map destination | Current migration status |
|---|---|---|
| Raw/local sign counts | `4/outputs/aggregate_statistics.csv` | Available only for coupled-logistic |
| Projection-vs-empirical counts | `4/outputs/aggregate_statistics.csv` | Available only for coupled-logistic |
| Projection-vs-numerical-refit counts | `4/outputs/aggregate_statistics.csv` | Available only for coupled-logistic |
| Second-order sign counts | `4/outputs/aggregate_statistics.csv` | Available only for coupled-logistic |
| Per-corpus sign checks | `4/outputs/seam_sign_checks.csv` | Available only for coupled-logistic |

## Approval Decision Needed

Please choose one before I execute experiment 4 migration:

| Option | Meaning | Consequence |
|---|---|---|
| A. Enforce decoupled-erf canon | Do not emit coupled-logistic seam theory as canonical | Experiment 4 is `BLOCKED`; legacy data may be archived only |
| B. Legacy-as-canonical with AUDIT | Emit historical coupled-logistic seam outputs into canonical `outputs/` | Phase 2 claim map can be satisfied, but model canon is inconsistent and manuscript must flag/update |
| C. Dual-output | Emit canonical `BLOCKED.md` plus legacy tables under `archive/legacy_coupled_logistic/` | Safest for audit package; current claims remain unsatisfied until a decoupled-erf seam derivation exists |

My recommendation is Option C unless the immediate goal is to close the v4 claim map mechanically. It preserves the historical numbers without laundering them into the decoupled-erf canon.
