# Experiment 7a Consolidation Plan: Canonical PMF Family

## Research Question

Within the Seam-Mandelbrot PMF family, what is the canonical held-out tradeoff among the free, fixed-`k`, fixed-`k,w`, soft-`k`, and soft-`k,w` variants under the train/test protocol used in the paper?

## Required Approval Gate

This experiment has five source bundles and a known split-fit/full-refit provenance risk. Per the migration rules, this plan must be approved before execution.

## Source Bundles Inspected

| Source bundle | Files inspected | Row count | Role |
|---|---|---:|---|
| `results/zipf_seam_mandelbrot_pmf` | `summary.json`, `seam_mandelbrot_table.csv`, `report.md` | 25 | Free PMF family baseline: Zipf, ZM, MOEZipf, free Seam-Mandelbrot |
| `results/zipf_seam_mandelbrot_regularized` | `summary.json`, `regularized_seam_table.csv`, `report.md` | 25 | Fixed-`k` and fixed-`k,w` variants compared against free seam |
| `results/zipf_seam_mandelbrot_softk_splitfit` | `summary.json`, `softk_table.csv`, `report.md` | 25 | Soft-`k` split-fit sweep and full-refit diagnostics |
| `results/zipf_seam_mandelbrot_softkw` | `summary.json`, `softkw_table.csv`, `report.md` | 25 | Soft-`k,w` extension and repeated soft-`k` values |
| `results/zipf_v4_verification` | `table_a_fourway_pmf.csv`, `report.md`, `pmf_winner_counts.png` | 25 | Final v4 four-way PMF verification table; no `summary.json` |

## Output Contract

Per `results/CONSOLIDATION_PROPOSAL_v3.md` and `results/MANUSCRIPT_CLAIM_TO_CSV_MAP_v4.md`, experiment 7a should emit:

| Canonical output | Purpose |
|---|---|
| `outputs/splitfit/table4_fourway.csv` | Table 4 body: Zipf/ZM/MOE/soft-`k` held-out average NLL, winner, `lambda_k` |
| `outputs/splitfit/pmf_variant_per_corpus.csv` | Per-corpus free/fixed/soft variant details |
| `outputs/splitfit/aggregate_statistics.csv` | Claim-facing Table 4 and PMF family aggregates |
| `outputs/fullrefit/fourway_per_corpus.csv` | Mirrored full-refit diagnostics, not canonical held-out claims |
| `outputs/fullrefit/aggregate_statistics.csv` | Full-refit diagnostic aggregates |
| `manifest.json` | Provenance, schemas, row counts |

The README must state clearly that canonical held-out manuscript metrics come only from `outputs/splitfit/`.

## Split-Fit vs Full-Refit Rule

Canonical held-out metrics:

- Use split-fit `test_avg_nll` fields only.
- Do not use full-refit training NLLs as held-out NLLs.
- Do not mix `best_lambda.selection.*` and `best_lambda.full_refit.*` in the same canonical held-out table.

Full-refit outputs:

- Preserve full-refit values as diagnostics in `outputs/fullrefit/`.
- Mark them non-canonical for Table 4 and held-out claims.

## Source-to-Output Mapping

### `outputs/splitfit/table4_fourway.csv`

Planned columns:

| Output column | Source |
|---|---|
| `slug`, `corpus`, `token_count`, `vocabulary_size` | `zipf_v4_verification/table_a_fourway_pmf.csv`, cross-checked against `zipf_seam_mandelbrot_pmf` |
| `zipf_test_avg_nll` | `zipf_seam_mandelbrot_pmf/seam_mandelbrot_table.csv` |
| `zm_test_avg_nll` | `zipf_seam_mandelbrot_pmf/seam_mandelbrot_table.csv` |
| `moe_test_avg_nll` | `zipf_seam_mandelbrot_pmf/seam_mandelbrot_table.csv` |
| `softk_test_avg_nll` | **Decision needed; see soft-`k` drift below** |
| `winner_family` | recompute from the four emitted NLL columns |
| `best_lambda_k` | soft-`k` source chosen for `softk_test_avg_nll` |

### `outputs/splitfit/pmf_variant_per_corpus.csv`

Planned columns:

| Variant family | Source |
|---|---|
| Free Seam-Mandelbrot | `zipf_seam_mandelbrot_pmf` and repeated in `zipf_seam_mandelbrot_regularized` |
| Fixed-`k` | `zipf_seam_mandelbrot_regularized` |
| Fixed-`k,w` | `zipf_seam_mandelbrot_regularized` |
| Soft-`k` | `zipf_seam_mandelbrot_softk_splitfit` or `zipf_seam_mandelbrot_softkw` depending on approved precedence |
| Soft-`k,w` | `zipf_seam_mandelbrot_softkw` |

### `outputs/splitfit/aggregate_statistics.csv`

Required claim-map rows:

| Metric row | Source / derivation |
|---|---|
| `winner_count_zipf` | recompute from `table4_fourway.csv` |
| `winner_count_zm` | recompute from `table4_fourway.csv` |
| `winner_count_moe` | recompute from `table4_fourway.csv` |
| `winner_count_softk` | recompute from `table4_fourway.csv` |
| `softk_beats_zipf_count` | `softk_test_avg_nll < zipf_test_avg_nll` |
| `softk_beats_zm_count` | `softk_test_avg_nll < zm_test_avg_nll` |
| `softk_beats_moe_count` | `softk_test_avg_nll < moe_test_avg_nll` |
| `median_softk_minus_zipf` | median of `softk - zipf` |
| `median_softk_minus_zm` | median of `softk - zm` |
| `median_softk_minus_moe` | median of `softk - moe` |
| `free_pmf_heldout_winner_count_zipf` | `zipf_seam_mandelbrot_pmf/counts/test_winner_counts` or recompute |
| `free_pmf_heldout_winner_count_zm` | `zipf_seam_mandelbrot_pmf/counts/test_winner_counts` or recompute |
| `free_pmf_heldout_winner_count_moe` | `zipf_seam_mandelbrot_pmf/counts/test_winner_counts` or recompute |
| `free_pmf_heldout_winner_count_seam` | `zipf_seam_mandelbrot_pmf/counts/test_winner_counts` or recompute |
| `softk_beats_free_count` | `zipf_seam_mandelbrot_softk_splitfit/counts/softk_beats_free_test`, unless soft-`k` precedence changes |
| `median_softk_minus_free` | `zipf_seam_mandelbrot_softk_splitfit/medians/delta_best_vs_free_test`, unless soft-`k` precedence changes |
| `softk_beats_fixedk_count` | `zipf_seam_mandelbrot_softk_splitfit/counts/softk_beats_fixedk_test` |
| `median_softk_minus_fixedk` | `zipf_seam_mandelbrot_softk_splitfit/medians/delta_best_vs_fixedk_test` |
| `softk_step2_help_count` | `zipf_seam_mandelbrot_softk_splitfit/counts/softk_step2_help_count` or recompute from approved soft-`k` source |
| `softk_step2_nonhelp_count` | `25 - softk_step2_help_count` |
| `softk_step2_gain_median` | `zipf_seam_mandelbrot_softk_splitfit/medians/softk_step2_gain` or approved soft-`k` source |
| `free_pmf_step2_help_count` | `zipf_seam_mandelbrot_pmf/counts/seam_step2_help_count` |
| `softk_parameter_count` | protocol constant: 6 |

## Soft-`k` Drift Requiring Approval

The final v4 verification table is not identical to `zipf_seam_mandelbrot_softk_splitfit/softk_table.csv` for three corpora:

| Corpus slug | `zipf_v4_verification` soft-`k` | `softk_splitfit` soft-`k` | Difference |
|---|---:|---:|---:|
| `critique_of_pure_reason` | 5.95290862285357 | 5.952965161503562 | -0.00005653865 |
| `jane_eyre` | 6.643433234761215 | 6.642788574812209 | +0.00064465995 |
| `king_james_bible` | 6.016015736044842 | 6.016425262838996 | -0.00040952679 |

`zipf_v4_verification/table_a_fourway_pmf.csv` matches:

- Zipf/ZM/MOE columns exactly from `zipf_seam_mandelbrot_pmf`.
- Soft-`k` values from `zipf_seam_mandelbrot_softkw/softkw_table.csv`'s repeated `softk_test_avg_nll` column for the three drifting rows.

This creates a provenance decision:

| Option | Meaning | Effect |
|---|---|---|
| A. Use `zipf_v4_verification` as canonical Table 4 target | Treat v4 verification as final hand/producer-assembled join | Preserves manuscript/checklist values: winner counts ZM 11, soft-k 10, MOE 4, Zipf 0; medians `-0.001949314592`, `-0.000584460494`, `-0.039357582849` |
| B. Use `zipf_seam_mandelbrot_softk_splitfit` as canonical soft-`k` source | Treat the named split-fit bundle as the authoritative soft-`k` source | Recompute Table 4 and flag DRIFT against v4 verification/manuscript for at least three rows |
| C. Use `zipf_seam_mandelbrot_softkw` repeated soft-`k` values | Treat soft-`k,w` bundle as later or corrected soft-`k` source for those rows | Same Table 4 soft-`k` values as v4 verification, with direct upstream source for the soft-`k` column |

My recommendation is **Option A with explicit upstream decomposition**:

- Emit `table4_fourway.csv` from `zipf_v4_verification/table_a_fourway_pmf.csv`.
- Add provenance columns or sidecar `table4_provenance.csv` stating that Zipf/ZM/MOE fields trace to `zipf_seam_mandelbrot_pmf`, while the soft-`k` field matches `zipf_seam_mandelbrot_softkw` for the three rows where `softk_splitfit` differs.
- Preserve `softk_splitfit` values in `pmf_variant_per_corpus.csv` with an audit column `softk_splitfit_vs_table4_delta`.
- Add an AUDIT note in README documenting the three-row drift.

This avoids reintroducing the split/full or hand-join opacity problem while still making the current paper's Table 4 cell values reproducible from one CSV.

## Precedence Rules If Option A Is Approved

1. `table4_fourway.csv` is emitted from `zipf_v4_verification/table_a_fourway_pmf.csv` as the canonical rendered Table 4 body.
2. `table4_provenance.csv` records source provenance per column and drift against `softk_splitfit`.
3. `pmf_variant_per_corpus.csv` is the union of all variant outputs:
   - free, Zipf, ZM, MOE from `zipf_seam_mandelbrot_pmf`
   - fixed-`k`, fixed-`k,w` from `zipf_seam_mandelbrot_regularized`
   - soft-`k` splitfit from both `softk_splitfit` and v4/softkw repeated source where applicable
   - soft-`k,w` from `zipf_seam_mandelbrot_softkw`
4. `aggregate_statistics.csv` Table 4 winner counts and medians are recomputed from the emitted `table4_fourway.csv`, not copied by hand.
5. Variant-specific regularization metrics are computed from their owning bundles, with source bundle named in a `source_bundle` column when useful.
6. `outputs/fullrefit/*` is populated only from `best_lambda.full_refit` and analogous saved diagnostics; it is clearly marked non-canonical for held-out claims.

## Union-vs-Select Rules

Use **union** for:

- variant diagnostics across free/fixed/soft/softkw
- step-2 help/gain diagnostics for each variant
- lambda and parameter fields
- full-refit diagnostics

Use **select** for:

- Table 4 canonical four-way held-out body
- Table 4 winner counts and median advantages
- manuscript checklist rows 644-650

Selection should follow the approved soft-`k` precedence option above.

## Known Non-Goals

- Do not import head-window results here; those belong to experiment 7b.
- Do not import metadata predictability of `lambda_k`; that belongs to 7c.
- Do not import asymmetric gate results; that belongs to 7d.
- Do not import hierarchical pooling results; that belongs to 7e.
- Do not delete or hide any pre-splitfit PMF sources; preserve them as provenance/archive where they contribute.

## Approval Needed

Please approve one soft-`k` precedence option before I execute migration:

- **A:** v4 verification table as canonical Table 4 target, with explicit provenance and soft-`k` drift audit.
- **B:** named `softk_splitfit` bundle as canonical, with DRIFT against v4/manuscript for three rows and possible aggregate shifts.
- **C:** `softkw` repeated soft-`k` values as canonical soft-`k` source, with v4 verification as a matching rendered table.

I recommend **A** because the paper/checklist Table 4 claims already cite the v4 verification values, and the migration can make that join explicit rather than implicit.
