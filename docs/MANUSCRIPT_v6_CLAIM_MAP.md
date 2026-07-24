# Manuscript v6 claim-to-source map (2026-07-22)

Rule inherited from the April program: every cited number lives in exactly one
canonical file. Paths relative to repo root. Bootstrap CIs: `docs/v6_bootstrap_cis.md`.

## §3.1 — residual reproducible; label scoring-conditioned
| claim | source |
|---|---|
| step-2 winners 25/25; terminal 25/25 (3.70–25.29%, mean 11.02%); both forms in beam | `experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv` |
| exact independent reproduction (Δc=0, winners 25/25) | `experiments/f1_fresh_reproduction/outputs/f1_summary.md` |
| free-amplitude winner map (14/25 change; exp→0; IS amplitudes 0.72–0.97 high-c) | `experiments/f1_fresh_reproduction/outputs/f1_per_corpus.csv` |
| weighted-objective flips (11/11 high-c IS→exp; stable 14/25) | `experiments/f9_legacy_reruns/outputs/f9_wls_winners.csv` |
| top-100 winner map (xpow 20/25) | `experiments/f3_lambda_zm_vs_moezipf/outputs/f3_summary.md` |
| widened-grammar span 17/17 R²>0.975 | canonical `experiments/10b_grammar_width_robustness/outputs` |
| identity e^(x−1)=e^(−0.95) r^(0.95/lnV) | algebraic; stated in `experiments/f1_.../README.md` |

## §3.2 — identity evidence
| claim | source |
|---|---|
| top-100-only c=5.3 RMSE 0.044; ablation behaviour | `experiments/2a_function_word_ablation/outputs/function_word_ablation.csv` |
| city control (33,535; c=100.4; 0.1206 vs 0.1225) | `experiments/1b_non_linguistic_control/outputs/city_population_control.csv` |
| surname control (diagnostic silent; ZM RMSE 0.031) | `experiments/f5_expanded_panel/outputs/f5a_panel.csv` |
| deterministic mixtures: small-gap seam; pure power law RMSE 0.0000; amplitude-gap suffices | `experiments/f11_loose_ends/outputs/f11_2b.csv` (+ canonical `2b_*`) |

## §3.3 — absorption, gate, width law
| claim | source |
|---|---|
| smooth beats ZM/piecewise 25/25; absorption 0/25 vs MOE 19/25 (+0.017); anti-IS MOE residual | `experiments/2c_mechanism_absorption_residual_space/outputs/aggregate_statistics.csv` |
| erf 24/25, arctan 1 (Dubliners), spread 677; logistic synthetic recovery 25/25 | `experiments/3e_gate_family_bic_sweep/outputs` |
| k identified ±1%; k~V β=0.758 CI [0.367,1.149] R²=0.386 | `experiments/f2_k_profile_likelihood/outputs/f2_summary.md` |
| s-law eq(3): β=1.003 CI [0.928,1.079] R²=0.967; s/V 0.0118 | `f2_.../f2_summary.md` + `experiments/f2b_s_profile/outputs/f2b_summary.md` |
| register points 0.0120–0.0122; language interiors 0.0068–0.0096; surnames 0.0266 | `experiments/f5_expanded_panel/outputs/f5b_gate_fits.csv` + `experiments/f5c_degenerate_refits/outputs/f5c_refits.csv` |
| composition split: single 0.0119 vs aggregate 0.0113, Welch p=0.19, boot CI [−0.0015,+0.0002] | computed from `f2_per_corpus.csv` (script in f12 README; classification list there) |
| forced mixing m=1..14: s/V stays ~0.012, k ≪ N_avg | `experiments/f12_forced_mixing/outputs/f12_mixtures.csv`, `f12_summary.md` |
| width = shared invariant: GA 0.0122, coreless Simon 0.0118, classic Simon 0.0101, surnames 0.0266 | `experiments/f14_ga_break/outputs/f14_ga_fits.csv`, `f14b_nulls.csv` |
| cross-language: corr(depth, s/V)=+0.94; deep median 0.0111; 11 languages | `experiments/f15_deep_multilingual/outputs/f15_deep_fits.csv`, `f15b_panel.csv` |
| calibration: bootstrap recovery median 1.032 IQR [1.006,1.062], 21/25 | `experiments/f16_nonparametric_width/outputs/f16d_bootstrap.csv` |
| model-light envelope: language wider than twins 16/25, median +0.011 | `experiments/f16_nonparametric_width/outputs/f16_widths.csv` |
| within-corpus depth walks: corr +0.93, collapse 0.37, 4 langs one curve | `experiments/f15_deep_multilingual/outputs/f15d_approach.csv` |
| matched-depth concats: 5 langs tok/V>=15 median 0.0122 vs EN 0.0118 | `experiments/f15_deep_multilingual/outputs/f15c_concat.csv` |

## §3.4 — c is sampling depth
| claim | source |
|---|---|
| collapse trajectories; 6/6 winner ladder | `experiments/f6_c_sampling_depth/outputs/f6_slices.csv` |
| binomial ≈ prefix | `f6_.../outputs/f6_summary.md` §C |
| matched-size panel (EN median 6.3 vs non-EN 0.00; W&P pair 10.1 vs 0.0) | `experiments/f8_matched_size_panel/outputs/f8_matched_panel.csv` |
| winner classifier 80% LOO | `experiments/f11_loose_ends/outputs/f11_summary.md` §D |

## §3.5 — generative account
| claim | source |
|---|---|
| PLN mixture fits (π_H median 2.15%, etc.) | `experiments/f4d_poisson_lognormal_mixture/outputs/f4d_mixture_fits.csv` |
| b prediction r=0.941, 4%; c scale corr 0.82; latent-quantile failure b≈8–12 | `experiments/f7_histogram_predicts_curve/` (run log in commit a588ab4; scripts reproduce) |
| downward trajectories + V(T) 1–2%; upward table (V/b/c) | `experiments/f6b_heavy_tail_extrapolation/outputs/f6b_trajectories.csv`, `f6b_summary.md` |
| identifiability: ΔNLL≤5 basins σ_T 1.9–4.0, n_total 20k–600k; ML ≠ best predictor | `experiments/f4f_identifiability_audit/outputs/f4f_landscape.csv` |
| histogram demands 2 comp only 3/10 broad-basin | `experiments/f4g_two_component_reaudit/outputs/f4g_summary.md` |
| Heaps β median 0.476; corr(β,1/b)=0.991 | `experiments/f11_loose_ends/outputs/f11_heaps.csv` |
| mixture contains width law: 16/25 at ML (median 1.02); 8/9 heal at pred-basin (median 1.01) | `experiments/f13_mixture_width_law/outputs/f13_summary.md`, `f13b_summary.md` |

## §3.6 — usage boundary
| claim | source |
|---|---|
| k_gate/k_lab median 3.10; corr 0.20 / −0.43; agreement 89.9–91.8 / 93.3–97.0 | `experiments/f10_gate_vs_labels/outputs/f10_per_corpus.csv`, `f10_summary.md` |
| labeled k_POS β=0.298 [0.076,0.519] | same |
| W&P→Bible transfer 0.1498 vs 0.1494 vs 0.1882 | `experiments/f9_legacy_reruns/outputs/f9_poly_transfer.csv` |
| labels + pilot | `data/annotations/` |

## §3.7 — λ-ZM
| claim | source |
|---|---|
| 25/25 BIC vs ZM+MOE (both fits); median −0.0157 vs MOE | `experiments/f3_lambda_zm_vs_moezipf/outputs/f3_english.csv` |
| median improvement 10.2% (2.6–23.0) | `experiments/f1_fresh_reproduction/outputs/f1_summary.md` |
| panel 10/10; languages 6.4–22.5% | `f3_.../f3_multilang.csv`, `f5a_panel.csv` |
| held-out 64/64 (EN 50/50 med +8.5%, ML 14/14 med +13.4%, worst +2.1%) | `experiments/f3_lambda_zm_vs_moezipf/outputs/f3b_heldout_summary.md`, `f3b_heldout.csv` |
| frozen-λ* equal-param: EN 25/25 (med +8.2%, ho 45/50, 95% retained); transfer 7/7 (ho 14/14, 99%; zh +19.4, ar +14.0, ru +12.7); λ* = 20.6 | `experiments/f19_parameter_free_lambda/outputs/f19_summary.md`, `f19_results.csv` |

## §3.7 — family bake-off
| claim | source |
|---|---|
| weight-class sweep: frozen beats all 3p (logn 0.1669, zm 0.1806), free beats all 4p (cubic 0.1638, cutoff 0.1761), only two-regime models above (hard5p 0.1490, canonical 9p 0.1094) | `experiments/f23_family_bakeoff/outputs/f23_summary.md`, `f23_results.csv` |

## §3.7 — likelihood-space transfer + scope
| claim | source |
|---|---|
| free λ-ZM PMF beats ZM 69/75, MOE 60/75 held-out NLL; λ_MLE median −7.35 scattered (no 2nd constant); frozen-at-MLE-LOO beats ZM 72/75, ties MOE 43/75 | `experiments/f22_pmf_likelihood/outputs/f22_summary.md`, `f22b_summary.md` |

## §3.3 — null calibration (two-sided instrument test)
| claim | source |
|---|---|
| seamless nulls: s/V median 0.152, IQR [0.064, 0.473], 2/75 in band, 27/75 pinned; paired reals 23/25 in band (median 0.0121) | `experiments/f18_null_width/outputs/f18_summary.md`, `f18_null.csv` |

## §4.2 — retired claims → killing experiments
1 k≈√V → `f2`; 2 WLS → `f9`; 3 σ_T universal → `f4f`; 4 morphology dial → `f4f` (A);
5 histogram 25/25 → `f4g`; 6 seam=grammar → `f10`; 7 √e → `f4d` vs degenerate-GMM note
(`f4b`/`f4d` READMEs); 8 convex conjugation → `docs/audits/gpt5_v1_audit_report.md`.

## §4.4 — neural-model observations
Companion-paper material (L-series: `L1_seam_in_embeddings/outputs/L1_summary.md`,
`L3_seam_in_learning_curves/outputs/L3_seam_fits_bounded.csv`), **not distributed
with the public Paper-1 repository** — §4.4 reports these as motivation only,
explicitly not as claims of this paper; they are substantiated in the companion
neural-models paper.

## Figures
fig1 residual (computed from pg100 inline) · fig2 ← f6_slices · fig3 ← f2_per_corpus
+ f5b · fig4 ← f7 pairs · fig5 ← f6b_trajectories · fig6 ← f11_heaps ·
fig7_fingerprints ← f2_per_corpus + f5b (registers, surnames) + f12_mixtures +
f15c_concat (tok/V ≥ 15 filled) + f14_ga_fits + f14b_nulls + x1 (belt, hardcoded
from x1_summary.md). Built by `paper/scripts/build_assets.py`.
