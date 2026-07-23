# legacy_engines — frozen originals

The scripts that produced the canonical results, copied verbatim from the archive repo
root. **Do not edit these** — they are provenance. Fresh work rewrites them into a
tested package under `src/seam/` (see ROADMAP §4).

Engine → canonical experiment map:

| Engine | Feeds |
|---|---|
| `eml_tree_port.py` | validated EML trainer port (foundation; audio-era origin) |
| `eml_zipf_experiment.py` | corpus load, tokenize, ZM fit, residuals (imported by search engines) |
| `eml_zipf_enriched_search.py` | 1a per-corpus enriched search (+1b via corpus swap) |
| `eml_zipf_discrete_search.py` | earlier deterministic search engine (10a depth ablation lineage) |
| `eml_zipf_formula_report.py` | formula → math-notation reporting |
| `eml_city_population_experiment.py` | 1b non-linguistic control |
| `zipf_analysis_common.py`, `zipf_eval_protocol_utils.py` | shared fitting/eval utilities |
| `zipf_correct_model_reranked_all.py`, `zipf_reranked_model_all_corpora_relaxed.py` | 3a smooth fits |
| `zipf_moezipf_comparison.py`, `zipf_continuous_piecewise.py`, `zipf_bic_comparison.py` | 3b model-family BIC |
| `zipf_kstat_scaling.py` | 3c k-scaling |
| `zipf_pos_all_corpora.py`, `zipf_pos_manual_v2.py` | 3d POS crossover |
| `zipf_multilang*.py` | 6 multilingual |
| `zipf_synthetic_mixture.py` | 2b mixtures |
| `zipf_function_word_test.py` | 2a ablation |
| `zipf_english_gap_verify.py` | 1d Euclidean gap |
| `zipf_guard_ablation.py`, `zipf_wls_test.py` | 1c robustness (E2 rerun base) |
| `zipf_simulation_recovery.py` | 5a recovery (E3 rerun base) |
| `zipf_smooth_parameter_sweep.py` | 5c sweep (E4 rerun base) |
| `zipf_lowc_manifold_analysis.py`, `zipf_phase_coordinate_analysis.py` | 5b manifold/span |
| `zipf_seam_{sign,projection,second_order}_theory.py` | 4 seam theory (E6 rerun base) |
| `zipf_step10_ablation.py`, `zipf_head_poly_{decomposition,transfer}.py` | 10a depth ablation |
| `zipf_widened_grammar_{diagnostic,extended}.py` | 10b grammar widening |
| `zipf_seam_mandelbrot_{pmf,regularized,softk,softkw}.py` | 7a PMF family (Paper 2) |
| `zipf_angle1_head_windows.py` … `zipf_angle4_hierk.py` | 7b–7e PMF diagnostics (Paper 2) |
| `zipf_angle6_bible_books.py` | 8 Bible decomposition |
| `zipf_hybrid_headtail.py`, `zipf_fully_ours_variants.py`, `zipf_hybrid_mechpenalty.py` | 9 likelihood frontier |

The decoupled gate-family engine (3e) lives separately in `src/s2_decoupled/`
(shared modules + `run_s2_v3_windows.py`; corpora stripped — point it at `data/zipf/`).
