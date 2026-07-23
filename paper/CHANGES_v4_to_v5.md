# Changes from MANUSCRIPT_DRAFT_v4 to MANUSCRIPT_DRAFT_v5

Date: 2026-04-18

This document records the Phase 3 manuscript migration from the coupled 8-parameter logistic smooth model used in `MANUSCRIPT_DRAFT_v4.md` to the decoupled 9-parameter erf-gate canonical model used in `MANUSCRIPT_DRAFT_v5.md`.

## Sources Checked

- Current manuscript: `MANUSCRIPT_DRAFT_v4.md`
- Drift tracker: `DRIFT_TRACKER.md`
- 1a enriched-search outputs: `experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv`
- 3a smooth-fit outputs: `experiments/3a_smooth_two_regime_fits/outputs/aggregate_statistics.csv`
- Decoupled five-gate sweep: `results/s2_v3_windows_full_outputs_2026-04-18/s2_v3_per_corpus_results.csv`
- Decoupled gate aggregate statistics: `results/s2_v3_windows_full_outputs_2026-04-18/s2_v3_aggregate_statistics.csv`
- Synthetic gate-recovery validation: `results/s2_synthetic_recovery_outputs_workers32_2026-04-18/outputs/synthetic_gate_recovery_aggregate_statistics.csv`
- Decoupled erf k-scaling recheck: `phase2_addon/k_scaling_decoupled_erf_recheck.csv`

Note: no local `LAB_NOTEBOOK.md` was found under `/Volumes/External2TB/emlexperiment`, `/Users/gregkara`, or `/mnt/user-data` during the pre-revision search, so the revision used the tracker and saved result files above.

## Major Scientific Updates

- Replaced the canonical smooth-rank model with a decoupled 9-parameter erf-gate two-regime ZM model: `(a1, b1, c1, a2, b2, c2, k, w_gate, w_tail)`.
- Updated Section 2.6 to define the erf gate, independent `w_gate`, and independent `w_tail` tail-coordinate smoothing scale.
- Added gate-family specificity as a new pillar inside Section 3.5.
- Reported the decoupled five-gate result: erf wins `24/25`, arctan wins `1/25` on Dubliners, logistic and algebraic win `0/25`, tanh calibration passes `25/25`, median BIC spread `677.37`, mean spread `864.66`.
- Added the synthetic logistic-gate recovery validation: logistic-generated data recovers logistic on `25/25`, erf wins `0/25`, and erf never beats the true generator.
- Weakened Section 3.6 k-scaling: decoupled erf gives `β = 0.7584`, 95% CI `[0.3457, 1.1711]`, `R² = 0.3858`, so the prior tight `0.5214` coupled-logistic claim is no longer canonical.

## Corrected Drift

- Replaced the unsupported Section 3.1 prose-only improvement range `0.15% to 2.0%` with terminal enriched-search improvements `3.70% to 25.29%`, mean `11.02%`, and `25/25` corpora improved.
- Clarified that the step-2 diagnostic itself improves only `12/25` corpora, with positive gains `0.33% to 2.50%` and all-corpus median `-0.23%`.
- Replaced the unsupported Section 3.1 `0.003` mean winner-runner gap with verified gaps: step-2 mean `0.00201`, terminal final-step mean `0.00075`.
- Corrected the small-width coupled diagnostic corpus list from `Aesop's Fables, Grimm's, Principia Ethica, Dubliners` to `Grimm's Fairy Tales, Aesop's Fables, Dubliners, Critique of Pure Reason`.

## Analytic Seam Section

- Re-derived the gate derivative terms for the canonical erf gate:
  - `τ(0) = 1/2`
  - `τ′(0) = 1/(sqrt(pi) * w_gate)`
  - `τ″(0) = 0`
  - `τ‴(0) = -2/(sqrt(pi) * w_gate^3)`
- Preserved the exact product-rule expansion and tangent-space projection formula.
- Removed the old coupled-logistic universal `(−, −, −)` sign-pattern claim as canonical evidence. The erf third-order term can oppose the old sign contribution depending on the fitted regime gap, so the sign-count validation must be rerun before being claimed for the decoupled erf model.

## Claims Intentionally Softened

- Low-c parameter-correlation claims involving "width" are now labeled as coupled-family diagnostics because no decoupled-erf parameter-correlation rerun was present in the supplied canonical outputs.
- PMF soft-k regularization still reports the historical `α = 0.521` PMF regularizer because those PMF experiments were not rerun under the decoupled erf k-scaling estimate.
- The manuscript no longer frames statistical transition-centre scaling as a precise physical law.

## Build Outputs

- New manuscript Markdown: `MANUSCRIPT_DRAFT_v5.md`
- New LaTeX builder: `build_manuscript_v5_latex.py`
- Generated LaTeX root copy: `MANUSCRIPT_DRAFT_v5.tex`
- Generated LaTeX working directory: `results/manuscript_v5_latex/`
- Generated PDF: `MANUSCRIPT_DRAFT_v5.pdf`
- Updated drift tracker: `DRIFT_TRACKER.md`

PDF compilation succeeded with bundled TinyTeX/XeLaTeX. The final build had layout warnings only; missing-glyph warnings from the first pass were resolved by replacing unsupported inline exponent glyphs with ASCII scientific notation in prose.
