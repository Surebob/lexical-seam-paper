# Changes from MANUSCRIPT_DRAFT_v5 to MANUSCRIPT_DRAFT_v5_1

Date: 2026-04-19

This is a surgical manuscript update documenting the T1 branch-structure test for the multilingual `x^x - sqrt(x)` winner.

## Source Checked

- T1 outputs: `phase2_addon/t1_branch_test/`
- Symbolic verification: `phase2_addon/t1_branch_test/t1_symbolic_verification.txt`
- Numerical verification: `phase2_addon/t1_branch_test/t1_numerical_verification.csv`
- Summary: `phase2_addon/t1_branch_test/t1_summary.json`
- Report: `phase2_addon/t1_branch_test/t1_report.md`

## Manuscript Edit

- Updated Section 3.10 to remove the claim that Russian, Mandarin, and Arabic select a non-Bregman expression.
- Added the T1 result: `f(x) = x^x - sqrt(x) - 0.5(x - 1)` satisfies `f(1) = 0`, `f'(1) = 0`, and `f''(x) > 0` on `[0.05, 1.0]`.
- Reframed the multilingual split as a within-Bregman-class branch variation after linear recentering, with the dual-coordinate interpretation of the `0.5` shift left open.
- Tightened the Section 3.10 conservative summary to distinguish the canonical IS/exponential branches from the broader linearly-recentered Bregman class.

## Scope

- No other scientific sections were revised.
- Abstract and Introduction did not contain an explicit non-Bregman claim requiring correction, so they were left unchanged.

## Build Outputs

- New manuscript Markdown: `MANUSCRIPT_DRAFT_v5_1.md`
- New LaTeX builder: `build_manuscript_v5_1_latex.py`
- Generated LaTeX root copy: `MANUSCRIPT_DRAFT_v5_1.tex`
- Generated PDF: `MANUSCRIPT_DRAFT_v5_1.pdf`
- Generated LaTeX working directory: `results/manuscript_v5_1_latex/`
