# f12 — composition robustness: is the seam a text-mixing artifact?

**Question.** Williams, Bagrow, Danforth & Dodds (PRE 91, 052811, 2015) argue
that two-regime structure in rank-frequency distributions is an artifact of
aggregating texts ("text mixing"), with a mixing-induced scaling break at rank
b ≈ N_avg (the mean per-constituent vocabulary) — and that the core-lexicon
reading of two regimes should be abandoned as corpus-relative. If they are right
about *our* object, the seam and its width law s ≈ 0.012·V (manuscript §3.3)
should depend on corpus composition.

**Two tests.**

## f12b — single-vs-aggregate split (observational, existing panel)

`scripts/run_f12b_composition_split.py` reclassifies the 25-corpus English panel
(from f2's canonical fits) by composition:

- **aggregate** (10): Shakespeare (≈37 plays), Complete Sherlock Holmes,
  Grimm's Fairy Tales, Aesop's Fables, Canterbury Tales, Arabian Nights,
  Federalist Papers (85 essays, 3 authors), King James Bible (66 books, many
  authors), Collected Poe, Dubliners (15 stories).
- **single** (15): the remaining single continuous works (one author, one
  narrative/argument).

Result (`outputs/f12b_composition_split.txt`): mean s/V 0.01193 (single) vs
0.01128 (aggregate); Welch t p = 0.190; Mann-Whitney p = 0.174; bootstrap 95%
CI on the mean difference [−0.00152, +0.00024]. The groups interleave fully
when sorted (lowest s/V is an aggregate, highest is a single work). The
insignificant trend points opposite to mixing-inflation.

## f12 — forced mixing (interventional)

`scripts/run_f12.py` takes 14 single continuous works (list in the script; no
collections), shuffles them once (seed 20260722), and fits — with one identical
erf-gate 9-parameter fitter (f5b's, k upper bound 5000) — every individual work
(m=1), 7 disjoint pairs (m=2), 3 disjoint quadruples (m=4), 2 disjoint
septuples (m=7), and the full concatenation (m=14). 27 fits.

Predictions in collision:
- **Williams:** the fitted break should migrate toward N_avg (~10⁴) and the
  break should grow more severe with m.
- **s-law:** s/V ≈ 0.012 at every m, with k a couple of orders below N_avg.

Result (run 2026-07-22, 27/27 fits, `outputs/f12_mixtures.csv`,
`outputs/f12_summary.md`): **the s-law wins on both axes.**

- median s/V by m = 1, 2, 4, 7, 14: **0.0119, 0.0127, 0.0128, 0.0123, 0.0120**
  — no trend across a fourteen-fold change in aggregation; the full 14-work,
  63,541-type concatenation sits dead-centre of the single-work distribution.
- median k/N_avg: 0.021 (m=1) → 0.079 (m=14) — the seam centre stays 1–2
  orders of magnitude below the Williams break location at every m; its mild
  growth tracks V (which keeps growing) rather than N_avg (which saturates).
- pooled regression over all 27 fits: s ∝ V^0.977, R² = 0.983 — the tightest
  width-law measurement in the project; forced mixtures land ON the line.

The m=1 fits also serve as an independent-reimplementation check of f2's s
values (they agree to the reported precision).

**Interpretation note.** Williams' break and our seam are different objects at
different scales: their break lives at rank ≈ N_avg; the seam sits at
k ~ 10²–10³ inside what their analysis treats as the unbroken Zipf/Simon head.
These tests show the seam and its width law are indifferent to composition —
whatever text mixing does to the deep tail, it neither creates nor moves the
seam.

**Verdict.** The seam is not a text-mixing artifact. Cited in manuscript §3.3
("Composition robustness") and §4.3.
