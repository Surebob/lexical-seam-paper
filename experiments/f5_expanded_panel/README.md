# f5 — expanded panel (2026-07-21)

**Question:** do the seam, the λ-ZM formula, and the s ≈ 1.2%·V law survive outside
classic English literature — modern registers, more languages, and a second
non-linguistic control?

Panel: Brown (balanced 1961 American), Cornell movie dialogs (conversational),
WikiText-103 1M-token slice (modern encyclopedic), US Census 2010 surnames
(control), + Italian/Finnish/Polish/Portuguese/German/Swedish (Gutenberg,
150k-token slices). `scripts/fetch_corpora.py`, `run_f5a.py` (light pipeline),
`run_f5b.py` (erf gate fits).

## f5a results (light pipeline)

- **λ-ZM BIC wins 10/10** (lifetime 42/42 across all corpora ever tested).
  Registers: Brown +14.1%, WikiText +9.9%, Cornell +5.8%. Languages +9.0–22.5%.
- **Seam diagnostic (unit-amplitude step-2) is language-specific**: fires on 8/9
  language corpora, stays **silent on the surname control** (step2_helps=False;
  its λ-ZM "improvement" is free-parameter polish on an already near-perfect fit,
  ZM RMSE 0.031). Second non-linguistic control in the bag.
- Modern registers reproduce the deep-sampling regime: Brown c=402, WikiText c=451,
  IS winners — the high-c phenomenon is not a Gutenberg artifact.
- Cornell dialogs: c=88, diagnostic silent — a mid-band register datapoint
  consistent with f6's sampling-depth ladder.
- New languages: all low-c with xpow/exp winners — as expected for 33–150k-token
  slices per f6 (size, not only morphology; matched-size redo queued).

## f5b results (s-law on the panel)

Interior (non-degenerate) fits — **the law holds across registers**:

| corpus | s/V |
|---|---:|
| Brown | 0.0122 |
| WikiText-1M | 0.0120 |
| Cornell dialogs | 0.0121 |
| (English books, f2 reference) | 0.0118 |
| Italian | 0.0093 |
| German | 0.0088 |
| Finnish | 0.0075 |
| Census surnames (control) | 0.0266 |

1961 balanced prose, 2016 encyclopedia, and film dialog land within ±3% of the
English-book constant. Languages sit somewhat lower (0.0075–0.0093) — possibly a
real language/morphology effect, possibly slice-size (their corpora are 33–150k
tokens); matched-size comparison required before any cross-language universality
claim. The surname control fits a different constant (different system — fine).

**Degenerate fits (excluded, flagged):** Portuguese, Swedish, Polish returned
s/V > 3 with **w_gate pinned at the 0.05 lower bound and w_tail at the 10.0 upper
bound** — a step-gate escape optimum on the shallowest corpora (Polish: 2.8 tokens
per type). These are optimizer pathologies, not measurements; rerun queued with
tightened w_tail bound / more starts / deeper slices. The combined-regression line
in `f5b_summary.md` includes them and is therefore not meaningful as printed.

## Verdict

- λ-ZM: fully general (registers, languages, even the control).
- Seam diagnostic: language-specific (two controls now silent).
- s-law: **register-independent for English at ≈1.2% of V**; cross-language value
  pending matched-size + non-degenerate refits.
