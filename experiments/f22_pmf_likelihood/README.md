# f22 / f22b — likelihood-space bake-off: lambda-ZM in MOEZipf's home arena (2026-07-24)

**Question:** the paper's comparisons are rank-curve LSQ. MOEZipf's native
objective is maximum likelihood on token counts. Does the seam correction
survive the arena change — and does the universal amplitude?

**Protocol (the old canon's):** per corpus, 80/20 binomial split of type
counts (x3, fixed seeds); ranks/support from train; every model truncated-
normalized over [1, V_train]; score = held-out NLL per test token. In PMF
form the amplitude is absorbed by the normalizer, so ZM, MOEZipf, and
frozen-lambda-ZM carry exactly TWO free parameters each; free lambda-ZM
carries three. `scripts/run_f22.py`, `run_f22b.py`.

## f22 — transplanting the rank-curve constant FAILS (instructively)

- free lambda-ZM beats ZM 69/75, beats MOEZipf 60/75 (median +0.0061 /
  +0.0017 nats/token; the one-extra-param BIC cost is ~3e-5 — negligible).
- frozen at the RANK-CURVE constant 20.6: loses 6/75 to ZM. **Constants are
  objective-relative** — MLE weights by token mass (head-dominated), LSQ
  weights every rank equally (tail-dominated); their optimal amplitudes
  differ.

## f22b — is there a second universal constant in likelihood space? NO

- lambda_MLE across 25 corpora: median **-7.35**, IQR [-11.3, -3.8], range
  [-18.9, +16.3] — different sign from the rank-curve amplitude, and far
  more dispersed. No clean second constant.
- Yet the SHAPE still wins: frozen at the LOO median, it beats ZM **72/75**
  at equal parameter count (median +0.0053 nats/token), and ties MOEZipf
  (43/75, median +0.0003 = statistical tie).

## Honest synthesis

1. The seam correction helps under every objective tested: rank-curve
   (dominant), token-likelihood (clear vs ZM, tie-to-win vs MOE).
2. The AMPLITUDE'S UNIVERSALITY is a rank-curve property. lambda* = 20.6 is
   a constant of the equal-per-rank lens; the token-mass lens has its own,
   noisier, sign-flipped amplitude. Manuscript section 3.7 now states this
   scope explicitly (added 2026-07-24).
3. MOEZipf remains genuinely competitive in its home arena at equal
   complexity — consistent with the April program's finding — while the
   free-amplitude seam form edges it 60/75.

Open thread (Paper 2): why the sign flip — what does the truncated-
normalized token-mass functional do to the seam term's optimal loading?

Outputs: `outputs/f22_results.csv`, `f22_summary.md`, `f22b_results.csv`,
`f22b_summary.md`.
