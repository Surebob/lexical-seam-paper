# f11 — loose ends (2026-07-21)

Four items, one sweep. `scripts/run_f11.py`.

## A. 2b faithful rerun (deterministic legacy protocol) — CLOSED, with a finding
- Pure single power law (α=1.5): ZM RMSE **0.0000**, no correction helps — the
  control passes exactly as the manuscript intends.
- Legacy-construction "same-exponent control" (α=1.5/1.5 but scales 10 vs 1):
  winner=IS, helps=True — because **an amplitude gap alone creates the seam**;
  equal exponents with unequal scales is still a two-population structure.
- v6 actions: reword §2.8/§3.4 (the control must be equal-scale copies, matching
  the archived α=0.8 control), and add the insight: the seam requires only
  rate-scale separation between populations, not an exponent gap. (Also resolves
  f9-A: April's deterministic protocol was sound; sampled-rate constructions are
  a different object.)

## B. Diachronic first look — essentially null
corr(year, s/V) = −0.34 on English originals (n=15, years approximate); no strong
six-century drift in seam structure — consistent with register invariance.
Reported as a scoped observation, not a claim.

## C. Heaps' law joins the framework
Per-corpus Heaps β from binomial thinning: median 0.476 (range 0.39–0.63), and
**corr(β, 1/b) = 0.991** — the classical Zipf–Heaps duality confirmed almost
exactly on our panel. Combined with f6b (which predicts V(T) to 1–2%), the
framework now contains both famous laws, one measured, one predicted.

## D. Winner-from-stats classifier
LOO logistic on (log tokens, hapax ratio, TTR): **80% accuracy (20/25)** predicting
the is/exp winner family with no curve fitting — the winner label is largely a
function of raw corpus statistics, consistent with the sampling-depth ladder.
