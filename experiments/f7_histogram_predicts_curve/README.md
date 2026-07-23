# f7 — the histogram predicts the curve (2026-07-21)

**Claim demonstrated:** the 5-parameter Poisson-lognormal mixture fitted ONLY to a
corpus's type-count histogram (f4d; never sees ranks) predicts the corpus's fitted
Zipf-Mandelbrot parameters through the observation model (draw N_total latent rates
from the mixture — N_total inferred from zero-truncation — Poisson counts, drop
zeros, rank, fit ZM):

- **b: corr(b_sim, b_actual) = 0.941, median |error| ≈ 4%** (median ratio 1.041).
  Exemplars: War and Peace 1.688 pred vs 1.681 actual; Ulysses 1.038 vs 1.045;
  Shakespeare 1.654 vs 1.724.
- **c: corr(log1p) = 0.822, median ratio 0.67** — correct scale (Arabian Nights
  6.6 vs 6.4; Iliad 119 vs 137; Shakespeare 315 vs 245). The 20x c-inflation seen
  with f4b's degenerate GMM is gone. Worst misses: Ulysses (33 vs 0.6), Critique
  (345 vs 58).
- Simulated V matches observed V per corpus (zero-truncation bookkeeping checks out).

**Why the naive version fails (kept as `naive_quantile_version.py`):** fitting ZM to
the *idealized latent* lognormal quantile curve gives b ≈ 8-12 with c pinned at V —
the latent curve plunges unboundedly in the deep tail. The observable curve is
censored at count 1: the sub-1 rate mass appears as the hapax pile instead of a
plunge. **The count floor is load-bearing for Zipf's observable shape**; any
sigma→b theory must go through the observation model.

**Interpretation:** b's cross-corpus variation (1.0-1.9) is not driven by sigma_T
(nearly constant ≈1.6) but by sampling depth and the head club — consistent with
the f6 ladder. Combined statement: type-population structure (two clubs, sigma_T ≈
1.6) + corpus size ⇒ the rank law's parameters, computed not fitted.

**Boundary:** prediction is at the fitted sampling depth. Cross-depth extrapolation
(predict full-corpus c from a slice) still fails with lognormal tails (f6 Part B)
→ f6b (heavier lower tail / dPlN) remains the open piece.

Reproduce: `scripts/run_f7.py` (reads f4d fits + f1 table; 2 reps/corpus, seeded).
