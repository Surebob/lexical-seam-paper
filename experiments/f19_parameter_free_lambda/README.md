# f19 — kill the parameter: frozen-λ (3-parameter) seam formula vs plain ZM (2026-07-23)

**Question (Greg's):** λ-ZM beats ZM, but it has one more parameter. Can the
parameter be eliminated — frozen at a universal value — and still win at
EXACTLY Zipf–Mandelbrot's parameter count?

**Setup.** f3's free fits showed (a) the exp generator g(x) = exp(x−1)−x wins
25/25, so the shape needs no per-corpus choice; (b) the fitted amplitude λ is
nearly universal (median 20.6, range [14.2, 24.4] across 25 corpora spanning
4k–29k vocabulary). `scripts/run_f19.py` therefore tests:

    ZM:        log f = a − b·log(r+c)               3 params/corpus
    λ*-ZM:     log f = a − b·log(r+c) + λ*·g(x)     3 params/corpus, λ* frozen

Anti-circularity protocols:
1. **LOO (English 25):** corpus i's λ* = median of the OTHER 24 corpora's
   free-fit amplitudes; (a, b, c) refit with the term frozen; full-data RMSE
   compared at equal parameter count.
2. **Held-out:** f3b protocol — binomial-thin 50/50 halves, fit both 3-param
   models on one half (λ still frozen), score fixed curves on the other, both
   folds.
3. **Cross-language transfer:** λ* = 20.61, the English median, applied
   UNCHANGED to the 7 canonical non-English corpora (Russian, Mandarin,
   Arabic, Latin, French, Spanish, Dutch) — those corpora contribute zero
   information to the constant.

## Result: the amplitude is a constant of the family, not a parameter

| panel | equal-param wins vs ZM | median impr | held-out folds | free-λ improvement retained |
|---|---:|---:|---:|---:|
| English 25 (LOO) | **25/25** | +8.2% | **45/50** | 95% |
| Non-English 7 (EN-frozen transfer) | **7/7** | +8.6% | **14/14** | 99% |

- Transfer highlights: Mandarin +19.4%, Arabic +14.0%, Russian +12.7% — with
  a constant fixed entirely by English books.
- Free-fit λ per language clusters at the same value: ru 22.8, zh 28.6,
  ar 23.2, la 16.0, fr 21.9, es 21.0, nl 21.1 (vs EN median 20.6).
- Worst cases stay positive: Aesop's Fables +0.9% (smallest corpus, V=4.4k);
  Latin +6.0% (lowest free λ, 16.0).

**Reading:** at Zipf–Mandelbrot's own per-corpus parameter count, the
two-population correction still beats it on 32/32 corpora and 59/64 held-out
fold-tests. The one "extra" number is a single global constant (λ* ≈ 20.6)
shared across languages and writing systems — i.e., the seam correction is
parameter-free relative to ZM, up to one universal constant of the family.

## Caveats
- λ* interacts with the fixed x-normalization (x = 0.05 + 0.95·log r/log V);
  the constant is defined relative to that convention.
- "Equal parameter count" means equal PER-CORPUS count; the global constant
  was estimated once, on English — the transfer panel is the test that this
  costs nothing.
- Rank-curve RMSE objective throughout (the paper's fitting protocol);
  BIC comparisons at equal p reduce to RMSE comparisons.

Outputs: `outputs/f19_results.csv`, `outputs/f19_summary.md`.
