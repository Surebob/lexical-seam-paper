# f23 — the completist family bake-off (2026-07-24)

**Question:** referees can name curve families we hadn't formally raced:
lognormal-type rank curves (the classic Zipf rival), flexible polynomials,
power law with exponential cutoff, Yule–Simon. Same LSQ-on-log-frequency
objective, paper BIC formula, 25 English corpora. `scripts/run_f23.py`.

## Result: the two-regime families sweep every weight class

| family | params | median RMSE | lzm_free better on |
|---|---:|---:|---:|
| smooth two-regime erf (canonical, from f18 paired fits) | 9 | **0.1094** | — (it IS the seam model) |
| hard two-regime break | 5 | 0.1490 | 5/25 |
| **lambda-ZM free** | 4 | **0.1582** | — |
| cubic in log r | 4 | 0.1638 | 25/25 |
| power law + exp cutoff | 4 | 0.1761 | 24/25 |
| **lambda-ZM frozen (20.6)** | 3 | **0.1595** | 23/25 |
| lognormal-type (quadratic) | 3 | 0.1669 | 25/25 |
| single ZM | 3 | 0.1806 | 25/25 |
| Yule–Simon shape | 2 | 0.2197 | 25/25 |

Readings:
1. **At every parameter count, the seam family wins its weight class.**
   3p: frozen lambda-ZM beats lognormal and ZM. 4p: free lambda-ZM beats the
   cubic (a MORE flexible polynomial) and the cutoff family. 5p+: the only
   thing that beats lambda-ZM is... an explicit two-regime model. Top: the
   canonical smooth two-regime at 0.1094.
2. **The meta-lesson: every family that acknowledges two regimes beats every
   family that does not, at every complexity level.** The hard 5-param break
   winning BIC over the 4-param minimal correction is not a threat to the
   paper — it is the two-population structure winning again through a
   different parameterization (and the paper's smooth-beats-hard claim is
   about the full two-regime level, where smooth erf at 0.1094 crushes hard
   at 0.1490).
3. No smooth single-regime family survives: the completist sweep is closed.

Caveats: zm_cutoff inherits c from the ZM fit (small gift to that family —
it still loses 24/25); BIC at these V is RMSE-dominated (param penalty is
second-order), stated for transparency.
