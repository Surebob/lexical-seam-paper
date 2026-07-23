# T3v2 Methodology Note

T3v2 replaces the fixed-exponent two-Pareto mixture from T3 with a c-generating two-regime rank-frequency construction.

Synthetic data:

`f(r) = A * (r + c_core)^(-b_head)` for `r <= K`

`f(r) = B * r^(-b_tail)` for `r > K`

`B` is chosen so the two pieces are continuous at `K`.

- V: `10000`
- b_head: `0.8`
- b_tail: `1.5`
- c_core: `10.0`
- K values: `[50, 100, 200, 500, 1000, 2000, 5000]`
- single-ZM fit: grid search over `c` plus affine least-squares solve for log-amplitude and slope at each `c`
- residual: `log(f_observed) - log(f_ZM_predicted)`
- SR coordinate: normalized log-rank `x = 0.05 + 0.95 log(r)/log(V)`
- SR grammar: live Section 2.4 grammar from `eml_zipf_enriched_search.py`: terminals `x, 1`; unary `neg, inv, sqr, sqrt, exp, log`; binary `eml, add, sub, mul, div, pow`
- depth: `2`; beam width: `50`; diversity weight: `0.35`
