# T3 Methodology Note

Synthetic data are rank-frequency curves generated from equal-support mixtures of two normalized discrete Pareto components over ranks `1..10000`.

- exponents: alpha1 `1.5`, alpha2 `1.3`
- sweep: `w1 = 0.05, 0.10, ..., 0.95`
- single-ZM fit: nonlinear least squares on log-frequencies, implemented as grid search over `c` with affine least-squares solve for intercept and slope at each `c`
- residual: `log(f_observed) - log(f_ZM_predicted)`
- SR coordinate: normalized log-rank `x = 0.05 + 0.95 log(r)/log(V)`
- grammar: live Section 2.4 grammar from `eml_zipf_enriched_search.py`: terminals `x, 1`; unary `neg, inv, sqr, sqrt, exp, log`; binary `eml, add, sub, mul, div, pow`
- depth: `2`; beam width: `50`; diversity weight: `0.35`; keep-all-until-step: `2`
