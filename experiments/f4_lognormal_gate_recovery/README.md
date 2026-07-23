# f4 — gate recovery on synthetic population mixtures (2026-07-21)

**Question (ROADMAP C1, empirical half):** do two-lognormal type-population mixtures
generate the empirical erf-gate preference?

`scripts/run_f4.py`: 16 synthetic datasets (4 two-lognormal configs, 3 two-Pareto
configs, 1 single-lognormal control; 2 seeds each), Poisson-sampled counts, all five
gates fit with the canonical 9-parameter decoupled model (20 starts each).

## Result: inconclusive on the naive hypothesis; strong hint from the control

- Two-lognormal (arbitrary configs): erf wins 3/8 (algebraic 3, arctan 1, logistic 1).
- Two-Pareto: erf 3/6 — no discrimination either.
- Gate BIC spreads (52–562, median ≈ 200) are mostly below the empirical median 677.
- **Single lognormal: erf wins 2/2 with the largest margins in the whole run**
  (ΔBIC −342 and −181 vs best other gate).

## Interpretation

The naive "two populations ⇒ erf" story is not supported by uncalibrated synthetics.
The single-lognormal result points at a cleaner hypothesis: **erf preference tracks
log-normality of the type-frequency distribution** (the erf is exactly the Gaussian
CDF, i.e. the lognormal rank curve is an erf-quantile object), with or without
two-ness. The discriminating experiment must therefore be calibrated to the real
corpora rather than hand-picked: see `f4b_calibrated_twins`.

## Caveats
- Configs were hand-chosen; token counts and f_max for the Pareto family are
  unrealistic (PA3 f_max ≈ 38% of tokens).
- 20 starts (vs canonical 100); adequate for spread-scale conclusions, not for
  sub-BIC-unit distinctions.
