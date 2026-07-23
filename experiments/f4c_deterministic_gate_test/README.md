# f4c — deterministic gate discrimination (2026-07-21)

**Question:** does the erf preference survive on noise-free curves, and does it
discriminate lognormal from power-law population structure? (Resolves the f4b
Poisson-confound hypothesis.)

Method: for 12 corpora, exact quantile rank curves from the corpus-calibrated
two-lognormal mixture and a two-Pareto mixture (head sized to the GMM head mass),
each with Poisson noise on/off; 5 gates × 16 starts. 240 fits. `scripts/run_f4c.py`.

## Results

Winner counts (independent gates):
| family | noise | erf | logistic | algebraic | arctan |
|---|---|---:|---:|---:|---:|
| lognormal_mix | off | 8/12 | 1 | 0 | 3 |
| lognormal_mix | on  | 11/12 | 1 | 0 | 0 |
| pareto_mix | off | 8/12 | 2 | 2 | 0 |
| pareto_mix | on  | 7/12 | 2 | 1 | 2 |

1. **The Poisson-confound hypothesis is refuted**: noise-free power-law curves still
   mostly prefer erf. The erf gate is NOT a unique signature of lognormal populations.
   Reading of the deviants (logistic wins on W&P/Quixote Pareto; algebraic on
   Moby/Ulysses Pareto with huge margins): the gate preference tracks how fast the
   head population *terminates* — Gaussian tails AND hard Pareto cutoffs both
   terminate fast and look erf-like; only genuinely exponential-tailed boundaries
   favor logistic.

2. **The family discrimination lives in the c-fingerprint, not the gate** (from
   f4b twins CSV):
   - lognormal twins: corr(twin c, real c) = **0.835**, step-2 winner match 7/12;
   - pareto twins: corr = **0.319**, match 4/12, and twin c collapses to ≈0–6
     regardless of parent — power-law populations cannot reproduce the empirical
     c structure at all.
   - Caveat: lognormal twins preserve c *ordering* but inflate its absolute scale
     ~20× (twin c 646–5897 vs real 1–245). Open calibration item (Poisson floor +
     GMM tail mismatch suspected).

## Final mechanism verdict (for manuscript v6)

The two-lognormal population model is supported as the generative mechanism by the
**conjunction** of: (a) sufficiency — calibrated lognormal twins reproduce the erf
preference (11–12/12 with noise); (b) fingerprint uniqueness — lognormal twins
transfer the parent's c-ordering and step-2 winners while power-law twins cannot;
(c) quantitative width prediction — w_gate ≈ √2·β·σ_H at r = 0.81 with Dubliners the
pre-registered outlier. The erf gate alone is NOT claimed as a lognormal signature
(f4c kills that shortcut); v6 §mechanism must present the conjunction, not the gate.

Follow-ups: fix the twin c-scale inflation; Laplace-log-rate control (exponential
boundary ⇒ logistic prediction); derive the 1.65 width factor.
