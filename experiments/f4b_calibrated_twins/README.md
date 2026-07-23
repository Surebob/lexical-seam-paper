# f4b — calibrated synthetic twins + P1 width check (2026-07-21)

**Question:** does each real corpus's own type-frequency structure reproduce the
erf-gate preference and the seam fingerprint (theory note predictions P1/P2)?

Method: for 12 representative corpora, fit a 2-component Gaussian mixture (EM) to the
type log-frequencies; simulate a **lognormal twin** (V rates from the mixture, Poisson
counts) and a **pareto twin** (two power-law populations, exponents from the f1
two-population fit); fit all 5 gates (16 starts) to each twin; record each twin's ZM c
and unit-amplitude step-2 winner. `scripts/run_f4b.py`, `scripts/p1_wgate_check.py`.

## Results

**Sufficiency: CONFIRMED.**
- Lognormal twins prefer **erf 12/12** (margins up to −729 Moby, −900 Ulysses).
- Twins inherit the seam fingerprint: **corr(twin c, real c) = 0.835**; step-2 winner
  matches the real corpus on 7/12.

**Quantitative width prediction (P1): SUPPORTED.**
- Derivation T3 predicts w_gate ≈ √2·β·σ_H (β = local |dlog r/dlog f| at the
  crossover; σ_H = GMM head-component width). Across 12 corpora:
  **corr(w_pred, w_fit) = 0.81** (log-corr 0.77), median ratio w_fit/w_pred = 1.65
  (systematic O(1) factor consistent with the mean-field approximation; to derive).
- The one strong outlier is **Dubliners** (ratio 0.32) — exactly the corpus the theory
  note pre-registered as deviant (sole empirical arctan winner, smallest gate spread).

**Specificity (P2): NOT ESTABLISHED.**
- Pareto twins ALSO prefer erf on 11/12 (exception: Shakespeare's pareto twin picks
  arctan with a huge 2806 BIC spread). Suspected confounds: fixed 130-type Pareto
  head (too small; transition sits at rank ~130 where Poisson noise smears the
  exponential tail toward Gaussian), rough level-matching, and possible
  zero-truncated-Poisson curvature at the hapax end that erf absorbs regardless of
  family.

## Follow-up (f4c, queued)
1. Noise-free deterministic rank curves per generative family (kills the Poisson
   confound).
2. Pareto twins with head size matched to the GMM head mass (hundreds of types).
3. P4: erf-vs-logistic margin vs head-component Gaussianity.
4. Derive the ~1.65 width factor from the interleaving correction.

## Status of the erf mechanism claim
"Two-lognormal structure ⇒ erf gate" is supported as *sufficient* and quantitatively
predictive (width, r=0.81), with c-fingerprint transfer; *uniqueness* against
power-law-population alternatives remains open pending f4c.
