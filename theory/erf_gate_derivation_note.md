# Why an erf gate? — derivation note (draft, 2026-07-21)

Status: working note (Claude/Fable + G.K.). Formalizes ROADMAP C1 in light of f4's
single-lognormal result. Not yet manuscript text.

## Setup

Let a corpus have V word types; type i has usage rate λ_i, observed count
~ Poisson(λ_i). The rank curve is f(r) = the r-th largest count. For large V the
deterministic skeleton is

    r(f) = V · S(f),        S = survival function of the type-rate distribution,

so the rank curve is the quantile function of the type-rate distribution:
f(r) = S⁻¹(r/V).

## Single lognormal population (Carroll's lognormal model)

If log λ ~ N(μ, σ²) then S(f) = Φ((μ − log f)/σ) and

    log f(r) = μ − σ · Φ⁻¹(r/V).                                   (T1)

The rank curve is a probit (Gaussian-quantile) object — all of its curvature is
erf-shaped by construction. This explains f4's control result: on single-lognormal
synthetic corpora the 9-parameter two-regime fit chooses the **erf** gate with the
largest margins in the run (ΔBIC −342, −181): the gate is being used to absorb
Gaussian-CDF curvature, and only the erf gate has exactly that shape.

## Two populations

Head population H (function words): N_H types, log λ ~ N(μ_H, σ_H²).
Tail population T (content words): N_T types, log λ ~ N(μ_T, σ_T²), μ_H > μ_T,
N_T ≫ N_H. Types above frequency f:

    R(f) = N_H·Φ((μ_H − log f)/σ_H) + N_T·Φ((μ_T − log f)/σ_T).

The head-population share among types ranked above f is

    π_H(f) = N_H·Φ_H(f) / R(f).

In the crossover window (where Φ_H transitions 0→1), total counts are dominated by
T, and the map u = log f ↦ log r = log R is smooth and monotone; to first order,
log r ≈ α − β·u with local slope β = −d log R/d u > 0. Substituting u ≈ (α − log r)/β
into Φ_H gives

    π_H(log r) ≈ Φ( (log k* − log r) / (β·σ_H) ),                  (T2)

i.e. the head share as a function of **log rank** is a Gaussian CDF —
π_H = [1 + erf((log k* − log r)/(√2·β·σ_H))]/2, exactly the model's erf gate with

    w_gate ≈ √2 · β · σ_H.                                          (T3)

Mean-field step: the fitted model blends the two branches' log-frequency
predictions with weight σ(log r); identifying that weight with the head share π_H
is an approximation (the true mixed rank curve interleaves types rather than
averaging predictions). This is the step to tighten for a manuscript-grade theorem.

## The discriminating prediction

The gate's *tail shape* inherits the head population's log-rate distribution tail:

- **lognormal head ⇒ Gaussian-tail gate (erf).**
- power-law (Pareto) head ⇒ exponential decay of Φ_H in u = log f ⇒
  **exponential-tail gate (logistic family)**.
- heavier/polynomial-tailed head ⇒ algebraic/arctan-like gates.

So the empirical result "erf wins BIC on 24/25 English corpora" (experiment 3e)
becomes a *measurement*: the function-word population's log-rate distribution has
Gaussian tails. This connects to the classical lognormal model of word frequencies
(Carroll 1967) and gives the erf preference a generative explanation instead of a
brute empirical one.

## Predictions & tests

- **P1** (quantitative): fitted w_gate ≈ √2·β·σ_H, with σ_H from a 2-component
  Gaussian mixture on the corpus's type log-frequencies and β the local log-rank /
  log-frequency slope at the crossover. → compute from f4b's GMM fits; correlate
  with canonical per-corpus w_gate. NOT YET RUN.
- **P2**: calibrated two-lognormal twins of real corpora prefer erf; calibrated
  two-Pareto twins do not. → **f4b (running)**.
- **P3**: single-lognormal synthetics prefer erf strongly. → **confirmed (f4)**.
- **P4**: corpora where the GMM head component is most cleanly Gaussian should show
  the largest erf-vs-logistic BIC margins. NOT YET RUN.

## Open items

- Tighten the mean-field step (T2) — either bound the interleaving error or derive
  the mixed rank curve directly as V·S_mix quantile and expand around the crossover.
- Dubliners (the single arctan win, smallest gate spread) predicts a head population
  with heavier-than-Gaussian log-rate tail — inspectable directly from its GMM fit.
- Poisson sampling adds a dispersion floor at small λ; check it does not distort the
  crossover-window argument (negative-binomial variant = B5 in the ledger).
