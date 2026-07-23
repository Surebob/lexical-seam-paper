# f1 — Fresh reproduction + proper re-experiments (2026-07-20)

**Research question:** does the canonical foundation (ZM fits, step-2 Bregman winners)
reproduce under a fully independent reimplementation — and what changes under
amplitude-fair scoring and a joint 4-parameter λ-ZM fit?

Independent code (no legacy imports): `scripts/run_f1.py`. Protocol matches canon
(tokenizer, Gutenberg strip, tie-break, c-grid ZM fit, x-normalization).
Environment: repo `.venv` (numpy 2.5.1, scipy 1.18.0), Python 3.12.10.

## Outputs
- `outputs/f1_per_corpus.csv` — one row per English corpus (25)
- `outputs/f1_summary.md` — aggregates

## Findings

1. **Exact reproduction.** 25/25 token counts identical to canon; max |Δc| = 0.000000;
   max |ΔRMSE| < 1e-8; unit-amplitude winner family matches canon 25/25. A continuous
   multi-start refit never beats the canonical c-grid by >1e-4 (grid optimum confirmed).

2. **The unit-amplitude winner convention is load-bearing.** With a free amplitude per
   generator (1-param LSQ on the residual), the winner map changes on 14/25 corpora:
   is 17, xpow 7, euclid 1, **exp 0**. The "exponential-Bregman dominates low-c" claim
   is an artifact of forcing amplitude 1; at fair amplitude the low-c side scatters
   across the near-degenerate span (consistent with §3.9's manifold story, but stronger
   — winner identity is convention-dependent even at fixed grammar and metric).
   Fitted IS amplitudes cluster at 0.72–0.97 on all 11 high-c corpora (unit amplitude
   was accidentally near-optimal there — this *explains* why the parameter-free search
   helped on high-c and hurt on low-c, where fitted amplitudes are ≈0).

3. **λ-ZM 4-parameter model** (`log f = a − b·log(r+c) + λ·g(x)`, joint fit, c-grid +
   3-col LSQ): with g = exp(x−1)−x it beats plain ZM on **25/25 by BIC**, RMSE
   improvement min 2.6% / median 10.2% / max 23.0% — comparable to the terminal
   10-step enriched search (mean 11.0%) with a closed-form formula.
   **Algebraic identity:** in the canonical normalization,
   `exp(x−1) = e^(−0.95) · r^(0.95/ln V)` — the "exponential Bregman generator" is
   exactly a shallow power law in rank. The λ-ZM-exp model is therefore a
   ZM-plus-second-power-law additive family: the step-2 discovery and the two-regime
   mechanism are the same object algebraically. (In the joint fit, c inflates ~10×
   and λ ≈ 14–24: the optimizer reorganizes ZM into a slowly-varying background plus
   the r^γ regime term.)

4. **A2/A3 reproduce:** corr(c, IS-winner) = 0.883; corr(exp_gap, c) = −0.689;
   corr(hapax, winner) = −0.653; corr(TTR, c) = −0.706.

## Companion result (B1 precursor, archive data)

Reproducing the decoupled-erf k-scaling from `s2_v3_per_fit_results.csv`:
β = 0.758, R² = 0.386 (matches canon). But residual(log k | log V) correlates with
log w_tail at **−0.957**, and log k ~ log V + log w_gate + log w_tail gives
**R² = 0.989** with the V-exponent at 0.98. The model couples k and w_tail through
its own tail coordinate (`s = max(1, k·w_tail)`); the wide CI that "weakened" the
k-scaling claim is substantially a parameterization degeneracy, with identified
combination ≈ `k·w_tail^1.1 ∝ V^0.98`. A profile-likelihood rerun (fixed-k grid)
is queued to settle whether a tight k law survives reparameterization.

## Caveats
- λ-ZM comparison vs MOEZipf not yet run (queued; same extra-parameter class).
- Free-λ comparison uses OLS on the full residual (tail-dominated); head-weighted
  variants queued.
- BIC formula matches the paper's (p·log n + n·log MSE), inheriting its caveats.
