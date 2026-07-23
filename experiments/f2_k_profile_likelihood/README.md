# f2 — k profile likelihood, decoupled erf model (2026-07-21)

**Question (ROADMAP B1):** is the wide decoupled-erf k-scaling CI an identifiability
artifact — and does any tight scaling law survive?

Method: per English corpus, k fixed on a 13-point geomspace grid over [20, 1000]
(+ the canonical 2026-04-18 k̂), remaining 8 parameters optimized (6 starts: warm,
3 jittered, 2 random; trf, bounds identical to the canonical sweep). BIC profile over
k → per-corpus identified interval; profile-minimum estimates regressed on V.
350 fits total. `scripts/run_f2.py`; outputs in `outputs/`.

## Results

1. **k is identified within corpus.** Median ΔBIC≤2 interval width: **0.007 log10**
   (≈ ±1% in k); 0/25 corpora have an interval wider than 0.5 log10. The B1
   "within-corpus flat ridge" hypothesis is **refuted** — the profile minima are sharp.
   (Loosest: Principia Ethica, width 0.28; Critique of Pure Reason is near-flat only
   at the ΔBIC≤10 level, spanning 20→376.)

2. **No tight law for k itself.** log k_prof ~ log V: β = 0.758, CI [0.367, 1.149],
   R² = 0.386 — identical to the canonical estimate. The April v5.1 walk-back stands
   *for k*: cross-corpus variability in the gate centre is real, not noise.

3. **The tight law lives in the product.**
   **log s_prof ~ log V: β = 1.003, 95% CI [0.928, 1.078], R² = 0.967**, where
   s = k·w_tail is the model's tail-coordinate smoothing scale. Equivalently
   **s ≈ 1.2% of V** (per-corpus ratio range 0.009–0.013). The earlier coupled-model
   "k ∝ V^0.52" tightness was the single-width parameterization compressing this real
   linear law into the wrong parameter; the decoupled model exposed k's true
   variability and moved the lawful quantity into s.

## Interpretation

Two distinct scales, previously conflated:
- the **gate centre k** (where head weighting crosses 1/2): corpus-specific,
  weakly V-lawful — plausibly tracks corpus-level lexical composition (open question;
  connects to exponent-gap/hapax organizers from the fresh-eyes ledger);
- the **crossover width in linear rank, s**: a fixed ~1.2% fraction of vocabulary,
  R² = 0.97 across six centuries of English. This is a new headline-grade empirical
  regularity, stronger than the retracted k ~ √V claim.

Cross-check: the archive-based multiple regression (f1 companion analysis)
log k ≈ 0.98·log V + 0.23·log w_gate − 1.13·log w_tail (R² = 0.989) is consistent —
rearranged it says s·(w_tail)^0.13 ≈ V^0.98.

## Caveats
- 6 starts per profile point (vs 100 in the canonical sweep): warm-started, and the
  profile minima agree with canonical k̂ where comparable, but very flat secondary
  basins could be under-explored on small corpora.
- s was measured at the k-profile minimum, not profiled directly; a dedicated
  s-profile (fix s, optimize rest) is the natural follow-up before publishing the law.
- V^1.0 vs V^0.93–1.08: the CI excludes 0.5 and 0.545 decisively but a mild
  sublinearity is not excluded.
