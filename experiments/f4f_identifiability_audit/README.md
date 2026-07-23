# f4f — identifiability audit for the latent mixture (2026-07-21)

**Question:** is the fitted tail spread σ_T (the "morphology dial", f4d/f4e/f8) a
well-defined measurement? Triggered by f6b's discovery that different parameter
basins fit one depth equally but predict differently. `scripts/run_f4f.py`.

## Verdict: the dial is DEAD as a latent-parameter claim

**A. Model-free check — the ordering REVERSES** (corr = −0.53 between
sd(log count | count ≥ 5) and the f8 fitted σ_T). The raw-count spread at matched
depth mostly measures tokens-per-type; morphology-rich languages have MORE types
per token → compressed observable counts → the model inflates latent σ_T to
compensate. The "dial" tracked sampling depth per type, not a latent usage spread.

**B. Likelihood landscape — genuinely unidentified from one depth.** Near-optimal
solutions (ΔNLL ≤ 5) span σ_T 1.9–4.0 and n_total 20k–600k on the same corpus.
The max-likelihood solution is usually NOT the best cross-size predictor (best
predictors have broader tails, σ_T ≈ 3.2–3.7) — retroactively explaining both
f4d's narrow-basin fits and f6b's success with broad basins.

**C. Joint multi-depth fitting — partial rescue only.** Some corpora pin σ_T
tightly (Shakespeare spread 0.04, P&P 0.05), others stay loose (Moby 1.9,
Ulysses 1.6); no clean universal or gradient emerges. (6 starts/corpus, 65k base —
indicative, not exhaustive.)

## Retired claims (moved to graveyard)
- "σ_T ≈ 1.6 universal for English" (f4d) — basin-conditioned.
- The cross-language σ_T gradient / morphology dial (f4e, f8) — depth-per-type
  artifact under a fixed protocol.
- Ulysses/Shakespeare "richer than German" ornaments — same artifact.

## What stands, and why
- All observable-level results are untouched: s-law, λ-ZM, the sampling-depth
  ladder, c-collapse, f6b's out-of-sample V/b/c predictions (prediction operates
  on observables and implicitly selects predictive basins — Part B confirms).
- f8's matched-size c results and the W&P translation pair (EN c=10.1 vs RU c=0.0)
  stand — ZM parameters are observables.
- f7's "histogram regenerates the curve" stands as a family-level statement (many
  basins regenerate the observables); it was never a parameter-value claim.

## New audit spawned (queued)
- **Re-test "two components demanded 25/25" (f4d) under broad-basin multi-start
  fitting** — the 1-vs-2-component BIC comparison used the same restricted
  protocol; v6 must not lean on it until re-verified. Also applies to the
  register-invariance absolute numbers (relative same-protocol comparison stands;
  absolute values do not).
