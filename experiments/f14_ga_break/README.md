# f14 — width behavior of rival generative models (the "break" that reversed)

**Original goal.** Show the Gerlach-Altmann (PRX 2013) core/non-core model
cannot produce the width law s ≈ 0.012·V — its fixed core was expected to pin
the crossover at fixed rank, so s/V should fall as V grows.

**What happened instead (f14, `scripts/run_f14.py`).** Implementing their
model from their Eqs. (3)-(5) with their published English constants
(Nc_max = b ≈ 7,900, γ = 1.77 → α = 0.77, p_c⁰ = 0.99) and fitting the
canonical erf-gate model to prefix snapshots (M = 0.2-3.2M tokens,
V = 16k-88k, 2 runs): **s/V = 0.0112-0.0125, slope of log(s/V) vs log V
= +0.06 — the width law, at our constant, from their model.** The prediction
behind the planned break was wrong: the fitted seam on GA curves is not the
core boundary and scales with V.

**The nulls (f14b, `scripts/run_f14b_nulls.py`).** What pins 0.012?

| process | s/V (median) |
|---|---:|
| GA reference (two-class, fixed core) | 0.0122 |
| pure Simon + decaying innovation, NO core | 0.0118 |
| real English, 25 corpora (f2) | 0.0118 |
| classic Simon (constant innovation; single-regime) | 0.0101 |
| US surnames (non-linguistic control, f5b) | 0.0266 |

1. **The core is irrelevant to the width** (0.0122 vs 0.0118): removing GA's
   distinguishing structure changes nothing. Consistent with f10 (the seam is
   not a membership boundary) and L3c (the frontier is not class membership).
2. **The width is partly generic**: even a genuinely single-regime process
   yields a stable fitted width fraction (0.0101) — the erf fit finds the
   finite-sample tail handover on any Zipfian curve. But the constants are
   **family-specific**: constant-innovation 0.0101, decaying-innovation
   ~0.0118-0.0122, surnames 0.0266. The width constant is a fingerprint of
   the growth family, not an operator artifact (an artifact would give one
   number everywhere).
3. **Natural language sits exactly on the decaying-innovation value** —
   0.0118 to three decimals. New, sharp datapoint: whatever generates real
   lexicons behaves, at the width level, like preferential attachment with
   decaying innovation (which LNRE/Heaps phenomenology independently
   supports).

**Consequences.**
- The width law's status STRENGTHENS as a law (robust invariant reproduced by
  every language-like generative account tested: PLN mixture at physical
  basin, GA, pure decaying-Simon) and WEAKENS as a mechanism discriminator
  (it cannot tell two-class from one-class growth).
- Manuscript corrections applied (§3.3 object paragraph, §4.3 GA block):
  drop "fixed core predicts fixed crossover rank" (empirically false for the
  fitted object), replace with the shared-invariant framing + the family-
  fingerprint table. What discriminates accounts is the boundary's *identity*
  (§3.6 usage vs grammar) and *depth behavior* (§3.4), not the width.
- `docs/PRIOR_ART_SWEEP.md` gains a dated addendum correcting its s-law row.
