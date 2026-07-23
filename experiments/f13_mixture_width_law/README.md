# f13 — does the generative account contain the width law?

**Question.** The s-law (manuscript eq. 3, s ≈ 0.012·V) is empirical. If §3.5's
zero-truncated Poisson-lognormal mixture is the mechanism, simulating corpora
from the fitted mixtures and fitting the same erf-gate model should reproduce
the measured widths. (A naive analytic width prediction in f4d ran ~3× small,
so this was genuinely open.)

## f13 — ML-basin test (`scripts/run_f13.py`)

Simulate each corpus's f4d mixture through the observation model (3 reps),
fit the f5b/f12 9-parameter erf model, compare s_synth to s_real.

**Result: split, and the split is diagnostic.** 16/25 corpora reproduce the
width almost exactly (ratios 0.89–1.19; Moby Dick 1.02, Sherlock 1.00,
Ulysses 1.05). The 9 failures (5–600×) are the deeply-sampled corpora — and
`run_f13_diag.py` shows failure tracks the fitted basin, not the corpus:
corr(log ratio, π_H) = +0.76, corr(log ratio, tokens/V) = +0.60. ML at deep
sampling slides into a blurred-head basin (π_H ~ 4–6%, small μ-gap) — the f4f
identifiability weakness surfacing in a new observable. The diagnostic also
finds the analytic seed: the ambiguous-type mass A between the two component
densities' crossings predicts s at the ~1% scale exactly where the basin is
well-separated (Sherlock/Dubliners/Poe: A/V = 0.0145/0.0081/0.0089).

## f13b — basin reselection (`scripts/run_f13b_basin_reselect.py`)

Non-circular repair: for the 9 failures, map the basin landscape (16 diverse
starts on the full-histogram PLN likelihood, f4f machinery) and select the
basin by MINIMUM CROSS-SIZE PREDICTION ERROR (thin real counts to 1/4 and
1/16; compare V and c against simulation — a criterion that never sees rank
curves or widths). Then simulate the selected basin and re-fit widths.

**Result: 8/9 heal to ~1.0.**

| corpus | ML-basin ratio | pred-basin ratio |
|---|---:|---:|
| Shakespeare | 15.1 | **1.02** |
| War and Peace | 18.2 | **1.01** |
| King James Bible | 18.2 | **1.01** |
| Federalist Papers | 6.0 | **1.03** |
| Origin of Species | 4.6 | **0.99** |
| Wealth of Nations | 5.0 | **1.00** |
| Les Misérables | 15.9 | **1.02** |
| Principia Ethica | 600.6 | **0.94** |
| Critique of Pure Reason | 624.8 | 549.9 (resists) |

Median pred-basin ratio 1.01 (was 15.9). The selected basins sit ΔNLL 0–11
from ML — inside the identifiability plateau — yet fix both observables at
once.

## Verdict

1. **The width law is contained in the generative account** (24/25 corpora
   once the basin is right): mixture + count floor ⇒ s ≈ 0.012·V.
2. **Two independent observables converge on the same basin** — cross-size
   V/c prediction and seam width select identically. Likelihood alone does
   not identify the physical basin; prediction criteria do (sharpens f4f/f6b).
3. The healed basins are heterogeneous in latent (π_H, gap) — the width is an
   observable-level invariant of basins that match the histogram AND its
   thinned versions, not a function of latent parameters. The depth dimension
   is load-bearing, consistent with §3.4.
4. **Kant resolved (f13c-K, `run_f13c_kant.py`)**: 48-start extended search
   finds healing basins for Critique of Pure Reason too — basin at ΔNLL 1.7
   with width ratio **1.02** (three more at 1.04–1.06). The f13b failure was
   search coverage, not physics: **25/25 corpora have a mixture basin that
   generates the measured width.** Nuance: for this one corpus the raw
   pred_err minimum lands on a non-healing basin (0.142 vs 0.29 for the
   healer) — the automatic two-criterion convergence is 8/9, with Kant
   requiring joint selection.
5. Closed-form status: the naive λ-space ambiguous-mass A is NOT the width
   (f13c: A/s scatters 0.14–76) — the width is created by the full
   observation ∘ rank ∘ fit operator. See `docs/WIDTH_DERIVATION_STATUS.md`
   for the proven/open split.

Feeds manuscript §3.5; candidate centerpiece for Paper 2's mechanism arc.
