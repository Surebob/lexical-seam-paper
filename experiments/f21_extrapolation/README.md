# f21 — extrapolation gauntlet: fit half the book, predict the whole curve (2026-07-24)

**Question:** lambda-ZM dominates at FITTING (42/42 BIC, 64/64 held-out,
32/32 equal-param). Does the advantage survive genuine FORECASTING —
chronological prefix half -> full curve, including tail ranks the fragment
never saw? `scripts/run_f21.py` (token-mass rescale; g's x frozen at
fit-time V; identical treatment for all models).

## Result: mixed — the edge survives but does not dominate

| band | free lambda-ZM beats ZM | frozen beats ZM | median improvement (free) |
|---|---:|---:|---:|
| full curve | 15/25 | 15/25 | +4.2% |
| head (top 100) | 14/25 | 11/25 | +40.1% |
| seen range | 15/25 | 15/25 | +4.9% |
| unseen tail | 14/25 | 15/25 | +2.8% |

**Honest verdict:** lambda-ZM is a better *descriptor at a given depth*, not
(yet) a depth *forecaster*. The reason is our own science: BOTH of the
model's depth-sensitive quantities drift between the half and the full
corpus — c collapses with depth (section 3.4) and lambda rises with depth
(f20: e.g. 16 -> 24 over the relevant range) — and naive extrapolation
freezes both at their half-corpus values. The improvement skew (median +40%
in the head where it wins, big losses where the drift bites) is exactly the
signature of un-modeled depth drift, not of a wrong correction shape.

**The obvious upgrade this begs (f21b, queued):** depth-AWARE extrapolation —
fit at depth T/2, then slide c and lambda along their measured depth curves
(sections 3.4 + f20) before predicting at depth T. The failure mode of naive
extrapolation is *predicted by the paper's own depth story*; correcting it
with the paper's own depth functions is the natural test of whether the
depth account is quantitatively right. If f21b wins big, the depth functions
graduate from descriptions to forecasting tools.

Scope note: the manuscript never claims formula-level extrapolation (the
section 3.5 out-of-sample claim belongs to the generative mixture, a
different object); nothing in the paper is touched by this result.
