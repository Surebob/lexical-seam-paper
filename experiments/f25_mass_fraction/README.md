# f25 — the club's token-mass share (2026-07-24)

**Question:** k (seam centre) is lawless in rank space. Is the seam's
location lawful in MASS space — does a constant fraction of all tokens sit
above it? `scripts/run_f25.py`; k from the f18 canonical real fits
(k = s/w_tail; the 2 pinned escapes excluded).

## Result: a soft third invariant — the club carries ~ two-thirds of usage

- M(k) across 23 corpora: **median 0.661**, IQR [0.608, 0.704], range
  [0.458, 0.796]; relative spread 0.15.
- Versus the rank coordinate k/V: median 0.0204, range [0.0072, 0.0520]
  (7x variation); relative spread 0.27.
- corr(M(k), log V) = +0.30 — mild size dependence.

**Reading (calibrated):** the seam centre is markedly more stable in
cumulative-mass coordinates (~2/3 of tokens above it, 1.7x range) than in
rank coordinates (7x range). Part of that compression is mechanical — the
cumulative mass curve is concave, so it compresses rank scatter — but not
all of it: a 7x rank wander maps to only 1.7x in mass. "Roughly two-thirds
of all usage is club usage" is a fair rough invariant; it is NOT a
precision constant at the level of s/V (±10%) or lambda (±25%), and k's own
estimation noise propagates into M(k).

Status: Paper-2 lead (mass-coordinate formulation of the seam; panel +
multilingual extension; does M(k) ride a depth curve like everything else?).
Not manuscript-facing for Paper 1 — the paper makes no location-lawfulness
claims.

Outputs: `outputs/f25_mass.csv`, `f25_summary.md`.
