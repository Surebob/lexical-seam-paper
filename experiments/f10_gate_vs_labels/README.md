# f10 — gate vs labels (2026-07-21)

**Question:** does the erf gate fitted to the rank curve (blind to identity) coincide
with the labeled function/content mixing curve from the EXP03 annotations?
10 label-complete corpora; instantaneous closed-class share (geometric window
r/1.4–r·1.4) vs the canonical gate σ(log r). `scripts/run_f10.py`
(v1 with cumulative-share bug superseded same day).

## Verdict: the STRONG identification fails — and is replaced by a finding

- corr(log k_gate, log k_lab) = 0.20; **median k_gate/k_lab = 3.10** — the
  statistical gate sits ~3× deeper in rank than the labeled grammar crossover.
- corr(w_gate, w_lab) = −0.43; median curve RMSE 0.23 (share units).
- Labels are reliable: 3-tagger pilot agreement 93–97% (closed/open), ~90–92%
  exact tag — the mismatch is real, not annotation noise.

**Reinterpretation (for v6):** the statistical head club ≠ the closed class. It
contains the grammatical core PLUS the most common content words — high-rate,
narrow-dispersion types that behave like function words without being them. The
seam is a **usage-dynamics boundary**, drawn ~3× deeper than the grammatical one.
Consistent with: f4d head sizes (20–820 > closed-class ~300 in large corpora),
f11-A (amplitude gap alone creates a seam — no grammatical distinction needed),
and the ablations (which manipulated the statistical club and still stand).
v6 drops any k_stat = k_POS identification; the two scales are related but
distinct objects (the paper's §3.6 already half-said this; now it's measured).

## Additional results
- **Labeled k_POS ~ V^0.298, CI [0.08, 0.52]** (n=10, mostly smaller corpora) —
  undercuts the spaCy-based 0.545 on this subsample; v6 must caveat the POS
  scaling claim and prefer trusted labels where available.
- Labeled ablation: unit step-2 helps 3/10 full → 4/10 content-only (no clean
  kill on this low-c-heavy subsample at unit amplitude; consistent with the
  amplitude findings — the removal test's proper form uses fitted amplitude).
- Inter-annotator pilot statistics preserved in outputs (reviewer armor).

## Status of the "flagship overlay figure" idea
Retired in its strong form. The v6 mechanism chapter instead gains a sharper
section: "the seam is a usage boundary, not a grammatical one — measured against
17k human-quality labels."
