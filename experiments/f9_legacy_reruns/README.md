# f9 — legacy reruns (2026-07-21)

Three blocked/unsupported legacy items addressed. `scripts/run_f9.py`.

## A. 2b same-exponent control (α=1.5/1.5) — STILL OPEN, with a new insight
My construction (two Pareto(1.5) rate populations, Poisson-sampled) produced a
HELPFUL IS correction on 3/3 reps with poor ZM fit (rmse ≈0.26) — contradicting the
manuscript's expectation for this control. Likely indicts the protocol, not the
claim: Poisson observation + the count-1 floor curves synthetic tails (f7), so a
sampled same-exponent mixture is NOT the clean single power law the deterministic
construction gives. **New confound finding: observation noise alone can generate
seam-like corrections on synthetic power laws.** Action: rerun 2b with the exact
legacy protocol (`src/legacy_engines/zipf_synthetic_mixture.py`) and re-examine
whether the April positive mixture results are robust to the observation model.

## B. E2 (WLS winner stability) — v5.1 CLAIM REFUTED
Fresh 4-generator comparison under three fit objectives: winners stable on only
**14/25**. Under rank-weighted and frequency-weighted objectives, **all 11 high-c
IS corpora flip to exp** — systematically, not noisily. Third independent proof
(after free-amplitude f1 and head-window f3) that winner identity is a property of
the scoring functional; the seam is the invariant. v6 §3.1 robustness paragraph
must be rewritten accordingly (the audit's E2 suspicion was correct).

## C. 10a transfers — UNBLOCKED, claim VINDICATED
- War&Peace → Bible: transfer RMSE 0.1498 vs Bible in-domain 0.1494 vs ZM 0.1882
  (essentially perfect transfer).
- Moby → Bible: 0.1799, still beats ZM.
- Shakespeare → W&P: 0.1513 vs in-domain 0.1458.
The universal high-c head shape (ledger A4) extends to the anthology corpus;
10a's BLOCKED.md can be retired with these rows.
