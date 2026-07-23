# f4d — Poisson-lognormal mixture fits + hunch re-tests (2026-07-21)

**Trigger:** f4b's GMM was degenerate — one component sat on the hapax spike
(1−π_H matched the hapax ratio to 3 decimals). This experiment fits the latent
type-rate distribution properly: 2-component lognormal mixture observed through
Poisson sampling with zero-truncation, MLE via Gauss-Hermite quadrature, head
constrained to the minority high-rate component. All 25 English corpora.
`scripts/run_f4d.py`.

## Headline result (new, for v6)

**A two-component latent structure is demanded by BIC on 25/25 corpora** — direct
distributional evidence for the two-population hypothesis, independent of rank
curves, symbolic search, or ablations:
- head component: median **π_H = 2.15%** of types (range 0.28–5.93%), N_H ≈ 20–820,
  **narrow** (σ_H median ≈ 0.63);
- tail component: broad, with **σ_T ≈ 1.6 nearly constant across all 25 corpora**
  (range 1.50–2.06) — a candidate new universal (six centuries, one dispersion).
- Character check: the small-N_H outliers are the empirically weird corpora
  (Dubliners N_H = 20 — the lone arctan corpus; Poe 38; Sherlock 46).

## Hunch outcomes (both resolved negative — recorded to prevent re-derivation)

1. **P1 with the true head width: FAILS** — corr(w_pred, w_fit) = 0.02. The earlier
   r = 0.81 was carried by the degenerate GMM's σ, which measured the **bulk**
   (non-hapax) dispersion. Empirical law restated: w_gate tracks √2·β·σ_bulk
   (r = 0.81), NOT the head-component width. The theory note's T2/T3 derivation must
   be rebuilt around the bulk/tail object (the gate width reflects how fast the
   rank curve's local slope transitions, governed by the broad component).
2. **√e factor: DEAD** — the 1.6491 ≈ √e match was an artifact of the degenerate σ;
   with corrected components the median ratio is 3.47 and no clean constant appears.
3. **s = N_H identity: DEAD** — corr 0.61, slope 0.36, s/N_H spans 0.19–3.3. The
   s-law remains "constant ≈1.2% quantile band of V" without a head-size reading.

## Follow-ups
- Rebuild the width derivation with the bulk dispersion as the driving scale
  (define σ_bulk operationally: sd of log-frequency over non-hapax types; verify
  r ≈ 0.81 reproduces under that exact definition).
- Regenerate calibrated twins from the PLN mixture (proper two-population twins)
  and re-run the f4b/f4c battery.
- Investigate σ_T ≈ 1.6 constancy (connection to the Zipf exponent b ≈ 1?) — if
  derivable, it links the lognormal bulk to Zipf's law itself.
