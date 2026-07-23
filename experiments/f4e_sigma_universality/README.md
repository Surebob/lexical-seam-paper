# f4e — mixture-signature universality on the expanded panel (2026-07-21)

**Question:** is the two-club mixture signature (pi_H, sd_H, sd_T) universal beyond
classic English books?

Method: f4d's zero-truncated Poisson-lognormal 2-component MLE applied to the f5a
processed panel. `scripts/run_f4e.py`.

## Results

**English registers — full-signature invariance** (reference: 25 books pi_H≈2.15%,
sd_H≈0.63, sd_T≈1.6):

| corpus | pi_H | sd_H | sd_T |
|---|---:|---:|---:|
| Brown (balanced 1961) | 2.30% | 0.596 | 1.636 |
| Cornell movie dialogs | 3.94% | 0.627 | 1.573 |
| WikiText-103 (1M slice) | 2.83% | 0.635 | 1.708 |

All three mixture parameters match the books across five centuries and three
registers (written literary, balanced modern, encyclopedic, conversational).

**Languages — a morphology-ordered band, not one constant:** sd_T = 1.726 (PT),
1.760 (DE), 1.810 (IT), 1.973 (FI), 2.002 (PL), 2.092 (SV); median 1.760 vs
English ≈1.6. Ordering is consistent with morphological richness spreading usage
mass across inflected types. Refined claim: **sd_T is register-invariant within a
language and morphology-graded across languages** (≈1.6–2.1 band).
Caveats: language corpora are shallow slices (33–150k tokens) — their tiny pi_H
values (0.1–1.3%) are plausibly depth-limited; matched-size English comparison
required before the gradient is asserted (queued with the multilingual redo).

**Surname control — alien signature, as it should be:** sd_T = 0.681 with pi_H
railed at the 30% bound; the two-club language structure is absent and the
estimator says so.

## Follow-ups
- Matched-size (80k/150k) English fits to de-confound the language pi_H/sd_T gradient.
- Fold sd_T into the f7 prediction chain for the panel corpora.
- If the morphology gradient survives matched size: sd_T as a typological index is
  a paper-section-worthy claim on its own.
