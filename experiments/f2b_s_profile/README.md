# f2b — direct s-profile (2026-07-21)

**Question:** does the f2 candidate law (s = k·w_tail ∝ V, ≈1.2% of V) survive when s
is profiled directly as a first-class parameter?

Method: per English corpus, s fixed on an 11-point geomspace grid of s/V ∈
[0.002, 0.08] (+ the f2 s_at_min point); free params (a1,b1,c1,a2,b2,c2,k,w_gate)
with the tail scale pinned to the fixed s; erf gate; 6 starts (warm from f2).
300 fits. `scripts/run_f2b.py`.

## Result — the law holds under direct measurement

- **log s_min ~ log V: β = 1.0033, 95% CI [0.9278, 1.0789], R² = 0.9672** —
  numerically identical to f2's indirect estimate.
- **median s/V at profile minimum: 0.0118** (≈ 1.2% of vocabulary).
- Profiles are sharp: the ΔBIC≤2 interval collapses to a single grid point on the
  median corpus (grid spacing ≈ 0.16 log10), i.e. s is identified to better than
  grid resolution. With V ≈ 10⁴ observations per corpus, BIC moves fast in s.

## Caveats
- Grid resolution bounds the measured CI width from below; a fine local grid around
  each minimum would give true per-corpus CIs (cosmetic for the law, which uses the
  minima).
- The f2 s_at_min point was injected into each grid; neighbors are independently fit,
  so the sharpness conclusion is not circular, but the regression minima coincide
  with f2 by construction where that point wins. The bracketing profile shape is the
  new information here.
- 6 starts per point (warm-started); very small corpora may under-explore.

Combined statement for the manuscript (f2 + f2b): *the decoupled-erf model's tail
crossover scale s is sharply identified per corpus and scales linearly with
vocabulary size, s ≈ 0.012·V (R² ≈ 0.97, exponent CI [0.93, 1.08]), while the gate
centre k is identified but not V-lawful.*
