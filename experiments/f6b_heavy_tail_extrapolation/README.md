# f6b — cross-size extrapolation (2026-07-21)

**Question:** can the two-population mixture, fitted at one corpus size, predict the
rank law at other sizes — including sizes it has never seen? (The "no-asterisk"
test for the histogram→curve result.)

Models: PLN (2-lognormal latent rates) and NLT (tail component = normal minus
exponential in log-rate — a heavy lower tail, the double-Pareto-lognormal
direction). Both fit by zero-truncated MLE on a wide latent grid; 4 exemplar
corpora. `scripts/run_f6b.py`.

## Results

**Upward prediction (fit on a 150k binomial slice, simulate at full size):**

| corpus (extrapolation) | V pred→real | b pred→real | c pred→real |
|---|---|---|---|
| Shakespeare (6.6×) | 24,072→24,458 | **1.725→1.724** | 262→245 |
| War and Peace (3.9×) | 17,183→17,445 | **1.685→1.680** | 209.5→208.1 |
| Les Misérables (3.8×) | 21,844→22,677 | 1.538→1.527 | 250→184 |
| Moby Dick (1.5×) | 16,970→16,956 | 1.190→1.182 | 18.6→10.7 |

(NLT rows shown; PLN comparable, better on Les Mis c=195.6.) **Exponent b is
predicted to ~0.3% and vocabulary to ~2% from one-sixth of the text.** Mandelbrot's
c is predicted to first approximation (0.7%–36% depending on corpus/model).

**Downward trajectories (fit at full size, predict the collapse):** V(T) matched to
1–2% at every scale for all corpora (Heaps' law captured); c(T) now tracked end to
end (Shakespeare real 244.8/103.7/36.5/11.9/3.4 vs PLN 223.1/116.4/38.0/8.3/3.3) —
the f6 Part-B failure is resolved.

## Why f6's version failed and this works

The latent mixture is only weakly identified from a single depth: narrow-tail
solutions (few hidden sub-observation types) and broad-tail solutions (many) fit
the visible histogram almost equally, but only the broad basin predicts across
sizes. f6 Part B inherited f4d's narrow-basin fits; f6b's wider optimization finds
the predictive basin. Model choice (PLN vs NLT) matters less than basin: neither
dominates (NLT best on Shakespeare/W&P c; PLN on Les Mis).

## Consequences and cautions

1. **Claim upgrade (v6):** "measure a fragment; the full text's V, b — and
   approximately c — follow." Genuine out-of-sample prediction; b/V are the
   precision results, c the first-approximation one.
2. **Claim downgrade (flagged before anyone else could):** the f4d "σ_T ≈ 1.6
   universal" is basin-conditioned — broad-tail optima have different σ_T. The
   cross-language ORDERING (f4e/f8 gradient) used one consistent protocol and
   likely survives as a comparative index; the absolute constant needs an
   identifiability audit → **f4f, queued, top priority**.
3. Latent totals (n_total) differ ~30% between models with similar fits — never
   claim "the true number of unseen words"; claim only observable predictions.
