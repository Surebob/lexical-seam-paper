# f18 — null-width calibration: the instrument on seamless data (2026-07-23)

**Question (Greg's, in its rawest form):** "maybe what we're seeing is just a
product of Zipf's itself" — i.e., is s/V ≈ 0.0118 simply what the canonical
9-parameter erf fitter returns on ANY Zipfian curve, seam or no seam?

**Design (the dual of f16d).** f16d planted the fitted two-regime truth and
showed the estimator recovers it (~1.03×). f18 plants the SEAMLESS null: for
each of the 25 English corpora, fit a single Zipf–Mandelbrot curve (a, b, c —
no gate, no second population), generate its exact expected curve at the
corpus's own token count, Poisson-sample it (3 resamples), and push the
samples through the IDENTICAL canonical 24-start fitter, alongside paired
fits of the real corpora with the same seeds. `scripts/run_f18.py`,
100 fits total, incremental CSV + resume.

## Result: the language band belongs to the data, not the operator

| | real corpora (25) | seamless nulls (75) |
|---|---|---|
| s/V median | **0.0121** | **0.1515** (12× the language value) |
| s/V IQR | [0.0114, 0.0127] | [0.0644, 0.4730] |
| s/V range | — | [0.0138, 8.17] (~3 decades) |
| in language band [0.009, 0.015] | **23/25** | **2/75** |
| width-bound pins (degeneracy flag) | 2/25 | 27/75 |

- On seamless data the gate parameters are unidentified; the fitter chases
  noise — widths scatter across three orders of magnitude and pin at bounds
  27/75 times (the detectable-degeneracy signature).
- The 2/25 real out-of-band fits (Aesop, Dubliners) are the same detectable
  failure mode: both width bounds pinned, s/V > 5 — a bad random-start draw,
  screened in the empirical pipeline by profile-likelihood cross-checks and
  multi-start search (both corpora sit at ≈0.012 in the canonical panel).
- The slope-collapse diagnostic (|b_head − b_tail| shrinking on nulls) did
  NOT discriminate (null median 0.787 vs real 0.601): degenerate null fits
  chase noise with wild slopes rather than collapsing to one. The
  discriminators are band occupancy, dispersion, and pin rate. Reported
  honestly as a wrong secondary prediction.

**Reading:** if the width law were an operator artifact, seamless Zipfian
data would reproduce it. It does not — by an order of magnitude in median
and by 23/25-vs-2/75 in band occupancy. Together with f16d (planted seams
recovered at ~1.03×), the instrument is now calibrated from BOTH sides:
it finds seams that exist, and does not find seams that don't.

Outputs: `outputs/f18_null.csv`, `outputs/f18_summary.md`.
