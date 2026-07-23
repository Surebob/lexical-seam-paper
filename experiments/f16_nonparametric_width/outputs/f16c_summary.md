# f16c — planted-truth calibration + FWHM instrument

## A. erf-fitter recovery on planted seams (n=72)
- log s_hat ~ log s*: slope 0.472, R2 0.085; median multiplicative bias s_hat/s* = 0.405
  - erf-planted: median s_hat/s* = 0.374 (n=36)
  - logistic-planted: median s_hat/s* = 0.428 (n=36)
- FWHM instrument vs planted s*: corr(log, log) = -0.543

## B. FWHM on real corpora vs single-regime twins
- English: median FWHM/V = 0.0043, median amplitude = 0.293 (n=25)
- twins:   median FWHM/V = 0.0063, median amplitude = 0.185 (n=25)
- amplitude separation (language/twin): x1.6