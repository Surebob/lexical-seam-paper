# f18 — null-width calibration (seamless single-ZM truth)

Paired design, identical fitter/seeds: 25 real fits, 75 null fits (3 Poisson resamples per corpus from the fitted single-ZM curve).

## Real corpora (positive control)
- s/V: median 0.0121, IQR [0.0114, 0.0127]
- in language band [0.009, 0.015]: 23/25
- median |b_head - b_tail|: 0.601
- width-bound pins: 2/25

## Seamless nulls (the test)
- s/V: median 0.1515, IQR [0.0644, 0.4730], range [0.0138, 8.1732]
- in language band [0.009, 0.015]: 2/75
- median |b_head - b_tail|: 0.787
- width-bound pins: 27/75

**Reading: the language band is a property of the DATA, not the fitter — seamless nulls do not reproduce the width law.**