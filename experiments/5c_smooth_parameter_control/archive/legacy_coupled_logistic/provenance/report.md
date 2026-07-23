# Smooth Parameter Sweep

Synthetic low-c smooth curves were generated from a continuity-normalized template using the median low-c English parameters, while varying:

- transition fraction `log(k)/log(V)`
- sigmoid width `w`
- tail slope `b2` (with head slope `b1` fixed at the low-c median)

- template vocab size: `10243`
- template medians: `a1=21.852757, b1=3.000000, c1=77.923884, c2=1000.000000`

## Winner Counts Across The Synthetic Grid

- full-RMSE winners: `{'is': 43, 'xpow': 58, 'exp': 115}`
- top-100 winners: `{'is': 89, 'exp': 88, 'xpow': 39}`

## Parameter Influence Rankings

- corr(theta, frac): `-0.539857`
- corr(theta, w): `0.305122`
- corr(theta, b2): `0.118858`
- standardized theta regression weights: `{'frac': -0.5398572685952798, 'w': 0.3051224792882102, 'b2': 0.11885807414466122}`

- corr(flip_lambda, frac): `0.392714`
- corr(flip_lambda, w): `0.554044`
- corr(flip_lambda, b2): `0.116274`
- standardized flip-lambda regression weights: `{'frac': 0.39271408483479714, 'w': 0.5540437051777384, 'b2': 0.1162736002471646}`

## Comparison To Empirical Low-c Corpora

- empirical low-c theta median: `105.902296`
- synthetic theta range: `[-145.359060, 172.837442]`
- empirical low-c flip-lambda median: `0.750000`
- synthetic flip-lambda range: `[0.000000, 1.050000]`

Interpretation: if one parameter dominates both theta and flip-lambda, it is the main low-c manifold control knob.

## Sample Rows

| frac | w | b2 | delta_b | zm_c | theta | winner full | winner top100 | flip lambda |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |
| 0.45 | 0.20 | 1.60 | 1.40 | 0.000 | -21.555 | `is` | `is` | 0.00 |
| 0.45 | 0.50 | 2.20 | 0.80 | 0.000 | -34.238 | `xpow` | `is` | 0.00 |
| 0.45 | 1.10 | 1.60 | 1.40 | 1.071 | -41.770 | `exp` | `xpow` | 0.40 |
| 0.45 | 1.40 | 2.20 | 0.80 | 53.037 | 75.134 | `exp` | `xpow` | 0.45 |
| 0.50 | 0.20 | 1.60 | 1.40 | 0.000 | -31.483 | `is` | `is` | 0.00 |
| 0.50 | 0.50 | 2.20 | 0.80 | 0.000 | -42.756 | `xpow` | `exp` | 0.00 |
| 0.50 | 1.10 | 1.60 | 1.40 | 23.436 | -18.799 | `exp` | `xpow` | 0.80 |
| 0.50 | 1.40 | 2.20 | 0.80 | 271.633 | -21.729 | `is` | `is` | 1.05 |
| 0.55 | 0.20 | 1.60 | 1.40 | 0.000 | -36.184 | `xpow` | `is` | 0.00 |
| 0.55 | 0.50 | 2.20 | 0.80 | 6.578 | -40.044 | `exp` | `exp` | 1.05 |
| 0.55 | 1.10 | 1.60 | 1.40 | 12.417 | -30.948 | `exp` | `xpow` | 0.95 |
| 0.55 | 1.40 | 2.20 | 0.80 | 157.585 | -23.485 | `is` | `is` | 1.05 |
| 0.60 | 0.20 | 1.60 | 1.40 | 0.000 | -40.250 | `xpow` | `is` | 0.00 |
| 0.60 | 0.50 | 2.20 | 0.80 | 40.397 | -33.548 | `exp` | `is` | 1.05 |
| 0.60 | 1.10 | 1.60 | 1.40 | 1.406 | -43.412 | `exp` | `exp` | 1.05 |
| 0.60 | 1.40 | 2.20 | 0.80 | 48.436 | -24.201 | `exp` | `is` | 1.05 |
| 0.65 | 0.20 | 1.60 | 1.40 | 0.000 | -47.753 | `xpow` | `is` | 0.00 |
| 0.65 | 0.50 | 2.20 | 0.80 | 19.546 | -41.721 | `exp` | `exp` | 1.05 |
| 0.65 | 1.10 | 1.60 | 1.40 | 0.473 | -46.303 | `xpow` | `exp` | 0.00 |
| 0.65 | 1.40 | 2.20 | 0.80 | 17.850 | -35.528 | `exp` | `exp` | 1.05 |
| 0.70 | 0.20 | 1.60 | 1.40 | 1.540 | -61.732 | `xpow` | `exp` | 0.00 |
| 0.70 | 0.50 | 2.20 | 0.80 | 11.339 | -67.824 | `exp` | `exp` | 1.05 |
| 0.70 | 1.10 | 1.60 | 1.40 | 0.518 | -49.626 | `xpow` | `exp` | 0.00 |
| 0.70 | 1.40 | 2.20 | 0.80 | 10.356 | -48.999 | `exp` | `exp` | 1.05 |
