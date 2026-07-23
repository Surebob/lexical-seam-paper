# T2 Coordinate-Alignment Verification

This addendum tests whether the T2 Poisson failure is primarily a coordinate-alignment issue. It compares a sorted-rank projection of the Poisson residual with a ZM residual expressed in a Mandelbrot-shift natural-coordinate proxy.

## Test 1: Poisson sorted-rank coordinate

The `(5, 20)` equal-weight Poisson mixture was fit by a single Poisson with lambda `12.5`. Counts were sorted by descending fitted single-Poisson probability and mapped to the manuscript normalized log-rank coordinate.

- winner: `div[pow[eml[x,1],eml[1,x]],eml[neg[x],mul[x,x]]]` ((pow(EML(x,1),EML(1,x))/EML((-x),(x*x))))
- residual RMSE: `23.6926`
- KL occurrence: `step 3, RMSE rank 6843, RMSE 23.6901`
- winner max |cosine vs ±KL|: `0.0754392`
- residual vs Poisson deviance cosine after sorting: `0.991815`
- residual vs KL on sorted-rank SR coordinate cosine: `0.0246038`
- top-20 Bregman-condition count: `0`

## Test 2: ZM residual in `c/(r+c)` coordinate

The Shakespeare single-ZM residual was paired with a linear-normalized `c/(r+c)` coordinate. This is a local coordinate induced by the Mandelbrot shift parameter; the output CSV also includes log-rank and log(r+c) reference coordinates.

- winner: `pow[mul[x,x],eml[1,x]]` (pow((x*x),EML(1,x)))
- residual RMSE: `0.18443`
- IS occurrence: `step 2, RMSE rank 385, RMSE 1.77941`
- winner max |cosine vs ±IS|: `0.000925182`
- residual vs IS on log-rank x cosine: `0.142046`
- residual vs IS on `c/(r+c)` x cosine: `0.0109909`
- residual vs IS on log(r+c) x cosine: `0.0196653`
- top-20 Bregman-condition count: `0`

### ZM coordinate controls

- log-rank control winner: `sub[sub[x,1],log[x]]`; IS occurrence: `step 2, RMSE rank 1, RMSE 0.182852`
- log(r+c) reference winner: `sub[sqrt[x],pow[x,x]]`; IS occurrence: `step 2, RMSE rank 29, RMSE 0.249814`

## Interpretation Matrix

Both predicted generators appeared under the tested coordinates.

Under this implementation, the Bregman generator recovery is coordinate-dependent rather than invariant. The Poisson deviance remains present in natural count/lambda space, but the sorted-rank projection does not make KL a competitive SR expression. Conversely, the Shakespeare ZM residual does not recover IS under the tested `c/(r+c)` coordinate, even though the manuscript log-rank coordinate is the one where IS is known to be recovered.
