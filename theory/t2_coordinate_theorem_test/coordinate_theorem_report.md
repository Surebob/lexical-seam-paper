# T2 Coordinate-Theorem Prediction Test

This primary test expresses the Poisson two-component mixture residual in affine-normalized count coordinate, `x = 0.05 + 0.95*k/k_max`, rather than log-rank or sorted-rank coordinate.

The `(5, 20)` equal-weight mixture is fit by a single Poisson with lambda `12.5`. The support `k=0..30` captures `0.993263` of mixture probability mass.

- residual RMSE: `3.26063`
- shifted KL occurrence `x log(x)-x+1`: `not found`
- unshifted negative-entropy occurrence `x log(x)-x`: `step 3, RMSE rank 7093, RMSE 3.71143`
- either predicted form in top 5 at depth 3: `False`
- either predicted form in top 20 at depth 3: `False`
- winner: `mul[pow[add[x,x],sub[x,1]],pow[add[x,x],add[x,x]]]` ((pow((x+x),(x-1))*pow((x+x),(x+x))))
- winner RMSE: `1.82744`
- winner max |cosine vs ± shifted KL|: `0.706705`
- best cosine candidate vs shifted KL: rank `3727`, cosine `0.999408`, RMSE `2.94114`
- winner satisfies Bregman boundary conditions: `False`

## Top 5 at depth 3

- rank `1`: `mul[pow[add[x,x],sub[x,1]],pow[add[x,x],add[x,x]]]` RMSE `1.82744`, cosine vs shifted KL `0.706705`
- rank `2`: `add[pow[add[x,x],eml[x,x]],mul[log[x],log[x]]]` RMSE `1.86919`, cosine vs shifted KL `0.702789`
- rank `3`: `add[pow[add[x,x],eml[x,x]],div[sub[1,x],add[x,x]]]` RMSE `1.87849`, cosine vs shifted KL `0.657056`
- rank `4`: `add[pow[add[x,x],eml[x,1]],div[sub[1,x],add[x,x]]]` RMSE `1.89541`, cosine vs shifted KL `0.692012`
- rank `5`: `add[mul[log[x],log[x]],pow[add[x,x],eml[x,1]]]` RMSE `1.90719`, cosine vs shifted KL `0.73504`

## Verdict

The primary theorem prediction fails under this affine count-coordinate test: neither the shifted KL generator nor the unshifted negative-entropy form appears in the top 5 at depth 3.
