# T2 Redesigned Methodology Note

This T2 restart uses the density-estimation frame rather than the rank-frequency frame. Each case generates `N=10000` iid samples from a two-component mixture, fits a single-family density by maximum likelihood, estimates the empirical log density on a grid, and runs the Section 2.4 deterministic search on `log p_empirical - log p_fitted`.

The SR coordinate is a monotone log-index coordinate over the density-evaluation grid: `x = 0.05 + 0.95*log(i)/log(M)` for grid index `i=1..M`. This preserves the manuscript search domain convention while acknowledging that density-estimation support is not word-rank.

The live Section 2.4 implementation includes unary operators `neg, inv, sqr, sqrt, exp, log` and binary operators `eml, add, sub, mul, div, pow`. Beam width is 50, semantic hashing rounds output vectors to ten decimals, `EXP_CLAMP=30`, and `VALUE_ABS_LIMIT=1e6`.

## Reachability

- `gaussian_euclidean` first reachable at step `2` as `mul[sub[1,x],sub[1,x]]` (((1-x)*(1-x))).
- `poisson_kl` first reachable at step `3` as `sub[add[1,neg[x]],mul[log[x],neg[x]]]` (((1+(-x))-(log(x)*(-x)))).
- `gamma_is` first reachable at step `2` as `sub[sub[x,1],log[x]]` (((x-1)-log(x))).
