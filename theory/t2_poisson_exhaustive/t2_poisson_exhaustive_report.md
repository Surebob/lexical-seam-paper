# T2 Poisson Exhaustive Follow-Up

This follow-up tests whether the redesigned T2 Poisson failure is caused by severity, tail truncation, misspecification choice, or a mismatch between the Poisson KL generator and the pointwise log-density residual given to symbolic regression.

## Test A: rate-ratio sweep

- ratio `1.2`: residual RMSE `1.36946`, winner `mul[eml[mul[x,x],eml[1,x]],eml[mul[x,x],eml[1,x]]]`, max |cosine vs ±KL| `0.0789`, KL rank `6072`.
- ratio `1.5`: residual RMSE `4.32701`, winner `div[add[add[x,x],mul[x,x]],eml[neg[x],mul[x,x]]]`, max |cosine vs ±KL| `0.0738`, KL rank `6376`.
- ratio `2.0`: residual RMSE `10.5558`, winner `div[eml[add[x,x],div[1,x]],eml[neg[x],mul[x,x]]]`, max |cosine vs ±KL| `0.0602`, KL rank `6569`.
- ratio `3.0`: residual RMSE `26.865`, winner `pow[add[eml[x,x],neg[x]],eml[add[x,x],1]]`, max |cosine vs ±KL| `0.1042`, KL rank `6775`.
- ratio `4.0`: residual RMSE `45.4431`, winner `pow[sub[add[x,x],log[x]],eml[add[x,x],add[x,x]]]`, max |cosine vs ±KL| `0.1014`, KL rank `not found`.

## Test B: truncation

The original `(5, 20)` mixture restricted to the 80% empirical-mass support closest to the fitted rate used `21` support points from `3` to `23`. Residual RMSE dropped from `23.6926` to `1.90317`. Winner: `add[pow[mul[x,x],eml[x,x]],sub[neg[x],log[x]]]`; max |cosine vs ±KL| `0.9348`; KL occurrence: `{'step': 3, 'rmse_rank': 1083, 'expr': 'add[add[1,neg[x]],div[log[x],div[1,x]]]', 'math': '((1+(-x))+(log(x)/(1/x)))', 'rmse': 1.792287182737333, 'in_top5_by_rmse': False, 'in_diversity_beam': False, 'diversity_rank': None}`.

## Test C: negative-binomial misspecification

NB(mean=10, dispersion=2) fitted by Poisson(lambda=10) gave residual RMSE `62.5361`. Winner: `pow[sub[add[x,x],log[x]],eml[add[x,x],1]]`; max |cosine vs ±KL| `0.0443`; KL occurrence: `{'step': 3, 'rmse_rank': 7051, 'expr': 'sub[add[1,neg[x]],mul[log[x],neg[x]]]', 'math': '((1+(-x))-(log(x)*(-x)))', 'rmse': 62.53185252307367, 'in_top5_by_rmse': False, 'in_diversity_beam': False, 'diversity_rank': None}`.

## Test D: theoretical expansion reference

For Poisson, the KL/Bregman generator is well-defined in mean-parameter/deviance space: `D(mu||lambda)=mu log(mu/lambda)-mu+lambda`. The SR residual in T2 is instead the pointwise log-density ratio `log q(k)-log p_lambda(k)`. The full distribution KL is the weighted sum of that log-ratio, not the same function as the per-count deviance. This means the predicted KL generator exists, but it is not algebraically identical to the object passed to SR.

Raw cosine between pointwise residual and Poisson deviance: `0.991815`. Centered cosine: `0.987129`. Raw cosine against KL on the SR coordinate: `0.073544`; centered cosine `-0.461781`.

## Verdict

No SR follow-up produced a clean KL winner, but Test D shows the empirical pointwise residual strongly agrees with the Poisson deviance in count/lambda space while strongly disagreeing with KL evaluated on the manuscript-style SR coordinate. This is a methodology-scope finding: the Poisson Bregman prediction is visible in the natural deviance coordinate, but the current density-frame SR encoding does not expose that coordinate to the grammar.
Final interpretation: `methodology_scope_theory_matches_deviance_coordinate_not_sr_coordinate`.
