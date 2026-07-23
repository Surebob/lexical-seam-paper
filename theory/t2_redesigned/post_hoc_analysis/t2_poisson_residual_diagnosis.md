# T2 Poisson Residual Diagnosis

The two Poisson mixture rates are `5.0` and `20.0` with equal weights; the fitted single-Poisson rate is `12.5606`. The component separation is `3.354` in sqrt-largest-lambda units, so this is a deliberately strong two-regime mismatch rather than a near-single-Poisson case.

The full-support log-density residual RMSE is `13.3334`. Restricting to the highest empirical-mass support points carrying about half the probability mass gives RMSE `2.61595`, while the remaining lower-mass support has RMSE `14.4855`. The largest absolute residual occurs at count `56` with residual `33.2975`. This indicates the large residual is not only random noise; the single-Poisson fit is structurally poor, with tail/log-ratio terms contributing strongly to the SR objective.

Interpretation caveat: because the Poisson residual scale is orders of magnitude larger than the Gaussian and Gamma residual scales, depth-3 compositional expressions can reduce RMSE by modeling gross tail structure rather than revealing a simple KL generator.
