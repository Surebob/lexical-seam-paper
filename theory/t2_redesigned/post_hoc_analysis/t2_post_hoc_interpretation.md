# T2 Post-Hoc Interpretation

This post-hoc audit does not modify the redesigned T2 outputs. It checks whether the winning expression is sign-equivalent to the predicted generator, whether any top-20 expressions satisfy numerical Bregman boundary conditions, and whether the Poisson miss is dominated by residual-scale pathology.

## Gaussian

- revised verdict: `BREGMAN_CLASS_PRESENT_BUT_PREDICTED_NOT_SELECTED`
- winner: `sub[sqrt[x],pow[x,x]]` ((sqrt(x)-pow(x,x)))
- winner sign-cosine verdict: `not_shape_equivalent`
- top-20 Bregman-class expression count: `3`
- predicted generator occurrence: step `2`, RMSE rank `16`, RMSE `0.218428`

## Poisson

- revised verdict: `NO_TOP20_BREGMAN_SUPPORT`
- winner: `div[pow[eml[x,1],eml[x,1]],pow[add[x,x],neg[x]]]` ((pow(EML(x,1),EML(x,1))/pow((x+x),(-x))))
- winner sign-cosine verdict: `not_shape_equivalent`
- top-20 Bregman-class expression count: `0`
- predicted generator occurrence: step `3`, RMSE rank `6436`, RMSE `13.316`

## Gamma

- revised verdict: `SIGN_CONVENTION_PASS`
- winner: `add[log[x],sub[1,x]]` ((log(x)+(1-x)))
- winner sign-cosine verdict: `equivalent_to_predicted_up_to_sign`
- top-20 Bregman-class expression count: `0`
- predicted generator occurrence: step `2`, RMSE rank `49`, RMSE `0.340191`

Summary: Gamma is upgraded from a raw fail to a sign-convention pass because the winner is exactly the negative IS form under the current residual direction. Gaussian remains a close but non-equivalent miss: Euclidean is present at rank 16, but the winner is not cosine-equivalent to Euclidean up to sign. Poisson remains the cleanest miss; KL is reachable but deeply ranked, and the residual diagnosis shows a severe two-regime mismatch with large tail/log-ratio structure.
