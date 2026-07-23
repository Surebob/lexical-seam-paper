# S2 Alternative Smooth Gate Families

- Alternative gates rerun with the historical reranked smooth-fit engine and optimizer settings.
- Logistic baseline BIC is computed from canonical 3a RMSE with the historical BIC formula `p * ln(n) + n * ln(MSE)` using `p = 8` and `n = vocab_size`.

## Winner counts

- logistic: `24`
- tanh: `0`
- erf: `1`
- algebraic: `0`

- median BIC spread: `495.826268932411`
- gates indistinguishable count (`spread < 10`): `0`

- Interpretation: logistic wins on 20+ corpora, so the specific functional form is empirically supported.
