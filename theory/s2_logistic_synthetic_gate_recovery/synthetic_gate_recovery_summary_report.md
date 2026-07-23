# Synthetic Logistic-Gate Recovery Check

- synthetic generator: exact decoupled logistic model using local S2 v3 logistic parameters
- corpora: `1`
- logistic wins: `1`
- erf wins: `0`
- algebraic wins: `0`
- arctan wins: `0`
- tanh calibration pass count: `0`
- erf beats logistic count: `0`
- median BIC(erf - logistic): `243143.97512281686`

Interpretation: erf does not win on logistic-generated synthetic data, so the empirical erf preference is not explained by a generic fitter artifact in this recovery check.
