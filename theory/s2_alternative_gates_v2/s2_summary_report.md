# S2 v2 Alternative Gate Fit Sweep

- Sanity check summary remains in `s2_pre_fit_sanity_check_v2.csv`.
- Full fit sweep uses five gates: logistic, tanh, erf, algebraic, arctan.
- Bounds widened only in `w`: `[0.05, 10.0]`; all other constraints and optimizer settings match the historical reranked smooth-model fit.

- tanh calibration pass count: `0` / `25`
- max |BIC_tanh - BIC_logistic|: `1296.979198424975`
- max |w_tanh/(2*w_logistic) - 1|: `10.358975492935`

## Tanh calibration failures

- `Complete Works of Shakespeare`
- `War and Peace`
- `Moby Dick`
- `King James Bible`
- `Federalist Papers`
- `Grimm's Fairy Tales`
- `Don Quixote`
- `Pride and Prejudice`
- `Canterbury Tales`
- `Arabian Nights (Vol 1)`
- `Aesop's Fables`
- `Complete Sherlock Holmes`
- `Jane Eyre`
- `Dubliners`
- `The Iliad`
- `Democracy in America`
- `Origin of Species`
- `Wealth of Nations`
- `Les Miserables`
- `Decline and Fall Vol 1`
- `Emile`
- `Ulysses`
- `Collected Poe`
- `Principia Ethica`
- `Critique of Pure Reason`

- Interpretation: more than 2 corpora fail tanh calibration, so the optimizer-correctness check does not pass and the gate-comparison results should not yet be interpreted.
