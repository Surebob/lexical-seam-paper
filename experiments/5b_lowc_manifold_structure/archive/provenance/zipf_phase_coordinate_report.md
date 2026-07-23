# Phase Coordinate Analysis

A weighted-loss sweep interpolates between full-curve RMSE (`lambda = 0`) and top-100-only RMSE (`lambda = 1`).

## High-vs-Low Winner Stability

- high-c/IS block full winners: `{'is': 11}`
- high-c/IS block top-100 winners: `{'is': 10, 'xpow': 1}`
- low-c/exp block full winners: `{'exp': 14}`
- low-c/exp block top-100 winners: `{'xpow': 11, 'is': 2, 'exp': 1}`

## Low-c Phase Coordinate

- median theta over low-c corpora: `105.902295746311` deg
- median head-200 fit R^2 in span{exp,xpow}: `0.850896336469`
- corr(theta, c): `-0.329788`
- corr(theta, transition_fraction): `0.389343`
- corr(theta, lambda_flip): `-0.552754`

## Link To Smooth Synthetic Recovery

- on the low-c family, smooth synthetic modal winner matches empirical full-RMSE winner only `0.071429` of the time
- but it matches the empirical top-100 winner `0.714286` of the time
- single-ZM control matches the empirical top-100 winner `0.571429` of the time

This is the key asymmetry: the high-c side is a stable IS phase, while the low-c side is a rotating manifold where the named winner depends strongly on the scoring functional.

## Per-Corpus Low-c Table

| corpus | c | transition frac | theta (deg) | span R^2 | full winner | top-100 winner | lambda flip exp->xpow |
| --- | ---: | ---: | ---: | ---: | --- | --- | ---: |
| Moby Dick | 10.498 | 0.528755 | 126.479 | 0.856805 | `exp` | `xpow` | 0.85 |
| Grimm's Fairy Tales | 45.983 | 0.672151 | 108.802 | 0.887399 | `exp` | `xpow` | 0.35 |
| Don Quixote | 65.664 | 0.528312 | 2.412 | 0.976825 | `exp` | `is` | 0.95 |
| Pride and Prejudice | 34.119 | 0.483433 | 15.528 | 0.962259 | `exp` | `xpow` | 0.70 |
| Canterbury Tales | 15.251 | 0.501032 | 91.399 | 0.914150 | `exp` | `xpow` | 0.85 |
| Arabian Nights (Vol 1) | 6.425 | 0.471460 | 123.298 | 0.537118 | `exp` | `xpow` | 0.75 |
| Aesop's Fables | 3.688 | 0.658638 | 133.686 | 0.494774 | `exp` | `xpow` | 0.70 |
| Complete Sherlock Holmes | 10.294 | 0.479984 | 95.587 | 0.737621 | `exp` | `xpow` | 0.75 |
| Jane Eyre | 16.308 | 0.498053 | 103.003 | 0.826654 | `exp` | `xpow` | 0.80 |
| Dubliners | 3.263 | 0.709681 | 137.991 | 0.226407 | `exp` | `xpow` | 0.65 |
| Emile | 51.481 | 0.515083 | 4.546 | 0.979460 | `exp` | `is` | 0.85 |
| Ulysses | 0.632 | 0.507735 | -52.375 | 0.585778 | `exp` | `exp` | never |
| Collected Poe | 6.958 | 0.463002 | 110.165 | 0.844988 | `exp` | `xpow` | 0.70 |
| Critique of Pure Reason | 57.958 | 0.688081 | 110.798 | 0.903924 | `exp` | `xpow` | 0.30 |
