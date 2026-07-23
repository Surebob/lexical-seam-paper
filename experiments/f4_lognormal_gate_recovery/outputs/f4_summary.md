# F4 gate recovery on synthetic population mixtures

| dataset | kind | winner | erf_bic - best_other | spread(min..max BIC gap) |
|---|---|---|---:|---:|
| LN1 s1 | lognormal | algebraic | 30.04 | 135.16 |
| LN1 s2 | lognormal | algebraic | 66.41 | 66.41 |
| LN2 s1 | lognormal | erf | -113.56 | 562.19 |
| LN2 s2 | lognormal | logistic | 34.30 | 263.27 |
| LN3 s1 | lognormal | algebraic | 325.96 | 325.96 |
| LN3 s2 | lognormal | arctan | 301.35 | 301.35 |
| LN4 s1 | lognormal | erf | -49.32 | 210.51 |
| LN4 s2 | lognormal | erf | -1.38 | 52.79 |
| PA1 s1 | pareto | erf | -26.45 | 125.77 |
| PA1 s2 | pareto | algebraic | 60.31 | 60.31 |
| PA2 s1 | pareto | erf | -25.49 | 72.55 |
| PA2 s2 | pareto | erf | -23.69 | 71.31 |
| PA3 s1 | pareto | arctan | 185.84 | 185.84 |
| PA3 s2 | pareto | arctan | 190.85 | 190.85 |
| SL1 s1 | single_lognormal | erf | -342.41 | 789.92 |
| SL1 s2 | single_lognormal | erf | -181.20 | 362.27 |

## Winner counts by generative family
- lognormal -> algebraic: 3
- lognormal -> arctan: 1
- lognormal -> erf: 3
- lognormal -> logistic: 1
- pareto -> algebraic: 1
- pareto -> arctan: 2
- pareto -> erf: 3
- single_lognormal -> erf: 2