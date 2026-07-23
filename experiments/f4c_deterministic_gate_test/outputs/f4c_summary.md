# F4c deterministic gate discrimination — summary

| corpus | family | noise | winner | erf margin | spread |
|---|---|---|---|---:|---:|
| Complete Works of Shakespeare | lognormal_mix | off | erf | -439.0 | 605.2 |
| Complete Works of Shakespeare | lognormal_mix | on | logistic | 220.3 | 220.3 |
| Complete Works of Shakespeare | pareto_mix | off | erf | -86.8 | 610.3 |
| Complete Works of Shakespeare | pareto_mix | on | erf | -81.9 | 698.4 |
| War and Peace | lognormal_mix | off | erf | -76.7 | 186.8 |
| War and Peace | lognormal_mix | on | erf | -17.6 | 126.4 |
| War and Peace | pareto_mix | off | logistic | 12.3 | 202.3 |
| War and Peace | pareto_mix | on | logistic | 26.3 | 212.4 |
| King James Bible | lognormal_mix | off | erf | -64.1 | 160.0 |
| King James Bible | lognormal_mix | on | erf | -113.4 | 451.3 |
| King James Bible | pareto_mix | off | erf | -77.5 | 331.6 |
| King James Bible | pareto_mix | on | erf | -72.8 | 391.7 |
| Federalist Papers | lognormal_mix | off | erf | -149.9 | 1269.8 |
| Federalist Papers | lognormal_mix | on | erf | -36.8 | 99.4 |
| Federalist Papers | pareto_mix | off | erf | -47.7 | 173.7 |
| Federalist Papers | pareto_mix | on | erf | -49.5 | 244.1 |
| Origin of Species | lognormal_mix | off | arctan | 19.4 | 605.3 |
| Origin of Species | lognormal_mix | on | erf | -32.0 | 98.6 |
| Origin of Species | pareto_mix | off | erf | -66.1 | 329.4 |
| Origin of Species | pareto_mix | on | erf | -60.5 | 322.8 |
| Wealth of Nations | lognormal_mix | off | erf | -61.2 | 164.1 |
| Wealth of Nations | lognormal_mix | on | erf | -74.1 | 268.3 |
| Wealth of Nations | pareto_mix | off | erf | -75.2 | 340.0 |
| Wealth of Nations | pareto_mix | on | erf | -80.3 | 426.1 |
| Moby Dick | lognormal_mix | off | erf | -674.7 | 986.4 |
| Moby Dick | lognormal_mix | on | erf | -536.9 | 934.6 |
| Moby Dick | pareto_mix | off | algebraic | 192.2 | 205.0 |
| Moby Dick | pareto_mix | on | algebraic | 129.5 | 129.5 |
| Pride and Prejudice | lognormal_mix | off | arctan | 20.9 | 102.8 |
| Pride and Prejudice | lognormal_mix | on | erf | -22.7 | 69.0 |
| Pride and Prejudice | pareto_mix | off | erf | -90.7 | 458.7 |
| Pride and Prejudice | pareto_mix | on | erf | -46.7 | 243.6 |
| Dubliners | lognormal_mix | off | arctan | 93.2 | 127.9 |
| Dubliners | lognormal_mix | on | erf | -24.0 | 124.0 |
| Dubliners | pareto_mix | off | erf | -43.2 | 198.0 |
| Dubliners | pareto_mix | on | arctan | 13.5 | 13.5 |
| Ulysses | lognormal_mix | off | erf | -1517.9 | 2334.2 |
| Ulysses | lognormal_mix | on | erf | -587.1 | 1759.8 |
| Ulysses | pareto_mix | off | algebraic | 1044.1 | 1044.1 |
| Ulysses | pareto_mix | on | arctan | 216.1 | 216.1 |
| Grimm's Fairy Tales | lognormal_mix | off | logistic | 609.7 | 609.7 |
| Grimm's Fairy Tales | lognormal_mix | on | erf | -27.5 | 93.2 |
| Grimm's Fairy Tales | pareto_mix | off | erf | -65.8 | 329.8 |
| Grimm's Fairy Tales | pareto_mix | on | erf | -35.6 | 177.1 |
| Don Quixote | lognormal_mix | off | erf | -495.9 | 844.4 |
| Don Quixote | lognormal_mix | on | erf | -3.6 | 17.2 |
| Don Quixote | pareto_mix | off | logistic | 17.6 | 173.3 |
| Don Quixote | pareto_mix | on | logistic | 27.6 | 114.1 |

## Winner counts (family, noise) -> gate
- lognormal_mix noise=off -> arctan: 3
- lognormal_mix noise=off -> erf: 8
- lognormal_mix noise=off -> logistic: 1
- lognormal_mix noise=on -> erf: 11
- lognormal_mix noise=on -> logistic: 1
- pareto_mix noise=off -> algebraic: 2
- pareto_mix noise=off -> erf: 8
- pareto_mix noise=off -> logistic: 2
- pareto_mix noise=on -> algebraic: 1
- pareto_mix noise=on -> arctan: 2
- pareto_mix noise=on -> erf: 7
- pareto_mix noise=on -> logistic: 2