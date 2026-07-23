# Simulation Recovery Test

Synthetic corpora were sampled from fitted rank-probability curves, then re-analyzed with single ZM plus exact step-2 search.

- replicates per corpus per generator: `10`
- corpora analyzed: `25`

## Overall Recovery

- empirical median winner-vs-Euclidean gap (full): `0.005957768526`
- empirical median winner-vs-Euclidean gap (top-100): `0.092759376913`

### Smooth-model synthetic corpora

- replicate exact winner match rate: `0.484000`
- majority winner match count: `12` / `25`
- replicate step-2 help rate: `0.588000`
- median winner-vs-Euclidean gap (full): `0.007826478737`
- median winner-vs-Euclidean gap (top-100): `0.037112483001`
- basis corr linear/quadratic/cubic: `{'linear_u': 0.892854405891926, 'quadratic_u2': 0.8794468873817035, 'cubic_u3': 0.806545507783293}`

### Single-ZM synthetic controls

- replicate exact winner match rate: `0.200000`
- majority winner match count: `5` / `25`
- replicate step-2 help rate: `0.000000`
- median winner-vs-Euclidean gap (full): `0.007714168828`
- median winner-vs-Euclidean gap (top-100): `0.095933590738`
- basis corr linear/quadratic/cubic: `{'linear_u': 0.37542454418261295, 'quadratic_u2': 0.48403455528452977, 'cubic_u3': 0.3417143411251593}`

## Per-Corpus Winner Recovery

| corpus | empirical winner | smooth modal (share) | smooth exact-match rate | smooth help rate | ZM modal (share) | ZM exact-match rate | ZM help rate |
| --- | --- | --- | ---: | ---: | --- | ---: | ---: |
| Complete Works of Shakespeare | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| War and Peace | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Moby Dick | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| King James Bible | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Federalist Papers | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Grimm's Fairy Tales | `eml[sub[x,1],eml[x,1]]` | `sub[sub[x,1],log[x]]` (1.00) | 0.00 | 0.90 | `eml[sub[x,1],eml[x,1]]` (1.00) | 1.00 | 0.00 |
| Don Quixote | `eml[sub[x,1],eml[x,1]]` | `sub[sub[x,1],log[x]]` (1.00) | 0.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 1.00 | 0.00 |
| Pride and Prejudice | `eml[sub[x,1],eml[x,1]]` | `eml[sub[x,1],eml[x,1]]` (1.00) | 1.00 | 0.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 1.00 | 0.00 |
| Canterbury Tales | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Arabian Nights (Vol 1) | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Aesop's Fables | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Complete Sherlock Holmes | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Jane Eyre | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Dubliners | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (0.90) | 0.10 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| The Iliad | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Democracy in America | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Origin of Species | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Wealth of Nations | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Les Miserables | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Decline and Fall Vol 1 | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Emile | `eml[sub[x,1],eml[x,1]]` | `sub[sub[x,1],log[x]]` (1.00) | 0.00 | 0.80 | `eml[sub[x,1],eml[x,1]]` (1.00) | 1.00 | 0.00 |
| Ulysses | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Collected Poe | `eml[sub[x,1],eml[x,1]]` | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 | `sub[sqrt[x],pow[x,x]]` (1.00) | 0.00 | 0.00 |
| Principia Ethica | `sub[sub[x,1],log[x]]` | `sub[sub[x,1],log[x]]` (1.00) | 1.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 0.00 | 0.00 |
| Critique of Pure Reason | `eml[sub[x,1],eml[x,1]]` | `sub[sub[x,1],log[x]]` (1.00) | 0.00 | 1.00 | `eml[sub[x,1],eml[x,1]]` (1.00) | 1.00 | 0.00 |
