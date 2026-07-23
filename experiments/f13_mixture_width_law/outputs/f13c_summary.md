# F13c — analytic ambiguous-mass A vs measured seam width s

| corpus | basin | A | s_real | A/s | A/V |
|---|---|---:|---:|---:|---:|
| Dubliners | ML | 10 | 65 | 0.14 | 0.0013 |
| Collected Poe | ML | 18 | 101 | 0.18 | 0.0019 |
| Complete Sherlock Holmes | ML | 44 | 86 | 0.51 | 0.0056 |
| Principia Ethica | pred-selected | 49 | 46 | 1.06 | 0.0128 |
| Canterbury Tales | ML | 271 | 199 | 1.36 | 0.0154 |
| Jane Eyre | ML | 219 | 148 | 1.48 | 0.0175 |
| Emile | ML | 234 | 132 | 1.77 | 0.0210 |
| Don Quixote | ML | 345 | 193 | 1.79 | 0.0222 |
| Democracy in America | ML | 249 | 116 | 2.15 | 0.0268 |
| Decline and Fall Vol 1 | ML | 466 | 181 | 2.57 | 0.0292 |
| Arabian Nights (Vol 1) | ML | 390 | 142 | 2.75 | 0.0281 |
| Complete Works of Shakespeare | pred-selected | 942 | 304 | 3.09 | 0.0385 |
| The Iliad | ML | 472 | 146 | 3.23 | 0.0384 |
| Aesop's Fables | ML | 144 | 44 | 3.29 | 0.0328 |
| War and Peace | pred-selected | 865 | 227 | 3.81 | 0.0496 |
| Moby Dick | ML | 768 | 193 | 3.99 | 0.0453 |
| Pride and Prejudice | ML | 336 | 79 | 4.25 | 0.0500 |
| Wealth of Nations | pred-selected | 553 | 124 | 4.45 | 0.0583 |
| Grimm's Fairy Tales | ML | 329 | 58 | 5.68 | 0.0694 |
| Ulysses | ML | 1902 | 274 | 6.94 | 0.0656 |
| King James Bible | pred-selected | 1977 | 163 | 12.11 | 0.1574 |
| Les Miserables | pred-selected | 4474 | 277 | 16.14 | 0.1973 |
| Federalist Papers | pred-selected | 7957 | 110 | 72.56 | 0.9234 |
| Origin of Species | pred-selected | 6857 | 90 | 76.42 | 0.9957 |

- corr(log A, log s_real) = 0.511  (n=24)
- median A/s = 3.16; median A/V = 0.0356 (empirical s/V: 0.0118)
- regression log s ~ log A: slope 0.167, R2 0.261