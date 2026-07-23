# F13 — does the fitted PLN mixture contain the width law?

| corpus | V_real | V_synth | s_real | s_synth | s/V real | s/V synth | ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| Complete Works of Shakespeare | 24458 | 24475 | 304.5 | 4594.3 | 0.0124 | 0.1877 | 15.09 |
| War and Peace | 17445 | 17482 | 227.2 | 4131.9 | 0.0130 | 0.2364 | 18.19 |
| Moby Dick | 16956 | 17117 | 192.6 | 196.0 | 0.0114 | 0.0115 | 1.02 |
| King James Bible | 12555 | 12551 | 163.3 | 2973.4 | 0.0130 | 0.2369 | 18.21 |
| Federalist Papers | 8617 | 8648 | 109.7 | 663.0 | 0.0127 | 0.0767 | 6.05 |
| Grimm's Fairy Tales | 4746 | 4751 | 58.0 | 60.4 | 0.0122 | 0.0127 | 1.04 |
| Don Quixote | 15574 | 15512 | 192.8 | 184.5 | 0.0124 | 0.0119 | 0.96 |
| Pride and Prejudice | 6729 | 6770 | 79.2 | 71.0 | 0.0118 | 0.0105 | 0.90 |
| Canterbury Tales | 17653 | 17700 | 198.7 | 184.7 | 0.0113 | 0.0104 | 0.93 |
| Arabian Nights (Vol 1) | 13883 | 13945 | 142.0 | 126.7 | 0.0102 | 0.0091 | 0.89 |
| Aesop's Fables | 4394 | 4330 | 43.8 | 45.3 | 0.0100 | 0.0105 | 1.03 |
| Complete Sherlock Holmes | 7817 | 7879 | 86.1 | 85.8 | 0.0110 | 0.0109 | 1.00 |
| Jane Eyre | 12541 | 12496 | 148.3 | 138.5 | 0.0118 | 0.0111 | 0.93 |
| Dubliners | 7157 | 7204 | 65.4 | 77.8 | 0.0091 | 0.0108 | 1.19 |
| The Iliad | 12300 | 12388 | 146.5 | 163.2 | 0.0119 | 0.0132 | 1.11 |
| Democracy in America | 9301 | 9309 | 116.1 | 105.1 | 0.0125 | 0.0113 | 0.91 |
| Origin of Species | 6887 | 6831 | 89.7 | 412.4 | 0.0130 | 0.0604 | 4.60 |
| Wealth of Nations | 9478 | 9442 | 124.4 | 618.9 | 0.0131 | 0.0655 | 4.98 |
| Les Miserables | 22677 | 22576 | 277.1 | 4417.4 | 0.0122 | 0.1957 | 15.94 |
| Decline and Fall Vol 1 | 15983 | 15960 | 181.4 | 182.7 | 0.0114 | 0.0114 | 1.01 |
| Emile | 11174 | 11167 | 132.1 | 125.0 | 0.0118 | 0.0112 | 0.95 |
| Ulysses | 28990 | 28949 | 274.1 | 286.9 | 0.0095 | 0.0099 | 1.05 |
| Collected Poe | 9312 | 9393 | 101.0 | 105.5 | 0.0108 | 0.0112 | 1.04 |
| Principia Ethica | 3814 | 3814 | 46.2 | 27722.9 | 0.0121 | 7.2687 | 600.64 |
| Critique of Pure Reason | 6537 | 6559 | 72.7 | 45427.3 | 0.0111 | 6.9260 | 624.78 |

- corr(log s_synth, log s_real) = 0.143
- median s_synth/s_real = 1.045
- median s_synth/V_synth = 0.0115  (real-corpus law: 0.0118)
- law within synthetic family: log s ~ log V slope = 0.272, R2 = 0.006

**Reading: NOT (fully) CONTAINED — the single-depth mixture does not reproduce the width law as measured; the law requires an ingredient beyond Section 3.5's generative account.**