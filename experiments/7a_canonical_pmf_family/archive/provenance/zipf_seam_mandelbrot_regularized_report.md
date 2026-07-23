# Regularized Seam-Mandelbrot PMF

- Prior structure came from the smooth reranked family, not from the free PMF fit.
- shared transition law: `k = V^0.521422287447`
- shared width constant for the fixed-kw variant: `w = 0.989992813029`

- held-out winner counts: `{"free": 5, "fixed_k": 3, "fixed_kw": 17}`
- train BIC winner counts: `{"free": 24, "fixed_k": 1}`
- fixed-k beats free on held-out avg NLL: `17` / `25`
- fixed-k,w beats free on held-out avg NLL: `17` / `25`
- median fixed-k minus free held-out avg NLL: `-0.000260615398`
- median fixed-k,w minus free held-out avg NLL: `-0.003798781709`
- free step-2 help count: `2` / `25`
- fixed-k step-2 help count: `5` / `25`
- fixed-k,w step-2 help count: `17` / `25`
- median free step-2 gain: `-0.001600784603`
- median fixed-k step-2 gain: `-0.001312566848`
- median fixed-k,w step-2 gain: `0.004855907600`

| corpus | winner (test) | free | fixed-k | fixed-k,w | free step2 | fixed-k step2 | fixed-k,w step2 | frac |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Dubliners | fixed_kw | 6.550685 | 6.517218 | 6.511116 | -0.002642 | 0.001100 | 0.000149 | 0.521 |
| Collected Poe | fixed_kw | 6.682486 | 6.671851 | 6.660085 | -0.002677 | -0.001499 | -0.000477 | 0.521 |
| Complete Sherlock Holmes | fixed_kw | 6.457765 | 6.453827 | 6.445628 | -0.002463 | -0.000684 | -0.000789 | 0.521 |
| Democracy in America | fixed_kw | 6.316106 | 6.312829 | 6.310677 | -0.001601 | -0.001350 | 0.022608 | 0.521 |
| Ulysses | fixed_kw | 7.356321 | 7.353073 | 7.336533 | -0.001190 | 0.002094 | 0.001912 | 0.521 |
| Arabian Nights (Vol 1) | fixed_kw | 6.689645 | 6.687969 | 6.677410 | -0.002300 | 0.001302 | 0.000215 | 0.521 |
| Aesop's Fables | fixed_kw | 6.232215 | 6.230610 | 6.222484 | -0.001388 | -0.000936 | -0.001966 | 0.521 |
| Emile | fixed_kw | 6.491234 | 6.490527 | 6.487436 | -0.000652 | -0.001289 | 0.031587 | 0.521 |
| Decline and Fall Vol 1 | fixed_kw | 6.718567 | 6.717888 | 6.712761 | 0.000424 | 0.000109 | 0.017640 | 0.521 |
| Origin of Species | fixed_k | 6.375772 | 6.375400 | 6.376035 | -0.002104 | -0.001752 | -0.000565 | 0.521 |
| The Iliad | fixed_kw | 6.880024 | 6.879755 | 6.873835 | -0.001455 | -0.001313 | -0.000248 | 0.521 |
| Les Miserables | fixed_k | 6.802747 | 6.802485 | 6.803990 | -0.000951 | -0.001151 | 0.084744 | 0.521 |
| War and Peace | fixed_k | 6.678952 | 6.678691 | 6.681961 | -0.002324 | -0.000734 | 0.092492 | 0.521 |
| Canterbury Tales | fixed_kw | 6.807097 | 6.806869 | 6.798762 | -0.001521 | -0.001097 | 0.004856 | 0.521 |
| Critique of Pure Reason | fixed_kw | 5.953364 | 5.953298 | 5.952325 | -0.002370 | -0.001915 | -0.000167 | 0.521 |
| Federalist Papers | fixed_kw | 6.200537 | 6.200501 | 6.197656 | -0.002349 | -0.001479 | 0.024521 | 0.521 |
| Principia Ethica | fixed_kw | 5.910708 | 5.910690 | 5.910474 | -0.003097 | -0.001636 | -0.001776 | 0.521 |
| Jane Eyre | fixed_kw | 6.644123 | 6.644193 | 6.634599 | -0.001414 | -0.001483 | 0.000003 | 0.521 |
| Complete Works of Shakespeare | free | 6.882207 | 6.882279 | 6.886567 | -0.001633 | -0.001149 | 0.089781 | 0.521 |
| Don Quixote | free | 6.496717 | 6.497166 | 6.498339 | -0.001711 | -0.001406 | 0.077182 | 0.521 |
| Pride and Prejudice | fixed_kw | 6.336766 | 6.337528 | 6.333796 | -0.002535 | -0.001328 | -0.001348 | 0.521 |
| Wealth of Nations | free | 6.155586 | 6.156624 | 6.156506 | -0.000617 | -0.001567 | 0.080519 | 0.521 |
| Moby Dick | fixed_kw | 6.948097 | 6.951886 | 6.942673 | -0.001394 | -0.001697 | 0.024447 | 0.521 |
| Grimm's Fairy Tales | free | 6.077594 | 6.082621 | 6.080534 | -0.000479 | -0.002937 | 0.011116 | 0.521 |
| King James Bible | free | 6.015869 | 6.021742 | 6.020711 | 0.097341 | 0.000071 | 0.139621 | 0.521 |
