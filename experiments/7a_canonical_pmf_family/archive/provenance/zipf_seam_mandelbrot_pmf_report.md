# Seam-Mandelbrot PMF Comparison

- Finite-support rank PMFs were fit by train likelihood on an 80/20 token split induced by binomially splitting each type count.
- Models: truncated Zipf, truncated Zipf-Mandelbrot, truncated MOEZipf, and the new finite-support Seam-Mandelbrot PMF.
- The Seam-Mandelbrot PMF uses a smooth transition between two ZM-like slopes and is normalized directly over ranks `1..V`.

- test winner counts: `{"seam": 8, "moe": 5, "zm": 12}`
- train BIC winner counts: `{"seam": 23, "moe": 2}`
- Seam beats MOE on held-out avg NLL: `13` / `25`
- Seam beats ZM on held-out avg NLL: `10` / `25`
- Seam residual step-2 help count: `2` / `25`
- median Seam minus MOE held-out avg NLL: `-0.000052229253`
- median Seam minus ZM held-out avg NLL: `0.000745891042`
- median Seam step-2 gain: `-0.001600784603`
- median Seam transition fraction: `0.640662027425`

| corpus | winner (test) | Zipf | ZM | MOE | Seam | winner (BIC) | Seam RMSE | MOE RMSE | Seam step-2 | gain | frac |
| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: |
| Wealth of Nations | seam | 6.195767 | 6.175560 | 6.161545 | 6.155586 | seam | 0.199497 | 0.483896 | eml[sub[x,1],eml[x,1]] | -0.000617 | 0.375 |
| Canterbury Tales | zm | 6.850582 | 6.803569 | 6.812279 | 6.807097 | seam | 0.124492 | 0.179158 | eml[sub[x,1],eml[x,1]] | -0.001521 | 0.417 |
| Jane Eyre | zm | 6.685154 | 6.639199 | 6.649151 | 6.644123 | seam | 0.132976 | 0.167017 | eml[sub[x,1],eml[x,1]] | -0.001414 | 0.868 |
| Ulysses | zm | 7.366644 | 7.351230 | 7.360040 | 7.356321 | seam | 0.125763 | 0.219711 | eml[sub[x,1],eml[x,1]] | -0.001190 | 0.785 |
| Complete Works of Shakespeare | seam | 6.987294 | 6.886706 | 6.885901 | 6.882207 | seam | 0.111380 | 0.262536 | eml[sub[x,1],eml[x,1]] | -0.001633 | 0.637 |
| The Iliad | seam | 6.895521 | 6.890807 | 6.883405 | 6.880024 | seam | 0.161066 | 0.336733 | eml[sub[x,1],eml[x,1]] | -0.001455 | 0.535 |
| Principia Ethica | seam | 5.978310 | 5.917477 | 5.913924 | 5.910708 | seam | 0.111742 | 0.305848 | eml[sub[x,1],eml[x,1]] | -0.003097 | 0.469 |
| Critique of Pure Reason | seam | 6.006653 | 5.967963 | 5.956428 | 5.953364 | seam | 0.114402 | 0.353394 | eml[sub[x,1],eml[x,1]] | -0.002370 | 0.350 |
| Decline and Fall Vol 1 | seam | 6.721171 | 6.721170 | 6.721115 | 6.718567 | seam | 0.183424 | 0.396220 | sub[pow[x,x],sqrt[x]] | 0.000424 | 0.652 |
| King James Bible | seam | 6.080495 | 6.035166 | 6.017965 | 6.015869 | seam | 0.346382 | 0.468982 | eml[neg[x],add[x,x]] | 0.097341 | 0.380 |
| Pride and Prejudice | zm | 6.413204 | 6.332858 | 6.338649 | 6.336766 | seam | 0.114128 | 0.160457 | eml[sub[x,1],eml[x,1]] | -0.002535 | 0.529 |
| Moby Dick | zm | 6.963919 | 6.947351 | 6.949106 | 6.948097 | moe | 0.116651 | 0.161470 | sub[pow[x,x],sqrt[x]] | -0.001394 | 0.411 |
| Grimm's Fairy Tales | seam | 6.128785 | 6.088524 | 6.077646 | 6.077594 | seam | 0.207186 | 0.229029 | eml[sub[x,1],eml[x,1]] | -0.000479 | 0.324 |
| Don Quixote | zm | 6.556167 | 6.496280 | 6.496303 | 6.496717 | seam | 0.146540 | 0.175435 | eml[sub[x,1],eml[x,1]] | -0.001711 | 0.820 |
| Complete Sherlock Holmes | zm | 6.493105 | 6.447377 | 6.456989 | 6.457765 | seam | 0.115579 | 0.185280 | eml[sub[x,1],eml[x,1]] | -0.002463 | 0.804 |
| Origin of Species | moe | 6.406564 | 6.383434 | 6.374942 | 6.375772 | seam | 0.141874 | 0.326233 | eml[sub[x,1],eml[x,1]] | -0.002104 | 0.561 |
| War and Peace | moe | 6.728315 | 6.685807 | 6.678028 | 6.678952 | seam | 0.117931 | 0.314596 | sub[pow[x,x],sqrt[x]] | -0.002324 | 0.580 |
| Emile | zm | 6.550809 | 6.487122 | 6.489192 | 6.491234 | seam | 0.161567 | 0.175332 | sub[pow[x,x],sqrt[x]] | -0.000652 | 0.704 |
| Les Miserables | moe | 6.827572 | 6.802444 | 6.800382 | 6.802747 | seam | 0.124938 | 0.276868 | eml[sub[x,1],eml[x,1]] | -0.000951 | 0.641 |
| Arabian Nights (Vol 1) | zm | 6.706086 | 6.685473 | 6.686671 | 6.689645 | seam | 0.116640 | 0.188905 | sub[pow[x,x],sqrt[x]] | -0.002300 | 0.821 |
| Federalist Papers | moe | 6.207761 | 6.199416 | 6.196934 | 6.200537 | seam | 0.112724 | 0.319281 | eml[sub[x,1],eml[x,1]] | -0.002349 | 0.787 |
| Democracy in America | moe | 6.318432 | 6.314584 | 6.311806 | 6.316106 | seam | 0.148497 | 0.330385 | eml[sub[x,1],eml[x,1]] | -0.001601 | 0.801 |
| Aesop's Fables | zm | 6.226453 | 6.221751 | 6.223659 | 6.232215 | moe | 0.157197 | 0.169042 | eml[sub[x,1],eml[x,1]] | -0.001388 | 0.773 |
| Collected Poe | zm | 6.668311 | 6.664039 | 6.666312 | 6.682486 | seam | 0.106380 | 0.172245 | eml[sub[x,1],eml[x,1]] | -0.002677 | 0.822 |
| Dubliners | zm | 6.541973 | 6.514009 | 6.521278 | 6.550685 | seam | 0.113508 | 0.204987 | eml[sub[x,1],eml[x,1]] | -0.002642 | 0.960 |
