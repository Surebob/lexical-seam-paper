# Low-c Manifold Analysis

These are the 14 English corpora whose empirical step-2 winner is the exponential Bregman generator `exp(x-1)-x`.

- median cosine between `exp(x-1)-x` and `x^x-sqrt(x)` over top-200 head coordinates: `0.965001758614`
- median head-200 R^2 of span{exp, xpow}: `0.835717144618`
- median head-200 R^2 of span{exp, IS}: `0.565392470978`
- median `(xpow - exp)` RMSE on full curve: `0.002566904494`
- median `(xpow - exp)` RMSE on top-100: `-0.074708721124`

## Winner Counts On The Same 14 Corpora

- full RMSE winners: `{'exp': 14}`
- top-50 RMSE winners: `{'xpow': 9, 'is': 4, 'exp': 1}`
- top-100 RMSE winners: `{'xpow': 11, 'is': 2, 'exp': 1}`
- top-200 RMSE winners: `{'xpow': 13, 'exp': 1}`

## Link To Simulation Recovery

- smooth synthetic exact winner match rate over all corpora: `0.484000`
- smooth synthetic exact winner match rate on the empirical low-c family is only `0.078571`
- smooth synthetic modal winners on the empirical low-c family: `{'sub[sqrt[x],pow[x,x]]': 9, 'sub[sub[x,1],log[x]]': 4, 'eml[sub[x,1],eml[x,1]]': 1}`

This says the low-c side is not a clean exp-vs-not-exp phase split. It is a near-degenerate head manifold where full-curve scoring picks `exp`, but head-focused scoring usually picks `x^x-sqrt(x)`.

## Per-Corpus Table

| corpus | c | cosine(exp,xpow) | full winner | top-100 winner | xpow-exp full | xpow-exp top100 | R^2 span{exp,xpow} | R^2 span{exp,IS} |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| Moby Dick | 10.498 | 0.968492 | `exp` | `xpow` | 0.003358 | -0.088661 | 0.663963 | 0.825199 |
| Grimm's Fairy Tales | 45.983 | 0.957765 | `exp` | `xpow` | 0.001552 | -0.076271 | 0.879598 | 0.699409 |
| Don Quixote | 65.664 | 0.967965 | `exp` | `is` | 0.002833 | -0.015572 | 0.904863 | 0.615724 |
| Pride and Prejudice | 34.119 | 0.961452 | `exp` | `xpow` | 0.002568 | -0.070584 | 0.881810 | 0.430398 |
| Canterbury Tales | 15.251 | 0.968734 | `exp` | `xpow` | 0.002746 | -0.088635 | 0.905551 | 0.605517 |
| Arabian Nights (Vol 1) | 6.425 | 0.967221 | `exp` | `xpow` | 0.001783 | -0.079092 | 0.500572 | 0.484593 |
| Aesop's Fables | 3.688 | 0.956844 | `exp` | `xpow` | 0.002369 | -0.057444 | 0.102236 | 0.362439 |
| Complete Sherlock Holmes | 10.294 | 0.962827 | `exp` | `xpow` | 0.002872 | -0.085935 | 0.730695 | 0.433116 |
| Jane Eyre | 16.308 | 0.966529 | `exp` | `xpow` | 0.003453 | -0.084927 | 0.826458 | 0.570921 |
| Dubliners | 3.263 | 0.962031 | `exp` | `xpow` | 0.001139 | -0.049677 | 0.043688 | 0.075547 |
| Emile | 51.481 | 0.965703 | `exp` | `is` | 0.002566 | -0.030985 | 0.899910 | 0.559864 |
| Ulysses | 0.632 | 0.971412 | `exp` | `exp` | 0.000560 | 0.018047 | 0.483741 | 0.541323 |
| Collected Poe | 6.958 | 0.964300 | `exp` | `xpow` | 0.003161 | -0.093314 | 0.844976 | 0.583867 |
| Critique of Pure Reason | 57.958 | 0.961174 | `exp` | `xpow` | 0.000998 | -0.073146 | 0.894128 | 0.710351 |
