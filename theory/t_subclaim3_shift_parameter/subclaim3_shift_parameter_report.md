# Subclaim 3 Shift-Parameter Fit

## Protocol

Fit `f(x) = x^x - sqrt(x) - λ(x - 1)` to the top-200 single-ZM log-frequency residual with λ as the only free parameter.
The normalized coordinate is the paper's log-rank coordinate over the full vocabulary, `x = 0.05 + 0.95 log(r) / log(V)`.
English targets are the canonical Section 3.1 low-c corpora with `c < 66` and step-2 winner `exp(x - 1) - x`; multilingual targets are Russian, Mandarin, and Arabic from the cached multilingual corpus run.

## Pre-flight

- English low-c targets loaded: `14`.
- Multilingual c≈0 targets loaded: `3`.
- Skipped targets: `0`.
- Minimum analytic/numerical `f''` on `[0.05, 1]`: `1.861240463098`.

## Aggregate Results

- N fitted corpora: `17`.
- λ mean / median / std: `-0.040620` / `-0.124009` / `0.271259`.
- λ range: `-0.319614` (Collected Poe) to `0.781753` (Romance of the Three Kingdoms (Chinese, Gutenberg 23950)).
- Pearson r(λ, single-ZM c): `-0.246179`.
- Pearson r(λ, V): `0.614078`.
- Bregman boundary conditions satisfied at fitted λ: `0` / `17`.

## Group Breakdown

- English low-c λ mean / median / std: `-0.143309` / `-0.157006` / `0.103730`; range `-0.319614` to `0.060771`.
- English low-c Pearson r(λ, c) / r(λ, V): `0.317200` / `0.326454`.
- Multilingual c≈0 λ mean / median / std: `0.438595` / `0.377042` / `0.316897`; range `0.156990` to `0.781753`.
- Multilingual c≈0 Pearson r(λ, V): `-0.977198`; r(λ, c) is undefined because all three fitted c values are zero.

## Interpretation Against Requested Framework

The fitted λ values vary systematically by the requested correlation criterion (`|r| > 0.5` for c or V).
Because `f'(1) = 0.5 - λ`, exact Bregman recentering requires λ = 0.5; the free residual fits should therefore be read as a data-fit diagnostic, not as automatically Bregman-valid generators.

## Per-Corpus Summary

| Corpus | Group | V | c | λ | R² | RMSE | Bregman conditions? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Moby Dick | english_low_c | 16956 | 10.497932 | -0.138679 | 0.579201 | 0.082832 | False |
| Grimm's Fairy Tales | english_low_c | 4746 | 45.982513 | -0.124009 | 0.412294 | 0.195731 | False |
| Don Quixote | english_low_c | 15574 | 65.663942 | 0.060771 | 0.319097 | 0.364535 | False |
| Pride and Prejudice | english_low_c | 6729 | 34.119421 | -0.194676 | 0.391453 | 0.186198 | False |
| Canterbury Tales | english_low_c | 17653 | 15.250535 | -0.242745 | 0.534988 | 0.117169 | False |
| Arabian Nights (Vol 1) | english_low_c | 13883 | 6.425403 | -0.232128 | 0.469278 | 0.078683 | False |
| Aesop's Fables | english_low_c | 4394 | 3.687693 | -0.076346 | -0.539780 | 0.091663 | False |
| Complete Sherlock Holmes | english_low_c | 7817 | 10.294480 | -0.196531 | 0.633764 | 0.069213 | False |
| Jane Eyre | english_low_c | 12541 | 16.307585 | -0.214242 | 0.472597 | 0.140484 | False |
| Dubliners | english_low_c | 7157 | 3.263282 | -0.095797 | -1.195603 | 0.107093 | False |
| Emile | english_low_c | 11174 | 51.480679 | -0.075106 | 0.314016 | 0.310805 | False |
| Ulysses | english_low_c | 28990 | 0.632415 | 0.018108 | -4.073463 | 0.126517 | False |
| Collected Poe | english_low_c | 9312 | 6.957528 | -0.319614 | 0.488129 | 0.112893 | False |
| Critique of Pure Reason | english_low_c | 6537 | 57.958304 | -0.175332 | 0.322100 | 0.263544 | False |
| War and Peace (Russian, Wikisource) | multilingual_c0 | 26011 | 0.000000 | 0.156990 | -3.210809 | 0.176921 | False |
| Romance of the Three Kingdoms (Chinese, Gutenberg 23950) | multilingual_c0 | 23483 | 0.000000 | 0.781753 | -3.188874 | 0.213243 | False |
| One Thousand and One Nights (Arabic, Wikisource) | multilingual_c0 | 25629 | 0.000000 | 0.377042 | -3.075605 | 0.208810 | False |
