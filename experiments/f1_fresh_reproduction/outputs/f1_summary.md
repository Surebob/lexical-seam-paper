# F1 fresh reproduction + proper re-experiments — summary

- corpora run: 25/25
- token-count exact matches vs canon: 25/25
- max |c - canon_c|: 0.000000
- max |zm_rmse - canon_rmse|: 0.00000000
- unit-amplitude winner-family matches: 25/25
- continuous refit beats grid by >1e-4: 0/25 (min delta -0.000001)

## Amplitude-fair (free-lambda) winner map
- is: 17/25
- exp: 0/25
- euclid: 1/25
- xpow: 7/25
- corpora where free winner != canonical family: 14

## lambda-ZM 4-parameter model (joint OLS fit)
- winner (is vs exp): is 0, exp 25
- RMSE improvement over ZM: min 2.581% median 10.219% max 23.010%
- BIC(4-param) beats BIC(ZM): 25/25
- corr(lzm_is_lam, c): 0.7487
- corr(lzm_is_lam, hapax_ratio): -0.4845
- corr(lzm_is_c minus zm_c shift): median 92.371

## Exponent gap (A2 rerun)
- corr(exp_gap, c): -0.6891
- corr(exp_gap, is_winner_flag): -0.5089
- corr(c, is_winner_flag): 0.8828

## Lexical confound (A3 rerun)
- corr(hapax, is_winner_flag): -0.6529
- corr(ttr, c): -0.7060
- corr(tokens, c): 0.5945  (log-log)