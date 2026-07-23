# T2 Coordinate-Theorem Follow-Up

This follow-up tests whether Poisson KL is present but subdominant, whether milder mixtures recover it, and whether a pure dispersion misspecification exposes it more cleanly.

## Test 1: high-cosine expression

The highest-cosine expression versus shifted KL is rank `3727` by RMSE: `div[mul[log[x],log[x]],div[1,sqrt[x]]]` (((log(x)*log(x))/(1/sqrt(x)))). Its cosine versus shifted KL is `0.999408` and its RMSE is `2.94114`.
It is not algebraically equivalent to shifted KL under affine fit: max residual `0.0361485`, R^2 `0.998724`.
Bregman boundary pass: `False`; integer-support discrete convexity pass: `True`.

## Test 2: KL component fit and remainder

Least-squares fit `a + b*(x log x - x + 1)` gives `a=1.23619`, `b=2.3468`, R^2 `0.0386756`, and RMSE `2.69208` versus original residual RMSE `3.26063`.
Remainder SR depth-2 winner: `mul[mul[x,x],eml[x,1]]` with RMSE `2.41543`.

## Test 3: mild mixtures in linear count coordinate

- rates `(10, 12)`: residual RMSE `0.171404`, winner `add[mul[mul[x,x],log[x]],eml[sub[x,1],add[1,x]]]`, predicted form in top 5 `False`, in top 20 `False`.
- rates `(10, 15)`: residual RMSE `0.757964`, winner `add[sub[sub[x,1],log[x]],pow[mul[x,x],eml[x,x]]]`, predicted form in top 5 `False`, in top 20 `False`.

## Test 4: negative-binomial dispersion misspecification

NB(mean=12, variance=24) fit by Poisson used NB size/dispersion `12`. Residual RMSE `1.492`; winner `sub[eml[mul[x,x],add[x,x]],eml[neg[x],pow[x,x]]]`; predicted form in top 5 `False`, in top 20 `False`.

## Synthesis

None of the follow-ups produced clean Poisson KL support under the SR selection rule. The high-cosine expression confirms KL-like shape is expressible, but it is not selected and is not algebraically equivalent to KL.
Final verdict: `no_clean_poisson_support_in_sr_selection`.
