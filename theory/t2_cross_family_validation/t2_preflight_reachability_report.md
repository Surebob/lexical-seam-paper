# T2 Pre-flight Reachability Check

Status: **STOP**

At least one predicted generator is not reachable at step <= 2 under the live Section 2.4 grammar.

Live Section 2.4 operator set used for this check:

- Unary: `['neg', 'inv', 'sqr', 'sqrt', 'exp', 'log']`
- Binary: `['eml', 'add', 'sub', 'mul', 'div', 'pow']`
- Note: the live implementation includes `eml` in addition to the manuscript base binary operators.

## Predicted Generator Reachability

### T2a_gaussian_euclidean

- family: `Gaussian`
- expected generator: `(1-x)^2 or (x-1)^2`
- reachable at step <= 2: `True`
- exact match: `mul[sub[1,x],sub[1,x]]` at step `2`; math `((1-x)*(1-x))`

### T2b_poisson_generalized_kl

- family: `Poisson`
- expected generator: `x*log(x)-x+1`
- reachable at step <= 2: `False`
- closest generated expression by RMSE: `add[neg[x],pow[x,x]]` at step `2`; RMSE to target generator `0.0403003`, max abs diff `0.06008`

### T2c_gamma_itakura_saito

- family: `Gamma`
- expected generator: `(x-1)-log(x)`
- reachable at step <= 2: `True`
- exact match: `sub[sub[x,1],log[x]]` at step `2`; math `((x-1)-log(x))`

## Stop Decision

The requested stop condition is triggered because the Poisson generalized-KL generator `x*log(x)-x+1` is not reachable at step 2 under the live Section 2.4 grammar. It requires composing `mul[x,log[x]]` with `sub[1,x]`, but `mul[x,log[x]]` is itself first generated at step 2, so the full KL expression would require a later step or a grammar macro.
