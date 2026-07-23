# F4d — Poisson-lognormal mixture fits and re-tests

- two-component demanded by BIC: 25/25
- head share pi_H: median 2.15% (range 0.28-5.93%)

## P1 redo (corrected sd_H)
- corr(w_pred, w_fit) = 0.0223 (excl Dubliners: -0.1031)
- median ratio = 3.4709; sqrt(e)=1.6487
- corr(log ratio, sd_h^2/2) = -0.5311
- corr(log ratio, sd_t^2/2) = -0.3920
- corr(log ratio, sd_t/sd_h) = 0.3187

## Hunch-1 redo (s vs corrected head size)
- corr(log s, log N_H) = 0.6149
- median s/N_H = 0.5509 (range 0.187-3.285)
- log s ~ log N_H slope = 0.362