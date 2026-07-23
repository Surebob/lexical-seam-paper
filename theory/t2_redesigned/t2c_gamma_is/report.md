# t2c_gamma_is

- family: `Gamma`
- predicted generator: `Itakura-Saito` (`(x-1)-log(x)`)
- search depth: `2`
- residual RMSE: `0.273726076375`
- residual lag-1 autocorrelation: `0.998435`
- residual second-difference sign-change fraction: `0.0167504`
- predicted occurrences: `[{'step': 2, 'rmse_rank': 49, 'expr': 'sub[sub[x,1],log[x]]', 'math': '((x-1)-log(x))', 'rmse': 0.3401911416181829, 'in_top5_by_rmse': False, 'in_diversity_beam': True, 'diversity_rank': 49}]`
- global best: `{'step': 2, 'expr': 'add[log[x],sub[1,x]]', 'math': '(log(x)+(1-x))', 'rmse': 0.2512283114259147, 'is_predicted_generator': False}`
- verdict: `FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5`

The predicted generator did not appear in the top-5 by RMSE or top-5 by diversity at its first reachable step under this density-estimation protocol.
