# t2a_gaussian_euclidean

- family: `Gaussian`
- predicted generator: `Euclidean` (`(1-x)^2`)
- search depth: `2`
- residual RMSE: `0.193482096611`
- residual lag-1 autocorrelation: `0.999778`
- residual second-difference sign-change fraction: `0.0134003`
- predicted occurrences: `[{'step': 2, 'rmse_rank': 16, 'expr': 'mul[sub[1,x],sub[1,x]]', 'math': '((1-x)*(1-x))', 'rmse': 0.2184282539251113, 'in_top5_by_rmse': False, 'in_diversity_beam': False, 'diversity_rank': None}]`
- global best: `{'step': 2, 'expr': 'sub[sqrt[x],pow[x,x]]', 'math': '(sqrt(x)-pow(x,x))', 'rmse': 0.18740059320489277, 'is_predicted_generator': False}`
- verdict: `FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5`

The predicted generator did not appear in the top-5 by RMSE or top-5 by diversity at its first reachable step under this density-estimation protocol.
