# t2b_poisson_kl

- family: `Poisson`
- predicted generator: `generalized KL` (`x*log(x)-x+1`)
- search depth: `3`
- residual RMSE: `13.3333615996`
- residual lag-1 autocorrelation: `0.998256`
- residual second-difference sign-change fraction: `0.185185`
- predicted occurrences: `[{'step': 3, 'rmse_rank': 6436, 'expr': 'sub[add[1,neg[x]],mul[log[x],neg[x]]]', 'math': '((1+(-x))-(log(x)*(-x)))', 'rmse': 13.316046615320841, 'in_top5_by_rmse': False, 'in_diversity_beam': False, 'diversity_rank': None}]`
- global best: `{'step': 3, 'expr': 'div[pow[eml[x,1],eml[x,1]],pow[add[x,x],neg[x]]]', 'math': '(pow(EML(x,1),EML(x,1))/pow((x+x),(-x)))', 'rmse': 4.987523814020072, 'is_predicted_generator': False}`
- verdict: `FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5`

The predicted generator did not appear in the top-5 by RMSE or top-5 by diversity at its first reachable step under this density-estimation protocol.
