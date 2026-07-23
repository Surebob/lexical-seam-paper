# T2 Redesigned Report

Density-frame T2 completed for `3` families. Predicted generators passed the top-5 criterion in `0/3` cases.

## t2a_gaussian_euclidean

- family: `Gaussian`
- predicted generator: `Euclidean`
- verdict: `FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5`
- global best: `{'step': 2, 'expr': 'sub[sqrt[x],pow[x,x]]', 'math': '(sqrt(x)-pow(x,x))', 'rmse': 0.18740059320489277, 'is_predicted_generator': False}`
- predicted occurrences: `[{'step': 2, 'rmse_rank': 16, 'expr': 'mul[sub[1,x],sub[1,x]]', 'math': '((1-x)*(1-x))', 'rmse': 0.2184282539251113, 'in_top5_by_rmse': False, 'in_diversity_beam': False, 'diversity_rank': None}]`

## t2b_poisson_kl

- family: `Poisson`
- predicted generator: `generalized KL`
- verdict: `FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5`
- global best: `{'step': 3, 'expr': 'div[pow[eml[x,1],eml[x,1]],pow[add[x,x],neg[x]]]', 'math': '(pow(EML(x,1),EML(x,1))/pow((x+x),(-x)))', 'rmse': 4.987523814020072, 'is_predicted_generator': False}`
- predicted occurrences: `[{'step': 3, 'rmse_rank': 6436, 'expr': 'sub[add[1,neg[x]],mul[log[x],neg[x]]]', 'math': '((1+(-x))-(log(x)*(-x)))', 'rmse': 13.316046615320841, 'in_top5_by_rmse': False, 'in_diversity_beam': False, 'diversity_rank': None}]`

## t2c_gamma_is

- family: `Gamma`
- predicted generator: `Itakura-Saito`
- verdict: `FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5`
- global best: `{'step': 2, 'expr': 'add[log[x],sub[1,x]]', 'math': '(log(x)+(1-x))', 'rmse': 0.2512283114259147, 'is_predicted_generator': False}`
- predicted occurrences: `[{'step': 2, 'rmse_rank': 49, 'expr': 'sub[sub[x,1],log[x]]', 'math': '((x-1)-log(x))', 'rmse': 0.3401911416181829, 'in_top5_by_rmse': False, 'in_diversity_beam': True, 'diversity_rank': 49}]`

The redesigned density-frame test does not support the cross-family generator-identity claim as stated.
