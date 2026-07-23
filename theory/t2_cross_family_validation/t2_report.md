# T2 Cross-Family Synthetic Validation

Depth-4 search completed for `3` synthetic families. Predicted generator wins: `0/3`.

## T2a Gaussian -> Euclidean

- predicted: `Euclidean`
- observed global best: `add[mul[div[1,eml[x,x]],log[eml[x,x]]],div[div[log[x],pow[x,x]],div[sqrt[x],add[x,x]]]]`
- predicted first generated: `{'step': 2, 'beam_rank': 66, 'expr': 'mul[sub[1,x],sub[1,x]]', 'math': '((1-x)*(1-x))', 'rmse': 1.3799198797255174}`
- predicted first retained top-50: `None`
- verdict: `FAIL_PREDICTED_GENERATED_BUT_RANKED_OUTSIDE_TOP50`

## T2b Poisson -> generalized KL

- predicted: `generalized KL`
- observed global best: `add[div[sub[mul[x,x],pow[x,x]],mul[pow[x,x],pow[x,x]]],sub[div[sqrt[x],eml[x,x]],add[pow[x,x],neg[x]]]]`
- predicted first generated: `None`
- predicted first retained top-50: `None`
- verdict: `FAIL_PREDICTED_NOT_IN_BEAM`

## T2c Gamma -> Itakura-Saito

- predicted: `Itakura-Saito`
- observed global best: `div[sub[mul[sub[1,x],eml[1,x]],eml[log[x],add[x,x]]],eml[add[neg[x],neg[x]],log[eml[x,x]]]]`
- predicted first generated: `{'step': 2, 'beam_rank': 53, 'expr': 'sub[sub[x,1],log[x]]', 'math': '((x-1)-log(x))', 'rmse': 1.4640393720385372}`
- predicted first retained top-50: `None`
- verdict: `FAIL_PREDICTED_GENERATED_BUT_RANKED_OUTSIDE_TOP50`

Zero or one synthetic family selected the predicted generator. This does not support the cross-family generalization as stated under this protocol.
