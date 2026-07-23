# F3b held-out summary

Protocol: Binomial(n,1/2) thin split (seed 20260723), fit on one half (f3 objective, generator selected on train), evaluate fixed model on the other half's rank curve at ranks 1..min(V_A,V_B). Two folds per corpus.

## English 25 (50 fold-tests)
- lambda-ZM beats ZM on held-out RMSE: 50/50 (median improvement 8.5%)
- lambda-ZM beats MOEZipf-lsq on held-out RMSE: 50/50
- worst case: Pride and Prejudice B>A (+2.1%)

## Multilingual 7 (14 fold-tests)
- lambda-ZM beats ZM on held-out RMSE: 14/14 (median improvement 13.4%)
- lambda-ZM beats MOEZipf-lsq on held-out RMSE: 14/14
- worst case: Latin Gallic Wars B>A (+10.8%)
