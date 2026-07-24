# f19 — frozen-lambda (3-parameter) seam formula vs plain ZM

g fixed to exp(x-1)-x (25/25 winner in f3); lambda frozen by LOO median (English) / global English median 20.61 (transfer). Equal parameter count everywhere: 3 vs 3.

## English 25 (LOO-frozen)
- full-data RMSE wins vs ZM (equal params): 25/25 (median improvement 8.2%)
- held-out fold wins: 45/50
- median share of the free-lambda improvement retained: 95%
- worst case: Aesop's Fables (+0.9%)

## Non-English 7 (English-frozen transfer)
- full-data RMSE wins vs ZM (equal params): 7/7 (median improvement 8.6%)
- held-out fold wins: 14/14
- median share of the free-lambda improvement retained: 99%
- worst case: Latin Gallic Wars (+6.0%)
