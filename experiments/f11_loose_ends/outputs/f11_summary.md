# F11 — loose ends

## A. Faithful 2b (deterministic legacy protocol)
- same_exponent_1.5_1.5: b=1.653 c=0.00 ZM rmse=0.1737 winner=is unit-step2 helps=True
- small_gap_1.5_1.3: b=1.426 c=0.00 ZM rmse=0.1361 winner=is unit-step2 helps=True
- pure_single_1.5: b=1.500 c=0.00 ZM rmse=0.0000 winner=exp unit-step2 helps=False

## B. Diachronic (first look)
- corr(year, s/V) all = -0.253
- corr(year, w_gate) all = -0.143
- corr(year, s/V) originals = -0.340
- corr(year, w_gate) originals = -0.193
- corr(year, hapax) originals = 0.267
  (n=25; originals n=15; years approximate, translations use translation year)

## C. Heaps beta: median 0.476 (range 0.391-0.626); corr(beta, 1/b) = 0.991

## D. Winner-from-stats classifier: LOO accuracy 0.800 (20/25) from (log tokens, hapax, TTR) — no curve fitting