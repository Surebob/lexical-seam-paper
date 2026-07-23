# Experiment B: Fully-Ours Variants

- B1 nested seam: early fixed mini-seam between two ZM components in the head, then a main seam into a tail ZM.
- B2 three-regime: three ZM components blended sequentially by two smooth gates.

## nested

- held-out wins vs MOE: `21` / 25
- held-out wins vs soft-k: `18` / 25
- held-out wins vs ZM: `16` / 25
- four-way winner counts: `{"zipf": 0, "zm": 9, "moe": 2, "nested": 14}`
- median minus MOE held-out avg NLL: `-0.002342149246`
- median minus ZM held-out avg NLL: `-0.002691216166`
- median minus soft-k held-out avg NLL: `-0.000850529644`
- step-2 help count: `11` / 25
- BIC wins vs two-regime seam: `0` / 25
- BIC wins vs hybrid head-tail: `13` / 25

## three_regime

- held-out wins vs MOE: `20` / 25
- held-out wins vs soft-k: `12` / 25
- held-out wins vs ZM: `15` / 25
- four-way winner counts: `{"zipf": 0, "zm": 9, "moe": 4, "three_regime": 12}`
- median minus MOE held-out avg NLL: `-0.001593793541`
- median minus ZM held-out avg NLL: `-0.001570814784`
- median minus soft-k held-out avg NLL: `0.000105268605`
- step-2 help count: `11` / 25
- BIC wins vs two-regime seam: `0` / 25
- BIC wins vs hybrid head-tail: `7` / 25
