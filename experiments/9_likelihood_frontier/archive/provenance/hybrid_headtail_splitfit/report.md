# Hybrid Head-Tail PMF

- Head component: normalized ZM law over ranks.
- Tail component: normalized MOEZipf law over ranks.
- The two are blended by a rank-dependent seam gate and then renormalized into one discrete PMF.
- The same soft-k lambda sweep is used to regularize seam location `k`.

## top50

- hybrid beats MOE: `15` / 25
- hybrid beats ZM: `0` / 25
- hybrid beats soft-k seam: `4` / 25
- hybrid beats Zipf: `16` / 25
- median hybrid minus MOE held-out avg NLL: `-0.002934838673`
- median hybrid minus ZM held-out avg NLL: `0.029812465488`
- median hybrid minus soft-k held-out avg NLL: `0.006710580172`
- winner counts: `{"zipf": 0, "zm": 23, "moe": 0, "softk": 2, "hybrid": 0}`
- hybrid step-2 help count: `0` / 25

## top100

- hybrid beats MOE: `13` / 25
- hybrid beats ZM: `0` / 25
- hybrid beats soft-k seam: `11` / 25
- hybrid beats Zipf: `22` / 25
- median hybrid minus MOE held-out avg NLL: `-0.002108418005`
- median hybrid minus ZM held-out avg NLL: `0.019677896763`
- median hybrid minus soft-k held-out avg NLL: `0.001103035332`
- winner counts: `{"zipf": 0, "zm": 25, "moe": 0, "softk": 0, "hybrid": 0}`
- hybrid step-2 help count: `0` / 25

## top200

- hybrid beats MOE: `10` / 25
- hybrid beats ZM: `5` / 25
- hybrid beats soft-k seam: `21` / 25
- hybrid beats Zipf: `23` / 25
- median hybrid minus MOE held-out avg NLL: `0.001740352172`
- median hybrid minus ZM held-out avg NLL: `0.012457448513`
- median hybrid minus soft-k held-out avg NLL: `-0.004083502788`
- winner counts: `{"zipf": 0, "zm": 14, "moe": 7, "softk": 1, "hybrid": 3}`
- hybrid step-2 help count: `0` / 25

## top500

- hybrid beats MOE: `9` / 25
- hybrid beats ZM: `17` / 25
- hybrid beats soft-k seam: `24` / 25
- hybrid beats Zipf: `24` / 25
- median hybrid minus MOE held-out avg NLL: `0.002133420227`
- median hybrid minus ZM held-out avg NLL: `-0.004547950531`
- median hybrid minus soft-k held-out avg NLL: `-0.005981033613`
- winner counts: `{"zipf": 0, "zm": 4, "moe": 13, "softk": 1, "hybrid": 7}`
- hybrid step-2 help count: `6` / 25

## top1000

- hybrid beats MOE: `11` / 25
- hybrid beats ZM: `23` / 25
- hybrid beats soft-k seam: `20` / 25
- hybrid beats Zipf: `25` / 25
- median hybrid minus MOE held-out avg NLL: `0.001427741186`
- median hybrid minus ZM held-out avg NLL: `-0.013584120913`
- median hybrid minus soft-k held-out avg NLL: `-0.002407572863`
- winner counts: `{"zipf": 0, "zm": 1, "moe": 13, "softk": 4, "hybrid": 7}`
- hybrid step-2 help count: `6` / 25

## full

- hybrid beats MOE: `21` / 25
- hybrid beats ZM: `15` / 25
- hybrid beats soft-k seam: `22` / 25
- hybrid beats Zipf: `25` / 25
- median hybrid minus MOE held-out avg NLL: `-0.003719569643`
- median hybrid minus ZM held-out avg NLL: `-0.001898584332`
- median hybrid minus soft-k held-out avg NLL: `-0.001677944397`
- winner counts: `{"zipf": 0, "zm": 9, "moe": 2, "softk": 1, "hybrid": 13}`
- hybrid step-2 help count: `9` / 25

