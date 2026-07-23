# Angle 1: Head-Window Held-out NLL

- Same 80/20 binomial type-count split as the Seam-Mandelbrot PMF runs.
- Held-out average NLL is restricted to train-rank windows `1-K` without refitting.
- Step-2 help counts are recomputed on the soft-k residual restricted to the same head window.

## top50

- soft-k beats MOE: `15` / 25
- soft-k beats ZM: `2` / 25
- soft-k beats Zipf: `17` / 25
- median soft-k minus MOE held-out avg NLL: `-0.010517713360`
- median soft-k minus ZM held-out avg NLL: `0.025947782146`
- winner counts: `{"zipf": 0, "zm": 23, "moe": 0, "softk": 2}`
- soft-k step-2 help count: `2` / 25

## top100

- soft-k beats MOE: `11` / 25
- soft-k beats ZM: `0` / 25
- soft-k beats Zipf: `21` / 25
- median soft-k minus MOE held-out avg NLL: `0.002818598264`
- median soft-k minus ZM held-out avg NLL: `0.022641031996`
- winner counts: `{"zipf": 0, "zm": 25, "moe": 0, "softk": 0}`
- soft-k step-2 help count: `6` / 25

## top200

- soft-k beats MOE: `7` / 25
- soft-k beats ZM: `3` / 25
- soft-k beats Zipf: `22` / 25
- median soft-k minus MOE held-out avg NLL: `0.008654673078`
- median soft-k minus ZM held-out avg NLL: `0.014337028806`
- winner counts: `{"zipf": 0, "zm": 15, "moe": 8, "softk": 2}`
- soft-k step-2 help count: `0` / 25

## top500

- soft-k beats MOE: `7` / 25
- soft-k beats ZM: `11` / 25
- soft-k beats Zipf: `24` / 25
- median soft-k minus MOE held-out avg NLL: `0.008138385900`
- median soft-k minus ZM held-out avg NLL: `0.004069542449`
- winner counts: `{"zipf": 0, "zm": 4, "moe": 14, "softk": 7}`
- soft-k step-2 help count: `1` / 25

## top1000

- soft-k beats MOE: `10` / 25
- soft-k beats ZM: `21` / 25
- soft-k beats Zipf: `25` / 25
- median soft-k minus MOE held-out avg NLL: `0.004404237552`
- median soft-k minus ZM held-out avg NLL: `-0.010661907584`
- winner counts: `{"zipf": 0, "zm": 1, "moe": 14, "softk": 10}`
- soft-k step-2 help count: `4` / 25

## full

- soft-k beats MOE: `17` / 25
- soft-k beats ZM: `13` / 25
- soft-k beats Zipf: `24` / 25
- median soft-k minus MOE held-out avg NLL: `-0.001866288732`
- median soft-k minus ZM held-out avg NLL: `-0.000584460494`
- winner counts: `{"zipf": 0, "zm": 11, "moe": 4, "softk": 10}`
- soft-k step-2 help count: `4` / 25

