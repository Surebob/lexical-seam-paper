# Angle 3: Asymmetric Transition Gate

- The Seam-Mandelbrot PMF is extended with a two-sided logistic gate: separate left and right widths `w_left`, `w_right` around the seam.
- The rest of the PMF is unchanged, and the same soft-k lambda sweep is used to select the per-corpus prior strength.
- median width ratio w_right / w_left: `1.000000`
- mean width ratio w_right / w_left: `1.027741`
- corpora with material asymmetry (>1.5x width imbalance): `5` / 25

## top50

- asymmetric seam beats MOE: `15` / 25
- asymmetric seam beats ZM: `1` / 25
- asymmetric seam beats symmetric soft-k: `7` / 25
- asymmetric seam beats Zipf: `17` / 25
- median asym minus MOE held-out avg NLL: `-0.008994300711`
- median asym minus ZM held-out avg NLL: `0.025941468092`
- median asym minus symmetric soft-k held-out avg NLL: `0.000529892979`
- winner counts: `{"zipf": 0, "zm": 23, "moe": 0, "softk": 2, "asym": 0}`
- asymmetric seam step-2 help count: `2` / 25

## top100

- asymmetric seam beats MOE: `12` / 25
- asymmetric seam beats ZM: `0` / 25
- asymmetric seam beats symmetric soft-k: `12` / 25
- asymmetric seam beats Zipf: `21` / 25
- median asym minus MOE held-out avg NLL: `0.001539325557`
- median asym minus ZM held-out avg NLL: `0.022729503295`
- median asym minus symmetric soft-k held-out avg NLL: `0.000000000000`
- winner counts: `{"zipf": 0, "zm": 25, "moe": 0, "softk": 0, "asym": 0}`
- asymmetric seam step-2 help count: `6` / 25

## top200

- asymmetric seam beats MOE: `7` / 25
- asymmetric seam beats ZM: `3` / 25
- asymmetric seam beats symmetric soft-k: `17` / 25
- asymmetric seam beats Zipf: `22` / 25
- median asym minus MOE held-out avg NLL: `0.006623906581`
- median asym minus ZM held-out avg NLL: `0.014432677128`
- median asym minus symmetric soft-k held-out avg NLL: `-0.000515980872`
- winner counts: `{"zipf": 0, "zm": 15, "moe": 8, "softk": 1, "asym": 1}`
- asymmetric seam step-2 help count: `0` / 25

## top500

- asymmetric seam beats MOE: `7` / 25
- asymmetric seam beats ZM: `12` / 25
- asymmetric seam beats symmetric soft-k: `21` / 25
- asymmetric seam beats Zipf: `24` / 25
- median asym minus MOE held-out avg NLL: `0.006885325630`
- median asym minus ZM held-out avg NLL: `0.001627104010`
- median asym minus symmetric soft-k held-out avg NLL: `-0.000622377330`
- winner counts: `{"zipf": 0, "zm": 4, "moe": 14, "softk": 4, "asym": 3}`
- asymmetric seam step-2 help count: `1` / 25

## top1000

- asymmetric seam beats MOE: `10` / 25
- asymmetric seam beats ZM: `21` / 25
- asymmetric seam beats symmetric soft-k: `20` / 25
- asymmetric seam beats Zipf: `25` / 25
- median asym minus MOE held-out avg NLL: `0.001149107852`
- median asym minus ZM held-out avg NLL: `-0.014224408322`
- median asym minus symmetric soft-k held-out avg NLL: `-0.000589362922`
- winner counts: `{"zipf": 0, "zm": 1, "moe": 14, "softk": 3, "asym": 7}`
- asymmetric seam step-2 help count: `6` / 25

## full

- asymmetric seam beats MOE: `19` / 25
- asymmetric seam beats ZM: `13` / 25
- asymmetric seam beats symmetric soft-k: `22` / 25
- asymmetric seam beats Zipf: `24` / 25
- median asym minus MOE held-out avg NLL: `-0.002064939772`
- median asym minus ZM held-out avg NLL: `-0.000204221086`
- median asym minus symmetric soft-k held-out avg NLL: `-0.000272496315`
- winner counts: `{"zipf": 0, "zm": 11, "moe": 4, "softk": 2, "asym": 8}`
- asymmetric seam step-2 help count: `6` / 25

