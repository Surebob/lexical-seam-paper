# Angle 4: Empirical-Bayes Hierarchical k Pooling

- Each corpus gets its own Seam-Mandelbrot PMF, but `log k_i` is pooled through a shared Gaussian prior.
- The prior center and spread are fit iteratively across corpora.
- final alpha: `0.504917434096`
- final sigma: `0.050000000000`

- iteration 1: alpha `0.504917386870`, sigma `0.050000000000`
- iteration 2: alpha `0.504917373901`, sigma `0.050000000000`
- iteration 3: alpha `0.504917380785`, sigma `0.050000000000`
- iteration 4: alpha `0.504917434096`, sigma `0.050000000000`

## top50

- hierarchical seam beats MOE: `15` / 25
- hierarchical seam beats ZM: `2` / 25
- hierarchical seam beats Zipf: `17` / 25
- median hier minus MOE held-out avg NLL: `-0.010704405594`
- median hier minus ZM held-out avg NLL: `0.021000591305`
- winner counts: `{"zipf": 0, "zm": 23, "moe": 0, "hier": 2}`
- hierarchical seam step-2 help count: `0` / 25

## top100

- hierarchical seam beats MOE: `11` / 25
- hierarchical seam beats ZM: `0` / 25
- hierarchical seam beats Zipf: `21` / 25
- median hier minus MOE held-out avg NLL: `0.007988222185`
- median hier minus ZM held-out avg NLL: `0.024910904508`
- winner counts: `{"zipf": 0, "zm": 25, "moe": 0, "hier": 0}`
- hierarchical seam step-2 help count: `8` / 25

## top200

- hierarchical seam beats MOE: `6` / 25
- hierarchical seam beats ZM: `0` / 25
- hierarchical seam beats Zipf: `22` / 25
- median hier minus MOE held-out avg NLL: `0.015036929433`
- median hier minus ZM held-out avg NLL: `0.018422752525`
- winner counts: `{"zipf": 0, "zm": 15, "moe": 10, "hier": 0}`
- hierarchical seam step-2 help count: `4` / 25

## top500

- hierarchical seam beats MOE: `4` / 25
- hierarchical seam beats ZM: `11` / 25
- hierarchical seam beats Zipf: `24` / 25
- median hier minus MOE held-out avg NLL: `0.011433201139`
- median hier minus ZM held-out avg NLL: `0.004793398880`
- winner counts: `{"zipf": 0, "zm": 4, "moe": 17, "hier": 4}`
- hierarchical seam step-2 help count: `1` / 25

## top1000

- hierarchical seam beats MOE: `9` / 25
- hierarchical seam beats ZM: `20` / 25
- hierarchical seam beats Zipf: `25` / 25
- median hier minus MOE held-out avg NLL: `0.005171093906`
- median hier minus ZM held-out avg NLL: `-0.007428225048`
- winner counts: `{"zipf": 0, "zm": 1, "moe": 15, "hier": 9}`
- hierarchical seam step-2 help count: `7` / 25

## full

- hierarchical seam beats MOE: `16` / 25
- hierarchical seam beats ZM: `16` / 25
- hierarchical seam beats Zipf: `24` / 25
- median hier minus MOE held-out avg NLL: `-0.001515605083`
- median hier minus ZM held-out avg NLL: `-0.000723616273`
- winner counts: `{"zipf": 0, "zm": 9, "moe": 5, "hier": 11}`
- hierarchical seam step-2 help count: `4` / 25

