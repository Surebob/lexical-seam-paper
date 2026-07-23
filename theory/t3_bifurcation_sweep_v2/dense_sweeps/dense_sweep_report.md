# T3v2 Dense Sweep Report

Dense sweeps refine the two transition windows observed in the primary T3v2 break-rank sweep.

## Sweep A: K=1000..2000 step 100

- classification: `sharp_discrete_transition_width_le_100`
- IS wins: `10` / `11`
- exp wins: `1` / `11`
- c range: `255.749` to `412.463`
- winner changes: `[{'from_K': 1000, 'to_K': 1100, 'from_c': 255.74927297999844, 'to_c': 278.2559402207126, 'from_winner': 'eml[sub[x,1],eml[x,1]]', 'to_winner': 'sub[sub[x,1],log[x]]', 'from_gap': 0.0007220762391961649, 'to_gap': -0.00040239259270427596}]`
- gap sign changes: `[{'from_K': 1000, 'to_K': 1100, 'from_gap': 0.0007220762391961649, 'to_gap': -0.00040239259270427596, 'linear_interpolated_zero_K': 1064.2148736106628}]`

## Sweep B: K=2000..5000 step 300

- classification: `moderately_broad_transition_width_le_500`
- IS wins: `8` / `11`
- exp wins: `3` / `11`
- c range: `89.8674` to `419.479`
- winner changes: `[{'from_K': 4100, 'to_K': 4400, 'from_c': 189.84024180728048, 'to_c': 149.06722725319557, 'from_winner': 'sub[sub[x,1],log[x]]', 'to_winner': 'eml[sub[x,1],eml[x,1]]', 'from_gap': -0.0011067242375614045, 'to_gap': 8.653760630086826e-05}]`
- gap sign changes: `[{'from_K': 4100, 'to_K': 4400, 'from_gap': -0.0011067242375614045, 'to_gap': 8.653760630086826e-05, 'linear_interpolated_zero_K': 4378.243432467236}]`

## Interpretation

IS is not a single isolated point; it occupies a contiguous window in the dense sweeps.
At least one transition is sharp at the current grid resolution.
At least one transition spans multiple grid intervals, so the signed RMSE gap should be read as a continuous crossing rather than a literal discontinuity.
