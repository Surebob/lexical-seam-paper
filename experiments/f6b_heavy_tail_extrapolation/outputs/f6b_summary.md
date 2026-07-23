# F6b — heavy-tail mixture and cross-size extrapolation

## Downward trajectories (median over reps)

| corpus | scale | V real | V pln | V nlt | c real | c pln | c nlt |
|---|---:|---:|---:|---:|---:|---:|---:|
| Complete Works of Shak | 1 | 24458 | 24480 | 24534 | 244.8 | 223.1 | 181.9 |
| Complete Works of Shak | 0.5 | 18836 | 18965 | 18856 | 103.7 | 116.4 | 67.0 |
| Complete Works of Shak | 0.25 | 14182 | 14260 | 14208 | 36.5 | 38.0 | 22.0 |
| Complete Works of Shak | 0.125 | 10302 | 10530 | 10426 | 11.9 | 8.3 | 10.4 |
| Complete Works of Shak | 0.0625 | 7426 | 7464 | 7252 | 3.4 | 3.3 | 3.4 |
| War and Peace | 1 | 17445 | 17502 | 17436 | 208.1 | 285.7 | 320.8 |
| War and Peace | 0.5 | 13607 | 13651 | 13662 | 74.8 | 128.5 | 118.8 |
| War and Peace | 0.25 | 10290 | 10237 | 10246 | 20.4 | 32.7 | 56.1 |
| War and Peace | 0.125 | 7476 | 7499 | 7514 | 5.0 | 12.9 | 21.8 |
| War and Peace | 0.0625 | 5302 | 5274 | 5280 | 0.9 | 4.1 | 6.7 |
| Moby Dick | 1 | 16956 | 17008 | 17072 | 10.7 | 14.1 | 12.2 |
| Moby Dick | 0.5 | 12268 | 12246 | 12312 | 1.3 | 2.2 | 1.7 |
| Moby Dick | 0.25 | 8513 | 8444 | 8522 | 0.0 | 0.5 | 0.0 |
| Moby Dick | 0.125 | 5624 | 5611 | 5702 | 0.0 | 0.0 | 0.0 |
| Moby Dick | 0.0625 | 3586 | 3608 | 3584 | 0.0 | 0.0 | 0.0 |
| Les Miserables | 1 | 22677 | 22488 | 22618 | 184.1 | 228.8 | 409.8 |
| Les Miserables | 0.5 | 17184 | 17184 | 17546 | 47.1 | 86.8 | 203.9 |
| Les Miserables | 0.25 | 12659 | 12514 | 12540 | 6.9 | 29.5 | 82.9 |
| Les Miserables | 0.125 | 9108 | 9040 | 8928 | 0.6 | 6.5 | 25.6 |
| Les Miserables | 0.0625 | 6170 | 6054 | 6139 | 0.0 | 2.1 | 1.5 |

## Upward prediction (fit on 150k slice, predict full corpus)

- Complete Works of Shakespeare: see console rows in run log / f6b_trajectories.csv (models *_upscaled_from_150k vs real full values)
- War and Peace: see console rows in run log / f6b_trajectories.csv (models *_upscaled_from_150k vs real full values)
- Moby Dick: see console rows in run log / f6b_trajectories.csv (models *_upscaled_from_150k vs real full values)
- Les Miserables: see console rows in run log / f6b_trajectories.csv (models *_upscaled_from_150k vs real full values)