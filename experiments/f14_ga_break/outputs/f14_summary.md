# F14 — width behavior of the Gerlach-Altmann model

- params: Nc_max=7900, alpha=0.77 (gamma=1.77), pc0=0.99, pnew0=0.1; 2 runs, prefix snapshots.
| run | M | V | s | s/V |
|---:|---:|---:|---:|---:|
| 0 | 0.2M | 16480 | 185.3 | 0.01125 |
| 0 | 0.4M | 25855 | 308.6 | 0.01194 |
| 0 | 0.8M | 39323 | 480.8 | 0.01223 |
| 0 | 1.6M | 58918 | 730.4 | 0.0124 |
| 0 | 3.2M | 87830 | 1097.6 | 0.0125 |
| 1 | 0.2M | 16481 | 184.8 | 0.01121 |
| 1 | 0.4M | 25825 | 305.3 | 0.01182 |
| 1 | 0.8M | 39220 | 479.7 | 0.01223 |
| 1 | 1.6M | 58938 | 730.3 | 0.01239 |
| 1 | 3.2M | 87808 | 1096.5 | 0.01249 |

- GA model: log(s/V) ~ log V slope = +0.062 (real corpora: ~0.00 by the s-law; f12 pooled slope-1 = -0.023)
- GA s/V range across the sweep: 0.0112 - 0.0125 (x1.1); real corpora hold 0.009-0.013 across V 3.8k-63k

**Reading: the GA model reproduces an approximately constant s/V over this range — the width law does NOT discriminate against it here; weaken the manuscript's contrast accordingly.**