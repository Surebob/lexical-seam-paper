# f3 — lambda-ZM vs MOEZipf + head maps + multilingual (2026-07-21)

**Question:** is the 4-parameter lambda-ZM model (f1 finding) competitive against the
published MOEZipf alternative under a fair comparison, and does it generalize beyond
English?

`scripts/run_f3.py`; outputs in `outputs/`.

## Results

**English 25 (rank-curve objective, 3-way BIC: ZM p=3, lambda-ZM p=4, MOE p=2):**
- BIC winner: **lambda-ZM 25/25**.
- lambda-ZM rank-RMSE beats MOEZipf on 25/25 whether MOE is fit by MLE (canonical
  protocol) or by direct least squares on the same rank-curve objective
  (median RMSE delta vs MOE-lsq: −0.0157).
- On several corpora MOE-lsq is worse than plain ZM (e.g. Bible 0.1896 vs 0.1882) —
  the single-regime skew family cannot buy rank-curve accuracy even when optimized
  for it, consistent with the mechanism-absorption story (canon 2c/3.7).

**Head-only (top-100) free-amplitude winner map:** xpow 20 / exp 4 / is 1 —
vs full-curve free-amplitude map is 17 / xpow 7 / euclid 1 (f1). Generator identity
is a projection of the scoring functional; the seam and the lambda-ZM fit are the
scoring-stable objects.

**Multilingual (7/7 improve):**
| corpus | V | ZM c | ZM RMSE | lambda-ZM improvement |
|---|---:|---:|---:|---:|
| Russian War and Peace | 25,414 | 0.11 | 0.1971 | 12.9% |
| Mandarin Three Kingdoms | 37,049 | 0.00 | 0.2451 | 21.6% |
| Arabic 1001 Nights | 25,016 | 0.08 | 0.2047 | 14.3% |
| Latin Gallic Wars | 11,126 | 3.69 | 0.1879 | 7.4% |
| French Les Miserables | 14,406 | 1.51 | 0.1760 | 6.4% |
| Spanish Don Quixote | 13,114 | 1.98 | 0.1917 | 8.6% |
| Dutch Max Havelaar | 13,772 | 1.35 | 0.1856 | 8.1% |

## Caveats
- Multilingual vocabulary sizes differ from canonical Table 3 (e.g. Mandarin 37,049
  vs 23,483): this run used a 150k-token slice + current jieba for Mandarin where the
  canonical run used ~80k tokens and an older jieba; Unicode tokenizer details may
  also differ for French/Russian. Protocol variant, not a discrepancy — improvements
  are robust to it, but Table-3-facing numbers should use the canonical protocol.
- MOE rank-curve prediction uses the infinite-support PMF evaluated on ranks 1..V
  (canonical 3b convention); no truncation correction.
- BIC inherits the paper's formula (p·log n + n·log MSE) and its caveats.
