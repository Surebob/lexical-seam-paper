# f4g — "two components demanded 25/25", re-audited (2026-07-21)

**Question:** does f4d's headline ("a 2-component latent mixture is demanded by BIC
on 25/25 corpora") survive when BOTH the 1-component and 2-component
Poisson-lognormal models get broad multi-start optimization? (f4f showed f4d's
restricted starts found a narrow basin.) 10 corpora × 8 diverse starts × 2 models.
`scripts/run_f4g.py`.

## Verdict: the pillar FALLS — 2-component wins only 3/10

| survives | falls |
|---|---|
| Shakespeare (ΔBIC 66), Moby Dick (38), Ulysses (0.1) | W&P, Bible, Federalist, P&P, Dubliners, Origin, Grimm (ΔBIC −10 to −23) |

A broad-tailed single lognormal explains the unordered count histogram about as
well as the two-component mixture on 7/10 corpora. f4d's 25/25 was a narrow-basin
artifact: restricted starts handicapped the 1-component model more than the
2-component one. Only the largest-V corpora retain histogram-level evidence for a
second component.

## Retired claim
- "Direct distributional two-population evidence, independent of rank curves"
  (f4d headline; also cited in the one-pager — both corrected this commit).

## What stands, and why this was always the weakest lens
The count histogram discards word identity and ordering; the ~300-type head club
is a tiny likelihood contribution among ~20k types. The two-vocabulary case rests
on identity- and structure-aware evidence, all untouched:
- function-word / content-word ablations (each population alone fits ZM cleanly;
  together they produce the seam) — canonical 2a;
- residual absorption asymmetry (smooth two-regime absorbs 25/25; MOEZipf leaves
  structure 19/25) — canonical 2c;
- the POS crossover, the erf-gate result (3e), the seam-width law (f2/f2b/f5b/f5c);
- the predicted c-collapse and out-of-sample V/b/c predictions (f6/f6b) — these are
  mixture-FAMILY observable results and never depended on the histogram demanding
  two components (note f6b's W&P fit predicted perfectly with π_H at a rail).

## v6 language
"The two-population structure is established by identity-aware and rank-structure
evidence; the count histogram alone resolves the second component only in the
largest corpora." The mechanism chapter loses an ornament and gains a spine.
