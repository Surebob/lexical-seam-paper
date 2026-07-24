# f20 — lambda: mechanism-generated, depth-dependent, family-scale (2026-07-23)

Follow-up to f19's universal-amplitude result. Three hypotheses tested; one
confirmed, one refuted, one weakened. `scripts/run_f20.py`.

## A. Is lambda mechanism-generated? YES (at family scale)

Free-fit lambda on the f14b simulated growth processes:

| process | median lam | range |
|---|---:|---|
| GA reference (published constants) | 24.1 | [22.6, 24.8] |
| coreless decaying Simon | 24.0 | [22.2, 25.6] |
| classic Simon (constant innovation) | 26.4 | [23.9, 28.3] |
| language (EN 25 / canonical 7) | 20.6 / 21.9 | [14.2, 28.6] |

Decaying-innovation growth generates amplitudes in the language ballpark
(~24) with no tuning — the amplitude, like the width, is produced by the
growth family. BUT classic Simon overlaps GA/decaying (26.4 vs 24.1,
ranges touching), so **lambda is NOT a sharp mechanism discriminator** —
same lesson f14 taught about the width: shared invariant of the family,
not a separator within it. The "2-D fingerprint" hope is weakened for
within-family separation (far-family systems — surnames, belt — untested).

## B. 12-language panel concats

median 26.4, IQR [25.1, 26.8], range [24.0, 35.5]; lambda-ZM improves fit
11/11. Panel lambdas sit above the EN-book value despite mostly shallower
depths — cross-language offsets exist (morphology/concat-mixture effects,
as with s/V), so lambda is lambda(language, depth) with modest offsets,
not one number.

## C. Is lambda depth-invariant? NO — hypothesis refuted

Binomial thinning, 4 deep corpora, fractions 1..1/16: lambda rises smoothly
with tokens-per-type in every corpus (Shakespeare 16.1→23.5 over tok/V 8→40;
Bible 14.2→23.8 over 12→63; within-corpus curves near-monotone, pooled corr
+0.63). The pre-registered sleeper hypothesis ("lambda is the deeper,
depth-invariant constant") is dead.

## The honest synthesis — and it's prettier than the hypothesis

Lambda has exactly the same structure as the width constant: **the
book-depth value of a slowly varying universal depth function.** Both of
the paper's constants (s/V ≈ 0.0118, lambda* ≈ 20.6) are values of depth
functions in the regime where books live (tok/V ~15–40) — one sampling-depth
dial behind both. f19's practical result is untouched: within the tested
regime the drift is small enough that ONE frozen value beats ZM on 32/32
corpora and 59/64 held-out folds. A harmony clause was added to §3.7 so the
manuscript states this before any referee does.

Outputs: `outputs/f20_results.csv`, `f20_summary.md`. Paper-2 leads: lambda
depth function alongside F(gamma, depth); far-family lambda (surnames, belt);
analytic lambda-from-geometry functional.
