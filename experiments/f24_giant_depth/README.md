# f24 — the giant depth ladder: 45,937 books, 2.84 billion tokens (2026-07-24)

**Question:** every constant in the paper is the book-depth value of a rising
depth function, and our data ended at tokens-per-type ~40 (English) / ~30
(other languages). Where do the functions GO? SPGC (Gerlach & Font-Clos's
Standardized Project Gutenberg Corpus; data/ gitignored, fetch via Zenodo
record 2422561) merged in a seeded-shuffled ladder — valid as a depth
amplifier by the f12 composition result — from 10 books to all 45,937
English books. `scripts/run_f24.py`: per rung, single-ZM (c), free
lambda-ZM (lambda), canonical 9-param erf (s/V) with giant-regime bounds
(`giant_bounds()`) and quadrature-weighted log-grid fits above 200k types
(full-vs-sub validated identical to 4 figures at rung 100).

## The ladder

| books | tokens | V | tok/V | c | lambda | s/V |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 0.63M | 22,255 | 28 | 232 | 22.3 | 0.01254 |
| 30 | 1.8M | 39,792 | 46 | 740 | 26.3 | 0.01316 |
| 100 | 5.9M | 72,992 | 80 | 1,571 | 27.1 | 0.01290 |
| 300 | 17.9M | 138,114 | 130 | 1,760 | 26.6 | 0.01173* |
| 1,000 | 58M | 322,883 | 180 | 714 | 34.5 | 0.00962 |
| 3,000 | 181M | 622,673 | 291 | 443 | 45.0 | 0.00980 |
| 10,000 | 613M | 1,337,407 | 458 | 168 | 58.5 | 0.00980 |
| 30,000 | 1.86B | 2,641,125 | 704 | 62 | 66.0 | 0.01019 |
| ALL | 2.84B | 3,491,525 | 815 | 38 | 69.5 | 0.01064 |

*rung 300: w_tail pinned at its bound (10.0) — the detectable escape
signature; treat its s/V as unreliable (flagged, not hidden).

## Headline: the width constant SURVIVES three billion tokens

Across a 150x growth in vocabulary (22k -> 3.5M types) and a 29x growth in
sampling depth, **s/V never leaves the band [0.0096, 0.0132]** — i.e. within
−19%/+12% of the book value 0.0118 — while every other quantity moves by one
to two orders of magnitude (c by 46x, k across basins by 500x, lambda by
3x). No runaway, no collapse: the 1.2%-scale width law holds from a
ten-book shelf to all of English Project Gutenberg. Fine structure: a hump
(~0.013 at tok/V 46–80), a shallow dip (~0.0096–0.0098 at 180–460), a gentle
rise to 0.0106 at the summit.

## Second finding: lambda does NOT saturate

lambda climbs smoothly 22 -> 70 across the ladder with no plateau — the
amplitude's depth function is strong, monotone, and still rising at 2.8B
tokens. This confirms f20/§3.7's scoping (lambda* = 20.6 is the book-depth
value) and hands Paper 2 a beautiful measured curve. Note the f26
code-vs-language amplitude contrast (55 vs 26) was at MATCHED depth and is
unaffected.

## Third finding: c turns over at aggregation scale

c rises through the book regime (232 -> 1,760 by 18M tokens) then COLLAPSES
back to 38 at 2.8B — while b falls 2.07 -> 1.62. The paper's §3.4 claim
(within-corpus thinning collapses c predictably) is untouched — this is a
DIFFERENT protocol (cross-scale aggregation), and the turnover is new,
undocumented structure: as thousands of heterogeneous books average out,
the apex flattening shrinks again. Open mechanism question for Paper 2.

## Caveats (read before quoting)

1. The ladder confounds depth with COMPOSITION: high rungs are massive
   heterogeneous mixtures (all genres/eras + OCR junk + proper-name floods),
   which inflate V mechanically and could deflate s/V — part of the dip may
   be junk-vocabulary inflation, not depth physics. A cleaned-vocabulary
   rerun (frequency floor / dictionary filter) is the natural control.
2. Sub-grid fits validated against full fits at V=73k only; validity at
   V=3.5M is extrapolated (quadrature-weighted, but unverified at scale).
3. One run, one shuffle order; no seed replication yet.
4. k/w basins switch across rungs (k: 331 -> 84k -> 34k) — s = k*w_tail is
   the stable observable, as everywhere else in the program.
5. Paper-1 claims are UNTOUCHED by design: the manuscript's depth statements
   are about the f15d within-corpus protocol and the cross-language panel at
   its depths; f24 is Paper-2 material (its opening figure).

Outputs: `outputs/f24_ladder.csv`, `f24_summary.md`.
