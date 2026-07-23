# f8 — matched-size cross-language panel (2026-07-21)

**Question:** what survives of the cross-language story when every corpus is measured
at the SAME sampling depth (first 65,000 tokens)? (f6 showed c is depth-dominated,
contaminating all unequal-size comparisons, including the v5.1 Table 3.)

21 corpora: 8 English classics + 13 other-language texts (Polish only 33k, flagged).
Per corpus at 65k: ZM (b, c), unit step-2 winner, PLN mixture signature.
`scripts/run_f8.py`; outputs in `outputs/`.

## Results

1. **The "c≈0 language family" claim mostly dies at matched depth.** English c
   collapses from [0.6–245] (full books) to [0–24] at 65k; non-English median 0.00
   (range 0–0.71). A small residual English-vs-other gap remains. Star exhibit —
   a natural experiment: **War and Peace at 65k in both languages: English
   translation c = 10.1, Russian original c = 0.0.** Same text, same depth; the
   difference isolates the language's morphology. (Don Quixote pair: EN 0.93 vs
   ES 0.02 — weaker but same direction.)

2. **The σ_T morphology dial survives and sharpens** (sorted, matched depth):
   ordinary English 1.68–1.76 < RU/PT/ES/SV 1.81–1.94 < DE/FR/LA/FI 2.0–2.08 <
   AR 2.37 < ZH 2.48. Within English, Shakespeare (1.94) and Ulysses (2.10) sit
   above ordinary prose — the dial flags exceptional vocabulary as well as
   morphology. (Mandarin value may be inflated by jieba single-character effects —
   tokenization caveat.)

3. **Winner ladder confirmed:** at 65k all corpora sit in the shallow regime
   (xpow 11/11 non-English, xpow/exp English) — consistent with f6.

4. **Caveat:** pi_H rails at its 0.1% floor for most non-English corpora at this
   depth — 65k under-resolves the head club; club-size comparisons need deeper
   matched slices (150k+ tier once deeper texts are fetched for all languages).

## Consequences for v6
- Replaces Table 3: all cross-corpus c claims must be depth-annotated.
- The multilingual section reframes from "low-c languages" to: (i) depth explains
  most of the old contrast; (ii) a small morphology-linked residual remains,
  cleanly shown by the translation pair; (iii) σ_T is the language-level dial.
