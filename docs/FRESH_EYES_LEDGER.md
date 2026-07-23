# Fresh-Eyes Ledger — 2026-07-20

Result of a skeptical re-mining of the FULL archive (all ~196 historical bundles,
including every run the manuscript never cited, plus recovered audit trail), done by
Claude (Fable) explicitly NOT trusting the April models' distillation. Three categories:
**(A)** findings that were walked back or never published and deserve resurrection,
**(B)** published conclusions I do not yet trust and will re-litigate, **(C)** genuinely
unexplored branches. Archive paths refer to the frozen `emlexperiment` repo.

---

## A. Walked-back / never-published findings worth resurrecting

**A1. The λ-ZM four-parameter formula.** `results/zipf_lambda_test`: fit
`log f = log a − b·log(r+c) + λ·[(x−1)−log x]` by frequency-weighted LS. Improves BOTH
weighted and unweighted RMSE on all 4 corpora tested (λ = 0.22 Shakespeare, 0.65 W&P,
0.70 Moby, 0.95 Bible). This is "the better Zipf formula" in its most parsimonious,
drop-in form — one interpretable seam-strength parameter — and it appears in NO
manuscript version. Run on all 25 English + 7 non-English; BIC vs MOEZipf.
Bonus puzzle: λ is NOT monotone in c (Shakespeare has the highest c but lowest λ) —
λ and c measure different aspects of the seam. Unexplored.

**A2. The exponent-gap organizer.** `results/zipf_exponent_gap`: fit separate power
laws to the function-word proxy (top-100) and re-ranked content population. Then
corr(c, IS-winner) = **0.883**, corr(gap, c) = −0.689. The mysterious "c ≈ 66–79
threshold" may reduce to the interpretable statement *the step-2 winner is organized by
the exponent gap between the two lexical populations*. Feeds threshold theory (T3v2)
and the erf derivation (C1). Never cited.

**A3. Lexical-statistics organizer / the "what is c?" question.**
`results/zipf_lexical_metrics`: hapax_ratio correlates with winner at r = −0.65, TTR
at −0.58; tokens correlate with c at +0.72. Two consequences: (i) a reviewer-killer
confound the paper never addresses — is c partly a corpus-length/richness artifact?
(Ulysses: 269k tokens but c = 0.63, so not purely — but v6 must treat this head-on);
(ii) an opportunity — predict the winner from raw lexical stats without fitting ZM at
all. Never cited.

**A4. Universal high-c head shape (zero-shot transfer).** `results/zipf_beta_kernel`:
a Beta-kernel residual model fitted ONLY on Shakespeare transfers zero-shot to the
other 6 high-c corpora at essentially in-domain accuracy (e.g. Federalist 0.17031 vs
in-domain 0.17007). The high-c head residual appears to be ONE universal curve in
normalized coordinates. Stronger than the published poly-transfer result. Caveat: the
Beta amplitude slammed its −10⁴ bound (broken parameterization — refit in log-amplitude
before trusting shapes). Never cited.

**A5. Chronological stability of the seam.** `results/zipf_chronological_split`:
freeze smooth-model params on the first half of a text, evaluate on the second half —
transfer gap ~1% on single-author corpora (Moby 0.9%, Shakespeare 1.1%) vs 80% on
Federalist (multi-author). The seam is a stable property of coherent authorship, and
transfer gap doubles as a corpus-coherence diagnostic (independent support for the
Bible/anthology scope condition). Never cited.

**A6. SPARC/RAR: an unused physics negative control.** `results/sparc_rar_eml*`: the
identical pipeline on 2,693 galaxy rotation points found NO accepted correction over
the McGaugh RAR baseline. Alongside city populations, this is a second, harder
specificity control ("the method finds nothing on a tight physical law"). Never cited.

**A7. The identity of the MOEZipf residual.** `results/zipf_breakthrough_probe` §3:
after MOEZipf subtraction the step-2 winner is consistently `(1−x)·log x` (11 corpora)
or `±[log x + (1−x)]` — the latter is the NEGATIVE of the IS generator. MOEZipf
over-bends the apex, leaving anti-IS residue. Sharp mechanistic sentence for v6 §3.7;
currently the paper only says "Bregman-shaped structure remains."

**A8. Hierarchical-k pooling is a POSITIVE result.** Experiment 7e data: empirical-
Bayes pooling of k beats canonical soft-k on 20/25 (median −0.000912). v5.1 asserts
the opposite (audit caught it). For Paper 2: the better regularizer is already built
and sitting unused.

**A9. log(x)² is a third valid Bregman generator** (`results/zipf_literature_review`
derivative checks) and shows up in MOE-residual searches. Worth one sentence when
describing the beam family.

---

## B. Published conclusions to re-litigate (do not trust yet)

**B1. "k ~ √V weakens under the decoupled model." — RESOLVED 2026-07-21 (f2), with a
twist.** Profile likelihood shows k IS identified within corpus (median ΔBIC≤2 interval
±1%; the flat-ridge hypothesis is refuted), and k itself genuinely has no tight V-law
(β=0.758, R²=0.386 — the April walk-back stands for k). BUT the tight law exists in
the product: **s = k·w_tail ∝ V^1.003, CI [0.93, 1.08], R² = 0.967** — the tail
crossover scale is ≈1.2% of vocabulary across all 25 corpora. The coupled-model
k~V^0.52 tightness was the single-width parameterization compressing this linear law
into the wrong parameter. See `experiments/f2_k_profile_likelihood/`. Follow-up: direct
s-profile before publishing; interpret k's residual variation via A2/A3 organizers.

**B2. erf-gate specificity magnitude.** Direction is well-supported (3e + logistic-
recovery falsification), but a median BIC spread of 677 from gate shape alone is
suspiciously enormous, and the falsification only tested logistic-generated truth.
Missing: erf-generated recovery, misspecified-truth battery, optimizer-convergence
asymmetry audit (= old EXP08 idea). Also only regime coefficients for erf were saved.

**B3. The c-threshold framing.** T3v2 showed winner flips are smooth RMSE-gap
crossings, and A2/A3 suggest c is a proxy for (exponent gap, hapax structure). v6
should state the organizer in those terms, not as a bare c-band.

**B4. Low-c span R² inconsistency.** 5b reports head-200 span R² median 0.836 (Aesop
0.49!) while 10b's widened-grammar table reports 17/17 at R² > 0.975. Both are true
under different weighting/protocols — but v6 must reconcile them explicitly or a
referee will.

**B5. Poisson noise in simulation recovery.** Real corpora are bursty (word
clustering); Poisson understates count dispersion. Rerun recovery with negative-
binomial noise to check the 58.8%/0% contrast survives.

**B6. All BIC counts** inherit n = V and no autocorrelation correction (audit noted).
Affects smooth-vs-MOE comparisons everywhere; consider row-subsampled or blocked
robustness checks.

**B7. Single non-linguistic control.** City populations (+ unused SPARC, A6) is thin
for the specificity claim; add 2–3 more Zipfian systems (surnames, firm sizes, city
data by country).

---

## C. Genuinely unexplored branches

**C1. Derive the erf gate.** If the two lexical populations have (approximately)
log-normal frequency distributions, the mixture fraction as a function of log-rank is
an erf by construction — turning 3e's empirical 24/25 into a predicted result. The
A2 (two-population exponents) and A3 (hapax/TTR) data are exactly the ingredients.
Highest theory value per unit effort.

**C2. Diachronic seam drift.** The 25 corpora span six centuries; seam parameters vs
composition date has never been plotted. A5's stability result makes this well-posed.

**C3. Heaps' law connection.** k ~ V^α is a Heaps-adjacent statement and hapax/TTR
correlations (A3) hint the seam and vocabulary growth are coupled. No formal link
attempted anywhere in the archive.

**C4. Winner-from-lexical-stats classifier.** Predict IS-vs-exp from (tokens, TTR,
hapax) without any ZM fit (from A3). Cheap, interpretable, and a strong robustness
exhibit.

**C5. Free the mini-seam.** Experiment 9's nested-seam fixed k₀, w₀ inside the head;
freeing them tests whether the head itself is multi-regime (n-regime generalization).

**C6. Cross-linguistic exponent gap.** Does A2's organizer hold on the 7 non-English
corpora? (Data already in `data/zipf_multilang*`.)

**C7. λ-vs-c decoupling** (from A1): what does λ measure that c doesn't? Likely seam
*amplitude* vs seam *position* — connects to the analytic seam expansion (a₁ coefficient).

---

## Source-coverage note

The Codex transcripts contain only run-execution sessions (verified: every other
"zipf" hit was `zipfile` noise). The original ideation conversations (Opus 4.7 "lab
notebook", GPT-5 Pro, Gemini) live in claude.ai / ChatGPT accounts, not on disk. The
orphan bundles above are the substantive record. Optional future step: export those
account conversations into `docs/archive_maps/` for completeness.
