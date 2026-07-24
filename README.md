# The Lexical Seam

**Public replication repository** for:

> Karapetyan, G. *The Lexical Seam: A Two-Population Account of Word
> Frequencies that Measures, Explains, and Predicts the Zipf–Mandelbrot Law.*
> Preprint, 2026. (`paper/MANUSCRIPT_v6.pdf`)

Count the words of any book, sort by frequency, and the Zipf–Mandelbrot law
appears — along with a systematic error nobody had diagnosed. This repository
contains every experiment behind the paper's account of that error: the
**lexical seam**, a measurable boundary between a small high-usage vocabulary
and the broad lexicon, whose rank-space width is a constant ≈1.2% of
vocabulary (s ≈ 0.012·V) across five centuries and four registers of English.

## What's here

- `paper/` — manuscript (Markdown source, LaTeX, PDF), figures, and the
  build pipeline (`scripts/md2tex.py`; compiles with
  [tectonic](https://tectonic-typesetting.github.io), not included).
- `experiments/` — 53 experiment directories, each with `scripts/`,
  `outputs/` (versioned CSVs), and a README stating its verdict. This
  includes the negative results: the paper's §4.2 documents eight retired
  claims, and the experiments that killed them are all here.
- `docs/MANUSCRIPT_v6_CLAIM_MAP.md` — **every quantitative claim in the
  paper mapped to its canonical CSV.** Start here to verify any number.
- `docs/PRIOR_ART_SWEEP.md` — the adversarial literature check each claim
  survived (with a dated correction the sweep itself later forced).
- `docs/v6_bootstrap_cis.md` — bootstrap confidence intervals for all
  headline numbers.
- `docs/lexical-seam-onepager.html` — a self-contained plain-language
  explainer (open in any browser).
- `data/` — the Project Gutenberg corpus panel (public domain), multilingual
  corpora, and 17k context-grounded word annotations. Large third-party
  corpora (Brown, WikiText-103, film dialogue, surname census) are fetched
  by scripts in `experiments/f5_expanded_panel/` rather than redistributed.
- `src/` — the shared fitting machinery (two-regime erf-gate model,
  Zipf–Mandelbrot baselines, corpus loaders).

## Reproducing

Python 3.12; `pip install -r requirements.txt`. Each experiment's README
states what it tests and where its outputs land; experiment scripts are
self-contained and deterministic (seeds in-file). The claim map is the
audit trail.

## Provenance

This work was carried out by an independent researcher in extensive
collaboration with Claude (Anthropic). Every claim traces to a reproducible
artifact in this repository; the git history documents the program,
including the claims that died along the way.
