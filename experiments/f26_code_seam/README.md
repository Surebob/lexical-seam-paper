# f26 — the code seam: programming languages join the family line (2026-07-24)

**Question:** does source code have a lexical seam, and whose fingerprint —
language's (0.0118-0.013 by depth), or its own (like surnames 0.0266 /
belt 0.0166)? `scripts/run_f26.py`; token = identifier word
([A-Za-z_][A-Za-z0-9_]*, lowercased, no snake/camel splitting — stated
choice); fits share f24's machinery incl. the giant-regime bounds.

## Result: code sits ON the language width line — with double the amplitude

| corpus | tokens | V | tok/V | s/V | lambda | k |
|---|---:|---:|---:|---:|---:|---:|
| Python (site-packages, 13,004 files) | 22.4M | 513,736 | 43.5 | **0.0132** | **54.5** | 6,836 |
| English books at matched depth (f24 rung, tok/V 46) | 1.8M | 39,792 | 46.2 | **0.0132** | 26.3 | 842 |
| JS/TS projects (485 files) | 431k | 18,178 | 23.7 | 0.505 = ESCAPE | 27.8 | — |

- **Width: identical to three decimals at matched depth.** Python source is
  a member of the decaying-innovation family carrying the language width
  invariant — the strongest cross-system extension of the width law yet
  (fig7's line now stretches from Aesop to CPython's dependency graph).
- **Amplitude: 2x language's.** lambda ~ 54.5 vs ~26 at the same depth —
  code's usage club (keywords, stdlib, idiom identifiers) dominates the
  apex far harder than function words dominate prose. So lambda, which
  f20 showed CANNOT separate siblings within the family (GA vs classic
  Simon), CAN separate prose from code: the (width, amplitude) pair is a
  two-coordinate fingerprint at family distance.
- JS corpus: too small/mixed (430k tokens, 2 game projects) — the erf fit
  escapes (s/V 0.50); recorded as insufficient-data, not evidence.

## The bound catch (methods note, applies beyond f26)

First run returned k = 5000.000 exactly — the canonical fitter's k ceiling,
a book-regime search prior that binds on 500k-type corpora. Bounds now
scale with the corpus (k <= V/4, c <= V/10, b <= 4; `giant_bounds()` in
f24's script, inherited here). Corrected fit: k = 6,836, rmse improved
0.136 -> 0.122. The same fix protects the f24 giant ladder, whose k
trajectory (331 -> 842 -> 1373 by rung 100) was headed for the ceiling.

## Implications

- **Paper 2 (Atlas):** a new system ON the line, first two-coordinate
  fingerprint separation; other languages' code (Rust/Go corpora), split
  vs unsplit identifier tokenization as robustness.
- **Paper 3 (tokenizers):** code models face the same seam with a STRONGER
  club — deep-coverage vocabulary allocation should matter at least as
  much for code as for prose; code-arm experiment is a natural P3 rung.

Outputs: `outputs/f26_results.csv`, `f26_summary.md`.
