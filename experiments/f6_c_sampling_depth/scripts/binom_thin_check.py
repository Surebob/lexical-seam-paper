"""Binomial thinning vs prefix slicing: which does the PLN prediction correspond to?"""
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

DATA = Path(r"C:\Users\Greg Kara\Desktop\temporary\lexical-seam\data\zipf")
TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
CORPORA = [("Shakespeare", "pg100.txt"), ("War and Peace", "pg2600.txt"), ("Moby Dick", "pg2701.txt")]
SLICES = [300_000, 150_000, 80_000, 40_000]
rng = np.random.default_rng(20260722)


def counts_of(fname):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    st, e = 0, len(text)
    for m in SM:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            st = nl + 1 if nl != -1 else i
            break
    for m in EMK:
        i = text.find(m)
        if i != -1:
            e = i
            break
    return np.array(sorted(Counter(TOKEN_RE.findall(text[st:e].lower())).values(), reverse=True), dtype=np.int64)


def fit_c(freqs):
    freqs = np.asarray(freqs, dtype=float)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 768)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(c))
    return best[1]


print(f"{'corpus':16} {'T':>7} {'c_binom(med of 3)':>18}")
for name, fname in CORPORA:
    counts = counts_of(fname)
    T0 = int(counts.sum())
    for T in SLICES:
        if T > T0:
            continue
        cs = []
        for rep in range(3):
            thin = rng.binomial(counts, T / T0)
            thin = np.sort(thin[thin > 0])[::-1]
            cs.append(fit_c(thin))
        print(f"{name:16} {T:7} {np.median(cs):18.1f}   (range {min(cs):.1f}-{max(cs):.1f})")
