"""P1 check: does fitted w_gate ~ sqrt(2)*beta*sigma_H (derivation T3)?
beta = local |d log r / d log f| around the fitted transition centre k;
sigma_H = GMM head-component sd of type log-frequencies (from f4b_twins.csv).
w_gate, k from the canonical erf fits (archive s2_v3_per_fit_results.csv).
"""
import csv
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(r"C:\Users\Greg Kara\Desktop\temporary\lexical-seam")
DATA = REPO / "data" / "zipf"
ARCH = Path(r"C:\Users\Greg Kara\Desktop\temporary\emlexperiment\results\s2_v3_windows_full_outputs_2026-04-18\s2_v3_per_fit_results.csv")
TWINS = REPO / "experiments" / "f4b_calibrated_twins" / "outputs" / "f4b_twins.csv"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EM = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
FILES = {
    "Complete Works of Shakespeare": "pg100.txt", "War and Peace": "pg2600.txt",
    "King James Bible": "pg10.txt", "Federalist Papers": "pg1404.txt",
    "Origin of Species": "pg1228.txt", "Wealth of Nations": "pg3300.txt",
    "Moby Dick": "pg2701.txt", "Pride and Prejudice": "pg1342.txt",
    "Dubliners": "pg2814.txt", "Ulysses": "pg4300.txt",
    "Grimm's Fairy Tales": "pg2591.txt", "Don Quixote": "pg996.txt",
}

erf_fit = {r["corpus"]: r for r in csv.DictReader(open(ARCH, newline="", encoding="utf-8")) if r["gate"] == "erf"}
gmm = {r["corpus"]: r for r in csv.DictReader(open(TWINS, newline="", encoding="utf-8")) if r["family"] == "lognormal_twin"}

print(f"{'corpus':32} {'k':>6} {'beta':>6} {'sd_H':>6} {'w_pred':>7} {'w_fit':>6} {'ratio':>6}")
preds, fits = [], []
for name, fname in FILES.items():
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    s, e = 0, len(text)
    for m in SM:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            s = nl + 1 if nl != -1 else i
            break
    for m in EM:
        i = text.find(m)
        if i != -1:
            e = i
            break
    counts = Counter(TOKEN_RE.findall(text[s:e].lower()))
    freqs = np.array(sorted(counts.values(), reverse=True), dtype=float)
    ranks = np.arange(1, len(freqs) + 1, dtype=float)
    k = float(erf_fit[name]["k"])
    w_fit = float(erf_fit[name]["w_gate"])
    sd_h = float(gmm[name]["gmm_sd_head"])
    lo, hi = int(max(1, k / 2)), int(min(len(freqs), k * 2))
    lr, lf = np.log(ranks[lo - 1:hi]), np.log(freqs[lo - 1:hi])
    # robust local slope: fit log r on log f in the window
    A = np.column_stack([np.ones_like(lf), lf])
    coef, *_ = np.linalg.lstsq(A, lr, rcond=None)
    beta = abs(float(coef[1]))
    w_pred = math.sqrt(2.0) * beta * sd_h
    preds.append(w_pred)
    fits.append(w_fit)
    print(f"{name[:31]:32} {k:6.0f} {beta:6.2f} {sd_h:6.3f} {w_pred:7.3f} {w_fit:6.3f} {w_fit / w_pred:6.3f}")

preds, fits = np.array(preds), np.array(fits)
print(f"\ncorr(w_pred, w_fit) = {np.corrcoef(preds, fits)[0, 1]:.4f}")
print(f"median ratio w_fit/w_pred = {np.median(fits / preds):.3f}")
print(f"corr(log w_pred, log w_fit) = {np.corrcoef(np.log(preds), np.log(fits))[0, 1]:.4f}")
