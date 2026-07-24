"""f23 — the completist bake-off: lambda-ZM vs every standard alternative family.

The paper beats ZM, MOEZipf, and hard/continuous piecewise. Referee-ammo
completeness: the other curve families people fit to rank-frequency data,
all on the same LSQ-on-log-frequency objective, BIC-compared (paper formula
p*log V + V*log rmse^2).

Families (per-corpus free params):
  zm          3  a - b log(r+c)
  lzm_free    4  + lam*g(x), lam free
  lzm_frozen  3  + 20.6*g(x)
  lognormal   3  quadratic in log r (lognormal-type rank curve)
  cubic       4  cubic in log r (flexible polynomial straw man)
  zm_cutoff   4  a - b log(r+c_zm) + d*r  (power law with exponential cutoff;
                 c inherited from the ZM fit — a small gift to this family)
  yule        2  a + lnGamma(r) - lnGamma(r+1+rho)  (Yule-Simon shape)
  hard_break  5  two independent power laws with a hard break at k (grid k)

Outputs: ../outputs/f23_results.csv, f23_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.special import gammaln

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
LAM_STAR = 20.6

CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"), ("War and Peace", "pg2600.txt"),
    ("Moby Dick", "pg2701.txt"), ("King James Bible", "pg10.txt"),
    ("Federalist Papers", "pg1404.txt"), ("Grimm's Fairy Tales", "pg2591.txt"),
    ("Don Quixote", "pg996.txt"), ("Pride and Prejudice", "pg1342.txt"),
    ("Canterbury Tales", "pg2383.txt"), ("Arabian Nights (Vol 1)", "pg3435.txt"),
    ("Aesop's Fables", "pg11339.txt"), ("Complete Sherlock Holmes", "pg1661.txt"),
    ("Jane Eyre", "pg1260.txt"), ("Dubliners", "pg2814.txt"),
    ("The Iliad", "pg6130.txt"), ("Democracy in America", "pg815.txt"),
    ("Origin of Species", "pg1228.txt"), ("Wealth of Nations", "pg3300.txt"),
    ("Les Miserables", "pg135.txt"), ("Decline and Fall Vol 1", "pg731.txt"),
    ("Emile", "pg5427.txt"), ("Ulysses", "pg4300.txt"),
    ("Collected Poe", "pg2147.txt"), ("Principia Ethica", "pg53430.txt"),
    ("Critique of Pure Reason", "pg4280.txt"),
]

PARAMS = {"zm": 3, "lzm_free": 4, "lzm_frozen": 3, "lognormal": 3, "cubic": 4,
          "zm_cutoff": 4, "yule": 2, "hard_break": 5}


def load_curve(fname):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    start, end = 0, len(text)
    for m in START_MARKERS:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            start = nl + 1 if nl != -1 else i
            break
    for m in END_MARKERS:
        i = text.find(m)
        if i != -1:
            end = i
            break
    counts = Counter(TOKEN_RE.findall(text[start:end].lower()))
    return np.array(sorted(counts.values(), reverse=True), dtype=np.float64)


def gx_of(ranks, V):
    x = 0.05 + 0.95 * np.log(ranks) / math.log(V)
    return np.exp(x - 1.0) - x


def c_grid_for(V):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, V, 2048)])


def lstsq_rmse(d, y):
    coef, *_ = np.linalg.lstsq(d, y, rcond=None)
    return float(np.sqrt(np.mean((d @ coef - y) ** 2))), coef


def fit_all(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    logf = np.log(freqs)
    logr = np.log(ranks)
    gx = gx_of(ranks, V)
    ones = np.ones_like(logf)
    out = {}

    best_zm, best_lf, best_lz, c_zm = None, None, None, 0.0
    for c in c_grid_for(V):
        lr = np.log(ranks + c)
        r2, _ = lstsq_rmse(np.column_stack([ones, lr]), logf)
        if best_zm is None or r2 < best_zm:
            best_zm, c_zm = r2, c
        r4, _ = lstsq_rmse(np.column_stack([ones, lr, gx]), logf)
        if best_lf is None or r4 < best_lf:
            best_lf = r4
        r3, _ = lstsq_rmse(np.column_stack([ones, lr]), logf - LAM_STAR * gx)
        if best_lz is None or r3 < best_lz:
            best_lz = r3
    out["zm"], out["lzm_free"], out["lzm_frozen"] = best_zm, best_lf, best_lz

    out["lognormal"], _ = lstsq_rmse(np.column_stack([ones, logr, logr ** 2]), logf)
    out["cubic"], _ = lstsq_rmse(np.column_stack([ones, logr, logr ** 2, logr ** 3]), logf)

    lr_zm = np.log(ranks + c_zm)
    out["zm_cutoff"], _ = lstsq_rmse(np.column_stack([ones, lr_zm, ranks]), logf)

    best_y = None
    for rho in np.geomspace(0.05, 10, 200):
        reg = gammaln(ranks) - gammaln(ranks + 1.0 + rho)
        a = float(np.mean(logf - reg))
        r = float(np.sqrt(np.mean((a + reg - logf) ** 2)))
        if best_y is None or r < best_y:
            best_y = r
    out["yule"] = best_y

    best_h = None
    for k in np.unique(np.geomspace(10, V / 2, 60).astype(int)):
        h = slice(0, k)
        t = slice(k, V)
        rh, _ = lstsq_rmse(np.column_stack([ones[h], logr[h]]), logf[h])
        rt, _ = lstsq_rmse(np.column_stack([ones[t], logr[t]]), logf[t])
        r = math.sqrt((rh ** 2 * k + rt ** 2 * (V - k)) / V)
        if best_h is None or r < best_h:
            best_h = r
    out["hard_break"] = best_h
    return out, V


def bic(V, p, rmse):
    return p * math.log(V) + V * math.log(max(rmse, 1e-12) ** 2)


def main():
    rows = []
    for name, fname in CORPORA:
        freqs = load_curve(fname)
        rmses, V = fit_all(freqs)
        bics = {k: bic(V, PARAMS[k], v) for k, v in rmses.items()}
        winner = min(bics, key=bics.get)
        row = {"corpus": name, "V": V, "bic_winner": winner}
        row.update({f"rmse_{k}": round(v, 4) for k, v in rmses.items()})
        rows.append(row)
        print(f"{name[:28]:29} winner={winner:11} lzm={rmses['lzm_free']:.4f} "
              f"logn={rmses['lognormal']:.4f} cub={rmses['cubic']:.4f} "
              f"cut={rmses['zm_cutoff']:.4f} yule={rmses['yule']:.4f} "
              f"hard={rmses['hard_break']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f23_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    winners = Counter(r["bic_winner"] for r in rows)
    lines = ["# f23 — family bake-off (BIC, paper formula)\n",
             f"- BIC winners: {dict(winners)}\n", "## median RMSE by family"]
    for k in PARAMS:
        med = float(np.median([r[f"rmse_{k}"] for r in rows]))
        beats = sum(1 for r in rows if r["rmse_lzm_free"] < r[f"rmse_{k}"]) if k != "lzm_free" else "-"
        lines.append(f"- {k} (p={PARAMS[k]}): median rmse {med:.4f}"
                     + (f"; lzm_free better on {beats}/25" if k != "lzm_free" else ""))
    (OUT / "f23_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
