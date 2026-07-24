"""f21 — the forecasting gauntlet: fit half the book, predict the whole curve.

f3/f3b/f19 proved lambda-ZM wins at FITTING. The practical question: does it
win at PREDICTING curves it hasn't seen? Protocol (chronological, realistic):

  1. tokenize in document order; split at 50% of tokens (prefix half);
  2. fit ZM (a,b,c), free lambda-ZM (a,b,c,lam), frozen lambda-ZM
     (a,b,c; lam=20.6) on the HALF curve;
  3. predict the FULL curve: model(r) + log(T_full/T_half) (token-mass
     rescale, same for all models; g's x frozen at fit-time V_half);
  4. score RMSE on log f over the full curve, split into head (r<=100),
     seen range (r<=V_half), and unseen tail (V_half < r <= V_full) — ranks
     the fragment never had.

Both models suffer the same depth systematics (c and lambda drift with
depth); the comparison is who lands closer AS USED, out of the box.

Outputs: ../outputs/f21_results.csv, f21_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

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


def load_tokens(fname):
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
    return TOKEN_RE.findall(text[start:end].lower())


def curve(tokens):
    return np.array(sorted(Counter(tokens).values(), reverse=True), dtype=np.float64)


def gx_of(ranks, V_norm):
    x = 0.05 + 0.95 * np.log(ranks) / math.log(V_norm)
    return np.exp(x - 1.0) - x


def c_grid_for(max_rank):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, max_rank, 2048)])


def fit(ranks, logf, gx=None, lam_fixed=None):
    """LSQ over the c-grid. gx None -> plain ZM. lam_fixed -> frozen term."""
    best = None
    for c in c_grid_for(ranks[-1]):
        lr = np.log(ranks + c)
        if gx is None:
            d = np.column_stack([np.ones_like(logf), lr])
            y = logf
        elif lam_fixed is not None:
            d = np.column_stack([np.ones_like(logf), lr])
            y = logf - lam_fixed * gx
        else:
            d = np.column_stack([np.ones_like(logf), lr, gx])
            y = logf
        coef, *_ = np.linalg.lstsq(d, y, rcond=None)
        mse = float(np.mean((d @ coef - y) ** 2))
        if best is None or mse < best[0]:
            lam = lam_fixed if lam_fixed is not None else (float(coef[2]) if gx is not None else 0.0)
            best = (mse, float(coef[0]), float(-coef[1]), float(c), lam)
    return {"a": best[1], "b": best[2], "c": best[3], "lam": best[4]}


def predict(m, ranks, V_norm, use_g):
    p = m["a"] - m["b"] * np.log(ranks + m["c"])
    if use_g:
        p = p + m["lam"] * gx_of(ranks, V_norm)
    return p


def main():
    rows = []
    for name, fname in CORPORA:
        toks = load_tokens(fname)
        half = toks[: len(toks) // 2]
        f_half, f_full = curve(half), curve(toks)
        Vh, Vf = len(f_half), len(f_full)
        Th, Tf = float(f_half.sum()), float(f_full.sum())
        shift = math.log(Tf / Th)
        rh = np.arange(1, Vh + 1, dtype=np.float64)
        lh = np.log(f_half)
        gh = gx_of(rh, Vh)
        models = {
            "zm": (fit(rh, lh), False),
            "lzm": (fit(rh, lh, gx=gh), True),
            "lzm_frozen": (fit(rh, lh, gx=gh, lam_fixed=LAM_STAR), True),
        }
        rf = np.arange(1, Vf + 1, dtype=np.float64)
        lf = np.log(f_full)
        row = {"corpus": name, "V_half": Vh, "V_full": Vf,
               "growth": round(Vf / Vh, 3)}
        for key, (m, use_g) in models.items():
            pred = predict(m, rf, Vh, use_g) + shift
            err = pred - lf
            row[f"{key}_full"] = round(float(np.sqrt(np.mean(err ** 2))), 4)
            row[f"{key}_head"] = round(float(np.sqrt(np.mean(err[:100] ** 2))), 4)
            row[f"{key}_seen"] = round(float(np.sqrt(np.mean(err[:Vh] ** 2))), 4)
            row[f"{key}_tail"] = round(float(np.sqrt(np.mean(err[Vh:] ** 2))), 4)
        rows.append(row)
        print(f"{name[:28]:29} full: zm={row['zm_full']:.4f} lzm={row['lzm_full']:.4f} "
              f"frozen={row['lzm_frozen_full']:.4f} | tail: zm={row['zm_tail']:.4f} "
              f"lzm={row['lzm_tail']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f21_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    lines = ["# f21 — extrapolation gauntlet (fit half, predict full)\n",
             "Chronological prefix split; token-mass rescale; x frozen at fit-time V.\n"]
    for band in ["full", "head", "seen", "tail"]:
        wz_free = sum(1 for r in rows if r[f"lzm_{band}"] < r[f"zm_{band}"])
        wz_froz = sum(1 for r in rows if r[f"lzm_frozen_{band}"] < r[f"zm_{band}"])
        med_impr = float(np.median([100 * (r[f"zm_{band}"] - r[f"lzm_{band}"]) / r[f"zm_{band}"] for r in rows]))
        lines.append(f"## {band}")
        lines.append(f"- free lambda-ZM beats ZM: {wz_free}/25 (median improvement {med_impr:+.1f}%)")
        lines.append(f"- frozen lambda*-ZM beats ZM: {wz_froz}/25")
    (OUT / "f21_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
