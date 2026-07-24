"""f21b — depth-AWARE forecasting: measure the dial inside the fragment,
then turn it forward.

f21 showed naive half->full prediction is mixed (15/25) because c and lambda
drift with depth. f21b uses ONLY information a forecaster would have — the
fragment itself: thin the half-corpus to a quarter (binomial), fit at both
depths, learn the per-corpus local drift per depth-doubling, then
extrapolate the parameters ONE doubling forward before predicting the full
curve. Both models get the same treatment (ZM extrapolates c; lambda-ZM
extrapolates c and lambda). Amplitude re-anchored by token mass as in f21.

If depth-aware lambda-ZM cleanly beats naive lambda-ZM and depth-aware ZM,
the paper's depth functions are quantitative forecasting tools, not just
descriptions.

Outputs: ../outputs/f21b_results.csv, f21b_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_f21 import (CORPORA, load_tokens, curve, gx_of, fit, predict)  # noqa: E402

OUT = HERE.parent / "outputs"
SEED = 20260736


def thin_counter(counts, rng):
    words = list(counts.keys())
    n = np.array([counts[w] for w in words], dtype=np.int64)
    a = rng.binomial(n, 0.5)
    return np.sort(a[a > 0])[::-1].astype(np.float64)


def fit_at(freqs, lam_free):
    V = len(freqs)
    r = np.arange(1, V + 1, dtype=np.float64)
    lf = np.log(freqs)
    g = gx_of(r, V)
    if lam_free:
        return fit(r, lf, gx=g), V
    return fit(r, lf), V


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for name, fname in CORPORA:
        toks = load_tokens(fname)
        half = toks[: len(toks) // 2]
        c_half = Counter(half)
        f_half = curve(half)
        f_qtr = thin_counter(c_half, rng)
        f_full = curve(toks)
        Vh, Vf = len(f_half), len(f_full)
        Th, Tf = float(f_half.sum()), float(f_full.sum())
        shift = math.log(Tf / Th)

        zm_h, _ = fit_at(f_half, False)
        zm_q, _ = fit_at(f_qtr, False)
        lz_h, _ = fit_at(f_half, True)
        lz_q, _ = fit_at(f_qtr, True)

        # one-doubling forward extrapolation of the drifting dials
        c_zm_pred = max(zm_h["c"] + (zm_h["c"] - zm_q["c"]), 0.0)
        c_lz_pred = max(lz_h["c"] + (lz_h["c"] - lz_q["c"]), 0.0)
        lam_pred = lz_h["lam"] + (lz_h["lam"] - lz_q["lam"])

        rf = np.arange(1, Vf + 1, dtype=np.float64)
        lf = np.log(f_full)

        def rmse_of(m, use_g):
            p = predict(m, rf, Vh, use_g) + shift
            return float(np.sqrt(np.mean((p - lf) ** 2)))

        row = {"corpus": name,
               "zm_naive": round(rmse_of(zm_h, False), 4),
               "zm_aware": round(rmse_of({**zm_h, "c": c_zm_pred}, False), 4),
               "lzm_naive": round(rmse_of(lz_h, True), 4),
               "lzm_aware": round(rmse_of({**lz_h, "c": c_lz_pred, "lam": lam_pred}, True), 4),
               "c_half": round(lz_h["c"], 1), "c_pred": round(c_lz_pred, 1),
               "lam_half": round(lz_h["lam"], 1), "lam_pred": round(lam_pred, 1)}
        rows.append(row)
        print(f"{name[:26]:27} lzm naive={row['lzm_naive']:.4f} aware={row['lzm_aware']:.4f} | "
              f"zm naive={row['zm_naive']:.4f} aware={row['zm_aware']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f21b_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    aware_beats_naive = sum(1 for r in rows if r["lzm_aware"] < r["lzm_naive"])
    aware_beats_zm_aware = sum(1 for r in rows if r["lzm_aware"] < r["zm_aware"])
    zm_aware_beats_naive = sum(1 for r in rows if r["zm_aware"] < r["zm_naive"])
    med = float(np.median([100 * (r["lzm_naive"] - r["lzm_aware"]) / r["lzm_naive"] for r in rows]))
    med2 = float(np.median([100 * (r["zm_aware"] - r["lzm_aware"]) / r["zm_aware"] for r in rows]))
    lines = ["# f21b — depth-aware forecasting (dial measured inside the fragment)\n",
             f"- depth-aware lambda-ZM beats naive lambda-ZM: {aware_beats_naive}/{n} "
             f"(median improvement {med:+.1f}%)",
             f"- depth-aware ZM beats naive ZM: {zm_aware_beats_naive}/{n}",
             f"- depth-aware lambda-ZM beats depth-aware ZM: {aware_beats_zm_aware}/{n} "
             f"(median {med2:+.1f}%)"]
    (OUT / "f21b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
