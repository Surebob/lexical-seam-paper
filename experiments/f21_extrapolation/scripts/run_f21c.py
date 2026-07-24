"""f21c — depth-aware forecasting, corrected extrapolation rule (FINAL attempt).

f21b extrapolated c linearly per depth-doubling and failed (12/25). That
rule contradicts the paper's own section 3.4: c moves along log-log
trajectories, so the pre-justified rule is log-space extrapolation:
log(c+1) advances by its measured per-doubling step; lambda stays linear
per doubling (f20's curves are near-linear in log depth). Declared as the
single refinement — if this does not clearly win, formula-level depth
forecasting is recorded as open and parked for Paper 2.

Outputs: ../outputs/f21c_results.csv, f21c_summary.md
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
from run_f21 import CORPORA, load_tokens, curve, fit, gx_of, predict  # noqa: E402
from run_f21b import thin_counter, fit_at  # noqa: E402

OUT = HERE.parent / "outputs"
SEED = 20260737


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for name, fname in CORPORA:
        toks = load_tokens(fname)
        half = toks[: len(toks) // 2]
        c_half_counts = Counter(half)
        f_half = curve(half)
        f_qtr = thin_counter(c_half_counts, rng)
        f_full = curve(toks)
        Vh, Vf = len(f_half), len(f_full)
        Th, Tf = float(f_half.sum()), float(f_full.sum())
        shift = math.log(Tf / Th)

        zm_h, _ = fit_at(f_half, False)
        zm_q, _ = fit_at(f_qtr, False)
        lz_h, _ = fit_at(f_half, True)
        lz_q, _ = fit_at(f_qtr, True)

        def c_log_extrap(c_h, c_q):
            step = math.log(c_h + 1.0) - math.log(c_q + 1.0)
            return max(math.exp(math.log(c_h + 1.0) + step) - 1.0, 0.0)

        c_zm_pred = c_log_extrap(zm_h["c"], zm_q["c"])
        c_lz_pred = c_log_extrap(lz_h["c"], lz_q["c"])
        lam_pred = lz_h["lam"] + (lz_h["lam"] - lz_q["lam"])

        rf = np.arange(1, Vf + 1, dtype=np.float64)
        lf = np.log(f_full)

        def rmse_of(m, use_g):
            p = predict(m, rf, Vh, use_g) + shift
            return float(np.sqrt(np.mean((p - lf) ** 2)))

        row = {"corpus": name,
               "zm_naive": round(rmse_of(zm_h, False), 4),
               "zm_logc": round(rmse_of({**zm_h, "c": c_zm_pred}, False), 4),
               "lzm_naive": round(rmse_of(lz_h, True), 4),
               "lzm_logc": round(rmse_of({**lz_h, "c": c_lz_pred, "lam": lam_pred}, True), 4)}
        rows.append(row)
        print(f"{name[:26]:27} lzm naive={row['lzm_naive']:.4f} logc={row['lzm_logc']:.4f} | "
              f"zm naive={row['zm_naive']:.4f} logc={row['zm_logc']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f21c_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    lz_wins = sum(1 for r in rows if r["lzm_logc"] < r["lzm_naive"])
    zm_wins = sum(1 for r in rows if r["zm_logc"] < r["zm_naive"])
    cross = sum(1 for r in rows if r["lzm_logc"] < r["zm_logc"])
    med = float(np.median([100 * (r["lzm_naive"] - r["lzm_logc"]) / r["lzm_naive"] for r in rows]))
    lines = ["# f21c — depth-aware forecasting, log-space c extrapolation (final attempt)\n",
             f"- log-c lambda-ZM beats naive lambda-ZM: {lz_wins}/{n} (median {med:+.1f}%)",
             f"- log-c ZM beats naive ZM: {zm_wins}/{n}",
             f"- log-c lambda-ZM beats log-c ZM: {cross}/{n}"]
    (OUT / "f21c_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
