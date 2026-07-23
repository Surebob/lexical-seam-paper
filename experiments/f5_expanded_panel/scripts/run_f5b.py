"""F5b — erf-gate 9-parameter fits on the expanded panel (from f5a's processed
frequency vectors): does the s ≈ 1.2%·V law hold beyond classic English books?

Fits corpora: brown, cornell_dialogs, wikitext_1M, 6 languages, census_surnames
(control). 24 starts each, canonical bounds. Outputs s = k*w_tail per corpus and
the panel s-vs-V regression combined with the f2 English points.

Outputs: ../outputs/f5b_gate_fits.csv, f5b_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PROC = REPO / "data" / "processed_panel"
OUT = HERE.parent / "outputs"
F2 = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"

SEED = 20260722
N_STARTS = 24
MAX_NFEV = 9000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 5000.0, 10.0, 10.0])
# note: k upper raised to 5000 because panel corpora have V up to 162k


def predict(ranks, log_ranks, params):
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    z = (log_ranks - math.log(k)) / w_gate
    sigma = 0.5 * (1.0 - sp_erf(z))
    head = a1 - b1 * np.log(ranks + max(c1, 0.0))
    scale = max(1.0, k * w_tail)
    zz = np.clip((ranks - k) / scale, -60.0, 60.0)
    tail_rank = 1.0 + scale * np.log1p(np.exp(zz))
    tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
    return sigma * head + (1.0 - sigma) * tail


def piecewise_anchor(ranks, logf):
    k = min(500, len(ranks) // 3)

    def zmfit(rr, yy):
        best = None
        for c in np.concatenate([[0.0], np.geomspace(1e-6, rr[-1], 384)]):
            d = np.column_stack([np.ones_like(yy), np.log(rr + c)])
            coef, *_ = np.linalg.lstsq(d, yy, rcond=None)
            mse = float(np.mean((d @ coef - yy) ** 2))
            if best is None or mse < best[0]:
                best = (mse, float(coef[0]), float(-coef[1]), float(c))
        return best[1:]

    a1, b1, c1 = zmfit(ranks[:k], logf[:k])
    rr = np.arange(1, len(ranks) - k + 1, dtype=np.float64)
    a2, b2, c2 = zmfit(rr, logf[k:])
    return [a1, min(max(b1, 0.5), 3.0), min(c1, 1000.0), a2, min(max(b2, 0.5), 3.0), min(c2, 1000.0)]


def fit_task(task):
    name, freqs_list = task
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + hash(name) % (2**31))
    anchor6 = piecewise_anchor(ranks, y)

    def resid(p):
        return predict(ranks, log_ranks, p) - y

    best_rmse, best_p = None, None
    for i in range(N_STARTS):
        if i < N_STARTS // 2:
            x0 = np.array(anchor6 + [0.0, 0.0, 0.0])
            x0[:6] += rng.normal(0.0, [6.0, 0.3, 90.0, 6.0, 0.3, 90.0])
            x0[6] = min(0.02 * V, 4000.0) * math.exp(rng.normal(0.0, 0.7))
            x0[7] = 0.5 * math.exp(rng.normal(0.0, 0.45))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.45))
        else:
            x0 = np.array([
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(5000))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
            ])
        x0 = np.clip(x0, LOWER + 1e-8, UPPER - 1e-8)
        try:
            sol = least_squares(resid, x0=x0, bounds=(LOWER, UPPER), method="trf", max_nfev=MAX_NFEV)
        except Exception:
            continue
        r = float(np.sqrt(np.mean(sol.fun**2)))
        if best_rmse is None or r < best_rmse:
            best_rmse, best_p = r, sol.x
    return {"corpus": name, "V": V, "rmse": best_rmse,
            "k": float(best_p[6]), "w_gate": float(best_p[7]), "w_tail": float(best_p[8]),
            "s": float(best_p[6] * best_p[8]), "s_over_V": float(best_p[6] * best_p[8] / V),
            "k_hit_upper": bool(abs(best_p[6] - UPPER[6]) < 1e-3)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = []
    for path in sorted(PROC.glob("*.npy")):
        freqs = np.load(path)
        tasks.append((path.stem, freqs.tolist()))
    print(f"{len(tasks)} corpora", flush=True)
    rows = []
    with Pool(processes=10) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['corpus']:18} V={res['V']:7} s={res['s']:8.1f} "
                  f"s/V={res['s_over_V']:.4f} rmse={res['rmse']:.4f}", flush=True)

    with open(OUT / "f5b_gate_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    f2rows = list(csv.DictReader(open(F2, newline="", encoding="utf-8")))
    eng_V = np.array([float(r["V"]) for r in f2rows])
    eng_s = np.array([float(r["s_at_min"]) for r in f2rows])
    lang_rows = [r for r in rows if r["corpus"] != "census_surnames"]
    pan_V = np.array([float(r["V"]) for r in lang_rows])
    pan_s = np.array([float(r["s"]) for r in lang_rows])
    allV = np.concatenate([eng_V, pan_V])
    allS = np.concatenate([eng_s, pan_s])
    lV, lS = np.log(allV), np.log(allS)
    X = np.column_stack([np.ones_like(lV), lV])
    b, *_ = np.linalg.lstsq(X, lS, rcond=None)
    ssr = float(np.sum((lS - X @ b) ** 2))
    r2 = 1 - ssr / float(np.sum((lS - lS.mean()) ** 2))
    se = math.sqrt(ssr / (len(lV) - 2) / float(np.sum((lV - lV.mean()) ** 2)))

    lines = ["# F5b — s-law on the expanded panel\n"]
    lines.append("| corpus | V | s | s/V | rmse |")
    lines.append("|---|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: r["V"]):
        lines.append(f"| {r['corpus']} | {r['V']} | {r['s']:.1f} | {r['s_over_V']:.4f} | {r['rmse']:.4f} |")
    lines.append(f"\n- panel (language corpora) median s/V: {float(np.median(pan_s / pan_V)):.4f}")
    lines.append(f"- combined regression (25 English f2 + {len(lang_rows)} panel): "
                 f"beta={b[1]:.4f} 95% CI [{b[1]-1.96*se:.4f},{b[1]+1.96*se:.4f}] R2={r2:.4f}")
    (OUT / "f5b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
