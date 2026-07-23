"""F5c — refit the three degenerate f5b gate fits (Portuguese, Swedish, Polish)
with the escape hatch closed: w_gate lower bound raised to 0.15, w_tail upper
bound tightened to 2.5, 32 starts. Also refit Finnish/German/Italian under the
same bounds as a consistency check.

Outputs: ../outputs/f5c_refits.csv, f5c_summary.md
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
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260726
N_STARTS = 32
MAX_NFEV = 9000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.15, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 3000.0, 10.0, 2.5])

TARGETS = ["lang_portuguese", "lang_swedish", "lang_polish", "lang_finnish", "lang_german", "lang_italian"]


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


def fit_task(job):
    name, freqs_list = job
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
            x0[6] = min(0.02 * V, 2000.0) * math.exp(rng.normal(0.0, 0.7))
            x0[7] = 1.0 * math.exp(rng.normal(0.0, 0.5))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.5))
        else:
            x0 = np.array([
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(3000))),
                math.exp(rng.uniform(math.log(0.15), math.log(10))),
                math.exp(rng.uniform(math.log(0.05), math.log(2.5))),
            ])
        x0 = np.clip(x0, LOWER + 1e-8, UPPER - 1e-8)
        try:
            sol = least_squares(resid, x0=x0, bounds=(LOWER, UPPER), method="trf", max_nfev=MAX_NFEV)
        except Exception:
            continue
        r = float(np.sqrt(np.mean(sol.fun**2)))
        if best_rmse is None or r < best_rmse:
            best_rmse, best_p = r, sol.x
    hits = {"w_gate_lo": abs(best_p[7] - LOWER[7]) < 1e-3, "w_tail_hi": abs(best_p[8] - UPPER[8]) < 1e-3}
    return {"corpus": name, "V": V, "rmse": best_rmse, "k": float(best_p[6]),
            "w_gate": float(best_p[7]), "w_tail": float(best_p[8]),
            "s": float(best_p[6] * best_p[8]), "s_over_V": float(best_p[6] * best_p[8] / V),
            "bound_hit": any(hits.values())}


def main():
    jobs = []
    for t in TARGETS:
        p = PROC / f"{t}.npy"
        if p.exists():
            jobs.append((t, np.load(p).tolist()))
    rows = []
    with Pool(processes=6) as pool:
        for res in pool.imap_unordered(fit_task, jobs, chunksize=1):
            rows.append(res)
            print(f"{res['corpus']:18} V={res['V']:6} s/V={res['s_over_V']:.4f} "
                  f"w_gate={res['w_gate']:.3f} w_tail={res['w_tail']:.3f} rmse={res['rmse']:.4f} "
                  f"bound_hit={res['bound_hit']}", flush=True)
    with open(OUT / "f5c_refits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    lines = ["# F5c — degenerate refits (tightened bounds, 32 starts)\n"]
    for r in sorted(rows, key=lambda r: r["corpus"]):
        lines.append(f"- {r['corpus']}: s/V={r['s_over_V']:.4f} (w_gate {r['w_gate']:.3f}, "
                     f"w_tail {r['w_tail']:.3f}, bound_hit={r['bound_hit']}, rmse {r['rmse']:.4f})")
    (OUT / "f5c_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
