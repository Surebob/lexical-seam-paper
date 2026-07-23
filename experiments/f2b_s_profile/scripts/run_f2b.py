"""F2b — Direct profile likelihood for s = k*w_tail (the f2 candidate law s ≈ 1.2% V).

For each English corpus, fix s on an 11-point grid of s/V in [0.002, 0.08]; free
parameters are (a1,b1,c1,a2,b2,c2,k,w_gate) with w_tail := s/k tied to the fixed s.
BIC profile over s gives the per-corpus identified interval for s; regression of the
profile-minimum s on V re-measures the law with s treated as a first-class parameter.

Outputs: ../outputs/f2b_profile_points.csv, f2b_per_corpus.csv, f2b_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
F2CSV = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

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

S_OVER_V = np.geomspace(0.002, 0.08, 11)
N_STARTS = 6
MAX_NFEV = 8000
SEED = 20260722

L8 = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05])
U8 = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 1000.0, 10.0])


@lru_cache(maxsize=32)
def load(fname: str):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    s, e = 0, len(text)
    for m in SM:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            s = nl + 1 if nl != -1 else i
            break
    for m in EMK:
        i = text.find(m)
        if i != -1:
            e = i
            break
    counts = Counter(TOKEN_RE.findall(text[s:e].lower()))
    freqs = np.array(sorted(counts.values(), reverse=True), dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    return ranks, np.log(ranks), np.log(freqs)


def predict_fixed_s(ranks, log_ranks, s_fixed, p8):
    a1, b1, c1, a2, b2, c2, k, w_gate = [float(v) for v in p8]
    z = (log_ranks - math.log(k)) / w_gate
    sigma = 0.5 * (1.0 - sp_erf(z))
    head = a1 - b1 * np.log(ranks + max(c1, 0.0))
    scale = max(1.0, s_fixed)
    zz = np.clip((ranks - k) / scale, -60.0, 60.0)
    tail_rank = 1.0 + scale * np.log1p(np.exp(zz))
    tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
    return sigma * head + (1.0 - sigma) * tail


def fit_task(task):
    name, fname, s_fixed, warm = task
    ranks, log_ranks, y = load(fname)
    n = len(ranks)
    rng = np.random.default_rng(SEED + hash((fname, round(math.log(s_fixed) * 1000))) % (2**31))

    def resid(p8):
        return predict_fixed_s(ranks, log_ranks, s_fixed, p8) - y

    starts = [np.clip(np.asarray(warm, dtype=np.float64), L8 + 1e-8, U8 - 1e-8)]
    for _ in range(3):
        j = starts[0].copy()
        j[:6] += rng.normal(0.0, [4.0, 0.2, 60.0, 4.0, 0.2, 60.0])
        j[6] *= math.exp(rng.normal(0.0, 0.3))
        j[7] *= math.exp(rng.normal(0.0, 0.4))
        starts.append(np.clip(j, L8 + 1e-8, U8 - 1e-8))
    for _ in range(2):
        r = np.array([
            rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
            rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
            math.exp(rng.uniform(math.log(20), math.log(1000))),
            math.exp(rng.uniform(math.log(0.05), math.log(10))),
        ])
        starts.append(np.clip(r, L8 + 1e-8, U8 - 1e-8))

    best_rmse, best_p = None, None
    for x0 in starts:
        try:
            sol = least_squares(resid, x0=x0, bounds=(L8, U8), method="trf", max_nfev=MAX_NFEV)
        except Exception:
            continue
        rr = float(np.sqrt(np.mean(sol.fun**2)))
        if best_rmse is None or rr < best_rmse:
            best_rmse, best_p = rr, sol.x
    bic = 9 * math.log(n) + n * math.log(best_rmse**2)
    return {"corpus": name, "V": n, "s": float(s_fixed), "s_over_V": float(s_fixed / n),
            "rmse": best_rmse, "bic": bic, "k": float(best_p[6]), "w_gate": float(best_p[7])}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2CSV, newline="", encoding="utf-8"))}
    tasks = []
    for name, fname in CORPORA:
        ranks, _, _ = load(fname)
        V = len(ranks)
        warm = [12.0, 1.0, 10.0, 15.0, 1.5, 100.0,
                float(f2[name]["k_profile_min"]), float(f2[name]["w_gate_at_min"])]
        s_grid = sorted(set(list(S_OVER_V * V) + [float(f2[name]["s_at_min"])]))
        for s in s_grid:
            tasks.append((name, fname, float(s), warm))

    print(f"{len(tasks)} tasks", flush=True)
    rows = []
    with Pool(processes=10) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['corpus'][:26]:27} s/V={res['s_over_V']:.4f} rmse={res['rmse']:.6f}", flush=True)

    with open(OUT / "f2b_profile_points.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    per = []
    for name, fname in CORPORA:
        pts = sorted((r for r in rows if r["corpus"] == name), key=lambda r: r["s"])
        bics = np.array([p["bic"] for p in pts])
        ss = np.array([p["s"] for p in pts])
        i0 = int(np.argmin(bics))
        d = bics - bics[i0]
        in2 = ss[d <= 2.0]
        per.append({"corpus": name, "V": pts[0]["V"], "s_min": float(ss[i0]),
                    "s_over_V_min": float(ss[i0] / pts[0]["V"]),
                    "s_lo_d2": float(in2.min()), "s_hi_d2": float(in2.max()),
                    "width_log10_d2": float(np.log10(in2.max() / in2.min()))})
    with open(OUT / "f2b_per_corpus.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(per[0].keys()))
        w.writeheader()
        w.writerows(per)

    lV = np.log(np.array([p["V"] for p in per], dtype=float))
    ls = np.log(np.array([p["s_min"] for p in per], dtype=float))
    X = np.column_stack([np.ones_like(lV), lV])
    b, *_ = np.linalg.lstsq(X, ls, rcond=None)
    ssr = float(np.sum((ls - X @ b) ** 2))
    r2 = 1 - ssr / float(np.sum((ls - ls.mean()) ** 2))
    se = math.sqrt(ssr / (len(lV) - 2) / float(np.sum((lV - lV.mean()) ** 2)))
    lines = ["# F2b direct s-profile — summary\n"]
    lines.append(f"- median identified s-interval width (dBIC<=2): {float(np.median([p['width_log10_d2'] for p in per])):.3f} log10")
    lines.append(f"- log s_min ~ log V: beta={b[1]:.4f} 95% CI [{b[1]-1.96*se:.4f},{b[1]+1.96*se:.4f}] R2={r2:.4f}")
    lines.append(f"- median s/V at profile minimum: {float(np.median([p['s_over_V_min'] for p in per])):.4f}")
    (OUT / "f2b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
