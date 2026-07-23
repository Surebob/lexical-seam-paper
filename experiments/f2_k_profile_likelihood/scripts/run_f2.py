"""F2 — Per-corpus profile likelihood for the transition centre k (decoupled erf model).

Question (ROADMAP B1): is the wide decoupled-erf k-scaling CI an identifiability
artifact? For each English corpus we FIX k on a 13-point grid over [20, 1000] and
optimize the remaining 8 parameters (warm-started from the canonical 2026-04-18 erf
fit + jittered/random restarts). The BIC profile over k gives a per-corpus identified
interval; regressing the profile-minimum k on V re-tests the scaling law with
degeneracy handled explicitly.

Model equations reproduced from src/s2_decoupled/shared (erf gate, softplus tail
coordinate with s = max(1, k*w_tail)); bounds identical to the canonical sweep.

Outputs: ../outputs/f2_profile_points.csv, f2_per_corpus.csv, f2_summary.md
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
ARCHIVE_FITS = Path(r"C:\Users\Greg Kara\Desktop\temporary\emlexperiment\results\s2_v3_windows_full_outputs_2026-04-18\s2_v3_per_fit_results.csv")

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

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

K_GRID = np.geomspace(20.0, 1000.0, 13)
N_JITTER = 3
N_RANDOM = 2
MAX_NFEV = 8000
BASE_SEED = 20260721

LOWER8 = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 0.05, 0.05])
UPPER8 = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 10.0, 10.0])


@lru_cache(maxsize=32)
def load_dataset(fname: str):
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
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    freqs = np.array([f for _, f in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    return ranks, np.log(ranks), np.log(freqs)


def predict(ranks, log_ranks, k, p8):
    a1, b1, c1, a2, b2, c2, w_gate, w_tail = [float(v) for v in p8]
    sigma = 0.5 * (1.0 - sp_erf((log_ranks - math.log(k)) / w_gate))
    head = a1 - b1 * np.log(ranks + max(c1, 0.0))
    scale = max(1.0, k * w_tail)
    z = np.clip((ranks - k) / scale, -60.0, 60.0)
    tail_rank = 1.0 + scale * np.log1p(np.exp(z))
    tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
    return sigma * head + (1.0 - sigma) * tail


def fit_fixed_k(task):
    name, fname, k, warm8 = task
    ranks, log_ranks, y = load_dataset(fname)
    n = len(ranks)
    rng = np.random.default_rng(BASE_SEED + hash((fname, round(math.log(k) * 1000))) % (2**31))

    def resid(p8):
        return predict(ranks, log_ranks, k, p8) - y

    starts = [np.clip(np.asarray(warm8, dtype=np.float64), LOWER8 + 1e-8, UPPER8 - 1e-8)]
    for _ in range(N_JITTER):
        j = starts[0].copy()
        j[:6] += rng.normal(0.0, [4.0, 0.2, 60.0, 4.0, 0.2, 60.0])
        j[6] *= math.exp(rng.normal(0.0, 0.4))
        j[7] *= math.exp(rng.normal(0.0, 0.4))
        starts.append(np.clip(j, LOWER8 + 1e-8, UPPER8 - 1e-8))
    for _ in range(N_RANDOM):
        r = np.empty(8)
        r[0] = rng.uniform(y.min() - 3, y.max() + 3)
        r[1] = rng.uniform(0.5, 3.0)
        r[2] = rng.uniform(0.0, 1000.0)
        r[3] = rng.uniform(y.min() - 3, y.max() + 3)
        r[4] = rng.uniform(0.5, 3.0)
        r[5] = rng.uniform(0.0, 1000.0)
        r[6] = math.exp(rng.uniform(math.log(0.05), math.log(10.0)))
        r[7] = math.exp(rng.uniform(math.log(0.05), math.log(10.0)))
        starts.append(np.clip(r, LOWER8 + 1e-8, UPPER8 - 1e-8))

    best_rmse, best_p = None, None
    for x0 in starts:
        try:
            sol = least_squares(resid, x0=x0, bounds=(LOWER8, UPPER8), method="trf", max_nfev=MAX_NFEV)
        except Exception:
            continue
        r = float(np.sqrt(np.mean(sol.fun**2)))
        if best_rmse is None or r < best_rmse:
            best_rmse, best_p = r, sol.x
    bic = 9 * math.log(n) + n * math.log(best_rmse**2)
    return {"corpus": name, "V": n, "k": float(k), "rmse": best_rmse, "bic": bic,
            "w_gate": float(best_p[6]), "w_tail": float(best_p[7]), "s": float(k * best_p[7])}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    warm = {}
    with open(ARCHIVE_FITS, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["gate"] == "erf":
                warm[row["corpus"]] = row

    tasks = []
    for name, fname in CORPORA:
        w = warm[name]
        # canonical params are a1,b1,c1,a2,b2,c2 absent from per-fit CSV -> derive warm8
        # from a piecewise anchor instead: use generic slopes and the canonical widths.
        warm8 = [12.0, 1.0, 10.0, 15.0, 1.5, 100.0, float(w["w_gate"]), float(w["w_tail"])]
        kgrid = sorted(set(list(K_GRID) + [float(w["k"])]))
        for k in kgrid:
            tasks.append((name, fname, float(k), warm8))

    print(f"{len(tasks)} tasks", flush=True)
    rows = []
    with Pool(processes=14) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_fixed_k, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['corpus'][:28]:29} k={res['k']:8.1f} rmse={res['rmse']:.6f}", flush=True)

    with open(OUT / "f2_profile_points.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(rows)

    per = []
    for name, fname in CORPORA:
        pts = sorted((r for r in rows if r["corpus"] == name), key=lambda r: r["k"])
        bics = np.array([p["bic"] for p in pts])
        ks = np.array([p["k"] for p in pts])
        i0 = int(np.argmin(bics))
        d = bics - bics[i0]
        in2 = ks[d <= 2.0]
        in10 = ks[d <= 10.0]
        per.append({
            "corpus": name, "V": pts[0]["V"], "k_profile_min": float(ks[i0]),
            "bic_min": float(bics[i0]),
            "k_lo_d2": float(in2.min()), "k_hi_d2": float(in2.max()),
            "width_log10_d2": float(np.log10(in2.max() / in2.min())),
            "k_lo_d10": float(in10.min()), "k_hi_d10": float(in10.max()),
            "width_log10_d10": float(np.log10(in10.max() / in10.min())),
            "s_at_min": float(pts[i0]["s"]), "w_tail_at_min": float(pts[i0]["w_tail"]),
            "w_gate_at_min": float(pts[i0]["w_gate"]),
        })

    with open(OUT / "f2_per_corpus.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(per[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(per)

    lV = np.log(np.array([p["V"] for p in per], dtype=float))
    lk = np.log(np.array([p["k_profile_min"] for p in per], dtype=float))
    ls = np.log(np.array([p["s_at_min"] for p in per], dtype=float))
    X = np.column_stack([np.ones_like(lV), lV])
    bk, *_ = np.linalg.lstsq(X, lk, rcond=None)
    rk = 1 - np.sum((lk - X @ bk) ** 2) / np.sum((lk - lk.mean()) ** 2)
    sek = math.sqrt(float(np.sum((lk - X @ bk) ** 2)) / (len(lV) - 2) / float(np.sum((lV - lV.mean()) ** 2)))
    bs, *_ = np.linalg.lstsq(X, ls, rcond=None)
    rs = 1 - np.sum((ls - X @ bs) ** 2) / np.sum((ls - ls.mean()) ** 2)
    ses = math.sqrt(float(np.sum((ls - X @ bs) ** 2)) / (len(lV) - 2) / float(np.sum((lV - lV.mean()) ** 2)))

    lines = ["# F2 k-profile likelihood — summary\n"]
    lines.append(f"- corpora: {len(per)}; k grid: {len(K_GRID)}+canonical per corpus; starts per point: {1+N_JITTER+N_RANDOM}")
    lines.append(f"- median identified k-interval width (dBIC<=2): {float(np.median([p['width_log10_d2'] for p in per])):.3f} log10 units")
    lines.append(f"- median identified k-interval width (dBIC<=10): {float(np.median([p['width_log10_d10'] for p in per])):.3f} log10 units")
    lines.append(f"- corpora with dBIC<=2 interval spanning >0.5 log10 (>3.2x in k): {sum(1 for p in per if p['width_log10_d2'] > 0.5)}/{len(per)}")
    lines.append(f"\n## Scaling regressions (profile-minimum estimates)")
    lines.append(f"- log k_prof ~ log V: beta={bk[1]:.4f} 95% CI [{bk[1]-1.96*sek:.4f},{bk[1]+1.96*sek:.4f}] R2={rk:.4f}")
    lines.append(f"- log s_prof ~ log V: beta={bs[1]:.4f} 95% CI [{bs[1]-1.96*ses:.4f},{bs[1]+1.96*ses:.4f}] R2={rs:.4f}")
    (OUT / "f2_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
