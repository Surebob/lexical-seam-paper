"""F4 — Does a two-LOGNORMAL type-population mixture generate the empirical erf-gate
preference? (Empirical half of the erf derivation, ROADMAP C1.)

Hypothesis: if the two lexical populations have log-normal frequency distributions,
the head-population share as a function of log-rank is approximately an erf, so the
decoupled 9-parameter fit should PREFER the erf gate on such synthetic corpora —
and should not prefer it on two-Pareto (power-law population) mixtures.

Datasets (2 seeds each):
  LN1..LN4  two-lognormal mixtures (varying gap/widths/head size)
  PA1..PA3  two-Pareto mixtures (power-law populations)
  SL1       single lognormal (degenerate control)

Fit: all 5 gates (logistic, tanh, erf, algebraic, arctan), 9 free params, 20 starts,
bounds identical to the canonical sweep. Winner per dataset = min RMSE (p equal).

Outputs: ../outputs/f4_fits.csv, f4_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"

LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 1000.0, 10.0, 10.0])
N_STARTS = 20
MAX_NFEV = 8000
BASE_SEED = 20260721

CONFIGS = []
for seed in (1, 2):
    CONFIGS += [
        ("LN1", "lognormal", dict(NH=250, NT=18000, muH=math.log(1500), sH=1.0, muT=math.log(2.5), sT=1.5), seed),
        ("LN2", "lognormal", dict(NH=400, NT=12000, muH=math.log(800), sH=1.2, muT=math.log(3.5), sT=1.8), seed),
        ("LN3", "lognormal", dict(NH=150, NT=25000, muH=math.log(3000), sH=0.8, muT=math.log(2.0), sT=1.3), seed),
        ("LN4", "lognormal", dict(NH=300, NT=9000, muH=math.log(600), sH=1.5, muT=math.log(4.0), sT=2.0), seed),
        ("PA1", "pareto", dict(NH=250, NT=18000, aH=1.05, mH=40.0, aT=1.25, mT=1.0), seed),
        ("PA2", "pareto", dict(NH=400, NT=12000, aH=0.95, mH=25.0, aT=1.45, mT=1.0), seed),
        ("PA3", "pareto", dict(NH=150, NT=25000, aH=1.15, mH=80.0, aT=1.10, mT=1.0), seed),
        ("SL1", "single_lognormal", dict(N=15000, mu=math.log(4.0), s=2.2), seed),
    ]


def gen_dataset(kind: str, p: dict, seed: int):
    rng = np.random.default_rng(BASE_SEED + seed * 1000 + hash(kind) % 997)
    if kind == "lognormal":
        lamH = np.exp(rng.normal(p["muH"], p["sH"], p["NH"]))
        lamT = np.exp(rng.normal(p["muT"], p["sT"], p["NT"]))
        lam = np.concatenate([lamH, lamT])
    elif kind == "pareto":
        lamH = p["mH"] * (1.0 + rng.pareto(p["aH"], p["NH"]))
        lamT = p["mT"] * (1.0 + rng.pareto(p["aT"], p["NT"]))
        lam = np.concatenate([lamH, lamT])
    else:
        lam = np.exp(rng.normal(p["mu"], p["s"], p["N"]))
    counts = rng.poisson(lam)
    counts = counts[counts > 0]
    freqs = np.sort(counts)[::-1].astype(np.float64)
    return freqs


def gate_sigma(gate: str, z: np.ndarray):
    if gate == "logistic":
        return 1.0 / (1.0 + np.exp(np.clip(z, -60, 60)))
    if gate == "tanh":
        return 0.5 * (1.0 - np.tanh(z))
    if gate == "erf":
        return 0.5 * (1.0 - sp_erf(z))
    if gate == "algebraic":
        return 0.5 * (1.0 - z / np.sqrt(1.0 + z * z))
    if gate == "arctan":
        return 0.5 - np.arctan(z) / math.pi
    raise ValueError(gate)


def predict(ranks, log_ranks, params, gate):
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    z = (log_ranks - math.log(k)) / w_gate
    sigma = gate_sigma(gate, z)
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
        for c in np.concatenate([[0.0], np.geomspace(1e-6, rr[-1], 512)]):
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
    cfg_name, kind, seed, gate, freqs_list = task
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(BASE_SEED + hash((cfg_name, seed, gate)) % (2**31))
    anchor6 = piecewise_anchor(ranks, y)

    def resid(p):
        return predict(ranks, log_ranks, p, gate) - y

    best_rmse, best_p = None, None
    for i in range(N_STARTS):
        if i < N_STARTS // 2:
            x0 = np.array(anchor6 + [0.0, 0.0, 0.0])
            x0[:6] += rng.normal(0.0, [6.0, 0.3, 90.0, 6.0, 0.3, 90.0])
            x0[6] = 400.0 + rng.normal(0.0, 150.0)
            x0[7] = 0.5 * math.exp(rng.normal(0.0, 0.45))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.45))
        else:
            x0 = np.array([
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(1000))),
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
    bic = 9 * math.log(V) + V * math.log(best_rmse**2)
    return {"config": cfg_name, "kind": kind, "seed": seed, "gate": gate, "V": V,
            "rmse": best_rmse, "bic": bic, "k": float(best_p[6]),
            "w_gate": float(best_p[7]), "w_tail": float(best_p[8])}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    datasets = {}
    for cfg_name, kind, p, seed in CONFIGS:
        freqs = gen_dataset(kind, p, seed)
        datasets[(cfg_name, seed)] = (kind, freqs)
        print(f"gen {cfg_name} seed{seed}: V={len(freqs)} tokens={int(freqs.sum())} fmax={int(freqs[0])}", flush=True)

    tasks = []
    for (cfg_name, seed), (kind, freqs) in datasets.items():
        for gate in ["logistic", "tanh", "erf", "algebraic", "arctan"]:
            tasks.append((cfg_name, kind, seed, gate, freqs.tolist()))

    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['config']} s{res['seed']} {res['gate']:9} rmse={res['rmse']:.5f}", flush=True)

    with open(OUT / "f4_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lines = ["# F4 gate recovery on synthetic population mixtures\n"]
    lines.append("| dataset | kind | winner | erf_bic - best_other | spread(min..max BIC gap) |")
    lines.append("|---|---|---|---:|---:|")
    win_counts = Counter()
    for (cfg_name, seed), (kind, _) in sorted(datasets.items()):
        sub = [r for r in rows if r["config"] == cfg_name and r["seed"] == seed]
        indep = [r for r in sub if r["gate"] != "tanh"]
        best = min(indep, key=lambda r: r["bic"])
        others = [r for r in indep if r["gate"] != "erf"]
        erf_row = next(r for r in indep if r["gate"] == "erf")
        gap = erf_row["bic"] - min(r["bic"] for r in others)
        spread = max(r["bic"] for r in indep) - min(r["bic"] for r in indep)
        win_counts[(kind, best["gate"])] += 1
        lines.append(f"| {cfg_name} s{seed} | {kind} | {best['gate']} | {gap:.2f} | {spread:.2f} |")
    lines.append("\n## Winner counts by generative family")
    for (kind, gate), nwin in sorted(win_counts.items()):
        lines.append(f"- {kind} -> {gate}: {nwin}")
    (OUT / "f4_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
