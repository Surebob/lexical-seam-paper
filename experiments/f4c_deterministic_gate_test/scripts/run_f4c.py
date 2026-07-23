"""F4c — Deterministic (noise-free) gate discrimination + properly-sized Pareto heads.

Fixes f4b's two confounds:
  1. NOISE-FREE curves: evaluate each generative family's rank curve as the exact
     quantile function of the type-rate distribution (no Poisson sampling). If the
     erf preference is driven by family SHAPE, it must appear here; if f4b's
     "Pareto twins pick erf too" was Poisson smearing, it must vanish here.
  2. Pareto heads sized to the GMM head mass (pi_H * V types, not a fixed 130).

For each of 12 corpora x {lognormal_mix, pareto_mix} x noise in {off, on}:
  build curve, fit 5 gates (16 starts). P4: correlate erf-vs-best-other margin with
  head-component Gaussianity proxy across corpora.

Outputs: ../outputs/f4c_fits.csv, f4c_summary.md
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
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
TWINS = REPO / "experiments" / "f4b_calibrated_twins" / "outputs" / "f4b_twins.csv"
F1CSV = REPO / "experiments" / "f1_fresh_reproduction" / "outputs" / "f1_per_corpus.csv"

SEED = 20260722
N_STARTS = 16
MAX_NFEV = 8000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 1000.0, 10.0, 10.0])


def quantile_curve_lognormal_mix(V, pi_h, mu_h, sd_h, mu_t, sd_t):
    """Exact rank curve: rate of the i-th largest type under the mixture."""
    q = (np.arange(1, V + 1) - 0.5) / V  # upper-tail probabilities
    # solve S(x) = q where S = pi_h*(1-Phi((x-mu_h)/sd_h)) + (1-pi_h)*(1-Phi((x-mu_t)/sd_t))
    lo = min(mu_h, mu_t) - 8 * max(sd_h, sd_t)
    hi = max(mu_h, mu_t) + 8 * max(sd_h, sd_t)
    xs = np.linspace(lo, hi, 4000)
    S = pi_h * norm.sf((xs - mu_h) / sd_h) + (1 - pi_h) * norm.sf((xs - mu_t) / sd_t)
    x_of_q = np.interp(q, S[::-1], xs[::-1])
    return np.exp(x_of_q)


def quantile_curve_pareto_mix(V, pi_h, a_h, m_h, a_t, m_t):
    q = (np.arange(1, V + 1) - 0.5) / V
    lo = math.log(min(m_h, m_t)) - 0.5
    hi = math.log(max(m_h, m_t)) + 14.0
    xs = np.linspace(lo, hi, 4000)
    lam = np.exp(xs)
    Sh = np.where(lam < m_h, 1.0, (lam / m_h) ** (-a_h))
    St = np.where(lam < m_t, 1.0, (lam / m_t) ** (-a_t))
    S = pi_h * Sh + (1 - pi_h) * St
    x_of_q = np.interp(q, S[::-1], xs[::-1])
    return np.exp(x_of_q)


def gate_sigma(gate, z):
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
    corpus, family, noise, gate, freqs_list = task
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + hash((corpus, family, noise, gate)) % (2**31))
    anchor6 = piecewise_anchor(ranks, y)

    def resid(p):
        return predict(ranks, log_ranks, p, gate) - y

    best_rmse = None
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
            best_rmse = r
    bic = 9 * math.log(V) + V * math.log(best_rmse**2)
    return {"corpus": corpus, "family": family, "noise": noise, "gate": gate, "V": V,
            "rmse": best_rmse, "bic": bic}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    twins = {r["corpus"]: r for r in csv.DictReader(open(TWINS, newline="", encoding="utf-8"))
             if r["family"] == "lognormal_twin"}
    f1 = {r["corpus"]: r for r in csv.DictReader(open(F1CSV, newline="", encoding="utf-8"))}
    rng = np.random.default_rng(SEED)

    tasks = []
    meta = []
    for corpus, t in twins.items():
        V = int(t["V"])
        pi_h = float(t["gmm_pi_head"]) if "gmm_pi_head" in t and t["gmm_pi_head"] else 0.02
        mu_h, sd_h = float(t["gmm_mu_head"]), float(t["gmm_sd_head"])
        mu_t, sd_t = float(t["gmm_mu_tail"]), float(t["gmm_sd_tail"])
        af = float(f1[corpus]["alpha_function"])
        ac = float(f1[corpus]["alpha_content"])
        # Pareto populations sized to the SAME head mass as the GMM; levels matched
        # so the head is m_h-times the tail floor with comparable dynamic range.
        a_h, a_t = 1.0 / max(af, 0.2), 1.0 / max(ac, 0.2)
        m_h = math.exp(mu_h - 1.0)
        m_t = max(math.exp(mu_t - 0.5), 0.5)

        ln_curve = quantile_curve_lognormal_mix(V, pi_h, mu_h, sd_h, mu_t, sd_t)
        pa_curve = quantile_curve_pareto_mix(V, pi_h, a_h, m_h, a_t, m_t)
        for fam, curve in [("lognormal_mix", ln_curve), ("pareto_mix", pa_curve)]:
            clean = np.maximum(curve, 1e-3)
            noisy = np.random.default_rng(SEED + hash((corpus, fam)) % 999983).poisson(clean)
            noisy = np.sort(noisy[noisy > 0])[::-1].astype(np.float64)
            for noise, fr in [("off", clean), ("on", noisy)]:
                meta.append({"corpus": corpus, "family": fam, "noise": noise, "V": len(fr)})
                for gate in ["logistic", "tanh", "erf", "algebraic", "arctan"]:
                    tasks.append((corpus, fam, noise, gate, fr.tolist()))

    print(f"{len(tasks)} tasks", flush=True)
    rows = []
    with Pool(processes=10) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['corpus'][:20]:21} {res['family']:14} n={res['noise']:3} "
                  f"{res['gate']:9} rmse={res['rmse']:.5f}", flush=True)

    with open(OUT / "f4c_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    from collections import Counter as C2
    wc = C2()
    lines = ["# F4c deterministic gate discrimination — summary\n"]
    lines.append("| corpus | family | noise | winner | erf margin | spread |")
    lines.append("|---|---|---|---|---:|---:|")
    for m in meta:
        sub = [r for r in rows if r["corpus"] == m["corpus"] and r["family"] == m["family"]
               and r["noise"] == m["noise"] and r["gate"] != "tanh"]
        if not sub:
            continue
        best = min(sub, key=lambda r: r["bic"])
        erf_b = next(r["bic"] for r in sub if r["gate"] == "erf")
        other = min(r["bic"] for r in sub if r["gate"] != "erf")
        spread = max(r["bic"] for r in sub) - min(r["bic"] for r in sub)
        wc[(m["family"], m["noise"], best["gate"])] += 1
        lines.append(f"| {m['corpus']} | {m['family']} | {m['noise']} | {best['gate']} | {erf_b - other:.1f} | {spread:.1f} |")
    lines.append("\n## Winner counts (family, noise) -> gate")
    for (fam, noise, gate), nwin in sorted(wc.items()):
        lines.append(f"- {fam} noise={noise} -> {gate}: {nwin}")
    (OUT / "f4c_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
