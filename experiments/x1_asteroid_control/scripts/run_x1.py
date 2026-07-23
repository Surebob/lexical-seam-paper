"""x1 — the asteroid belt as a non-linguistic Zipf control.

Asteroid sizes follow a famous power-law family (the collisional cascade,
Dohnanyi 1969) with known physics-driven breaks (e.g., near D ~ 100 km).
Question, in the spirit of the city/surname controls: Zipf-like, yes — but
does the LANGUAGE seam signature appear? Prediction from the paper: two-regime
structure may exist (collisional physics has its own breaks), but the
language fingerprint — s/V ~ 0.012 — should NOT (surnames gave 0.0266).

Data: Minor Planet Center MPCORB catalog; absolute magnitude H (fixed-width
cols 9-13). log-diameter is affine in H (D = 1329/sqrt(albedo) * 10^(-H/5)),
so the rank-size curve shape is albedo-independent. "Frequency" := diameter,
exactly analogous to city population.

Outputs: ../outputs/x1_summary.md, x1_asteroids.png
"""
from __future__ import annotations

import gzip
import io
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"
URLS = [
    "https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz",
    "https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz",
]
ALBEDO = 0.14
N_FIT = 20000  # log-spaced rank subsample for fitting


def fetch_H():
    import requests
    last = None
    for url in URLS:
        try:
            print(f"downloading {url} ...", flush=True)
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            raw = gzip.decompress(r.content)
            break
        except Exception as e:
            last = e
            continue
    else:
        raise RuntimeError(f"could not fetch MPCORB: {last}")
    H = []
    started = False
    for line in io.TextIOWrapper(io.BytesIO(raw), encoding="ascii", errors="ignore"):
        if not started:
            if line.startswith("----"):
                started = True
            continue
        s = line[8:13].strip()
        if s:
            try:
                H.append(float(s))
            except ValueError:
                pass
    return np.asarray(H, dtype=np.float64)


def fit_zm(logf, ranks):
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, ranks[-1], 384)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        rmse = float(np.sqrt(np.mean((d @ coef - logf) ** 2)))
        if best is None or rmse < best[0]:
            best = (rmse, float(-coef[1]), float(c))
    return best


def fit_erf(ranks, logf, V, seed=20260726, n_starts=24):
    from scipy.optimize import least_squares
    from scipy.special import erf as sp_erf
    log_ranks = np.log(ranks)
    LOWER = np.array([-100.0, 0.01, 0.0, -100.0, 0.01, 0.0, 20.0, 0.05, 0.05])
    UPPER = np.array([100.0, 30.0, 1e6, 100.0, 30.0, 1e6, V * 0.9, 10.0, 10.0])

    def predict(p):
        a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in p]
        z = (log_ranks - math.log(k)) / w_gate
        sigma = 0.5 * (1.0 - sp_erf(z))
        head = a1 - b1 * np.log(ranks + max(c1, 0.0))
        scale = max(1.0, k * w_tail)
        zz = np.clip((ranks - k) / scale, -60.0, 60.0)
        tail_rank = 1.0 + scale * np.log1p(np.exp(zz))
        tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
        return sigma * head + (1.0 - sigma) * tail

    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_starts):
        x0 = np.array([
            rng.uniform(logf.min() - 3, logf.max() + 3), rng.uniform(0.05, 5.0), rng.uniform(0, 1000),
            rng.uniform(logf.min() - 3, logf.max() + 3), rng.uniform(0.05, 5.0), rng.uniform(0, 1000),
            math.exp(rng.uniform(math.log(20), math.log(V * 0.5))),
            math.exp(rng.uniform(math.log(0.05), math.log(10))),
            math.exp(rng.uniform(math.log(0.05), math.log(10))),
        ])
        x0 = np.clip(x0, LOWER + 1e-9, UPPER - 1e-9)
        try:
            sol = least_squares(lambda p: predict(p) - logf, x0=x0,
                                bounds=(LOWER, UPPER), method="trf", max_nfev=9000)
        except Exception:
            continue
        rmse = float(np.sqrt(np.mean(sol.fun ** 2)))
        if best is None or rmse < best[0]:
            best = (rmse, sol.x, predict(sol.x))
    return best


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    H = fetch_H()
    n = len(H)
    logD = math.log(1329.0 / math.sqrt(ALBEDO)) - (H / 5.0) * math.log(10.0)
    logD = np.sort(logD)[::-1]
    V = len(logD)
    print(f"{V:,} asteroids with H", flush=True)

    idx = np.unique(np.clip(np.geomspace(1, V, N_FIT).astype(int) - 1, 0, V - 1))
    ranks = (idx + 1).astype(float)
    y = logD[idx]

    rmse_zm, b_zm, c_zm = fit_zm(y, ranks)
    rmse_erf, p_erf, pred = fit_erf(ranks, y, V)
    k, w_gate, w_tail = float(p_erf[6]), float(p_erf[7]), float(p_erf[8])
    s = k * w_tail
    n_pts = len(ranks)
    bic_zm = 3 * math.log(n_pts) + n_pts * math.log(max(rmse_zm ** 2, 1e-18))
    bic_erf = 9 * math.log(n_pts) + n_pts * math.log(max(rmse_erf ** 2, 1e-18))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5), dpi=150)
    ax.plot(np.log10(ranks), y / math.log(10), ".", ms=2, color="#666", label=f"{V:,} asteroids (MPC)")
    d = np.column_stack([np.ones_like(y), np.log(ranks + c_zm)])
    coef, *_ = np.linalg.lstsq(d, y, rcond=None)
    ax.plot(np.log10(ranks), (d @ coef) / math.log(10), "-", lw=1.2, color="#1f77b4",
            label=f"Zipf–Mandelbrot (RMSE {rmse_zm:.3f})")
    ax.plot(np.log10(ranks), pred / math.log(10), "-", lw=1.2, color="#d62728",
            label=f"two-regime erf (RMSE {rmse_erf:.3f})")
    ax.set_xlabel("log10 rank")
    ax.set_ylabel("log10 diameter (km, albedo 0.14)")
    ax.set_title("The asteroid belt under the lexical-seam instrument")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "x1_asteroids.png")

    lines = ["# x1 — asteroid-belt control\n"]
    lines.append(f"- {V:,} asteroids (MPCORB, absolute magnitude H; D via albedo {ALBEDO})")
    lines.append(f"- single ZM: b = {b_zm:.3f}, c = {c_zm:.1f}, RMSE {rmse_zm:.4f}")
    lines.append(f"- two-regime erf: RMSE {rmse_erf:.4f}, dBIC (ZM - erf) = {bic_zm - bic_erf:.0f}")
    lines.append(f"- fitted transition: k = {k:,.0f}, w_gate = {w_gate:.2f}, w_tail = {w_tail:.3f}")
    lines.append(f"- **s = k*w_tail = {s:,.0f};  s/V = {s / V:.4f}**  (language: 0.0118; surnames: 0.0266)")
    lines.append("")
    if abs(s / V - 0.0118) < 0.003:
        lines.append("Reading: s/V lands NEAR the language constant — unexpected; investigate before quoting.")
    else:
        lines.append("Reading: Zipf-family structure yes, but the LANGUAGE fingerprint no — "
                     "the belt fits its own constant, like every non-linguistic system tested.")
    (OUT / "x1_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines[1:]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
