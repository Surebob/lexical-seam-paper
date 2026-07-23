"""F4e — is sigma_T ~ 1.6 universal? PLN mixture fits on the expanded panel
(modern English registers, new languages, surname control) using f4d's estimator.

Reads processed frequency vectors from data/processed_panel/*.npy (built by f5a).
Outputs: ../outputs/f4e_panel_mixtures.csv, f4e_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PROC = REPO / "data" / "processed_panel"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

GH_X, GH_W = hermgauss(48)
LOG_GH_W = np.log(GH_W)


def pln_log_pn(ns, mu, sd):
    loglam = mu + math.sqrt(2.0) * sd * GH_X
    lam = np.exp(np.clip(loglam, -30, 30))
    n = ns[:, None]
    ll = n * loglam[None, :] - lam[None, :] - gammaln(n + 1.0)
    return logsumexp(ll + LOG_GH_W[None, :], axis=1) - 0.5 * math.log(math.pi)


def mixture_nll(theta, uniq, mult):
    pi_h = 0.001 + (0.30 - 0.001) / (1.0 + math.exp(-theta[0]))
    mu_t = theta[1]
    sd_t = math.exp(theta[2])
    mu_h = mu_t + 0.5 + math.exp(theta[3])
    sd_h = math.exp(theta[4])
    if not (0.03 < sd_t < 6.0 and 0.03 < sd_h < 6.0):
        return 1e12
    lp_h = pln_log_pn(uniq, mu_h, sd_h)
    lp_t = pln_log_pn(uniq, mu_t, sd_t)
    lp_mix = np.logaddexp(math.log(pi_h) + lp_h, math.log(1 - pi_h) + lp_t)
    z = np.array([0.0])
    lp0 = np.logaddexp(math.log(pi_h) + pln_log_pn(z, mu_h, sd_h),
                       math.log(1 - pi_h) + pln_log_pn(z, mu_t, sd_t))[0]
    return -float(np.dot(mult, lp_mix - math.log1p(-min(math.exp(lp0), 1 - 1e-12))))


def fit(freqs):
    uniq, mult = np.unique(freqs, return_counts=True)
    mult = mult.astype(np.float64)
    best = None
    for x0 in [
        [-3.0, 0.3, math.log(1.2), math.log(4.0), math.log(1.0)],
        [-2.0, 0.0, math.log(1.0), math.log(6.0), math.log(1.4)],
        [-4.0, 0.5, math.log(1.5), math.log(3.0), math.log(0.7)],
        [-1.5, -0.3, math.log(0.8), math.log(5.0), math.log(1.8)],
    ]:
        try:
            res = minimize(mixture_nll, x0=x0, args=(uniq, mult), method="Nelder-Mead",
                           options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-7})
        except Exception:
            continue
        if best is None or res.fun < best.fun:
            best = res
    t = best.x
    return {
        "pi_h": 0.001 + (0.30 - 0.001) / (1.0 + math.exp(-t[0])),
        "mu_t": float(t[1]), "sd_t": float(math.exp(t[2])),
        "mu_h": float(t[1]) + 0.5 + float(math.exp(t[3])), "sd_h": float(math.exp(t[4])),
        "nll": float(best.fun),
    }


def main():
    rows = []
    for path in sorted(PROC.glob("*.npy")):
        freqs = np.load(path).astype(np.float64)
        r = fit(freqs)
        rows.append({"corpus": path.stem, "V": len(freqs), "tokens": int(freqs.sum()), **r})
        print(f"{path.stem:18} V={len(freqs):7} pi_H={100*r['pi_h']:5.2f}% sd_H={r['sd_h']:5.2f} "
              f"sd_T={r['sd_t']:5.2f} mu_T={r['mu_t']:6.2f}", flush=True)

    with open(OUT / "f4e_panel_mixtures.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lang = [r for r in rows if r["corpus"] != "census_surnames"]
    lines = ["# F4e — sigma_T universality on the expanded panel\n"]
    lines.append(f"- language corpora sd_T: median {float(np.median([r['sd_t'] for r in lang])):.3f} "
                 f"(range {min(r['sd_t'] for r in lang):.3f}-{max(r['sd_t'] for r in lang):.3f}); "
                 f"English-books reference (f4d): ~1.6 (1.50-2.06)")
    for r in rows:
        lines.append(f"- {r['corpus']}: sd_T={r['sd_t']:.3f} sd_H={r['sd_h']:.3f} pi_H={100*r['pi_h']:.2f}%")
    (OUT / "f4e_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
