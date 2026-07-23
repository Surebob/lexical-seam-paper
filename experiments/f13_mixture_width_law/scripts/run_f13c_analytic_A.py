"""F13c — the analytic width: ambiguous-mass A at the physical basins.

f13/f13b established computationally that the mixture generates the width law
at the physical basin. This script tests the ANALYTIC reduction: define, from
mixture parameters alone (no simulation, no fitting),

    A = N_total * Integral over the handover band of
        [pi_H phi_H(x) + (1-pi_H) phi_T(x)] * (1 - exp(-e^x)) dx

where the handover band is where the head's posterior share crosses from 90%
to 10%, and (1 - exp(-e^x)) is the count-floor observation probability. A
counts the expected number of OBSERVED types in the zone where the two
populations genuinely mix — the candidate analytic seam width.

Basins: f4d ML params for the 16 corpora where ML already reproduces the
width; f13b's pred_err-selected basins for the 8 healed corpora. (Kant runs
separately in run_f13c_kant.py.)

Test: A vs measured s (corr, median ratio), A/V vs 0.0118, and the s ~ A
regression. If A ~ s with a stable O(1) factor, the derivation of the 1.2%%
constant reduces to explaining the fitted-basin regularities that set A/V.

Outputs: ../outputs/f13c_analytic.csv, f13c_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import logsumexp
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
F4D = REPO / "experiments" / "f4d_poisson_lognormal_mixture" / "outputs" / "f4d_mixture_fits.csv"
F2 = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"
F13B = OUT / "f13b_landscape.csv"

NLL_SLACK = 50.0
GH_X, GH_W = hermgauss(48)

HEALED = {"Complete Works of Shakespeare", "War and Peace", "King James Bible",
          "Federalist Papers", "Origin of Species", "Wealth of Nations",
          "Les Miserables", "Principia Ethica"}
RESIST = {"Critique of Pure Reason"}


def pln_p0(mu, sd):
    lam = np.exp(np.clip(mu + math.sqrt(2.0) * sd * GH_X, -30, 30))
    return float(np.exp(logsumexp(-lam + np.log(GH_W)) - 0.5 * math.log(math.pi)))


def analytic_A(pi_h, mu_h, sd_h, mu_t, sd_t, V_obs):
    p0 = pi_h * pln_p0(mu_h, sd_h) + (1 - pi_h) * pln_p0(mu_t, sd_t)
    n_total = V_obs / max(1e-9, 1 - p0)
    xs = np.linspace(min(mu_t - 7 * sd_t, mu_h - 6 * sd_h),
                     max(mu_t + 7 * sd_t, mu_h + 6 * sd_h), 40001)
    dx = xs[1] - xs[0]
    lg = (math.log(pi_h) + norm.logpdf(xs, mu_h, sd_h)) - \
         (math.log1p(-pi_h) + norm.logpdf(xs, mu_t, sd_t))
    share = 1.0 / (1.0 + np.exp(-np.clip(lg, -500, 500)))
    band = (share > 0.10) & (share < 0.90)
    if not band.any():
        return float("nan"), n_total, float("nan"), float("nan")
    dens = pi_h * norm.pdf(xs, mu_h, sd_h) + (1 - pi_h) * norm.pdf(xs, mu_t, sd_t)
    p_obs = 1.0 - np.exp(-np.exp(np.clip(xs, -30, 30)))
    A = n_total * float(np.sum(dens[band] * p_obs[band]) * dx)
    lo, hi = float(xs[band].min()), float(xs[band].max())
    return A, n_total, lo, hi


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f4d = {r["corpus"]: r for r in csv.DictReader(open(F4D, newline="", encoding="utf-8"))}
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    land = list(csv.DictReader(open(F13B, newline="", encoding="utf-8"))) if F13B.exists() else []

    rows = []
    for name, r4 in f4d.items():
        V = float(r4["V"])
        s_real = float(f2[name]["s_at_min"])
        if name in RESIST:
            basin_src = "resists"
            continue
        if name in HEALED and land:
            sub = [x for x in land if x["corpus"] == name]
            mn = min(float(x["nll"]) for x in sub)
            adm = [x for x in sub if float(x["nll"]) <= mn + NLL_SLACK]
            pick = min(adm, key=lambda x: float(x["pred_err"]))
            pi_h, mu_h, sd_h = float(pick["pi_h"]), float(pick["mu_h"]), float(pick["sd_h"])
            mu_t, sd_t = float(pick["mu_t"]), float(pick["sd_t"])
            basin_src = "pred-selected"
        else:
            pi_h, mu_h, sd_h = float(r4["pi_h"]), float(r4["mu_h"]), float(r4["sd_h"])
            mu_t, sd_t = float(r4["mu_t"]), float(r4["sd_t"])
            basin_src = "ML"
        A, n_total, lo, hi = analytic_A(pi_h, mu_h, sd_h, mu_t, sd_t, V)
        rows.append({"corpus": name, "basin": basin_src, "V": V, "s_real": s_real,
                     "A": round(A, 1), "A_over_s": round(A / s_real, 3),
                     "A_over_V": round(A / V, 5),
                     "band_lo": round(lo, 2), "band_hi": round(hi, 2),
                     "pi_h": round(pi_h, 5), "gap": round(mu_h - mu_t, 3)})

    with open(OUT / "f13c_analytic.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    A_ = np.array([r["A"] for r in rows])
    s_ = np.array([r["s_real"] for r in rows])
    V_ = np.array([r["V"] for r in rows])
    ok = np.isfinite(A_) & (A_ > 0)
    corr = float(np.corrcoef(np.log(A_[ok]), np.log(s_[ok]))[0, 1])
    med_ratio = float(np.median(A_[ok] / s_[ok]))
    med_AV = float(np.median(A_[ok] / V_[ok]))
    # regression log s ~ log A
    X = np.column_stack([np.ones(ok.sum()), np.log(A_[ok])])
    b, *_ = np.linalg.lstsq(X, np.log(s_[ok]), rcond=None)
    r2 = 1 - float(np.sum((np.log(s_[ok]) - X @ b) ** 2)) / float(np.sum((np.log(s_[ok]) - np.log(s_[ok]).mean()) ** 2))

    lines = ["# F13c — analytic ambiguous-mass A vs measured seam width s\n"]
    lines.append("| corpus | basin | A | s_real | A/s | A/V |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: r["A_over_s"]):
        lines.append(f"| {r['corpus']} | {r['basin']} | {r['A']:.0f} | {r['s_real']:.0f} | "
                     f"{r['A_over_s']:.2f} | {r['A_over_V']:.4f} |")
    lines.append("")
    lines.append(f"- corr(log A, log s_real) = {corr:.3f}  (n={int(ok.sum())})")
    lines.append(f"- median A/s = {med_ratio:.2f}; median A/V = {med_AV:.4f} (empirical s/V: 0.0118)")
    lines.append(f"- regression log s ~ log A: slope {b[1]:.3f}, R2 {r2:.3f}")
    (OUT / "f13c_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote f13c outputs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
