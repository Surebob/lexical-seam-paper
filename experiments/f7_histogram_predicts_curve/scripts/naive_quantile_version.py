"""Does sigma_T predict the fitted Zipf-Mandelbrot exponent b?

Theory: for a lognormal type-rate population, the rank curve is
ln f(r) = mu - sigma * PhiInv(r/V); fitting ZM to that curve yields b = g(sigma),
a V-(nearly-)independent function of sigma alone. Compute g numerically, then test
b_pred = g(sigma_T) against each corpus's actually fitted zm_b (f1) using the
f4d-fitted sigma_T. Also handle the head: real corpora have the head club, so also
compute g on a two-component curve using each corpus's full PLN params.
"""
import csv
import math
from pathlib import Path

import numpy as np
from scipy.stats import norm

REPO = Path(r"C:\Users\Greg Kara\Desktop\temporary\lexical-seam")
F1 = REPO / "experiments" / "f1_fresh_reproduction" / "outputs" / "f1_per_corpus.csv"
F4D = REPO / "experiments" / "f4d_poisson_lognormal_mixture" / "outputs" / "f4d_mixture_fits.csv"

f1 = {r["corpus"]: r for r in csv.DictReader(open(F1, newline="", encoding="utf-8"))}
f4d = list(csv.DictReader(open(F4D, newline="", encoding="utf-8")))


def fit_zm_b(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 512)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(-coef[1]), float(c))
    return best[1], best[2]


def lognormal_curve(V, mu, sigma):
    q = (np.arange(1, V + 1) - 0.5) / V
    return np.exp(mu - sigma * norm.ppf(q))


def mixture_curve(V, pi_h, mu_h, sd_h, mu_t, sd_t):
    q = (np.arange(1, V + 1) - 0.5) / V
    lo = min(mu_h - 8 * sd_h, mu_t - 8 * sd_t)
    hi = max(mu_h + 8 * sd_h, mu_t + 8 * sd_t)
    xs = np.linspace(lo, hi, 4000)
    S = pi_h * norm.sf((xs - mu_h) / sd_h) + (1 - pi_h) * norm.sf((xs - mu_t) / sd_t)
    x_of_q = np.interp(q, S[::-1], xs[::-1])
    return np.exp(x_of_q)


print("=== g(sigma): ZM exponent fitted to a pure lognormal quantile curve (V=15000) ===")
grid = np.arange(0.8, 2.61, 0.1)
gvals = []
for s in grid:
    b, c = fit_zm_b(lognormal_curve(15000, 5.0, s))
    gvals.append(b)
    print(f"  sigma={s:.1f}  ->  b_ZM={b:.3f}  (c={c:.1f})")

print("\n=== per-corpus prediction ===")
print(f"{'corpus':32} {'sd_T':>5} {'b_pred(tail)':>12} {'b_pred(mix)':>11} {'b_actual':>8}")
preds_t, preds_m, actuals = [], [], []
for r in f4d:
    name = r["corpus"]
    sd_t = float(r["sd_t"])
    V = int(r["V"])
    b_actual = float(f1[name]["zm_b"])
    b_pred_tail = float(np.interp(sd_t, grid, gvals))
    curve = mixture_curve(V, float(r["pi_h"]), float(r["mu_h"]), float(r["sd_h"]),
                          float(r["mu_t"]), sd_t)
    b_pred_mix, _ = fit_zm_b(curve)
    preds_t.append(b_pred_tail)
    preds_m.append(b_pred_mix)
    actuals.append(b_actual)
    print(f"{name[:31]:32} {sd_t:5.2f} {b_pred_tail:12.3f} {b_pred_mix:11.3f} {b_actual:8.3f}")

pt, pm, ac = np.array(preds_t), np.array(preds_m), np.array(actuals)
print(f"\ncorr(b_pred_tail, b_actual) = {np.corrcoef(pt, ac)[0,1]:.4f}   median ratio {np.median(ac/pt):.3f}")
print(f"corr(b_pred_mix,  b_actual) = {np.corrcoef(pm, ac)[0,1]:.4f}   median ratio {np.median(ac/pm):.3f}")
print(f"median b_actual = {np.median(ac):.3f}; median sd_T = {np.median([float(r['sd_t']) for r in f4d]):.3f}")
print(f"implied sigma for b=1 (invert g): {float(np.interp(1.0, gvals, grid)):.3f}")
