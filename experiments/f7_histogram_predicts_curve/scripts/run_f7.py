"""sigma->b, corrected: through the observation model (Poisson + count>=1 floor).
Per corpus: simulate counts from the f4d PLN mixture (N_total via zero-truncation),
fit ZM to simulated counts, compare (b_sim, c_sim) with the corpus's actual (b, c).
"""
import csv
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import logsumexp

REPO = Path(r"C:\Users\Greg Kara\Desktop\temporary\lexical-seam")
F1 = REPO / "experiments" / "f1_fresh_reproduction" / "outputs" / "f1_per_corpus.csv"
F4D = REPO / "experiments" / "f4d_poisson_lognormal_mixture" / "outputs" / "f4d_mixture_fits.csv"

f1 = {r["corpus"]: r for r in csv.DictReader(open(F1, newline="", encoding="utf-8"))}
f4d = list(csv.DictReader(open(F4D, newline="", encoding="utf-8")))
GH_X, GH_W = hermgauss(48)
rng = np.random.default_rng(20260723)


def pln_p0(mu, sd):
    lam = np.exp(np.clip(mu + math.sqrt(2.0) * sd * GH_X, -30, 30))
    return float(np.exp(logsumexp(-lam + np.log(GH_W)) - 0.5 * math.log(math.pi)))


def fit_zm(freqs):
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


print(f"{'corpus':32} {'b_sim':>6} {'b_act':>6} {'c_sim':>8} {'c_act':>8} {'V_sim':>6} {'V_act':>6}")
bs, ba, cs, ca = [], [], [], []
for r in f4d:
    name = r["corpus"]
    pi_h, mu_h, sd_h = float(r["pi_h"]), float(r["mu_h"]), float(r["sd_h"])
    mu_t, sd_t = float(r["mu_t"]), float(r["sd_t"])
    V_obs = int(r["V"])
    p0 = pi_h * pln_p0(mu_h, sd_h) + (1 - pi_h) * pln_p0(mu_t, sd_t)
    n_total = int(round(V_obs / max(1e-9, 1 - p0)))
    b_reps, c_reps, v_reps = [], [], []
    for rep in range(2):
        comp = rng.random(n_total) < pi_h
        loglam = np.where(comp, rng.normal(mu_h, sd_h, n_total), rng.normal(mu_t, sd_t, n_total))
        counts = rng.poisson(np.exp(np.clip(loglam, -30, 30)))
        counts = counts[counts > 0]
        freqs = np.sort(counts)[::-1].astype(float)
        b, c = fit_zm(freqs)
        b_reps.append(b)
        c_reps.append(c)
        v_reps.append(len(freqs))
    b_sim, c_sim, v_sim = float(np.median(b_reps)), float(np.median(c_reps)), int(np.median(v_reps))
    b_act, c_act = float(f1[name]["zm_b"]), float(f1[name]["zm_c"])
    bs.append(b_sim); ba.append(b_act); cs.append(c_sim); ca.append(c_act)
    print(f"{name[:31]:32} {b_sim:6.3f} {b_act:6.3f} {c_sim:8.1f} {c_act:8.1f} {v_sim:6} {V_obs:6}")

bs, ba, cs, ca = map(np.array, (bs, ba, cs, ca))
print(f"\ncorr(b_sim, b_act) = {np.corrcoef(bs, ba)[0,1]:.4f}   median ratio b_act/b_sim = {np.median(ba/bs):.3f}")
print(f"corr(log1p c_sim, log1p c_act) = {np.corrcoef(np.log1p(cs), np.log1p(ca))[0,1]:.4f}   median ratio = {np.median((ca+1)/(cs+1)):.3f}")
print(f"median b_sim = {np.median(bs):.3f}  median b_act = {np.median(ba):.3f}")
