"""F13 diagnostic — what separates the corpora where the mixture reproduces the
width from those where it fails, and does an analytic overlap predictor exist?

A. Failure correlates: s_synth/s_real ratio vs sampling depth (tokens/V), vs
   fitted ZM c, vs mixture params (pi_h, mu_h - mu_t gap, sd_t).
B. Analytic seed: for each corpus compute, from the f4d mixture alone, the
   "ambiguous mass" A = expected number of types whose latent log-rate lies
   between the two points where pi_h*phi_H = (1-pi_h)*phi_T (the density
   crossings). Compare A to the measured s. If A tracks s on the 16 successes,
   the derivation target becomes: show A/V ~ 0.012 from the fitted parameter
   regularities.

Output: ../outputs/f13_diag.md
"""
from __future__ import annotations

import csv
import math
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
F13 = OUT / "f13_width_fits.csv"

GH_X, GH_W = hermgauss(48)


def pln_p0(mu, sd):
    lam = np.exp(np.clip(mu + math.sqrt(2.0) * sd * GH_X, -30, 30))
    return float(np.exp(logsumexp(-lam + np.log(GH_W)) - 0.5 * math.log(math.pi)))


def crossings(pi_h, mu_h, sd_h, mu_t, sd_t):
    """log-rate points where head and tail component densities are equal."""
    xs = np.linspace(min(mu_t - 6 * sd_t, mu_h - 6 * sd_h),
                     max(mu_t + 6 * sd_t, mu_h + 6 * sd_h), 20001)
    d = (math.log(pi_h) + norm.logpdf(xs, mu_h, sd_h)) - \
        (math.log1p(-pi_h) + norm.logpdf(xs, mu_t, sd_t))
    sign = np.sign(d)
    idx = np.where(np.diff(sign) != 0)[0]
    return [float(xs[i]) for i in idx]


def main():
    f4d = list(csv.DictReader(open(F4D, newline="", encoding="utf-8")))
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    f13rows = list(csv.DictReader(open(F13, newline="", encoding="utf-8")))

    lines = ["# F13 diagnostic\n"]
    lines.append("| corpus | ratio | tok/V | pi_h | mu_h-mu_t | sd_t | crossings | A (amb. types) | s_real | A/s_real | A/V |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

    ratios, depths, gaps, sds, pis = [], [], [], [], []
    As, Ss, Vs = [], [], []
    for r in f4d:
        name = r["corpus"]
        pi_h, mu_h, sd_h = float(r["pi_h"]), float(r["mu_h"]), float(r["sd_h"])
        mu_t, sd_t = float(r["mu_t"]), float(r["sd_t"])
        V = float(r["V"]); toks = float(r["tokens"])
        sub = [float(x["s_synth"]) for x in f13rows if x["corpus"] == name]
        s_synth = float(np.median(sub)) if sub else float("nan")
        s_real = float(f2[name]["s_at_min"])
        ratio = s_synth / s_real

        # ambiguous mass between density crossings, as expected OBSERVED types
        cr = crossings(pi_h, mu_h, sd_h, mu_t, sd_t)
        if len(cr) >= 2:
            lo, hi = cr[-2], cr[-1]  # the crossing pair bracketing the handover
        elif len(cr) == 1:
            # single crossing: measure where the minority share is in [10%, 90%]
            lo = hi = cr[0]
        p0 = pi_h * pln_p0(mu_h, sd_h) + (1 - pi_h) * pln_p0(mu_t, sd_t)
        n_total = V / max(1e-9, 1 - p0)
        if len(cr) >= 1 and hi > lo:
            mass = (pi_h * (norm.cdf(hi, mu_h, sd_h) - norm.cdf(lo, mu_h, sd_h)) +
                    (1 - pi_h) * (norm.cdf(hi, mu_t, sd_t) - norm.cdf(lo, mu_t, sd_t)))
        else:
            # width of the share transition: share_H(x) in [0.1, 0.9]
            xs = np.linspace(mu_t - 6 * sd_t, mu_h + 6 * sd_h, 20001)
            lg = (math.log(pi_h) + norm.logpdf(xs, mu_h, sd_h)) - \
                 (math.log1p(-pi_h) + norm.logpdf(xs, mu_t, sd_t))
            share = 1.0 / (1.0 + np.exp(-lg))
            band = (share > 0.1) & (share < 0.9)
            if band.any():
                lo, hi = float(xs[band].min()), float(xs[band].max())
                mass = (pi_h * (norm.cdf(hi, mu_h, sd_h) - norm.cdf(lo, mu_h, sd_h)) +
                        (1 - pi_h) * (norm.cdf(hi, mu_t, sd_t) - norm.cdf(lo, mu_t, sd_t)))
            else:
                mass = float("nan")
        # only count OBSERVABLE types (weight by 1 - P(n=0 | lambda)) -- approximate
        # with the observation probability at the band centre
        mid = 0.5 * (lo + hi)
        p_obs_mid = 1.0 - math.exp(-math.exp(min(mid, 30)))
        A = n_total * mass * p_obs_mid

        ratios.append(ratio); depths.append(toks / V); gaps.append(mu_h - mu_t)
        sds.append(sd_t); pis.append(pi_h)
        if np.isfinite(A) and A > 0:
            As.append(A); Ss.append(s_real); Vs.append(V)

        lines.append(f"| {name[:28]} | {ratio:.2f} | {toks/V:.1f} | {pi_h:.4f} | "
                     f"{mu_h-mu_t:.2f} | {sd_t:.2f} | {len(cr)} | {A:.0f} | {s_real:.0f} | "
                     f"{A/s_real:.2f} | {A/V:.4f} |")

    ratios = np.array(ratios); depths = np.array(depths); gaps = np.array(gaps)
    sds = np.array(sds); pis = np.array(pis)
    lr = np.log(ratios)
    lines.append("")
    lines.append("## Failure correlates (corr with log width-ratio)")
    for nm, v in [("log tokens/V", np.log(depths)), ("mu_h - mu_t", gaps),
                  ("sd_t", sds), ("pi_h", pis)]:
        lines.append(f"- corr(log ratio, {nm}) = {float(np.corrcoef(lr, v)[0,1]):+.3f}")

    if len(As) > 3:
        As_, Ss_, Vs_ = map(np.array, (As, Ss, Vs))
        lines.append("")
        lines.append("## Analytic overlap predictor A vs measured s")
        lines.append(f"- corr(log A, log s_real) = {float(np.corrcoef(np.log(As_), np.log(Ss_))[0,1]):.3f}")
        lines.append(f"- median A/s_real = {float(np.median(As_/Ss_)):.3f}")
        lines.append(f"- median A/V = {float(np.median(As_/Vs_)):.4f}  (target 0.0118)")

    (OUT / "f13_diag.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote f13_diag.md")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
