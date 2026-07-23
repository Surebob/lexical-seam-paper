"""F6b — heavy-tailed mixture and cross-size extrapolation.

Why: f6 showed the 2-lognormal (PLN) mixture predicts the c(T) collapse only
qualitatively — predicted c ~140 where thinned reality gives ~15. Hypothesis: the
latent rate distribution's LOWER tail is power-law (double Pareto-lognormal, Reed),
not lognormal; that tail controls how vocabulary appears/disappears with corpus
size, i.e. Heaps' law, and through it the apex-vs-tail balance that sets c.

Model NLT ("normal-minus-exponential tail component"):
  head:  x ~ N(mu_H, sd_H)                                  (unchanged)
  tail:  x = N(mu_T, sd_T) - Exp(rate beta)                 (heavy lower tail;
         density f(x) = beta*exp(beta*(x-mu)+beta^2 s^2/2)*Phi((mu-x)/s - beta*s))
  observed count ~ Poisson(e^x), zero-truncated.

Parts:
  A. Diagnostic: V(T) growth curves — binomially thinned REAL counts vs PLN sim
     vs NLT sim (fit both at full depth) for 4 exemplar corpora.
  B. Downward trajectories: c(T) and b(T) predicted by each model vs thinned real.
  C. Upward test: fit both models on a 150k binomial slice ONLY, predict full-size
     (V, b, c). The hard, honest prediction.

Outputs: ../outputs/f6b_fits.csv, f6b_trajectories.csv, f6b_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"),
    ("War and Peace", "pg2600.txt"),
    ("Moby Dick", "pg2701.txt"),
    ("Les Miserables", "pg135.txt"),
]
SCALES = [1.0, 0.5, 0.25, 0.125, 0.0625]
SEED = 20260724
rng_global = np.random.default_rng(SEED)

# latent-grid for NLT integrals
XGRID_N = 900


def counts_of(fname):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    st, e = 0, len(text)
    for m in SM:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            st = nl + 1 if nl != -1 else i
            break
    for m in EMK:
        i = text.find(m)
        if i != -1:
            e = i
            break
    return np.array(sorted(Counter(TOKEN_RE.findall(text[st:e].lower())).values(), reverse=True), dtype=np.int64)


def fit_zm(freqs):
    freqs = np.asarray(freqs, dtype=float)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 640)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(-coef[1]), float(c))
    return best[1], best[2]


def nlt_logpdf(x, mu, sd, beta):
    """log density of N(mu, sd) - Exp(beta)."""
    z = (mu - x) / sd - beta * sd
    return (math.log(beta) + beta * (x - mu) + 0.5 * beta * beta * sd * sd
            + norm.logcdf(z))


def make_grid(mu_t, sd_t, beta, mu_h, sd_h):
    lo = min(mu_t - 6 * sd_t - 18.0 / max(beta, 0.15), mu_h - 6 * sd_h)
    hi = max(mu_t + 6 * sd_t, mu_h + 6 * sd_h)
    return np.linspace(lo, hi, XGRID_N)


def comp_log_pn(uniq, xs, logpdf_xs):
    """log P(n) for counts uniq under latent-x density given on grid xs."""
    dx = xs[1] - xs[0]
    lam = np.exp(np.clip(xs, -35, 35))
    n = uniq[:, None]
    ll = n * xs[None, :] - lam[None, :] - gammaln(n + 1.0)
    return logsumexp(ll + logpdf_xs[None, :] + math.log(dx), axis=1)


def model_log_pn(uniq, params, kind):
    pi_h, mu_h, sd_h, mu_t, sd_t, beta = params
    if kind == "pln":
        xs = make_grid(mu_t, sd_t, 10.0, mu_h, sd_h)
        lp_t = norm.logpdf(xs, mu_t, sd_t)
    else:
        xs = make_grid(mu_t, sd_t, beta, mu_h, sd_h)
        lp_t = nlt_logpdf(xs, mu_t, sd_t, beta)
    lp_h = norm.logpdf(xs, mu_h, sd_h)
    both = np.concatenate([uniq, [0.0]])
    l_h = comp_log_pn(both, xs, lp_h)
    l_t = comp_log_pn(both, xs, lp_t)
    lp = np.logaddexp(math.log(pi_h) + l_h, math.log(1 - pi_h) + l_t)
    return lp[:-1], lp[-1]


def nll(theta, uniq, mult, kind):
    pi_h = 0.001 + 0.299 / (1.0 + math.exp(-theta[0]))
    mu_t = theta[1]
    sd_t = math.exp(theta[2])
    mu_h = mu_t + 0.5 + math.exp(theta[3])
    sd_h = math.exp(theta[4])
    beta = 0.1 + math.exp(theta[5]) if kind == "nlt" else 10.0
    if not (0.05 < sd_t < 6 and 0.05 < sd_h < 6 and beta < 25):
        return 1e12
    lp, lp0 = model_log_pn(uniq, (pi_h, mu_h, sd_h, mu_t, sd_t, beta), kind)
    if not np.all(np.isfinite(lp)):
        return 1e12
    return -float(np.dot(mult, lp - math.log1p(-min(math.exp(lp0), 1 - 1e-12))))


def unpack(theta, kind):
    pi_h = 0.001 + 0.299 / (1.0 + math.exp(-theta[0]))
    mu_t = float(theta[1])
    sd_t = float(math.exp(theta[2]))
    mu_h = mu_t + 0.5 + float(math.exp(theta[3]))
    sd_h = float(math.exp(theta[4]))
    beta = 0.1 + float(math.exp(theta[5])) if kind == "nlt" else 10.0
    return pi_h, mu_h, sd_h, mu_t, sd_t, beta


def fit_model(freqs, kind):
    uniq, mult = np.unique(freqs.astype(float), return_counts=True)
    mult = mult.astype(float)
    starts = [
        [-3.0, 0.3, math.log(1.2), math.log(4.0), math.log(0.8), math.log(0.7)],
        [-2.5, 1.0, math.log(0.9), math.log(5.0), math.log(0.6), math.log(1.5)],
        [-3.5, -0.5, math.log(1.6), math.log(4.5), math.log(1.0), math.log(0.4)],
    ]
    best = None
    for x0 in starts:
        try:
            r = minimize(nll, x0=x0, args=(uniq, mult, kind), method="Nelder-Mead",
                         options={"maxiter": 4000, "xatol": 1e-5, "fatol": 1e-7})
        except Exception:
            continue
        if best is None or r.fun < best.fun:
            best = r
    params = unpack(best.x, kind)
    # zero-truncation -> total latent types
    _, lp0 = model_log_pn(np.array([1.0]), params, kind)
    p0 = min(math.exp(lp0), 1 - 1e-12)
    n_total = len(freqs) / (1 - p0)
    return {"params": params, "nll": float(best.fun), "p0": p0, "n_total": n_total}


def simulate(fitres, kind, scale, rng):
    pi_h, mu_h, sd_h, mu_t, sd_t, beta = fitres["params"]
    n_total = int(round(fitres["n_total"]))
    comp = rng.random(n_total) < pi_h
    x = np.where(comp, rng.normal(mu_h, sd_h, n_total), rng.normal(mu_t, sd_t, n_total))
    if kind == "nlt":
        x = np.where(comp, x, x - rng.exponential(1.0 / beta, n_total))
    counts = rng.poisson(np.exp(np.clip(x, -35, 35)) * scale)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(float)


def traj(freqs_or_fit, kind, rng, is_fit=False):
    out = []
    for s in SCALES:
        if is_fit:
            reps = [simulate(freqs_or_fit, kind, s, rng) for _ in range(2)]
        else:
            reps = []
            for _ in range(2):
                thin = rng.binomial(freqs_or_fit, s)
                thin = np.sort(thin[thin > 0])[::-1].astype(float)
                reps.append(thin)
        Vs, bs, cs = [], [], []
        for fr in reps:
            b, c = fit_zm(fr)
            Vs.append(len(fr)); bs.append(b); cs.append(c)
        out.append({"scale": s, "V": float(np.median(Vs)), "b": float(np.median(bs)), "c": float(np.median(cs))})
    return out


def main():
    fit_rows, traj_rows = [], []
    lines = ["# F6b — heavy-tail mixture and cross-size extrapolation\n"]
    for name, fname in CORPORA:
        counts = counts_of(fname)
        T0 = int(counts.sum())
        rng = np.random.default_rng(SEED + hash(name) % 99991)
        print(f"=== {name} (T={T0}, V={len(counts)}) ===", flush=True)

        real = traj(counts, "real", rng)
        for r in real:
            traj_rows.append({"corpus": name, "model": "real_thinned", **r})

        fits = {}
        for kind in ["pln", "nlt"]:
            f = fit_model(counts, kind)
            fits[kind] = f
            pi_h, mu_h, sd_h, mu_t, sd_t, beta = f["params"]
            fit_rows.append({"corpus": name, "fit_depth": "full", "kind": kind,
                             "pi_h": pi_h, "mu_h": mu_h, "sd_h": sd_h, "mu_t": mu_t,
                             "sd_t": sd_t, "beta": beta, "nll": f["nll"],
                             "n_total": f["n_total"], "aic_like": 2 * f["nll"] + 2 * (5 if kind == "pln" else 6)})
            print(f"  {kind}: pi_H={100*pi_h:.2f}% sd_H={sd_h:.2f} mu_T={mu_t:.2f} sd_T={sd_t:.2f} "
                  f"beta={beta:.2f} nll={f['nll']:.0f} n_total={f['n_total']:.0f}", flush=True)
            for r in traj(f, kind, rng, is_fit=True):
                traj_rows.append({"corpus": name, "model": f"{kind}_fullfit", **r})

        # Part C: fit on a 150k binomial slice only, predict full size
        p150 = 150000 / T0
        slice_counts = rng.binomial(counts, p150)
        slice_counts = np.sort(slice_counts[slice_counts > 0])[::-1]
        b_full_real, c_full_real = fit_zm(counts.astype(float))
        for kind in ["pln", "nlt"]:
            fs = fit_model(slice_counts.astype(float), kind)
            # upscale: rates were fitted at slice depth; scale up by 1/p150
            ups = [simulate(fs, kind, 1.0 / p150, rng) for _ in range(2)]
            Vp = float(np.median([len(u) for u in ups]))
            bs, cs = zip(*[fit_zm(u) for u in ups])
            fit_rows.append({"corpus": name, "fit_depth": "150k_slice", "kind": kind,
                             **{k: v for k, v in zip(["pi_h", "mu_h", "sd_h", "mu_t", "sd_t", "beta"], fs["params"])},
                             "nll": fs["nll"], "n_total": fs["n_total"],
                             "aic_like": 2 * fs["nll"] + 2 * (5 if kind == "pln" else 6)})
            traj_rows.append({"corpus": name, "model": f"{kind}_upscaled_from_150k",
                              "scale": 1.0, "V": Vp, "b": float(np.median(bs)), "c": float(np.median(cs))})
            print(f"  UP {kind}: from 150k slice -> predict full: V={Vp:.0f} (real {len(counts)}), "
                  f"b={np.median(bs):.3f} (real {b_full_real:.3f}), c={np.median(cs):.1f} (real {c_full_real:.1f})", flush=True)

    with open(OUT / "f6b_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fit_rows[0].keys()))
        w.writeheader(); w.writerows(fit_rows)
    with open(OUT / "f6b_trajectories.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(traj_rows[0].keys()))
        w.writeheader(); w.writerows(traj_rows)

    # summary comparison
    lines.append("## Downward trajectories (median over reps)\n")
    lines.append("| corpus | scale | V real | V pln | V nlt | c real | c pln | c nlt |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for name, _ in CORPORA:
        for s in SCALES:
            g = lambda model, key: next((r[key] for r in traj_rows
                                         if r["corpus"] == name and r["model"] == model and abs(r["scale"] - s) < 1e-9), float("nan"))
            lines.append(f"| {name[:22]} | {s:g} | {g('real_thinned','V'):.0f} | {g('pln_fullfit','V'):.0f} | "
                         f"{g('nlt_fullfit','V'):.0f} | {g('real_thinned','c'):.1f} | {g('pln_fullfit','c'):.1f} | {g('nlt_fullfit','c'):.1f} |")
    lines.append("\n## Upward prediction (fit on 150k slice, predict full corpus)\n")
    for name, fname in CORPORA:
        lines.append(f"- {name}: see console rows in run log / f6b_trajectories.csv "
                     f"(models *_upscaled_from_150k vs real full values)")
    (OUT / "f6b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
