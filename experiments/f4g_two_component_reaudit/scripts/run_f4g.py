"""F4g — re-audit of "two components demanded 25/25" (f4d) under broad-basin fitting.

f4f showed the latent PLN mixture is weakly identified from one depth and that
f4d's restricted starts found a narrow basin. Question: does the 2-vs-1-component
BIC preference survive when BOTH models get broad multi-start optimization?
(If a broad-tailed single lognormal matches the 2-component fit, v6's "direct
distributional two-population evidence" must be retired or reframed.)

10 corpora x {1-comp, 2-comp} x 8 diverse starts, latent-grid zero-truncated MLE.
Outputs: ../outputs/f4g_results.csv, f4g_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from multiprocessing import Pool
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
SEED = 20260726
XGRID_N = 640

CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"), ("War and Peace", "pg2600.txt"),
    ("Moby Dick", "pg2701.txt"), ("King James Bible", "pg10.txt"),
    ("Federalist Papers", "pg1404.txt"), ("Pride and Prejudice", "pg1342.txt"),
    ("Ulysses", "pg4300.txt"), ("Dubliners", "pg2814.txt"),
    ("Origin of Species", "pg1228.txt"), ("Grimm's Fairy Tales", "pg2591.txt"),
]


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
    return np.array(sorted(Counter(TOKEN_RE.findall(text[st:e].lower())).values(), reverse=True), dtype=np.float64)


def log_pn_grid(uniq, xs, logpdf):
    dx = xs[1] - xs[0]
    lam = np.exp(np.clip(xs, -35, 35))
    n = uniq[:, None]
    ll = n * xs[None, :] - lam[None, :] - gammaln(n + 1.0)
    return logsumexp(ll + logpdf[None, :] + math.log(dx), axis=1)


def nll_1comp(theta, uniq, mult):
    mu, sd = theta[0], math.exp(theta[1])
    if not (0.05 < sd < 6.0 and -14 < mu < 10):
        return 1e12
    xs = np.linspace(min(mu - 8 * sd, -16), mu + 8 * sd, XGRID_N)
    lp = norm.logpdf(xs, mu, sd)
    both = np.concatenate([uniq, [0.0]])
    l = log_pn_grid(both, xs, lp)
    p0 = min(math.exp(l[-1]), 1 - 1e-12)
    return -float(np.dot(mult, l[:-1] - math.log1p(-p0)))


def nll_2comp(theta, uniq, mult):
    pi_h = 0.001 + 0.299 / (1.0 + math.exp(-theta[0]))
    mu_t = theta[1]
    sd_t = math.exp(theta[2])
    mu_h = mu_t + 0.5 + math.exp(theta[3])
    sd_h = math.exp(theta[4])
    if not (0.05 < sd_t < 6.0 and 0.05 < sd_h < 6.0 and -14 < mu_t < 10):
        return 1e12
    lo = min(mu_t - 8 * sd_t, mu_h - 6 * sd_h, -16)
    hi = max(mu_t + 8 * sd_t, mu_h + 6 * sd_h)
    xs = np.linspace(lo, hi, XGRID_N)
    lp_t = norm.logpdf(xs, mu_t, sd_t)
    lp_h = norm.logpdf(xs, mu_h, sd_h)
    both = np.concatenate([uniq, [0.0]])
    l_t = log_pn_grid(both, xs, lp_t)
    l_h = log_pn_grid(both, xs, lp_h)
    lp = np.logaddexp(math.log(pi_h) + l_h, math.log(1 - pi_h) + l_t)
    p0 = min(math.exp(lp[-1]), 1 - 1e-12)
    return -float(np.dot(mult, lp[:-1] - math.log1p(-p0)))


def fit_task(job):
    name, counts_list, model, x0 = job
    counts = np.asarray(counts_list, dtype=np.float64)
    uniq, mult = np.unique(counts, return_counts=True)
    mult = mult.astype(float)
    fn = nll_1comp if model == "one" else nll_2comp
    try:
        r = minimize(fn, x0=x0, args=(uniq, mult), method="Nelder-Mead",
                     options={"maxiter": 2400, "xatol": 1e-5, "fatol": 1e-7})
    except Exception:
        return None
    return {"corpus": name, "model": model, "nll": float(r.fun),
            "theta": ";".join(f"{v:.4f}" for v in r.x)}


def main():
    rng = np.random.default_rng(SEED)
    jobs = []
    Vs = {}
    for name, fname in CORPORA:
        counts = counts_of(fname)
        Vs[name] = len(counts)
        for i in range(8):
            jobs.append((name, counts.tolist(), "one",
                         [rng.uniform(-7, 2), math.log(rng.uniform(0.6, 4.5))]))
            jobs.append((name, counts.tolist(), "two",
                         [rng.uniform(-5, -1), rng.uniform(-7, 1.5), math.log(rng.uniform(0.6, 4.0)),
                          math.log(rng.uniform(2.0, 9.0)), math.log(rng.uniform(0.4, 2.2))]))
    print(f"{len(jobs)} fits", flush=True)
    rows = []
    with Pool(processes=12) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, jobs, chunksize=1), 1):
            if res:
                rows.append(res)
                if i % 10 == 0:
                    print(f"[{i}/{len(jobs)}]", flush=True)

    with open(OUT / "f4g_results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    lines = ["# F4g — two-component demand, re-audited under broad-basin fitting\n"]
    lines.append("| corpus | best NLL 1-comp | best NLL 2-comp | ΔNLL (1−2) | ΔBIC (1−2) | 2-comp wins |")
    lines.append("|---|---:|---:|---:|---:|:---:|")
    wins = 0
    for name, _ in CORPORA:
        n1 = min((r["nll"] for r in rows if r["corpus"] == name and r["model"] == "one"), default=float("nan"))
        n2 = min((r["nll"] for r in rows if r["corpus"] == name and r["model"] == "two"), default=float("nan"))
        V = Vs[name]
        bic1 = 2 * math.log(V) + 2 * n1
        bic2 = 5 * math.log(V) + 2 * n2
        win = bic2 < bic1
        wins += int(win)
        lines.append(f"| {name} | {n1:.1f} | {n2:.1f} | {n1-n2:.1f} | {bic1-bic2:.1f} | {'YES' if win else 'no'} |")
    lines.append(f"\n**2-component wins under broad-basin fitting: {wins}/{len(CORPORA)}**")
    (OUT / "f4g_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
