"""F13b — do the width-law failures heal when the mixture basin is selected by
an INDEPENDENT criterion (cross-size prediction) instead of raw ML?

f13 found the fitted PLN mixture reproduces the measured seam width almost
exactly on 16/25 corpora but fails by 5-600x on 9 — and the failures are the
deeply-sampled corpora whose ML basin has a large, blurred head (pi_h ~ 4-6%,
small mu-gap), the same population where f4f showed the latent mixture is
weakly identified and the ML solution is often not the best cross-size
predictor.

Test (non-circular): for each failure corpus, map the basin landscape (16
diverse Nelder-Mead starts on the full-histogram PLN likelihood, f4f
machinery), select the basin by MINIMUM CROSS-SIZE PREDICTION ERROR (thin the
real counts to 1/4 and 1/16, compare V and c against simulation — f4f's
pred_err, which never sees rank curves or widths), then simulate that basin
forward and fit the same 9-parameter erf-gate model as f13. If the
pred-selected basin's synthetic width lands near the real width, the width law
is contained in the generative account at the physical basin, and width
containment converges with cross-size prediction as a basin criterion.

Outputs: ../outputs/f13b_landscape.csv, f13b_summary.md
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
from scipy.optimize import least_squares, minimize
from scipy.special import erf as sp_erf, gammaln, logsumexp
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
ZIPF = REPO / "data" / "zipf"
F2 = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"

SEED = 20260724
XGRID_N = 900
N_STARTS_LANDSCAPE = 16
NLL_SLACK = 50.0          # basins admitted to selection: nll <= min_nll + slack
N_REPS = 3
N_STARTS_ERF = 24
MAX_NFEV = 9000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 5000.0, 10.0, 10.0])

# the nine f13 failure corpora (median s_synth/s_real > 2 under the ML basin)
FAILURES = [
    ("Complete Works of Shakespeare", "pg100.txt"),
    ("War and Peace", "pg2600.txt"),
    ("King James Bible", "pg10.txt"),
    ("Federalist Papers", "pg1404.txt"),
    ("Origin of Species", "pg1228.txt"),
    ("Wealth of Nations", "pg3300.txt"),
    ("Les Miserables", "pg135.txt"),
    ("Principia Ethica", "pg53430.txt"),
    ("Critique of Pure Reason", "pg4280.txt"),
]

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]


def counts_for(path: Path) -> np.ndarray:
    text = path.read_text(encoding="utf-8", errors="ignore")
    start, end = 0, len(text)
    for m in START_MARKERS:
        i = text.find(m)
        if i != -1:
            j = text.find("\n", i)
            start = j + 1 if j != -1 else i
            break
    for m in END_MARKERS:
        i = text.find(m)
        if i != -1:
            end = i
            break
    c = Counter(TOKEN_RE.findall(text[start:end].lower()))
    return np.array(sorted(c.values(), reverse=True), dtype=np.int64)


# ---------- PLN machinery (f4f) ----------

def grid_for(mu_t, sd_t, mu_h, sd_h):
    lo = min(mu_t - 7 * sd_t, mu_h - 6 * sd_h, -14.0)
    hi = max(mu_t + 7 * sd_t, mu_h + 6 * sd_h)
    return np.linspace(lo, hi, XGRID_N)


def log_pn_grid(uniq, xs, logpdf, scale=1.0):
    dx = xs[1] - xs[0]
    lam = np.exp(np.clip(xs, -35, 35)) * scale
    n = uniq[:, None]
    ll = n * np.log(np.maximum(lam[None, :], 1e-300)) - lam[None, :] - gammaln(n + 1.0)
    return logsumexp(ll + logpdf[None, :] + math.log(dx), axis=1)


def pln_nll(theta, uniq, mult):
    pi_h = 0.001 + 0.299 / (1.0 + math.exp(-theta[0]))
    mu_t = theta[1]
    sd_t = math.exp(theta[2])
    mu_h = mu_t + 0.5 + math.exp(theta[3])
    sd_h = math.exp(theta[4])
    if not (0.05 < sd_t < 5.0 and 0.05 < sd_h < 5.0 and -12 < mu_t < 10):
        return 1e12
    xs = grid_for(mu_t, sd_t, mu_h, sd_h)
    lp_t = norm.logpdf(xs, mu_t, sd_t)
    lp_h = norm.logpdf(xs, mu_h, sd_h)
    both = np.concatenate([uniq, [0.0]])
    l_h = log_pn_grid(both, xs, lp_h)
    l_t = log_pn_grid(both, xs, lp_t)
    lp = np.logaddexp(math.log(pi_h) + l_h, math.log(1 - pi_h) + l_t)
    lp0 = min(math.exp(lp[-1]), 1 - 1e-12)
    return -float(np.dot(mult, lp[:-1] - math.log1p(-lp0)))


def unpack(theta):
    pi_h = 0.001 + 0.299 / (1.0 + math.exp(-theta[0]))
    mu_t = float(theta[1])
    sd_t = float(math.exp(theta[2]))
    mu_h = mu_t + 0.5 + float(math.exp(theta[3]))
    sd_h = float(math.exp(theta[4]))
    return pi_h, mu_h, sd_h, mu_t, sd_t


def n_total_of(params, V_obs):
    pi_h, mu_h, sd_h, mu_t, sd_t = params
    xs = grid_for(mu_t, sd_t, mu_h, sd_h)
    lp_t = norm.logpdf(xs, mu_t, sd_t)
    lp_h = norm.logpdf(xs, mu_h, sd_h)
    l0 = np.logaddexp(math.log(pi_h) + log_pn_grid(np.array([0.0]), xs, lp_h),
                      math.log(1 - pi_h) + log_pn_grid(np.array([0.0]), xs, lp_t))[0]
    p0 = min(math.exp(l0), 1 - 1e-12)
    return V_obs / (1 - p0)


def simulate(params, n_total, scale, rng):
    pi_h, mu_h, sd_h, mu_t, sd_t = params
    n = int(round(n_total))
    comp = rng.random(n) < pi_h
    x = np.where(comp, rng.normal(mu_h, sd_h, n), rng.normal(mu_t, sd_t, n))
    counts = rng.poisson(np.exp(np.clip(x, -35, 35)) * scale)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(float)


def fit_zm_c(freqs):
    freqs = np.asarray(freqs, dtype=float)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 512)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(c))
    return best[1]


def landscape_task(job):
    name, counts_list, x0 = job
    counts = np.asarray(counts_list, dtype=np.int64)
    uniq, mult = np.unique(counts.astype(float), return_counts=True)
    try:
        r = minimize(pln_nll, x0=x0, args=(uniq, mult.astype(float)), method="Nelder-Mead",
                     options={"maxiter": 2600, "xatol": 1e-5, "fatol": 1e-7})
    except Exception:
        return None
    params = unpack(r.x)
    nt = n_total_of(params, len(counts))
    rng = np.random.default_rng(SEED + int(abs(r.fun)) % 99991)
    err = 0.0
    for s in (0.25, 0.0625):
        thin = rng.binomial(counts, s)
        thin = np.sort(thin[thin > 0])[::-1].astype(float)
        sim = simulate(params, nt, s, rng)
        cV = abs(math.log(max(len(sim), 1) / len(thin)))
        c_r = fit_zm_c(thin)
        c_s = fit_zm_c(sim)
        cC = abs(math.log1p(c_s) - math.log1p(c_r))
        err += cV + 0.5 * cC
    return {"corpus": name, "nll": float(r.fun), "pi_h": params[0], "mu_h": params[1],
            "sd_h": params[2], "mu_t": params[3], "sd_t": params[4],
            "n_total": nt, "pred_err": float(err)}


# ---------- erf-gate fitter (f13/f5b) ----------

def predict(ranks, log_ranks, params):
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    z = (log_ranks - math.log(k)) / w_gate
    sigma = 0.5 * (1.0 - sp_erf(z))
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


def erf_task(job):
    name, rep, freqs_list = job
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + (hash((name, rep)) % (2**31)))
    anchor6 = piecewise_anchor(ranks, y)

    def resid(p):
        return predict(ranks, log_ranks, p) - y

    best_rmse, best_p = None, None
    for i in range(N_STARTS_ERF):
        if i < N_STARTS_ERF // 2:
            x0 = np.array(anchor6 + [0.0, 0.0, 0.0])
            x0[:6] += rng.normal(0.0, [6.0, 0.3, 90.0, 6.0, 0.3, 90.0])
            x0[6] = min(0.02 * V, 4000.0) * math.exp(rng.normal(0.0, 0.7))
            x0[7] = 0.5 * math.exp(rng.normal(0.0, 0.45))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.45))
        else:
            x0 = np.array([
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(5000))),
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
    return {"corpus": name, "rep": rep, "V_synth": V,
            "s_synth": float(best_p[6] * best_p[8]), "rmse": best_rmse}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    rng = np.random.default_rng(SEED)

    print("tokenizing failure corpora...", flush=True)
    corpus_counts = {}
    for name, fn in FAILURES:
        corpus_counts[name] = counts_for(ZIPF / fn)
        print(f"  {name[:36]:36} V={len(corpus_counts[name]):6,}", flush=True)

    # ---- landscape ----
    jobs = []
    for name, _ in FAILURES:
        for _ in range(N_STARTS_LANDSCAPE):
            x0 = [rng.uniform(-5, -1),
                  rng.uniform(-6.5, 1.0),
                  math.log(rng.uniform(0.6, 4.0)),
                  math.log(rng.uniform(2.0, 8.0)),
                  math.log(rng.uniform(0.4, 2.0))]
            jobs.append((name, corpus_counts[name].tolist(), x0))
    print(f"{len(jobs)} landscape fits queued", flush=True)
    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(landscape_task, jobs, chunksize=1), 1):
            if res:
                rows.append(res)
                print(f"[{i}/{len(jobs)}] {res['corpus'][:24]:24} nll={res['nll']:.0f} "
                      f"pi={res['pi_h']:.4f} gap={res['mu_h']-res['mu_t']:.2f} "
                      f"pred_err={res['pred_err']:.3f}", flush=True)
    with open(OUT / "f13b_landscape.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # ---- basin selection + resimulation ----
    erf_jobs = []
    chosen = {}
    for name, _ in FAILURES:
        sub = [r for r in rows if r["corpus"] == name]
        if not sub:
            continue
        min_nll = min(r["nll"] for r in sub)
        adm = [r for r in sub if r["nll"] <= min_nll + NLL_SLACK]
        pick = min(adm, key=lambda r: r["pred_err"])
        ml = min(sub, key=lambda r: r["nll"])
        chosen[name] = {"pick": pick, "ml": ml, "min_nll": min_nll}
        params = (pick["pi_h"], pick["mu_h"], pick["sd_h"], pick["mu_t"], pick["sd_t"])
        nt = pick["n_total"]
        for rep in range(N_REPS):
            freqs = simulate(params, nt, 1.0, rng)
            erf_jobs.append((name, rep, freqs.tolist()))
    print(f"\n{len(erf_jobs)} erf fits queued on pred-selected basins", flush=True)
    efits = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, erf_jobs, chunksize=1), 1):
            efits.append(res)
            print(f"[{i}/{len(erf_jobs)}] {res['corpus'][:24]:24} rep{res['rep']} "
                  f"V={res['V_synth']:6} s={res['s_synth']:8.1f}", flush=True)

    # ---- summary ----
    ml_ratio_f13 = {  # from f13_summary.md (ML basin, for contrast)
        "Complete Works of Shakespeare": 15.09, "War and Peace": 18.19,
        "King James Bible": 18.21, "Federalist Papers": 6.05,
        "Origin of Species": 4.60, "Wealth of Nations": 4.98,
        "Les Miserables": 15.94, "Principia Ethica": 600.64,
        "Critique of Pure Reason": 624.78,
    }
    lines = ["# F13b — basin reselection by cross-size prediction: do the widths heal?\n"]
    lines.append("| corpus | ML ratio (f13) | pred-basin ratio | pred pi_h | pred gap | dNLL(pred-ML) | pred_err pred/ML |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    healed = 0
    ratios = []
    for name, _ in FAILURES:
        if name not in chosen:
            continue
        pick, ml = chosen[name]["pick"], chosen[name]["ml"]
        sub = [e for e in efits if e["corpus"] == name]
        s_med = float(np.median([e["s_synth"] for e in sub]))
        s_real = float(f2[name]["s_at_min"])
        ratio = s_med / s_real
        ratios.append(ratio)
        if 0.5 <= ratio <= 2.0:
            healed += 1
        lines.append(f"| {name} | {ml_ratio_f13.get(name, float('nan')):.1f} | {ratio:.2f} | "
                     f"{pick['pi_h']:.4f} | {pick['mu_h']-pick['mu_t']:.2f} | "
                     f"{pick['nll']-ml['nll']:.1f} | {pick['pred_err']:.2f}/{ml['pred_err']:.2f} |")
    lines.append("")
    lines.append(f"- healed (ratio in [0.5, 2]): {healed}/{len(ratios)}")
    lines.append(f"- median pred-basin ratio: {float(np.median(ratios)):.2f} "
                 f"(ML-basin median over these 9: {float(np.median(list(ml_ratio_f13.values()))):.1f})")
    lines.append("")
    if healed >= 6:
        lines.append("**Reading: the width law IS contained in the generative account at the "
                     "basin selected by cross-size prediction — two independent observables "
                     "(cross-size V/c and seam width) converge on the same basin, and the "
                     "asymptotic derivation of the 1.2% constant is justified.**")
    else:
        lines.append("**Reading: basin reselection does not (fully) heal the widths — the "
                     "width law requires an ingredient beyond the single-depth PLN mixture "
                     "under any admissible basin.**")
    (OUT / "f13b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
