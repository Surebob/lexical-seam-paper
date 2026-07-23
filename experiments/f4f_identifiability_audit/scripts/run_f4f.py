"""F4f — identifiability audit for the mixture's tail parameters (sigma_T dial).

Trigger (f6b): narrow-tail and broad-tail PLN solutions fit a single-depth count
histogram comparably but predict cross-size behavior differently, making the f4d
"sigma_T ~ 1.6 universal" suspect in absolute terms.

Parts:
  A. MODEL-FREE dial: sd(log count | count >= 5) and IQR-based spread for the f8
     matched-size panel — does the cross-language ordering reproduce with no model?
  B. LANDSCAPE: for 5 English corpora, PLN fits from 16 diverse starts; keep all
     solutions with dNLL <= 5; record (sd_T, mu_T, n_total, NLL) and each
     solution's downscale prediction error (binomial-thin real counts to 1/4 and
     1/16; compare simulated V and c). Is the likelihood flat across basins? Does
     prediction single out a basin? Does NLL itself prefer the predictive basin?
  C. JOINT MULTI-DEPTH FIT: fit one parameter set to histograms at three depths
     simultaneously (scales 1, 1/4, 1/16 of the real counts via binomial thinning;
     model side is exact: rates scale). Report start-dispersion of sd_T under the
     joint fit (identified?) and the joint-fit sd_T ordering for 8 corpora
     (4 English + RU/ZH/FI/IT at matched 65k base).

Outputs: ../outputs/f4f_modelfree.csv, f4f_landscape.csv, f4f_joint.csv, f4f_summary.md
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
ML1 = REPO / "data" / "zipf_multilang"
LX = REPO / "data" / "langext"
F8CSV = REPO / "experiments" / "f8_matched_size_panel" / "outputs" / "f8_matched_panel.csv"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

EN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
SEED = 20260725
XGRID_N = 800

EN_CORPORA = [
    ("Shakespeare", DATA / "pg100.txt"), ("War and Peace", DATA / "pg2600.txt"),
    ("Moby Dick", DATA / "pg2701.txt"), ("Pride and Prejudice", DATA / "pg1342.txt"),
    ("Federalist Papers", DATA / "pg1404.txt"),
]
JOINT_PANEL = [
    ("EN Shakespeare@65k", DATA / "pg100.txt", "en"),
    ("EN Moby@65k", DATA / "pg2701.txt", "en"),
    ("EN P&P@65k", DATA / "pg1342.txt", "en"),
    ("EN Ulysses@65k", DATA / "pg4300.txt", "en"),
    ("RU WarPeace@65k", ML1 / "russian_war_and_peace" / "combined_clean.txt", "uni"),
    ("ZH ThreeKingdoms@65k", ML1 / "mandarin_three_kingdoms" / "combined_clean.txt", "jieba"),
    ("FI Seitseman@65k", LX / "finnish.txt", "uni"),
    ("IT Promessi@65k", LX / "italian.txt", "uni"),
]


def strip_g(text):
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
    return text[st:e]


def counts_from(path, mode, cap=None):
    text = strip_g(path.read_text(encoding="utf-8", errors="ignore"))
    if mode == "en":
        toks = EN_RE.findall(text.lower())
    elif mode == "jieba":
        import jieba
        toks = [t for t in jieba.cut(text) if UNI_RE.fullmatch(t)]
    else:
        toks = UNI_RE.findall(text.lower())
    if cap and len(toks) > cap:
        toks = toks[:cap]
    return np.array(sorted(Counter(toks).values(), reverse=True), dtype=np.int64)


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


# ---------- shared PLN machinery (latent-grid likelihood) ----------

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


def pln_nll_at_scales(theta, datasets):
    """datasets: list of (uniq, mult, scale). One parameter set, joint likelihood."""
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
    total = 0.0
    for uniq, mult, scale in datasets:
        both = np.concatenate([uniq, [0.0]])
        l_h = log_pn_grid(both, xs, lp_h, scale)
        l_t = log_pn_grid(both, xs, lp_t, scale)
        lp = np.logaddexp(math.log(pi_h) + l_h, math.log(1 - pi_h) + l_t)
        lp0 = min(math.exp(lp[-1]), 1 - 1e-12)
        total += -float(np.dot(mult, lp[:-1] - math.log1p(-lp0)))
    return total


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


def landscape_task(job):
    name, counts_list, x0 = job
    counts = np.asarray(counts_list, dtype=np.int64)
    uniq, mult = np.unique(counts.astype(float), return_counts=True)
    datasets = [(uniq, mult.astype(float), 1.0)]
    try:
        r = minimize(pln_nll_at_scales, x0=x0, args=(datasets,), method="Nelder-Mead",
                     options={"maxiter": 2600, "xatol": 1e-5, "fatol": 1e-7})
    except Exception:
        return None
    params = unpack(r.x)
    nt = n_total_of(params, len(counts))
    rng = np.random.default_rng(SEED + int(abs(r.fun)) % 99991)
    # prediction score: thin real counts vs simulate at 1/4 and 1/16
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
    return {"corpus": name, "nll": float(r.fun), "pi_h": params[0], "sd_h": params[2],
            "mu_t": params[3], "sd_t": params[4], "n_total": nt, "pred_err": float(err)}


def main():
    rng = np.random.default_rng(SEED)

    # ---------- Part A: model-free dial ----------
    a_rows = []
    f8 = {r["corpus"]: r for r in csv.DictReader(open(F8CSV, newline="", encoding="utf-8"))}
    jobsA = [("EN Shakespeare", DATA / "pg100.txt", "en"), ("EN War and Peace (tr)", DATA / "pg2600.txt", "en"),
             ("EN Moby Dick", DATA / "pg2701.txt", "en"), ("EN Pride and Prejudice", DATA / "pg1342.txt", "en"),
             ("EN Ulysses", DATA / "pg4300.txt", "en"), ("EN Federalist", DATA / "pg1404.txt", "en"),
             ("RU War and Peace", ML1 / "russian_war_and_peace" / "combined_clean.txt", "uni"),
             ("ZH Three Kingdoms", ML1 / "mandarin_three_kingdoms" / "combined_clean.txt", "jieba"),
             ("AR 1001 Nights", ML1 / "arabic_1001_nights" / "combined_clean.txt", "uni"),
             ("FI Seitseman veljesta", LX / "finnish.txt", "uni"),
             ("IT Promessi Sposi", LX / "italian.txt", "uni"),
             ("DE Wahlverwandtschaften", LX / "german.txt", "uni"),
             ("SV Gosta Berling", LX / "swedish.txt", "uni"),
             ("PT Dom Casmurro", LX / "portuguese.txt", "uni")]
    for name, path, mode in jobsA:
        counts = counts_from(path, mode, cap=65000)
        yc = np.log(counts[counts >= 5].astype(float))
        sd5 = float(np.std(yc, ddof=1))
        q = np.quantile(yc, [0.25, 0.75])
        f8row = f8.get(name) or f8.get(name.replace("EN ", "EN "))
        sdt_f8 = float(f8row["sd_t"]) if f8row else float("nan")
        a_rows.append({"corpus": name, "sd_logcount_ge5": sd5, "iqr_logcount_ge5": float(q[1] - q[0]),
                       "n_ge5": int(np.sum(counts >= 5)), "sd_t_f8": sdt_f8})
        print(f"A {name:26} sd(logN|N>=5)={sd5:.3f} f8_sd_T={sdt_f8:.3f}", flush=True)
    with open(OUT / "f4f_modelfree.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(a_rows[0].keys())); w.writeheader(); w.writerows(a_rows)
    ok = [r for r in a_rows if np.isfinite(r["sd_t_f8"])]
    corrA = float(np.corrcoef([r["sd_logcount_ge5"] for r in ok], [r["sd_t_f8"] for r in ok])[0, 1]) if len(ok) > 3 else float("nan")

    # ---------- Part B: landscape ----------
    jobsB = []
    for name, path in EN_CORPORA:
        counts = counts_from(path, "en")
        for i in range(16):
            x0 = [rng.uniform(-5, -1),                       # pi_h logit
                  rng.uniform(-6.5, 1.0),                    # mu_t
                  math.log(rng.uniform(0.6, 4.0)),           # sd_t
                  math.log(rng.uniform(2.0, 8.0)),           # mu gap
                  math.log(rng.uniform(0.4, 2.0))]           # sd_h
            jobsB.append((name, counts.tolist(), x0))
    b_rows = []
    with Pool(processes=12) as pool:
        for i, res in enumerate(pool.imap_unordered(landscape_task, jobsB, chunksize=1), 1):
            if res:
                b_rows.append(res)
                print(f"B [{i}/{len(jobsB)}] {res['corpus'][:18]:19} nll={res['nll']:.0f} sd_T={res['sd_t']:.2f} "
                      f"mu_T={res['mu_t']:+.2f} n_tot={res['n_total']:.0f} pred_err={res['pred_err']:.3f}", flush=True)
    with open(OUT / "f4f_landscape.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(b_rows[0].keys())); w.writeheader(); w.writerows(b_rows)

    # ---------- Part C: joint multi-depth fit ----------
    c_rows = []
    for name, path, mode in JOINT_PANEL:
        counts = counts_from(path, mode, cap=65000)
        sets = []
        for s in (1.0, 0.25, 0.0625):
            thin = counts if s == 1.0 else rng.binomial(counts, s)
            thin = thin[thin > 0].astype(float)
            uq, ml_ = np.unique(thin, return_counts=True)
            sets.append((uq, ml_.astype(float), s))
        sols = []
        for i in range(6):
            x0 = [rng.uniform(-5, -1), rng.uniform(-6, 1), math.log(rng.uniform(0.6, 3.5)),
                  math.log(rng.uniform(2.0, 8.0)), math.log(rng.uniform(0.4, 2.0))]
            try:
                r = minimize(pln_nll_at_scales, x0=x0, args=(sets,), method="Nelder-Mead",
                             options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-7})
                sols.append((float(r.fun), unpack(r.x)))
            except Exception:
                continue
        sols.sort(key=lambda t: t[0])
        best_nll = sols[0][0]
        near = [p for f_, p in sols if f_ - best_nll <= 2.0]
        sd_ts = [p[4] for p in near]
        c_rows.append({"corpus": name, "joint_sd_t": sols[0][1][4], "joint_mu_t": sols[0][1][3],
                       "joint_pi_h": sols[0][1][0], "n_near_solutions": len(near),
                       "sd_t_spread_near": float(max(sd_ts) - min(sd_ts)) if len(sd_ts) > 1 else 0.0,
                       "best_nll": best_nll})
        print(f"C {name:22} joint sd_T={sols[0][1][4]:.3f} (near-solution spread "
              f"{c_rows[-1]['sd_t_spread_near']:.3f}, n={len(near)})", flush=True)
    with open(OUT / "f4f_joint.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(c_rows[0].keys())); w.writeheader(); w.writerows(c_rows)

    # ---------- summary ----------
    lines = ["# F4f — identifiability audit\n"]
    lines.append(f"## A. Model-free dial: corr(sd(logN|N>=5), f8 sd_T) = {corrA:.3f}")
    for r in sorted(a_rows, key=lambda r: r["sd_logcount_ge5"]):
        lines.append(f"- {r['corpus']}: sd_logcount_ge5={r['sd_logcount_ge5']:.3f} (f8 sd_T={r['sd_t_f8']:.3f})")
    lines.append("\n## B. Landscape (per corpus: solutions within dNLL<=5 of best)")
    for name, _ in EN_CORPORA:
        sub = [r for r in b_rows if r["corpus"] == name]
        if not sub:
            continue
        best = min(r["nll"] for r in sub)
        near = [r for r in sub if r["nll"] - best <= 5.0]
        sd_rng = (min(r["sd_t"] for r in near), max(r["sd_t"] for r in near))
        nt_rng = (min(r["n_total"] for r in near), max(r["n_total"] for r in near))
        bestpred = min(near, key=lambda r: r["pred_err"])
        bestnll = min(near, key=lambda r: r["nll"])
        lines.append(f"- {name}: {len(near)} near-solutions; sd_T range {sd_rng[0]:.2f}-{sd_rng[1]:.2f}; "
                     f"n_total range {nt_rng[0]:.0f}-{nt_rng[1]:.0f}; "
                     f"best-NLL sd_T={bestnll['sd_t']:.2f} (pred_err {bestnll['pred_err']:.3f}); "
                     f"best-pred sd_T={bestpred['sd_t']:.2f} (pred_err {bestpred['pred_err']:.3f})")
    lines.append("\n## C. Joint multi-depth fits")
    for r in c_rows:
        lines.append(f"- {r['corpus']}: joint sd_T={r['joint_sd_t']:.3f}, near-solution spread {r['sd_t_spread_near']:.3f}")
    (OUT / "f4f_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
