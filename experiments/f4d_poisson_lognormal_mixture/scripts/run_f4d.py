"""F4d — Proper 2-component Poisson-lognormal mixture fits (all 25 English corpora)
and downstream re-tests.

Fixes f4b's degenerate GMM (which put one component on the hapax spike): here the
latent type-rate distribution is a 2-component lognormal mixture observed through
Poisson sampling with zero-truncation, fit by MLE (Gauss-Hermite quadrature over the
latent rate, aggregated over unique counts). Head component constrained to be the
minority high-rate population (pi_H <= 0.30, mu_H >= mu_T + 0.5).

Downstream:
  A. Corrected population parameters (pi_H, mu_H, sd_H, mu_T, sd_T) per corpus.
  B. P1 redo: w_gate ~ sqrt(2)*beta*sd_H with the corrected sd_H; factor analysis
     including the sqrt(e) candidate and e^(sd^2/2) variants.
  C. Hunch-1 redo: s (f2) vs corrected head size N_H = pi_H*V (non-circular now).

Outputs: ../outputs/f4d_mixture_fits.csv, f4d_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
F2 = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"
ARCH = Path(r"C:\Users\Greg Kara\Desktop\temporary\emlexperiment\results\s2_v3_windows_full_outputs_2026-04-18\s2_v3_per_fit_results.csv")

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"), ("War and Peace", "pg2600.txt"),
    ("Moby Dick", "pg2701.txt"), ("King James Bible", "pg10.txt"),
    ("Federalist Papers", "pg1404.txt"), ("Grimm's Fairy Tales", "pg2591.txt"),
    ("Don Quixote", "pg996.txt"), ("Pride and Prejudice", "pg1342.txt"),
    ("Canterbury Tales", "pg2383.txt"), ("Arabian Nights (Vol 1)", "pg3435.txt"),
    ("Aesop's Fables", "pg11339.txt"), ("Complete Sherlock Holmes", "pg1661.txt"),
    ("Jane Eyre", "pg1260.txt"), ("Dubliners", "pg2814.txt"),
    ("The Iliad", "pg6130.txt"), ("Democracy in America", "pg815.txt"),
    ("Origin of Species", "pg1228.txt"), ("Wealth of Nations", "pg3300.txt"),
    ("Les Miserables", "pg135.txt"), ("Decline and Fall Vol 1", "pg731.txt"),
    ("Emile", "pg5427.txt"), ("Ulysses", "pg4300.txt"),
    ("Collected Poe", "pg2147.txt"), ("Principia Ethica", "pg53430.txt"),
    ("Critique of Pure Reason", "pg4280.txt"),
]

GH_X, GH_W = hermgauss(48)
LOG_GH_W = np.log(GH_W)


def load_freqs(fname):
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
    counts = Counter(TOKEN_RE.findall(text[st:e].lower()))
    return np.array(sorted(counts.values(), reverse=True), dtype=np.float64)


def pln_log_pn(ns, mu, sd):
    """log P(N=n) under Poisson-lognormal via Gauss-Hermite. ns: array of counts."""
    loglam = mu + math.sqrt(2.0) * sd * GH_X          # (Q,)
    lam = np.exp(np.clip(loglam, -30, 30))
    # log Poisson(n | lam): n*loglam - lam - gammaln(n+1); shape (len(ns), Q)
    n = ns[:, None]
    ll = n * loglam[None, :] - lam[None, :] - gammaln(n + 1.0)
    return logsumexp(ll + LOG_GH_W[None, :], axis=1) - 0.5 * math.log(math.pi)


def mixture_nll(theta, uniq, mult):
    lpi = theta[0]
    pi_h = 1.0 / (1.0 + math.exp(-lpi))
    pi_h = 0.001 + (0.30 - 0.001) * pi_h            # head share in [0.001, 0.30]
    mu_t = theta[1]
    sd_t = math.exp(theta[2])
    dmu = 0.5 + math.exp(theta[3])                   # mu_h = mu_t + dmu >= +0.5
    mu_h = mu_t + dmu
    sd_h = math.exp(theta[4])
    if not (0.03 < sd_t < 6.0 and 0.03 < sd_h < 6.0):
        return 1e12
    lp_h = pln_log_pn(uniq, mu_h, sd_h)
    lp_t = pln_log_pn(uniq, mu_t, sd_t)
    lp_mix = np.logaddexp(math.log(pi_h) + lp_h, math.log(1 - pi_h) + lp_t)
    z = np.array([0.0])
    lp0 = np.logaddexp(math.log(pi_h) + pln_log_pn(z, mu_h, sd_h),
                       math.log(1 - pi_h) + pln_log_pn(z, mu_t, sd_t))[0]
    log_trunc = math.log1p(-min(math.exp(lp0), 1 - 1e-12))
    return -float(np.dot(mult, lp_mix - log_trunc))


def fit_pln_mixture(freqs):
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
    pi_h = 0.001 + (0.30 - 0.001) / (1.0 + math.exp(-t[0]))
    mu_t, sd_t = float(t[1]), float(math.exp(t[2]))
    mu_h = mu_t + 0.5 + math.exp(t[3])
    sd_h = float(math.exp(t[4]))
    # single-component PLN for comparison (is 2 components even demanded?)
    def nll1(th):
        lp = pln_log_pn(uniq, th[0], math.exp(th[1]))
        lp0 = pln_log_pn(np.array([0.0]), th[0], math.exp(th[1]))[0]
        return -float(np.dot(mult, lp - math.log1p(-min(math.exp(lp0), 1 - 1e-12))))
    b1 = None
    for x0 in [[1.0, math.log(1.5)], [0.0, math.log(2.0)]]:
        r = minimize(nll1, x0=x0, method="Nelder-Mead", options={"maxiter": 2000})
        if b1 is None or r.fun < b1.fun:
            b1 = r
    V = int(len(freqs))
    bic2 = 5 * math.log(V) + 2 * best.fun
    bic1 = 2 * math.log(V) + 2 * b1.fun
    return {"pi_h": pi_h, "mu_h": mu_h, "sd_h": sd_h, "mu_t": mu_t, "sd_t": sd_t,
            "nll2": float(best.fun), "nll1": float(b1.fun),
            "bic2": bic2, "bic1": bic1, "two_comp_wins": bic2 < bic1}


def main():
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    erf_fit = {r["corpus"]: r for r in csv.DictReader(open(ARCH, newline="", encoding="utf-8")) if r["gate"] == "erf"}

    rows = []
    for name, fname in CORPORA:
        freqs = load_freqs(fname)
        V = len(freqs)
        fit = fit_pln_mixture(freqs)
        k = float(erf_fit[name]["k"])
        w_fit = float(erf_fit[name]["w_gate"])
        ranks = np.arange(1, V + 1, dtype=float)
        lo, hi = int(max(1, k / 2)), int(min(V, k * 2))
        lr, lf = np.log(ranks[lo - 1:hi]), np.log(freqs[lo - 1:hi])
        A = np.column_stack([np.ones_like(lf), lf])
        beta = abs(float(np.linalg.lstsq(A, lr, rcond=None)[0][1]))
        w_pred = math.sqrt(2.0) * beta * fit["sd_h"]
        s = float(f2[name]["s_at_min"])
        rows.append({"corpus": name, "V": V, "tokens": int(freqs.sum()),
                     **{k2: (round(v, 6) if isinstance(v, float) else v) for k2, v in fit.items()},
                     "N_H": fit["pi_h"] * V, "beta": beta, "w_fit": w_fit,
                     "w_pred": w_pred, "ratio": w_fit / w_pred, "s_f2": s,
                     "s_over_NH": s / (fit["pi_h"] * V)})
        r = rows[-1]
        print(f"{name[:30]:31} pi_H={100*r['pi_h']:5.2f}% mu_H={r['mu_h']:6.2f} sd_H={r['sd_h']:5.2f} "
              f"sd_T={r['sd_t']:5.2f} 2comp={'Y' if r['two_comp_wins'] else 'N'} "
              f"N_H={r['N_H']:7.0f} s/N_H={r['s_over_NH']:6.3f} ratio={r['ratio']:5.2f}", flush=True)

    with open(OUT / "f4d_mixture_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    def arr(k2):
        return np.array([r[k2] for r in rows], dtype=float)

    def cc(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    ratio = arr("ratio")
    sdh, sdt = arr("sd_h"), arr("sd_t")
    s_vals, nh = arr("s_f2"), arr("N_H")
    notdub = [i for i, r in enumerate(rows) if r["corpus"] != "Dubliners"]
    lines = ["# F4d — Poisson-lognormal mixture fits and re-tests\n"]
    lines.append(f"- two-component demanded by BIC: {sum(1 for r in rows if r['two_comp_wins'])}/25")
    lines.append(f"- head share pi_H: median {float(np.median(arr('pi_h')))*100:.2f}% "
                 f"(range {arr('pi_h').min()*100:.2f}-{arr('pi_h').max()*100:.2f}%)")
    lines.append(f"\n## P1 redo (corrected sd_H)")
    lines.append(f"- corr(w_pred, w_fit) = {cc(arr('w_pred'), arr('w_fit')):.4f} "
                 f"(excl Dubliners: {cc(arr('w_pred')[notdub], arr('w_fit')[notdub]):.4f})")
    lines.append(f"- median ratio = {float(np.median(ratio)):.4f}; sqrt(e)={math.sqrt(math.e):.4f}")
    lines.append(f"- corr(log ratio, sd_h^2/2) = {cc(np.log(np.abs(ratio)), sdh**2 / 2):.4f}")
    lines.append(f"- corr(log ratio, sd_t^2/2) = {cc(np.log(np.abs(ratio)), sdt**2 / 2):.4f}")
    lines.append(f"- corr(log ratio, sd_t/sd_h) = {cc(np.log(np.abs(ratio)), sdt / sdh):.4f}")
    lines.append(f"\n## Hunch-1 redo (s vs corrected head size)")
    lines.append(f"- corr(log s, log N_H) = {cc(np.log(s_vals), np.log(nh)):.4f}")
    lines.append(f"- median s/N_H = {float(np.median(s_vals / nh)):.4f} "
                 f"(range {(s_vals/nh).min():.3f}-{(s_vals/nh).max():.3f})")
    lines.append(f"- log s ~ log N_H slope = "
                 f"{float(np.linalg.lstsq(np.column_stack([np.ones(len(nh)), np.log(nh)]), np.log(s_vals), rcond=None)[0][1]):.3f}")
    (OUT / "f4d_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
