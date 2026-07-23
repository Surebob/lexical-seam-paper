"""F6 — c as a sampling-depth dial: slice trajectories, winner flips, and the
PLN-thinning prediction of c(T).

Parts:
  A. For 6 English corpora x slices {full, 300k, 150k, 80k, 40k} tokens:
     ZM fit (a,b,c), unit-amplitude step-2 winner, free-amplitude winner.
     -> does c collapse and does the IS/exp winner flip with size on the SAME text?
  B. PLN thinning prediction: from the f4d full-corpus mixture fit, predict the
     c(T) trajectory with NO new parameters — simulate counts ~ Poisson(lambda*T/T0)
     over the latent mixture (including unseen types via zero-truncation mass),
     fit ZM c on each simulated slice, compare to the empirical trajectory.

Outputs: ../outputs/f6_slices.csv, f6_pln_prediction.csv, f6_summary.md
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
from scipy.special import gammaln, logsumexp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
F4D = REPO / "experiments" / "f4d_poisson_lognormal_mixture" / "outputs" / "f4d_mixture_fits.csv"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"), ("War and Peace", "pg2600.txt"),
    ("King James Bible", "pg10.txt"), ("Les Miserables", "pg135.txt"),
    ("Don Quixote", "pg996.txt"), ("Moby Dick", "pg2701.txt"),
]
SLICES = [None, 300_000, 150_000, 80_000, 40_000]
SEED = 20260722
GH_X, GH_W = hermgauss(48)

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


def tokens_of(fname):
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
    return TOKEN_RE.findall(text[st:e].lower())


def zm_and_winners(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    logr = np.log(ranks)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 1024)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(c), d @ coef)
    _, c, pred = best
    resid = logf - pred
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
    unit = {k: float(np.sqrt(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
    free = {}
    for k, g in GENS.items():
        gx = g(x)
        lam = float(np.dot(resid, gx) / np.dot(gx, gx))
        free[k] = float(np.sqrt(np.mean((resid - lam * gx) ** 2)))
    return c, min(unit, key=unit.get), min(free, key=free.get), math.sqrt(best[0]), V


def zm_c_only(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 768)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(c))
    return best[1]


def pln_p0(mu, sd):
    loglam = mu + math.sqrt(2.0) * sd * GH_X
    lam = np.exp(np.clip(loglam, -30, 30))
    return float(np.exp(logsumexp(-lam + np.log(GH_W)) - 0.5 * math.log(math.pi)))


def simulate_pln_at_scale(fit, V_obs, scale, rng):
    """Simulate counts from the fitted mixture with all rates scaled by `scale`.
    Total latent type count inferred from zero-truncation at scale=1."""
    pi_h, mu_h, sd_h = fit["pi_h"], fit["mu_h"], fit["sd_h"]
    mu_t, sd_t = fit["mu_t"], fit["sd_t"]
    p0 = pi_h * pln_p0(mu_h, sd_h) + (1 - pi_h) * pln_p0(mu_t, sd_t)
    n_total = int(round(V_obs / max(1e-9, 1.0 - p0)))
    comp = rng.random(n_total) < pi_h
    loglam = np.where(comp, rng.normal(mu_h, sd_h, n_total), rng.normal(mu_t, sd_t, n_total))
    lam = np.exp(np.clip(loglam, -30, 30)) * scale
    counts = rng.poisson(lam)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(np.float64)


def main():
    f4d = {r["corpus"]: {k: float(r[k]) for k in ("pi_h", "mu_h", "sd_h", "mu_t", "sd_t")}
           for r in csv.DictReader(open(F4D, newline="", encoding="utf-8"))}

    rows = []
    fulltoks = {}
    for name, fname in CORPORA:
        toks = tokens_of(fname)
        fulltoks[name] = len(toks)
        for s in SLICES:
            tt = toks if s is None else toks[:s]
            if s is not None and len(toks) < s:
                continue
            freqs = np.array(sorted(Counter(tt).values(), reverse=True), dtype=float)
            c, uw, fw, rmse, V = zm_and_winners(freqs)
            rows.append({"corpus": name, "tokens": len(tt), "slice": "full" if s is None else s,
                         "V": V, "zm_c": c, "zm_rmse": rmse, "unit_winner": uw, "free_winner": fw})
            print(f"{name[:26]:27} T={len(tt):7} V={V:6} c={c:7.1f} unit={uw:6} free={fw}", flush=True)

    with open(OUT / "f6_slices.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Part B: PLN thinning prediction of c(T) for 3 exemplars
    rng = np.random.default_rng(SEED)
    pred_rows = []
    for name in ["Complete Works of Shakespeare", "War and Peace", "Moby Dick"]:
        fit = f4d[name]
        V_obs = next(r["V"] for r in rows if r["corpus"] == name and r["slice"] == "full")
        T0 = fulltoks[name]
        for s in SLICES:
            T = T0 if s is None else s
            if T > T0:
                continue
            cs = []
            for rep in range(3):
                freqs = simulate_pln_at_scale(fit, V_obs, T / T0, rng)
                cs.append(zm_c_only(freqs))
            c_emp = next((r["zm_c"] for r in rows if r["corpus"] == name and
                          (r["slice"] == ("full" if s is None else s))), float("nan"))
            pred_rows.append({"corpus": name, "T": T, "c_pred_median": float(np.median(cs)),
                              "c_pred_min": float(np.min(cs)), "c_pred_max": float(np.max(cs)),
                              "c_empirical": c_emp})
            print(f"PLN {name[:24]:25} T={T:7} c_pred={np.median(cs):7.1f} "
                  f"[{np.min(cs):.1f},{np.max(cs):.1f}] c_emp={c_emp:7.1f}", flush=True)

    with open(OUT / "f6_pln_prediction.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(pred_rows[0].keys()))
        w.writeheader()
        w.writerows(pred_rows)

    lines = ["# F6 — c as sampling depth\n"]
    lines.append("## A. Empirical slice trajectories (c and winners vs T)")
    for name, _ in CORPORA:
        sub = [r for r in rows if r["corpus"] == name]
        lines.append(f"- {name}: " + " | ".join(
            f"T={r['tokens']//1000}k c={r['zm_c']:.1f} {r['unit_winner']}" for r in sub))
    flips = sum(1 for name, _ in CORPORA
                if len({r["unit_winner"] for r in rows if r["corpus"] == name}) > 1)
    lines.append(f"\n- corpora whose unit winner CHANGES with slice size: {flips}/{len(CORPORA)}")
    lines.append("\n## B. PLN thinning prediction of c(T) (no new parameters)")
    for r in pred_rows:
        lines.append(f"- {r['corpus']} T={r['T']}: pred {r['c_pred_median']:.1f} "
                     f"[{r['c_pred_min']:.1f},{r['c_pred_max']:.1f}] vs empirical {r['c_empirical']:.1f}")
    lp = np.log1p(np.array([r["c_pred_median"] for r in pred_rows]))
    le = np.log1p(np.array([r["c_empirical"] for r in pred_rows]))
    lines.append(f"\n- corr(log1p c_pred, log1p c_emp) = {float(np.corrcoef(lp, le)[0, 1]):.4f}")
    (OUT / "f6_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
