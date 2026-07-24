"""f22b — is the LIKELIHOOD-space amplitude also universal?

f22 found: free lambda-ZM wins in MOE's arena (69/75 vs ZM, 60/75 vs MOE),
but the LSQ constant 20.6 transplants badly (6/75) — the MLE objective is
token-mass-weighted, so its optimal amplitude is a different number.
Objective-relative constants: measure the MLE lambda per corpus, test its
universality, freeze at the LOO median, refight at 2 params.

Protocol identical to f22 (80/20 binomial x3, train support, truncated
normalization); pass 1 records lambda_MLE per corpus (fit on FULL counts,
one fit per corpus, for the freeze pool); pass 2 fights frozen-at-LOO
lambda** vs ZM and MOE on the same splits as f22 (same seed stream).

Outputs: ../outputs/f22b_results.csv, f22b_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_f22 import (CORPORA, N_SPLITS, SEED, load_counts, make_gx,  # noqa: E402
                     logsumexp, nll_per_token)
from scipy.optimize import minimize  # noqa: E402
from scipy.special import zeta as hurwitz  # noqa: E402

OUT = HERE.parent / "outputs"


def fit_lzm_pmf(train, lam_fixed=None):
    V = len(train)
    ranks, gx = make_gx(V)
    Ttr = float(train.sum())

    def log_p(theta):
        b, c = math.exp(theta[0]), math.exp(theta[1]) - 1e-9
        lam = lam_fixed if lam_fixed is not None else theta[2]
        u = -b * np.log(ranks + c) + lam * gx
        if not np.all(np.isfinite(u)):
            return None
        return u - logsumexp(u)

    def obj(theta):
        lp = log_p(theta)
        if lp is None:
            return 1e12
        return -float(np.dot(train, lp)) / Ttr

    if lam_fixed is None:
        starts = [[math.log(1.1), math.log(3.0), 8.0], [math.log(0.9), math.log(30.0), 3.0],
                  [math.log(1.3), math.log(1.0), 15.0], [math.log(1.0), math.log(10.0), 0.0],
                  [math.log(1.1), math.log(3.0), 25.0]]
    else:
        starts = [[math.log(1.1), math.log(3.0)], [math.log(0.9), math.log(30.0)],
                  [math.log(1.3), math.log(1.0)], [math.log(1.0), math.log(100.0)]]
    best = None
    for x0 in starts:
        try:
            res = minimize(obj, x0=x0, method="Nelder-Mead",
                           options={"maxiter": 1500, "xatol": 1e-6, "fatol": 1e-9})
        except Exception:
            continue
        if best is None or res.fun < best.fun:
            best = res
    lam = lam_fixed if lam_fixed is not None else float(best.x[2])
    return best, lam, log_p(best.x)


def fit_simple_pmf(train, kind):
    V = len(train)
    ranks, _ = make_gx(V)
    Ttr = float(train.sum())

    def log_p(theta):
        if kind == "zm":
            b, c = math.exp(theta[0]), math.exp(theta[1]) - 1e-9
            u = -b * np.log(ranks + c)
        else:
            alpha, beta = 1.0 + math.exp(theta[0]), math.exp(theta[1])
            bb = 1.0 - beta
            z = float(hurwitz(alpha, 1.0))
            num = np.power(ranks, -alpha) * beta * z
            den = (z - bb * hurwitz(alpha, ranks)) * (z - bb * hurwitz(alpha, ranks + 1.0))
            with np.errstate(divide="ignore", invalid="ignore"):
                u = np.log(num) - np.log(den)
        if not np.all(np.isfinite(u)):
            return None
        return u - logsumexp(u)

    def obj(theta):
        lp = log_p(theta)
        return 1e12 if lp is None else -float(np.dot(train, lp)) / Ttr

    starts = ([[math.log(0.1), math.log(1.0)], [math.log(0.6), math.log(0.3)],
               [math.log(0.6), math.log(3.0)], [math.log(0.2), math.log(8.0)]]
              if kind == "moe" else
              [[math.log(1.1), math.log(3.0)], [math.log(0.9), math.log(30.0)],
               [math.log(1.3), math.log(1.0)], [math.log(1.0), math.log(100.0)]])
    best = None
    for x0 in starts:
        try:
            res = minimize(obj, x0=x0, method="Nelder-Mead",
                           options={"maxiter": 1200, "xatol": 1e-6, "fatol": 1e-9})
        except Exception:
            continue
        if best is None or res.fun < best.fun:
            best = res
    return log_p(best.x)


def main():
    print("pass 1: MLE lambda per corpus (full counts)...", flush=True)
    lams = {}
    full = {}
    for name, fname in CORPORA:
        counts = load_counts(fname)
        n = np.array(sorted(counts.values(), reverse=True), dtype=np.float64)
        full[name] = counts
        _, lam, _ = fit_lzm_pmf(n)
        lams[name] = lam
        print(f"  {name[:30]:31} lam_mle={lam:7.2f}", flush=True)
    vals = np.array(list(lams.values()))
    print(f"lam_MLE: median {np.median(vals):.2f}, IQR [{np.percentile(vals,25):.2f}, "
          f"{np.percentile(vals,75):.2f}], range [{vals.min():.2f}, {vals.max():.2f}]", flush=True)

    print("pass 2: frozen-at-LOO fights (same splits as f22)...", flush=True)
    rng = np.random.default_rng(SEED)
    rows = []
    names = [n for n, _ in CORPORA]
    for name, fname in CORPORA:
        counts = full[name]
        words = list(counts.keys())
        n = np.array([counts[wd] for wd in words], dtype=np.int64)
        loo = float(np.median([lams[m] for m in names if m != name]))
        for split in range(N_SPLITS):
            tr = rng.binomial(n, 0.8)
            te = n - tr
            keep = tr > 0
            order = np.argsort(-tr[keep], kind="stable")
            train = tr[keep][order].astype(np.float64)
            test = te[keep][order].astype(np.float64)
            _, _, lp_froz = fit_lzm_pmf(train, lam_fixed=loo)
            lp_zm = fit_simple_pmf(train, "zm")
            lp_moe = fit_simple_pmf(train, "moe")
            row = {"corpus": name, "split": split, "lam_mle": round(lams[name], 2),
                   "lam_loo": round(loo, 2),
                   "nll_frozen2": round(nll_per_token(lp_froz, test), 5),
                   "nll_zm": round(nll_per_token(lp_zm, test), 5),
                   "nll_moe": round(nll_per_token(lp_moe, test), 5)}
            rows.append(row)
            print(f"  {name[:26]:27} s{split} frozen2={row['nll_frozen2']:.4f} "
                  f"zm={row['nll_zm']:.4f} moe={row['nll_moe']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f22b_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    nt = len(rows)
    w_zm = sum(1 for r in rows if r["nll_frozen2"] < r["nll_zm"])
    w_moe = sum(1 for r in rows if r["nll_frozen2"] < r["nll_moe"])
    lines = ["# f22b — likelihood-space universal amplitude\n",
             f"- lam_MLE across 25 corpora: median {np.median(vals):.2f}, IQR "
             f"[{np.percentile(vals,25):.2f}, {np.percentile(vals,75):.2f}], "
             f"range [{vals.min():.2f}, {vals.max():.2f}]",
             f"- frozen-at-LOO (2 params) beats ZM (2 params): {w_zm}/{nt} "
             f"(median delta {float(np.median([r['nll_zm']-r['nll_frozen2'] for r in rows])):+.5f})",
             f"- frozen-at-LOO beats MOEZipf (2 params): {w_moe}/{nt} "
             f"(median delta {float(np.median([r['nll_moe']-r['nll_frozen2'] for r in rows])):+.5f})"]
    (OUT / "f22b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
