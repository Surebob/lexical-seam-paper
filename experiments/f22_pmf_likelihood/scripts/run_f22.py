"""f22 — likelihood-space bake-off: lambda-ZM as a PMF, on MOEZipf's home turf.

The paper's comparisons are rank-curve LSQ. MOEZipf's native arena is
maximum likelihood on token counts. Here every model becomes a proper PMF
over the train rank support and fights by held-out negative log-likelihood:

  ZM PMF          p(r) prop (r+c)^-b                    2 params (b, c)
  frozen lzm PMF  p(r) prop exp(-b log(r+c) + 20.6 g)   2 params (b, c)
  free lzm PMF    + lam free                             3 params
  MOEZipf PMF     truncated-renormalized over [1, V]     2 params (alpha, beta)

Protocol (the old canon's): per corpus, 80/20 binomial split of each type's
count (3 splits, fixed seeds); rank map and support from the TRAIN counts;
test counts projected onto train types; all models normalized over
[1, V_train]; score = held-out NLL per test token. Identical treatment for
every model; amplitude lives in the normalizer, so ZM/frozen/MOE all carry
exactly TWO free parameters — the equal-complexity fight, in MOE's arena.

Outputs: ../outputs/f22_results.csv, f22_summary.md
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
from scipy.special import zeta as hurwitz

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
LAM_STAR = 20.6
SEED = 20260735
N_SPLITS = 3

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


def load_counts(fname):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    start, end = 0, len(text)
    for m in START_MARKERS:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            start = nl + 1 if nl != -1 else i
            break
    for m in END_MARKERS:
        i = text.find(m)
        if i != -1:
            end = i
            break
    return Counter(TOKEN_RE.findall(text[start:end].lower()))


def logsumexp(v):
    m = np.max(v)
    return m + math.log(float(np.sum(np.exp(v - m))))


def nll_per_token(log_p, test_counts):
    T = float(test_counts.sum())
    return -float(np.dot(test_counts, log_p)) / T


def make_gx(V):
    ranks = np.arange(1, V + 1, dtype=np.float64)
    x = 0.05 + 0.95 * np.log(ranks) / math.log(V)
    return ranks, np.exp(x - 1.0) - x


def fit_pmf(train, test, kind):
    """Fit PMF by MLE on train counts (ranks = train order), score test NLL/token."""
    V = len(train)
    ranks, gx = make_gx(V)
    lr_base = np.log(ranks)

    def log_p(theta):
        if kind == "zm":
            b, c = math.exp(theta[0]), math.exp(theta[1]) - 1e-9
            u = -b * np.log(ranks + c)
        elif kind == "lzm_frozen":
            b, c = math.exp(theta[0]), math.exp(theta[1]) - 1e-9
            u = -b * np.log(ranks + c) + LAM_STAR * gx
        elif kind == "lzm_free":
            b, c, lam = math.exp(theta[0]), math.exp(theta[1]) - 1e-9, theta[2]
            u = -b * np.log(ranks + c) + lam * gx
        elif kind == "moe":
            alpha, beta = 1.0 + math.exp(theta[0]), math.exp(theta[1])
            bb = 1.0 - beta
            z = float(hurwitz(alpha, 1.0))
            hz_k = hurwitz(alpha, ranks)
            hz_k1 = hurwitz(alpha, ranks + 1.0)
            num = np.power(ranks, -alpha) * beta * z
            den = (z - bb * hz_k) * (z - bb * hz_k1)
            with np.errstate(divide="ignore", invalid="ignore"):
                u = np.log(num) - np.log(den)
        if not np.all(np.isfinite(u)):
            return None
        return u - logsumexp(u)

    Ttr = float(train.sum())

    def obj(theta):
        lp = log_p(theta)
        if lp is None:
            return 1e12
        return -float(np.dot(train, lp)) / Ttr

    if kind == "moe":
        starts = [[math.log(0.1), math.log(1.0)], [math.log(0.6), math.log(0.3)],
                  [math.log(0.6), math.log(3.0)], [math.log(0.2), math.log(8.0)],
                  [math.log(1.0), math.log(1.0)]]
    elif kind == "lzm_free":
        starts = [[math.log(1.1), math.log(3.0), 20.0], [math.log(0.9), math.log(30.0), 15.0],
                  [math.log(1.3), math.log(1.0), 25.0], [math.log(1.0), math.log(10.0), 10.0]]
    else:
        starts = [[math.log(1.1), math.log(3.0)], [math.log(0.9), math.log(30.0)],
                  [math.log(1.3), math.log(1.0)], [math.log(1.0), math.log(100.0)]]
    best = None
    for x0 in starts:
        try:
            res = minimize(obj, x0=x0, method="Nelder-Mead",
                           options={"maxiter": 1200, "xatol": 1e-6, "fatol": 1e-9})
        except Exception:
            continue
        if best is None or res.fun < best.fun:
            best = res
    lp = log_p(best.x)
    return nll_per_token(lp, test)


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for name, fname in CORPORA:
        counts = load_counts(fname)
        words = list(counts.keys())
        n = np.array([counts[wd] for wd in words], dtype=np.int64)
        for split in range(N_SPLITS):
            tr = rng.binomial(n, 0.8)
            te = n - tr
            keep = tr > 0
            tr_k, te_k = tr[keep], te[keep]
            order = np.argsort(-tr_k, kind="stable")
            train = tr_k[order].astype(np.float64)
            test = te_k[order].astype(np.float64)
            res = {k: fit_pmf(train, test, k) for k in ["zm", "lzm_frozen", "lzm_free", "moe"]}
            unseen_frac = float(te[~keep].sum()) / max(float(te.sum() + te[~keep].sum()), 1.0)
            row = {"corpus": name, "split": split, "V_train": len(train),
                   "unseen_mass": round(unseen_frac, 4)}
            row.update({f"nll_{k}": round(v, 5) for k, v in res.items()})
            rows.append(row)
            print(f"{name[:26]:27} s{split} zm={res['zm']:.4f} frozen={res['lzm_frozen']:.4f} "
                  f"free={res['lzm_free']:.4f} moe={res['moe']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f22_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    lines = ["# f22 — likelihood-space PMF bake-off (held-out NLL per token)\n",
             f"80/20 binomial splits x{N_SPLITS}; support/ranks from train; all models "
             "normalized over [1, V_train]. ZM, frozen lambda-ZM, MOEZipf carry 2 free "
             "params each; free lambda-ZM carries 3.\n"]
    nt = len(rows)
    for a, b in [("lzm_frozen", "zm"), ("lzm_frozen", "moe"), ("lzm_free", "zm"),
                 ("lzm_free", "moe"), ("lzm_free", "lzm_frozen")]:
        wins = sum(1 for r in rows if r[f"nll_{a}"] < r[f"nll_{b}"])
        med = float(np.median([r[f"nll_{b}"] - r[f"nll_{a}"] for r in rows]))
        lines.append(f"- {a} beats {b}: {wins}/{nt} (median NLL delta {med:+.5f} nats/token)")
    (OUT / "f22_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
