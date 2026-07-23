"""F3b — held-out check for the lambda-ZM BIC claim (manuscript section 3.7).

The 42/42 record is in-sample BIC. Referee-proof version: binomial-thin each
corpus into two same-depth halves (the paper's own section-3.4 equivalence:
per-type Binomial(n, 1/2) split == random token split), fit ZM / lambda-ZM /
MOEZipf-lsq on half A's rank curve with the exact f3 protocol, then evaluate
each fitted model AS A FIXED FUNCTION OF RANK on half B's rank curve
(ranks 1..min(V_A, V_B), RMSE on log frequency). Both fold directions.

Model selection inside lambda-ZM (exp vs is generator) is made on TRAIN rmse
only; the selected variant is then scored on the held-out half.

Corpora: English 25 (the headline "25/25" claim) + canonical multilingual 7.
Outputs: ../outputs/f3b_heldout.csv, f3b_heldout_summary.md
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
ML1 = REPO / "data" / "zipf_multilang"
ML2 = REPO / "data" / "zipf_multilang_romance"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

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

MULTILANG = [
    ("Russian War and Peace", ML1 / "russian_war_and_peace" / "combined_clean.txt", "uni"),
    ("Mandarin Three Kingdoms", ML1 / "mandarin_three_kingdoms" / "combined_clean.txt", "jieba"),
    ("Arabic 1001 Nights", ML1 / "arabic_1001_nights" / "combined_clean.txt", "uni"),
    ("Latin Gallic Wars", ML1 / "latin_gallic_wars" / "combined_clean.txt", "uni"),
    ("French Les Miserables", ML2 / "french_les_miserables" / "combined_clean.txt", "uni"),
    ("Spanish Don Quixote", ML2 / "spanish_don_quixote" / "combined_clean.txt", "uni"),
    ("Dutch Max Havelaar", ML2 / "dutch_max_havelaar" / "combined_clean.txt", "uni"),
]

SLICE_TOKENS = 150_000
SEED = 20260723


def load_english_counts(fname: str) -> Counter:
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


def load_multilang_counts(path: Path, mode: str) -> Counter:
    text = path.read_text(encoding="utf-8", errors="ignore")
    if mode == "jieba":
        import jieba
        toks = [t for t in jieba.cut(text) if UNI_RE.fullmatch(t)]
    else:
        toks = UNI_RE.findall(text.lower())
    if len(toks) > SLICE_TOKENS:
        toks = toks[:SLICE_TOKENS]
    return Counter(toks)


def thin_split(counts: Counter, rng) -> tuple[np.ndarray, np.ndarray]:
    """Binomial(n, 1/2) per type == random per-token assignment. Returns the two
    half-corpus rank-frequency curves (descending, zero-count types dropped)."""
    words = list(counts.keys())
    n = np.array([counts[w] for w in words], dtype=np.int64)
    a = rng.binomial(n, 0.5)
    b = n - a
    fa = np.sort(a[a > 0])[::-1].astype(np.float64)
    fb = np.sort(b[b > 0])[::-1].astype(np.float64)
    return fa, fb


def c_grid_for(max_rank: float):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, max_rank, 2048)])


def fit_zm(ranks, logf):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c))
    return {"a": best[1], "b": best[2], "c": best[3], "rmse": math.sqrt(best[0])}


GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
}


def xcoord(ranks, V_norm):
    """f3's normalized log-rank coordinate, frozen at the FIT vocabulary size so
    the fitted model is a fixed function of rank at evaluation time."""
    logr = np.log(ranks)
    return 0.05 + 0.95 * logr / math.log(V_norm)


def fit_lambda_zm(ranks, logf, gx):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c), gx])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c), float(coef[2]))
    return {"a": best[1], "b": best[2], "c": best[3], "lam": best[4], "rmse": math.sqrt(best[0])}


def moe_log_pmf(alpha, beta, ranks):
    bb = 1.0 - beta
    z = float(hurwitz(alpha, 1.0))
    hz_k = hurwitz(alpha, ranks)
    hz_k1 = hurwitz(alpha, ranks + 1.0)
    num = np.power(ranks, -alpha) * beta * z
    den = (z - bb * hz_k) * (z - bb * hz_k1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(num) - np.log(den)


def fit_moe_lsq(ranks, freqs, logf):
    n_tok = float(freqs.sum())

    def unpack(t):
        return 1.0 + math.exp(t[0]), math.exp(t[1])

    def obj(t):
        alpha, beta = unpack(t)
        lp = moe_log_pmf(alpha, beta, ranks)
        if not np.all(np.isfinite(lp)):
            return 1e12
        return float(np.mean((math.log(n_tok) + lp - logf) ** 2))

    best = None
    for a0, b0 in [(0.1, 1.0), (0.6, 0.3), (0.6, 3.0), (0.2, 8.0), (1.0, 1.0)]:
        try:
            res = minimize(obj, x0=[math.log(a0), math.log(b0)], method="Nelder-Mead",
                           options={"maxiter": 800, "xatol": 1e-6, "fatol": 1e-10})
        except Exception:
            continue
        if best is None or res.fun < best.fun:
            best = res
    alpha, beta = unpack(best.x)
    return {"alpha": alpha, "beta": beta, "amp": math.log(n_tok)}


def fold(train: np.ndarray, test: np.ndarray):
    """Fit on train curve, evaluate fixed functions on test curve."""
    Vtr, Vte = len(train), len(test)
    r_tr = np.arange(1, Vtr + 1, dtype=np.float64)
    lf_tr = np.log(train)
    x_tr = xcoord(r_tr, Vtr)

    zm = fit_zm(r_tr, lf_tr)
    lzms = {k: fit_lambda_zm(r_tr, lf_tr, g(x_tr)) for k, g in GENS.items()}
    lkey = min(lzms, key=lambda k: lzms[k]["rmse"])  # selected on TRAIN only
    lzm = lzms[lkey]
    moe = fit_moe_lsq(r_tr, train, lf_tr)

    m = min(Vtr, Vte)
    r_ev = np.arange(1, m + 1, dtype=np.float64)
    lf_ev = np.log(test[:m])
    x_ev = xcoord(r_ev, Vtr)  # frozen at fit-time normalization

    pred_zm = zm["a"] - zm["b"] * np.log(r_ev + zm["c"])
    pred_lzm = lzm["a"] - lzm["b"] * np.log(r_ev + lzm["c"]) + lzm["lam"] * GENS[lkey](x_ev)
    pred_moe = moe["amp"] + moe_log_pmf(moe["alpha"], moe["beta"], r_ev)

    def rmse(p):
        return float(np.sqrt(np.mean((p - lf_ev) ** 2)))

    return {
        "V_train": Vtr, "V_test": Vte, "lzm_gen": lkey,
        "zm_train_rmse": zm["rmse"], "lzm_train_rmse": lzm["rmse"],
        "zm_ho": rmse(pred_zm), "lzm_ho": rmse(pred_lzm), "moe_ho": rmse(pred_moe),
    }


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    jobs = [(n, ("en", f)) for n, f in CORPORA] + [(n, ("ml", (p, m))) for n, p, m in MULTILANG]
    for name, (kind, spec) in jobs:
        try:
            counts = load_english_counts(spec) if kind == "en" else load_multilang_counts(*spec)
        except Exception as e:
            print(f"SKIP {name}: {e}", flush=True)
            continue
        fa, fb = thin_split(counts, rng)
        for tag, (tr, te) in [("A>B", (fa, fb)), ("B>A", (fb, fa))]:
            r = fold(tr, te)
            r.update({"corpus": name, "panel": kind, "fold": tag,
                      "lzm_impr_pct": 100 * (r["zm_ho"] - r["lzm_ho"]) / r["zm_ho"]})
            rows.append(r)
            print(f"{name[:28]:29} {tag}  zm_ho={r['zm_ho']:.4f}  lzm_ho={r['lzm_ho']:.4f}  "
                  f"moe_ho={r['moe_ho']:.4f}  impr={r['lzm_impr_pct']:+.1f}%", flush=True)

    cols = ["corpus", "panel", "fold", "V_train", "V_test", "lzm_gen",
            "zm_train_rmse", "lzm_train_rmse", "zm_ho", "lzm_ho", "moe_ho", "lzm_impr_pct"]
    with open(OUT / "f3b_heldout.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    en = [r for r in rows if r["panel"] == "en"]
    ml = [r for r in rows if r["panel"] == "ml"]
    lines = ["# F3b held-out summary\n",
             f"Protocol: Binomial(n,1/2) thin split (seed {SEED}), fit on one half "
             f"(f3 objective, generator selected on train), evaluate fixed model on the "
             f"other half's rank curve at ranks 1..min(V_A,V_B). Two folds per corpus.\n"]
    for label, sub, nf in [("English 25", en, len(en)), ("Multilingual 7", ml, len(ml))]:
        if not sub:
            continue
        w_zm = sum(1 for r in sub if r["lzm_ho"] < r["zm_ho"])
        w_moe = sum(1 for r in sub if r["lzm_ho"] < r["moe_ho"])
        med = float(np.median([r["lzm_impr_pct"] for r in sub]))
        lines.append(f"## {label} ({nf} fold-tests)")
        lines.append(f"- lambda-ZM beats ZM on held-out RMSE: {w_zm}/{nf} (median improvement {med:.1f}%)")
        lines.append(f"- lambda-ZM beats MOEZipf-lsq on held-out RMSE: {w_moe}/{nf}")
        worst = min(sub, key=lambda r: r["lzm_impr_pct"])
        lines.append(f"- worst case: {worst['corpus']} {worst['fold']} ({worst['lzm_impr_pct']:+.1f}%)\n")
    (OUT / "f3b_heldout_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)


if __name__ == "__main__":
    sys.exit(main())
