"""F3 — lambda-ZM vs MOEZipf (fair, both objectives) + head-weighted winner maps
+ multilingual lambda-ZM.

Parts:
  A. English 25: plain ZM (3p), lambda-ZM with g=exp (4p), MOEZipf (2p) fit two ways
     (MLE on token counts, and direct rank-curve LSQ). BIC comparison across families.
  B. Head-only (top-100) free-amplitude winner maps for {IS, exp, euclid, xpow}.
  C. Multilingual 7: tokenize combined_clean.txt (unicode regex; jieba for Mandarin;
     leading 150k-token slice where longer), ZM + lambda-ZM(is/exp).

Outputs: ../outputs/f3_english.csv, f3_multilang.csv, f3_summary.md
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


def freqs_from_counts(counts: Counter):
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return np.array([f for _, f in ranked], dtype=np.float64)


def load_english(fname: str):
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
    return freqs_from_counts(Counter(TOKEN_RE.findall(text[start:end].lower())))


def load_multilang(path: Path, mode: str):
    text = path.read_text(encoding="utf-8", errors="ignore")
    if mode == "jieba":
        import jieba
        toks = [t for t in jieba.cut(text) if UNI_RE.fullmatch(t)]
    else:
        toks = UNI_RE.findall(text.lower())
    if len(toks) > SLICE_TOKENS:
        toks = toks[:SLICE_TOKENS]
    return freqs_from_counts(Counter(toks))


def c_grid_for(max_rank: float):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, max_rank, 2048)])


def fit_zm(ranks, logf):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        pred = d @ coef
        mse = float(np.mean((pred - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c))
    return {"a": best[1], "b": best[2], "c": best[3], "rmse": math.sqrt(best[0])}


def fit_lambda_zm(ranks, logf, gx):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c), gx])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        pred = d @ coef
        mse = float(np.mean((pred - logf) ** 2))
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


def fit_moe(ranks, freqs, logf, mode: str):
    n_tok = float(freqs.sum())

    def unpack(t):
        return 1.0 + math.exp(t[0]), math.exp(t[1])

    def nll(t):
        alpha, beta = unpack(t)
        lp = moe_log_pmf(alpha, beta, ranks)
        if not np.all(np.isfinite(lp)):
            return 1e12
        return -float(np.dot(freqs, lp))

    def rank_rmse_obj(t):
        alpha, beta = unpack(t)
        lp = moe_log_pmf(alpha, beta, ranks)
        if not np.all(np.isfinite(lp)):
            return 1e12
        pred = math.log(n_tok) + lp
        return float(np.mean((pred - logf) ** 2))

    obj = nll if mode == "mle" else rank_rmse_obj
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
    lp = moe_log_pmf(alpha, beta, ranks)
    pred = math.log(n_tok) + lp
    return {"alpha": alpha, "beta": beta, "rmse": float(np.sqrt(np.mean((pred - logf) ** 2)))}


def bic(n, p, rmse):
    return p * math.log(n) + n * math.log(rmse**2)


GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


def analyze(name, freqs, is_english=True):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    logf = np.log(freqs)
    logr = np.log(ranks)
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])

    zm = fit_zm(ranks, logf)
    lzm_exp = fit_lambda_zm(ranks, logf, GENS["exp"](x))
    lzm_is = fit_lambda_zm(ranks, logf, GENS["is"](x))
    row = {
        "corpus": name, "V": V, "tokens": int(freqs.sum()),
        "zm_c": zm["c"], "zm_rmse": zm["rmse"],
        "lzm_exp_rmse": lzm_exp["rmse"], "lzm_exp_lam": lzm_exp["lam"], "lzm_exp_c": lzm_exp["c"],
        "lzm_is_rmse": lzm_is["rmse"], "lzm_is_lam": lzm_is["lam"],
        "lzm_impr_pct": 100 * (zm["rmse"] - min(lzm_exp["rmse"], lzm_is["rmse"])) / zm["rmse"],
    }
    if is_english:
        moe_mle = fit_moe(ranks, freqs, logf, "mle")
        moe_lsq = fit_moe(ranks, freqs, logf, "lsq")
        b_zm = bic(V, 3, zm["rmse"])
        b_lzm = bic(V, 4, min(lzm_exp["rmse"], lzm_is["rmse"]))
        b_moe = bic(V, 2, moe_lsq["rmse"])
        winner = min([("zm", b_zm), ("lzm", b_lzm), ("moe", b_moe)], key=lambda t: t[1])[0]
        # head-only (top-100) free-amplitude winner
        pred_zm = zm["a"] - zm["b"] * np.log(ranks + zm["c"])
        resid = logf - pred_zm
        h = slice(0, min(100, V))
        head_winner, head_rmse = None, None
        full_winner, full_best = None, None
        for kname, g in GENS.items():
            gx = g(x)
            lam_h = float(np.dot(resid[h], gx[h]) / np.dot(gx[h], gx[h]))
            rh = float(np.sqrt(np.mean((resid[h] - lam_h * gx[h]) ** 2)))
            if head_rmse is None or rh < head_rmse:
                head_rmse, head_winner = rh, kname
            lam_f = float(np.dot(resid, gx) / np.dot(gx, gx))
            rf = float(np.sqrt(np.mean((resid - lam_f * gx) ** 2)))
            if full_best is None or rf < full_best:
                full_best, full_winner = rf, kname
        row.update({
            "moe_mle_rmse": moe_mle["rmse"], "moe_mle_alpha": moe_mle["alpha"], "moe_mle_beta": moe_mle["beta"],
            "moe_lsq_rmse": moe_lsq["rmse"], "moe_lsq_alpha": moe_lsq["alpha"], "moe_lsq_beta": moe_lsq["beta"],
            "bic_zm": b_zm, "bic_lzm": b_lzm, "bic_moe": b_moe, "bic_winner": winner,
            "lzm_minus_moe_rmse": min(lzm_exp["rmse"], lzm_is["rmse"]) - moe_lsq["rmse"],
            "head100_free_winner": head_winner, "full_free_winner": full_winner,
        })
    return row


def main():
    eng = []
    for name, fname in CORPORA:
        freqs = load_english(fname)
        eng.append(analyze(name, freqs, is_english=True))
        r = eng[-1]
        print(f"{name[:30]:31} zm={r['zm_rmse']:.4f} lzm={min(r['lzm_exp_rmse'], r['lzm_is_rmse']):.4f} "
              f"moe_lsq={r['moe_lsq_rmse']:.4f} bic_winner={r['bic_winner']} head100={r['head100_free_winner']}", flush=True)

    ml = []
    for name, path, mode in MULTILANG:
        try:
            freqs = load_multilang(path, mode)
        except Exception as e:
            print(f"SKIP {name}: {e}", flush=True)
            continue
        ml.append(analyze(name, freqs, is_english=False))
        r = ml[-1]
        print(f"{name[:30]:31} V={r['V']} c={r['zm_c']:.2f} zm={r['zm_rmse']:.4f} "
              f"lzm_impr={r['lzm_impr_pct']:.2f}%", flush=True)

    for rows, fn in [(eng, "f3_english.csv"), (ml, "f3_multilang.csv")]:
        cols = sorted({k for r in rows for k in r}, key=lambda k: (k != "corpus", k))
        with open(OUT / fn, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)

    wc = Counter(r["bic_winner"] for r in eng)
    hw = Counter(r["head100_free_winner"] for r in eng)
    fw = Counter(r["full_free_winner"] for r in eng)
    lines = ["# F3 summary\n"]
    lines.append(f"## English 25 — 3-way BIC (zm p=3, lambda-ZM p=4, MOE-lsq p=2)")
    lines.append(f"- BIC winners: {dict(wc)}")
    lines.append(f"- lambda-ZM rank-RMSE beats MOE(lsq): {sum(1 for r in eng if r['lzm_minus_moe_rmse'] < 0)}/25 "
                 f"(median delta {float(np.median([r['lzm_minus_moe_rmse'] for r in eng])):.5f})")
    lines.append(f"- lambda-ZM beats MOE(mle-fit): {sum(1 for r in eng if min(r['lzm_exp_rmse'], r['lzm_is_rmse']) < r['moe_mle_rmse'])}/25")
    lines.append(f"- median MOE(lsq) alpha={float(np.median([r['moe_lsq_alpha'] for r in eng])):.3f} "
                 f"beta={float(np.median([r['moe_lsq_beta'] for r in eng])):.3f}")
    lines.append(f"\n## Head-only (top-100) free-amplitude winner map: {dict(hw)}")
    lines.append(f"## Full free-amplitude winner map (check vs f1): {dict(fw)}")
    lines.append(f"\n## Multilingual lambda-ZM ({len(ml)}/7)")
    for r in ml:
        lines.append(f"- {r['corpus']}: V={r['V']} c={r['zm_c']:.2f} zm_rmse={r['zm_rmse']:.4f} "
                     f"impr={r['lzm_impr_pct']:.2f}% lam_exp={r['lzm_exp_lam']:.2f}")
    (OUT / "f3_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nwrote outputs", flush=True)


if __name__ == "__main__":
    sys.exit(main())
