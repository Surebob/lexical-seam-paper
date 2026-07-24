"""f19 — kill the parameter: does lambda-ZM still beat ZM with lambda FROZEN?

f3's free-fit lambda (exp generator, 25/25 winner) is nearly universal across
the English 25: median 20.6, range [14.2, 24.4]. If lambda can be frozen at a
universal value, the seam correction costs ZERO extra per-corpus parameters:
  ZM:            log f = a - b log(r+c)                  (3 params)
  lambda*-ZM:    log f = a - b log(r+c) + L* g(x)        (3 params, L* frozen)
with g(x) = exp(x-1) - x, x = 0.05 + 0.95 log(r)/log(V).

Protocols (all anti-circularity):
  1. LOO train fit: for corpus i, L*_i = median of free-fit lambdas of the
     OTHER 24 corpora; fit (a,b,c) with the term frozen; compare full-data
     RMSE at equal parameter count (same p => lower RMSE wins outright).
  2. Held-out (f3b protocol): binomial-thin 50/50 (seed fixed), fit both
     3-param models on half A (lambda still LOO-frozen), score fixed curves
     on half B, both folds.
  3. Cross-language transfer: L* = median of ALL 25 English free lambdas,
     applied to the 7 canonical non-English corpora (150k-token slices).
Free-lambda (4p) fits are also run as the ceiling.

Outputs: ../outputs/f19_results.csv, f19_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

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
SEED = 20260732


def load_english_counts(fname):
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


def load_multilang_counts(path, mode):
    text = path.read_text(encoding="utf-8", errors="ignore")
    if mode == "jieba":
        import jieba
        toks = [t for t in jieba.cut(text) if UNI_RE.fullmatch(t)]
    else:
        toks = UNI_RE.findall(text.lower())
    if len(toks) > SLICE_TOKENS:
        toks = toks[:SLICE_TOKENS]
    return Counter(toks)


def curve(counts):
    return np.array(sorted(counts.values(), reverse=True), dtype=np.float64)


def gx_of(ranks, V_norm):
    return np.exp(0.05 + 0.95 * np.log(ranks) / math.log(V_norm) - 1.0) - (
        0.05 + 0.95 * np.log(ranks) / math.log(V_norm))


def c_grid_for(max_rank):
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


def fit_lzm_free(ranks, logf, gx):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c), gx])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c), float(coef[2]))
    return {"a": best[1], "b": best[2], "c": best[3], "lam": best[4], "rmse": math.sqrt(best[0])}


def fit_lzm_frozen(ranks, logf, gx, lam_star):
    """3-param fit: lambda frozen, term moved to the target side."""
    y = logf - lam_star * gx
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(y), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, y, rcond=None)
        mse = float(np.mean((d @ coef - y) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c))
    return {"a": best[1], "b": best[2], "c": best[3], "lam": lam_star, "rmse": math.sqrt(best[0])}


def thin_split(counts, rng):
    words = list(counts.keys())
    n = np.array([counts[w] for w in words], dtype=np.int64)
    a = rng.binomial(n, 0.5)
    b = n - a
    fa = np.sort(a[a > 0])[::-1].astype(np.float64)
    fb = np.sort(b[b > 0])[::-1].astype(np.float64)
    return fa, fb


def heldout_rmse(fitres, has_lam, V_fit, test):
    m = min(V_fit, len(test))
    r = np.arange(1, m + 1, dtype=np.float64)
    pred = fitres["a"] - fitres["b"] * np.log(r + fitres["c"])
    if has_lam:
        pred = pred + fitres["lam"] * gx_of(r, V_fit)
    return float(np.sqrt(np.mean((pred - np.log(test[:m])) ** 2)))


def main():
    rng = np.random.default_rng(SEED)
    print("pass 1: free-lambda fits (English 25)...", flush=True)
    eng = []
    for name, fname in CORPORA:
        counts = load_english_counts(fname)
        freqs = curve(counts)
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=np.float64)
        logf = np.log(freqs)
        gx = gx_of(ranks, V)
        free = fit_lzm_free(ranks, logf, gx)
        eng.append({"name": name, "counts": counts, "freqs": freqs,
                    "lam_free": free["lam"], "rmse_free": free["rmse"]})
        print(f"  {name[:30]:31} lam_free={free['lam']:.2f}", flush=True)

    lams = np.array([e["lam_free"] for e in eng])
    rows = []
    print("pass 2: LOO frozen-lambda vs ZM (train + held-out)...", flush=True)
    for i, e in enumerate(eng):
        loo = float(np.median(np.delete(lams, i)))
        freqs = e["freqs"]
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=np.float64)
        logf = np.log(freqs)
        gx = gx_of(ranks, V)
        zm = fit_zm(ranks, logf)
        froz = fit_lzm_frozen(ranks, logf, gx, loo)
        fa, fb = thin_split(e["counts"], rng)
        ho = {}
        for tag, (tr, te) in [("AB", (fa, fb)), ("BA", (fb, fa))]:
            Vtr = len(tr)
            rtr = np.arange(1, Vtr + 1, dtype=np.float64)
            ltr = np.log(tr)
            gtr = gx_of(rtr, Vtr)
            z = fit_zm(rtr, ltr)
            f = fit_lzm_frozen(rtr, ltr, gtr, loo)
            ho[f"zm_{tag}"] = heldout_rmse(z, False, Vtr, te)
            ho[f"fz_{tag}"] = heldout_rmse(f, True, Vtr, te)
        row = {"corpus": e["name"], "panel": "en", "V": V, "lam_free": round(e["lam_free"], 2),
               "lam_frozen": round(loo, 2), "rmse_zm": round(zm["rmse"], 4),
               "rmse_frozen": round(froz["rmse"], 4), "rmse_free": round(e["rmse_free"], 4),
               "train_win": froz["rmse"] < zm["rmse"],
               "impr_pct": round(100 * (zm["rmse"] - froz["rmse"]) / zm["rmse"], 2),
               "retained_pct": round(100 * (zm["rmse"] - froz["rmse"]) / max(zm["rmse"] - e["rmse_free"], 1e-9), 1),
               "ho_zm_AB": round(ho["zm_AB"], 4), "ho_fz_AB": round(ho["fz_AB"], 4),
               "ho_zm_BA": round(ho["zm_BA"], 4), "ho_fz_BA": round(ho["fz_BA"], 4),
               "ho_wins": sum(1 for t in ("AB", "BA") if ho[f"fz_{t}"] < ho[f"zm_{t}"])}
        rows.append(row)
        print(f"  {e['name'][:28]:29} frozen={row['rmse_frozen']:.4f} zm={row['rmse_zm']:.4f} "
              f"impr={row['impr_pct']:+.1f}% retained={row['retained_pct']:.0f}% ho={row['ho_wins']}/2", flush=True)

    lam_global = float(np.median(lams))
    print(f"pass 3: cross-language transfer (lam* = {lam_global:.2f} from English)...", flush=True)
    for name, path, mode in MULTILANG:
        try:
            counts = load_multilang_counts(path, mode)
        except Exception as ex:
            print(f"  SKIP {name}: {ex}", flush=True)
            continue
        freqs = curve(counts)
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=np.float64)
        logf = np.log(freqs)
        gx = gx_of(ranks, V)
        zm = fit_zm(ranks, logf)
        froz = fit_lzm_frozen(ranks, logf, gx, lam_global)
        free = fit_lzm_free(ranks, logf, gx)
        fa, fb = thin_split(counts, rng)
        ho = {}
        for tag, (tr, te) in [("AB", (fa, fb)), ("BA", (fb, fa))]:
            Vtr = len(tr)
            rtr = np.arange(1, Vtr + 1, dtype=np.float64)
            ltr = np.log(tr)
            z = fit_zm(rtr, ltr)
            f = fit_lzm_frozen(rtr, ltr, gx_of(rtr, Vtr), lam_global)
            ho[f"zm_{tag}"] = heldout_rmse(z, False, Vtr, te)
            ho[f"fz_{tag}"] = heldout_rmse(f, True, Vtr, te)
        row = {"corpus": name, "panel": "ml", "V": V, "lam_free": round(free["lam"], 2),
               "lam_frozen": round(lam_global, 2), "rmse_zm": round(zm["rmse"], 4),
               "rmse_frozen": round(froz["rmse"], 4), "rmse_free": round(free["rmse"], 4),
               "train_win": froz["rmse"] < zm["rmse"],
               "impr_pct": round(100 * (zm["rmse"] - froz["rmse"]) / zm["rmse"], 2),
               "retained_pct": round(100 * (zm["rmse"] - froz["rmse"]) / max(zm["rmse"] - free["rmse"], 1e-9), 1),
               "ho_zm_AB": round(ho["zm_AB"], 4), "ho_fz_AB": round(ho["fz_AB"], 4),
               "ho_zm_BA": round(ho["zm_BA"], 4), "ho_fz_BA": round(ho["fz_BA"], 4),
               "ho_wins": sum(1 for t in ("AB", "BA") if ho[f"fz_{t}"] < ho[f"zm_{t}"])}
        rows.append(row)
        print(f"  {name[:28]:29} frozen={row['rmse_frozen']:.4f} zm={row['rmse_zm']:.4f} "
              f"impr={row['impr_pct']:+.1f}% (their lam_free={row['lam_free']:.1f}) ho={row['ho_wins']}/2", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f19_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    en = [r for r in rows if r["panel"] == "en"]
    ml = [r for r in rows if r["panel"] == "ml"]
    lines = ["# f19 — frozen-lambda (3-parameter) seam formula vs plain ZM\n",
             f"g fixed to exp(x-1)-x (25/25 winner in f3); lambda frozen by LOO median "
             f"(English) / global English median {lam_global:.2f} (transfer). Equal parameter "
             f"count everywhere: 3 vs 3.\n"]
    for label, sub in [("English 25 (LOO-frozen)", en), ("Non-English 7 (English-frozen transfer)", ml)]:
        if not sub:
            continue
        tw = sum(1 for r in sub if r["train_win"])
        how = sum(r["ho_wins"] for r in sub)
        med_i = float(np.median([r["impr_pct"] for r in sub]))
        med_r = float(np.median([r["retained_pct"] for r in sub]))
        lines.append(f"## {label}")
        lines.append(f"- full-data RMSE wins vs ZM (equal params): {tw}/{len(sub)} (median improvement {med_i:.1f}%)")
        lines.append(f"- held-out fold wins: {how}/{2*len(sub)}")
        lines.append(f"- median share of the free-lambda improvement retained: {med_r:.0f}%")
        worst = min(sub, key=lambda r: r["impr_pct"])
        lines.append(f"- worst case: {worst['corpus']} ({worst['impr_pct']:+.1f}%)\n")
    (OUT / "f19_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
