"""F1 — Fresh reproduction + proper re-experiments (2026-07-20).

Independent reimplementation (no legacy-engine imports). Protocol spec matches canon:
  tokens  = [a-z]+(?:'[a-z]+)?  on lowercased text, Gutenberg boilerplate stripped
  ranking = count desc, tie-break word ascending
  ZM fit  = grid over c in {0} + geomspace(1e-6, max_rank, 2048), affine LSQ in log space
  x       = 0.05 + 0.95 * minmax(log rank)
  step-2 composite = RMSE(residual - g(x)) at unit amplitude (canonical convention)

Fresh additions beyond reproduction:
  1. Continuous ZM refit (multi-start least_squares) auditing the grid optimum.
  2. Free-amplitude generator comparison (lambda fit per generator) — amplitude-fair
     winner map, vs the canonical forced-amplitude-1 convention.
  3. lambda-ZM 4-parameter joint fit (a, b, c, lambda·g) via c-grid + 3-col LSQ,
     OLS and frequency-weighted (WLS), for g in {IS, EXP}. BIC vs plain ZM.
  4. Two-population exponent gap (top-100 function proxy vs re-ranked content tail).
  5. Lexical stats: TTR, hapax ratio.

Outputs: ../outputs/f1_per_corpus.csv, ../outputs/f1_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
CANON_TABLE1 = REPO / "experiments" / "1a_per_corpus_enriched_search" / "outputs" / "table1_per_corpus.csv"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
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


def load_ranked_freqs(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
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
    counts = Counter(TOKEN_RE.findall(text[start:end].lower()))
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    freqs = np.array([f for _, f in ranked], dtype=np.float64)
    return freqs, int(freqs.sum())


def affine(z, y, w=None):
    d = np.column_stack([np.ones_like(z), z])
    if w is not None:
        sw = np.sqrt(w)
        coef, *_ = np.linalg.lstsq(d * sw[:, None], y * sw, rcond=None)
    else:
        coef, *_ = np.linalg.lstsq(d, y, rcond=None)
    pred = d @ coef
    return coef, pred


def c_grid_for(max_rank: float):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, max_rank, 2048)])


def fit_zm_grid(ranks, logf):
    best = None
    for c in c_grid_for(ranks[-1]):
        coef, pred = affine(np.log(ranks + c), logf)
        mse = float(np.mean((pred - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c), pred)
    mse, a, b, c, pred = best
    return {"a": a, "b": b, "c": c, "rmse": math.sqrt(mse), "pred": pred}


def fit_zm_continuous(ranks, logf, c0_list):
    def resid(p):
        a, b, logc = p
        return a - b * np.log(ranks + np.exp(logc)) - logf

    best = None
    for c0 in c0_list:
        try:
            sol = least_squares(resid, x0=[logf[0], 1.0, math.log(max(c0, 1e-6))], method="lm", max_nfev=4000)
        except Exception:
            continue
        rmse = math.sqrt(float(np.mean(sol.fun ** 2)))
        if best is None or rmse < best[0]:
            best = (rmse, float(sol.x[0]), float(sol.x[1]), float(math.exp(sol.x[2])))
    return {"rmse": best[0], "a": best[1], "b": best[2], "c": best[3]}


GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


def lambda_zm_fit(ranks, logf, g_of_x, w=None):
    """Joint 4-param fit: logf ~ a - b*log(r+c) + lam*g(x) via c-grid + 3-col LSQ."""
    gx = g_of_x
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c), gx])
        if w is not None:
            sw = np.sqrt(w)
            coef, *_ = np.linalg.lstsq(d * sw[:, None], logf * sw, rcond=None)
        else:
            coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        pred = d @ coef
        mse = float(np.mean((pred - logf) ** 2))  # always report unweighted mse too
        if w is not None:
            wmse = float(np.average((pred - logf) ** 2, weights=w))
            key = wmse
        else:
            wmse = mse
            key = mse
        if best is None or key < best[0]:
            best = (key, mse, wmse, float(coef[0]), float(-coef[1]), float(c), float(coef[2]))
    _, mse, wmse, a, b, c, lam = best
    return {"a": a, "b": b, "c": c, "lam": lam, "rmse": math.sqrt(mse), "wrmse": math.sqrt(wmse)}


def bic(n, p, mse):
    return p * math.log(n) + n * math.log(mse)


def main():
    canon = {}
    with open(CANON_TABLE1, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            canon[row["corpus"]] = row

    rows = []
    for name, fname in CORPORA:
        path = DATA / fname
        if not path.exists():
            print(f"MISSING corpus file: {name} -> {fname}", flush=True)
            continue
        freqs, tokens = load_ranked_freqs(path)
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=np.float64)
        logf = np.log(freqs)
        logr = np.log(ranks)

        zm = fit_zm_grid(ranks, logf)
        zmc = fit_zm_continuous(ranks, logf, c0_list=[0.01, 1.0, 10.0, 100.0, zm["c"] + 1e-6, 1000.0])
        resid = logf - zm["pred"]

        x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])

        unit_rmse = {k: math.sqrt(float(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
        unit_winner = min(unit_rmse, key=unit_rmse.get)

        free = {}
        for k, g in GENS.items():
            gx = g(x)
            lam = float(np.dot(resid, gx) / np.dot(gx, gx))
            free[k] = (math.sqrt(float(np.mean((resid - lam * gx) ** 2))), lam)
        free_winner = min(free, key=lambda k: free[k][0])

        lam_is = lambda_zm_fit(ranks, logf, GENS["is"](x))
        lam_exp = lambda_zm_fit(ranks, logf, GENS["exp"](x))
        lam_is_w = lambda_zm_fit(ranks, logf, GENS["is"](x), w=freqs)
        joint_winner = "is" if lam_is["rmse"] <= lam_exp["rmse"] else "exp"
        jw = lam_is if joint_winner == "is" else lam_exp
        bic_zm = bic(V, 3, zm["rmse"] ** 2)
        bic_lzm = bic(V, 4, jw["rmse"] ** 2)

        n_f = min(100, V)
        af = -np.polyfit(logr[:n_f], logf[:n_f], 1)[0]
        tail = freqs[n_f:]
        tr = np.arange(1, len(tail) + 1, dtype=np.float64)
        ac = -np.polyfit(np.log(tr), np.log(tail), 1)[0]

        hapax = float(np.sum(freqs == 1.0)) / V
        ttr = V / tokens

        cn = canon.get(name, {})
        c_canon = float(cn.get("zm_c", "nan"))
        fam_canon = cn.get("winner_family", "?")
        rows.append({
            "corpus": name, "tokens": tokens, "V": V,
            "zm_a": zm["a"], "zm_b": zm["b"], "zm_c": zm["c"], "zm_rmse": zm["rmse"],
            "canon_c": c_canon, "dc": zm["c"] - c_canon,
            "canon_tokens": int(cn.get("tokens", -1)), "canon_rmse": float(cn.get("zm_rmse", "nan")),
            "cont_rmse": zmc["rmse"], "cont_minus_grid": zmc["rmse"] - zm["rmse"], "cont_c": zmc["c"],
            "unit_winner": unit_winner, "canon_family": fam_canon,
            "unit_match": unit_winner == fam_canon,
            "unit_is_rmse": unit_rmse["is"], "unit_exp_rmse": unit_rmse["exp"],
            "unit_euclid_rmse": unit_rmse["euclid"], "unit_xpow_rmse": unit_rmse["xpow"],
            "free_winner": free_winner,
            "free_is_rmse": free[
                "is"][0], "free_is_lam": free["is"][1],
            "free_exp_rmse": free["exp"][0], "free_exp_lam": free["exp"][1],
            "free_euclid_rmse": free["euclid"][0], "free_xpow_rmse": free["xpow"][0],
            "lzm_winner": joint_winner,
            "lzm_is_rmse": lam_is["rmse"], "lzm_is_lam": lam_is["lam"], "lzm_is_c": lam_is["c"],
            "lzm_exp_rmse": lam_exp["rmse"], "lzm_exp_lam": lam_exp["lam"], "lzm_exp_c": lam_exp["c"],
            "lzm_impr_pct": 100.0 * (zm["rmse"] - jw["rmse"]) / zm["rmse"],
            "lzm_wls_lam": lam_is_w["lam"], "lzm_wls_wrmse": lam_is_w["wrmse"],
            "bic_zm": bic_zm, "bic_lzm": bic_lzm, "lzm_bic_wins": bic_lzm < bic_zm,
            "alpha_function": af, "alpha_content": ac, "exp_gap": af - ac,
            "hapax_ratio": hapax, "ttr": ttr,
        })
        print(f"done {name}: c={zm['c']:.2f} (canon {c_canon:.2f}) unit={unit_winner}/{fam_canon} "
              f"free={free_winner} lzm={joint_winner} lam={jw['lam']:.3f} impr={rows[-1]['lzm_impr_pct']:.2f}%", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f1_per_corpus.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=cols)
        wcsv.writeheader()
        wcsv.writerows(rows)

    def arr(key):
        return np.array([r[key] for r in rows], dtype=np.float64)

    def corr(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    is_flag = np.array([1.0 if r["canon_family"] == "is" else 0.0 for r in rows])
    c_all = arr("zm_c")
    lines = []
    lines.append("# F1 fresh reproduction + proper re-experiments — summary\n")
    lines.append(f"- corpora run: {len(rows)}/25")
    lines.append(f"- token-count exact matches vs canon: {sum(1 for r in rows if r['tokens'] == r['canon_tokens'])}/{len(rows)}")
    lines.append(f"- max |c - canon_c|: {float(np.max(np.abs(arr('dc')))):.6f}")
    lines.append(f"- max |zm_rmse - canon_rmse|: {float(np.max(np.abs(arr('zm_rmse') - arr('canon_rmse')))):.8f}")
    lines.append(f"- unit-amplitude winner-family matches: {sum(1 for r in rows if r['unit_match'])}/{len(rows)}")
    lines.append(f"- continuous refit beats grid by >1e-4: {sum(1 for r in rows if r['cont_minus_grid'] < -1e-4)}/{len(rows)}"
                 f" (min delta {float(np.min(arr('cont_minus_grid'))):.6f})\n")
    fw = [r["free_winner"] for r in rows]
    lines.append("## Amplitude-fair (free-lambda) winner map")
    for k in GENS:
        lines.append(f"- {k}: {fw.count(k)}/{len(rows)}")
    lines.append(f"- corpora where free winner != canonical family: "
                 f"{sum(1 for r in rows if r['free_winner'] != r['canon_family'])}\n")
    lines.append("## lambda-ZM 4-parameter model (joint OLS fit)")
    lines.append(f"- winner (is vs exp): is {sum(1 for r in rows if r['lzm_winner']=='is')}, "
                 f"exp {sum(1 for r in rows if r['lzm_winner']=='exp')}")
    lines.append(f"- RMSE improvement over ZM: min {float(np.min(arr('lzm_impr_pct'))):.3f}% "
                 f"median {float(np.median(arr('lzm_impr_pct'))):.3f}% max {float(np.max(arr('lzm_impr_pct'))):.3f}%")
    lines.append(f"- BIC(4-param) beats BIC(ZM): {sum(1 for r in rows if r['lzm_bic_wins'])}/{len(rows)}")
    lines.append(f"- corr(lzm_is_lam, c): {corr(arr('lzm_is_lam'), c_all):.4f}")
    lines.append(f"- corr(lzm_is_lam, hapax_ratio): {corr(arr('lzm_is_lam'), arr('hapax_ratio')):.4f}")
    lines.append(f"- corr(lzm_is_c minus zm_c shift): median {float(np.median(arr('lzm_is_c') - c_all)):.3f}\n")
    lines.append("## Exponent gap (A2 rerun)")
    lines.append(f"- corr(exp_gap, c): {corr(arr('exp_gap'), c_all):.4f}")
    lines.append(f"- corr(exp_gap, is_winner_flag): {corr(arr('exp_gap'), is_flag):.4f}")
    lines.append(f"- corr(c, is_winner_flag): {corr(c_all, is_flag):.4f}\n")
    lines.append("## Lexical confound (A3 rerun)")
    lines.append(f"- corr(hapax, is_winner_flag): {corr(arr('hapax_ratio'), is_flag):.4f}")
    lines.append(f"- corr(ttr, c): {corr(arr('ttr'), c_all):.4f}")
    lines.append(f"- corr(tokens, c): {corr(np.log(arr('tokens')), np.log1p(c_all)):.4f}  (log-log)")
    (OUT / "f1_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nwrote", OUT / "f1_per_corpus.csv", "and f1_summary.md", flush=True)


if __name__ == "__main__":
    sys.exit(main())
