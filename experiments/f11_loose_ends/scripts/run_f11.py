"""F11 — four loose ends in one sweep.

  A. 2b FAITHFUL rerun (legacy deterministic protocol, ported from
     src/legacy_engines/zipf_synthetic_mixture.py): power-law frequency VECTORS
     (no sampling). Same-exponent control alpha=1.5/1.5 with the legacy sizes
     (A: 50 types @ scale 10, B: 10000 @ 1) + replication of the positive
     small-gap mixture (1.5/1.3). Family winner + does unit step-2 help.
  B. Diachronic: English composition/translation year vs seam quantities
     (s/V from f2, w_gate, c, b). First-ever look.
  C. Heaps' beta per corpus (binomial thinning, 5 scales) + relation to b.
  D. Winner-from-stats classifier: predict canonical is/exp family from
     (log tokens, hapax ratio, TTR), leave-one-out logistic regression.

Outputs: ../outputs/f11_2b.csv, f11_diachronic.csv, f11_heaps.csv, f11_summary.md
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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
F1 = REPO / "experiments" / "f1_fresh_reproduction" / "outputs" / "f1_per_corpus.csv"
F2 = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"
ARCH = Path(r"C:\Users\Greg Kara\Desktop\temporary\emlexperiment\results\s2_v3_windows_full_outputs_2026-04-18\s2_v3_per_fit_results.csv")
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260727

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}

# English-text year (composition, or translation year for translated corpora)
YEARS = {
    "Canterbury Tales": 1400, "Complete Works of Shakespeare": 1600, "King James Bible": 1611,
    "Wealth of Nations": 1776, "Decline and Fall Vol 1": 1776, "Federalist Papers": 1788,
    "Pride and Prejudice": 1813, "Democracy in America": 1838, "Collected Poe": 1845,
    "Jane Eyre": 1847, "Moby Dick": 1851, "Origin of Species": 1859,
    "Critique of Pure Reason": 1855, "Don Quixote": 1885, "Grimm's Fairy Tales": 1884,
    "Les Miserables": 1887, "Complete Sherlock Holmes": 1892, "The Iliad": 1898,
    "Arabian Nights (Vol 1)": 1898, "Principia Ethica": 1903, "War and Peace": 1904,
    "Emile": 1911, "Aesop's Fables": 1912, "Dubliners": 1914, "Ulysses": 1922,
}
TRANSLATED = {"Democracy in America", "Critique of Pure Reason", "Don Quixote",
              "Grimm's Fairy Tales", "Les Miserables", "The Iliad", "Arabian Nights (Vol 1)",
              "War and Peace", "Emile", "Aesop's Fables"}

CORPFILES = {
    "Complete Works of Shakespeare": "pg100.txt", "War and Peace": "pg2600.txt",
    "Moby Dick": "pg2701.txt", "King James Bible": "pg10.txt",
    "Federalist Papers": "pg1404.txt", "Grimm's Fairy Tales": "pg2591.txt",
    "Don Quixote": "pg996.txt", "Pride and Prejudice": "pg1342.txt",
    "Canterbury Tales": "pg2383.txt", "Arabian Nights (Vol 1)": "pg3435.txt",
    "Aesop's Fables": "pg11339.txt", "Complete Sherlock Holmes": "pg1661.txt",
    "Jane Eyre": "pg1260.txt", "Dubliners": "pg2814.txt",
    "The Iliad": "pg6130.txt", "Democracy in America": "pg815.txt",
    "Origin of Species": "pg1228.txt", "Wealth of Nations": "pg3300.txt",
    "Les Miserables": "pg135.txt", "Decline and Fall Vol 1": "pg731.txt",
    "Emile": "pg5427.txt", "Ulysses": "pg4300.txt",
    "Collected Poe": "pg2147.txt", "Principia Ethica": "pg53430.txt",
    "Critique of Pure Reason": "pg4280.txt",
}


def counts_of(fname):
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
    return np.array(sorted(Counter(TOKEN_RE.findall(text[st:e].lower())).values(), reverse=True), dtype=np.int64)


def fit_zm(freqs):
    freqs = np.asarray(freqs, dtype=float)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 768)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, d @ coef, float(-coef[1]), float(c))
    return best[1], math.sqrt(best[0]), best[2], best[3]


def family_check(freqs):
    pred, rmse, b, c = fit_zm(freqs)
    logf = np.log(np.asarray(freqs, dtype=float))
    resid = logf - pred
    logr = np.log(np.arange(1, len(freqs) + 1, dtype=float))
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
    scores = {k: float(np.sqrt(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
    w = min(scores, key=scores.get)
    return {"zm_b": b, "zm_c": c, "zm_rmse": rmse, "winner": w,
            "winner_rmse": scores[w], "helps": scores[w] < rmse}


def power_law(size, alpha, scale):
    return scale / np.power(np.arange(1, size + 1, dtype=float), alpha)


def main():
    rng = np.random.default_rng(SEED)

    # ---------- A. faithful 2b ----------
    a_rows = []
    cases = [
        ("same_exponent_1.5_1.5", np.concatenate([power_law(50, 1.5, 10.0), power_law(10000, 1.5, 1.0)])),
        ("small_gap_1.5_1.3", np.concatenate([power_law(50, 1.5, 10.0), power_law(10000, 1.3, 1.0)])),
        ("pure_single_1.5", power_law(10000, 1.5, 1.0)),
    ]
    for name, freqs in cases:
        freqs = np.sort(freqs)[::-1]
        r = family_check(freqs)
        a_rows.append({"case": name, "V": len(freqs), **r})
        print(f"A {name:24} b={r['zm_b']:.3f} c={r['zm_c']:.2f} rmse={r['zm_rmse']:.4f} "
              f"winner={r['winner']} helps={r['helps']}", flush=True)
    with open(OUT / "f11_2b.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(a_rows[0].keys())); w.writeheader(); w.writerows(a_rows)

    # ---------- B. diachronic ----------
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    erf_fit = {r["corpus"]: r for r in csv.DictReader(open(ARCH, newline="", encoding="utf-8")) if r["gate"] == "erf"}
    f1 = {r["corpus"]: r for r in csv.DictReader(open(F1, newline="", encoding="utf-8"))}
    b_rows = []
    for name, yr in YEARS.items():
        if name not in f2:
            continue
        b_rows.append({"corpus": name, "year": yr, "translated": name in TRANSLATED,
                       "s_over_V": float(f2[name]["s_at_min"]) / float(f2[name]["V"]),
                       "w_gate": float(erf_fit[name]["w_gate"]),
                       "zm_c": float(f1[name]["zm_c"]), "zm_b": float(f1[name]["zm_b"]),
                       "hapax": float(f1[name]["hapax_ratio"])})
    with open(OUT / "f11_diachronic.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(b_rows[0].keys())); w.writeheader(); w.writerows(b_rows)

    def cc(a, b):
        return float(np.corrcoef(a, b)[0, 1])

    yrs = np.array([r["year"] for r in b_rows], dtype=float)
    sV = np.array([r["s_over_V"] for r in b_rows])
    wg = np.array([r["w_gate"] for r in b_rows])
    orig = [i for i, r in enumerate(b_rows) if not r["translated"]]
    dia = {
        "corr(year, s/V) all": cc(yrs, sV),
        "corr(year, w_gate) all": cc(yrs, wg),
        "corr(year, s/V) originals": cc(yrs[orig], sV[orig]),
        "corr(year, w_gate) originals": cc(yrs[orig], wg[orig]),
        "corr(year, hapax) originals": cc(yrs[orig], np.array([b_rows[i]["hapax"] for i in orig])),
    }
    for k, v in dia.items():
        print(f"B {k} = {v:.3f}", flush=True)

    # ---------- C. Heaps ----------
    c_rows = []
    for name, fname in CORPFILES.items():
        counts = counts_of(fname)
        T0 = counts.sum()
        Ts, Vs = [], []
        for s in (1.0, 0.5, 0.25, 0.125, 0.0625):
            thin = counts if s == 1.0 else rng.binomial(counts, s)
            Vs.append(int((thin > 0).sum()))
            Ts.append(float(T0 * s))
        X = np.column_stack([np.ones(5), np.log(Ts)])
        beta = float(np.linalg.lstsq(X, np.log(Vs), rcond=None)[0][1])
        c_rows.append({"corpus": name, "heaps_beta": beta, "zm_b": float(f1[name]["zm_b"])})
        print(f"C {name[:28]:29} heaps_beta={beta:.3f} b={float(f1[name]['zm_b']):.3f}", flush=True)
    with open(OUT / "f11_heaps.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(c_rows[0].keys())); w.writeheader(); w.writerows(c_rows)
    hb = np.array([r["heaps_beta"] for r in c_rows])
    bb = np.array([r["zm_b"] for r in c_rows])
    heaps_corr = cc(hb, 1.0 / bb)

    # ---------- D. winner classifier (LOO logistic) ----------
    feats, ys = [], []
    for name in CORPFILES:
        r = f1[name]
        feats.append([math.log(float(r["tokens"])), float(r["hapax_ratio"]), float(r["ttr"])])
        ys.append(1.0 if r["canon_family"] == "is" else 0.0)
    Xf = np.array(feats)
    Xf = (Xf - Xf.mean(0)) / Xf.std(0)
    Xf = np.column_stack([np.ones(len(ys)), Xf])
    yv = np.array(ys)

    def nll_log(wv, Xa, ya):
        z = Xa @ wv
        return float(np.sum(np.log1p(np.exp(-np.abs(z))) + np.maximum(z, 0) - ya * z) + 1e-3 * np.sum(wv[1:]**2))

    correct = 0
    for i in range(len(yv)):
        mask = np.ones(len(yv), bool); mask[i] = False
        res = minimize(nll_log, x0=np.zeros(Xf.shape[1]), args=(Xf[mask], yv[mask]), method="BFGS")
        p = 1.0 / (1.0 + math.exp(-float(Xf[i] @ res.x)))
        correct += int((p > 0.5) == (yv[i] > 0.5))
    acc = correct / len(yv)
    print(f"D winner classifier LOO accuracy = {acc:.3f} ({correct}/{len(yv)})", flush=True)

    lines = ["# F11 — loose ends\n"]
    lines.append("## A. Faithful 2b (deterministic legacy protocol)")
    for r in a_rows:
        lines.append(f"- {r['case']}: b={r['zm_b']:.3f} c={r['zm_c']:.2f} ZM rmse={r['zm_rmse']:.4f} "
                     f"winner={r['winner']} unit-step2 helps={r['helps']}")
    lines.append("\n## B. Diachronic (first look)")
    for k, v in dia.items():
        lines.append(f"- {k} = {v:.3f}")
    lines.append(f"  (n={len(b_rows)}; originals n={len(orig)}; years approximate, translations use translation year)")
    lines.append(f"\n## C. Heaps beta: median {float(np.median(hb)):.3f} "
                 f"(range {hb.min():.3f}-{hb.max():.3f}); corr(beta, 1/b) = {heaps_corr:.3f}")
    lines.append(f"\n## D. Winner-from-stats classifier: LOO accuracy {acc:.3f} "
                 f"({correct}/{len(yv)}) from (log tokens, hapax, TTR) — no curve fitting")
    (OUT / "f11_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
