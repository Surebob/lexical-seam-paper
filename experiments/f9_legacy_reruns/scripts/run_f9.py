"""F9 — three legacy blocked items in one sweep.

  A. 2b unblock: the same-exponent synthetic control at alpha=1.5/1.5 (manuscript
     cites it; only an alpha=0.8/0.8 run was ever saved). Claim to verify: a pure
     same-exponent mixture is a clean single power law — good ZM fit, NO helpful
     step-2 correction.
  B. E2-lite (WLS winner stability): for all 25 English corpora, refit ZM under
     OLS / rank-weighted / frequency-weighted objectives and re-run the 4-generator
     family comparison (unit amplitude) on each residual. Does winner identity
     change under weighting? (v5.1 claimed "no on all 25" from a bundle that never
     actually re-ran the comparison.)
  C. 10a unblock: the two missing poly-transfer rows (War&Peace→Bible, Moby→Bible):
     fit degree-5 polynomial to source corpus residual-after-step2, apply to
     target, report RMSE vs target's own ZM baseline.

Outputs: ../outputs/f9_2b_control.csv, f9_wls_winners.csv, f9_poly_transfer.csv,
         f9_summary.md
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
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20260726

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

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


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
    return np.array(sorted(Counter(TOKEN_RE.findall(text[st:e].lower())).values(), reverse=True), dtype=np.float64)


def fit_zm(ranks, logf, w=None):
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, ranks[-1], 1024)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        if w is None:
            coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        else:
            sw = np.sqrt(w)
            coef, *_ = np.linalg.lstsq(d * sw[:, None], logf * sw, rcond=None)
        pred = d @ coef
        if w is None:
            key = float(np.mean((pred - logf) ** 2))
        else:
            key = float(np.average((pred - logf) ** 2, weights=w))
        if best is None or key < best[0]:
            best = (key, pred, float(-coef[1]), float(c))
    return best[1], best[2], best[3]


def x_of(V):
    logr = np.log(np.arange(1, V + 1, dtype=float))
    return 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])


def family_winner(resid, x, w=None):
    scores = {}
    for k, g in GENS.items():
        d = resid - g(x)
        scores[k] = float(np.mean(d * d)) if w is None else float(np.average(d * d, weights=w))
    return min(scores, key=scores.get)


def main():
    rng = np.random.default_rng(SEED)

    # ---------- A. 2b same-exponent control (alpha = 1.5/1.5) ----------
    a_rows = []
    for rep in range(3):
        lamA = 1.0 + rng.pareto(1.5, 12000)
        lamB = 1.0 + rng.pareto(1.5, 12000)
        counts = rng.poisson(np.concatenate([lamA, lamB]) * 3.0)
        counts = counts[counts > 0]
        freqs = np.sort(counts)[::-1].astype(float)
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=float)
        logf = np.log(freqs)
        pred, b, c = fit_zm(ranks, logf)
        resid = logf - pred
        x = x_of(V)
        rmse = float(np.sqrt(np.mean(resid**2)))
        w = family_winner(resid, x)
        helps = float(np.sqrt(np.mean((resid - GENS[w](x)) ** 2))) < rmse
        a_rows.append({"rep": rep, "V": V, "zm_b": b, "zm_c": c, "zm_rmse": rmse,
                       "winner": w, "step2_helps": helps})
        print(f"A 2b control rep{rep}: V={V} b={b:.3f} c={c:.2f} rmse={rmse:.4f} "
              f"winner={w} helps={helps}", flush=True)
    with open(OUT / "f9_2b_control.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(a_rows[0].keys())); wcsv.writeheader(); wcsv.writerows(a_rows)

    # ---------- B. E2-lite WLS winner stability ----------
    b_rows = []
    for name, fname in CORPORA:
        freqs = load_freqs(fname)
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=float)
        logf = np.log(freqs)
        x = x_of(V)
        winners = {}
        for wname, wvec in [("ols", None), ("rank_w", 1.0 / ranks), ("freq_w", freqs)]:
            pred, b, c = fit_zm(ranks, logf, wvec)
            resid = logf - pred
            winners[wname] = family_winner(resid, x, wvec)
        stable = len(set(winners.values())) == 1
        b_rows.append({"corpus": name, **winners, "stable": stable})
        print(f"B {name[:30]:31} ols={winners['ols']:6} rank_w={winners['rank_w']:6} "
              f"freq_w={winners['freq_w']:6} stable={stable}", flush=True)
    with open(OUT / "f9_wls_winners.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(b_rows[0].keys())); wcsv.writeheader(); wcsv.writerows(b_rows)

    # ---------- C. 10a poly-transfer rows ----------
    c_rows = []
    def resid_after_step2(fname, gen="is"):
        freqs = load_freqs(fname)
        V = len(freqs)
        ranks = np.arange(1, V + 1, dtype=float)
        logf = np.log(freqs)
        pred, _, _ = fit_zm(ranks, logf)
        x = x_of(V)
        return logf - pred - GENS[gen](x), x, float(np.sqrt(np.mean((logf - pred) ** 2)))

    for src_name, src_file, tgt_name, tgt_file in [
        ("War and Peace", "pg2600.txt", "King James Bible", "pg10.txt"),
        ("Moby Dick", "pg2701.txt", "King James Bible", "pg10.txt"),
        ("Complete Works of Shakespeare", "pg100.txt", "War and Peace", "pg2600.txt"),
    ]:
        r_src, x_src, _ = resid_after_step2(src_file)
        coefs = np.polyfit(x_src, r_src, 5)
        r_tgt, x_tgt, zm_rmse_tgt = resid_after_step2(tgt_file)
        transfer_rmse = float(np.sqrt(np.mean((r_tgt - np.polyval(coefs, x_tgt)) ** 2)))
        indom = np.polyfit(x_tgt, r_tgt, 5)
        indom_rmse = float(np.sqrt(np.mean((r_tgt - np.polyval(indom, x_tgt)) ** 2)))
        c_rows.append({"source": src_name, "target": tgt_name, "zm_rmse_target": zm_rmse_tgt,
                       "transfer_rmse": transfer_rmse, "indomain_rmse": indom_rmse,
                       "transfer_beats_zm": transfer_rmse < zm_rmse_tgt})
        print(f"C {src_name[:18]} -> {tgt_name[:18]}: transfer {transfer_rmse:.4f} "
              f"(in-domain {indom_rmse:.4f}, ZM {zm_rmse_tgt:.4f})", flush=True)
    with open(OUT / "f9_poly_transfer.csv", "w", newline="", encoding="utf-8") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(c_rows[0].keys())); wcsv.writeheader(); wcsv.writerows(c_rows)

    n_stable = sum(1 for r in b_rows if r["stable"])
    lines = ["# F9 — legacy reruns\n"]
    lines.append(f"## A. 2b same-exponent control (alpha=1.5/1.5, 3 reps)")
    for r in a_rows:
        lines.append(f"- rep{r['rep']}: ZM rmse {r['zm_rmse']:.4f}, winner {r['winner']}, "
                     f"unit-amplitude helps: {r['step2_helps']}")
    lines.append(f"\n## B. E2-lite: winner stability under weighting — stable on {n_stable}/25")
    ch = [r for r in b_rows if not r["stable"]]
    for r in ch:
        lines.append(f"- CHANGES: {r['corpus']}: ols={r['ols']} rank_w={r['rank_w']} freq_w={r['freq_w']}")
    lines.append(f"\n## C. 10a transfers")
    for r in c_rows:
        lines.append(f"- {r['source']} → {r['target']}: transfer {r['transfer_rmse']:.4f} vs "
                     f"ZM {r['zm_rmse_target']:.4f} vs in-domain {r['indomain_rmse']:.4f} "
                     f"(beats ZM: {r['transfer_beats_zm']})")
    (OUT / "f9_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
