"""F12 — forced text-mixing test (adversarial check of Williams et al. 2015).

Williams, Bagrow, Danforth & Dodds (PRE 91, 052811, 2015) argue that two-regime
structure in rank-frequency curves is an artifact of AGGREGATING texts: mixing
produces an effective decay of word introduction with a scaling break at rank
b ~ N_avg (the mean per-text vocabulary), so the break should be corpus-relative
and grow more severe with aggregation.

If the lexical seam were that mixing break, then forcibly concatenating m
independent single-author works (m = 1, 2, 4, 7, 14) should (a) move the fitted
seam location k toward N_avg (~10-15k) and (b) break the width law s = k*w_tail
away from s ~ 0.012*V. If instead the seam is a property of usage concentration
that every text carries, s/V should stay on the 0.012 line at every m and k
should stay orders of magnitude below N_avg.

Design: 14 single continuous works (one author, one narrative/argument — no
collections), fixed seeded shuffle, disjoint groups per m:
  m=1: 14 fits, m=2: 7 disjoint pairs, m=4: 3 disjoint quadruples,
  m=7: 2 disjoint septuples, m=14: 1 fit. 27 fits, identical fitter throughout
(erf gate, 9-param decoupled model, same bounds as f5b with k upper 5000).

Outputs: ../outputs/f12_mixtures.csv, f12_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ZIPF = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"

SEED = 20260722
N_STARTS = 24
MAX_NFEV = 9000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 5000.0, 10.0, 10.0])

# Single continuous works only (no collections/anthologies/multi-book scriptures).
SINGLES = [
    ("moby_dick", "pg2701.txt"),
    ("war_and_peace", "pg2600.txt"),
    ("don_quixote", "pg996.txt"),
    ("pride_and_prejudice", "pg1342.txt"),
    ("jane_eyre", "pg1260.txt"),
    ("the_iliad", "pg6130.txt"),
    ("democracy_in_america", "pg815.txt"),
    ("origin_of_species", "pg1228.txt"),
    ("wealth_of_nations", "pg3300.txt"),
    ("les_miserables", "pg135.txt"),
    ("decline_and_fall_vol1", "pg731.txt"),
    ("emile", "pg5427.txt"),
    ("ulysses", "pg4300.txt"),
    ("critique_of_pure_reason", "pg4280.txt"),
]

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]


def tokens_for(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    start, end = 0, len(text)
    for m in START_MARKERS:
        i = text.find(m)
        if i != -1:
            j = text.find("\n", i)
            start = j + 1 if j != -1 else i
            break
    for m in END_MARKERS:
        i = text.find(m)
        if i != -1:
            end = i
            break
    return TOKEN_RE.findall(text[start:end].lower())


def predict(ranks, log_ranks, params):
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    z = (log_ranks - math.log(k)) / w_gate
    sigma = 0.5 * (1.0 - sp_erf(z))
    head = a1 - b1 * np.log(ranks + max(c1, 0.0))
    scale = max(1.0, k * w_tail)
    zz = np.clip((ranks - k) / scale, -60.0, 60.0)
    tail_rank = 1.0 + scale * np.log1p(np.exp(zz))
    tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
    return sigma * head + (1.0 - sigma) * tail


def piecewise_anchor(ranks, logf):
    k = min(500, len(ranks) // 3)

    def zmfit(rr, yy):
        best = None
        for c in np.concatenate([[0.0], np.geomspace(1e-6, rr[-1], 384)]):
            d = np.column_stack([np.ones_like(yy), np.log(rr + c)])
            coef, *_ = np.linalg.lstsq(d, yy, rcond=None)
            mse = float(np.mean((d @ coef - yy) ** 2))
            if best is None or mse < best[0]:
                best = (mse, float(coef[0]), float(-coef[1]), float(c))
        return best[1:]

    a1, b1, c1 = zmfit(ranks[:k], logf[:k])
    rr = np.arange(1, len(ranks) - k + 1, dtype=np.float64)
    a2, b2, c2 = zmfit(rr, logf[k:])
    return [a1, min(max(b1, 0.5), 3.0), min(c1, 1000.0), a2, min(max(b2, 0.5), 3.0), min(c2, 1000.0)]


def fit_task(task):
    name, m, n_avg, freqs_list = task
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + hash(name) % (2**31))
    anchor6 = piecewise_anchor(ranks, y)

    def resid(p):
        return predict(ranks, log_ranks, p) - y

    best_rmse, best_p = None, None
    for i in range(N_STARTS):
        if i < N_STARTS // 2:
            x0 = np.array(anchor6 + [0.0, 0.0, 0.0])
            x0[:6] += rng.normal(0.0, [6.0, 0.3, 90.0, 6.0, 0.3, 90.0])
            x0[6] = min(0.02 * V, 4000.0) * math.exp(rng.normal(0.0, 0.7))
            x0[7] = 0.5 * math.exp(rng.normal(0.0, 0.45))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.45))
        else:
            x0 = np.array([
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(5000))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
            ])
        x0 = np.clip(x0, LOWER + 1e-8, UPPER - 1e-8)
        try:
            sol = least_squares(resid, x0=x0, bounds=(LOWER, UPPER), method="trf", max_nfev=MAX_NFEV)
        except Exception:
            continue
        r = float(np.sqrt(np.mean(sol.fun**2)))
        if best_rmse is None or r < best_rmse:
            best_rmse, best_p = r, sol.x
    return {"mixture": name, "m": m, "V": V, "N_avg": round(n_avg, 1),
            "rmse": round(best_rmse, 5),
            "k": round(float(best_p[6]), 2), "w_gate": round(float(best_p[7]), 4),
            "w_tail": round(float(best_p[8]), 4),
            "s": round(float(best_p[6] * best_p[8]), 2),
            "s_over_V": round(float(best_p[6] * best_p[8] / V), 5),
            "k_over_Navg": round(float(best_p[6]) / n_avg, 4),
            "k_hit_upper": bool(abs(best_p[6] - UPPER[6]) < 1e-3)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # Tokenize all singles once; record per-text vocab for N_avg.
    tok = {}
    vocab = {}
    for slug, fn in SINGLES:
        path = ZIPF / fn
        if not path.exists():
            print(f"MISSING {fn} — aborting", flush=True)
            return 1
        t = tokens_for(path)
        tok[slug] = t
        vocab[slug] = len(set(t))
        print(f"loaded {slug:24} tokens={len(t):9,} V={vocab[slug]:6,}", flush=True)

    order = [s for s, _ in SINGLES]
    rng = np.random.default_rng(SEED)
    rng.shuffle(order)

    def groups_of(m):
        n_groups = len(order) // m
        return [order[i * m:(i + 1) * m] for i in range(n_groups)]

    tasks = []
    for m in (1, 2, 4, 7, 14):
        for gi, group in enumerate(groups_of(m)):
            merged = Counter()
            for slug in group:
                merged.update(tok[slug])
            freqs = sorted(merged.values(), reverse=True)
            n_avg = float(np.mean([vocab[s] for s in group]))
            name = group[0] if m == 1 else f"mix{m}_{chr(97 + gi)}"
            tasks.append((name, m, n_avg, freqs))
    print(f"\n{len(tasks)} fits queued", flush=True)

    rows = []
    with Pool(processes=10) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['mixture']:24} m={res['m']:2} V={res['V']:6} "
                  f"k={res['k']:8.1f} s={res['s']:8.1f} s/V={res['s_over_V']:.4f} "
                  f"k/Navg={res['k_over_Navg']:.3f}", flush=True)

    rows.sort(key=lambda r: (r["m"], r["mixture"]))
    with open(OUT / "f12_mixtures.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lines = ["# F12 — forced text-mixing test (Williams et al. 2015 adversarial check)\n"]
    lines.append("| mixture | m | V | N_avg | k | s | s/V | k/N_avg | rmse |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        lines.append(f"| {r['mixture']} | {r['m']} | {r['V']} | {r['N_avg']} | {r['k']} | "
                     f"{r['s']} | {r['s_over_V']} | {r['k_over_Navg']} | {r['rmse']} |")
    for m in (1, 2, 4, 7, 14):
        sub = [r for r in rows if r["m"] == m]
        sv = np.array([r["s_over_V"] for r in sub])
        kn = np.array([r["k_over_Navg"] for r in sub])
        lines.append(f"\n- m={m}: n={len(sub)}, median s/V={float(np.median(sv)):.4f} "
                     f"(range {sv.min():.4f}-{sv.max():.4f}), median k/N_avg={float(np.median(kn)):.3f}")
    # s-vs-V regression across all mixtures
    lV = np.log([r["V"] for r in rows])
    lS = np.log([r["s"] for r in rows])
    X = np.column_stack([np.ones_like(lV), lV])
    b, *_ = np.linalg.lstsq(X, lS, rcond=None)
    ssr = float(np.sum((lS - X @ b) ** 2))
    r2 = 1 - ssr / float(np.sum((lS - lS.mean()) ** 2))
    lines.append(f"\n- s ~ V regression across all 27 mixtures: beta={b[1]:.4f}, R2={r2:.4f}")
    (OUT / "f12_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
