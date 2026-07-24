"""f24 — the giant depth ladder: SPGC (Gerlach & Font-Clos), ~3B tokens.

Where do the depth functions go? Every constant in the paper (c, s/V, lambda)
is the book-depth value of a rising depth function, and our data ends at
tokens-per-type ~40-63. SPGC's 55,905 books let us merge (valid by f12's
composition result) a ladder to tok/V ~hundreds and watch s/V, lambda, and c
at depths no one has measured.

Ladder: English books in a seeded-shuffled order, cumulative merges at
rungs of [10, 30, 100, 300, 1000, 3000, 10000, 30000, ALL] books. Per rung:
  - single ZM fit (b, c)        -> c(depth)
  - free lambda-ZM fit          -> lambda(depth)
  - canonical 9-param erf fit   -> s(depth), s/V(depth)   [24 starts]

Scale handling: above V_FULL_MAX ranks, fits run on a log-spaced rank
subsample with trapezoid quadrature weights (sqrt-w rows), which preserves
the paper's equal-per-rank objective to quadrature accuracy. The first rung
where V <= V_FULL_MAX runs BOTH full and subsampled fits to measure the
subsampling bias in situ.

Incremental CSV + rung-level resume. Outputs: ../outputs/f24_ladder.csv,
f24_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
import zipfile
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SPGC = REPO / "data" / "spgc"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
import run_f13b_basin_reselect as fitmod  # noqa: E402

SEED = 20260734
RUNGS = [10, 30, 100, 300, 1000, 3000, 10000, 30000, -1]  # -1 = all English
V_FULL_MAX = 200_000
SUB_POINTS = 120_000
CSV_PATH = OUT / "f24_ladder.csv"
FIELDS = ["rung_books", "mode", "tokens", "V", "tok_per_V", "zm_b", "zm_c",
          "zm_rmse", "lam", "lzm_rmse", "s", "s_over_V", "erf_rmse",
          "wg", "wt", "k"]


def rank_subsample(V):
    """Log-spaced unique integer ranks + trapezoid quadrature weights."""
    if V <= V_FULL_MAX:
        r = np.arange(1, V + 1, dtype=np.float64)
        return r, np.ones(V)
    r = np.unique(np.geomspace(1, V, SUB_POINTS).astype(np.int64)).astype(np.float64)
    edges = np.empty(len(r) + 1)
    edges[1:-1] = (r[:-1] + r[1:]) / 2
    edges[0] = 0.5
    edges[-1] = V + 0.5
    w = np.diff(edges)
    return r, w


def fit_zm_w(ranks, logf, sw, V):
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 1024)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)]) * sw[:, None]
        y = logf * sw
        coef, *_ = np.linalg.lstsq(d, y, rcond=None)
        mse = float(np.sum((d @ coef - y) ** 2) / np.sum(sw ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c))
    return {"a": best[1], "b": best[2], "c": best[3], "rmse": math.sqrt(best[0])}


def fit_lzm_w(ranks, logf, sw, V):
    x = 0.05 + 0.95 * np.log(ranks) / math.log(V)
    gx = np.exp(x - 1.0) - x
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 1024)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c), gx]) * sw[:, None]
        y = logf * sw
        coef, *_ = np.linalg.lstsq(d, y, rcond=None)
        mse = float(np.sum((d @ coef - y) ** 2) / np.sum(sw ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[2]))
    return {"lam": best[1], "rmse": math.sqrt(best[0])}


def giant_bounds(V):
    """The canonical bounds are book-regime search priors (k <= 5000,
    c <= 1000, b <= 3) that never bind on books but strangle giant corpora
    (f26 caught k pinned at exactly 5000.0 on a 514k-type corpus; the ladder's
    k trajectory 331->842->1373 was headed for the ceiling). A binding bound
    is a measurement artifact by definition, so above book scale the bounds
    scale with the corpus: k <= max(5000, V/4), c <= max(1000, V/10),
    b <= 4. Documented protocol adaptation; the small rungs (k << 5000) are
    unaffected and were kept."""
    lower = fitmod.LOWER.copy()
    upper = fitmod.UPPER.copy()
    upper[2] = upper[5] = max(1000.0, V / 10.0)
    upper[1] = upper[4] = 4.0
    upper[6] = max(5000.0, V / 4.0)
    return lower, upper


def _erf_start(args):
    x0, ranks, log_ranks, y, sw, lower, upper = args
    try:
        sol = least_squares(lambda p: (fitmod.predict(ranks, log_ranks, p) - y) * sw,
                            x0=x0, bounds=(lower, upper),
                            method="trf", max_nfev=9000)
        r = float(np.sqrt(np.sum(sol.fun ** 2) / np.sum(sw ** 2)))
        return r, sol.x
    except Exception:
        return None


def fit_erf_w(ranks, logf, sw, V, seed_salt, procs=8):
    log_ranks = np.log(ranks)
    rng = np.random.default_rng(SEED + seed_salt)
    anchor6 = fitmod.piecewise_anchor(ranks, logf)
    lower, upper = giant_bounds(V)
    k_cap = float(upper[6])
    starts = []
    for i in range(24):
        if i < 12:
            x0 = np.array(anchor6 + [0.0, 0.0, 0.0])
            x0[:6] += rng.normal(0.0, [6.0, 0.3, 90.0, 6.0, 0.3, 90.0])
            x0[6] = min(0.02 * V, 0.8 * k_cap) * math.exp(rng.normal(0.0, 0.7))
            x0[7] = 0.5 * math.exp(rng.normal(0.0, 0.45))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.45))
        else:
            x0 = np.array([
                rng.uniform(logf.min() - 3, logf.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(logf.min() - 3, logf.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(max(k_cap / 2, 40)))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
            ])
        starts.append((np.clip(x0, lower + 1e-8, upper - 1e-8),
                       ranks, log_ranks, logf, sw, lower, upper))
    best = None
    with Pool(processes=procs) as pool:
        for res in pool.imap_unordered(_erf_start, starts, chunksize=1):
            if res is None:
                continue
            if best is None or res[0] < best[0]:
                best = res
    r, p = best
    return {"rmse": r, "k": float(p[6]), "wg": float(p[7]), "wt": float(p[8]),
            "s": float(p[6] * p[8])}


def english_ids():
    ids = []
    with open(SPGC / "SPGC-metadata-2018-07-18.csv", encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            if "'en'" in (row.get("language") or "") and row.get("type") == "Text":
                ids.append(row["id"])
    return ids


def fit_rung(rung_books, counter, tokens, mode, w, fh, done):
    key = (str(rung_books), mode)
    if key in done:
        return
    freqs = np.array(sorted(counter.values(), reverse=True), dtype=np.float64)
    V = len(freqs)
    if mode == "sub":
        r_idx = np.unique(np.geomspace(1, V, SUB_POINTS).astype(np.int64))
        ranks = r_idx.astype(np.float64)
        edges = np.empty(len(ranks) + 1)
        edges[1:-1] = (ranks[:-1] + ranks[1:]) / 2
        edges[0], edges[-1] = 0.5, V + 0.5
        wts = np.diff(edges)
        logf = np.log(freqs[r_idx - 1])
    else:
        ranks = np.arange(1, V + 1, dtype=np.float64)
        wts = np.ones(V)
        logf = np.log(freqs)
    sw = np.sqrt(wts / wts.mean())
    zm = fit_zm_w(ranks, logf, sw, V)
    lzm = fit_lzm_w(ranks, logf, sw, V)
    erf = fit_erf_w(ranks, logf, sw, V, seed_salt=rung_books if rung_books > 0 else 99999)
    row = {"rung_books": rung_books, "mode": mode, "tokens": int(tokens), "V": V,
           "tok_per_V": round(tokens / V, 1), "zm_b": round(zm["b"], 4),
           "zm_c": round(zm["c"], 2), "zm_rmse": round(zm["rmse"], 4),
           "lam": round(lzm["lam"], 2), "lzm_rmse": round(lzm["rmse"], 4),
           "s": round(erf["s"], 1), "s_over_V": round(erf["s"] / V, 5),
           "erf_rmse": round(erf["rmse"], 4), "wg": round(erf["wg"], 4),
           "wt": round(erf["wt"], 4), "k": round(erf["k"], 1)}
    w.writerow(row)
    fh.flush()
    print(f"rung {rung_books:>6} [{mode}] tokens={int(tokens):,} V={V:,} tok/V={row['tok_per_V']} "
          f"c={row['zm_c']} lam={row['lam']} s/V={row['s_over_V']}", flush=True)


def main():
    done = set()
    if CSV_PATH.exists():
        with open(CSV_PATH, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                done.add((r["rung_books"], r["mode"]))
        print(f"resume: {len(done)} rung-fits present", flush=True)

    ids = english_ids()
    rng = np.random.default_rng(SEED)
    rng.shuffle(ids)
    print(f"english books: {len(ids)}", flush=True)

    z = zipfile.ZipFile(SPGC / "SPGC-counts-2018-07-18.zip")
    members = set(z.namelist())

    new_file = not CSV_PATH.exists()
    fh = open(CSV_PATH, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new_file:
        w.writeheader(); fh.flush()

    rungs = [r if r > 0 else len(ids) for r in RUNGS]
    counter = Counter()
    tokens = 0
    merged = 0
    validated_full = False
    for target in rungs:
        while merged < min(target, len(ids)):
            name = f"SPGC-counts-2018-07-18/{ids[merged]}_counts.txt"
            merged += 1
            if name not in members:
                continue
            try:
                for line in z.read(name).decode("utf-8", "ignore").splitlines():
                    parts = line.split("\t")
                    if len(parts) == 2:
                        try:
                            n = int(parts[1])
                        except ValueError:
                            continue
                        counter[parts[0]] += n
                        tokens += n
            except Exception as e:
                print(f"skip {name}: {e}", flush=True)
        V = len(counter)
        label = target if target < len(ids) else -1
        if V <= V_FULL_MAX:
            fit_rung(label, counter, tokens, "full", w, fh, done)
            if not validated_full and V > 50_000:
                fit_rung(label, counter, tokens, "sub", w, fh, done)
                validated_full = True
        else:
            fit_rung(label, counter, tokens, "sub", w, fh, done)
        if merged >= len(ids):
            break
    fh.close()

    rows = list(csv.DictReader(open(CSV_PATH, encoding="utf-8")))
    lines = ["# f24 — giant depth ladder (SPGC, English)\n",
             "| books | mode | tokens | V | tok/V | c | lambda | s/V |",
             "|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in sorted(rows, key=lambda r: (int(r["rung_books"]) if int(r["rung_books"]) > 0 else 10**9, r["mode"])):
        lines.append(f"| {r['rung_books']} | {r['mode']} | {int(r['tokens']):,} | {int(r['V']):,} | "
                     f"{r['tok_per_V']} | {r['zm_c']} | {r['lam']} | {r['s_over_V']} |")
    (OUT / "f24_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
