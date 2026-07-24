"""f18 — NULL-width calibration: what does the canonical width instrument report
on data that has NO seam?

The dual of f16d. f16d planted the fitted two-regime truth and showed the
estimator recovers s (~1.03x). f18 plants the SEAMLESS null: for each of the
25 corpora, fit a single Zipf-Mandelbrot curve (a, b, c — no gate, no second
population), generate its exact expected curve at the corpus's own token count,
Poisson-sample it (3 resamples), and push the samples through the identical
canonical 24-start 9-parameter erf fitter.

The question this kills (or confirms): "maybe s/V ~ 0.0118 is just what the
fitter returns on any Zipfian curve" — i.e. the width law as operator artifact.

Predictions if the law is real:
  - null s/V should be dispersed / degenerate / bound-pinned, NOT clustered
    in the language band around 0.0118;
  - on nulls the two fitted regimes should collapse toward one slope
    (|b_head - b_tail| small) since there is no second population to find;
  - paired real fits (same fitter, same seeds, same pipeline, run here for
    direct comparison) should cluster at 0.0118 with distinct slopes.

Incremental CSV writes + resume (background-process-death protection).
Outputs: ../outputs/f18_null.csv, f18_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
ZIPF = REPO / "data" / "zipf"
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
import run_f13b_basin_reselect as fitmod  # noqa: E402  (canonical erf machinery)
sys.path.insert(0, str(REPO / "src" / "s2_decoupled"))
from shared.corpus_loader import build_zipf_dataset  # noqa: E402
from shared.fit_config import SEARCHED_CORPORA  # noqa: E402

SEED = 20260731
N_NULL = 3
CSV_PATH = OUT / "f18_null.csv"
FIELDS = ["kind", "corpus", "rep", "V", "s_hat", "s_over_V", "b_head", "b_tail",
          "d_slope", "w_gate", "w_tail", "wg_pinned", "wt_pinned", "rmse"]


def fit_single_zm(freqs):
    """Seamless null truth: least-squares single ZM on the log rank curve
    (the f3 protocol: c-grid + linear lstsq for a, b)."""
    freqs = np.asarray(freqs, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([np.array([0.0]), np.geomspace(1e-6, V, 2048)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c))
    return {"a": best[1], "b": best[2], "c": best[3]}


def full_fit(freqs, seed_salt):
    """Canonical 24-start erf fit (identical to f16d)."""
    from scipy.optimize import least_squares
    freqs = np.asarray(freqs, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + seed_salt)
    anchor6 = fitmod.piecewise_anchor(ranks, y)
    best_rmse, best_p = None, None
    for i in range(24):
        if i < 12:
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
        x0 = np.clip(x0, fitmod.LOWER + 1e-8, fitmod.UPPER - 1e-8)
        try:
            sol = least_squares(lambda p: fitmod.predict(ranks, log_ranks, p) - y,
                                x0=x0, bounds=(fitmod.LOWER, fitmod.UPPER),
                                method="trf", max_nfev=9000)
        except Exception:
            continue
        r = float(np.sqrt(np.mean(sol.fun ** 2)))
        if best_rmse is None or r < best_rmse:
            best_rmse, best_p = r, sol.x
    return best_p, best_rmse


def near(v, bound, tol=1e-3):
    return abs(v - bound) <= tol * max(abs(bound), 1.0)


def task(job):
    kind, name, rep, freqs_list, salt = job
    freqs = np.asarray(freqs_list)
    p, rmse = full_fit(freqs, salt)
    V = len(freqs)
    wg, wt = float(p[7]), float(p[8])
    s_hat = float(p[6]) * wt
    lo_wg, hi_wg = float(fitmod.LOWER[7]), float(fitmod.UPPER[7])
    lo_wt, hi_wt = float(fitmod.LOWER[8]), float(fitmod.UPPER[8])
    return {
        "kind": kind, "corpus": name, "rep": rep, "V": V,
        "s_hat": round(s_hat, 1), "s_over_V": round(s_hat / V, 5),
        "b_head": round(float(p[1]), 3), "b_tail": round(float(p[4]), 3),
        "d_slope": round(abs(float(p[1]) - float(p[4])), 3),
        "w_gate": round(wg, 4), "w_tail": round(wt, 4),
        "wg_pinned": near(wg, lo_wg) or near(wg, hi_wg),
        "wt_pinned": near(wt, lo_wt) or near(wt, hi_wt),
        "rmse": round(rmse, 5),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    done = set()
    if CSV_PATH.exists():
        with open(CSV_PATH, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                done.add((r["kind"], r["corpus"], int(r["rep"])))
        print(f"resume: {len(done)} rows already present", flush=True)

    jobs = []
    salt = 0
    for spec in SEARCHED_CORPORA:
        ds = build_zipf_dataset(ZIPF / spec["filename"])
        freqs = np.asarray(ds["freqs"], dtype=np.float64)
        name = spec["slug"]
        salt += 1
        if ("real", name, 0) not in done:
            jobs.append(("real", name, 0, freqs.tolist(), salt))
        zm = fit_single_zm(freqs)
        V = len(freqs)
        T = float(freqs.sum())
        ranks = np.arange(1, V + 1, dtype=np.float64)
        f_model = np.exp(zm["a"] - zm["b"] * np.log(ranks + zm["c"]))
        f_model *= T / f_model.sum()
        for b in range(N_NULL):
            salt += 1
            if ("null", name, b) in done:
                rng.poisson(f_model)  # burn draw to keep stream aligned on resume
                continue
            counts = rng.poisson(f_model)
            counts = counts[counts > 0]
            null = np.sort(counts)[::-1].astype(float)
            jobs.append(("null", name, b, null.tolist(), salt))
        print(f"prepared {name} (zm b={zm['b']:.3f} c={zm['c']:.1f})", flush=True)

    print(f"{len(jobs)} fits to run...", flush=True)
    new_file = not CSV_PATH.exists()
    with open(CSV_PATH, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()
            f.flush()
        with Pool(processes=8) as pool:
            for i, res in enumerate(pool.imap_unordered(task, jobs, chunksize=1), 1):
                w.writerow(res)
                f.flush()
                print(f"[{i}/{len(jobs)}] {res['kind']:4} {res['corpus'][:24]:24} "
                      f"s/V={res['s_over_V']:.4f} d_slope={res['d_slope']:.3f} "
                      f"pin={res['wg_pinned'] or res['wt_pinned']}", flush=True)

    rows = list(csv.DictReader(open(CSV_PATH, encoding="utf-8")))
    for r in rows:
        r["s_over_V"] = float(r["s_over_V"]); r["d_slope"] = float(r["d_slope"])
        r["pin"] = (r["wg_pinned"] == "True") or (r["wt_pinned"] == "True")
    real = [r for r in rows if r["kind"] == "real"]
    null = [r for r in rows if r["kind"] == "null"]
    band = (0.009, 0.015)
    in_band = lambda rs: sum(1 for r in rs if band[0] <= r["s_over_V"] <= band[1])
    sv_r = np.array([r["s_over_V"] for r in real])
    sv_n = np.array([r["s_over_V"] for r in null])
    lines = ["# f18 — null-width calibration (seamless single-ZM truth)\n",
             f"Paired design, identical fitter/seeds: {len(real)} real fits, {len(null)} null fits "
             f"({N_NULL} Poisson resamples per corpus from the fitted single-ZM curve).\n",
             f"## Real corpora (positive control)",
             f"- s/V: median {np.median(sv_r):.4f}, IQR [{np.percentile(sv_r,25):.4f}, {np.percentile(sv_r,75):.4f}]",
             f"- in language band [{band[0]}, {band[1]}]: {in_band(real)}/{len(real)}",
             f"- median |b_head - b_tail|: {np.median([r['d_slope'] for r in real]):.3f}",
             f"- width-bound pins: {sum(1 for r in real if r['pin'])}/{len(real)}\n",
             f"## Seamless nulls (the test)",
             f"- s/V: median {np.median(sv_n):.4f}, IQR [{np.percentile(sv_n,25):.4f}, {np.percentile(sv_n,75):.4f}], "
             f"range [{sv_n.min():.4f}, {sv_n.max():.4f}]",
             f"- in language band [{band[0]}, {band[1]}]: {in_band(null)}/{len(null)}",
             f"- median |b_head - b_tail|: {np.median([r['d_slope'] for r in null]):.3f}",
             f"- width-bound pins: {sum(1 for r in null if r['pin'])}/{len(null)}\n"]
    frac_null_band = in_band(null) / max(len(null), 1)
    frac_real_band = in_band(real) / max(len(real), 1)
    if frac_null_band < 0.3 and frac_real_band > 0.7:
        lines.append("**Reading: the language band is a property of the DATA, not the fitter — "
                     "seamless nulls do not reproduce the width law.**")
    elif frac_null_band >= 0.5:
        lines.append("**Reading: WARNING — nulls cluster in the language band; the width law "
                     "may be partially operator-induced. Escalate before submission.**")
    else:
        lines.append("**Reading: intermediate — inspect distributions before claiming either way.**")
    (OUT / "f18_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
