"""F13 — does the fitted PLN mixture CONTAIN the width law?

The s-law (manuscript eq. 3, s ~ 0.012*V) is empirical. If the generative
account of Section 3.5 is the mechanism, the width law should be derivable —
and the first, decisive check is computational: push each corpus's fitted
5-parameter zero-truncated Poisson-lognormal mixture (f4d) through the
observation model to a synthetic rank-frequency curve, fit the SAME 9-parameter
erf-gate model used on real corpora (f5b/f12 fitter), and ask whether the
synthetic seam width s_synth = k*w_tail reproduces:

  (a) the per-corpus measured widths  (corr, median ratio s_synth/s_real)
  (b) the LAW itself                  (log s_synth ~ log V_synth slope ~ 1,
                                       median s_synth/V_synth ~ 0.012)

Pre-registered readings:
  - CONTAINED:  slope ~1, median s/V in [0.008, 0.016], strong per-corpus corr.
    -> an asymptotic derivation of 0.012 from PLN overlap + count floor is
       justified; the width law is a consequence of the generative account.
  - NOT CONTAINED: systematic offset (cf. f4d's naive w_pred, ~3x small) or no
    V-scaling. -> the width law needs an ingredient beyond the single-depth
    mixture; it stands as an independent empirical regularity (and the
    manuscript's separation of Sections 3.3 and 3.5 is vindicated).

Simulation matches f7 exactly: N_total from zero-truncation, component draw,
lognormal rate, Poisson count, drop zeros, sort. 3 replicates per corpus.

Outputs: ../outputs/f13_width_fits.csv, f13_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf, logsumexp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
F4D = REPO / "experiments" / "f4d_poisson_lognormal_mixture" / "outputs" / "f4d_mixture_fits.csv"
F2 = REPO / "experiments" / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv"

SEED = 20260723
N_REPS = 3
N_STARTS = 24
MAX_NFEV = 9000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 5000.0, 10.0, 10.0])

GH_X, GH_W = hermgauss(48)


def pln_p0(mu, sd):
    lam = np.exp(np.clip(mu + math.sqrt(2.0) * sd * GH_X, -30, 30))
    return float(np.exp(logsumexp(-lam + np.log(GH_W)) - 0.5 * math.log(math.pi)))


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
    name, rep, freqs_list = task
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + (hash((name, rep)) % (2**31)))
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
    return {"corpus": name, "rep": rep, "V_synth": V, "rmse": round(best_rmse, 5),
            "k": round(float(best_p[6]), 2), "w_gate": round(float(best_p[7]), 4),
            "w_tail": round(float(best_p[8]), 4),
            "s_synth": round(float(best_p[6] * best_p[8]), 2),
            "k_hit_upper": bool(abs(best_p[6] - UPPER[6]) < 1e-3)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f4d = list(csv.DictReader(open(F4D, newline="", encoding="utf-8")))
    f2 = {r["corpus"]: r for r in csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    rng = np.random.default_rng(SEED)

    tasks = []
    for r in f4d:
        name = r["corpus"]
        pi_h, mu_h, sd_h = float(r["pi_h"]), float(r["mu_h"]), float(r["sd_h"])
        mu_t, sd_t = float(r["mu_t"]), float(r["sd_t"])
        V_obs = int(r["V"])
        p0 = pi_h * pln_p0(mu_h, sd_h) + (1 - pi_h) * pln_p0(mu_t, sd_t)
        n_total = int(round(V_obs / max(1e-9, 1 - p0)))
        for rep in range(N_REPS):
            comp = rng.random(n_total) < pi_h
            loglam = np.where(comp, rng.normal(mu_h, sd_h, n_total), rng.normal(mu_t, sd_t, n_total))
            counts = rng.poisson(np.exp(np.clip(loglam, -30, 30)))
            counts = counts[counts > 0]
            freqs = np.sort(counts)[::-1].astype(float)
            tasks.append((name, rep, freqs.tolist()))
    print(f"{len(tasks)} synthetic fits queued ({len(f4d)} corpora x {N_REPS} reps)", flush=True)

    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['corpus'][:28]:28} rep{res['rep']} V={res['V_synth']:6} "
                  f"k={res['k']:7.1f} s={res['s_synth']:8.1f} s/V={res['s_synth']/res['V_synth']:.4f}",
                  flush=True)

    rows.sort(key=lambda r: (r["corpus"], r["rep"]))
    with open(OUT / "f13_width_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # per-corpus medians vs real widths
    agg = []
    for r4 in f4d:
        name = r4["corpus"]
        sub = [r for r in rows if r["corpus"] == name]
        if not sub:
            continue
        s_med = float(np.median([r["s_synth"] for r in sub]))
        v_med = float(np.median([r["V_synth"] for r in sub]))
        s_real = float(f2[name]["s_at_min"])
        v_real = float(f2[name]["V"])
        agg.append({"corpus": name, "V_real": v_real, "V_synth_med": v_med,
                    "s_real": s_real, "s_synth_med": s_med,
                    "sv_real": s_real / v_real, "sv_synth": s_med / v_med,
                    "ratio": s_med / s_real})

    sr = np.array([a["s_real"] for a in agg])
    ss = np.array([a["s_synth_med"] for a in agg])
    vs = np.array([a["V_synth_med"] for a in agg])
    corr = float(np.corrcoef(np.log(ss), np.log(sr))[0, 1])
    med_ratio = float(np.median(ss / sr))
    med_sv = float(np.median(ss / vs))
    lV, lS = np.log(vs), np.log(ss)
    X = np.column_stack([np.ones_like(lV), lV])
    b, *_ = np.linalg.lstsq(X, lS, rcond=None)
    ssr_ = float(np.sum((lS - X @ b) ** 2))
    r2 = 1 - ssr_ / float(np.sum((lS - lS.mean()) ** 2))

    lines = ["# F13 — does the fitted PLN mixture contain the width law?\n"]
    lines.append("| corpus | V_real | V_synth | s_real | s_synth | s/V real | s/V synth | ratio |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for a in agg:
        lines.append(f"| {a['corpus']} | {a['V_real']:.0f} | {a['V_synth_med']:.0f} | "
                     f"{a['s_real']:.1f} | {a['s_synth_med']:.1f} | {a['sv_real']:.4f} | "
                     f"{a['sv_synth']:.4f} | {a['ratio']:.2f} |")
    lines.append("")
    lines.append(f"- corr(log s_synth, log s_real) = {corr:.3f}")
    lines.append(f"- median s_synth/s_real = {med_ratio:.3f}")
    lines.append(f"- median s_synth/V_synth = {med_sv:.4f}  (real-corpus law: 0.0118)")
    lines.append(f"- law within synthetic family: log s ~ log V slope = {b[1]:.3f}, R2 = {r2:.3f}")
    lines.append("")
    if 0.008 <= med_sv <= 0.016 and 0.8 <= b[1] <= 1.2:
        lines.append("**Reading: CONTAINED — the mixture + observation model reproduce the "
                     "width law; an asymptotic derivation of the 1.2% constant is justified.**")
    else:
        lines.append("**Reading: NOT (fully) CONTAINED — the single-depth mixture does not "
                     "reproduce the width law as measured; the law requires an ingredient "
                     "beyond Section 3.5's generative account.**")
    (OUT / "f13_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
