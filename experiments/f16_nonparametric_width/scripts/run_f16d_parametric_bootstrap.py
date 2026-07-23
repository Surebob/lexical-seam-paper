"""f16d — parametric-bootstrap calibration of the width instrument, in-regime.

f16c's planted grid failed recovery — forensics attribute this to the rig
(planted gate widths 5x outside the empirical regime, gate/tail confounding
by construction, order-statistics structure contaminating the twins). The
decisive in-regime test is the parametric bootstrap:

  For each of the 25 English corpora:
    1. fit the canonical 9-parameter erf model -> truth vector theta_i,
       s_true = k * w_tail;
    2. generate the model's exact expected frequency curve from theta_i
       (src decoupled_prediction — the very code the paper uses);
    3. Poisson-sample token counts at the corpus's own scale, drop zeros
       (3 resamples);
    4. refit each resample with the same fitter -> s_hat;
    5. recovery ratio s_hat / s_true, per corpus and pooled.

If ratios concentrate near 1, the width estimator is calibrated where it is
used, and the s-law's absolute value stands. Systematic bias -> report and
correct. Wide scatter -> the width claim must be restated. Either way it
goes in the paper's calibration appendix.

Outputs: ../outputs/f16d_bootstrap.csv, f16d_summary.md
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

SEED = 20260730
N_BOOT = 3


def full_fit(freqs, seed_salt):
    """Canonical 24-start erf fit returning the full 9-parameter vector."""
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


def task(job):
    kind, name, freqs_list, salt = job
    p, rmse = full_fit(np.asarray(freqs_list), salt)
    return {"kind": kind, "name": name, "params": [float(v) for v in p],
            "rmse": rmse}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    print("fitting truth vectors (25 corpora)...", flush=True)
    truth_jobs = []
    freqs_by = {}
    for i, spec in enumerate(SEARCHED_CORPORA):
        ds = build_zipf_dataset(ZIPF / spec["filename"])
        freqs_by[spec["slug"]] = ds["freqs"]
        truth_jobs.append(("truth", spec["slug"], ds["freqs"].tolist(), i))
    truths = {}
    with Pool(processes=8) as pool:
        for res in pool.imap_unordered(task, truth_jobs, chunksize=1):
            truths[res["name"]] = res["params"]
            p = res["params"]
            print(f"truth {res['name'][:28]:28} s_true={p[6]*p[8]:8.1f}", flush=True)

    print("generating parametric resamples...", flush=True)
    boot_jobs = []
    s_true = {}
    for name, p in truths.items():
        freqs = freqs_by[name]
        V = len(freqs)
        T = float(np.sum(freqs))
        ranks = np.arange(1, V + 1, dtype=np.float64)
        log_ranks = np.log(ranks)
        y_model = fitmod.predict(ranks, log_ranks, np.asarray(p))
        f_model = np.exp(y_model)
        f_model *= T / f_model.sum()
        s_true[name] = p[6] * p[8]
        for b in range(N_BOOT):
            counts = rng.poisson(f_model)
            counts = counts[counts > 0]
            boot = np.sort(counts)[::-1].astype(float)
            boot_jobs.append(("boot", f"{name}|{b}", boot.tolist(), 1000 + b))

    print(f"{len(boot_jobs)} bootstrap refits...", flush=True)
    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(task, boot_jobs, chunksize=1), 1):
            name, b = res["name"].split("|")
            p = res["params"]
            s_hat = p[6] * p[8]
            rows.append({"corpus": name, "boot": int(b),
                         "s_true": round(s_true[name], 1), "s_hat": round(s_hat, 1),
                         "ratio": round(s_hat / s_true[name], 3)})
            print(f"[{i}/{len(boot_jobs)}] {name[:26]:26} b{b} "
                  f"ratio={s_hat/s_true[name]:.3f}", flush=True)

    rows.sort(key=lambda r: (r["corpus"], r["boot"]))
    with open(OUT / "f16d_bootstrap.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    ratios = np.array([r["ratio"] for r in rows])
    per_corpus = {}
    for r in rows:
        per_corpus.setdefault(r["corpus"], []).append(r["ratio"])
    med_per = {k: float(np.median(v)) for k, v in per_corpus.items()}

    lines = ["# f16d — parametric-bootstrap calibration (in-regime)\n"]
    lines.append(f"- {len(rows)} refits ({len(per_corpus)} corpora x {N_BOOT} resamples)")
    lines.append(f"- pooled recovery ratio s_hat/s_true: median {float(np.median(ratios)):.3f}, "
                 f"IQR [{float(np.percentile(ratios,25)):.3f}, {float(np.percentile(ratios,75)):.3f}]")
    lines.append(f"- corpora with median ratio in [0.8, 1.25]: "
                 f"{sum(1 for v in med_per.values() if 0.8 <= v <= 1.25)}/{len(med_per)}")
    lines.append(f"- worst corpora: " + ", ".join(
        f"{k} ({v:.2f})" for k, v in sorted(med_per.items(), key=lambda kv: abs(math.log(max(kv[1], 1e-9))))[-3:]))
    lines.append("")
    med = float(np.median(ratios))
    if 0.85 <= med <= 1.18 and float(np.percentile(ratios, 75)) - float(np.percentile(ratios, 25)) < 0.5:
        lines.append("**Reading: CALIBRATED — the estimator recovers the width in-regime; "
                     "the s-law's absolute value stands as measured. (f16c's grid failure "
                     "was the rig: out-of-regime plants + gate/tail confound.)**")
    else:
        lines.append("**Reading: NOT calibrated in-regime — bias/scatter must be propagated "
                     "into the manuscript's width claims before submission.**")
    (OUT / "f16d_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
