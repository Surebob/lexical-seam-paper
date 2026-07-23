"""F13c-K — extended basin search for the one f13b holdout, Critique of Pure
Reason (pred-basin width ratio 550 where 8/9 peers healed to ~1.0).

48 Nelder-Mead starts over WIDER ranges than f13b's 16 (pi logit down to -7,
mu_t to -9, sd_t up to 4.5, gap 1-9, sd_h 0.3-2.5). Every admissible basin
(nll <= min + 80) is width-tested (2 simulation reps, erf fit), so we learn
whether ANY basin heals Kant — not just the pred_err winner. Reports the
(pred_err, width-ratio) landscape.

Outputs: ../outputs/f13c_kant.csv, f13c_kant_summary.md
"""
from __future__ import annotations

import csv
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_f13b_basin_reselect import (  # noqa: E402
    counts_for, landscape_task, simulate, erf_task, ZIPF, F2, OUT, NLL_SLACK,
)

SEED = 20260725
N_STARTS = 48
CORPUS = ("Critique of Pure Reason", "pg4280.txt")


def main():
    import csv as _csv
    f2 = {r["corpus"]: r for r in _csv.DictReader(open(F2, newline="", encoding="utf-8"))}
    rng = np.random.default_rng(SEED)
    name, fn = CORPUS
    counts = counts_for(ZIPF / fn)
    print(f"{name}: V={len(counts):,}", flush=True)

    jobs = []
    for _ in range(N_STARTS):
        x0 = [rng.uniform(-7, 0),
              rng.uniform(-9.0, 2.0),
              math.log(rng.uniform(0.5, 4.5)),
              math.log(rng.uniform(1.0, 9.0)),
              math.log(rng.uniform(0.3, 2.5))]
        jobs.append((name, counts.tolist(), x0))
    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(landscape_task, jobs, chunksize=1), 1):
            if res:
                rows.append(res)
                print(f"[{i}/{len(jobs)}] nll={res['nll']:.0f} pi={res['pi_h']:.4f} "
                      f"gap={res['mu_h']-res['mu_t']:.2f} pred_err={res['pred_err']:.3f}", flush=True)

    mn = min(r["nll"] for r in rows)
    adm = [r for r in rows if r["nll"] <= mn + 80.0]
    # dedupe near-identical basins by (pi, gap, sd_t) rounding
    seen, basins = set(), []
    for r in sorted(adm, key=lambda r: r["pred_err"]):
        key = (round(r["pi_h"], 3), round(r["mu_h"] - r["mu_t"], 1), round(r["sd_t"], 1))
        if key in seen:
            continue
        seen.add(key)
        basins.append(r)
    basins = basins[:10]
    print(f"\n{len(basins)} distinct admissible basins to width-test", flush=True)

    erf_jobs = []
    for bi, b in enumerate(basins):
        params = (b["pi_h"], b["mu_h"], b["sd_h"], b["mu_t"], b["sd_t"])
        for rep in range(2):
            freqs = simulate(params, b["n_total"], 1.0, rng)
            erf_jobs.append((f"basin{bi}", rep, freqs.tolist()))
    efits = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, erf_jobs, chunksize=1), 1):
            efits.append(res)
            print(f"[{i}/{len(erf_jobs)}] {res['corpus']} rep{res['rep']} s={res['s_synth']:.1f}", flush=True)

    s_real = float(f2[name]["s_at_min"])
    out_rows = []
    for bi, b in enumerate(basins):
        sub = [e["s_synth"] for e in efits if e["corpus"] == f"basin{bi}"]
        ratio = float(np.median(sub)) / s_real
        out_rows.append({"basin": bi, "nll": round(b["nll"], 1), "dnll": round(b["nll"] - mn, 1),
                         "pi_h": round(b["pi_h"], 5), "gap": round(b["mu_h"] - b["mu_t"], 2),
                         "sd_t": round(b["sd_t"], 2), "pred_err": round(b["pred_err"], 3),
                         "width_ratio": round(ratio, 2)})
    with open(OUT / "f13c_kant.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)

    lines = ["# F13c-K — Kant extended basin search\n"]
    lines.append(f"- s_real = {s_real:.1f}; {len(rows)} fits, {len(basins)} distinct admissible basins")
    lines.append("| basin | dNLL | pi_h | gap | sd_t | pred_err | width ratio |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for r in sorted(out_rows, key=lambda r: r["width_ratio"]):
        lines.append(f"| {r['basin']} | {r['dnll']} | {r['pi_h']} | {r['gap']} | {r['sd_t']} | "
                     f"{r['pred_err']} | {r['width_ratio']} |")
    best = min(out_rows, key=lambda r: abs(math.log(max(r["width_ratio"], 1e-9))))
    lines.append("")
    if 0.5 <= best["width_ratio"] <= 2.0:
        lines.append(f"**A healing basin EXISTS (basin {best['basin']}, ratio {best['width_ratio']}, "
                     f"dNLL {best['dnll']}, pred_err {best['pred_err']}) — Kant joins the other 24; "
                     "the earlier failure was search coverage, not physics.**")
    else:
        lines.append(f"**No admissible basin heals Kant (best ratio {best['width_ratio']}). The "
                     "corpus genuinely deviates from the two-lognormal account — an honest "
                     "outlier to report.**")
    (OUT / "f13c_kant_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
