"""F14 — the Gerlach-Altmann model cannot produce the width law.

GA (PRX 3, 021006, 2013) generative model, implemented from their Eqs. (3)-(5):
  - each step: new word-type with prob p_new, else repeat an existing type with
    prob proportional to its past count (uniform draw over past tokens);
  - a new type is a core word with prob p_c = p_c0 (0.99) while N_c < Nc_max,
    else noncore;
  - after each new NONCORE type: p_new <- p_new * (1 - alpha/(N_noncore + s)),
    s = Nc_max  (their Eq. 4; asymptotically p_new ~ N^-alpha, their Eq. 5).
Their correspondence: b = Nc_max (core size), gamma = alpha + 1, with (b,
gamma) language-dependent constants — English crossover rank ~ 7.9k.

Test: grow two independent realizations to 3.2M tokens; snapshot the count
vector at 0.2/0.4/0.8/1.6/3.2M (prefix growth, same protocol as our f6
slicing); fit the SAME 9-parameter erf-gate model; track k(M), s(M), s/V(M).

Discriminator: real corpora hold s/V ~ 0.012 across V (3.8k-63k, f2/f5/f12).
A fixed-core model pins its crossover near Nc_max: as V grows with M, its
width cannot stay a constant fraction of V. If s/V drifts systematically with
V (and k tracks Nc_max rather than V), the GA model class is incompatible
with the measured width law.

Outputs: ../outputs/f14_ga_fits.csv, f14_summary.md
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
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
from run_f13b_basin_reselect import erf_task  # noqa: E402  (same fitter as f5b/f12/f13)

SEED = 20260725
NC_MAX = 7900        # English core size (their b; crossover rank ~7.9k)
ALPHA = 0.77         # gamma = 1.77
PC0 = 0.99
PNEW0 = 0.10
SNAPSHOTS = [200_000, 400_000, 800_000, 1_600_000, 3_200_000]
N_RUNS = 2


def simulate_ga(m_max, rng):
    """Return list of (M, counts_dict_snapshot) at SNAPSHOTS."""
    token_types = np.empty(m_max, dtype=np.int32)
    counts = {}
    p_new = PNEW0
    n_types = 0
    n_core = 0
    n_noncore = 0
    snaps = []
    snap_i = 0
    for m in range(m_max):
        if n_types == 0 or rng.random() < p_new:
            # new word-type
            t = n_types
            n_types += 1
            if n_core < NC_MAX and rng.random() < PC0:
                n_core += 1
            else:
                n_noncore += 1
                p_new *= (1.0 - ALPHA / (n_noncore + NC_MAX))
            counts[t] = 1
        else:
            # repeat proportional to past counts: copy a uniform past token
            t = int(token_types[rng.integers(0, m)])
            counts[t] += 1
        token_types[m] = t
        if snap_i < len(SNAPSHOTS) and (m + 1) == SNAPSHOTS[snap_i]:
            freqs = np.sort(np.fromiter(counts.values(), dtype=np.int64))[::-1]
            snaps.append((m + 1, freqs.astype(float).copy(), n_types, p_new))
            snap_i += 1
    return snaps


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = []
    meta = {}
    for run in range(N_RUNS):
        rng = np.random.default_rng(SEED + run)
        print(f"simulating GA run {run} to {SNAPSHOTS[-1]/1e6:.1f}M tokens...", flush=True)
        for M, freqs, n_types, p_new in simulate_ga(SNAPSHOTS[-1], rng):
            name = f"run{run}_M{M//1000}k"
            meta[name] = {"run": run, "M": M, "V": len(freqs), "p_new_end": p_new}
            jobs.append((name, 0, freqs.tolist()))
            print(f"  run{run} M={M/1e6:.1f}M V={len(freqs):,} p_new={p_new:.4f}", flush=True)

    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            m = meta[res["corpus"]]
            rows.append({"name": res["corpus"], "run": m["run"], "M": m["M"], "V": res["V_synth"],
                         "s": round(res["s_synth"], 1),
                         "s_over_V": round(res["s_synth"] / res["V_synth"], 5),
                         "rmse": round(res["rmse"], 5)})
            print(f"[{i}/{len(jobs)}] {res['corpus']:16} V={res['V_synth']:7,} "
                  f"s={res['s_synth']:9.1f} s/V={res['s_synth']/res['V_synth']:.4f}", flush=True)

    rows.sort(key=lambda r: (r["run"], r["M"]))
    with open(OUT / "f14_ga_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    lines = ["# F14 — width behavior of the Gerlach-Altmann model\n"]
    lines.append(f"- params: Nc_max={NC_MAX}, alpha={ALPHA} (gamma=1.77), pc0={PC0}, pnew0={PNEW0}; "
                 f"{N_RUNS} runs, prefix snapshots.")
    lines.append("| run | M | V | s | s/V |")
    lines.append("|---:|---:|---:|---:|---:|")
    for r in rows:
        lines.append(f"| {r['run']} | {r['M']/1e6:.1f}M | {r['V']} | {r['s']} | {r['s_over_V']} |")
    # trend of s/V with V (per run pooled)
    lV = np.log([r["V"] for r in rows])
    lsv = np.log([max(r["s_over_V"], 1e-9) for r in rows])
    X = np.column_stack([np.ones_like(lV), lV])
    b, *_ = np.linalg.lstsq(X, lsv, rcond=None)
    lines.append("")
    lines.append(f"- GA model: log(s/V) ~ log V slope = {b[1]:+.3f} "
                 f"(real corpora: ~0.00 by the s-law; f12 pooled slope-1 = -0.023)")
    sv = [r["s_over_V"] for r in rows]
    lines.append(f"- GA s/V range across the sweep: {min(sv):.4f} - {max(sv):.4f} "
                 f"(x{max(sv)/max(min(sv),1e-9):.1f}); real corpora hold 0.009-0.013 across V 3.8k-63k")
    lines.append("")
    if abs(b[1]) > 0.15 or max(sv) / max(min(sv), 1e-9) > 2.0:
        lines.append("**Reading: the GA core/non-core model does NOT hold s/V constant — its "
                     "transition width is tied to the fixed core scale, drifting as V grows. "
                     "The measured width law s = 0.012 V is outside this model class as "
                     "parameterized by its own fitted constants.**")
    else:
        lines.append("**Reading: the GA model reproduces an approximately constant s/V over this "
                     "range — the width law does NOT discriminate against it here; weaken the "
                     "manuscript's contrast accordingly.**")
    (OUT / "f14_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
