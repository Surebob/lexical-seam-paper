"""F14b — is the 1.2% constant specific to two-class structure, or generic?

f14 found the GA core/noncore model reproduces s/V ~ 0.012. Two nulls decide
what that means:

  pure_simon: pc0 = 0 (NO core class; single word class, p_new decaying with
      vocabulary exactly as in GA). If this also gives ~0.012, the constant is
      generic to decaying-innovation preferential attachment, and the
      two-CLASS structure is not what pins it.
  classic_simon: alpha = 0 (constant p_new — the original Simon model, a
      single clean Zipf regime). The erf fitter on a genuinely single-regime
      curve should NOT produce the law (degenerate gate or a different
      constant); if it does, 0.012 is partly a fitting-operator artifact.

Same simulation and fitter as f14; snapshots at 0.4M and 1.6M tokens, 2 runs.

Outputs: ../outputs/f14b_nulls.csv, f14b_summary.md
"""
from __future__ import annotations

import csv
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
from run_f13b_basin_reselect import erf_task  # noqa: E402

SEED = 20260726
SNAPSHOTS = [400_000, 1_600_000]
N_RUNS = 2

VARIANTS = {
    # name: (NC_MAX, ALPHA, PC0, PNEW0)
    "ga_reference": (7900, 0.77, 0.99, 0.10),
    "pure_simon_decay": (7900, 0.77, 0.00, 0.10),   # no core; same decay rule
    "classic_simon": (7900, 0.00, 0.99, 0.05),      # constant p_new (alpha=0)
}


def simulate(nc_max, alpha, pc0, pnew0, m_max, rng):
    token_types = np.empty(m_max, dtype=np.int32)
    counts = {}
    p_new = pnew0
    n_types = 0
    n_core = 0
    n_noncore = 0
    snaps = []
    si = 0
    for m in range(m_max):
        if n_types == 0 or rng.random() < p_new:
            t = n_types
            n_types += 1
            if n_core < nc_max and rng.random() < pc0:
                n_core += 1
            else:
                n_noncore += 1
                if alpha > 0:
                    p_new *= (1.0 - alpha / (n_noncore + nc_max))
            counts[t] = 1
        else:
            t = int(token_types[rng.integers(0, m)])
            counts[t] += 1
        token_types[m] = t
        if si < len(SNAPSHOTS) and (m + 1) == SNAPSHOTS[si]:
            freqs = np.sort(np.fromiter(counts.values(), dtype=np.int64))[::-1]
            snaps.append((m + 1, freqs.astype(float)))
            si += 1
    return snaps


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs, meta = [], {}
    for vname, (nc, al, pc, pn) in VARIANTS.items():
        for run in range(N_RUNS):
            rng = np.random.default_rng(SEED + run)
            for M, freqs in simulate(nc, al, pc, pn, SNAPSHOTS[-1], rng):
                name = f"{vname}_r{run}_M{M//1000}k"
                meta[name] = {"variant": vname, "run": run, "M": M, "V": len(freqs)}
                jobs.append((name, 0, freqs.tolist()))
                print(f"{name:36} V={len(freqs):,}", flush=True)

    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            m = meta[res["corpus"]]
            rows.append({"variant": m["variant"], "run": m["run"], "M": m["M"],
                         "V": res["V_synth"], "s": round(res["s_synth"], 1),
                         "s_over_V": round(res["s_synth"] / res["V_synth"], 5)})
            print(f"[{i}/{len(jobs)}] {res['corpus']:36} s/V={res['s_synth']/res['V_synth']:.4f}",
                  flush=True)

    rows.sort(key=lambda r: (r["variant"], r["run"], r["M"]))
    with open(OUT / "f14b_nulls.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    lines = ["# F14b — nulls: what pins the 1.2%?\n"]
    lines.append("| variant | run | M | V | s | s/V |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for r in rows:
        lines.append(f"| {r['variant']} | {r['run']} | {r['M']/1e6:.1f}M | {r['V']} | "
                     f"{r['s']} | {r['s_over_V']} |")
    lines.append("")
    for vname in VARIANTS:
        sv = [r["s_over_V"] for r in rows if r["variant"] == vname]
        lines.append(f"- {vname}: s/V {min(sv):.4f}-{max(sv):.4f} (median {np.median(sv):.4f})")
    (OUT / "f14b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
