"""f16c — planted-truth calibration of the width instruments.

The gap this closes: no width estimator in the program (including the
canonical erf fitter) has ever been validated against corpora with a KNOWN
transition width. 5a validated residual-signal reproduction; 3e validated
gate-family identification; neither validated width RECOVERY.

Part A — planted seams. Generate two-regime rank curves from the decoupled
model with known (k*, w_gate*, w_tail*), hence known s* = k*.w_tail*, at
realistic head/tail parameters; Poisson-sample token counts at controlled
depth; drop zeros. Grid: V* in {8k, 16k, 32k} x s*/V* in {0.006, 0.012,
0.024} x depth (tokens/type) in {10, 20} x planted gate in {erf, logistic}
(model-mismatch arm), 2 replicates = 72 corpora. Fit the canonical erf model
-> s_hat. Recovery regression s_hat vs s* (slope, R2, bias) is the
calibration of the paper's headline instrument.

Part B — the FWHM instrument (fitter-free v3). The seam is, definitionally,
the structured residual left by the global 3-parameter ZM fit (manuscript
Fig. 1). Measure its width directly: residual vs log-rank, binned + smoothed,
find the dominant interior extremum, take the full-width-at-half-maximum.
Only the field-standard ZM baseline plus a smoothing bandwidth (sensitivity
grid) — no gate model. Run on planted corpora (does FWHM track s*?), on the
25 English corpora, and on matched single-regime ZM twins (amplitude should
collapse to the noise floor: the presence/absence contrast).

Outputs: ../outputs/f16c_planted.csv, f16c_real.csv, f16c_summary.md
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
from run_f13b_basin_reselect import erf_task  # noqa: E402
sys.path.insert(0, str(REPO / "src" / "s2_decoupled"))
from shared.corpus_loader import build_zipf_dataset  # noqa: E402
from shared.fit_config import SEARCHED_CORPORA  # noqa: E402

SEED = 20260729
V_GRID = [8000, 16000, 32000]
SV_GRID = [0.006, 0.012, 0.024]
DEPTH_GRID = [10, 20]
GATES = ["erf", "logistic"]
N_REP = 2


def gate_fn(name, z):
    from scipy.special import erf as sp_erf
    if name == "erf":
        return 0.5 * (1.0 - sp_erf(z))
    return 1.0 / (1.0 + np.exp(1.7 * z))  # logistic matched to erf slope scale


def planted_curve(V, s_frac, gate, depth, rng):
    """Expected freqs from the decoupled model with known seam scale."""
    ranks = np.arange(1, V + 1, dtype=float)
    k = 0.02 * V * math.exp(rng.normal(0.0, 0.15))
    w_tail = s_frac * V / k
    w_gate = 0.6
    b1, c1 = 0.95 + rng.normal(0, 0.05), 8.0
    b2, c2 = 1.55 + rng.normal(0, 0.08), 40.0
    z = (np.log(ranks) - math.log(k)) / w_gate
    sig = gate_fn(gate, z)
    scale = max(1.0, k * w_tail)
    zz = np.clip((ranks - k) / scale, -60, 60)
    rho = 1.0 + scale * np.log1p(np.exp(zz))
    y = sig * (-b1 * np.log(ranks + c1)) + (1 - sig) * (
        -b2 * np.log(rho + c2) - b1 * math.log(k + c1) + b2 * math.log(1 + scale * math.log(2) + c2))
    f = np.exp(y - y.max())
    f *= depth * V / f.sum()
    counts = rng.poisson(f)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(float), k * w_tail


def fit_global_zm(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 256)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        resid = logf - d @ coef
        mse = float(np.mean(resid ** 2))
        if best is None or mse < best[0]:
            best = (mse, resid)
    return best[1]


def fwhm_width(freqs, n_bins=60, smooth_bins=(3, 5)):
    """FWHM of the dominant interior structure of the global-ZM residual.
    Returns (median width fraction over smoothing grid, amplitude)."""
    V = len(freqs)
    if V < 3000:
        return float("nan"), float("nan")
    resid = fit_global_zm(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    hapax_start = int(np.searchsorted(-freqs, -2.5))  # ranks with count < 3 excluded
    lr = np.log(ranks)
    edges = np.linspace(0, lr[-1], n_bins + 1)
    idx = np.clip(np.digitize(lr, edges) - 1, 0, n_bins - 1)
    vals, amps = [], []
    for sb in smooth_bins:
        bx, by = [], []
        for b in range(n_bins):
            m = (idx == b) & (ranks < max(hapax_start, 50))
            if m.sum() >= 3:
                bx.append(lr[m].mean())
                by.append(resid[m].mean())
        if len(by) < 8:
            continue
        by_s = np.convolve(by, np.ones(sb) / sb, mode="same")
        interior = slice(2, len(by_s) - 2)
        seg = by_s[interior]
        bxs = np.asarray(bx)[interior]
        i_ext = int(np.argmax(np.abs(seg)))
        amp = float(abs(seg[i_ext]))
        half = amp / 2.0
        sgn = np.sign(seg[i_ext])
        lo = i_ext
        while lo > 0 and sgn * seg[lo] >= half:
            lo -= 1
        hi = i_ext
        while hi < len(seg) - 1 and sgn * seg[hi] >= half:
            hi += 1
        r_lo, r_hi = math.exp(bxs[lo]), math.exp(bxs[hi])
        vals.append((r_hi - r_lo) / V)
        amps.append(amp)
    if not vals:
        return float("nan"), float("nan")
    return float(np.median(vals)), float(np.median(amps))


def zm_twin(freqs, rng):
    V = len(freqs)
    T = int(freqs.sum())
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 128)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(-coef[1]), float(c))
    _, b, c = best
    p = (ranks + c) ** (-b)
    p /= p.sum()
    counts = rng.multinomial(T, p)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(float)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    # ---- Part A: planted corpora ----
    jobs, meta = [], {}
    for V in V_GRID:
        for sv in SV_GRID:
            for depth in DEPTH_GRID:
                for gate in GATES:
                    for rep in range(N_REP):
                        freqs, s_true = planted_curve(V, sv, gate, depth, rng)
                        name = f"pl_V{V//1000}k_sv{sv}_d{depth}_{gate}_r{rep}"
                        w_np, amp = fwhm_width(freqs)
                        meta[name] = {"V_star": V, "s_frac_star": sv, "depth": depth,
                                      "gate": gate, "rep": rep, "V_real": len(freqs),
                                      "s_star": s_true, "fwhm": w_np, "amp": amp}
                        jobs.append((name, 0, freqs.tolist()))
    print(f"{len(jobs)} planted fits queued", flush=True)
    prow = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            m = meta[res["corpus"]]
            prow.append({**{k: (round(v, 5) if isinstance(v, float) else v) for k, v in m.items()},
                         "s_hat": round(res["s_synth"], 1),
                         "s_hat_frac": round(res["s_synth"] / m["V_real"], 5)})
            print(f"[{i}/{len(jobs)}] {res['corpus'][:34]:34} s*={m['s_star']:.0f} "
                  f"s^={res['s_synth']:.0f}", flush=True)
    with open(OUT / "f16c_planted.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(prow[0].keys()))
        w.writeheader(); w.writerows(prow)

    # ---- Part B: real corpora + twins under FWHM ----
    rrows = []
    for spec in SEARCHED_CORPORA:
        ds = build_zipf_dataset(ZIPF / spec["filename"])
        freqs = ds["freqs"]
        w_np, amp = fwhm_width(freqs)
        rrows.append({"corpus": spec["slug"], "group": "english", "V": len(freqs),
                      "fwhm_frac": round(w_np, 5) if np.isfinite(w_np) else "",
                      "amplitude": round(amp, 4) if np.isfinite(amp) else ""})
        tw_w, tw_a = [], []
        for _ in range(3):
            tf = zm_twin(freqs, rng)
            w2, a2 = fwhm_width(tf)
            if np.isfinite(w2):
                tw_w.append(w2); tw_a.append(a2)
        if tw_w:
            rrows.append({"corpus": spec["slug"] + "_twin", "group": "twin",
                          "V": len(freqs),
                          "fwhm_frac": round(float(np.median(tw_w)), 5),
                          "amplitude": round(float(np.median(tw_a)), 4)})
        print(f"real {spec['slug'][:28]:28} fwhm={w_np:.4f} amp={amp:.3f}", flush=True)
    with open(OUT / "f16c_real.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rrows[0].keys()))
        w.writeheader(); w.writerows(rrows)

    # ---- summary ----
    s_star = np.array([r["s_star"] for r in prow])
    s_hat = np.array([r["s_hat"] for r in prow])
    ok = (s_star > 0) & (s_hat > 0)
    lx, ly = np.log(s_star[ok]), np.log(s_hat[ok])
    X = np.column_stack([np.ones_like(lx), lx])
    b, *_ = np.linalg.lstsq(X, ly, rcond=None)
    r2 = 1 - float(np.sum((ly - X @ b) ** 2)) / float(np.sum((ly - ly.mean()) ** 2))
    ratio = np.exp(np.median(ly - lx))
    fw = np.array([r["fwhm"] for r in prow])
    okf = ok & np.isfinite(fw) & (fw > 0)
    corr_f = float(np.corrcoef(np.log(s_star[okf]), np.log(fw[okf]))[0, 1]) if okf.sum() > 5 else float("nan")

    en_amp = [float(r["amplitude"]) for r in rrows if r["group"] == "english" and r["amplitude"] != ""]
    tw_amp = [float(r["amplitude"]) for r in rrows if r["group"] == "twin" and r["amplitude"] != ""]
    en_fw = [float(r["fwhm_frac"]) for r in rrows if r["group"] == "english" and r["fwhm_frac"] != ""]
    tw_fw = [float(r["fwhm_frac"]) for r in rrows if r["group"] == "twin" and r["fwhm_frac"] != ""]

    lines = ["# f16c — planted-truth calibration + FWHM instrument\n"]
    lines.append("## A. erf-fitter recovery on planted seams (n=%d)" % ok.sum())
    lines.append(f"- log s_hat ~ log s*: slope {b[1]:.3f}, R2 {r2:.3f}; "
                 f"median multiplicative bias s_hat/s* = {ratio:.3f}")
    for gate in GATES:
        m = ok & np.array([r["gate"] == gate for r in prow])
        if m.sum():
            rr = np.exp(np.median(np.log(s_hat[m] / s_star[m])))
            lines.append(f"  - {gate}-planted: median s_hat/s* = {rr:.3f} (n={m.sum()})")
    lines.append(f"- FWHM instrument vs planted s*: corr(log, log) = {corr_f:.3f}")
    lines.append("")
    lines.append("## B. FWHM on real corpora vs single-regime twins")
    lines.append(f"- English: median FWHM/V = {np.median(en_fw):.4f}, "
                 f"median amplitude = {np.median(en_amp):.3f} (n={len(en_fw)})")
    lines.append(f"- twins:   median FWHM/V = {np.median(tw_fw):.4f}, "
                 f"median amplitude = {np.median(tw_amp):.3f} (n={len(tw_fw)})")
    amp_sep = float(np.median(en_amp) / max(np.median(tw_amp), 1e-9))
    lines.append(f"- amplitude separation (language/twin): x{amp_sep:.1f}")
    (OUT / "f16c_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
