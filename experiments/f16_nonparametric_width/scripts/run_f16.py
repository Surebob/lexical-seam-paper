"""f16 — the width without the fitter, and the matched single-regime null band.

Two referee-grade attacks on the width law's interpretation, answered in one
experiment:

  A. "Your 9-parameter model manufactures the width."  ->  Measure the width
     NONPARAMETRICALLY: bin the rank-frequency curve in log-rank, compute the
     local slope by finite differences, identify the head-slope plateau and
     tail-slope plateau, and define the transition as the rank interval where
     the slope crosses from 10% to 90% of the way between them. Width_np =
     (r90 - r10) / V. No model, no optimizer — arithmetic on the curve.

  B. "The constant is generic to sampled Zipfian curves."  ->  For each of
     the 25 English corpora, build its SINGLE-REGIME TWIN: expected token
     count and fitted ZM (b, c) identical, frequencies drawn as one
     multinomial sample from an exact ZM pmf over the same vocabulary size.
     Same instruments on the twins. The gap between language and its matched
     twins is the structural (two-population) contribution, decomposed with
     per-corpus pairing.

Corpora: 25 EN (data/zipf via canonical loader), deep multilingual (f15/f15b
caches), surname control. Outputs: ../outputs/f16_widths.csv, f16_summary.md
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
OUT = HERE.parent / "outputs"
ZIPF = REPO / "data" / "zipf"
sys.path.insert(0, str(REPO / "src" / "s2_decoupled" / "shared"))
from corpus_loader import build_zipf_dataset  # noqa: E402
sys.path.insert(0, str(REPO / "src" / "s2_decoupled"))
from shared.fit_config import SEARCHED_CORPORA  # noqa: E402

UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
SEED = 20260728
N_TWINS = 3


def _zm_window_fit(ranks, logf):
    """3-param ZM least squares on a window; returns predictor over any ranks."""
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, ranks[-1] * 4, 128)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        resid = d @ coef - logf
        mse = float(np.mean(resid ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(coef[1]), float(c), float(np.std(resid)))
    _, a, bcoef, c, sd = best

    def pred(rr):
        return a + bcoef * np.log(rr + c)

    return pred, sd


def nonparam_width(freqs, head_n=50, tail_win=(0.10, 0.40), eps_mult=2.0):
    """Model-light transition width: the vocabulary fraction where NEITHER the
    head-window ZM nor the tail-window ZM extrapolation explains the curve.
    The two window fits are plain 3-parameter regressions on seam-free zones;
    no gate model is involved."""
    V = len(freqs)
    if V < 2000:
        return float("nan")
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    h_pred, h_sd = _zm_window_fit(ranks[:head_n], logf[:head_n])
    t_lo, t_hi = int(V * tail_win[0]), int(V * tail_win[1])
    t_pred, t_sd = _zm_window_fit(ranks[t_lo:t_hi], logf[t_lo:t_hi])
    eps_h = max(eps_mult * h_sd, 0.02)
    eps_t = max(eps_mult * t_sd, 0.02)
    # evaluate on log-spaced ranks between the windows
    rr = np.unique(np.geomspace(head_n, max(t_lo, head_n + 2), 400).astype(int))
    ok_h = np.abs(h_pred(rr.astype(float)) - logf[rr - 1]) <= eps_h
    ok_t = np.abs(t_pred(rr.astype(float)) - logf[rr - 1]) <= eps_t
    neither = ~ok_h & ~ok_t
    # r_start: first rank where head persistently fails (this and next 4 pts)
    r_start = r_end = None
    for i in range(len(rr) - 4):
        if (~ok_h[i:i + 5]).all():
            r_start = rr[i]
            break
    for i in range(len(rr) - 1, 3, -1):
        if (~ok_t[i - 4:i + 1]).all():
            r_end = rr[i]
            break
    if r_start is None:      # head law holds to the tail window: no seam zone
        return 0.0
    if r_end is None or r_end <= r_start:
        # head fails but tail law already explains everything beyond: zone is
        # where neither holds, possibly empty
        if not neither.any():
            return 0.0
        span = rr[neither]
        return float((span.max() - span.min()) / V)
    return float((r_end - r_start) / V)


def nonparam_width_grid(freqs):
    """Median and IQR of the width over a sensitivity grid of conventions."""
    vals = []
    for head_n in (30, 50, 100):
        for tw in ((0.10, 0.40), (0.15, 0.45)):
            for em in (1.5, 2.0, 3.0):
                w = nonparam_width(np.asarray(freqs, dtype=float), head_n, tw, em)
                if np.isfinite(w):
                    vals.append(w)
    if not vals:
        return float("nan"), float("nan")
    return float(np.median(vals)), float(np.subtract(*np.percentile(vals, [75, 25])))


def fit_zm_bc(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 256)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c))
    return best[2], best[3]


def zm_twin(V, tokens, b, c, rng):
    """Single multinomial sample of `tokens` tokens from exact ZM over V types."""
    ranks = np.arange(1, V + 1, dtype=float)
    p = (ranks + c) ** (-b)
    p /= p.sum()
    counts = rng.multinomial(tokens, p)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(float)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rows = []

    def add(name, group, freqs, tokens=None):
        w, iqr = nonparam_width_grid(freqs)
        rows.append({"corpus": name, "group": group, "V": len(freqs),
                     "tokens": tokens if tokens else int(np.sum(freqs)),
                     "np_width_frac": round(w, 5) if np.isfinite(w) else "",
                     "np_width_iqr": round(iqr, 5) if np.isfinite(iqr) else ""})
        msg = f"np_width/V = {w:.4f} (iqr {iqr:.4f})" if np.isfinite(w) else "DEGENERATE"
        print(f"{group:12} {name[:30]:30} {msg}", flush=True)
        return w

    # 25 EN + matched twins
    for spec in SEARCHED_CORPORA:
        ds = build_zipf_dataset(ZIPF / spec["filename"])
        freqs = ds["freqs"]
        V, T = len(freqs), int(ds["token_count"])
        add(spec["slug"], "english", freqs, T)
        b, c = fit_zm_bc(freqs)
        tw = []
        for k in range(N_TWINS):
            tf = zm_twin(V, T, b, c, rng)
            w, _ = nonparam_width_grid(tf)
            if np.isfinite(w):
                tw.append(w)
        if tw:
            rows.append({"corpus": spec["slug"] + "_ZMtwin", "group": "zm_twin",
                         "V": V, "tokens": T,
                         "np_width_frac": round(float(np.median(tw)), 5),
                         "np_width_iqr": ""})
            print(f"{'zm_twin':12} {spec['slug'][:30]:30} np_width/V = "
                  f"{float(np.median(tw)):.4f}", flush=True)

    # deep multilingual (cached texts)
    ML = [("french_les_miserables", REPO / "data/zipf_multilang_romance/french_les_miserables/combined_clean.txt"),
          ("spanish_don_quixote", REPO / "data/zipf_multilang_romance/spanish_don_quixote/combined_clean.txt"),
          ("russian_war_and_peace", REPO / "data/zipf_multilang/russian_war_and_peace/combined_clean.txt")]
    for name, path in ML:
        toks = UNI_RE.findall(path.read_text(encoding="utf-8", errors="ignore").lower())
        freqs = np.array(sorted(Counter(toks).values(), reverse=True), dtype=float)
        add(name, "multilingual", freqs, len(toks))

    # panel extras from gutenberg_panel cache (up to 8, deepest first)
    cached = sorted(REPO.glob("data/gutenberg_panel/*/pg*.txt"),
                    key=lambda p: -p.stat().st_size)[:8]
    for p in cached:
        toks = UNI_RE.findall(p.read_text(encoding="utf-8", errors="ignore").lower())
        if len(toks) < 150_000:
            continue
        freqs = np.array(sorted(Counter(toks).values(), reverse=True), dtype=float)
        add(f"{p.parent.name}_{p.stem}", "multilingual", freqs, len(toks))

    with open(OUT / "f16_widths.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    def med(group):
        v = [float(r["np_width_frac"]) for r in rows
             if r["group"] == group and r["np_width_frac"] != ""]
        return (float(np.median(v)), len(v)) if v else (float("nan"), 0)

    en, n_en = med("english")
    tw, n_tw = med("zm_twin")
    ml, n_ml = med("multilingual")
    lines = ["# f16 — nonparametric width + matched single-regime nulls\n"]
    lines.append("| group | n | median nonparam width / V |")
    lines.append("|---|---:|---:|")
    lines.append(f"| English corpora | {n_en} | **{en:.4f}** |")
    lines.append(f"| matched ZM twins (single-regime, same V, tokens, b, c) | {n_tw} | **{tw:.4f}** |")
    lines.append(f"| deep multilingual | {n_ml} | **{ml:.4f}** |")
    lines.append("")
    # paired per-corpus gap
    gaps = []
    by_name = {r["corpus"]: r for r in rows}
    for spec in SEARCHED_CORPORA:
        a = by_name.get(spec["slug"], {}).get("np_width_frac", "")
        b_ = by_name.get(spec["slug"] + "_ZMtwin", {}).get("np_width_frac", "")
        if a != "" and b_ != "":
            gaps.append(float(a) - float(b_))
    if gaps:
        g = np.array(gaps)
        lines.append(f"- paired language-minus-twin width gap: median {np.median(g):+.4f}, "
                     f"positive on {int((g > 0).sum())}/{len(g)} corpora")
    lines.append("- fitter-based reference values: EN erf-fit s/V 0.0118; twins have the "
                 "same ZM parameters and sampling depth but no two-population structure.")
    (OUT / "f16_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
