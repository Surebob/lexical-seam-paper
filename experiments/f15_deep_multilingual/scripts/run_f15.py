"""f15 — the s-law at FULL DEPTH in four new languages.

The manuscript's §4.5 concedes that cross-language s-law values were measured
only at shallow interior depths (f5c: 0.0068-0.0096 for IT/FI/SV/DE/PL/PT at
~100-200k tokens). This experiment fits the width law at full depth on four
corpora already in the repository — none of which contributed an s-law value
before:

  French   — Les Misérables (original), ~0.5M words
  Spanish  — Don Quixote (original), ~0.4M words
  Russian  — War and Peace (original), ~0.4M words (Cyrillic)
  Dutch    — Max Havelaar, ~0.13M words

Tokenization: f8's Unicode word regex ([^\\W\\d_]+ on lowercased text), the
same convention behind the manuscript's matched-size panel. Fitter: the
canonical f5b/f12 9-parameter erf-gate (imported), 3 fit replicates per
corpus (fit noise only; deterministic data). ZM baseline (b, c) for context.

Outputs: ../outputs/f15_deep_fits.csv, f15_summary.md
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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
from run_f13b_basin_reselect import erf_task  # noqa: E402  (canonical fitter)

UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)

CORPORA = [
    ("french_les_miserables", REPO / "data/zipf_multilang_romance/french_les_miserables/combined_clean.txt"),
    ("spanish_don_quixote", REPO / "data/zipf_multilang_romance/spanish_don_quixote/combined_clean.txt"),
    ("russian_war_and_peace", REPO / "data/zipf_multilang/russian_war_and_peace/combined_clean.txt"),
    ("dutch_max_havelaar", REPO / "data/zipf_multilang_romance/dutch_max_havelaar/combined_clean.txt"),
]
N_REPS = 3


def fit_zm(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 512)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(-coef[1]), float(c))
    return best[1], best[2]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs, meta = [], {}
    for name, path in CORPORA:
        text = path.read_text(encoding="utf-8", errors="ignore")
        toks = UNI_RE.findall(text.lower())
        counts = Counter(toks)
        freqs = np.array(sorted(counts.values(), reverse=True), dtype=float)
        b, c = fit_zm(freqs)
        meta[name] = {"tokens": len(toks), "V": len(freqs), "zm_b": b, "zm_c": c}
        print(f"{name:24} tokens={len(toks):,} V={len(freqs):,} b={b:.3f} c={c:.1f}", flush=True)
        for rep in range(N_REPS):
            jobs.append((name, rep, freqs.tolist()))

    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(jobs)}] {res['corpus']:24} rep{res['rep']} "
                  f"s={res['s_synth']:8.1f} s/V={res['s_synth']/res['V_synth']:.4f}", flush=True)

    out_rows = []
    for name, path in CORPORA:
        m = meta[name]
        sub = [r["s_synth"] for r in rows if r["corpus"] == name]
        s_med = float(np.median(sub))
        out_rows.append({"corpus": name, "tokens": m["tokens"], "V": m["V"],
                         "zm_b": round(m["zm_b"], 3), "zm_c": round(m["zm_c"], 1),
                         "s": round(s_med, 1), "s_over_V": round(s_med / m["V"], 5),
                         "s_reps": ";".join(f"{v:.1f}" for v in sub)})
    with open(OUT / "f15_deep_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)

    lines = ["# f15 — the s-law at full depth in four new languages\n"]
    lines.append("| corpus | tokens | V | ZM b | ZM c | s | s/V |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for r in out_rows:
        lines.append(f"| {r['corpus']} | {r['tokens']:,} | {r['V']:,} | {r['zm_b']} | "
                     f"{r['zm_c']} | {r['s']} | **{r['s_over_V']:.4f}** |")
    sv = [r["s_over_V"] for r in out_rows]
    lines.append("")
    lines.append(f"- median s/V across the four deep non-English corpora: {float(np.median(sv)):.4f}")
    lines.append("- English reference: 0.0118 (25 corpora); registers 0.0120-0.0122; "
                 "shallow non-English interiors (f5c): 0.0068-0.0096; surnames 0.0266.")
    (OUT / "f15_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
