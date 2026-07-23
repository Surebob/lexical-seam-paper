"""f15d — the approach law: does s/V(depth) collapse onto one universal curve?

f15b found corr(tokens-per-type, s/V) = +0.94 across languages; f15c pushes
languages to matched depth. The remaining question is the SHAPE: slice large
corpora across languages to multiple sampling depths (binomial thinning of
the count vector, validated depth-equivalent to prefix slicing in f6),
measure s/V at each, and test whether all languages ride ONE curve in
(tokens-per-type, s/V) space.

Corpora: three deep English + French/Spanish/Russian deep sets. Depth
fractions {0.1, 0.2, 0.4, 0.7, 1.0}, 2 thinning replicates below 1.0.
Collapse metric: within-bin cross-corpus spread vs total spread after binning
by realized tokens-per-type.

Outputs: ../outputs/f15d_approach.csv, f15d_summary.md
"""
from __future__ import annotations

import csv
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
from run_f13b_basin_reselect import erf_task, counts_for, ZIPF  # noqa: E402

UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
SEED = 20260801
FRACS = [0.1, 0.2, 0.4, 0.7, 1.0]
N_REP = 2

EN = [("en_shakespeare", "pg100.txt"), ("en_war_and_peace", "pg2600.txt"),
      ("en_les_miserables", "pg135.txt")]
ML = [("fr_deep", [REPO / "data/zipf_multilang_romance/french_les_miserables/combined_clean.txt"]
       + sorted((REPO / "data/gutenberg_panel/fr").glob("pg*.txt"))),
      ("es_deep", [REPO / "data/zipf_multilang_romance/spanish_don_quixote/combined_clean.txt"]
       + sorted((REPO / "data/gutenberg_panel/es").glob("pg*.txt"))),
      ("ru_deep", [REPO / "data/zipf_multilang/russian_war_and_peace/combined_clean.txt"])]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    base = {}
    for name, fn in EN:
        base[name] = counts_for(ZIPF / fn).astype(np.int64)
    for name, files in ML:
        toks = []
        for f in files:
            toks.extend(UNI_RE.findall(f.read_text(encoding="utf-8", errors="ignore").lower()))
        base[name] = np.array(sorted(Counter(toks).values(), reverse=True), dtype=np.int64)

    jobs, meta = [], {}
    for name, counts in base.items():
        T = int(counts.sum())
        for frac in FRACS:
            reps = 1 if frac == 1.0 else N_REP
            for rep in range(reps):
                if frac == 1.0:
                    thin = counts.astype(float)
                else:
                    thin = rng.binomial(counts, frac)
                    thin = np.sort(thin[thin > 0])[::-1].astype(float)
                key = f"{name}|{frac}|{rep}"
                meta[key] = {"corpus": name, "frac": frac, "rep": rep,
                             "tokens": int(thin.sum()), "V": len(thin)}
                jobs.append((key, 0, thin.tolist()))
        print(f"{name:18} T={T:,} V={len(counts):,}", flush=True)

    print(f"{len(jobs)} fits queued", flush=True)
    rows = []
    with Pool(processes=6) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            m = meta[res["corpus"]]
            rows.append({**m, "s": round(res["s_synth"], 1),
                         "s_over_V": round(res["s_synth"] / m["V"], 5),
                         "tok_per_V": round(m["tokens"] / m["V"], 2)})
            print(f"[{i}/{len(jobs)}] {m['corpus']:18} frac={m['frac']} "
                  f"tok/V={m['tokens']/m['V']:.1f} s/V={res['s_synth']/m['V']:.4f}", flush=True)

    rows.sort(key=lambda r: (r["corpus"], r["frac"], r["rep"]))
    with open(OUT / "f15d_approach.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # collapse metric: bin by log tok/V; compare within-bin cross-corpus spread
    # of s/V against total spread
    d = np.array([r["tok_per_V"] for r in rows])
    sv = np.array([r["s_over_V"] for r in rows])
    corp = np.array([r["corpus"] for r in rows])
    bins = np.quantile(np.log(d), np.linspace(0, 1, 6))
    b = np.clip(np.digitize(np.log(d), bins) - 1, 0, 4)
    within = []
    for k in range(5):
        m = b == k
        if m.sum() >= 3 and len(set(corp[m])) >= 2:
            within.append(float(np.std(sv[m])))
    collapse = float(np.mean(within)) / float(np.std(sv)) if within else float("nan")

    lines = ["# f15d — the approach law: s/V vs sampling depth across languages\n"]
    lines.append("| corpus | frac | tokens | V | tok/V | s/V |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for r in rows:
        lines.append(f"| {r['corpus']} | {r['frac']} | {r['tokens']:,} | {r['V']:,} | "
                     f"{r['tok_per_V']} | **{r['s_over_V']:.4f}** |")
    lines.append("")
    lines.append(f"- collapse ratio (mean within-depth-bin cross-corpus std / total std): "
                 f"{collapse:.3f}  (<< 1 means the curves collapse onto one function)")
    lines.append(f"- corr(log tok/V, s/V) pooled: "
                 f"{float(np.corrcoef(np.log(d), sv)[0,1]):+.3f}")
    (OUT / "f15d_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
