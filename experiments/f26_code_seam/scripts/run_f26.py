"""f26 — does source code have a lexical seam, and whose fingerprint?

Programming languages have their own usage club (keywords, stdlib names,
common identifiers) and long tail (project-specific names). Measure the
seam observables on code-token rank curves:

  corpus A: every .py file under the repo venv's site-packages (a large,
            diverse, real-world Python corpus already on disk);
  corpus B: Greg's JS/TS projects (haul + planetfall source trees, excluding
            node_modules), if present.

Token = identifier-ish word: [A-Za-z_][A-Za-z0-9_]* lowercased, split on
nothing else (raw code words; snake/camel NOT split — stated choice).
Fits: single ZM (b, c), free lambda-ZM (lambda), canonical 9-param erf
(s, s/V) with the f24 weighted log-grid above 200k types.

Question: does code land on the language line (s/V ~ 0.0118 at matched
depth, lambda ~ 20) or carry its own fingerprint (like surnames 0.0266 /
belt 0.0166)? Either answer feeds Paper 2 (atlas) and Paper 3 (code-model
tokenizers inherit whatever boundary code has).

Outputs: ../outputs/f26_results.csv, f26_summary.md
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
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "experiments" / "f24_giant_depth" / "scripts"))
from run_f24 import fit_zm_w, fit_lzm_w, fit_erf_w, SUB_POINTS, V_FULL_MAX  # noqa: E402

WORD = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
DESKTOP = Path("C:/Users/Greg Kara/Desktop/temporary")

CORPORA = {
    "python_sitepackages": {
        "roots": [REPO / ".venv" / "Lib" / "site-packages"],
        "exts": {".py"},
        "exclude": set(),
    },
    "js_projects": {
        "roots": [DESKTOP / "haul", DESKTOP / "planetfall"],
        "exts": {".js", ".ts", ".mjs", ".jsx", ".tsx"},
        "exclude": {"node_modules", "dist", "build", ".git"},
    },
}


def count_corpus(spec):
    counts = Counter()
    files = 0
    for root in spec["roots"]:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if p.suffix.lower() not in spec["exts"]:
                continue
            if any(part in spec["exclude"] for part in p.parts):
                continue
            try:
                text = p.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            files += 1
            for w in WORD.findall(text):
                counts[w.lower()] += 1
    return counts, files


def main():
    rows = []
    for name, spec in CORPORA.items():
        counts, files = count_corpus(spec)
        if not counts:
            print(f"skip {name}: no files found", flush=True)
            continue
        freqs = np.array(sorted(counts.values(), reverse=True), dtype=np.float64)
        T = float(freqs.sum())
        V = len(freqs)
        print(f"{name}: files={files:,} tokens={int(T):,} V={V:,} tok/V={T/V:.1f}", flush=True)
        if V > V_FULL_MAX:
            r_idx = np.unique(np.geomspace(1, V, SUB_POINTS).astype(np.int64))
            ranks = r_idx.astype(np.float64)
            edges = np.empty(len(ranks) + 1)
            edges[1:-1] = (ranks[:-1] + ranks[1:]) / 2
            edges[0], edges[-1] = 0.5, V + 0.5
            wts = np.diff(edges)
            logf = np.log(freqs[r_idx - 1])
        else:
            ranks = np.arange(1, V + 1, dtype=np.float64)
            wts = np.ones(V)
            logf = np.log(freqs)
        sw = np.sqrt(wts / wts.mean())
        zm = fit_zm_w(ranks, logf, sw, V)
        lzm = fit_lzm_w(ranks, logf, sw, V)
        erf = fit_erf_w(ranks, logf, sw, V, seed_salt=hash(name) % 1000, procs=4)
        row = {"corpus": name, "files": files, "tokens": int(T), "V": V,
               "tok_per_V": round(T / V, 1), "zm_b": round(zm["b"], 4),
               "zm_c": round(zm["c"], 2), "zm_rmse": round(zm["rmse"], 4),
               "lam": round(lzm["lam"], 2), "lzm_rmse": round(lzm["rmse"], 4),
               "k": round(erf["k"], 1), "s": round(erf["s"], 1),
               "s_over_V": round(erf["s"] / V, 5), "erf_rmse": round(erf["rmse"], 4)}
        rows.append(row)
        print(f"  b={row['zm_b']} c={row['zm_c']} lam={row['lam']} "
              f"s/V={row['s_over_V']} (erf rmse {row['erf_rmse']})", flush=True)

    if not rows:
        print("no corpora produced results", flush=True)
        return 1
    cols = list(rows[0].keys())
    with open(OUT / "f26_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    lines = ["# f26 — the code seam\n",
             "Reference fingerprints: language s/V 0.0118 (lam ~20-26 by depth), "
             "surnames 0.0266, belt 0.0166, classic Simon 0.0101.\n",
             "| corpus | tokens | V | tok/V | b | c | lambda | s/V |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['corpus']} | {r['tokens']:,} | {r['V']:,} | {r['tok_per_V']} | "
                     f"{r['zm_b']} | {r['zm_c']} | {r['lam']} | {r['s_over_V']} |")
    lines.append("\nCaveats: identifier tokenization (no snake/camel splitting — a raw "
                 "'code word' definition); site-packages mixes library code + vendored "
                 "assets; depth differs from book regime — compare via the depth curves, "
                 "not raw constants.")
    (OUT / "f26_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
