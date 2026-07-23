"""f15c — the width law at MATCHED depth: within-language concatenation.

Closes the manuscript's final stated limitation ("deepest sampling regimes
unprobed outside English"). Justification for the design is the paper's own
f12 result: corpus composition does not move the width constant (median s/V
flat from single novels to 14-work concatenations), so within-language
concatenation is a valid depth amplifier.

Per language: concatenate the token streams of all cached panel texts plus
the pre-existing deep corpora (French += Les Miserables, Spanish += Don
Quixote, Dutch += Max Havelaar; Russian = W&P alone, already deep). Fit the
canonical erf model (3 replicates). Compare s/V at the resulting depths with
the English reference 0.0118 (tok/V ~ 19).

Outputs: ../outputs/f15c_concat.csv, f15c_summary.md
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
PANEL = REPO / "data" / "gutenberg_panel"
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
from run_f13b_basin_reselect import erf_task  # noqa: E402

UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
N_REPS = 3
EXTRA = {
    "fr": [REPO / "data/zipf_multilang_romance/french_les_miserables/combined_clean.txt"],
    "es": [REPO / "data/zipf_multilang_romance/spanish_don_quixote/combined_clean.txt"],
    "nl": [REPO / "data/zipf_multilang_romance/dutch_max_havelaar/combined_clean.txt"],
    "ru": [REPO / "data/zipf_multilang/russian_war_and_peace/combined_clean.txt"],
}
START_M = ["*** START OF THE PROJECT GUTENBERG", "***START OF THE PROJECT GUTENBERG"]
END_M = ["*** END OF THE PROJECT GUTENBERG", "***END OF THE PROJECT GUTENBERG"]


def strip_bp(text):
    s, e = 0, len(text)
    for m in START_M:
        i = text.find(m)
        if i != -1:
            j = text.find("\n", i)
            s = j + 1 if j != -1 else i
            break
    for m in END_M:
        i = text.find(m)
        if i != -1:
            e = i
            break
    return text[s:e]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    langs = sorted({p.name for p in PANEL.iterdir() if p.is_dir()} | set(EXTRA))
    jobs, meta = [], {}
    for lang in langs:
        toks = []
        for f in sorted((PANEL / lang).glob("pg*.txt")) if (PANEL / lang).exists() else []:
            toks.extend(UNI_RE.findall(strip_bp(f.read_text(encoding="utf-8", errors="ignore")).lower()))
        for f in EXTRA.get(lang, []):
            if f.exists():
                toks.extend(UNI_RE.findall(f.read_text(encoding="utf-8", errors="ignore").lower()))
        if len(toks) < 250_000:
            print(f"skip {lang}: only {len(toks):,} tokens", flush=True)
            continue
        counts = Counter(toks)
        freqs = np.array(sorted(counts.values(), reverse=True), dtype=float)
        name = f"{lang}_concat"
        meta[name] = {"lang": lang, "tokens": len(toks), "V": len(freqs)}
        print(f"{name:12} tokens={len(toks):,} V={len(freqs):,} tok/V={len(toks)/len(freqs):.1f}",
              flush=True)
        for rep in range(N_REPS):
            jobs.append((name, rep, freqs.tolist()))

    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(jobs)}] {res['corpus']:12} rep{res['rep']} "
                  f"s/V={res['s_synth']/res['V_synth']:.4f}", flush=True)

    out_rows = []
    for name, m in meta.items():
        sub = [r["s_synth"] for r in rows if r["corpus"] == name]
        if not sub:
            continue
        s_med = float(np.median(sub))
        out_rows.append({"corpus": name, "lang": m["lang"], "tokens": m["tokens"],
                         "V": m["V"], "tok_per_V": round(m["tokens"] / m["V"], 1),
                         "s": round(s_med, 1), "s_over_V": round(s_med / m["V"], 5)})
    out_rows.sort(key=lambda r: -r["tok_per_V"])
    with open(OUT / "f15c_concat.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)

    lines = ["# f15c — within-language deep concatenation (matched depth)\n"]
    lines.append("Justified by f12 (composition does not move the constant).\n")
    lines.append("| language | tokens | V | tok/V | s/V |")
    lines.append("|---|---:|---:|---:|---:|")
    for r in out_rows:
        lines.append(f"| {r['lang']} | {r['tokens']:,} | {r['V']:,} | {r['tok_per_V']} | "
                     f"**{r['s_over_V']:.4f}** |")
    deep = [r["s_over_V"] for r in out_rows if r["tok_per_V"] >= 15]
    lines.append("")
    if deep:
        lines.append(f"- languages at EN-like depth (tok/V >= 15): n={len(deep)}, "
                     f"median s/V = {float(np.median(deep)):.4f} (EN reference 0.0118)")
    lines.append(f"- all concatenations: median {float(np.median([r['s_over_V'] for r in out_rows])):.4f}")
    (OUT / "f15c_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
