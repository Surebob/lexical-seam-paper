"""Reconstruct f15b panel results from the run log (the orphaned process died
one fit short of writing outputs; all 57 completed fits and per-corpus
metadata are in the log). The one missing replicate (pt_pg31552 rep0) leaves
that corpus with a single-rep median, noted."""
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent.parent / "outputs"
log = (OUT / "run_log_f15b.txt").read_text(encoding="utf-8", errors="ignore")

meta = {}
for m in re.finditer(r"^\s{2}(\w+_pg\d+)\s+tokens=([\d,]+) V=([\d,]+) b=([\d.]+)", log, re.M):
    name = m.group(1)
    meta[name] = {"tokens": int(m.group(2).replace(",", "")),
                  "V": int(m.group(3).replace(",", "")),
                  "zm_b": float(m.group(4))}

fits = defaultdict(list)
for m in re.finditer(r"\[\d+/58\]\s+(\w+_pg\d+)\s+rep\d\s+s/V=([\d.]+)", log):
    fits[m.group(1)].append(float(m.group(2)))

rows = []
for name, svs in fits.items():
    mm = meta.get(name)
    if not mm:
        continue
    sv = float(np.median(svs))
    rows.append({"corpus": name, "lang": name.split("_")[0], "tokens": mm["tokens"],
                 "V": mm["V"], "zm_b": mm["zm_b"], "s": round(sv * mm["V"], 1),
                 "s_over_V": round(sv, 5), "n_reps": len(svs)})
rows.sort(key=lambda r: (r["lang"], -r["tokens"]))
with open(OUT / "f15b_panel.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

lines = ["# f15b — cross-language width-law panel (full depth; reconstructed from run log)\n"]
lines.append("| corpus | lang | tokens | V | tok/V | s/V | reps |")
lines.append("|---|---|---:|---:|---:|---:|---:|")
for r in rows:
    lines.append(f"| {r['corpus']} | {r['lang']} | {r['tokens']:,} | {r['V']:,} | "
                 f"{r['tokens']/r['V']:.1f} | **{r['s_over_V']:.4f}** | {r['n_reps']} |")
lines.append("")
for lang in sorted({r["lang"] for r in rows}):
    sub = [r for r in rows if r["lang"] == lang]
    sv = [r["s_over_V"] for r in sub]
    dep = [r["tokens"] / r["V"] for r in sub]
    lines.append(f"- {lang}: n={len(sv)}, median s/V {float(np.median(sv)):.4f}, "
                 f"median tok/V {float(np.median(dep)):.1f}")
sv_all = np.array([r["s_over_V"] for r in rows])
dep_all = np.array([r["tokens"] / r["V"] for r in rows])
corr = float(np.corrcoef(np.log(dep_all), sv_all)[0, 1])
deep = sv_all[dep_all >= 12]
shallow = sv_all[dep_all < 12]
lines.append("")
lines.append(f"- panel: {len(rows)} corpora, {len(set(r['lang'] for r in rows))} languages; "
             f"median s/V {float(np.median(sv_all)):.4f}")
lines.append(f"- corr(log tokens-per-type, s/V) = {corr:+.3f}")
lines.append(f"- deep corpora (tok/V >= 12): n={len(deep)}, median s/V "
             f"{float(np.median(deep)):.4f}  |  shallow (tok/V < 12): n={len(shallow)}, "
             f"median {float(np.median(shallow)):.4f}")
lines.append("- EN reference 0.0118 (median tok/V ~ 19). One corpus (pt_pg31552) has a "
             "single replicate (process died before its second).")
(OUT / "f15b_summary.md").write_text("\n".join(lines), encoding="utf-8")
print(f"reconstructed {len(rows)} corpora")
