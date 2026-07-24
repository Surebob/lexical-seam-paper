"""f25 — is there a THIRD constant: the club's token-mass share?

The seam's location k is lawless in V; its width s is lawful. Untested: the
fraction of ALL tokens carried by words above the seam centre, M(k). If M(k)
is stable across corpora while k itself wanders, language commits a fixed
share of usage to the club — a mass constant to sit beside the width and
amplitude constants.

k per corpus: the canonical 9-param fits from f18 (kind=real), excluding the
two detectable bound-pinned escapes. Sensitivity: mass also reported at k/2
and 2k. Cheap arithmetic on existing artifacts.

Outputs: ../outputs/f25_mass.csv, f25_summary.md
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "src" / "s2_decoupled"))
from shared.corpus_loader import build_zipf_dataset  # noqa: E402
from shared.fit_config import SEARCHED_CORPORA  # noqa: E402

F18 = REPO / "experiments" / "f18_null_width" / "outputs" / "f18_null.csv"


def main():
    kmap = {}
    for r in csv.DictReader(open(F18, encoding="utf-8")):
        if r["kind"] == "real" and r["wg_pinned"] != "True" and r["wt_pinned"] != "True":
            kmap[r["corpus"]] = float(r["s_hat"]) / float(r["w_tail"])
    rows = []
    for spec in SEARCHED_CORPORA:
        slug = spec["slug"]
        if slug not in kmap:
            print(f"skip {slug} (pinned escape in f18)", flush=True)
            continue
        ds = build_zipf_dataset(REPO / "data" / "zipf" / spec["filename"])
        freqs = np.asarray(ds["freqs"], dtype=np.float64)
        T = float(freqs.sum())
        V = len(freqs)
        k = kmap[slug]
        cum = np.cumsum(freqs)

        def mass(at):
            i = int(min(max(round(at), 1), V))
            return float(cum[i - 1] / T)

        row = {"corpus": slug, "V": V, "tokens": int(T), "k": round(k, 1),
               "k_over_V": round(k / V, 4), "M_k": round(mass(k), 4),
               "M_half_k": round(mass(k / 2), 4), "M_2k": round(mass(2 * k), 4)}
        rows.append(row)
        print(f"{slug[:28]:29} k={k:7.1f} k/V={row['k_over_V']:.4f} M(k)={row['M_k']:.4f}", flush=True)

    cols = list(rows[0].keys())
    with open(OUT / "f25_mass.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    mk = np.array([r["M_k"] for r in rows])
    kv = np.array([r["k_over_V"] for r in rows])
    lines = ["# f25 — club token-mass share at the seam centre\n",
             f"- corpora: {len(rows)} (2 f18-pinned escapes excluded)",
             f"- M(k): median {np.median(mk):.4f}, IQR [{np.percentile(mk,25):.4f}, "
             f"{np.percentile(mk,75):.4f}], range [{mk.min():.4f}, {mk.max():.4f}]",
             f"- rel spread (IQR/median): {(np.percentile(mk,75)-np.percentile(mk,25))/np.median(mk):.2f}",
             f"- compare k/V (the lawless location): median {np.median(kv):.4f}, "
             f"range [{kv.min():.4f}, {kv.max():.4f}], rel spread "
             f"{(np.percentile(kv,75)-np.percentile(kv,25))/np.median(kv):.2f}",
             f"- corr(M_k, log V): {np.corrcoef(mk, np.log([r['V'] for r in rows]))[0,1]:+.3f}"]
    m = float(np.median(mk))
    iqr_rel = (np.percentile(mk, 75) - np.percentile(mk, 25)) / m
    if iqr_rel < 0.15:
        lines.append(f"\n**Reading: M(k) ~ {m:.2f} is tight across corpora — a candidate "
                     "third constant (club mass share). Panel/language extension warranted.**")
    else:
        lines.append("\n**Reading: M(k) varies substantially — the mass share is not a "
                     "constant; k's wandering carries real mass variation with it.**")
    (OUT / "f25_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote outputs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
