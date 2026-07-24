"""f20 — is the lambda amplitude a second fingerprint axis, and is it deeper
than the width?

f19 found lambda* ~ 20.6 universal across 32 language corpora. Three questions
the finding begs:

  A. MECHANISM: do the simulated growth processes generate it? Run the f14b
     variants (GA reference w/ published constants, coreless decaying Simon,
     classic constant-innovation Simon) and fit free lambda on their rank
     curves. If GA/decaying land near ~20 while classic Simon lands elsewhere,
     lambda is mechanism-generated and (s/V, lambda) becomes a 2-D family
     fingerprint.
  B. BREADTH: fit free lambda on the 12-language Gutenberg panel concats
     (languages the f19 transfer never touched).
  C. DEPTH: binomial-thin 4 deep corpora to {1, 1/2, 1/4, 1/8, 1/16} x 2 reps
     and fit lambda at each depth. s/V rides a rising depth curve (f15d);
     if lambda holds flat, the amplitude is the depth-invariant constant.

All fits: the f3/f19 protocol — g(x) = exp(x-1)-x, x = 0.05+0.95 log r/log V,
c-grid + lstsq. Incremental CSV. Outputs: ../outputs/f20_results.csv,
f20_summary.md
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
DATA = REPO / "data" / "zipf"
PANEL = REPO / "data" / "gutenberg_panel"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "experiments" / "f14_ga_break" / "scripts"))
from run_f14b_nulls import simulate, VARIANTS, SNAPSHOTS  # noqa: E402

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
SEED = 20260733
CSV_PATH = OUT / "f20_results.csv"
FIELDS = ["part", "name", "detail", "tokens", "V", "tok_per_V", "lam", "rmse_zm", "rmse_lzm", "impr_pct"]

DEPTH_CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"),
    ("War and Peace", "pg2600.txt"),
    ("Les Miserables", "pg135.txt"),
    ("King James Bible", "pg10.txt"),
]
FRACTIONS = [1.0, 0.5, 0.25, 0.125, 0.0625]


def strip_gutenberg(text):
    start, end = 0, len(text)
    for m in START_MARKERS:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            start = nl + 1 if nl != -1 else i
            break
    for m in END_MARKERS:
        i = text.find(m)
        if i != -1:
            end = i
            break
    return text[start:end]


def gx_of(ranks, V_norm):
    x = 0.05 + 0.95 * np.log(ranks) / math.log(V_norm)
    return np.exp(x - 1.0) - x


def c_grid_for(max_rank):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, max_rank, 2048)])


def fit_pair(freqs):
    """ZM (3p) and free lambda-ZM (4p) on one curve; returns lam + both rmses."""
    freqs = np.asarray(freqs, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    logf = np.log(freqs)
    gx = gx_of(ranks, V)
    best_z, best_l = None, None
    for c in c_grid_for(float(V)):
        lr = np.log(ranks + c)
        d2 = np.column_stack([np.ones_like(logf), lr])
        coef, *_ = np.linalg.lstsq(d2, logf, rcond=None)
        mse = float(np.mean((d2 @ coef - logf) ** 2))
        if best_z is None or mse < best_z:
            best_z = mse
        d3 = np.column_stack([np.ones_like(logf), lr, gx])
        coef3, *_ = np.linalg.lstsq(d3, logf, rcond=None)
        mse3 = float(np.mean((d3 @ coef3 - logf) ** 2))
        if best_l is None or mse3 < best_l[0]:
            best_l = (mse3, float(coef3[2]))
    rz, rl = math.sqrt(best_z), math.sqrt(best_l[0])
    return {"lam": best_l[1], "rmse_zm": rz, "rmse_lzm": rl,
            "impr_pct": 100 * (rz - rl) / rz}


def emit(writer, fh, done, part, name, detail, freqs):
    key = (part, name, detail)
    if key in done:
        return
    freqs = np.asarray(freqs, dtype=np.float64)
    T = float(freqs.sum())
    V = len(freqs)
    r = fit_pair(freqs)
    row = {"part": part, "name": name, "detail": detail, "tokens": int(T), "V": V,
           "tok_per_V": round(T / V, 1), "lam": round(r["lam"], 2),
           "rmse_zm": round(r["rmse_zm"], 4), "rmse_lzm": round(r["rmse_lzm"], 4),
           "impr_pct": round(r["impr_pct"], 2)}
    writer.writerow(row)
    fh.flush()
    print(f"{part} {name[:26]:27} {detail:12} lam={row['lam']:7.2f} "
          f"impr={row['impr_pct']:+5.1f}% tok/V={row['tok_per_V']}", flush=True)


def main():
    rng = np.random.default_rng(SEED)
    done = set()
    if CSV_PATH.exists():
        with open(CSV_PATH, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                done.add((r["part"], r["name"], r["detail"]))
        print(f"resume: {len(done)} rows present", flush=True)
    new_file = not CSV_PATH.exists()
    fh = open(CSV_PATH, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new_file:
        w.writeheader(); fh.flush()

    print("== part A: simulated mechanisms ==", flush=True)
    for vname, (nc, al, pc, pn) in VARIANTS.items():
        for run in range(2):
            srng = np.random.default_rng(SEED + 100 * run + hash(vname) % 50)
            for M, freqs in simulate(nc, al, pc, pn, SNAPSHOTS[-1], srng):
                emit(w, fh, done, "sim", vname, f"run{run}_M{M//1000}k", freqs)

    print("== part B: 12-language panel concats ==", flush=True)
    if PANEL.exists():
        for lang_dir in sorted(PANEL.iterdir()):
            if not lang_dir.is_dir():
                continue
            toks = []
            for pg in sorted(lang_dir.glob("pg*.txt")):
                text = strip_gutenberg(pg.read_text(encoding="utf-8", errors="ignore"))
                toks.extend(UNI_RE.findall(text.lower()))
            if len(toks) < 50_000:
                print(f"skip {lang_dir.name}: only {len(toks)} tokens", flush=True)
                continue
            freqs = np.array(sorted(Counter(toks).values(), reverse=True), dtype=np.float64)
            emit(w, fh, done, "panel", lang_dir.name, "concat", freqs)

    print("== part C: lambda vs depth (binomial thinning) ==", flush=True)
    for name, fname in DEPTH_CORPORA:
        text = strip_gutenberg((DATA / fname).read_text(encoding="utf-8", errors="ignore"))
        counts = Counter(TOKEN_RE.findall(text.lower()))
        words = list(counts.keys())
        n_full = np.array([counts[wd] for wd in words], dtype=np.int64)
        for frac in FRACTIONS:
            reps = 1 if frac == 1.0 else 2
            for rep in range(reps):
                if frac == 1.0:
                    sub = n_full
                else:
                    sub = rng.binomial(n_full, frac)
                freqs = np.sort(sub[sub > 0])[::-1].astype(np.float64)
                emit(w, fh, done, "depth", name, f"f{frac}_r{rep}", freqs)
    fh.close()

    rows = list(csv.DictReader(open(CSV_PATH, encoding="utf-8")))
    for r in rows:
        r["lam"] = float(r["lam"]); r["tok_per_V"] = float(r["tok_per_V"])
        r["impr_pct"] = float(r["impr_pct"])
    lines = ["# f20 — lambda as fingerprint axis / depth invariance\n",
             "Reference: language free-fit lambda median 20.6 (EN 25, f3), 21.9 (7 languages, f19).\n"]
    sims = [r for r in rows if r["part"] == "sim"]
    if sims:
        lines.append("## A. Simulated mechanisms (free lambda)")
        for vname in VARIANTS:
            ls = [r["lam"] for r in sims if r["name"] == vname]
            if ls:
                lines.append(f"- {vname}: median lam {np.median(ls):.1f} (n={len(ls)}, "
                             f"range [{min(ls):.1f}, {max(ls):.1f}])")
        lines.append("")
    panel = [r for r in rows if r["part"] == "panel"]
    if panel:
        ls = [r["lam"] for r in panel]
        lines.append("## B. 12-language panel concats")
        lines.append(f"- median lam {np.median(ls):.1f}, IQR [{np.percentile(ls,25):.1f}, "
                     f"{np.percentile(ls,75):.1f}], range [{min(ls):.1f}, {max(ls):.1f}]")
        lines.append(f"- corr(lam, tok/V) across panel: "
                     f"{np.corrcoef([r['tok_per_V'] for r in panel], ls)[0,1]:+.3f}")
        lines.append(f"- lambda-ZM improves fit on {sum(1 for r in panel if r['impr_pct'] > 0)}/{len(panel)}")
        lines.append("")
    depth = [r for r in rows if r["part"] == "depth"]
    if depth:
        lines.append("## C. Lambda vs depth (within-corpus thinning)")
        for name, _ in DEPTH_CORPORA:
            sub = sorted([r for r in depth if r["name"] == name], key=lambda r: -r["tok_per_V"])
            if sub:
                lines.append(f"- {name}: " + "  ".join(
                    f"tok/V {r['tok_per_V']:.0f}: lam {r['lam']:.1f}" for r in sub))
        alld = [(r["tok_per_V"], r["lam"]) for r in depth]
        if len(alld) > 3:
            tv, lm = zip(*alld)
            lines.append(f"- pooled corr(lam, tok/V): {np.corrcoef(tv, lm)[0,1]:+.3f}")
    (OUT / "f20_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
