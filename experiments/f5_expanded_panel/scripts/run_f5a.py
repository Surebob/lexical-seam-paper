"""F5a — light pipeline over the expanded panel: ZM, lambda-ZM, unit step-2 winners,
lexical stats. Registers: Brown (balanced 1961), Cornell movie dialogs (speech-like),
WikiText-103 1M-token slice (modern encyclopedic). Languages: IT/FI/PL/PT/DE/SV
(150k-token slices). Control: US Census 2010 surnames (frequency table direct).

Saves processed frequency vectors to data/processed_panel/*.npy for f5b (gate fits).
Outputs: ../outputs/f5a_panel.csv, f5a_summary.md
"""
from __future__ import annotations

import csv
import io
import math
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MODERN = REPO / "data" / "modern"
LANG = REPO / "data" / "langext"
PROC = REPO / "data" / "processed_panel"
PROC.mkdir(parents=True, exist_ok=True)
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

EN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
SLICE = 150_000
WIKI_TOKENS = 1_000_000

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


def freqs_from_tokens(tokens):
    c = Counter(tokens)
    ranked = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))
    return np.array([f for _, f in ranked], dtype=np.float64)


def load_brown():
    toks = []
    with zipfile.ZipFile(MODERN / "brown.zip") as z:
        for nm in sorted(z.namelist()):
            base = nm.rsplit("/", 1)[-1]
            if not re.fullmatch(r"c[a-r]\d\d", base):
                continue
            text = z.read(nm).decode("utf-8", errors="ignore")
            for wt in text.split():
                w = wt.rsplit("/", 1)[0].lower()
                if EN_RE.fullmatch(w):
                    toks.append(w)
    return freqs_from_tokens(toks), len(toks)


def load_cornell():
    toks = []
    with zipfile.ZipFile(MODERN / "cornell_movie_dialogs.zip") as z:
        nm = next(n for n in z.namelist() if n.endswith("movie_lines.txt"))
        text = z.read(nm).decode("iso-8859-1", errors="ignore")
    for line in text.splitlines():
        parts = line.split(" +++$+++ ")
        if len(parts) >= 5:
            toks.extend(EN_RE.findall(parts[-1].lower()))
    return freqs_from_tokens(toks), len(toks)


def load_wikitext():
    toks = []
    with zipfile.ZipFile(MODERN / "wikitext-103-raw-v1.zip") as z:
        nm = next(n for n in z.namelist() if n.endswith("wiki.train.raw"))
        with z.open(nm) as f:
            buf = io.TextIOWrapper(f, encoding="utf-8", errors="ignore")
            while len(toks) < WIKI_TOKENS:
                chunk = buf.read(1 << 20)
                if not chunk:
                    break
                toks.extend(EN_RE.findall(chunk.lower()))
    toks = toks[:WIKI_TOKENS]
    return freqs_from_tokens(toks), len(toks)


def load_census():
    with zipfile.ZipFile(MODERN / "census_surnames.zip") as z:
        nm = next(n for n in z.namelist() if n.lower().endswith(".csv"))
        text = z.read(nm).decode("utf-8", errors="ignore")
    counts = []
    rdr = csv.DictReader(io.StringIO(text))
    ckey = next(k for k in rdr.fieldnames if k.lower() == "count")
    for row in rdr:
        try:
            counts.append(float(row[ckey]))
        except Exception:
            continue
    freqs = np.array(sorted(counts, reverse=True), dtype=np.float64)
    freqs = freqs[freqs > 0]
    return freqs, int(freqs.sum())


def load_gutenberg_lang(path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    st, e = 0, len(text)
    for m in SM:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            st = nl + 1 if nl != -1 else i
            break
    for m in EMK:
        i = text.find(m)
        if i != -1:
            e = i
            break
    toks = UNI_RE.findall(text[st:e].lower())
    if len(toks) > SLICE:
        toks = toks[:SLICE]
    return freqs_from_tokens(toks), len(toks)


def c_grid_for(mr):
    return np.concatenate([np.array([0.0]), np.geomspace(1e-6, mr, 1536)])


def fit_zm(ranks, logf):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c), d @ coef)
    return {"a": best[1], "b": best[2], "c": best[3], "rmse": math.sqrt(best[0]), "pred": best[4]}


def fit_lzm(ranks, logf, gx):
    best = None
    for c in c_grid_for(ranks[-1]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c), gx])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[2]), float(c))
    return {"rmse": math.sqrt(best[0]), "lam": best[1], "c": best[2]}


def analyze(name, freqs, tokens):
    V = len(freqs)
    np.save(PROC / f"{name}.npy", freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    logf = np.log(freqs)
    logr = np.log(ranks)
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
    zm = fit_zm(ranks, logf)
    resid = logf - zm["pred"]
    unit = {k: math.sqrt(float(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
    lz_exp = fit_lzm(ranks, logf, GENS["exp"](x))
    lz_is = fit_lzm(ranks, logf, GENS["is"](x))
    lz = min(lz_exp["rmse"], lz_is["rmse"])
    bic_zm = 3 * math.log(V) + V * math.log(zm["rmse"] ** 2)
    bic_lzm = 4 * math.log(V) + V * math.log(lz**2)
    hapax = float(np.sum(freqs == 1.0)) / V
    return {
        "corpus": name, "V": V, "tokens": tokens,
        "zm_c": zm["c"], "zm_b": zm["b"], "zm_rmse": zm["rmse"],
        "unit_winner": min(unit, key=unit.get),
        "step2_helps": unit[min(unit, key=unit.get)] < zm["rmse"],
        "lzm_rmse": lz, "lzm_impr_pct": 100 * (zm["rmse"] - lz) / zm["rmse"],
        "lzm_lam_exp": lz_exp["lam"], "lzm_bic_wins": bic_lzm < bic_zm,
        "hapax_ratio": hapax, "ttr": V / max(tokens, 1),
    }


def main():
    panel = []
    for name, loader in [("brown", load_brown), ("cornell_dialogs", load_cornell),
                         ("wikitext_1M", load_wikitext), ("census_surnames", load_census)]:
        try:
            freqs, tokens = loader()
            panel.append(analyze(name, freqs, tokens))
            r = panel[-1]
            print(f"{name:18} V={r['V']:7} tok={r['tokens']:9} c={r['zm_c']:9.2f} zm={r['zm_rmse']:.4f} "
                  f"lzm_impr={r['lzm_impr_pct']:6.2f}% bic_win={r['lzm_bic_wins']} w={r['unit_winner']}", flush=True)
        except Exception as ex:
            print(f"FAIL {name}: {ex!r}", flush=True)
    for path in sorted(LANG.glob("*.txt")):
        name = f"lang_{path.stem}"
        try:
            freqs, tokens = load_gutenberg_lang(path)
            panel.append(analyze(name, freqs, tokens))
            r = panel[-1]
            print(f"{name:18} V={r['V']:7} tok={r['tokens']:9} c={r['zm_c']:9.2f} zm={r['zm_rmse']:.4f} "
                  f"lzm_impr={r['lzm_impr_pct']:6.2f}% bic_win={r['lzm_bic_wins']} w={r['unit_winner']}", flush=True)
        except Exception as ex:
            print(f"FAIL {name}: {ex!r}", flush=True)

    with open(OUT / "f5a_panel.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(panel[0].keys()))
        w.writeheader()
        w.writerows(panel)

    lines = ["# F5a expanded panel — summary\n"]
    lang_rows = [r for r in panel if r["corpus"].startswith("lang_")]
    reg_rows = [r for r in panel if r["corpus"] in ("brown", "cornell_dialogs", "wikitext_1M")]
    ctrl = [r for r in panel if r["corpus"] == "census_surnames"]
    lines.append(f"- lambda-ZM BIC wins: registers {sum(1 for r in reg_rows if r['lzm_bic_wins'])}/{len(reg_rows)}, "
                 f"languages {sum(1 for r in lang_rows if r['lzm_bic_wins'])}/{len(lang_rows)}, "
                 f"surname control {sum(1 for r in ctrl if r['lzm_bic_wins'])}/{len(ctrl)}")
    lines.append(f"- median lambda-ZM improvement: registers "
                 f"{float(np.median([r['lzm_impr_pct'] for r in reg_rows])):.2f}%, languages "
                 f"{float(np.median([r['lzm_impr_pct'] for r in lang_rows])):.2f}%")
    for r in panel:
        lines.append(f"- {r['corpus']}: V={r['V']} c={r['zm_c']:.2f} impr={r['lzm_impr_pct']:.2f}% "
                     f"winner={r['unit_winner']} hapax={r['hapax_ratio']:.3f}")
    (OUT / "f5a_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
