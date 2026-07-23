"""f15b — the definitive cross-language width-law panel.

Supersedes f15's four-corpus check with a properly powered panel: for each of
~10 languages, fetch the most popular long texts from Project Gutenberg (via
the gutendex catalog API), keep up to 3 texts with >= 150k word tokens, and
fit the canonical erf-gate model at full depth. Non-Latin scripts with
whitespace word boundaries (Cyrillic, Greek) are in scope via the Unicode
tokenizer; Chinese/Arabic segmentation issues remain out of scope as in the
manuscript.

Every download is cached under data/gutenberg_panel/<lang>/ with a
sources.csv (id, title, author, url) for provenance. Deterministic fitter
(f5b/f12/f13), 2 replicates per corpus.

Outputs: ../outputs/f15b_panel.csv, f15b_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OUT = HERE.parent / "outputs"
CACHE = REPO / "data" / "gutenberg_panel"
sys.path.insert(0, str(REPO / "experiments" / "f13_mixture_width_law" / "scripts"))
from run_f13b_basin_reselect import erf_task  # noqa: E402

UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
LANGS = ["fr", "de", "es", "it", "pt", "nl", "fi", "sv", "da", "hu", "pl", "el"]
MIN_TOKENS = 150_000
MAX_PER_LANG = 3
N_REPS = 2
START_M = ["*** START OF THE PROJECT GUTENBERG", "***START OF THE PROJECT GUTENBERG"]
END_M = ["*** END OF THE PROJECT GUTENBERG", "***END OF THE PROJECT GUTENBERG"]


def strip_boilerplate(text):
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


def plain_url(formats):
    for k, v in formats.items():
        if k.startswith("text/plain") and not v.endswith(".zip"):
            return v
    return None


def fetch_language(lang):
    """Return list of (corpus_name, tokens_list) for up to MAX_PER_LANG texts."""
    lang_dir = CACHE / lang
    lang_dir.mkdir(parents=True, exist_ok=True)
    src_path = lang_dir / "sources.csv"
    sources = []
    kept = []

    # reuse cache first
    for f in sorted(lang_dir.glob("pg*.txt")):
        text = strip_boilerplate(f.read_text(encoding="utf-8", errors="ignore"))
        toks = UNI_RE.findall(text.lower())
        if len(toks) >= MIN_TOKENS:
            kept.append((f"{lang}_{f.stem}", toks))
        if len(kept) >= MAX_PER_LANG:
            return kept

    page_url = f"https://gutendex.com/books?languages={lang}&sort=popular"
    seen_pages = 0
    while page_url and len(kept) < MAX_PER_LANG and seen_pages < 6:
        try:
            r = requests.get(page_url, timeout=60)
            r.raise_for_status()
            data = r.json()
        except Exception as ex:
            print(f"  [{lang}] catalog fetch failed: {ex}", flush=True)
            break
        seen_pages += 1
        page_url = data.get("next")
        for book in data.get("results", []):
            if len(kept) >= MAX_PER_LANG:
                break
            bid = book["id"]
            fpath = lang_dir / f"pg{bid}.txt"
            if fpath.exists():
                continue
            url = plain_url(book.get("formats", {}))
            if not url:
                continue
            try:
                t = requests.get(url, timeout=120)
                t.raise_for_status()
                raw = t.text
            except Exception:
                continue
            time.sleep(1.0)
            text = strip_boilerplate(raw)
            toks = UNI_RE.findall(text.lower())
            if len(toks) < MIN_TOKENS:
                continue
            fpath.write_text(raw, encoding="utf-8")
            sources.append({"id": bid, "title": book.get("title", "")[:80],
                            "authors": ";".join(a.get("name", "") for a in book.get("authors", []))[:80],
                            "url": url, "tokens": len(toks)})
            kept.append((f"{lang}_pg{bid}", toks))
            safe_title = book.get("title", "")[:40].encode("ascii", "replace").decode()
            print(f"  [{lang}] kept pg{bid} '{safe_title}' "
                  f"({len(toks):,} tokens)", flush=True)
    if sources:
        exists = src_path.exists()
        with open(src_path, "a", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["id", "title", "authors", "url", "tokens"])
            if not exists:
                w.writeheader()
            w.writerows(sources)
    return kept


def fit_zm(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 384)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(-coef[1]), float(c))
    return best[1], best[2]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs, meta = [], {}
    for lang in LANGS:
        print(f"language {lang}...", flush=True)
        for name, toks in fetch_language(lang):
            counts = Counter(toks)
            freqs = np.array(sorted(counts.values(), reverse=True), dtype=float)
            if len(freqs) < 3000:
                continue
            b, c = fit_zm(freqs)
            meta[name] = {"lang": lang, "tokens": len(toks), "V": len(freqs),
                          "zm_b": b, "zm_c": c}
            print(f"  {name:22} tokens={len(toks):,} V={len(freqs):,} b={b:.3f}", flush=True)
            for rep in range(N_REPS):
                jobs.append((name, rep, freqs.tolist()))

    print(f"\n{len(jobs)} fits queued over {len(meta)} corpora", flush=True)
    rows = []
    with Pool(processes=8) as pool:
        for i, res in enumerate(pool.imap_unordered(erf_task, jobs, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(jobs)}] {res['corpus']:22} rep{res['rep']} "
                  f"s/V={res['s_synth']/res['V_synth']:.4f}", flush=True)

    out_rows = []
    for name, m in meta.items():
        sub = [r["s_synth"] for r in rows if r["corpus"] == name]
        if not sub:
            continue
        s_med = float(np.median(sub))
        out_rows.append({"corpus": name, "lang": m["lang"], "tokens": m["tokens"],
                         "V": m["V"], "zm_b": round(m["zm_b"], 3), "zm_c": round(m["zm_c"], 1),
                         "s": round(s_med, 1), "s_over_V": round(s_med / m["V"], 5)})
    out_rows.sort(key=lambda r: (r["lang"], -r["tokens"]))
    with open(OUT / "f15b_panel.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)

    lines = ["# f15b — cross-language width-law panel (full depth)\n"]
    lines.append("| corpus | lang | tokens | V | ZM b | s | s/V |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for r in out_rows:
        lines.append(f"| {r['corpus']} | {r['lang']} | {r['tokens']:,} | {r['V']:,} | "
                     f"{r['zm_b']} | {r['s']} | **{r['s_over_V']:.4f}** |")
    lines.append("")
    for lang in sorted({r["lang"] for r in out_rows}):
        sv = [r["s_over_V"] for r in out_rows if r["lang"] == lang]
        lines.append(f"- {lang}: n={len(sv)}, median s/V = {float(np.median(sv)):.4f}")
    all_sv = [r["s_over_V"] for r in out_rows]
    lV = np.log([r["V"] for r in out_rows]); lS = np.log([r["s"] for r in out_rows])
    X = np.column_stack([np.ones_like(lV), lV])
    b, *_ = np.linalg.lstsq(X, lS, rcond=None)
    r2 = 1 - float(np.sum((lS - X @ b) ** 2)) / float(np.sum((lS - lS.mean()) ** 2))
    lines.append("")
    lines.append(f"- panel: {len(out_rows)} corpora, median s/V = {float(np.median(all_sv)):.4f} "
                 f"(EN reference 0.0118)")
    lines.append(f"- pooled s ~ V: slope {b[1]:.3f}, R2 {r2:.3f}")
    (OUT / "f15b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("summary written", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
