"""F5 fetch — expanded corpus panel: modern English registers, additional languages,
and a surname control. All sources are public research corpora; raw files land in
data/modern/ and data/langext/ (gitignored), with a manifest recording URL + size.

Modern English:
  - Brown corpus (nltk_data mirror)          data/modern/brown.zip
  - Cornell movie-dialogs corpus             data/modern/cornell_movie_dialogs.zip
  - WikiText-103 raw (CC-SA)                 data/modern/wikitext-103-raw-v1.zip
Control:
  - US Census 2010 surnames                  data/modern/census_surnames.zip
Languages (Project Gutenberg, candidate IDs validated by language markers):
  Italian, Finnish, Polish, Portuguese, German, Swedish
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
MODERN = REPO / "data" / "modern"
LANG = REPO / "data" / "langext"
MODERN.mkdir(parents=True, exist_ok=True)
LANG.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "lexical-seam-research/1.0 (corpus fetch for academic replication)"}

MODERN_FILES = [
    ("brown.zip", ["https://raw.githubusercontent.com/nltk/nltk_data/gh-pages/packages/corpora/brown.zip"]),
    ("cornell_movie_dialogs.zip", ["https://www.cs.cornell.edu/~cristian/data/cornell_movie_dialogs_corpus.zip"]),
    ("wikitext-103-raw-v1.zip", [
        "https://wikitext.smerity.com/wikitext-103-raw-v1.zip",
        "https://s3.amazonaws.com/research.metamind.io/wikitext/wikitext-103-raw-v1.zip",
    ]),
    ("census_surnames.zip", ["https://www2.census.gov/topics/genealogy/2010surnames/names.zip"]),
]

GUTEN = [
    ("italian", [45334], " che "),
    ("finnish", [11940], " ja "),
    ("polish", [31536], " nie "),
    ("portuguese", [55752], " que "),
    ("german", [2403, 6343, 34811], " und "),
    ("swedish", [39147, 30078, 5518], " och "),
]


def fetch(url: str, dest: Path, min_bytes=10_000) -> bool:
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        ok = dest.stat().st_size >= min_bytes
        if not ok:
            dest.unlink(missing_ok=True)
        return ok
    except Exception as e:
        print(f"  FAIL {url}: {e}", flush=True)
        dest.unlink(missing_ok=True)
        return False


def main():
    manifest = {}
    for name, urls in MODERN_FILES:
        dest = MODERN / name
        if dest.exists() and dest.stat().st_size > 10_000:
            print(f"have {name} ({dest.stat().st_size/1e6:.1f} MB)", flush=True)
            manifest[name] = {"status": "cached", "bytes": dest.stat().st_size}
            continue
        got = False
        for url in urls:
            print(f"fetching {name} <- {url}", flush=True)
            if fetch(url, dest):
                print(f"  ok ({dest.stat().st_size/1e6:.1f} MB)", flush=True)
                manifest[name] = {"status": "ok", "url": url, "bytes": dest.stat().st_size}
                got = True
                break
        if not got:
            manifest[name] = {"status": "FAILED"}

    for lang, ids, marker in GUTEN:
        dest = LANG / f"{lang}.txt"
        if dest.exists() and dest.stat().st_size > 100_000:
            print(f"have {lang} ({dest.stat().st_size/1e3:.0f} KB)", flush=True)
            manifest[f"lang_{lang}"] = {"status": "cached", "bytes": dest.stat().st_size}
            continue
        chosen = None
        for gid in ids:
            url = f"https://www.gutenberg.org/cache/epub/{gid}/pg{gid}.txt"
            tmp = LANG / f"_try_{lang}_{gid}.txt"
            print(f"trying {lang} pg{gid}", flush=True)
            if not fetch(url, tmp, min_bytes=100_000):
                continue
            sample = tmp.read_text(encoding="utf-8", errors="ignore")[:200_000].lower()
            if marker in sample:
                tmp.rename(dest)
                chosen = {"status": "ok", "gutenberg_id": gid, "bytes": dest.stat().st_size}
                print(f"  accepted pg{gid} ({dest.stat().st_size/1e3:.0f} KB)", flush=True)
                break
            print(f"  rejected pg{gid} (marker {marker!r} not found)", flush=True)
            tmp.unlink(missing_ok=True)
        manifest[f"lang_{lang}"] = chosen or {"status": "FAILED", "tried": ids}

    (MODERN.parent / "expanded_panel_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    n_ok = sum(1 for v in manifest.values() if v.get("status") in ("ok", "cached"))
    print(f"\nmanifest written: {n_ok}/{len(manifest)} items ok", flush=True)


if __name__ == "__main__":
    sys.exit(main())
