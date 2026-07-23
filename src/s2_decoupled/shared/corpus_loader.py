from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import numpy as np


TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = [
    "*** START OF THE PROJECT GUTENBERG EBOOK",
    "*** START OF THIS PROJECT GUTENBERG EBOOK",
]
END_MARKERS = [
    "*** END OF THE PROJECT GUTENBERG EBOOK",
    "*** END OF THIS PROJECT GUTENBERG EBOOK",
]


def strip_gutenberg_boilerplate(text: str) -> str:
    start = 0
    end = len(text)
    for marker in START_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            line_end = text.find("\n", idx)
            start = line_end + 1 if line_end != -1 else idx
            break
    for marker in END_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            end = idx
            break
    return text[start:end]


def tokenize_text(text: str) -> list[str]:
    return TOKEN_RE.findall(text.lower())


def build_zipf_dataset(corpus_path: Path) -> dict:
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = strip_gutenberg_boilerplate(raw_text)
    counts = Counter(tokenize_text(clean_text))
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    return {
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": np.log(ranks),
        "log_freq": np.log(freqs),
        "token_count": int(sum(counts.values())),
        "unique_words": len(freqs),
    }

