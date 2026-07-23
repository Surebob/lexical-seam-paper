from __future__ import annotations

import argparse
import csv
import gc
import html
from html.parser import HTMLParser
import importlib.util
import json
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path

import jieba
import numpy as np
import regex


ROOT = Path("/Volumes/External2TB/emlexperiment")
DATA_DIR = ROOT / "data" / "zipf_multilang"
OUTDIR = ROOT / "results" / "zipf_multilang"

COMMON_PATH = ROOT / "zipf_analysis_common.py"
ENRICHED_PATH = ROOT / "eml_zipf_enriched_search.py"
MOE_PATH = ROOT / "zipf_moezipf_comparison.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
REDUCTION_PATH = ROOT / "zipf_model_reduction.py"
CORRECT_MODEL_PATH = ROOT / "zipf_correct_model.py"

BEAM_WIDTH = 50
MAX_STEPS = 10
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
CONSTANT_VARIANCE_THRESHOLD = 1e-10
RERANKED_STARTS = 30
RERANKED_MAX_NFEV = 8000
REDUCED_STARTS = 30

BREGMAN_EXPR_HIGH = "sub[sub[x,1],log[x]]"
BREGMAN_EXPR_LOW = "eml[sub[x,1],eml[x,1]]"

UNICODE_WORD_RE = regex.compile(r"[\p{L}\p{M}]+(?:['’\-][\p{L}\p{M}]+)*", flags=regex.VERSION1)
HAN_RE = regex.compile(r"\p{Script=Han}", flags=regex.VERSION1)


CORPORA = [
    {
        "slug": "russian_war_and_peace",
        "language": "Russian",
        "name": "War and Peace (Russian, Wikisource)",
        "tokenizer": "unicode_words",
        "source_type": "wikisource_render",
        "max_tokens": 150000,
        "urls": [
            "https://ru.wikisource.org/wiki/%D0%92%D0%BE%D0%B9%D0%BD%D0%B0_%D0%B8_%D0%BC%D0%B8%D1%80_(%D0%A2%D0%BE%D0%BB%D1%81%D1%82%D0%BE%D0%B9)/%D0%A2%D0%BE%D0%BC_1?action=render",
            "https://ru.wikisource.org/wiki/%D0%92%D0%BE%D0%B9%D0%BD%D0%B0_%D0%B8_%D0%BC%D0%B8%D1%80_(%D0%A2%D0%BE%D0%BB%D1%81%D1%82%D0%BE%D0%B9)/%D0%A2%D0%BE%D0%BC_2?action=render",
            "https://ru.wikisource.org/wiki/%D0%92%D0%BE%D0%B9%D0%BD%D0%B0_%D0%B8_%D0%BC%D0%B8%D1%80_(%D0%A2%D0%BE%D0%BB%D1%81%D1%82%D0%BE%D0%B9)/%D0%A2%D0%BE%D0%BC_3?action=render",
            "https://ru.wikisource.org/wiki/%D0%92%D0%BE%D0%B9%D0%BD%D0%B0_%D0%B8_%D0%BC%D0%B8%D1%80_(%D0%A2%D0%BE%D0%BB%D1%81%D1%82%D0%BE%D0%B9)/%D0%A2%D0%BE%D0%BC_4?action=render",
            "https://ru.wikisource.org/wiki/%D0%92%D0%BE%D0%B9%D0%BD%D0%B0_%D0%B8_%D0%BC%D0%B8%D1%80_(%D0%A2%D0%BE%D0%BB%D1%81%D1%82%D0%BE%D0%B9)/%D0%AD%D0%BF%D0%B8%D0%BB%D0%BE%D0%B3?action=render",
        ],
    },
    {
        "slug": "mandarin_three_kingdoms",
        "language": "Mandarin",
        "name": "Romance of the Three Kingdoms (Chinese, Gutenberg 23950)",
        "tokenizer": "jieba",
        "source_type": "gutenberg_text",
        "max_tokens": 80000,
        "urls": [
            "https://www.gutenberg.org/cache/epub/23950/pg23950.txt",
        ],
    },
    {
        "slug": "arabic_1001_nights",
        "language": "Arabic",
        "name": "One Thousand and One Nights (Arabic, Wikisource)",
        "tokenizer": "unicode_words",
        "source_type": "wikisource_render",
        "max_tokens": 150000,
        "urls": [
            "https://ar.wikisource.org/wiki/%D8%A3%D9%84%D9%81_%D9%84%D9%8A%D9%84%D8%A9_%D9%88%D9%84%D9%8A%D9%84%D8%A9/%D8%A7%D9%84%D8%AC%D8%B2%D8%A1_%D8%A7%D9%84%D8%A3%D9%88%D9%84?action=render",
            "https://ar.wikisource.org/wiki/%D8%A3%D9%84%D9%81_%D9%84%D9%8A%D9%84%D8%A9_%D9%88%D9%84%D9%8A%D9%84%D8%A9/%D8%A7%D9%84%D8%AC%D8%B2%D8%A1_%D8%A7%D9%84%D8%AB%D8%A7%D9%86%D9%8A?action=render",
            "https://ar.wikisource.org/wiki/%D8%A3%D9%84%D9%81_%D9%84%D9%8A%D9%84%D8%A9_%D9%88%D9%84%D9%8A%D9%84%D8%A9/%D8%A7%D9%84%D8%AC%D8%B2%D8%A1_%D8%A7%D9%84%D8%AB%D8%A7%D9%84%D8%AB?action=render",
            "https://ar.wikisource.org/wiki/%D8%A3%D9%84%D9%81_%D9%84%D9%8A%D9%84%D8%A9_%D9%88%D9%84%D9%8A%D9%84%D8%A9/%D8%A7%D9%84%D8%AC%D8%B2%D8%A1_%D8%A7%D9%84%D8%B1%D8%A7%D8%A8%D8%B9?action=render",
        ],
    },
    {
        "slug": "latin_gallic_wars",
        "language": "Latin",
        "name": "De Bello Gallico (Latin, Gutenberg 218 + 18837)",
        "tokenizer": "unicode_words",
        "source_type": "gutenberg_text",
        "max_tokens": None,
        "urls": [
            "https://www.gutenberg.org/cache/epub/218/pg218.txt",
            "https://www.gutenberg.org/cache/epub/18837/pg18837.txt",
        ],
    },
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_multilang_common")
enriched = load_module(ENRICHED_PATH, "zipf_multilang_enriched")
moe = load_module(MOE_PATH, "zipf_multilang_moe")
reranked = load_module(RERANKED_PATH, "zipf_multilang_reranked")
reduction = load_module(REDUCTION_PATH, "zipf_multilang_reduction")
correct_model = load_module(CORRECT_MODEL_PATH, "zipf_multilang_correct_model")


def parse_args():
    parser = argparse.ArgumentParser(description="Multilingual Zipf extension")
    parser.add_argument("--only-slug", type=str, default=None)
    return parser.parse_args()


class TextExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)

    def text(self) -> str:
        return " ".join(self.parts)


def download_text(url: str, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return dest.read_text(encoding="utf-8", errors="ignore")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        raw = response.read()
        charset = response.headers.get_content_charset() or "utf-8"
    text = raw.decode(charset, errors="ignore")
    dest.write_text(text, encoding="utf-8")
    return text


def strip_rendered_html(html_text: str) -> str:
    stripped = re.sub(r"(?is)<(script|style).*?>.*?</\\1>", " ", html_text)
    parser = TextExtractor()
    parser.feed(stripped)
    text = html.unescape(parser.text())
    text = re.sub(r"\s+", " ", text)
    return text


def download_and_clean(spec: dict) -> str:
    parts = []
    source_dir = DATA_DIR / spec["slug"]
    for idx, url in enumerate(spec["urls"], start=1):
        ext = ".html" if spec["source_type"] == "wikisource_render" else ".txt"
        raw_text = download_text(url, source_dir / f"source_{idx}{ext}")
        if spec["source_type"] == "wikisource_render":
            cleaned = strip_rendered_html(raw_text)
        else:
            cleaned = common.strip_gutenberg_boilerplate(raw_text)
        parts.append(cleaned)
    combined = "\n".join(parts)
    (source_dir / "combined_clean.txt").write_text(combined, encoding="utf-8")
    (source_dir / "sources.json").write_text(json.dumps(spec["urls"], indent=2, ensure_ascii=False), encoding="utf-8")
    return combined


def tokenize_unicode_words(text: str) -> list[str]:
    return [tok.casefold() for tok in UNICODE_WORD_RE.findall(text)]


def tokenize_jieba_words(text: str) -> list[str]:
    tokens = []
    for token in jieba.cut(text, cut_all=False):
        token = token.strip()
        if not token:
            continue
        if not HAN_RE.search(token):
            continue
        token = regex.sub(r"[^\p{Script=Han}]+", "", token)
        if token:
            tokens.append(token)
    return tokens


def tokenize_text(spec: dict, text: str) -> list[str]:
    if spec["tokenizer"] == "jieba":
        return tokenize_jieba_words(text)
    return tokenize_unicode_words(text)


def build_dataset(tokens: list[str]) -> dict:
    counts = Counter(tokens)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)
    return {
        "tokens": tokens,
        "counts": counts,
        "ranked": ranked,
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": log_rank,
        "log_freq": log_freq,
        "unique_words": int(len(freqs)),
        "token_count": int(len(tokens)),
    }


def rmse_bic(rmse_value: float, n_obs: int, p: int) -> float:
    mse = max(float(rmse_value) ** 2, 1e-300)
    return float(p * math.log(n_obs) + n_obs * math.log(mse))


def run_step2_search(dataset: dict) -> dict:
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    zm_fit = common.fit_zipf_mandelbrot(dataset["ranks"], dataset["log_freq"])
    target = dataset["log_freq"] - zm_fit["prediction"]
    search = enriched.run_search(
        x,
        target,
        BEAM_WIDTH,
        MAX_STEPS,
        KEEP_ALL_UNTIL_STEP,
        DIVERSITY_WEIGHT,
        CONSTANT_VARIANCE_THRESHOLD,
    )
    step2 = next((step for step in search["steps"] if step["step"] == 2), None)
    step2_top = step2["top_candidates"] if step2 else []
    exprs = [row["expr"] for row in step2_top]
    return {
        "zm_fit": zm_fit,
        "step2_winner": step2_top[0]["expr"] if step2_top else None,
        "step2_rmse": step2_top[0]["rmse"] if step2_top else None,
        "step2_top5": step2_top[:5],
        "bregman_in_beam": bool(BREGMAN_EXPR_HIGH in exprs or BREGMAN_EXPR_LOW in exprs),
        "high_bregman_in_beam": bool(BREGMAN_EXPR_HIGH in exprs),
        "low_bregman_in_beam": bool(BREGMAN_EXPR_LOW in exprs),
    }


def fit_baselines(dataset: dict, idx: int) -> dict:
    freqs = dataset["freqs"].astype(np.int64)
    n_obs = int(dataset["unique_words"])
    observed_max = int(freqs.max())

    zm_fit = common.fit_zipf_mandelbrot(dataset["ranks"], dataset["log_freq"])
    zm_rmse = float(zm_fit["rmse"])
    zm_bic = rmse_bic(zm_rmse, n_obs, p=3)

    zipf_mle = moe.fit_zipf_mle(freqs)
    moe_fit = moe.fit_moe_mle(freqs, zipf_mle)
    moe_pred = moe.rank_curve_from_tail(
        lambda x: moe.moe_tail_ge(moe_fit["alpha"], moe_fit["beta"], x),
        n_obs,
        observed_max,
    )
    moe_rmse = float(moe.rmse_log_rank(freqs, moe_pred))
    moe_bic = rmse_bic(moe_rmse, n_obs, p=2)

    reduced = reduction.fit_reduced_model(dataset, fix_c2_zero=False, n_starts=REDUCED_STARTS, seed=20270400 + idx)
    reduced_summary = reduction.summarize_fit(dataset, reduced)
    sqrt_v_rmse = float(reduced_summary["rmse"])
    sqrt_v_bic = rmse_bic(sqrt_v_rmse, n_obs, p=7)

    reranked_best, _ = reranked.fit_reranked_model(dataset, n_starts=RERANKED_STARTS, max_nfev=RERANKED_MAX_NFEV)
    reranked_rmse = float(reranked_best["rmse"])
    reranked_bic = rmse_bic(reranked_rmse, n_obs, p=8)
    reranked_params = reranked.summarize_params(reranked_best["params"], dataset["unique_words"])

    return {
        "zm_fit": zm_fit,
        "zm_rmse": zm_rmse,
        "zm_bic": zm_bic,
        "moe_fit": moe_fit,
        "moe_rmse": moe_rmse,
        "moe_bic": moe_bic,
        "sqrt_v_rmse": sqrt_v_rmse,
        "sqrt_v_bic": sqrt_v_bic,
        "sqrt_v_params": reduced_summary["params"],
        "reranked_rmse": reranked_rmse,
        "reranked_bic": reranked_bic,
        "reranked_params": reranked_params,
    }


def artifact_note(spec: dict, top_tokens: list[tuple[str, int]]) -> str:
    words = [word for word, _ in top_tokens]
    if spec["language"] == "Mandarin":
        if any(len(word) == 1 for word in words[:5]):
            return "Top Chinese tokens are mostly single characters; segmentation may be too fine."
    if any("wikipedia" in word.lower() for word in words[:10]):
        return "Wikipedia/Wikisource metadata leaked into the top tokens."
    return "Top tokens look lexical rather than markup-heavy."


def write_partial(summary: dict) -> None:
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")


def write_csv(rows: list[dict]) -> None:
    with (OUTDIR / "multilang_table.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "language",
                "corpus",
                "token_count",
                "raw_token_count",
                "vocab_size",
                "single_zm_rmse",
                "moezipf_rmse",
                "smooth_7param_rmse",
                "smooth_8param_rmse",
                "step2_winner",
                "bregman_in_beam",
                "transition_fraction",
                "artifact_note",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_report(rows: list[dict], failures: list[dict]) -> str:
    lines = [
        "# Multilingual Zipf Extension",
        "",
        "| Language | Corpus | V | Single ZM RMSE | 8-param RMSE | Step-2 winner | Bregman in beam? | log(k)/log(V) |",
        "| --- | --- | ---: | ---: | ---: | --- | --- | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['language']} | {row['corpus']} | {row['vocab_size']} | {row['single_zm_rmse']:.12f} | {row['smooth_8param_rmse']:.12f} | {row['step2_winner']} | {row['bregman_in_beam']} | {row['transition_fraction']:.6f} |"
        )
    lines.extend(["", "## Per-corpus notes", ""])
    for row in rows:
        lines.extend(
            [
                f"### {row['language']} — {row['corpus']}",
                "",
                f"- tokens: `{row['token_count']}`",
                f"- raw tokens before any truncation: `{row['raw_token_count']}`",
                f"- vocab: `{row['vocab_size']}`",
                f"- top tokens: `{row['top_tokens']}`",
                f"- artifact note: {row['artifact_note']}",
                f"- step-2 top 5: `{row['step2_top5']}`",
                f"- single ZM RMSE/BIC: `{row['single_zm_rmse']:.12f}` / `{row['single_zm_bic']:.6f}`",
                f"- MOEZipf RMSE/BIC: `{row['moezipf_rmse']:.12f}` / `{row['moezipf_bic']:.6f}`",
                f"- smooth 7-param sqrt(V) RMSE/BIC: `{row['smooth_7param_rmse']:.12f}` / `{row['smooth_7param_bic']:.6f}`",
                f"- smooth 8-param RMSE/BIC: `{row['smooth_8param_rmse']:.12f}` / `{row['smooth_8param_bic']:.6f}`",
                f"- smooth 8-param beats single ZM: `{row['smooth_8param_rmse'] < row['single_zm_rmse']}`",
                f"- transition fraction: `{row['transition_fraction']:.6f}`",
                f"- sources: `{row['sources']}`",
                "",
            ]
        )
    if failures:
        lines.extend(["## Failures", ""])
        for failure in failures:
            lines.extend(
                [
                    f"### {failure['language']} — {failure['corpus']}",
                    "",
                    f"- error: `{failure['error']}`",
                    f"- sources: `{failure['sources']}`",
                    "",
                ]
            )
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    failures = []
    corpus_list = [spec for spec in CORPORA if args.only_slug in (None, spec["slug"])]
    for idx, spec in enumerate(corpus_list, start=1):
        print(f"[multilang] processing {spec['slug']}")
        try:
            cleaned_text = download_and_clean(spec)
            raw_tokens = tokenize_text(spec, cleaned_text)
            max_tokens = spec.get("max_tokens")
            tokens = raw_tokens[:max_tokens] if max_tokens else raw_tokens
            dataset = build_dataset(tokens)
            search_info = run_step2_search(dataset)
            model_info = fit_baselines(dataset, idx)

            top_tokens = dataset["ranked"][:15]
            row = {
                "slug": spec["slug"],
                "language": spec["language"],
                "corpus": spec["name"],
                "sources": spec["urls"],
                "token_count": dataset["token_count"],
                "raw_token_count": len(raw_tokens),
                "vocab_size": dataset["unique_words"],
                "top_tokens": top_tokens,
                "artifact_note": artifact_note(spec, top_tokens)
                + (
                    f" Truncated to first {max_tokens} tokens for feasibility."
                    if max_tokens and len(raw_tokens) > max_tokens
                    else ""
                ),
                "single_zm_rmse": model_info["zm_rmse"],
                "single_zm_bic": model_info["zm_bic"],
                "moezipf_rmse": model_info["moe_rmse"],
                "moezipf_bic": model_info["moe_bic"],
                "smooth_7param_rmse": model_info["sqrt_v_rmse"],
                "smooth_7param_bic": model_info["sqrt_v_bic"],
                "smooth_8param_rmse": model_info["reranked_rmse"],
                "smooth_8param_bic": model_info["reranked_bic"],
                "transition_fraction": float(model_info["reranked_params"]["transition_fraction"]),
                "step2_winner": search_info["step2_winner"],
                "step2_rmse": search_info["step2_rmse"],
                "step2_top5": search_info["step2_top5"],
                "bregman_in_beam": search_info["bregman_in_beam"],
                "high_bregman_in_beam": search_info["high_bregman_in_beam"],
                "low_bregman_in_beam": search_info["low_bregman_in_beam"],
            }
            rows.append(row)
        except Exception as exc:
            failures.append(
                {
                    "slug": spec["slug"],
                    "language": spec["language"],
                    "corpus": spec["name"],
                    "sources": spec["urls"],
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
        finally:
            gc.collect()
        summary = {"rows": rows, "failures": failures}
        write_partial(summary)
        write_csv(rows)
        (OUTDIR / "report.md").write_text(build_report(rows, failures), encoding="utf-8")


if __name__ == "__main__":
    main()
