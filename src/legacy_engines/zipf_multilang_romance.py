from __future__ import annotations

import argparse
import gc
import importlib.util
import json
from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
BASE_PATH = ROOT / "zipf_multilang.py"
DATA_DIR = ROOT / "data" / "zipf_multilang_romance"
OUTDIR = ROOT / "results" / "zipf_multilang_romance"
PREVIOUS_SUMMARY_PATH = ROOT / "results" / "zipf_multilang" / "summary.json"


CORPORA = [
    {
        "slug": "french_les_miserables",
        "language": "French",
        "name": "Les Misérables (French, Gutenberg)",
        "tokenizer": "unicode_words",
        "source_type": "gutenberg_text",
        "max_tokens": 150000,
        "urls": [
            "https://www.gutenberg.org/ebooks/17489.txt.utf-8",
            "https://www.gutenberg.org/ebooks/17493.txt.utf-8",
            "https://www.gutenberg.org/ebooks/17494.txt.utf-8",
            "https://www.gutenberg.org/ebooks/17518.txt.utf-8",
            "https://www.gutenberg.org/ebooks/17519.txt.utf-8",
        ],
    },
    {
        "slug": "spanish_don_quixote",
        "language": "Spanish",
        "name": "Don Quijote (Spanish, Gutenberg)",
        "tokenizer": "unicode_words",
        "source_type": "gutenberg_text",
        "max_tokens": 150000,
        "urls": [
            "https://www.gutenberg.org/ebooks/2000.txt.utf-8",
        ],
    },
    {
        "slug": "dutch_max_havelaar",
        "language": "Dutch",
        "name": "Max Havelaar (Dutch, Gutenberg)",
        "tokenizer": "unicode_words",
        "source_type": "gutenberg_text",
        "max_tokens": 150000,
        "urls": [
            "https://www.gutenberg.org/ebooks/11024.txt.utf-8",
        ],
    },
]

ORDER = [
    "russian_war_and_peace",
    "mandarin_three_kingdoms",
    "arabic_1001_nights",
    "latin_gallic_wars",
    "french_les_miserables",
    "spanish_don_quixote",
    "dutch_max_havelaar",
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


base = load_module(BASE_PATH, "zipf_multilang_romance_base")
base.DATA_DIR = DATA_DIR


def parse_args():
    parser = argparse.ArgumentParser(description="Romance/Germanic multilingual Zipf extension")
    parser.add_argument("--only-slug", type=str, default=None)
    return parser.parse_args()


def load_previous_rows():
    if not PREVIOUS_SUMMARY_PATH.exists():
        return [], []
    obj = json.loads(PREVIOUS_SUMMARY_PATH.read_text(encoding="utf-8"))
    return obj.get("rows", []), obj.get("failures", [])


def load_current_rows():
    rows = []
    failures = []
    per_dir = OUTDIR / "per_corpus"
    if not per_dir.exists():
        return rows, failures
    for path in sorted(per_dir.glob("*_summary.json")):
        obj = json.loads(path.read_text(encoding="utf-8"))
        rows.extend(obj.get("rows", []))
        failures.extend(obj.get("failures", []))
    return rows, failures


def ordered_rows(rows):
    by_slug = {row["slug"]: row for row in rows}
    return [by_slug[slug] for slug in ORDER if slug in by_slug]


def write_combined_outputs():
    prev_rows, prev_failures = load_previous_rows()
    new_rows, new_failures = load_current_rows()
    rows = ordered_rows(prev_rows + new_rows)
    failures = prev_failures + new_failures
    summary = {"rows": rows, "failures": failures}
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    base.OUTDIR = OUTDIR
    base.write_csv(rows)
    (OUTDIR / "report.md").write_text(base.build_report(rows, failures), encoding="utf-8")


def process_corpus(spec: dict, idx: int):
    cleaned_text = base.download_and_clean(spec)
    raw_tokens = base.tokenize_text(spec, cleaned_text)
    max_tokens = spec.get("max_tokens")
    tokens = raw_tokens[:max_tokens] if max_tokens else raw_tokens
    dataset = base.build_dataset(tokens)
    search_info = base.run_step2_search(dataset)
    model_info = base.fit_baselines(dataset, idx)
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
        "artifact_note": base.artifact_note(spec, top_tokens)
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
    return row


def main():
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "per_corpus").mkdir(parents=True, exist_ok=True)
    corpus_list = [spec for spec in CORPORA if args.only_slug in (None, spec["slug"])]
    for idx, spec in enumerate(corpus_list, start=1):
        try:
            row = process_corpus(spec, idx)
            payload = {"rows": [row], "failures": []}
        except Exception as exc:
            payload = {
                "rows": [],
                "failures": [
                    {
                        "slug": spec["slug"],
                        "language": spec["language"],
                        "corpus": spec["name"],
                        "sources": spec["urls"],
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                ],
            }
        (OUTDIR / "per_corpus" / f"{spec['slug']}_summary.json").write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        gc.collect()
        write_combined_outputs()


if __name__ == "__main__":
    main()
