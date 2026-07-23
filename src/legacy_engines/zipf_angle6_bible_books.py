from __future__ import annotations

import csv
import importlib.util
import json
import math
from collections import Counter
from pathlib import Path
import re

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
BASE_PATH = ROOT / "zipf_seam_mandelbrot_pmf.py"
SOFTK_PATH = ROOT / "zipf_seam_mandelbrot_softk.py"
OUTDIR = ROOT / "results" / "zipf_angle6_bible_books"
BIBLE_PATH = ROOT / "data" / "zipf" / "pg10.txt"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"

TRAIN_STARTS = 5
FULL_STARTS = 3
TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")

BOOK_HEADINGS = [
    ("The First Book of Moses: Called Genesis", "Genesis"),
    ("The Second Book of Moses: Called Exodus", "Exodus"),
    ("The Third Book of Moses: Called Leviticus", "Leviticus"),
    ("The Fourth Book of Moses: Called Numbers", "Numbers"),
    ("The Fifth Book of Moses: Called Deuteronomy", "Deuteronomy"),
    ("The Book of Joshua", "Joshua"),
    ("The Book of Judges", "Judges"),
    ("The Book of Ruth", "Ruth"),
    ("The First Book of Samuel", "1 Samuel"),
    ("The Second Book of Samuel", "2 Samuel"),
    ("The First Book of the Kings", "1 Kings"),
    ("The Second Book of the Kings", "2 Kings"),
    ("The First Book of the Chronicles", "1 Chronicles"),
    ("The Second Book of the Chronicles", "2 Chronicles"),
    ("Ezra", "Ezra"),
    ("The Book of Nehemiah", "Nehemiah"),
    ("The Book of Esther", "Esther"),
    ("The Book of Job", "Job"),
    ("The Book of Psalms", "Psalms"),
    ("The Proverbs", "Proverbs"),
    ("Ecclesiastes", "Ecclesiastes"),
    ("The Song of Solomon", "Song of Solomon"),
    ("The Book of the Prophet Isaiah", "Isaiah"),
    ("The Book of the Prophet Jeremiah", "Jeremiah"),
    ("The Lamentations of Jeremiah", "Lamentations"),
    ("The Book of the Prophet Ezekiel", "Ezekiel"),
    ("The Book of Daniel", "Daniel"),
    ("Hosea", "Hosea"),
    ("Joel", "Joel"),
    ("Amos", "Amos"),
    ("Obadiah", "Obadiah"),
    ("Jonah", "Jonah"),
    ("Micah", "Micah"),
    ("Nahum", "Nahum"),
    ("Habakkuk", "Habakkuk"),
    ("Zephaniah", "Zephaniah"),
    ("Haggai", "Haggai"),
    ("Zechariah", "Zechariah"),
    ("Malachi", "Malachi"),
    ("The Gospel According to Saint Matthew", "Matthew"),
    ("The Gospel According to Saint Mark", "Mark"),
    ("The Gospel According to Saint Luke", "Luke"),
    ("The Gospel According to Saint John", "John"),
    ("The Acts of the Apostles", "Acts"),
    ("The Epistle of Paul the Apostle to the Romans", "Romans"),
    ("The First Epistle of Paul the Apostle to the Corinthians", "1 Corinthians"),
    ("The Second Epistle of Paul the Apostle to the Corinthians", "2 Corinthians"),
    ("The Epistle of Paul the Apostle to the Galatians", "Galatians"),
    ("The Epistle of Paul the Apostle to the Ephesians", "Ephesians"),
    ("The Epistle of Paul the Apostle to the Philippians", "Philippians"),
    ("The Epistle of Paul the Apostle to the Colossians", "Colossians"),
    ("The First Epistle of Paul the Apostle to the Thessalonians", "1 Thessalonians"),
    ("The Second Epistle of Paul the Apostle to the Thessalonians", "2 Thessalonians"),
    ("The First Epistle of Paul the Apostle to Timothy", "1 Timothy"),
    ("The Second Epistle of Paul the Apostle to Timothy", "2 Timothy"),
    ("The Epistle of Paul the Apostle to Titus", "Titus"),
    ("The Epistle of Paul the Apostle to Philemon", "Philemon"),
    ("The Epistle of Paul the Apostle to the Hebrews", "Hebrews"),
    ("The General Epistle of James", "James"),
    ("The First Epistle General of Peter", "1 Peter"),
    ("The Second General Epistle of Peter", "2 Peter"),
    ("The First Epistle General of John", "1 John"),
    ("The Second Epistle General of John", "2 John"),
    ("The Third Epistle General of John", "3 John"),
    ("The General Epistle of Jude", "Jude"),
    ("The Revelation of Saint John the Divine", "Revelation"),
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_angle6_common")
base = load_module(BASE_PATH, "zipf_angle6_base")
softk_mod = load_module(SOFTK_PATH, "zipf_angle6_softk")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def tokenize_text(text: str):
    return TOKEN_RE.findall(text.lower())


def build_dataset_from_text(text: str) -> dict:
    tokens = tokenize_text(text)
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
        "token_count": len(tokens),
        "unique_words": len(freqs),
    }


def extract_bible_books() -> list[dict]:
    text = BIBLE_PATH.read_text(encoding="utf-8", errors="ignore")
    text = common.strip_gutenberg_boilerplate(text)
    lines = [line.rstrip("\r") for line in text.splitlines()]
    heading_to_index = {}
    for heading, short in BOOK_HEADINGS:
        indices = [i for i, line in enumerate(lines) if line.strip() == heading]
        if not indices:
            raise ValueError(f"Could not find heading: {heading}")
        heading_to_index[heading] = max(indices)
    ordered = sorted(((heading_to_index[heading], heading, short) for heading, short in BOOK_HEADINGS), key=lambda item: item[0])
    books = []
    for idx, (start_idx, heading, short) in enumerate(ordered):
        end_idx = ordered[idx + 1][0] if idx + 1 < len(ordered) else len(lines)
        segment = "\n".join(lines[start_idx:end_idx]).strip()
        books.append({"heading": heading, "name": short, "text": segment})
    return books


def minimal_base_row(params: dict) -> dict:
    return {"models": {"seam": {"params": params}}}


def load_whole_bible_softk():
    rows = json.loads(SOFTK_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    for row in rows:
        if row["slug"] == "king_james_bible":
            return row
    raise KeyError("king_james_bible not found in soft-k summary")


def analyze_book(book: dict, book_index: int) -> dict:
    dataset = build_dataset_from_text(book["text"])
    split = base.split_counts_by_train_rank(dataset, seed=base.SPLIT_SEED + 50000 + book_index, train_fraction=base.TRAIN_FRACTION)
    ranks = np.arange(1, len(split["train_counts"]) + 1, dtype=np.float64)

    zipf_fit = base.fit_zipf_rank(split["train_counts"], ranks)
    zm_fit = base.fit_zm_rank(split["train_counts"], ranks, dataset_hint=dataset)
    moe_fit = base.fit_moe_rank(split["train_counts"], ranks, zipf_fit)
    free_fit = base.fit_seam_rank(split["train_counts"], ranks, vocab_size=split["vocab_size"], reranked_row=None, seed=base.SPLIT_SEED + 51000 + book_index)
    base_row = minimal_base_row(free_fit["params"])

    lambda_rows = []
    prev_params = None
    for lambda_index, lam in enumerate(softk_mod.LAMBDAS, start=1):
        fit = softk_mod.fit_softk(
            split["train_counts"],
            ranks,
            split["vocab_size"],
            reranked_row=None,
            base_row=base_row,
            lam=lam,
            seed=base.SPLIT_SEED + 52000 + book_index * 10 + lambda_index,
            n_starts=TRAIN_STARTS,
            lambda_prev_params=prev_params,
        )
        prev_params = fit["params_vec"]
        evals = softk_mod.evaluate_fit(fit, split["test_counts"])
        lambda_rows.append({"lambda": float(lam), "fit": fit, "test_avg_nll": float(evals["test_avg_nll"]), "test_loglike": float(evals["test_loglike"])})

    best_lambda_row = min(lambda_rows, key=lambda row: row["test_avg_nll"])
    full = softk_mod.full_context_softk(
        dataset,
        reranked_row=None,
        base_row=base_row,
        lam=float(best_lambda_row["lambda"]),
        seed=base.SPLIT_SEED + 53000 + book_index,
        n_starts=FULL_STARTS,
        lambda_prev_params=best_lambda_row["fit"]["params_vec"],
    )
    step2 = base.exact_step2_on_residual(dataset, full["prediction_log"])

    zipf_eval = base.evaluate_test_fit(zipf_fit, split["test_counts"])
    zm_eval = base.evaluate_test_fit(zm_fit, split["test_counts"])
    moe_eval = base.evaluate_test_fit(moe_fit, split["test_counts"])

    return {
        "name": book["name"],
        "heading": book["heading"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "test_token_count": int(np.sum(split["test_counts"])),
        "softk_lambda": float(best_lambda_row["lambda"]),
        "softk_test_loglike": float(best_lambda_row["test_loglike"]),
        "softk_test_avg_nll": float(best_lambda_row["test_avg_nll"]),
        "zm_test_avg_nll": float(zm_eval["test_avg_nll"]),
        "moe_test_avg_nll": float(moe_eval["test_avg_nll"]),
        "zipf_test_avg_nll": float(zipf_eval["test_avg_nll"]),
        "softk_minus_zm": float(best_lambda_row["test_avg_nll"] - zm_eval["test_avg_nll"]),
        "softk_minus_moe": float(best_lambda_row["test_avg_nll"] - moe_eval["test_avg_nll"]),
        "step2_gain": float(step2["gain"]),
        "step2_helps": bool(step2["helps"]),
        "step2_expr": step2["expr"],
        "softk_rmse": float(full["rmse"]),
    }


def build_summary(rows: list[dict], whole_bible_softk: dict) -> dict:
    total_test_loglike = float(sum(row["softk_test_loglike"] for row in rows))
    total_test_tokens = int(sum(row["test_token_count"] for row in rows))
    return {
        "rows": rows,
        "counts": {
            "n_books": len(rows),
            "step2_help_count": int(sum(row["step2_helps"] for row in rows)),
            "softk_beats_zm": int(sum(row["softk_minus_zm"] < 0.0 for row in rows)),
            "softk_beats_moe": int(sum(row["softk_minus_moe"] < 0.0 for row in rows)),
        },
        "medians": {
            "softk_minus_zm": float(np.median([row["softk_minus_zm"] for row in rows])),
            "softk_minus_moe": float(np.median([row["softk_minus_moe"] for row in rows])),
            "step2_gain": float(np.median([row["step2_gain"] for row in rows])),
        },
        "aggregate": {
            "whole_bible_softk_test_avg_nll": float(whole_bible_softk["best_lambda"]["test_avg_nll"]),
            "bookwise_softk_total_test_loglike": total_test_loglike,
            "bookwise_softk_test_token_count": total_test_tokens,
            "bookwise_softk_test_avg_nll": float(-total_test_loglike / total_test_tokens),
        },
    }


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "name",
        "token_count",
        "vocab_size",
        "softk_lambda",
        "softk_test_avg_nll",
        "zm_test_avg_nll",
        "moe_test_avg_nll",
        "softk_minus_zm",
        "softk_minus_moe",
        "step2_gain",
        "step2_help",
        "step2_expr",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "name": row["name"],
                    "token_count": row["token_count"],
                    "vocab_size": row["vocab_size"],
                    "softk_lambda": row["softk_lambda"],
                    "softk_test_avg_nll": row["softk_test_avg_nll"],
                    "zm_test_avg_nll": row["zm_test_avg_nll"],
                    "moe_test_avg_nll": row["moe_test_avg_nll"],
                    "softk_minus_zm": row["softk_minus_zm"],
                    "softk_minus_moe": row["softk_minus_moe"],
                    "step2_gain": row["step2_gain"],
                    "step2_help": int(row["step2_helps"]),
                    "step2_expr": row["step2_expr"],
                }
            )


def build_report(summary: dict, aggregate_test_tokens: int):
    lines = [
        "# Angle 6: King James Bible by Book",
        "",
        "- Each canonical book was fit independently with the soft-k Seam-Mandelbrot PMF.",
        f"- books analyzed: `{summary['counts']['n_books']}`",
        f"- per-book step-2 help count: `{summary['counts']['step2_help_count']}` / `{summary['counts']['n_books']}`",
        f"- per-book soft-k beats ZM on held-out avg NLL: `{summary['counts']['softk_beats_zm']}` / `{summary['counts']['n_books']}`",
        f"- per-book soft-k beats MOE on held-out avg NLL: `{summary['counts']['softk_beats_moe']}` / `{summary['counts']['n_books']}`",
        f"- median per-book soft-k minus ZM held-out avg NLL: `{summary['medians']['softk_minus_zm']:.12f}`",
        f"- median per-book soft-k minus MOE held-out avg NLL: `{summary['medians']['softk_minus_moe']:.12f}`",
        f"- median per-book step-2 gain: `{summary['medians']['step2_gain']:.12f}`",
        "",
        f"- whole-Bible single-fit soft-k held-out avg NLL: `{summary['aggregate']['whole_bible_softk_test_avg_nll']:.12f}`",
        f"- aggregate per-book soft-k held-out avg NLL: `{summary['aggregate']['bookwise_softk_test_avg_nll']:.12f}`",
        "",
        "| book | soft-k - ZM | soft-k - MOE | step-2 gain | step-2 expr |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for row in sorted(summary["rows"], key=lambda item: item["softk_minus_moe"]):
        lines.append(
            f"| {row['name']} | {row['softk_minus_zm']:.6f} | {row['softk_minus_moe']:.6f} | {row['step2_gain']:.6f} | {row['step2_expr']} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    books = extract_bible_books()
    whole_bible_softk = load_whole_bible_softk()
    rows = []
    for book_index, book in enumerate(books, start=1):
        row = analyze_book(book, book_index)
        rows.append(row)
    summary = build_summary(rows, whole_bible_softk)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary, summary["aggregate"]["bookwise_softk_test_token_count"]), encoding="utf-8")
    write_csv(rows, OUTDIR / "bible_books_table.csv")


if __name__ == "__main__":
    main()
