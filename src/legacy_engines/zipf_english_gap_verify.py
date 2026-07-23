from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

import importlib.util


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_english_gap_verify"
COMMON_PATH = ROOT / "zipf_analysis_common.py"
MULTILANG_VERIFY_SUMMARY_PATH = ROOT / "results" / "zipf_multilang_verify" / "summary.json"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_english_gap_verify_common")


def high_bregman(x: np.ndarray):
    return (x - 1.0) - np.log(x)


def low_bregman(x: np.ndarray):
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def euclidean(x: np.ndarray):
    return (1.0 - x) ** 2


def formula_rmse(values: np.ndarray, target: np.ndarray) -> float:
    return common.rmse(target, values)


def analyze_corpus(spec: dict):
    summary = common.load_enriched_summary(spec)
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    zm_pred = common.zm_prediction(summary, dataset["ranks"])
    target = dataset["log_freq"] - zm_pred

    winner = common.get_step2_candidate(summary)
    winner_rmse = float(winner["rmse"])
    high_rmse = formula_rmse(high_bregman(x), target)
    low_rmse = formula_rmse(low_bregman(x), target)
    euclid_rmse = formula_rmse(euclidean(x), target)

    return {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "winner_expr": winner["expr"],
        "winner_rmse": winner_rmse,
        "high_bregman_rmse": high_rmse,
        "low_bregman_rmse": low_rmse,
        "euclidean_rmse": euclid_rmse,
        "gap_vs_euclidean": euclid_rmse - winner_rmse,
        "gap_vs_is_bregman": high_rmse - winner_rmse,
        "gap_vs_exp_bregman": low_rmse - winner_rmse,
    }


def load_multilang_medians():
    obj = json.loads(MULTILANG_VERIFY_SUMMARY_PATH.read_text(encoding="utf-8"))
    corpora = obj["corpora"]
    return {
        "median_gap_vs_euclidean": float(np.median([row["gap_to_euclidean"] for row in corpora])),
        "median_gap_vs_is_bregman": float(np.median([row["gap_to_high_bregman"] for row in corpora])),
    }


def write_csv(rows):
    with (OUTDIR / "english_gap_table.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "slug",
                "corpus",
                "winner_expr",
                "winner_rmse",
                "gap_vs_euclidean",
                "gap_vs_is_bregman",
                "gap_vs_exp_bregman",
                "euclidean_rmse",
                "high_bregman_rmse",
                "low_bregman_rmse",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def build_report(rows, medians, multilang_medians):
    lines = [
        "# English Gap Verification",
        "",
        "| Corpus | Winner | winner RMSE | gap vs Euclidean | gap vs IS-Bregman | gap vs exp-Bregman |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['corpus']} | {row['winner_expr']} | {row['winner_rmse']:.12f} | {row['gap_vs_euclidean']:.12f} | {row['gap_vs_is_bregman']:.12f} | {row['gap_vs_exp_bregman']:.12f} |"
        )

    lines.extend(
        [
            "",
            "## Summary Statistics",
            "",
            f"- median winner-vs-Euclidean gap across 25 English corpora: `{medians['median_gap_vs_euclidean']:.12f}`",
            f"- median winner-vs-IS-Bregman gap across 25 English corpora: `{medians['median_gap_vs_is_bregman']:.12f}`",
            f"- median winner-vs-exp-Bregman gap across 25 English corpora: `{medians['median_gap_vs_exp_bregman']:.12f}`",
            "",
            "## Multilingual Comparison",
            "",
            f"- multilingual 7-corpus median winner-vs-Euclidean gap: `{multilang_medians['median_gap_vs_euclidean']:.12f}`",
            f"- multilingual 7-corpus median winner-vs-IS-Bregman gap: `{multilang_medians['median_gap_vs_is_bregman']:.12f}`",
        ]
    )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = [analyze_corpus(spec) for spec in common.SEARCHED_CORPORA]
    rows.sort(key=lambda row: row["corpus"])
    medians = {
        "median_gap_vs_euclidean": float(np.median([row["gap_vs_euclidean"] for row in rows])),
        "median_gap_vs_is_bregman": float(np.median([row["gap_vs_is_bregman"] for row in rows])),
        "median_gap_vs_exp_bregman": float(np.median([row["gap_vs_exp_bregman"] for row in rows])),
    }
    multilang_medians = load_multilang_medians()
    summary = {
        "rows": rows,
        "medians": medians,
        "multilang_medians": multilang_medians,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(rows)
    (OUTDIR / "report.md").write_text(build_report(rows, medians, multilang_medians), encoding="utf-8")


if __name__ == "__main__":
    main()
