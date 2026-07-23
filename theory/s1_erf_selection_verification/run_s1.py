from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "s1_erf_selection_verification"
WIDENED_SUMMARY_PATH = ROOT / "results" / "zipf_widened_grammar_extended" / "summary.json"
SMOOTH_FIT_PATH = ROOT / "experiments" / "3a_smooth_two_regime_fits" / "outputs" / "smooth_fit_per_corpus.csv"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import zipf_step10_ablation as step10
import zipf_widened_grammar_diagnostic as widened


def sigma_curve(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = np.clip((np.log(ranks) - math.log(k)) / w, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def rank_grid_x(vocab_size: int) -> tuple[np.ndarray, np.ndarray]:
    ranks = np.arange(1, vocab_size + 1, dtype=np.float64)
    x = 0.05 + 0.95 * np.log(ranks) / math.log(vocab_size)
    return ranks, x


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(a, b) / denom)


def classify_expr(expr: str) -> str:
    expr_lower = expr.lower()
    if (
        "erf[" in expr_lower
        or "tanh[" in expr_lower
        or "logistic" in expr_lower
        or "div[1,add[1,exp[" in expr_lower
        or "div[1,add[exp[" in expr_lower
    ):
        return "sigmoid-family"
    if expr in {"sub[sub[x,1],log[x]]", "eml[sub[x,1],eml[x,1]]"}:
        return "bregman-family"
    if expr in {"mul[sub[1,x],sub[1,x]]", "mul[sub[1,x],sub[x,1]]"}:
        return "polynomial"
    return "other"


def load_smooth_rows() -> dict[str, dict[str, str]]:
    with SMOOTH_FIT_PATH.open(newline="", encoding="utf-8") as handle:
        return {row["slug"]: row for row in csv.DictReader(handle)}


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    widened.install_extended_grammar()

    summary = json.loads(WIDENED_SUMMARY_PATH.read_text(encoding="utf-8"))
    smooth_rows = load_smooth_rows()
    english_rows = [row for row in summary["rows"] if row["language"] == "English"]

    top1_rows: list[dict[str, str]] = []
    top10_rows: list[dict[str, str]] = []
    top10_counter: Counter[str] = Counter()
    top1_counter: Counter[str] = Counter()
    cosines: list[float] = []
    abs_cosines: list[float] = []

    for row in english_rows:
        slug = row["slug"]
        smooth = smooth_rows[slug]
        vocab_size = int(float(smooth["vocab_size"]))
        k = float(smooth["k"])
        w = float(smooth["w"])
        ranks, x_values = rank_grid_x(vocab_size)
        sigma = sigma_curve(ranks, k, w)

        for item in row["top10"]:
            expr = item["expr"]
            family = classify_expr(expr)
            top10_counter[family] += 1
            top10_rows.append(
                {
                    "slug": slug,
                    "corpus": row["name"],
                    "beam_rank": str(item["rank"]),
                    "expr": expr,
                    "math": item["math"],
                    "rmse": repr(float(item["rmse"])),
                    "uses_new_operator": str(bool(item["uses_new_operator"])),
                    "symbolic_family": family,
                }
            )
            if int(item["rank"]) == 1:
                values = step10.eval_node(step10.parse_expr(expr), x_values)
                cosine = cosine_similarity(sigma, values)
                abs_cosine = abs(cosine)
                cosines.append(cosine)
                abs_cosines.append(abs_cosine)
                top1_counter[family] += 1
                top1_rows.append(
                    {
                        "slug": slug,
                        "corpus": row["name"],
                        "zm_c": repr(float(row["zm_c"])),
                        "vocab_size": str(vocab_size),
                        "k": repr(k),
                        "w": repr(w),
                        "top1_expr": expr,
                        "top1_math": item["math"],
                        "top1_symbolic_family": family,
                        "top1_rmse": repr(float(item["rmse"])),
                        "sigma_winner_cosine": repr(cosine),
                        "sigma_winner_abs_cosine": repr(abs_cosine),
                    }
                )

    family_rows = []
    for scope, counter in [("top1_only", top1_counter), ("top10_all_entries", top10_counter)]:
        total = sum(counter.values())
        for family, count in sorted(counter.items()):
            family_rows.append(
                {
                    "scope": scope,
                    "symbolic_family": family,
                    "count": str(count),
                    "proportion": repr(count / total if total else float("nan")),
                }
            )

    summary_rows = [
        {
            "metric_name": "english_low_c_corpus_count",
            "value": str(len(english_rows)),
            "notes": "English low-c corpora in results/zipf_widened_grammar_extended.",
        },
        {
            "metric_name": "sigma_winner_cosine_min",
            "value": repr(float(np.min(cosines))),
            "notes": "Minimum raw cosine between smooth sigma_r and widened top-1 winner.",
        },
        {
            "metric_name": "sigma_winner_cosine_max",
            "value": repr(float(np.max(cosines))),
            "notes": "Maximum raw cosine between smooth sigma_r and widened top-1 winner.",
        },
        {
            "metric_name": "sigma_winner_cosine_mean",
            "value": repr(float(np.mean(cosines))),
            "notes": "Mean raw cosine between smooth sigma_r and widened top-1 winner.",
        },
        {
            "metric_name": "sigma_winner_cosine_median",
            "value": repr(float(np.median(cosines))),
            "notes": "Median raw cosine between smooth sigma_r and widened top-1 winner.",
        },
        {
            "metric_name": "sigma_winner_abs_cosine_mean",
            "value": repr(float(np.mean(abs_cosines))),
            "notes": "Mean absolute cosine between smooth sigma_r and widened top-1 winner.",
        },
        {
            "metric_name": "sigma_winner_abs_cosine_median",
            "value": repr(float(np.median(abs_cosines))),
            "notes": "Median absolute cosine between smooth sigma_r and widened top-1 winner.",
        },
    ]

    with (OUTDIR / "s1_top1_sigma_cosine_per_corpus.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "slug",
                "corpus",
                "zm_c",
                "vocab_size",
                "k",
                "w",
                "top1_expr",
                "top1_math",
                "top1_symbolic_family",
                "top1_rmse",
                "sigma_winner_cosine",
                "sigma_winner_abs_cosine",
            ],
        )
        writer.writeheader()
        writer.writerows(top1_rows)

    with (OUTDIR / "s1_top10_family_classification.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "slug",
                "corpus",
                "beam_rank",
                "expr",
                "math",
                "rmse",
                "uses_new_operator",
                "symbolic_family",
            ],
        )
        writer.writeheader()
        writer.writerows(top10_rows)

    with (OUTDIR / "s1_family_distribution.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scope", "symbolic_family", "count", "proportion"])
        writer.writeheader()
        writer.writerows(family_rows)

    with (OUTDIR / "s1_summary_statistics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric_name", "value", "notes"])
        writer.writeheader()
        writer.writerows(summary_rows)

    report_lines = [
        "# S1 erf Selection Verification",
        "",
        f"- English low-c corpora analyzed: `{len(english_rows)}`",
        f"- raw cosine min/max/mean/median: `{np.min(cosines):.12f}`, `{np.max(cosines):.12f}`, `{np.mean(cosines):.12f}`, `{np.median(cosines):.12f}`",
        f"- absolute cosine mean/median: `{np.mean(abs_cosines):.12f}`, `{np.median(abs_cosines):.12f}`",
        "",
        "## Top-1 symbolic-family counts",
        "",
    ]
    for family, count in sorted(top1_counter.items()):
        report_lines.append(f"- `{family}`: `{count}`")
    report_lines.extend(["", "## Top-10 symbolic-family counts", ""])
    for family, count in sorted(top10_counter.items()):
        report_lines.append(f"- `{family}`: `{count}`")
    (OUTDIR / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
