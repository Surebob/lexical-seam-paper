from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t, ttest_1samp

SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
ROOT = EXPERIMENT_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

import source_config as cfg


def forced_scaling_model(x: np.ndarray, alpha: float) -> np.ndarray:
    return alpha * x


def load_points() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with cfg.SMOOTH_FIT_CSV.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            vocab_size = float(row["vocab_size"])
            k_stat = float(row["k"])
            rows.append(
                {
                    "slug": row["slug"],
                    "corpus": row["corpus"],
                    "vocabulary_size": vocab_size,
                    "k_stat": k_stat,
                    "log_vocabulary_size": float(math.log(vocab_size)),
                    "log_k_stat": float(math.log(k_stat)),
                    "alpha_stat_per_corpus": float(math.log(k_stat) / math.log(vocab_size)),
                }
            )
    return sorted(rows, key=lambda item: str(item["slug"]))


def fit_forced_alpha(rows: list[dict[str, float | str]]) -> dict[str, float]:
    x = np.array([float(row["log_vocabulary_size"]) for row in rows], dtype=float)
    y = np.array([float(row["log_k_stat"]) for row in rows], dtype=float)
    popt, pcov = curve_fit(forced_scaling_model, x, y, p0=(cfg.HALF_REFERENCE,), maxfev=10000)
    alpha = float(popt[0])
    alpha_se = float(math.sqrt(pcov[0, 0]))
    df = len(rows) - 1
    tcrit = float(t.ppf(0.975, df))
    return {
        "forced_alpha": alpha,
        "forced_alpha_se": alpha_se,
        "forced_alpha_ci_low": alpha - tcrit * alpha_se,
        "forced_alpha_ci_high": alpha + tcrit * alpha_se,
        "forced_alpha_df": float(df),
    }


def summarize_distribution(rows: list[dict[str, float | str]]) -> dict[str, float]:
    values = np.array([float(row["alpha_stat_per_corpus"]) for row in rows], dtype=float)
    mean_alpha = float(np.mean(values))
    median_alpha = float(np.median(values))
    se = float(np.std(values, ddof=1) / math.sqrt(len(values)))
    df = len(values) - 1
    tcrit = float(t.ppf(0.975, df))
    test = ttest_1samp(values, popmean=cfg.HALF_REFERENCE)
    return {
        "mean_alpha_per_corpus": mean_alpha,
        "median_alpha_per_corpus": median_alpha,
        "mean_alpha_ci_low": mean_alpha - tcrit * se,
        "mean_alpha_ci_high": mean_alpha + tcrit * se,
        "alpha_stat_vs_half_t_statistic": float(test.statistic),
        "alpha_stat_vs_half_pvalue": float(test.pvalue),
        "forced_alpha_vs_half_pvalue": float(test.pvalue),
        "alpha_stat_vs_half_df": float(df),
    }


def build_aggregate_rows(rows: list[dict[str, float | str]]) -> list[dict[str, str]]:
    forced = fit_forced_alpha(rows)
    dist = summarize_distribution(rows)
    aggregate_rows = [
        {
            "metric_name": "english_corpus_count",
            "value": str(len(rows)),
            "display_format": "integer",
            "notes": "Canonical English corpus count inherited from experiment 3a.",
        },
        {
            "metric_name": "forced_alpha",
            "value": repr(forced["forced_alpha"]),
            "display_format": "decimal_4",
            "notes": "Forced fit of log(k_stat) = alpha * log(V).",
        },
        {
            "metric_name": "forced_alpha_se",
            "value": repr(forced["forced_alpha_se"]),
            "display_format": "decimal_6",
            "notes": "SE from curve_fit covariance.",
        },
        {
            "metric_name": "forced_alpha_ci_low",
            "value": repr(forced["forced_alpha_ci_low"]),
            "display_format": "decimal_4",
            "notes": cfg.CI_METHOD,
        },
        {
            "metric_name": "forced_alpha_ci_high",
            "value": repr(forced["forced_alpha_ci_high"]),
            "display_format": "decimal_4",
            "notes": cfg.CI_METHOD,
        },
        {
            "metric_name": "forced_alpha_df",
            "value": str(int(forced["forced_alpha_df"])),
            "display_format": "integer",
            "notes": "Degrees of freedom used for the forced-fit alpha CI.",
        },
        {
            "metric_name": "mean_alpha_per_corpus",
            "value": repr(dist["mean_alpha_per_corpus"]),
            "display_format": "decimal_6",
            "notes": "Mean of per-corpus alpha_stat = log(k_stat)/log(V).",
        },
        {
            "metric_name": "median_alpha_per_corpus",
            "value": repr(dist["median_alpha_per_corpus"]),
            "display_format": "decimal_6",
            "notes": "Median of per-corpus alpha_stat = log(k_stat)/log(V).",
        },
        {
            "metric_name": "mean_alpha_ci_low",
            "value": repr(dist["mean_alpha_ci_low"]),
            "display_format": "decimal_6",
            "notes": "Student-t 95% CI for the mean per-corpus alpha_stat value.",
        },
        {
            "metric_name": "mean_alpha_ci_high",
            "value": repr(dist["mean_alpha_ci_high"]),
            "display_format": "decimal_6",
            "notes": "Student-t 95% CI for the mean per-corpus alpha_stat value.",
        },
        {
            "metric_name": "alpha_stat_vs_half_t_statistic",
            "value": repr(dist["alpha_stat_vs_half_t_statistic"]),
            "display_format": "decimal_6",
            "notes": f"Two-sided one-sample t statistic for alpha_stat vs {cfg.HALF_REFERENCE}.",
        },
        {
            "metric_name": "alpha_stat_vs_half_pvalue",
            "value": repr(dist["alpha_stat_vs_half_pvalue"]),
            "display_format": "decimal_6",
            "notes": f"Two-sided one-sample p-value for mean(alpha_stat) vs {cfg.HALF_REFERENCE}.",
        },
        {
            "metric_name": "forced_alpha_vs_half_pvalue",
            "value": repr(dist["forced_alpha_vs_half_pvalue"]),
            "display_format": "decimal_6",
            "notes": "Compatibility row for the manuscript claim map; numerically identical to alpha_stat_vs_half_pvalue and based on the per-corpus one-sample t-test, not a separate test on the forced-fit alpha.",
        },
        {
            "metric_name": "alpha_stat_vs_half_df",
            "value": str(int(dist["alpha_stat_vs_half_df"])),
            "display_format": "integer",
            "notes": "Degrees of freedom for the one-sample t-test on alpha_stat values.",
        },
    ]
    return aggregate_rows


def write_points_csv(rows: list[dict[str, float | str]], path: Path) -> None:
    fieldnames = [
        "slug",
        "corpus",
        "vocabulary_size",
        "k_stat",
        "log_vocabulary_size",
        "log_k_stat",
        "alpha_stat_per_corpus",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_aggregate_csv(rows: list[dict[str, str]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric_name", "value", "display_format", "notes"])
        writer.writeheader()
        writer.writerows(rows)


def load_historical_summary() -> dict:
    return json.loads(cfg.HISTORICAL_SUMMARY_JSON.read_text(encoding="utf-8"))


def compare_points(
    rows: list[dict[str, float | str]],
    historical_summary: dict,
    path: Path,
) -> tuple[float, int]:
    historical_rows = {row["slug"]: row for row in historical_summary["rows"]}
    field_map = {
        "corpus": "name",
        "vocabulary_size": "V",
        "k_stat": "k_stat",
        "log_vocabulary_size": "log_V",
        "log_k_stat": "log_k_stat",
        "alpha_stat_per_corpus": "alpha_stat_per_corpus",
    }
    max_abs = 0.0
    text_mismatches = 0
    diff_rows: list[dict[str, str]] = []
    for row in rows:
        slug = str(row["slug"])
        historical = historical_rows[slug]
        for new_field, historical_field in field_map.items():
            new_value = row[new_field]
            historical_value = historical[historical_field]
            diff_value = ""
            matches = False
            if isinstance(new_value, str):
                matches = str(new_value) == str(historical_value)
                if not matches:
                    text_mismatches += 1
            else:
                numeric_diff = abs(float(new_value) - float(historical_value))
                diff_value = repr(numeric_diff)
                max_abs = max(max_abs, numeric_diff)
                matches = numeric_diff <= 1e-12
            diff_rows.append(
                {
                    "slug": slug,
                    "field": new_field,
                    "new_value": str(new_value),
                    "historical_value": str(historical_value),
                    "abs_diff": diff_value,
                    "matches": str(matches),
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["slug", "field", "new_value", "historical_value", "abs_diff", "matches"],
        )
        writer.writeheader()
        writer.writerows(diff_rows)
    return max_abs, text_mismatches


def compare_aggregates(
    aggregate_rows: list[dict[str, str]],
    historical_summary: dict,
    path: Path,
) -> tuple[float, int]:
    historical_map = {
        "forced_alpha": historical_summary["kstat_fit"]["alpha"],
        "forced_alpha_se": historical_summary["kstat_fit"]["alpha_se"],
        "forced_alpha_ci_low": historical_summary["kstat_fit"]["alpha_ci_95"][0],
        "forced_alpha_ci_high": historical_summary["kstat_fit"]["alpha_ci_95"][1],
        "forced_alpha_df": historical_summary["kstat_fit"]["df"],
        "mean_alpha_per_corpus": historical_summary["kstat_distribution"]["mean_alpha"],
        "median_alpha_per_corpus": historical_summary["kstat_distribution"]["median_alpha"],
        "mean_alpha_ci_low": historical_summary["kstat_distribution"]["mean_alpha_ci_95"][0],
        "mean_alpha_ci_high": historical_summary["kstat_distribution"]["mean_alpha_ci_95"][1],
        "alpha_stat_vs_half_t_statistic": historical_summary["kstat_distribution"]["t_statistic_vs_0_5"],
        "alpha_stat_vs_half_pvalue": historical_summary["kstat_distribution"]["p_value_vs_0_5"],
        "forced_alpha_vs_half_pvalue": historical_summary["kstat_distribution"]["p_value_vs_0_5"],
        "alpha_stat_vs_half_df": historical_summary["kstat_distribution"]["df"],
    }
    max_abs = 0.0
    text_mismatches = 0
    diff_rows: list[dict[str, str]] = []
    for row in aggregate_rows:
        metric = row["metric_name"]
        if metric not in historical_map:
            continue
        new_value = row["value"]
        historical_value = historical_map[metric]
        diff_value = ""
        matches = False
        try:
            numeric_diff = abs(float(new_value) - float(historical_value))
            diff_value = repr(numeric_diff)
            max_abs = max(max_abs, numeric_diff)
            matches = numeric_diff <= 1e-12
        except ValueError:
            matches = new_value == str(historical_value)
            if not matches:
                text_mismatches += 1
        diff_rows.append(
            {
                "metric_name": metric,
                "new_value": new_value,
                "historical_value": str(historical_value),
                "abs_diff": diff_value,
                "matches": str(matches),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["metric_name", "new_value", "historical_value", "abs_diff", "matches"],
        )
        writer.writeheader()
        writer.writerows(diff_rows)
    return max_abs, text_mismatches


def write_diff_summary(
    point_count: int,
    point_max_abs: float,
    point_text_mismatches: int,
    aggregate_count: int,
    aggregate_max_abs: float,
    aggregate_text_mismatches: int,
    path: Path,
) -> None:
    rows = [
        {
            "metric_name": "kstat_scaling_point_row_count",
            "value": str(point_count),
            "notes": "Expected 25 per-corpus k_stat scaling rows.",
        },
        {
            "metric_name": "point_max_abs_numeric_diff",
            "value": repr(point_max_abs),
            "notes": "Maximum absolute numeric diff across overlapping per-corpus fields.",
        },
        {
            "metric_name": "point_text_mismatch_count",
            "value": str(point_text_mismatches),
            "notes": "Count of non-numeric mismatches against the historical k_stat summary rows.",
        },
        {
            "metric_name": "aggregate_row_count",
            "value": str(aggregate_count),
            "notes": "Count of overlapping aggregate rows compared to the historical bundle.",
        },
        {
            "metric_name": "aggregate_max_abs_numeric_diff",
            "value": repr(aggregate_max_abs),
            "notes": "Maximum absolute numeric diff across overlapping aggregate metrics.",
        },
        {
            "metric_name": "aggregate_text_mismatch_count",
            "value": str(aggregate_text_mismatches),
            "notes": "Count of non-numeric mismatches in aggregate metrics.",
        },
        {
            "metric_name": "deferred_historical_sections",
            "value": "2",
            "notes": "Historical POS forced-fit reference and historical k_stat-vs-POS comparison are intentionally deferred to experiments 3d and J2 under the consolidated design.",
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric_name", "value", "notes"])
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(output_dir: Path, point_rows: list[dict[str, float | str]], aggregate_rows: list[dict[str, str]]) -> None:
    manifest = {
        "experiment_id": "3c_kstat_scaling",
        "research_question": "How does the transition centre inferred from the canonical smooth fit scale with vocabulary size across the 25 English corpora?",
        "upstream_dependencies": [
            str(cfg.SMOOTH_FIT_CSV),
            str(cfg.HISTORICAL_SUMMARY_JSON),
            str(cfg.HISTORICAL_REPORT_MD),
        ],
        "output_files": [
            "kstat_scaling_points.csv",
            "aggregate_statistics.csv",
            "historical_point_diff.csv",
            "historical_aggregate_diff.csv",
            "historical_diff_summary.csv",
        ],
        "protocol": {
            "forced_model": cfg.FORCED_SCALING_MODEL,
            "half_reference": cfg.HALF_REFERENCE,
            "ci_method": cfg.CI_METHOD,
            "pvalue_method": "Two-sided one-sample Student-t test on per-corpus alpha_stat values against 0.5",
        },
        "row_counts": {
            "kstat_scaling_points": len(point_rows),
            "aggregate_statistics": len(aggregate_rows),
        },
        "deferred_to_j2": [
            "historical POS forced-fit reference",
            "historical k_stat-vs-k_POS comparison block",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=cfg.OUTPUT_DIR,
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    point_rows = load_points()
    aggregate_rows = build_aggregate_rows(point_rows)
    historical_summary = load_historical_summary()

    write_points_csv(point_rows, output_dir / "kstat_scaling_points.csv")
    write_aggregate_csv(aggregate_rows, output_dir / "aggregate_statistics.csv")

    point_max_abs, point_text_mismatches = compare_points(
        point_rows,
        historical_summary,
        output_dir / "historical_point_diff.csv",
    )
    aggregate_max_abs, aggregate_text_mismatches = compare_aggregates(
        aggregate_rows,
        historical_summary,
        output_dir / "historical_aggregate_diff.csv",
    )
    compared_aggregate_count = sum(
        1
        for row in aggregate_rows
        if row["metric_name"]
        in {
            "forced_alpha",
            "forced_alpha_se",
            "forced_alpha_ci_low",
            "forced_alpha_ci_high",
            "forced_alpha_df",
            "mean_alpha_per_corpus",
            "median_alpha_per_corpus",
            "mean_alpha_ci_low",
            "mean_alpha_ci_high",
            "alpha_stat_vs_half_t_statistic",
            "alpha_stat_vs_half_pvalue",
            "forced_alpha_vs_half_pvalue",
            "alpha_stat_vs_half_df",
        }
    )
    write_diff_summary(
        point_count=len(point_rows),
        point_max_abs=point_max_abs,
        point_text_mismatches=point_text_mismatches,
        aggregate_count=compared_aggregate_count,
        aggregate_max_abs=aggregate_max_abs,
        aggregate_text_mismatches=aggregate_text_mismatches,
        path=output_dir / "historical_diff_summary.csv",
    )
    write_manifest(output_dir, point_rows, aggregate_rows)


if __name__ == "__main__":
    main()
