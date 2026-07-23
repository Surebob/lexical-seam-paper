#!/usr/bin/env python3
"""Build J2 scaling comparison from canonical decoupled-erf k scaling and POS scaling."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import sys

from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
JOIN_DIR = SCRIPT_DIR.parents[0]
ROOT = SCRIPT_DIR.parents[2]
sys.path.insert(0, str(JOIN_DIR))

from source_config import OUTPUTS, UPSTREAMS  # noqa: E402


def read_metric_csv(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing upstream aggregate: {path}")
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if "metric_name" not in (reader.fieldnames or []):
            raise ValueError(f"Expected metric_name column in {path}")
        return {row["metric_name"]: row for row in reader if row.get("metric_name")}


def read_table_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing upstream table: {path}")
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def float_value(row: dict[str, str], key: str = "value") -> float:
    return float(row[key])


def compute_erf_k_scaling(per_fit_path: Path, per_corpus_path: Path) -> dict[str, float]:
    per_fit = read_table_csv(per_fit_path)
    per_corpus = read_table_csv(per_corpus_path)

    vocab_by_slug = {row["slug"]: float(row["vocab_size"]) for row in per_corpus}
    points: list[tuple[float, float]] = []
    for row in per_fit:
        if row.get("gate_family") != "erf":
            continue
        slug = row["slug"]
        if slug not in vocab_by_slug:
            raise ValueError(f"Missing per-corpus vocab_size for {slug}")
        points.append((math.log(vocab_by_slug[slug]), math.log(float(row["k"]))))

    if len(points) < 3:
        raise ValueError(f"Need at least 3 erf k-scaling points; found {len(points)}")

    n = len(points)
    xs = [x for x, _ in points]
    ys = [y for _, y in points]
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    sxy = sum((x - xbar) * (y - ybar) for x, y in points)
    beta = sxy / sxx
    alpha = ybar - beta * xbar
    residuals = [y - (alpha + beta * x) for x, y in points]
    rss = sum(r * r for r in residuals)
    tss = sum((y - ybar) ** 2 for y in ys)
    df = float(n - 2)
    beta_se = math.sqrt((rss / df) / sxx)
    tcrit = stats.t.ppf(0.975, df)
    return {
        "n": float(n),
        "alpha": alpha,
        "beta": beta,
        "beta_se": beta_se,
        "beta_ci_low": beta - tcrit * beta_se,
        "beta_ci_high": beta + tcrit * beta_se,
        "df": df,
        "r2": 1.0 - rss / tss,
        "rss": rss,
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["metric_name", "value", "display_format", "notes"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS["summary"].parent)
    args = parser.parse_args()

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = (out_dir / "scaling_comparison_summary.csv").resolve()
    manifest_path = (out_dir / "manifest.json").resolve()

    erf_scaling = compute_erf_k_scaling(UPSTREAMS["3e_erf_fits"], UPSTREAMS["3e_per_corpus"])
    pos = read_metric_csv(UPSTREAMS["3d_pos"])

    beta = erf_scaling["beta"]
    beta_se = erf_scaling["beta_se"]
    beta_df = erf_scaling["df"]

    pos_alpha = float_value(pos["forced_alpha"])
    pos_ci_low = float_value(pos["forced_alpha_ci_low"])
    pos_ci_high = float_value(pos["forced_alpha_ci_high"])

    # Recover the POS standard error from the saved t interval. The 25-corpus
    # forced-alpha fit uses df = n - 1 = 24 in experiment 3d.
    pos_df = 24.0
    pos_tcrit = stats.t.ppf(0.975, pos_df)
    pos_se = (pos_ci_high - pos_ci_low) / (2.0 * pos_tcrit)

    diff = pos_alpha - beta
    diff_se = math.sqrt(beta_se**2 + pos_se**2)
    welch_num = diff_se**4
    welch_den = (beta_se**4 / beta_df) + (pos_se**4 / pos_df)
    welch_df = welch_num / welch_den
    tcrit = stats.t.ppf(0.975, welch_df)
    ci_low = diff - tcrit * diff_se
    ci_high = diff + tcrit * diff_se
    pvalue = 2.0 * stats.t.sf(abs(diff / diff_se), welch_df)

    rows = [
        {
            "metric_name": "alpha_kstat",
            "value": repr(beta),
            "display_format": "decimal_3",
            "notes": "Canonical decoupled-erf OLS beta recomputed from experiment 3e erf k values.",
        },
        {
            "metric_name": "alpha_pos",
            "value": repr(pos_alpha),
            "display_format": "decimal_3",
            "notes": "POS crossover forced alpha from experiment 3d.",
        },
        {
            "metric_name": "alpha_difference_pos_minus_kstat",
            "value": repr(diff),
            "display_format": "signed_decimal_3",
            "notes": "POS alpha minus decoupled-erf k-stat beta.",
        },
        {
            "metric_name": "alpha_difference_ci_low",
            "value": repr(float(ci_low)),
            "display_format": "signed_decimal_3",
            "notes": "Lower 95% Welch CI for POS alpha minus decoupled-erf beta.",
        },
        {
            "metric_name": "alpha_difference_ci_high",
            "value": repr(float(ci_high)),
            "display_format": "signed_decimal_3",
            "notes": "Upper 95% Welch CI for POS alpha minus decoupled-erf beta.",
        },
        {
            "metric_name": "alpha_difference_pvalue",
            "value": repr(float(pvalue)),
            "display_format": "decimal_3",
            "notes": "Two-sided Welch t-test for difference between POS alpha and decoupled-erf beta.",
        },
    ]

    write_csv(summary_path, rows)

    manifest = {
        "join_id": "J2_scaling_comparison",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve().relative_to(ROOT)),
        "outputs": {
            "scaling_comparison_summary.csv": {
                "path": str(summary_path.relative_to(ROOT)),
                "rows": len(rows),
                "schema": ["metric_name", "value", "display_format", "notes"],
            }
        },
        "upstreams": {name: str(path.relative_to(ROOT)) for name, path in UPSTREAMS.items()},
        "method": {
            "kstat_source": "OLS beta from experiment 3e erf k values joined to experiment 3e vocab sizes.",
            "kstat_ols": {
                "n": erf_scaling["n"],
                "intercept_alpha": erf_scaling["alpha"],
                "beta": erf_scaling["beta"],
                "beta_se": erf_scaling["beta_se"],
                "beta_ci_low": erf_scaling["beta_ci_low"],
                "beta_ci_high": erf_scaling["beta_ci_high"],
                "df": erf_scaling["df"],
                "r2": erf_scaling["r2"],
                "rss": erf_scaling["rss"],
            },
            "ci": "Welch-propagated 95% Student-t CI using decoupled beta SE and POS alpha SE.",
            "p_value": "Two-sided Welch t-test.",
            "welch_df": welch_df,
            "pos_alpha_se_from_ci": pos_se,
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
