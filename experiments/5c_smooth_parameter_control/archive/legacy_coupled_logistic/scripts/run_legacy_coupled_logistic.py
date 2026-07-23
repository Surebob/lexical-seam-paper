"""Reconstruct legacy coupled-logistic smooth-parameter sweep archive outputs."""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def find_repo_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / "results").is_dir() and (candidate / "experiments").is_dir():
            return candidate
    raise RuntimeError(f"Could not find repository root from {start}")


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = find_repo_root(SCRIPT_PATH)
ARCHIVE_DIR = SCRIPT_PATH.parents[1]
OUTPUT_DIR = ARCHIVE_DIR / "outputs"
SOURCE_BUNDLE = REPO_ROOT / "results" / "zipf_smooth_parameter_sweep"


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None and rows:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames or [])
        writer.writeheader()
        writer.writerows(rows)


def metric(name: str, value: Any, display: str, notes: str) -> dict[str, Any]:
    return {"metric_name": name, "value": value, "display_format": display, "notes": notes}


def main() -> None:
    data = json.loads((SOURCE_BUNDLE / "summary.json").read_text())
    summary = data["summary"]

    rows = []
    for row in data["rows"]:
        rows.append(
            {
                "transition_fraction": row["frac"],
                "sigmoid_width": row["w"],
                "tail_slope_b2": row["b2"],
                "tail_head_slope_contrast": row["delta_b"],
                "transition_k": row["k"],
                "single_zm_c": row["zm_c"],
                "theta_deg": row["theta_deg"],
                "angle_r2": row["angle_r2"],
                "winner_full": row["winner_full"],
                "winner_top100": row["winner_top100"],
                "flip_lambda": row["flip_lambda"],
                "source_model": "legacy_coupled_logistic",
            }
        )

    corr_rows = [
        {
            "parameter": "transition_fraction",
            "corr_theta": summary["corr_theta_frac"],
            "corr_lambda_flip": summary["corr_flip_frac"],
            "theta_regression_weight": summary["theta_regression"]["frac"],
            "flip_regression_weight": summary["flip_regression"]["frac"],
            "notes": "Legacy coupled-logistic sweep parameter.",
        },
        {
            "parameter": "sigmoid_width",
            "corr_theta": summary["corr_theta_w"],
            "corr_lambda_flip": summary["corr_flip_w"],
            "theta_regression_weight": summary["theta_regression"]["w"],
            "flip_regression_weight": summary["flip_regression"]["w"],
            "notes": "Legacy coupled-logistic width; not decoupled-erf w_gate/w_tail.",
        },
        {
            "parameter": "tail_head_slope_contrast",
            "corr_theta": summary["corr_theta_b2"],
            "corr_lambda_flip": summary["corr_flip_b2"],
            "theta_regression_weight": summary["theta_regression"]["b2"],
            "flip_regression_weight": summary["flip_regression"]["b2"],
            "notes": "Implemented as varying b2 with b1 fixed at the low-c median.",
        },
    ]

    aggregate_rows = [
        metric("synthetic_row_count", summary["n_rows"], "integer", "Number of legacy synthetic sweep configurations."),
        metric("synthetic_full_winner_count_is", summary["full_winner_counts"].get("is", 0), "count", "Legacy full-RMSE winner count."),
        metric("synthetic_full_winner_count_xpow", summary["full_winner_counts"].get("xpow", 0), "count", "Legacy full-RMSE winner count."),
        metric("synthetic_full_winner_count_exp", summary["full_winner_counts"].get("exp", 0), "count", "Legacy full-RMSE winner count."),
        metric("synthetic_top100_winner_count_is", summary["top100_winner_counts"].get("is", 0), "count", "Legacy top-100 winner count."),
        metric("synthetic_top100_winner_count_exp", summary["top100_winner_counts"].get("exp", 0), "count", "Legacy top-100 winner count."),
        metric("synthetic_top100_winner_count_xpow", summary["top100_winner_counts"].get("xpow", 0), "count", "Legacy top-100 winner count."),
        metric("theta_range_min", summary["theta_range"][0], "degrees", "Minimum theta in legacy sweep."),
        metric("theta_range_max", summary["theta_range"][1], "degrees", "Maximum theta in legacy sweep."),
        metric("flip_lambda_range_min", summary["flip_range"][0], "lambda", "Minimum flip lambda in legacy sweep."),
        metric("flip_lambda_range_max", summary["flip_range"][1], "lambda", "Maximum flip lambda in legacy sweep."),
    ]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "parameter_sweep_rows.csv", rows)
    write_csv(OUTPUT_DIR / "parameter_sweep_correlations.csv", corr_rows)
    write_csv(OUTPUT_DIR / "aggregate_statistics.csv", aggregate_rows, ["metric_name", "value", "display_format", "notes"])

    provenance = ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for src in [SOURCE_BUNDLE / "summary.json", SOURCE_BUNDLE / "report.md"]:
        shutil.copy2(src, provenance / src.name)

    manifest = {
        "archive_id": "legacy_coupled_logistic",
        "parent_experiment": "5c_smooth_parameter_control",
        "status": "first_class_legacy_prior_canon",
        "model": "coupled_logistic_parameter_sweep",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundle": str(SOURCE_BUNDLE.relative_to(REPO_ROOT)),
        "outputs": {
            "parameter_sweep_rows.csv": {"rows": len(rows)},
            "parameter_sweep_correlations.csv": {"rows": len(corr_rows)},
            "aggregate_statistics.csv": {"rows": len(aggregate_rows)},
        },
    }
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote legacy 5c archive: {len(rows)} sweep rows, {len(corr_rows)} correlation rows")


if __name__ == "__main__":
    main()

