from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

import source_config as cfg


GAP_FIELDS = [
    "slug",
    "corpus",
    "winner_expr",
    "winner_rmse",
    "euclidean_rmse",
    "is_bregman_rmse",
    "exp_bregman_rmse",
    "gap_vs_euclidean",
    "gap_vs_is_bregman",
    "gap_vs_exp_bregman",
    "winner_is_euclidean",
    "max_gap_to_core_bregman_family",
    "core_family_gap_below_0p01",
]

AGG_FIELDS = ["metric_name", "value", "display_format", "notes"]


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def metric(name: str, value: object, fmt: str, notes: str) -> dict[str, object]:
    return {"metric_name": name, "value": value, "display_format": fmt, "notes": notes}


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(cfg.GAP_SUMMARY.read_text(encoding="utf-8"))
    rows: list[dict[str, object]] = []
    euclidean_exprs = {"mul[sub[1,x],sub[1,x]]", "mul[sub[x,1],sub[x,1]]"}

    for item in summary["rows"]:
        core_gap = max(abs(float(item["gap_vs_is_bregman"])), abs(float(item["gap_vs_exp_bregman"])))
        winner_is_euclidean = item["winner_expr"] in euclidean_exprs
        rows.append(
            {
                "slug": item["slug"],
                "corpus": item["corpus"],
                "winner_expr": item["winner_expr"],
                "winner_rmse": item["winner_rmse"],
                "euclidean_rmse": item["euclidean_rmse"],
                "is_bregman_rmse": item["high_bregman_rmse"],
                "exp_bregman_rmse": item["low_bregman_rmse"],
                "gap_vs_euclidean": item["gap_vs_euclidean"],
                "gap_vs_is_bregman": item["gap_vs_is_bregman"],
                "gap_vs_exp_bregman": item["gap_vs_exp_bregman"],
                "winner_is_euclidean": str(winner_is_euclidean),
                "max_gap_to_core_bregman_family": core_gap,
                "core_family_gap_below_0p01": str(core_gap < 0.01),
            }
        )

    euclidean_winner_count = sum(row["winner_is_euclidean"] == "True" for row in rows)
    max_core_gap = max(float(row["max_gap_to_core_bregman_family"]) for row in rows)
    below_count = sum(row["core_family_gap_below_0p01"] == "True" for row in rows)

    aggregate_rows = [
        metric(
            "median_winner_minus_euclidean_gap_full_rmse",
            summary["medians"]["median_gap_vs_euclidean"],
            "float",
            "Median full-RMSE gap between the winning step-2 expression and Euclidean Bregman candidate.",
        ),
        metric(
            "euclidean_step2_winner_count",
            euclidean_winner_count,
            "X/25 count",
            "Count of corpora where the Euclidean Bregman candidate is the step-2 winner.",
        ),
        metric(
            "max_core_generator_gap_full_rmse_typical_band",
            max_core_gap,
            "float rendered as below 0.01",
            "Maximum full-RMSE gap between the winner and the better of the IS/exp Bregman core family.",
        ),
        metric(
            "winner_family_gap_below_0p01_count",
            below_count,
            "X/25 count",
            "Count of corpora where the winner is within 0.01 full-RMSE of the core IS/exp Bregman family.",
        ),
    ]

    write_csv(output_dir / "gap_analysis_per_corpus.csv", rows, GAP_FIELDS)
    write_csv(output_dir / "aggregate_statistics.csv", aggregate_rows, AGG_FIELDS)
    manifest = {
        "experiment_id": "1d_winner_vs_euclidean_gap_analysis",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(path.relative_to(cfg.ROOT)) for path in cfg.SOURCE_BUNDLES],
        "outputs": {
            "gap_analysis_per_corpus.csv": {"rows": len(rows), "schema": GAP_FIELDS},
            "aggregate_statistics.csv": {"rows": len(aggregate_rows), "schema": AGG_FIELDS},
        },
        "notes": ["No symbolic search was rerun; values are normalized from zipf_english_gap_verify."],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps({"source_bundles": manifest["source_bundles"]}, indent=2),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()

