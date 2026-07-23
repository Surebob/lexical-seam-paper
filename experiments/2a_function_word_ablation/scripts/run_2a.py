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


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(cfg.FUNCTION_WORD_SUMMARY_JSON.read_text(encoding="utf-8"))
    rows: list[dict[str, object]] = []
    for item in summary["results"]:
        rows.append(
            {
                "ablation": cfg.ABLATION_LABELS[item["case"]],
                "historical_case": item["case"],
                "token_count": item["token_count"],
                "unique_words": item["unique_words"],
                "zm_a": item["zm_a"],
                "zm_b": item["zm_b"],
                "zm_c": item["zm_c"],
                "zm_rmse": item["zm_rmse"],
                "step2_winner": item["step2_winner"],
                "step2_winner_display": item["step2_math"],
                "step2_rmse": item["step2_rmse"],
                "step2_helpful": str(item["step2_helps"]),
                "rmse_delta_step2_minus_zm": float(item["step2_rmse"]) - float(item["zm_rmse"]),
            }
        )
    rows.sort(key=lambda row: ["remove_top_100", "remove_top_50", "top_100_only"].index(str(row["ablation"])))
    fields = [
        "ablation",
        "historical_case",
        "token_count",
        "unique_words",
        "zm_a",
        "zm_b",
        "zm_c",
        "zm_rmse",
        "step2_winner",
        "step2_winner_display",
        "step2_rmse",
        "step2_helpful",
        "rmse_delta_step2_minus_zm",
    ]
    write_csv(output_dir / "function_word_ablation.csv", rows, fields)
    write_csv(output_dir / "aggregate_statistics.csv", [], ["metric_name", "value", "display_format", "notes"])
    manifest = {
        "experiment_id": "2a_function_word_ablation",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(cfg.FUNCTION_WORD_DIR.relative_to(cfg.ROOT))],
        "outputs": {
            "function_word_ablation.csv": {"rows": len(rows), "schema": fields},
            "aggregate_statistics.csv": {"rows": 0, "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "notes": ["No claim-map aggregate rows are assigned to 2a; claims map directly to function_word_ablation.csv cells."],
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
