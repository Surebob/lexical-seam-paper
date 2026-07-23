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


def find_candidate(step_summary: list[dict], expr: str) -> dict[str, object] | None:
    for step in step_summary:
        for rank, candidate in enumerate(step.get("top_candidates", []), start=1):
            if candidate.get("expr") == expr:
                return {"rank": rank, **candidate}
    return None


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    step2 = json.loads(cfg.CITY_STEP2_SUMMARY.read_text(encoding="utf-8"))
    full = json.loads(cfg.CITY_FULL_SUMMARY.read_text(encoding="utf-8"))

    corpus = step2["corpus"]
    zm = step2["zm_baseline"]
    zm_search = step2["zm_search"]
    full_zm_search = full["zm_search"]
    exp_expr = "eml[sub[x,1],eml[x,1]]"
    is_expr = "sub[sub[x,1],log[x]]"
    exp_candidate = find_candidate(zm_search["step_summary"], exp_expr)
    is_candidate = find_candidate(zm_search["step_summary"], is_expr)

    row = {
        "dataset": "world_city_populations",
        "city_count": corpus["unique_words"],
        "vocab_size": corpus["unique_words"],
        "population_sum": corpus["token_count"],
        "zm_a": zm["a"],
        "zm_b": zm["b"],
        "zm_c": zm["c"],
        "zm_rmse": zm["rmse_full"],
        "step2_winner_expr": zm_search["best"]["expr"],
        "step2_winner_math": zm_search["best"]["math"],
        "step2_composite_rmse": zm_search["best"]["composite_rmse"],
        "step2_helpful": str(float(zm_search["best"]["composite_rmse"]) < float(zm["rmse_full"])),
        "step2_delta_vs_zm": float(zm_search["best"]["composite_rmse"]) - float(zm["rmse_full"]),
        "is_bregman_rank_in_step2_beam": is_candidate["rank"] if is_candidate else "",
        "is_bregman_rmse": is_candidate["rmse"] if is_candidate else "",
        "exp_bregman_rank_in_step2_beam": exp_candidate["rank"] if exp_candidate else "",
        "exp_bregman_rmse": exp_candidate["rmse"] if exp_candidate else "",
        "terminal_winner_expr": full_zm_search["best"]["expr"],
        "terminal_winner_math": full_zm_search["best"]["math"],
        "terminal_composite_rmse": full_zm_search["best"]["composite_rmse"],
        "terminal_helpful": str(float(full_zm_search["best"]["composite_rmse"]) < float(zm["rmse_full"])),
    }
    fields = [
        "dataset",
        "city_count",
        "vocab_size",
        "population_sum",
        "zm_a",
        "zm_b",
        "zm_c",
        "zm_rmse",
        "step2_winner_expr",
        "step2_winner_math",
        "step2_composite_rmse",
        "step2_helpful",
        "step2_delta_vs_zm",
        "is_bregman_rank_in_step2_beam",
        "is_bregman_rmse",
        "exp_bregman_rank_in_step2_beam",
        "exp_bregman_rmse",
        "terminal_winner_expr",
        "terminal_winner_math",
        "terminal_composite_rmse",
        "terminal_helpful",
    ]
    write_csv(output_dir / "city_population_control.csv", [row], fields)
    write_csv(output_dir / "aggregate_statistics.csv", [], ["metric_name", "value", "display_format", "notes"])

    manifest = {
        "experiment_id": "1b_non_linguistic_control",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(path.relative_to(cfg.ROOT)) for path in cfg.SOURCE_BUNDLES],
        "outputs": {
            "city_population_control.csv": {"rows": 1, "schema": fields},
            "aggregate_statistics.csv": {"rows": 0, "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "notes": ["No fresh city-data download or symbolic search was run; values are normalized from saved summaries."],
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

