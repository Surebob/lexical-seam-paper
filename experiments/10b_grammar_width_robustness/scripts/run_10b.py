from __future__ import annotations

import csv
import json
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

import source_config as cfg


def load_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def metric(name: str, value: object, display_format: str, notes: str) -> dict[str, object]:
    return {
        "metric_name": name,
        "value": value,
        "display_format": display_format,
        "notes": notes,
    }


def build_diagnostic_rows(summary: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    per_mode: list[dict[str, object]] = []
    step_rows: list[dict[str, object]] = []
    top10_rows: list[dict[str, object]] = []
    for item in summary["results"]:
        per_mode.append(
            {
                "slug": item["slug"],
                "corpus": item["name"],
                "mode": item["mode"],
                "zm_c": item["zm_c"],
                "original_step2_expression": item["original_step2_expr"],
                "original_step2_math": item["original_step2_math"],
                "original_step2_rmse": item["original_step2_rmse"],
                "widened_step2_expression": item["widened_step2_expr"],
                "widened_step2_math": item["widened_step2_math"],
                "widened_step2_rmse": item["widened_step2_rmse"],
                "widened_step2_matches_original": item["widened_step2_matches_original"],
                "widened_step2_is_bregman": item["widened_step2_is_bregman"],
                "any_winner_uses_new_operator": item["any_winner_uses_new_operator"],
                "any_top10_uses_new_operator": item["any_top10_uses_new_operator"],
                "canonical_source": "results/zipf_widened_grammar_diagnostic/summary.json",
            }
        )
        for step in item["step_rows"]:
            step_rows.append(
                {
                    "slug": item["slug"],
                    "corpus": item["name"],
                    "mode": item["mode"],
                    "step": step["step"],
                    "winner_expression": step["winner_expr"],
                    "winner_math": step["winner_math"],
                    "winner_rmse": step["winner_rmse"],
                    "winner_uses_new_operator": step["winner_uses_new_operator"],
                    "top10_has_new_operator": step["top10_has_new_operator"],
                    "canonical_source": "results/zipf_widened_grammar_diagnostic/summary.json",
                }
            )
            for candidate in step["top10"]:
                top10_rows.append(
                    {
                        "slug": item["slug"],
                        "corpus": item["name"],
                        "mode": item["mode"],
                        "step": step["step"],
                        "rank": candidate["rank"],
                        "expression": candidate["expr"],
                        "math": candidate["math"],
                        "rmse": candidate["rmse"],
                        "uses_new_operator": candidate["uses_new_operator"],
                        "canonical_source": "results/zipf_widened_grammar_diagnostic/summary.json",
                    }
                )
    return per_mode, step_rows, top10_rows


def build_lowc_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in summary["rows"]:
        manifold = item["manifold"]
        rows.append(
            {
                "slug": item["slug"],
                "corpus": item["name"],
                "language": item["language"],
                "zm_c": item["zm_c"],
                "selection_note": item["selection_note"],
                "original_step2_expression": item["original_step2_expr"],
                "original_step2_math": item["original_step2_math"],
                "widened_step2_expression": item["widened_step2_expr"],
                "widened_step2_math": item["widened_step2_math"],
                "widened_step2_rmse": item["widened_step2_rmse"],
                "widened_matches_original": item["widened_matches_original"],
                "widened_is_bregman": item["widened_is_bregman"],
                "cos_vs_exp": manifold["cos_vs_exp"],
                "cos_vs_xx": manifold["cos_vs_xx"],
                "cos_vs_is": manifold["cos_vs_is"],
                "span_r2": manifold["span_r2"],
                "span_coeffs_json": json.dumps(manifold["span_coeffs"]),
                "manifold_verdict": manifold["verdict"],
                "top10_json": json.dumps(item["top10"], sort_keys=True),
                "canonical_source": "results/zipf_widened_grammar_extended/summary.json",
            }
        )
    return rows


def build_aggregates(
    diagnostic_summary: dict[str, object],
    per_mode: list[dict[str, object]],
    lowc_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    zm_rows = [row for row in per_mode if row["mode"] == "zm_residual"]
    post_rows = [row for row in per_mode if row["mode"] == "post_is_residual"]
    high_c_rows = [row for row in zm_rows if row["slug"] in {"shakespeare", "war_and_peace"}]
    english_lowc = [row for row in lowc_rows if row["language"] == "English"]
    all_lowc = lowc_rows
    winner_counts = Counter(str(row["widened_step2_expression"]) for row in all_lowc)
    english_winner_counts = Counter(str(row["widened_step2_expression"]) for row in english_lowc)
    return [
        metric("beam_width", diagnostic_summary["beam_width"], "integer protocol constant", "Widened grammar deterministic beam width."),
        metric("max_steps", diagnostic_summary["max_steps"], "integer protocol constant", "Maximum search depth in the original widened diagnostic."),
        metric("keep_all_until_step", diagnostic_summary["keep_all_until_step"], "integer protocol constant", "Search keeps all expressions until this step."),
        metric("new_unary_operator_count", len(diagnostic_summary["new_unary_ops"]), "integer count", "Number of added unary operators."),
        metric("new_unary_operators", json.dumps(diagnostic_summary["new_unary_ops"]), "JSON string list", "Added unary operators."),
        metric("diagnostic_mode_count", len(per_mode), "integer count", "Rows in the high-c/deeper-step diagnostic."),
        metric("diagnostic_zm_mode_count", len(zm_rows), "integer count", "Raw ZM residual rows in the diagnostic."),
        metric("diagnostic_post_is_mode_count", len(post_rows), "integer count", "Post-IS residual rows in the diagnostic."),
        metric("high_c_zm_step2_survives_count", sum(bool(row["widened_step2_matches_original"]) for row in high_c_rows), "X/2 count", "Shakespeare and War and Peace retain the original IS-Bregman step-2 winner under widened grammar."),
        metric("high_c_zm_step2_prior_bregman_count", sum(bool(row["widened_step2_is_bregman"]) for row in high_c_rows), "X/2 count", "High-c exemplar rows whose widened step-2 winner is the prior Bregman winner."),
        metric("lowc_total_count", len(all_lowc), "integer count", "Low-c and multilingual rows in the extended step-2 manifold check."),
        metric("lowc_english_count", len(english_lowc), "integer count", "English low-c rows in the extended step-2 manifold check."),
        metric("lowc_all_manifold_yes_count", sum(row["manifold_verdict"] == "yes" for row in all_lowc), "X/17 count", "Rows with manifold verdict yes."),
        metric("lowc_all_span_r2_gt_0_975_count", sum(float(row["span_r2"]) > 0.975 for row in all_lowc), "X/17 count", "Rows with weighted centered span R^2 above 0.975."),
        metric("lowc_all_min_span_r2", min(float(row["span_r2"]) for row in all_lowc), "R^2", "Minimum low-c span R^2."),
        metric("lowc_english_erf_sin_flip_count", sum(("erf" in str(row["widened_step2_expression"]) and "sin" in str(row["widened_step2_expression"])) for row in english_lowc), "X/14 count", "English low-c rows whose widened step-2 winner is erf-sin or sin-erf."),
        metric("lowc_all_erf_sin_winner_count", sum(("erf" in str(row["widened_step2_expression"]) and "sin" in str(row["widened_step2_expression"])) for row in all_lowc), "X/17 count", "All low-c/multilingual rows whose widened winner is erf-sin or sin-erf."),
        metric("lowc_all_winner_counts", json.dumps(dict(sorted(winner_counts.items()))), "JSON count map", "Widened step-2 winner-expression counts across all low-c/multilingual rows."),
        metric("lowc_english_winner_counts", json.dumps(dict(sorted(english_winner_counts.items()))), "JSON count map", "Widened step-2 winner-expression counts across English low-c rows."),
    ]


def archive_sources() -> None:
    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps({name: str(path.relative_to(cfg.ROOT)) for name, path in cfg.SOURCE_BUNDLES.items()}, indent=2) + "\n",
        encoding="utf-8",
    )
    provenance = cfg.ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for name, bundle in cfg.SOURCE_BUNDLES.items():
        dest = provenance / name
        dest.mkdir(parents=True, exist_ok=True)
        for src in sorted(bundle.iterdir()):
            if src.is_file():
                shutil.copy2(src, dest / src.name)


def write_manifest(outputs: dict[str, list[dict[str, object]]], schemas: dict[str, list[str]]) -> None:
    manifest = {
        "experiment_id": "10b_grammar_width_robustness",
        "status": "complete_from_historical_bundles",
        "migration_type": "historical_bundle_consolidation",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundles": {name: str(path.relative_to(cfg.ROOT)) for name, path in cfg.SOURCE_BUNDLES.items()},
        "outputs": {
            filename: {"rows": len(rows), "schema": schemas[filename]}
            for filename, rows in outputs.items()
        },
        "claim_map_rows_satisfied": [
            "high-c widened grammar step-2 robustness on Shakespeare and War and Peace",
            "low-c widened grammar erf/sin flip count",
            "low-c manifold membership under widened grammar",
        ],
    }
    cfg.OUTPUTS["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    diagnostic_summary = load_json(cfg.SOURCE_FILES["diagnostic_summary"])
    extended_summary = load_json(cfg.SOURCE_FILES["extended_summary"])
    per_mode, diagnostic_steps, diagnostic_top10 = build_diagnostic_rows(diagnostic_summary)
    lowc_rows = build_lowc_rows(extended_summary)
    aggregate_rows = build_aggregates(diagnostic_summary, per_mode, lowc_rows)

    schemas = {
        "widened_diagnostic_per_mode.csv": [
            "slug",
            "corpus",
            "mode",
            "zm_c",
            "original_step2_expression",
            "original_step2_math",
            "original_step2_rmse",
            "widened_step2_expression",
            "widened_step2_math",
            "widened_step2_rmse",
            "widened_step2_matches_original",
            "widened_step2_is_bregman",
            "any_winner_uses_new_operator",
            "any_top10_uses_new_operator",
            "canonical_source",
        ],
        "widened_diagnostic_step_winners.csv": [
            "slug",
            "corpus",
            "mode",
            "step",
            "winner_expression",
            "winner_math",
            "winner_rmse",
            "winner_uses_new_operator",
            "top10_has_new_operator",
            "canonical_source",
        ],
        "widened_diagnostic_top10.csv": [
            "slug",
            "corpus",
            "mode",
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "uses_new_operator",
            "canonical_source",
        ],
        "widened_lowc_manifold_summary.csv": [
            "slug",
            "corpus",
            "language",
            "zm_c",
            "selection_note",
            "original_step2_expression",
            "original_step2_math",
            "widened_step2_expression",
            "widened_step2_math",
            "widened_step2_rmse",
            "widened_matches_original",
            "widened_is_bregman",
            "cos_vs_exp",
            "cos_vs_xx",
            "cos_vs_is",
            "span_r2",
            "span_coeffs_json",
            "manifold_verdict",
            "top10_json",
            "canonical_source",
        ],
        "aggregate_statistics.csv": ["metric_name", "value", "display_format", "notes"],
    }
    outputs = {
        "widened_diagnostic_per_mode.csv": per_mode,
        "widened_diagnostic_step_winners.csv": diagnostic_steps,
        "widened_diagnostic_top10.csv": diagnostic_top10,
        "widened_lowc_manifold_summary.csv": lowc_rows,
        "aggregate_statistics.csv": aggregate_rows,
    }
    write_csv(cfg.OUTPUTS["diagnostic"], per_mode, schemas["widened_diagnostic_per_mode.csv"])
    write_csv(cfg.OUTPUTS["diagnostic_steps"], diagnostic_steps, schemas["widened_diagnostic_step_winners.csv"])
    write_csv(cfg.OUTPUTS["diagnostic_top10"], diagnostic_top10, schemas["widened_diagnostic_top10.csv"])
    write_csv(cfg.OUTPUTS["lowc_summary"], lowc_rows, schemas["widened_lowc_manifold_summary.csv"])
    write_csv(cfg.OUTPUTS["aggregate"], aggregate_rows, schemas["aggregate_statistics.csv"])
    archive_sources()
    write_manifest(outputs, schemas)
    print(
        "Wrote 10b outputs: "
        f"{len(per_mode)} diagnostic rows, "
        f"{len(diagnostic_steps)} step-winner rows, "
        f"{len(diagnostic_top10)} top10 rows, "
        f"{len(lowc_rows)} low-c rows, "
        f"{len(aggregate_rows)} aggregate rows"
    )


if __name__ == "__main__":
    main()
