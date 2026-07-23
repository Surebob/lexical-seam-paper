from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

import source_config as cfg


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def f(value: object) -> float:
    return float(value)


def metric(name: str, value: object, display_format: str, notes: str) -> dict[str, object]:
    return {"metric_name": name, "value": value, "display_format": display_format, "notes": notes}


def build_per_book(summary: dict[str, object]) -> list[dict[str, object]]:
    fields = [
        "book",
        "heading",
        "token_count",
        "vocabulary_size",
        "test_token_count",
        "softk_lambda",
        "zipf_test_avg_nll",
        "zm_test_avg_nll",
        "moe_test_avg_nll",
        "softk_test_avg_nll",
        "softk_minus_zm",
        "softk_minus_moe",
        "step2_gain",
        "step2_helps",
        "step2_expression",
        "softk_rmse",
    ]
    rows: list[dict[str, object]] = []
    for item in summary["rows"]:
        rows.append(
            {
                "book": item["name"],
                "heading": item["heading"],
                "token_count": item["token_count"],
                "vocabulary_size": item["vocab_size"],
                "test_token_count": item["test_token_count"],
                "softk_lambda": item["softk_lambda"],
                "zipf_test_avg_nll": item["zipf_test_avg_nll"],
                "zm_test_avg_nll": item["zm_test_avg_nll"],
                "moe_test_avg_nll": item["moe_test_avg_nll"],
                "softk_test_avg_nll": item["softk_test_avg_nll"],
                "softk_minus_zm": item["softk_minus_zm"],
                "softk_minus_moe": item["softk_minus_moe"],
                "step2_gain": item["step2_gain"],
                "step2_helps": item["step2_helps"],
                "step2_expression": item["step2_expr"],
                "softk_rmse": item["softk_rmse"],
            }
        )
    return rows


def load_whole_bible_context() -> dict[str, float]:
    rows = read_csv(cfg.PMF_VARIANT_TABLE)
    by_slug = {row["slug"]: row for row in rows}
    bible = by_slug["king_james_bible"]
    other24 = [f(row["softk_splitfit_step2_gain"]) for row in rows if row["slug"] != "king_james_bible"]
    return {
        "whole_bible_step2_gain": f(bible["softk_splitfit_step2_gain"]),
        "median_other24_softk_step2_gain": statistics.median(other24),
    }


def build_table5(summary: dict[str, object], context: dict[str, float]) -> list[dict[str, object]]:
    counts = summary["counts"]
    medians = summary["medians"]
    aggregate = summary["aggregate"]
    whole = aggregate["whole_bible_softk_test_avg_nll"]
    per_book = aggregate["bookwise_softk_test_avg_nll"]
    rows = [
        ("books_analyzed", "Books analyzed", f"{counts['n_books']}", counts["n_books"], "integer"),
        ("per_book_step2_help_count", "Per-book step-2 helps", f"{counts['step2_help_count']}/{counts['n_books']}", counts["step2_help_count"], "X/66 ratio"),
        ("per_book_softk_beats_zm_count", "Per-book soft-k beats ZM", f"{counts['softk_beats_zm']}/{counts['n_books']}", counts["softk_beats_zm"], "X/66 ratio"),
        ("per_book_softk_beats_moe_count", "Per-book soft-k beats MOE", f"{counts['softk_beats_moe']}/{counts['n_books']}", counts["softk_beats_moe"], "X/66 ratio"),
        ("median_per_book_softk_minus_zm", "Median per-book soft-k minus ZM", f"{medians['softk_minus_zm']:.6g}", medians["softk_minus_zm"], "signed decimal"),
        ("median_per_book_softk_minus_moe", "Median per-book soft-k minus MOE", f"{medians['softk_minus_moe']:.6g}", medians["softk_minus_moe"], "signed decimal"),
        ("median_per_book_step2_gain", "Median per-book step-2 gain", f"{medians['step2_gain']:.6g}", medians["step2_gain"], "signed decimal"),
        ("whole_bible_singlefit_softk_nll", "Whole-Bible single-fit soft-k held-out NLL", f"{whole:.3f}", whole, "average NLL"),
        ("aggregate_per_book_softk_nll", "Aggregate per-book soft-k held-out NLL", f"{per_book:.3f}", per_book, "average NLL"),
        ("improvement_from_decomposition", "Whole-fit minus aggregate per-book NLL", f"{whole - per_book:.3f}", whole - per_book, "average NLL delta"),
        ("per_book_step2_nonhelp_count", "Per-book step-2 non-help count", f"{counts['n_books'] - counts['step2_help_count']}/{counts['n_books']}", counts["n_books"] - counts["step2_help_count"], "X/66 ratio"),
        ("whole_bible_step2_gain", "Whole-Bible soft-k step-2 gain", f"{context['whole_bible_step2_gain']:.6g}", context["whole_bible_step2_gain"], "signed decimal"),
        ("median_other24_softk_step2_gain", "Median soft-k step-2 gain across other 24 English corpora", f"{context['median_other24_softk_step2_gain']:.6g}", context["median_other24_softk_step2_gain"], "signed decimal"),
    ]
    return [
        {
            "metric_name": metric_name,
            "metric_label": label,
            "rendered_value": rendered,
            "value": value,
            "display_format": display_format,
        }
        for metric_name, label, rendered, value, display_format in rows
    ]


def build_aggregates(table5: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        metric(row["metric_name"], row["value"], row["display_format"], f"Table 5 row: {row['metric_label']}.")
        for row in table5
    ]


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(cfg.SOURCE_SUMMARY.read_text(encoding="utf-8"))
    per_book = build_per_book(summary)
    context = load_whole_bible_context()
    table5 = build_table5(summary, context)
    aggregates = build_aggregates(table5)

    per_book_fields = [
        "book",
        "heading",
        "token_count",
        "vocabulary_size",
        "test_token_count",
        "softk_lambda",
        "zipf_test_avg_nll",
        "zm_test_avg_nll",
        "moe_test_avg_nll",
        "softk_test_avg_nll",
        "softk_minus_zm",
        "softk_minus_moe",
        "step2_gain",
        "step2_helps",
        "step2_expression",
        "softk_rmse",
    ]
    write_csv(output_dir / "bible_per_book.csv", per_book, per_book_fields)
    write_csv(output_dir / "table5_bible_summary.csv", table5, ["metric_name", "metric_label", "rendered_value", "value", "display_format"])
    write_csv(output_dir / "aggregate_statistics.csv", aggregates, ["metric_name", "value", "display_format", "notes"])

    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps(
            {
                "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
                "cross_experiment_dependency": str(cfg.PMF_VARIANT_TABLE.relative_to(cfg.ROOT)),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    provenance = cfg.ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for src in [cfg.SOURCE_TABLE, cfg.SOURCE_SUMMARY, cfg.SOURCE_REPORT, cfg.SOURCE_FIGURE]:
        if src.exists():
            shutil.copy2(src, provenance / src.name)

    manifest = {
        "experiment_id": "8_bible_per_book_decomposition",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
        "cross_experiment_dependencies": [str(cfg.PMF_VARIANT_TABLE.relative_to(cfg.ROOT))],
        "status": "complete",
        "outputs": {
            "bible_per_book.csv": {"rows": len(per_book), "schema": per_book_fields},
            "table5_bible_summary.csv": {"rows": len(table5)},
            "aggregate_statistics.csv": {"rows": len(aggregates), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "claim_map_rows_satisfied": [
            "line 417-432 / Table 5 Bible decomposition summary",
            "line 544-546 scope-restatement aggregate NLL and non-help count",
            "appendix checklist Bible rows",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote 8 outputs: {len(per_book)} per-book rows, {len(table5)} Table 5 rows")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()
