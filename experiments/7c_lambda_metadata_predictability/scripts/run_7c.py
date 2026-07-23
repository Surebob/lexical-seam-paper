from __future__ import annotations

import argparse
import csv
import json
import shutil
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


def metric(name: str, value: object, display_format: str, notes: str) -> dict[str, object]:
    return {
        "metric_name": name,
        "value": value,
        "display_format": display_format,
        "notes": notes,
    }


def build_per_corpus(summary: dict[str, object]) -> list[dict[str, object]]:
    fields = [
        "slug",
        "name",
        "lambda_k",
        "log10_lambda_k",
        "token_count",
        "vocab_size",
        "structure",
        "structure_code",
        "unit_count",
        "author_count",
        "era_span_years",
        "log_unit_count",
        "log_author_count",
        "log_era_span_years",
        "heterogeneity_score",
    ]
    return [{field: row.get(field, "") for field in fields} for row in summary["rows"]]


def build_correlation_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    correlations = summary["correlations"]
    for predictor in cfg.PREDICTOR_ORDER:
        item = correlations[predictor]
        rows.append(
            {
                "predictor": predictor,
                "pearson_r": item["pearson"],
                "spearman_r": item["spearman"],
                "abs_pearson_r": abs(item["pearson"]),
                "abs_spearman_r": abs(item["spearman"]),
                "n_corpora": len(summary["rows"]),
                "response": "log10_lambda_k",
                "source": "results/zipf_angle2_lambda_metadata/summary.json",
            }
        )
    return rows


def build_structure_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for structure, item in summary["by_structure"].items():
        rows.append(
            {
                "structure": structure,
                "count": item["count"],
                "median_lambda_k": item["median_lambda_k"],
                "median_log10_lambda_k": item["median_log10_lambda_k"],
            }
        )
    return rows


def build_aggregates(summary: dict[str, object], correlations: list[dict[str, object]]) -> list[dict[str, object]]:
    by_predictor = {row["predictor"]: row for row in correlations}
    manuscript_corrs = [by_predictor[name]["pearson_r"] for name in cfg.MANUSCRIPT_METADATA_PREDICTORS]
    all_corrs = [row["pearson_r"] for row in correlations]
    lambdas = [row["lambda_k"] for row in summary["rows"]]
    rows = [
        metric(
            "correlation_pearson_log_unit_count",
            by_predictor["log_unit_count"]["pearson_r"],
            "Pearson r",
            "Correlation between log10(lambda_k) and log unit count.",
        ),
        metric(
            "correlation_pearson_log_author_count",
            by_predictor["log_author_count"]["pearson_r"],
            "Pearson r",
            "Correlation between log10(lambda_k) and log author count.",
        ),
        metric(
            "correlation_pearson_log_era_span_years",
            by_predictor["log_era_span_years"]["pearson_r"],
            "Pearson r",
            "Correlation between log10(lambda_k) and log era span.",
        ),
        metric(
            "correlation_pearson_heterogeneity_score",
            by_predictor["heterogeneity_score"]["pearson_r"],
            "Pearson r",
            "Correlation between log10(lambda_k) and composite heterogeneity score.",
        ),
        metric(
            "max_abs_metadata_correlation",
            max(abs(value) for value in manuscript_corrs),
            "absolute Pearson r",
            "Maximum absolute Pearson correlation among manuscript-listed metadata predictors.",
        ),
        metric(
            "max_abs_all_available_correlation",
            max(abs(value) for value in all_corrs),
            "absolute Pearson r",
            "Maximum absolute Pearson correlation among all saved predictors, including structure code and corpus-size variables.",
        ),
        metric("lambda_k_min", min(lambdas), "scientific notation", "Minimum lambda_k in the metadata table."),
        metric("lambda_k_max", max(lambdas), "scientific notation", "Maximum lambda_k in the metadata table."),
        metric("corpus_count", len(summary["rows"]), "integer", "Number of English corpora in the lambda-metadata table."),
    ]
    return rows


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(cfg.SOURCE_SUMMARY.read_text(encoding="utf-8"))
    per_corpus = build_per_corpus(summary)
    correlations = build_correlation_rows(summary)
    structures = build_structure_rows(summary)
    aggregates = build_aggregates(summary, correlations)

    per_corpus_fields = [
        "slug",
        "name",
        "lambda_k",
        "log10_lambda_k",
        "token_count",
        "vocab_size",
        "structure",
        "structure_code",
        "unit_count",
        "author_count",
        "era_span_years",
        "log_unit_count",
        "log_author_count",
        "log_era_span_years",
        "heterogeneity_score",
    ]
    write_csv(output_dir / "lambda_metadata_per_corpus.csv", per_corpus, per_corpus_fields)
    write_csv(
        output_dir / "lambda_metadata_summary.csv",
        correlations,
        ["predictor", "pearson_r", "spearman_r", "abs_pearson_r", "abs_spearman_r", "n_corpora", "response", "source"],
    )
    write_csv(output_dir / "lambda_by_structure.csv", structures, ["structure", "count", "median_lambda_k", "median_log10_lambda_k"])
    write_csv(output_dir / "aggregate_statistics.csv", aggregates, ["metric_name", "value", "display_format", "notes"])

    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps({"source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))]}, indent=2),
        encoding="utf-8",
    )
    provenance = cfg.ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for src in [cfg.SOURCE_TABLE, cfg.SOURCE_SUMMARY, cfg.SOURCE_REPORT]:
        shutil.copy2(src, provenance / src.name)

    manifest = {
        "experiment_id": "7c_lambda_metadata_predictability",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
        "status": "complete",
        "outputs": {
            "lambda_metadata_per_corpus.csv": {"rows": len(per_corpus), "schema": per_corpus_fields},
            "lambda_metadata_summary.csv": {"rows": len(correlations)},
            "lambda_by_structure.csv": {"rows": len(structures)},
            "aggregate_statistics.csv": {"rows": len(aggregates), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "claim_map_rows_satisfied": [
            "line 504 max_abs_metadata_correlation",
            "line 536 lambda_k range restatement is also present here, though 7a remains the canonical Table 4 lambda source",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote 7c outputs: {len(per_corpus)} corpus rows, {len(correlations)} correlation rows, {len(aggregates)} aggregate rows")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()
