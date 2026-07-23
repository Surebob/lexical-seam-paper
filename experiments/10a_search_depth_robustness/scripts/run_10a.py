from __future__ import annotations

import csv
import json
import shutil
import sys
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


def polyfit_map(polyfits: list[dict[str, object]]) -> dict[int, dict[str, object]]:
    return {int(item["degree"]): item for item in polyfits}


def build_step10_rows(source_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in source_rows:
        fits = polyfit_map(item["polyfits"])
        rows.append(
            {
                "slug": item["slug"],
                "corpus": item["name"],
                "step2_expression": item["step2_expr"],
                "step10_expression": item["step10_expr"],
                "grid_min": item["grid_min"],
                "grid_max": item["grid_max"],
                "difference_min": item["difference_min"],
                "difference_max": item["difference_max"],
                "difference_mean_abs": item["difference_mean_abs"],
                "difference_std": item["difference_std"],
                "best_poly_degree": item["best_poly_degree"],
                "best_poly_r2": item["best_poly_r2"],
                "polyfit_degree3_r2": fits[3]["r2"],
                "polyfit_degree4_r2": fits[4]["r2"],
                "polyfit_degree5_r2": fits[5]["r2"],
                "polyfit_degree3_coefficients_json": json.dumps(fits[3]["coefficients"]),
                "polyfit_degree4_coefficients_json": json.dumps(fits[4]["coefficients"]),
                "polyfit_degree5_coefficients_json": json.dumps(fits[5]["coefficients"]),
                "plot_path": item["plot_path"],
                "canonical_source": "results/zipf_step10_ablation/summary.json",
            }
        )
    return rows


def build_decomposition_rows(source_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in source_rows:
        models = polyfit_map(item["poly_models"])
        rows.append(
            {
                "slug": item["slug"],
                "corpus": item["name"],
                "zm_rmse": item["zm_rmse"],
                "step2_rmse": item["step2_rmse"],
                "step10_rmse": item["step10_rmse"],
                "poly3_rmse": models[3]["rmse"],
                "poly4_rmse": models[4]["rmse"],
                "poly5_rmse": models[5]["rmse"],
                "poly3_coefficients_json": json.dumps(models[3]["coefficients"]),
                "poly4_coefficients_json": json.dumps(models[4]["coefficients"]),
                "poly5_coefficients_json": json.dumps(models[5]["coefficients"]),
                "poly5_minus_step10_rmse": float(models[5]["rmse"]) - float(item["step10_rmse"]),
                "canonical_source": "results/zipf_head_poly_decomposition/summary.json",
            }
        )
    return rows


def build_transfer_rows(source: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in source["transfers"]:
        rows.append(
            {
                "source_corpus": item["source"],
                "target_corpus": item["target"],
                "source_poly5_coefficients_json": json.dumps(item["source_poly5_coefficients"]),
                "target_zm_rmse": item["target_zm_rmse"],
                "target_step2_rmse": item["target_step2_rmse"],
                "target_in_domain_poly5_rmse": item["target_in_domain_poly5_rmse"],
                "target_step10_rmse": item["target_step10_rmse"],
                "zero_shot_rmse": item["zero_shot_rmse"],
                "zero_shot_minus_in_domain_poly5": float(item["zero_shot_rmse"]) - float(item["target_in_domain_poly5_rmse"]),
                "zero_shot_minus_step10": float(item["zero_shot_rmse"]) - float(item["target_step10_rmse"]),
                "canonical_source": "results/zipf_head_poly_transfer/summary.json",
            }
        )
    return rows


def validate_overlap(decomp_rows: list[dict[str, object]], transfer_source: dict[str, object]) -> None:
    by_slug = {row["slug"]: row for row in decomp_rows}
    for slug, record in transfer_source["records"].items():
        row = by_slug[slug]
        checks = [
            ("zm_rmse", record["zm_rmse"], row["zm_rmse"]),
            ("step2_rmse", record["step2_rmse"], row["step2_rmse"]),
            ("step10_rmse", record["step10_rmse"], row["step10_rmse"]),
            ("poly5_rmse", record["poly5_rmse"], row["poly5_rmse"]),
        ]
        for name, left, right in checks:
            if abs(float(left) - float(right)) > 1e-12:
                raise RuntimeError(f"Overlap drift for {slug} {name}: {left} vs {right}")
        if json.dumps(record["poly5_coefficients"]) != row["poly5_coefficients_json"]:
            raise RuntimeError(f"Overlap drift for {slug} poly5 coefficients")


def build_aggregate_rows(
    step10_rows: list[dict[str, object]],
    decomp_rows: list[dict[str, object]],
    transfer_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    step10 = {row["slug"]: row for row in step10_rows}
    decomp = {row["slug"]: row for row in decomp_rows}
    transfer = {(row["source_corpus"], row["target_corpus"]): row for row in transfer_rows}
    return [
        metric("shakespeare_step10_minus_step2_poly5_r2", step10["shakespeare"]["polyfit_degree5_r2"], "R^2", "Degree-5 polynomial fit to step10 minus step2 on the 1000-point grid."),
        metric("war_and_peace_step10_minus_step2_poly5_r2", step10["war_and_peace"]["polyfit_degree5_r2"], "R^2", "Degree-5 polynomial fit to step10 minus step2 on the 1000-point grid."),
        metric("shakespeare_poly5_direct_rmse", decomp["shakespeare"]["poly5_rmse"], "RMSE", "Direct ZM + step-2 + poly5 residual-decomposition RMSE."),
        metric("shakespeare_step10_rmse", decomp["shakespeare"]["step10_rmse"], "RMSE", "Step-10 monster RMSE."),
        metric("war_and_peace_poly5_direct_rmse", decomp["war_and_peace"]["poly5_rmse"], "RMSE", "Direct ZM + step-2 + poly5 residual-decomposition RMSE."),
        metric("war_and_peace_step10_rmse", decomp["war_and_peace"]["step10_rmse"], "RMSE", "Step-10 monster RMSE."),
        metric("shakespeare_to_war_and_peace_zero_shot_rmse", transfer[("Shakespeare", "War and Peace")]["zero_shot_rmse"], "RMSE", "Zero-shot degree-5 polynomial transfer from Shakespeare to War and Peace."),
        metric("war_and_peace_to_shakespeare_zero_shot_rmse", transfer[("War and Peace", "Shakespeare")]["zero_shot_rmse"], "RMSE", "Zero-shot degree-5 polynomial transfer from War and Peace to Shakespeare."),
        metric("saved_transfer_count", len(transfer_rows), "integer count", "Number of transfer rows saved in the historical bundle."),
        metric("missing_transfer_claims_count", len(cfg.MISSING_TRANSFER_CLAIMS), "integer count", "Transfer claims present in manuscript prose but absent from saved outputs."),
    ]


def write_blocked() -> None:
    path = cfg.EXPERIMENT_DIR / "BLOCKED.md"
    missing = "\n".join(f"- {item}" for item in cfg.MISSING_TRANSFER_CLAIMS)
    path.write_text(
        "# Experiment 10a Partial Block\n\n"
        "Canonical migration completed for the saved source bundles, but the current manuscript "
        "contains transfer claims that are not reconstructible from the inspected historical outputs.\n\n"
        "Missing saved transfer rows:\n\n"
        f"{missing}\n\n"
        "The saved `zipf_head_poly_transfer` bundle contains only Shakespeare -> War and Peace and "
        "War and Peace -> Shakespeare. Do not fill these missing rows from manuscript prose. They require "
        "either manuscript revision or an explicitly authorized fresh run.\n",
        encoding="utf-8",
    )


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
    shutil.copy2(cfg.EXPERIMENT_DIR / "CONSOLIDATION_PLAN.md", provenance / "CONSOLIDATION_PLAN.md")


def write_manifest(outputs: dict[str, list[dict[str, object]]], schemas: dict[str, list[str]]) -> None:
    manifest = {
        "experiment_id": "10a_search_depth_robustness",
        "status": "partial_blocked_missing_transfer_rows",
        "migration_type": "multi_bundle_consolidation",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundles": {name: str(path.relative_to(cfg.ROOT)) for name, path in cfg.SOURCE_BUNDLES.items()},
        "outputs": {
            filename: {"rows": len(rows), "schema": schemas[filename]}
            for filename, rows in outputs.items()
        },
        "blocked_items": cfg.MISSING_TRANSFER_CLAIMS,
        "claim_map_rows_satisfied": [
            "step-10 minus step-2 degree-5 polynomial R^2 for Shakespeare and War and Peace",
            "direct polynomial decomposition RMSE comparisons for Shakespeare and War and Peace",
            "two saved zero-shot transfer rows",
        ],
        "claim_map_rows_unsatisfied": [
            "War and Peace -> Bible transfer",
            "Moby Dick -> Bible transfer",
        ],
    }
    cfg.OUTPUTS["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    step10_source = load_json(cfg.SOURCE_FILES["step10_summary"])
    decomp_source = load_json(cfg.SOURCE_FILES["decomposition_summary"])
    transfer_source = load_json(cfg.SOURCE_FILES["transfer_summary"])

    step10_rows = build_step10_rows(step10_source)
    decomp_rows = build_decomposition_rows(decomp_source)
    transfer_rows = build_transfer_rows(transfer_source)
    validate_overlap(decomp_rows, transfer_source)
    aggregate_rows = build_aggregate_rows(step10_rows, decomp_rows, transfer_rows)

    schemas = {
        "step10_ablation_per_corpus.csv": [
            "slug",
            "corpus",
            "step2_expression",
            "step10_expression",
            "grid_min",
            "grid_max",
            "difference_min",
            "difference_max",
            "difference_mean_abs",
            "difference_std",
            "best_poly_degree",
            "best_poly_r2",
            "polyfit_degree3_r2",
            "polyfit_degree4_r2",
            "polyfit_degree5_r2",
            "polyfit_degree3_coefficients_json",
            "polyfit_degree4_coefficients_json",
            "polyfit_degree5_coefficients_json",
            "plot_path",
            "canonical_source",
        ],
        "polynomial_decomposition_per_corpus.csv": [
            "slug",
            "corpus",
            "zm_rmse",
            "step2_rmse",
            "step10_rmse",
            "poly3_rmse",
            "poly4_rmse",
            "poly5_rmse",
            "poly3_coefficients_json",
            "poly4_coefficients_json",
            "poly5_coefficients_json",
            "poly5_minus_step10_rmse",
            "canonical_source",
        ],
        "poly_transfer.csv": [
            "source_corpus",
            "target_corpus",
            "source_poly5_coefficients_json",
            "target_zm_rmse",
            "target_step2_rmse",
            "target_in_domain_poly5_rmse",
            "target_step10_rmse",
            "zero_shot_rmse",
            "zero_shot_minus_in_domain_poly5",
            "zero_shot_minus_step10",
            "canonical_source",
        ],
        "aggregate_statistics.csv": ["metric_name", "value", "display_format", "notes"],
    }
    outputs = {
        "step10_ablation_per_corpus.csv": step10_rows,
        "polynomial_decomposition_per_corpus.csv": decomp_rows,
        "poly_transfer.csv": transfer_rows,
        "aggregate_statistics.csv": aggregate_rows,
    }

    write_csv(cfg.OUTPUTS["step10"], step10_rows, schemas["step10_ablation_per_corpus.csv"])
    write_csv(cfg.OUTPUTS["decomposition"], decomp_rows, schemas["polynomial_decomposition_per_corpus.csv"])
    write_csv(cfg.OUTPUTS["transfer"], transfer_rows, schemas["poly_transfer.csv"])
    write_csv(cfg.OUTPUTS["aggregate"], aggregate_rows, schemas["aggregate_statistics.csv"])
    write_blocked()
    archive_sources()
    write_manifest(outputs, schemas)
    print(
        "Wrote 10a outputs: "
        f"{len(step10_rows)} step10 rows, "
        f"{len(decomp_rows)} decomposition rows, "
        f"{len(transfer_rows)} transfer rows, "
        f"{len(aggregate_rows)} aggregate rows; partial BLOCKED for missing transfer rows"
    )


if __name__ == "__main__":
    main()
