from __future__ import annotations

import argparse
import csv
import json
import statistics
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


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def aggregate_rows(breakthrough: dict, simulation: dict) -> list[dict[str, str]]:
    residual = breakthrough["residual_summary"]
    head = breakthrough["head_summary"]
    phase = breakthrough["phase_summary"]
    basis_rows = breakthrough["english_basis_rows"]
    sim = simulation["summary"]
    smooth_hurt_count = sum(float(row["smooth_improvement"]) < 0.0 for row in residual["rows"])
    smooth_nonhelp_count = sum(float(row["smooth_improvement"]) <= cfg.HELPFUL_STEP2_GAIN_THRESHOLD for row in residual["rows"])
    rows: list[tuple[str, object, str, str]] = [
        ("moe_residual_step2_help_count", residual["moe_count_gt_1e3"], "integer", "MOE residual improvement > 0.001."),
        ("moe_residual_step2_gain_median", residual["moe_median_improvement"], "signed_decimal_3", "Median MOE residual step-2 improvement."),
        ("smooth_residual_step2_gain_median", residual["smooth_median_improvement"], "signed_decimal_3", "Median smooth residual step-2 improvement."),
        ("smooth_residual_step2_hurt_count", smooth_hurt_count, "integer", "Count where smooth residual step-2 improvement is negative."),
        ("smooth_residual_step2_help_count", residual["smooth_count_gt_1e3"], "integer", "Smooth residual improvement > 0.001."),
        ("smooth_residual_step2_nonhelp_count", smooth_nonhelp_count, "integer", "Count where smooth residual step-2 does not exceed the helpful threshold."),
        ("head_basis_r2_median_top200", statistics.median(float(row["head_r2"]) for row in basis_rows), "decimal_6", "Median top-200 cubic head-basis R^2."),
        ("slant_vs_c_correlation", phase["corr_slant_c"], "decimal_6", "Correlation between head-basis slant and c."),
        ("median_winner_minus_euclidean_gap_full_rmse", head["median_full_euclidean_gap"], "decimal_6", "Median formula winner-vs-Euclidean full-RMSE gap from breakthrough probe."),
        ("median_winner_minus_euclidean_gap_top100_rmse", head["median_top100_euclidean_gap"], "decimal_6", "Median formula winner-vs-Euclidean top-100 RMSE gap from breakthrough probe."),
        ("empirical_euclidean_gap_full_rmse_median", sim["empirical_median_gap_full"], "decimal_6", "Auxiliary value imported from zipf_simulation_recovery to satisfy v4 no-builder-join mapping."),
        ("empirical_euclidean_gap_top100_rmse_median", sim["empirical_median_gap_top100"], "decimal_6", "Auxiliary value imported from zipf_simulation_recovery to satisfy v4 no-builder-join mapping."),
        ("smooth_synthetic_euclidean_gap_full_rmse_median", sim["smooth"]["replicate_median_winner_vs_euclidean_gap_full"], "decimal_6", "Auxiliary value imported from zipf_simulation_recovery to satisfy v4 no-builder-join mapping."),
        ("smooth_synthetic_euclidean_gap_top100_rmse_median", sim["smooth"]["replicate_median_winner_vs_euclidean_gap_top100"], "decimal_6", "Auxiliary value imported from zipf_simulation_recovery to satisfy v4 no-builder-join mapping."),
        ("zm_control_euclidean_gap_top100_rmse_median", sim["single_zm"]["replicate_median_winner_vs_euclidean_gap_top100"], "decimal_6", "Auxiliary value imported from zipf_simulation_recovery to satisfy v4 no-builder-join mapping."),
    ]
    return [
        {"metric_name": name, "value": repr(value), "display_format": fmt, "notes": notes}
        for name, value, fmt, notes in rows
    ]


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    breakthrough = load_json(cfg.BREAKTHROUGH_SUMMARY_JSON)
    simulation = load_json(cfg.SIMULATION_RECOVERY_SUMMARY_JSON)
    residual_rows = breakthrough["residual_summary"]["rows"]
    formula_rows = breakthrough["english_formula_rows"]
    basis_rows = breakthrough["english_basis_rows"]
    aggregate = aggregate_rows(breakthrough, simulation)
    residual_fields = [
        "slug",
        "corpus",
        "moe_base_rmse",
        "moe_winner_expr",
        "moe_step2_residual_rmse",
        "moe_step2_composite_rmse",
        "moe_improvement",
        "smooth_base_rmse",
        "smooth_winner_expr",
        "smooth_step2_residual_rmse",
        "smooth_step2_composite_rmse",
        "smooth_improvement",
    ]
    formula_fields = [
        "slug",
        "corpus",
        "winner_expr",
        "winner_family",
        "c",
        "full_best_formula",
        "top50_best_formula",
        "top100_best_formula",
        "top200_best_formula",
        "full_gap_euclidean",
        "top50_gap_euclidean",
        "top100_gap_euclidean",
        "top200_gap_euclidean",
    ]
    basis_fields = [
        "slug",
        "corpus",
        "family",
        "c",
        "linear_u",
        "quadratic_u2",
        "cubic_u3",
        "head_rmse",
        "head_r2",
        "head_k",
    ]
    write_csv(output_dir / "residual_search_per_corpus.csv", residual_rows, residual_fields)
    write_csv(output_dir / "head_formula_competition.csv", formula_rows, formula_fields)
    write_csv(output_dir / "head_basis_per_corpus.csv", basis_rows, basis_fields)
    write_csv(output_dir / "aggregate_statistics.csv", aggregate, ["metric_name", "value", "display_format", "notes"])
    manifest = {
        "experiment_id": "2c_mechanism_absorption_residual_space",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [
            str(cfg.BREAKTHROUGH_DIR.relative_to(cfg.ROOT)),
            str(cfg.SIMULATION_RECOVERY_DIR.relative_to(cfg.ROOT)),
        ],
        "outputs": {
            "residual_search_per_corpus.csv": {"rows": len(residual_rows), "schema": residual_fields},
            "head_formula_competition.csv": {"rows": len(formula_rows), "schema": formula_fields},
            "head_basis_per_corpus.csv": {"rows": len(basis_rows), "schema": basis_fields},
            "aggregate_statistics.csv": {"rows": len(aggregate), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "audit_flags": [
            "Five Euclidean-gap comparison aggregates are imported from zipf_simulation_recovery because v4 maps them to 2c to avoid LaTeX builder joins.",
        ],
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
