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


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_synthetic_rows(summary: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in summary["results"]:
        meta = cfg.MIXTURE_METADATA.get(item["name"], {})
        zm = item["zm_baseline"]
        step2 = item["step2"]
        rows.append(
            {
                "mixture_id": meta.get("mixture_id", item["name"]),
                "historical_name": item["name"],
                "alpha_1": meta.get("alpha_1", ""),
                "alpha_2": meta.get("alpha_2", ""),
                "weight_1": meta.get("weight_1", ""),
                "claimed_by_v4": str(meta.get("claimed_by_v4", False)),
                "token_equivalent": item["token_equivalent"],
                "vocab_equivalent": item["vocab_equivalent"],
                "zm_a": zm["a"],
                "zm_b": zm["b"],
                "zm_c": zm["c"],
                "zm_rmse": zm["rmse"],
                "step2_winner": step2["winner"],
                "step2_winner_display": step2["math"],
                "step2_rmse": step2["rmse"],
                "step2_helps": str(step2["helps"]),
                "top5_json": json.dumps(step2["top5"], ensure_ascii=False),
            }
        )
    order = {
        "small_gap": 0,
        "medium_gap": 1,
        "large_gap": 2,
        "historical_same_exponent_control_alpha0p8": 3,
    }
    rows.sort(key=lambda row: order.get(str(row["mixture_id"]), 99))
    return rows


def build_planted_rows(breakthrough: dict) -> list[dict[str, object]]:
    s = breakthrough["synthetic_summary"]
    return [
        {
            "mixture_id": "planted_small_gap",
            "mixture_name": "planted_small_gap",
            "source_bundle": str(cfg.BREAKTHROUGH_DIR.relative_to(cfg.ROOT)),
            "original_step2_winner": s["original_step2_winner"],
            "original_step2_rmse": s["original_step2_rmse"],
            "zm_baseline_rmse": s["zm_rmse"],
            "zm_step2_rmse": s["original_step2_rmse"],
            "zm_rmse": s["zm_rmse"],
            "zm_plus_step2_rmse": s["original_step2_rmse"],
            "moezipf_alpha": s["moe_alpha"],
            "moezipf_beta": s["moe_beta"],
            "moezipf_baseline_rmse": s["moe_base_rmse"],
            "moezipf_step2_rmse": s["moe_step2_composite_rmse"],
            "moe_rmse": s["moe_base_rmse"],
            "moe_plus_step2_rmse": s["moe_step2_composite_rmse"],
            "moezipf_step2_winner": s["moe_step2_winner"],
            "moezipf_step2_winner_math": "log(x)^2" if s["moe_step2_winner"] == "mul[log[x],log[x]]" else s["moe_step2_winner"],
            "step2_winner_display": "log(x)^2" if s["moe_step2_winner"] == "mul[log[x],log[x]]" else s["moe_step2_winner"],
            "moezipf_step2_improvement": s["moe_step2_improvement"],
            "integer_scale_factor": s["integer_scale_factor"],
        }
    ]


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    synthetic = load_json(cfg.SYNTHETIC_MIXTURE_SUMMARY_JSON)
    breakthrough = load_json(cfg.BREAKTHROUGH_SUMMARY_JSON)
    synthetic_rows = build_synthetic_rows(synthetic)
    planted_rows = build_planted_rows(breakthrough)
    aggregate_rows = [
        {
            "metric_name": "synthetic_mixture_configuration_count",
            "value": "4",
            "display_format": "integer",
            "notes": "Four historical non-component rows are present, but the same-exponent control is alpha=0.8/0.8, not the manuscript/map's claimed alpha=1.5/1.5; see BLOCKED.md.",
        }
    ]
    synthetic_fields = [
        "mixture_id",
        "historical_name",
        "alpha_1",
        "alpha_2",
        "weight_1",
        "claimed_by_v4",
        "token_equivalent",
        "vocab_equivalent",
        "zm_a",
        "zm_b",
        "zm_c",
        "zm_rmse",
        "step2_winner",
        "step2_winner_display",
        "step2_rmse",
        "step2_helps",
        "top5_json",
    ]
    planted_fields = [
        "mixture_id",
        "mixture_name",
        "source_bundle",
        "original_step2_winner",
        "original_step2_rmse",
        "zm_baseline_rmse",
        "zm_step2_rmse",
        "zm_rmse",
        "zm_plus_step2_rmse",
        "moezipf_alpha",
        "moezipf_beta",
        "moezipf_baseline_rmse",
        "moezipf_step2_rmse",
        "moe_rmse",
        "moe_plus_step2_rmse",
        "moezipf_step2_winner",
        "moezipf_step2_winner_math",
        "step2_winner_display",
        "moezipf_step2_improvement",
        "integer_scale_factor",
    ]
    write_csv(output_dir / "synthetic_mixture_runs.csv", synthetic_rows, synthetic_fields)
    write_csv(output_dir / "planted_mixture_runs.csv", planted_rows, planted_fields)
    write_csv(output_dir / "aggregate_statistics.csv", aggregate_rows, ["metric_name", "value", "display_format", "notes"])
    blocked = (
        "# Experiment 2b BLOCKED Items\n\n"
        "- The v4 claim map and manuscript refer to `mixture_id=same_exponent_control` with `alpha_1 = alpha_2 = 1.5`.\n"
        "- The saved historical bundle `results/zipf_synthetic_mixture` instead contains `control_two_copies_alpha0.8`, represented here as `historical_same_exponent_control_alpha0p8` with `alpha_1 = alpha_2 = 0.8`.\n"
        "- No saved output in the inspected source bundles reconstructs a `1.5/1.5` same-exponent control without rerunning synthetic generation/search.\n"
    )
    (EXPERIMENT_DIR / "BLOCKED.md").write_text(blocked, encoding="utf-8")
    manifest = {
        "experiment_id": "2b_synthetic_two_regime_mixtures",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [
            str(cfg.SYNTHETIC_MIXTURE_DIR.relative_to(cfg.ROOT)),
            str(cfg.BREAKTHROUGH_DIR.relative_to(cfg.ROOT)),
        ],
        "outputs": {
            "synthetic_mixture_runs.csv": {"rows": len(synthetic_rows), "schema": synthetic_fields},
            "planted_mixture_runs.csv": {"rows": len(planted_rows), "schema": planted_fields},
            "aggregate_statistics.csv": {"rows": len(aggregate_rows), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "blocked": ["claimed alpha=1.5/1.5 same-exponent control is not present in source bundles"],
        "audit_flags": ["line 277 planted-mixture values are migrated from zipf_breakthrough_probe into 2b to repair the narrative provenance mismatch"],
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
