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


def width_imbalance(width_ratio: float) -> float:
    if width_ratio <= 0:
        return float("inf")
    return max(width_ratio, 1.0 / width_ratio)


def build_per_corpus(table_rows: list[dict[str, str]], summary: dict[str, object]) -> list[dict[str, object]]:
    by_slug = {row["slug"]: row for row in summary["rows"]}
    out: list[dict[str, object]] = []
    for row in table_rows:
        item = by_slug[row["slug"]]
        params = item["best_lambda"]["selection"]["params"]
        width_ratio_value = f(row["width_ratio"])
        out.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "token_count": item["token_count"],
                "vocabulary_size": item["vocab_size"],
                "lambda_k": f(row["lambda"]),
                "w_left": f(row["w_left"]),
                "w_right": f(row["w_right"]),
                "width_ratio_w_right_over_w_left": width_ratio_value,
                "width_imbalance_ratio": width_imbalance(width_ratio_value),
                "material_asymmetry": width_imbalance(width_ratio_value) > 1.5,
                "k": params["k"],
                "transition_fraction": params["transition_fraction"],
                "asymmetric_full_avg_nll": f(row["full_asym"]),
                "softk_full_avg_nll": f(row["full_softk"]),
                "moe_full_avg_nll": f(row["full_moe"]),
                "zm_full_avg_nll": f(row["full_zm"]),
                "asymmetric_minus_softk": f(row["full_asym"]) - f(row["full_softk"]),
                "asymmetric_minus_moe": f(row["full_asym"]) - f(row["full_moe"]),
                "asymmetric_minus_zm": f(row["full_asym"]) - f(row["full_zm"]),
                "full_step2_help": bool(int(row["full_step2_help"])),
                "full_step2_gain": f(row["full_step2_gain"]),
                "canonical_source": "results/zipf_angle3_asymmetric_gate_splitfit/asymmetric_gate_table.csv",
            }
        )
    return out


def build_aggregates(summary: dict[str, object]) -> list[dict[str, object]]:
    full = summary["cutoffs"]["full"]
    rows = [
        metric("median_width_ratio", summary["width_ratio"]["median"], "ratio", "Median w_right / w_left across 25 corpus-level asymmetric fits."),
        metric("mean_width_ratio", summary["width_ratio"]["mean"], "ratio", "Mean w_right / w_left across 25 corpus-level asymmetric fits."),
        metric("material_asymmetry_count", summary["width_ratio"]["material_asymmetry_count"], "X/25 count", "Count with greater than 1.5x width imbalance."),
        metric("asymmetric_beats_softk_count", full["asym_beats_softk"], "X/25 count", "Full held-out NLL count where asymmetric gate beats symmetric soft-k."),
        metric("asymmetric_beats_moe_count", full["asym_beats_moe"], "X/25 count", "Full held-out NLL count where asymmetric gate beats MOEZipf."),
        metric("asymmetric_beats_zm_count", full["asym_beats_zm"], "X/25 count", "Full held-out NLL count where asymmetric gate beats ZM."),
        metric("median_asymmetric_minus_softk", full["median_asym_minus_softk"], "signed average NLL delta", "Median full held-out NLL delta, asymmetric minus symmetric soft-k."),
        metric("median_asymmetric_minus_moe", full["median_asym_minus_moe"], "signed average NLL delta", "Median full held-out NLL delta, asymmetric minus MOEZipf."),
        metric("median_asymmetric_minus_zm", full["median_asym_minus_zm"], "signed average NLL delta", "Median full held-out NLL delta, asymmetric minus ZM."),
        metric("asymmetric_step2_help_count", full["step2_help_count"], "X/25 count", "Count where step-2 search helps on full asymmetric-gate residuals."),
    ]
    for model in cfg.MODEL_ORDER:
        rows.append(
            metric(
                f"full_winner_count_{model}",
                full["winner_counts"].get(model, 0),
                "X/25 count",
                "Full held-out winner count with asymmetric gate included as a fifth candidate.",
            )
        )
    return rows


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    table_rows = read_csv(cfg.SOURCE_TABLE)
    summary = json.loads(cfg.SOURCE_SUMMARY.read_text(encoding="utf-8"))
    per_corpus = build_per_corpus(table_rows, summary)
    aggregates = build_aggregates(summary)

    per_corpus_fields = [
        "slug",
        "corpus",
        "token_count",
        "vocabulary_size",
        "lambda_k",
        "w_left",
        "w_right",
        "width_ratio_w_right_over_w_left",
        "width_imbalance_ratio",
        "material_asymmetry",
        "k",
        "transition_fraction",
        "asymmetric_full_avg_nll",
        "softk_full_avg_nll",
        "moe_full_avg_nll",
        "zm_full_avg_nll",
        "asymmetric_minus_softk",
        "asymmetric_minus_moe",
        "asymmetric_minus_zm",
        "full_step2_help",
        "full_step2_gain",
        "canonical_source",
    ]
    write_csv(output_dir / "asymmetric_gate_per_corpus.csv", per_corpus, per_corpus_fields)
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
        "experiment_id": "7d_asymmetric_gate_diagnostic",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
        "status": "complete",
        "outputs": {
            "asymmetric_gate_per_corpus.csv": {"rows": len(per_corpus), "schema": per_corpus_fields},
            "aggregate_statistics.csv": {"rows": len(aggregates), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "claim_map_rows_satisfied": [
            "line 413 / 504 asymmetric gate width ratio, material asymmetry count, held-out delta, and step-2 help count",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote 7d outputs: {len(per_corpus)} corpus rows, {len(aggregates)} aggregate rows")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()
