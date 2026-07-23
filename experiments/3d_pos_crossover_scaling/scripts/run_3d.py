from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
ROOT = EXPERIMENT_DIR.parents[1]
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


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def bool_string(value: object) -> str:
    if isinstance(value, str):
        return value
    return "True" if bool(value) else "False"


def build_pos_points(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        out.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "vocabulary_size": row["V"],
                "sqrt_v": row["sqrt_V"],
                "k_crossover": row["k_crossover"],
                "alpha_per_corpus": row["alpha_per_corpus"],
                "log_vocabulary_size": row["log_V"],
                "log_k_crossover": row["log_k"],
                "crossover_deviation_from_sqrt_v": row["crossover_deviation_from_sqrt_V"],
                "censored_at_top_n": bool_string(row["censored_at_top_n"]),
            }
        )
    return out


def build_manual_validation(
    manual_rows: list[dict[str, str]],
    pos_points_by_slug: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in manual_rows:
        v = float(row["V"])
        k = float(row["k_crossover"])
        alpha = math.log(k) / math.log(v)
        spacy = pos_points_by_slug.get(row["slug"], {})
        spacy_k = float(spacy["k_crossover"]) if spacy else math.nan
        out.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "source": row["source"],
                "vocabulary_size": repr(v),
                "manual_k_crossover": repr(k),
                "manual_alpha_per_corpus": repr(alpha),
                "spacy_k_crossover": "" if math.isnan(spacy_k) else repr(spacy_k),
                "abs_manual_minus_spacy_rank_gap": "" if math.isnan(spacy_k) else repr(abs(k - spacy_k)),
                "has_top300_manual_csv": str((cfg.POS_MANUAL_DIR / f"{row['slug']}_top300_manual_pos.csv").exists()),
            }
        )
    return out


def build_aggregate_rows(pos_summary: dict, manual_rows: list[dict[str, object]]) -> list[dict[str, str]]:
    forced = pos_summary["forced_fit"]
    alpha_stats = pos_summary["alpha_stats"]
    max_gap = max(float(row["abs_manual_minus_spacy_rank_gap"]) for row in manual_rows if row["abs_manual_minus_spacy_rank_gap"] != "")
    rows = [
        ("top_pos_window", cfg.TOP_POS_WINDOW, "integer", "Protocol constant: automatic POS tagging window."),
        ("crossover_fraction", cfg.CROSSOVER_FRACTION, "decimal_2", "Protocol constant: closed-class fraction threshold."),
        ("manual_validation_window", cfg.MANUAL_VALIDATION_WINDOW, "integer", "Protocol constant: manual validation window."),
        ("manual_validation_corpus_count", len(manual_rows), "integer", "Four alpha points exist, but only three top300 manual CSV files are present; see README AUDIT."),
        ("max_manual_vs_spacy_rank_gap", max_gap, "integer", "Computed from available manual alpha points vs automatic POS crossover ranks; manuscript/map mention 15, data gives 10."),
        ("forced_alpha", forced["alpha"], "decimal_4", "Forced fit k_POS = V^alpha across 25 automatic POS points."),
        ("forced_alpha_ci_low", forced["alpha_ci_95"][0], "decimal_4", cfg.CI_METHOD),
        ("forced_alpha_ci_high", forced["alpha_ci_95"][1], "decimal_4", cfg.CI_METHOD),
        ("forced_alpha_vs_half_pvalue", alpha_stats["p_value"], "scientific", cfg.PVALUE_METHOD),
    ]
    return [
        {"metric_name": str(name), "value": repr(value), "display_format": fmt, "notes": notes}
        for name, value, fmt, notes in rows
    ]


def write_manifest(output_dir: Path, pos_rows: list[dict[str, object]], manual_rows: list[dict[str, object]], aggregate_rows: list[dict[str, str]]) -> None:
    manifest = {
        "experiment_id": "3d_pos_crossover_scaling",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [
            str(cfg.POS_ALL_DIR.relative_to(cfg.ROOT)),
            str(cfg.POS_MANUAL_DIR.relative_to(cfg.ROOT)),
        ],
        "outputs": {
            "pos_scaling_points.csv": {
                "rows": len(pos_rows),
                "schema": list(pos_rows[0].keys()) if pos_rows else [],
            },
            "manual_validation.csv": {
                "rows": len(manual_rows),
                "schema": list(manual_rows[0].keys()) if manual_rows else [],
            },
            "aggregate_statistics.csv": {
                "rows": len(aggregate_rows),
                "schema": ["metric_name", "value", "display_format", "notes"],
            },
        },
        "audit_flags": [
            "Manual validation has four alpha points but only three top300 manual POS CSV files.",
            "Claim map/manuscript mention max rank gap 15; available data gives 10.",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps({"source_bundles": manifest["source_bundles"]}, indent=2),
        encoding="utf-8",
    )


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pos_summary = load_json(cfg.POS_ALL_SUMMARY_JSON)
    pos_points = build_pos_points(read_csv(cfg.POS_ALL_POINTS_CSV))
    pos_points_by_slug = {str(row["slug"]): row for row in pos_points}
    manual_rows = build_manual_validation(read_csv(cfg.POS_MANUAL_POINTS_CSV), pos_points_by_slug)
    aggregate_rows = build_aggregate_rows(pos_summary, manual_rows)

    write_csv(
        output_dir / "pos_scaling_points.csv",
        pos_points,
        [
            "slug",
            "corpus",
            "vocabulary_size",
            "sqrt_v",
            "k_crossover",
            "alpha_per_corpus",
            "log_vocabulary_size",
            "log_k_crossover",
            "crossover_deviation_from_sqrt_v",
            "censored_at_top_n",
        ],
    )
    write_csv(
        output_dir / "manual_validation.csv",
        manual_rows,
        [
            "slug",
            "corpus",
            "source",
            "vocabulary_size",
            "manual_k_crossover",
            "manual_alpha_per_corpus",
            "spacy_k_crossover",
            "abs_manual_minus_spacy_rank_gap",
            "has_top300_manual_csv",
        ],
    )
    write_csv(output_dir / "aggregate_statistics.csv", aggregate_rows, ["metric_name", "value", "display_format", "notes"])
    write_manifest(output_dir, pos_points, manual_rows, aggregate_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()
