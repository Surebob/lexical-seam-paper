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


FIELDS = [
    "ablation_family",
    "slug",
    "corpus",
    "setting",
    "x_low",
    "exp_clamp",
    "value_abs_limit",
    "weighting",
    "token_count",
    "vocabulary_size",
    "zm_c",
    "zm_rmse",
    "step2_winner_expr",
    "step2_winner_math",
    "step2_rmse",
    "step2_delta_vs_zm",
    "step2_helpful",
    "canonical_step2_winner_expr",
    "matches_canonical_step2_winner",
    "notes",
]


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def load_canonical_winners() -> dict[str, str]:
    out: dict[str, str] = {}
    with cfg.ONE_A_TABLE.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            out[row["corpus"]] = row["step2_winner_expression"]
    aliases = {
        "Shakespeare": "Complete Works of Shakespeare",
        "bible": "King James Bible",
    }
    for alias, canonical in aliases.items():
        if canonical in out:
            out[alias] = out[canonical]
    return out


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    canonical = load_canonical_winners()
    rows: list[dict[str, object]] = []

    for item in json.loads(cfg.BOUNDARY_SUMMARY.read_text(encoding="utf-8")):
        canonical_expr = canonical.get(item["corpus"], "")
        rows.append(
            {
                "ablation_family": "boundary_x_low",
                "slug": item["corpus_slug"],
                "corpus": item["corpus"],
                "setting": f"x_low={item['x_low']}",
                "x_low": item["x_low"],
                "zm_c": item["c"],
                "zm_rmse": item["zm_baseline_rmse"],
                "step2_winner_expr": item["step2_expr"],
                "step2_winner_math": item["step2_math"],
                "step2_rmse": item["step2_rmse"],
                "step2_delta_vs_zm": float(item["step2_rmse"]) - float(item["zm_baseline_rmse"]),
                "step2_helpful": str(item["step2_helps"]),
                "canonical_step2_winner_expr": canonical_expr,
                "matches_canonical_step2_winner": str(item["step2_expr"] == canonical_expr) if canonical_expr else "",
                "notes": "Boundary-normalization x_low ablation.",
            }
        )

    for item in json.loads(cfg.GUARD_SUMMARY.read_text(encoding="utf-8")):
        canonical_expr = canonical.get(item["corpus_name"], "")
        rows.append(
            {
                "ablation_family": "numeric_guard",
                "slug": item["corpus_slug"],
                "corpus": item["corpus_name"],
                "setting": item["profile"],
                "exp_clamp": item["exp_clamp"],
                "value_abs_limit": item["value_abs_limit"],
                "zm_c": item["zm_c"],
                "zm_rmse": item["zm_rmse"],
                "step2_winner_expr": item["step2_expr"],
                "step2_winner_math": item["step2_math"],
                "step2_rmse": item["step2_rmse"],
                "step2_delta_vs_zm": float(item["step2_rmse"]) - float(item["zm_rmse"]),
                "step2_helpful": str(float(item["step2_rmse"]) < float(item["zm_rmse"])),
                "canonical_step2_winner_expr": canonical_expr,
                "matches_canonical_step2_winner": str(item["step2_expr"] == canonical_expr) if canonical_expr else "",
                "notes": "Exp-clamp and value-absolute-limit guard ablation.",
            }
        )

    for corpus in json.loads(cfg.WLS_SUMMARY.read_text(encoding="utf-8")):
        canonical_expr = canonical.get(corpus["name"], "")
        for method in corpus["methods"]:
            rows.append(
                {
                    "ablation_family": "weighted_least_squares",
                    "slug": corpus["slug"],
                    "corpus": corpus["name"],
                    "setting": method["slug"],
                    "weighting": method["name"],
                    "token_count": corpus["token_count"],
                    "vocabulary_size": corpus["unique_words"],
                    "zm_c": method["c"],
                    "zm_rmse": method["zm_rmse_unweighted"],
                    "step2_winner_expr": "sub[sub[x,1],log[x]]",
                    "step2_winner_math": "((x-1)-log(x))",
                    "step2_rmse": method["zm_plus_step2_rmse_unweighted"],
                    "step2_delta_vs_zm": method["step2_delta"],
                    "step2_helpful": str(float(method["step2_delta"]) < 0.0),
                    "canonical_step2_winner_expr": canonical_expr,
                    "matches_canonical_step2_winner": str(canonical_expr == "sub[sub[x,1],log[x]]") if canonical_expr else "",
                    "notes": "Historical WLS script evaluated the canonical IS correction; it did not save a fresh search winner per weighting method.",
                }
            )

    write_csv(output_dir / "robustness_ablation_by_corpus.csv", rows, FIELDS)
    write_csv(output_dir / "aggregate_statistics.csv", [], ["metric_name", "value", "display_format", "notes"])
    manifest = {
        "experiment_id": "1c_search_robustness_ablations",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(path.relative_to(cfg.ROOT)) for path in cfg.SOURCE_BUNDLES],
        "outputs": {
            "robustness_ablation_by_corpus.csv": {"rows": len(rows), "schema": FIELDS},
            "aggregate_statistics.csv": {"rows": 0, "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "audit_flags": [
            "Historical WLS output evaluates a fixed IS correction and does not store a freshly searched step-2 winner per weighted fit.",
            "Ablation scope is diagnostic subsets, not all 25 English corpora.",
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

