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


AGG_FIELDS = ["metric_name", "value", "display_format", "notes"]


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def metric(name: str, value: object, fmt: str, notes: str) -> dict[str, object]:
    return {"metric_name": name, "value": value, "display_format": fmt, "notes": notes}


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    bic_obj = load_json(cfg.BIC_COMPARISON_SUMMARY)
    moe_obj = load_json(cfg.MOEZIPF_COMPARISON_SUMMARY)
    cont_obj = load_json(cfg.CONTINUOUS_PIECEWISE_SUMMARY)
    sqrt_obj = load_json(cfg.SQRT_V_SUMMARY)

    bic_rows = {row["slug"]: row for row in bic_obj["rows"]}
    moe_rows = {row["slug"]: row for row in moe_obj["rows"]}
    cont_rows = {row["slug"]: row for row in cont_obj["rows"]}
    sqrt_rows = {row["slug"]: row for row in sqrt_obj["rows"]}

    table_rows: list[dict[str, object]] = []
    rmse_rows: list[dict[str, object]] = []
    winner_counts = {
        "single_zm": 0,
        "moezipf": 0,
        "hard_piecewise": 0,
        "continuous_piecewise": 0,
        "smooth_ksqrtv": 0,
        "smooth_freek": 0,
    }

    for slug in sorted(bic_rows):
        bic = bic_rows[slug]
        moe = moe_rows[slug]
        cont = cont_rows[slug]
        sqrt = sqrt_rows[slug]
        bics = {
            "single_zm": float(bic["single_zm_bic"]),
            "moezipf": float(cont["moezipf_bic"]),
            "hard_piecewise": float(bic["piecewise_k500_bic"]),
            "continuous_piecewise": float(cont["continuous_piecewise_bic"]),
            "smooth_ksqrtv": float(bic["reranked_7param_sqrtv_bic"]),
            "smooth_freek": float(bic["reranked_8param_bic"]),
        }
        winner = min(bics, key=bics.get)
        winner_counts[winner] += 1

        table_rows.append(
            {
                "slug": slug,
                "corpus": bic["name"],
                "vocabulary_size": bic["vocab_size"],
                "single_zm_bic": bics["single_zm"],
                "moe_bic": bics["moezipf"],
                "hard_piecewise_bic": bics["hard_piecewise"],
                "continuous_piecewise_bic": bics["continuous_piecewise"],
                "smooth_ksqrtv_bic": bics["smooth_ksqrtv"],
                "smooth_freek_bic": bics["smooth_freek"],
                "winner_family": winner,
                "legacy_v4_zipf_rank_bic": cont["zipf_bic"],
                "legacy_v4_zipf_rank_bic_delta_vs_single_zm_bic": float(cont["zipf_bic"]) - bics["single_zm"],
            }
        )
        rmse_rows.append(
            {
                "slug": slug,
                "corpus": bic["name"],
                "vocabulary_size": bic["vocab_size"],
                "zipf_rank_rmse_legacy": moe["zipf_rmse"],
                "moezipf_rank_rmse": moe["moe_rmse"],
                "single_zm_rank_rmse": sqrt["single_zm_rmse"],
                "hard_piecewise_rank_rmse": "",
                "continuous_piecewise_rank_rmse": cont["continuous_piecewise_rmse"],
                "smooth_ksqrtv_rank_rmse": sqrt["sqrt_v_rmse"],
                "smooth_free_k_rank_rmse": sqrt["reranked_8param_rmse"],
                "sqrtv_minus_free_k_rmse": sqrt["delta_vs_8param"],
            }
        )

    smooth_total = winner_counts["smooth_ksqrtv"] + winner_counts["smooth_freek"]
    piecewise_loses_count = 0
    for row in table_rows:
        smooth_best = min(float(row["smooth_ksqrtv_bic"]), float(row["smooth_freek_bic"]))
        if float(row["hard_piecewise_bic"]) > smooth_best and float(row["continuous_piecewise_bic"]) > smooth_best:
            piecewise_loses_count += 1

    sqrt_gaps = [float(row["delta_vs_8param"]) for row in sqrt_rows.values()]
    aggregate_rows = [
        metric(
            "piecewise_loses_to_smooth_on_all_25",
            piecewise_loses_count,
            "X/25 count",
            "Count of corpora where both hard and continuous piecewise BIC are worse than the better smooth family.",
        ),
        metric(
            "winner_count_smooth_family_total",
            smooth_total,
            "X/25 count",
            "BIC winner count for smooth families combined under true single-ZM baseline convention.",
        ),
        metric(
            "moezipf_bic_winner_count",
            winner_counts["moezipf"],
            "X/25 count",
            "MOEZipf winner count after replacing the v4 mislabeled pure-Zipf baseline with true single-ZM BIC.",
        ),
        metric(
            "winner_count_moezipf",
            winner_counts["moezipf"],
            "X/25 count",
            "Alias for manuscript discussion restatements.",
        ),
        metric(
            "median_rmse_cost_sqrtv_vs_free_k",
            statistics.median(sqrt_gaps),
            "float",
            "Median RMSE cost of fixed k=sqrt(V) relative to free-k smooth fit.",
        ),
    ]

    table_fields = [
        "slug",
        "corpus",
        "vocabulary_size",
        "single_zm_bic",
        "moe_bic",
        "hard_piecewise_bic",
        "continuous_piecewise_bic",
        "smooth_ksqrtv_bic",
        "smooth_freek_bic",
        "winner_family",
        "legacy_v4_zipf_rank_bic",
        "legacy_v4_zipf_rank_bic_delta_vs_single_zm_bic",
    ]
    rmse_fields = [
        "slug",
        "corpus",
        "vocabulary_size",
        "zipf_rank_rmse_legacy",
        "moezipf_rank_rmse",
        "single_zm_rank_rmse",
        "hard_piecewise_rank_rmse",
        "continuous_piecewise_rank_rmse",
        "smooth_ksqrtv_rank_rmse",
        "smooth_free_k_rank_rmse",
        "sqrtv_minus_free_k_rmse",
    ]
    write_csv(output_dir / "table2_model_family.csv", table_rows, table_fields)
    write_csv(output_dir / "model_family_rmse_per_corpus.csv", rmse_rows, rmse_fields)
    write_csv(output_dir / "aggregate_statistics.csv", aggregate_rows, AGG_FIELDS)

    manifest = {
        "experiment_id": "3b_model_family_bic_comparison",
        "migration_type": "historical_bundle_consolidation",
        "baseline_convention": "true_3_parameter_single_zm_bic",
        "source_bundles": [str(path.relative_to(cfg.ROOT)) for path in cfg.SOURCE_BUNDLES],
        "outputs": {
            "table2_model_family.csv": {"rows": len(table_rows), "schema": table_fields},
            "model_family_rmse_per_corpus.csv": {"rows": len(rmse_rows), "schema": rmse_fields},
            "aggregate_statistics.csv": {"rows": len(aggregate_rows), "schema": AGG_FIELDS},
        },
        "audit_flags": [
            "v4 Table 2 rendered zipf_continuous_piecewise.zipf_bic under a ZM label; canonical 3b corrects this to true single_zm_bic.",
            "v5.1 Table 2 no longer uses historical 3b; it uses S2 v3 decoupled gate-family outputs.",
            "All manuscript citations of single-ZM BIC must be re-verified against table2_model_family.csv.",
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
