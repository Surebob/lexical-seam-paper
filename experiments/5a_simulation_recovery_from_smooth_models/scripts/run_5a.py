"""Migrate historical simulation-recovery outputs into experiment 5a."""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))
from source_config import EXPERIMENT_DIR, EXPERIMENT_ID, OUTPUTS, SOURCE_FILES  # noqa: E402


IS_WINNER = "sub[sub[x,1],log[x]]"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None and rows:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames or [])
        writer.writeheader()
        writer.writerows(rows)


def metric(name: str, value: object, display: str, notes: str) -> dict[str, object]:
    return {"metric_name": name, "value": value, "display_format": display, "notes": notes}


def main() -> None:
    summary = json.loads(SOURCE_FILES["summary"].read_text())
    source_rows = read_csv(SOURCE_FILES["table"])

    per_rows: list[dict[str, object]] = []
    for row in source_rows:
        per_rows.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "empirical_winner": row["empirical_winner"],
                "empirical_zm_c": row["empirical_zm_c"],
                "empirical_linear": row["empirical_linear"],
                "empirical_quadratic": row["empirical_quadratic"],
                "empirical_cubic": row["empirical_cubic"],
                "smooth_modal_winner": row["smooth_modal_winner"],
                "smooth_modal_share": row["smooth_modal_share"],
                "smooth_exact_match_rate": row["smooth_exact_match_rate"],
                "smooth_help_rate": row["smooth_help_rate"],
                "smooth_median_gap_full": row["smooth_median_gap_full"],
                "smooth_median_gap_top100": row["smooth_median_gap_top100"],
                "single_zm_control_modal_winner": row["single_zm_modal_winner"],
                "single_zm_control_modal_share": row["single_zm_modal_share"],
                "single_zm_control_exact_match_rate": row["single_zm_exact_match_rate"],
                "single_zm_control_help_rate": row["single_zm_help_rate"],
                "single_zm_control_median_gap_full": row["single_zm_median_gap_full"],
                "single_zm_control_median_gap_top100": row["single_zm_median_gap_top100"],
            }
        )

    high_rows = [row for row in source_rows if row["empirical_winner"] == IS_WINNER]
    single_help_count = sum(1 for row in source_rows if float(row["single_zm_help_rate"]) > 0)
    high_exact_count = sum(1 for row in high_rows if float(row["smooth_exact_match_rate"]) == 1.0)

    s = summary["summary"]
    smooth = s["smooth"]
    single = s["single_zm"]
    aggregate_rows = [
        metric("english_corpus_count", s["n_corpora"], "integer", "Source corpus count."),
        metric("simulation_replicates_per_corpus", s["n_replicates"], "integer", "Historical replicate count per corpus."),
        metric("empirical_median_gap_full", s["empirical_median_gap_full"], "decimal", "Empirical winner-vs-Euclidean full-RMSE gap."),
        metric("empirical_median_gap_top100", s["empirical_median_gap_top100"], "decimal", "Empirical winner-vs-Euclidean top-100 RMSE gap."),
        metric("smooth_exact_winner_match_rate", smooth["replicate_exact_winner_match_rate"], "proportion", "Replicate-level exact winner match rate for smooth-generated synthetic corpora."),
        metric("smooth_step2_help_rate", smooth["replicate_step2_help_rate"], "proportion", "Replicate-level step-2 help rate for smooth-generated synthetic corpora."),
        metric("smooth_majority_winner_match_count", smooth["majority_winner_match_count"], "X/25 count", "Corpus-level majority winner match count for smooth-generated synthetic corpora."),
        metric("smooth_median_winner_vs_euclidean_gap_full", smooth["replicate_median_winner_vs_euclidean_gap_full"], "decimal", "Replicate-level median gap, smooth synthetic, full RMSE."),
        metric("smooth_median_winner_vs_euclidean_gap_top100", smooth["replicate_median_winner_vs_euclidean_gap_top100"], "decimal", "Replicate-level median gap, smooth synthetic, top-100 RMSE."),
        metric("single_zm_control_exact_winner_match_rate", single["replicate_exact_winner_match_rate"], "proportion", "Replicate-level exact winner match rate for single-ZM control synthetic corpora."),
        metric("single_zm_control_step2_help_rate", single["replicate_step2_help_rate"], "proportion", "Replicate-level step-2 help rate for single-ZM control synthetic corpora."),
        metric("single_zm_control_step2_help_count", single_help_count, "X/25 count", "Corpus count with any positive single-ZM-control help rate."),
        metric("single_zm_control_majority_winner_match_count", single["majority_winner_match_count"], "X/25 count", "Corpus-level majority winner match count for single-ZM controls."),
        metric("single_zm_control_median_winner_vs_euclidean_gap_full", single["replicate_median_winner_vs_euclidean_gap_full"], "decimal", "Replicate-level median gap, single-ZM control, full RMSE."),
        metric("single_zm_control_median_winner_vs_euclidean_gap_top100", single["replicate_median_winner_vs_euclidean_gap_top100"], "decimal", "Replicate-level median gap, single-ZM control, top-100 RMSE."),
        metric("basis_corr_linear_smooth", smooth["basis_correlations"]["linear_u"], "correlation", "Correlation of empirical and smooth-generated linear head-basis coefficient."),
        metric("basis_corr_quadratic_smooth", smooth["basis_correlations"]["quadratic_u2"], "correlation", "Correlation of empirical and smooth-generated quadratic head-basis coefficient."),
        metric("basis_corr_cubic_smooth", smooth["basis_correlations"]["cubic_u3"], "correlation", "Correlation of empirical and smooth-generated cubic head-basis coefficient."),
        metric("basis_corr_linear_zm", single["basis_correlations"]["linear_u"], "correlation", "Correlation of empirical and single-ZM-control linear head-basis coefficient."),
        metric("basis_corr_quadratic_zm", single["basis_correlations"]["quadratic_u2"], "correlation", "Correlation of empirical and single-ZM-control quadratic head-basis coefficient."),
        metric("basis_corr_cubic_zm", single["basis_correlations"]["cubic_u3"], "correlation", "Correlation of empirical and single-ZM-control cubic head-basis coefficient."),
        metric("high_c_total_count", len(high_rows), "X/25 count", "Empirical IS-winner high-c block size."),
        metric("high_c_block_size", len(high_rows), "X/25 count", "Alias for high-c IS-winner block size."),
        metric("high_c_exact_match_count", high_exact_count, "X/25 count", "Smooth-generated exact winner matches in high-c IS block."),
        metric("high_c_exact_match_rate_smooth", high_exact_count / len(high_rows), "proportion", "Smooth-generated exact winner match rate in high-c IS block."),
        metric("high_c_help_rate_smooth", sum(float(row["smooth_help_rate"]) for row in high_rows) / len(high_rows), "proportion", "Smooth-generated step-2 help rate in high-c IS block."),
        metric("high_c_help_rate_zm_control", sum(float(row["single_zm_help_rate"]) for row in high_rows) / len(high_rows), "proportion", "Single-ZM-control step-2 help rate in high-c IS block."),
    ]

    write_csv(OUTPUTS["per_corpus"], per_rows)
    write_csv(OUTPUTS["aggregate"], aggregate_rows, ["metric_name", "value", "display_format", "notes"])

    archive = EXPERIMENT_DIR / "archive" / "provenance"
    archive.mkdir(parents=True, exist_ok=True)
    for key, src in SOURCE_FILES.items():
        shutil.copy2(src, archive / src.name)

    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "status": "complete_from_historical_bundle",
        "source_bundle": str(SOURCE_FILES["summary"].parents[0].relative_to(EXPERIMENT_DIR.parents[1])),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "simulation_recovery_per_corpus.csv": {"rows": len(per_rows)},
            "aggregate_statistics.csv": {"rows": len(aggregate_rows)},
        },
        "audit": "Historical smooth-generation source is preserved as canonical Phase 2 migration data; if strict decoupled-erf recovery is required, this experiment should be rerun rather than inferred.",
    }
    OUTPUTS["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote 5a outputs: {len(per_rows)} per-corpus rows, {len(aggregate_rows)} aggregate rows")


if __name__ == "__main__":
    main()

