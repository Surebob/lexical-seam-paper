"""Migrate low-c manifold and phase-coordinate outputs into experiment 5b."""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))
from source_config import EXPERIMENT_DIR, EXPERIMENT_ID, OUTPUTS, SOURCE_BUNDLES  # noqa: E402


EXPR_TO_FAMILY = {
    "sub[sub[x,1],log[x]]": "is",
    "eml[sub[x,1],eml[x,1]]": "exp",
    "sub[sqrt[x],pow[x,x]]": "xpow",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


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


def family(expr: str) -> str:
    return EXPR_TO_FAMILY.get(expr, expr)


def main() -> None:
    lowc = load_json(SOURCE_BUNDLES["lowc_manifold"] / "summary.json")
    phase = load_json(SOURCE_BUNDLES["phase_coordinate"] / "summary.json")
    sim_rows = read_csv(SOURCE_BUNDLES["simulation_recovery"] / "simulation_recovery_table.csv")

    low_rows = []
    for row in lowc["rows"]:
        low_rows.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "single_zm_c": row["c"],
                "full_winner": row["full_winner"],
                "top50_winner": row["top50_winner"],
                "top100_winner": row["top100_winner"],
                "top200_winner": row["top200_winner"],
                "cos_exp_xpow_top200": row["cos_exp_xpow_top200"],
                "r2_span_exp_xpow_top200": row["r2_span_expxpow_top200"],
                "r2_span_exp_is_top200": row["r2_span_expis_top200"],
                "xpow_minus_exp_full_rmse": row["xpow_minus_exp_full"],
                "xpow_minus_exp_top100_rmse": row["xpow_minus_exp_top100"],
                "full_scores_json": json.dumps(row["full_scores"], separators=(",", ":")),
                "top100_scores_json": json.dumps(row["top100_scores"], separators=(",", ":")),
            }
        )

    phase_rows = []
    for row in phase["rows"]:
        sweep_rows = row["sweep"]["rows"]
        lambda_1 = next((item for item in sweep_rows if abs(item["lambda"] - 1.0) < 1e-12), sweep_rows[-1])
        phase_rows.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "single_zm_c": row["c"],
                "transition_fraction": row["transition_fraction"],
                "empirical_winner": row["empirical_winner"],
                "theta_deg": row["angle"]["theta_deg"],
                "angle_r2": row["angle"]["r2"],
                "beta_exp": row["angle"]["beta_exp"],
                "beta_xpow": row["angle"]["beta_xpow"],
                "winner_full": row["winner_full"],
                "winner_top100": row["winner_top100"],
                "flip_lambda_exp_to_xpow": row.get("flip_lambda_exp_to_xpow", ""),
                "scores_lambda_1_json": json.dumps(lambda_1["scores"], separators=(",", ":")),
            }
        )

    phase_top100 = {row["slug"]: row["winner_top100"] for row in phase["low_rows"]}
    low_sim = [row for row in sim_rows if row["slug"] in phase_top100]
    smooth_matches = [family(row["smooth_modal_winner"]) == phase_top100[row["slug"]] for row in low_sim]
    zm_matches = [family(row["single_zm_modal_winner"]) == phase_top100[row["slug"]] for row in low_sim]

    summary = lowc["summary"]
    phase_summary = phase["summary"]
    aggregate_rows = [
        metric("lowc_corpus_count", summary["n_lowc_corpora"], "integer", "Low-c empirical exp-winner corpus count."),
        metric("lowc_median_cosine_exp_vs_xpow_top200", summary["median_cos_exp_xpow_top200"], "decimal", "Median cosine between exp and x^x-sqrt(x) on top-200 head coordinates."),
        metric("lowc_median_span_r2_exp_xpow_top200", summary["median_r2_span_expxpow_top200"], "proportion", "Median top-200 R^2 for span{exp,xpow}."),
        metric("lowc_median_span_r2_exp_is_top200", summary["median_r2_span_expis_top200"], "proportion", "Median top-200 R^2 for span{exp,IS}."),
        metric("lowc_median_delta_xpow_minus_exp_full_rmse", summary["median_xpow_minus_exp_full"], "signed decimal", "Median xpow-exp RMSE difference on full curve."),
        metric("lowc_median_delta_xpow_minus_exp_top100_rmse", summary["median_xpow_minus_exp_top100"], "signed decimal", "Median xpow-exp RMSE difference on top-100 window."),
        metric("lowc_full_rmse_exp_winner_count", summary["full_winner_counts"].get("exp", 0), "X/14 count", "Low-c full-RMSE exp winner count."),
        metric("lowc_top50_xpow_winner_count", summary["top50_winner_counts"].get("xpow", 0), "X/14 count", "Low-c top-50 xpow winner count."),
        metric("lowc_top100_xpow_winner_count", summary["top100_winner_counts"].get("xpow", 0), "X/14 count", "Low-c top-100 xpow winner count."),
        metric("lowc_top200_xpow_winner_count", summary["top200_winner_counts"].get("xpow", 0), "X/14 count", "Low-c top-200 xpow winner count."),
        metric("lowc_smooth_modal_vs_empirical_top100_match_count", sum(smooth_matches), "X/14 count", "Smooth synthetic modal winner matches empirical top-100 winner."),
        metric("lowc_smooth_modal_vs_empirical_top100_match_rate", sum(smooth_matches) / len(smooth_matches), "proportion", "Smooth synthetic modal winner match rate against empirical top-100 winner."),
        metric("lowc_zm_modal_vs_empirical_top100_match_count", sum(zm_matches), "X/14 count", "Single-ZM control modal winner matches empirical top-100 winner."),
        metric("lowc_zm_modal_vs_empirical_top100_match_rate", sum(zm_matches) / len(zm_matches), "proportion", "Single-ZM control modal winner match rate against empirical top-100 winner."),
        metric("phase_high_full_is_winner_count", phase_summary["high_full_counts"].get("is", 0), "X/11 count", "High-c full-RMSE IS winner count."),
        metric("phase_high_top100_is_winner_count", phase_summary["high_top100_counts"].get("is", 0), "X/11 count", "High-c top-100 IS winner count."),
        metric("phase_corr_theta_c", phase_summary["corr_theta_c"], "correlation", "Correlation between theta and single-ZM c."),
        metric("phase_corr_theta_transition_fraction", phase_summary["corr_theta_transition_fraction"], "correlation", "Correlation between theta and transition fraction."),
        metric("phase_corr_theta_lambda_flip", phase_summary["corr_theta_lambda_flip"], "correlation", "Correlation between theta and lambda flip."),
        metric("lowc_median_theta_deg", median(row["angle"]["theta_deg"] for row in phase["low_rows"]), "degrees", "Median low-c phase angle."),
    ]

    write_csv(OUTPUTS["lowc_per_corpus"], low_rows)
    write_csv(OUTPUTS["phase_per_corpus"], phase_rows)
    write_csv(OUTPUTS["aggregate"], aggregate_rows, ["metric_name", "value", "display_format", "notes"])

    archive = EXPERIMENT_DIR / "archive" / "provenance"
    archive.mkdir(parents=True, exist_ok=True)
    for bundle in SOURCE_BUNDLES.values():
        for src in [bundle / "summary.json", bundle / "report.md"]:
            if src.exists():
                shutil.copy2(src, archive / f"{bundle.name}_{src.name}")

    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "status": "complete_from_historical_bundles",
        "source_bundles": {name: str(path.relative_to(EXPERIMENT_DIR.parents[1])) for name, path in SOURCE_BUNDLES.items()},
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "lowc_manifold_per_corpus.csv": {"rows": len(low_rows)},
            "phase_coordinate_per_corpus.csv": {"rows": len(phase_rows)},
            "aggregate_statistics.csv": {"rows": len(aggregate_rows)},
        },
    }
    OUTPUTS["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote 5b outputs: {len(low_rows)} low-c rows, {len(phase_rows)} phase rows, {len(aggregate_rows)} aggregate rows")


if __name__ == "__main__":
    main()

