"""Reconstruct legacy coupled-logistic seam-theory archive outputs.

This is a migration/consolidation script only. It reads the original historical
summary files and emits first-class archive tables under the canonical
experiment directory. It does not modify historical results and does not run
fresh seam-theory computation.
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def find_repo_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / "results").is_dir() and (candidate / "experiments").is_dir():
            return candidate
    raise RuntimeError(f"Could not find repository root from {start}")


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = find_repo_root(SCRIPT_PATH)
ARCHIVE_DIR = SCRIPT_PATH.parents[1]
OUTPUT_DIR = ARCHIVE_DIR / "outputs"

SOURCE_BUNDLES = {
    "sign": REPO_ROOT / "results" / "zipf_seam_sign_theory" / "summary.json",
    "projection": REPO_ROOT / "results" / "zipf_seam_projection_theory" / "summary.json",
    "second_order": REPO_ROOT / "results" / "zipf_seam_second_order_theory" / "summary.json",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def as_json_cell(value: Any) -> str:
    return json.dumps(value, separators=(",", ":"))


def bools_all(values: list[bool]) -> bool:
    return all(bool(v) for v in values)


def metric_row(metric_name: str, value: Any, notes: str, display_format: str = "X/25 count") -> dict[str, Any]:
    return {
        "metric_name": metric_name,
        "value": value,
        "display_format": display_format,
        "notes": notes,
    }


def build_per_corpus_rows(sign: dict[str, Any], proj: dict[str, Any], second: dict[str, Any]) -> list[dict[str, Any]]:
    sign_rows = {row["slug"]: row for row in sign["rows"]}
    proj_rows = {row["slug"]: row for row in proj["rows"]}
    second_rows = {row["slug"]: row for row in second["rows"]}
    slugs = sorted(sign_rows)
    if set(slugs) != set(proj_rows) or set(slugs) != set(second_rows):
        raise ValueError("Source seam-theory bundles do not contain matching corpus slug sets.")

    rows: list[dict[str, Any]] = []
    for slug in slugs:
        s = sign_rows[slug]
        p = proj_rows[slug]
        q = second_rows[slug]
        rows.append(
            {
                "slug": slug,
                "corpus": s["corpus"],
                "family": s["family"],
                "single_zm_c": s["c"],
                "local_signs": as_json_cell(s["local_signs"]),
                "empirical_signs_head10": as_json_cell(s["empirical_signs_10"]),
                "empirical_signs_head50": as_json_cell(s["empirical_signs_50"]),
                "smooth_signs_head50": as_json_cell(s["smooth_signs_50"]),
                "seam_signs_head50": as_json_cell(s["seam_signs_50"]),
                "projected_signs_head50": as_json_cell(p["projected_signs"]),
                "second_order_signs_head50": as_json_cell(q["second_order_signs"]),
                "quadratic_surrogate_signs_head50": as_json_cell(q["quadratic_surrogate_signs"]),
                "match_local_vs_emp10": as_json_cell(s["match_local_vs_emp10"]),
                "match_local_vs_emp50": as_json_cell(s["match_local_vs_emp50"]),
                "match_local_vs_emp10_all_terms": bools_all(s["match_local_vs_emp10"]),
                "match_local_vs_emp50_all_terms": bools_all(s["match_local_vs_emp50"]),
                "match_smooth_vs_emp50": as_json_cell(s["match_smooth_vs_emp50"]),
                "match_smooth_vs_emp50_all_terms": bools_all(s["match_smooth_vs_emp50"]),
                "match_projected_vs_emp50": as_json_cell(p["match_proj_vs_emp"]),
                "match_projected_vs_emp50_all_terms": bools_all(p["match_proj_vs_emp"]),
                "match_projected_vs_smooth50": as_json_cell(p["match_proj_vs_smooth"]),
                "match_projected_vs_smooth50_all_terms": bools_all(p["match_proj_vs_smooth"]),
                "match_second_order_vs_emp50": as_json_cell(q["match_second_vs_emp"]),
                "match_second_order_vs_emp50_all_terms": bools_all(q["match_second_vs_emp"]),
                "match_second_order_vs_smooth50": as_json_cell(q["match_second_vs_smooth"]),
                "match_second_order_vs_smooth50_all_terms": bools_all(q["match_second_vs_smooth"]),
                "match_quadratic_surrogate_vs_emp50": as_json_cell(q["match_quad_vs_emp"]),
                "match_quadratic_surrogate_vs_emp50_all_terms": bools_all(q["match_quad_vs_emp"]),
                "match_quadratic_surrogate_vs_smooth50": as_json_cell(q["match_quad_vs_smooth"]),
                "match_quadratic_surrogate_vs_smooth50_all_terms": bools_all(q["match_quad_vs_smooth"]),
                "source_model": "legacy_coupled_logistic",
            }
        )
    return rows


def build_aggregate_rows(sign: dict[str, Any], proj: dict[str, Any], second: dict[str, Any]) -> list[dict[str, Any]]:
    s = sign["summary"]
    p = proj["summary"]
    q = second["summary"]
    rows: list[dict[str, Any]] = []
    note = "Legacy coupled-logistic seam-theory aggregate; preserved as prior canonical research output."

    rows.extend(
        [
            metric_row("raw_sign_prediction_match_count_internal", s["local_all_negative_count"], note),
            metric_row("internal_predicted_sign_match_count", s["local_all_negative_count"], note),
            metric_row("raw_seam_term_match_head10_linear", s["local_match_counts_head10"][0], note),
            metric_row("raw_seam_term_match_head10_quadratic", s["local_match_counts_head10"][1], note),
            metric_row("raw_seam_term_match_head10_cubic", s["local_match_counts_head10"][2], note),
            metric_row("raw_seam_term_match_head50_linear", s["local_match_counts_head50"][0], note),
            metric_row("raw_seam_term_match_head50_quadratic", s["local_match_counts_head50"][1], note),
            metric_row("raw_seam_term_match_head50_cubic", s["local_match_counts_head50"][2], note),
            metric_row("raw_seam_full_sign_match_head10", s["local_full_match_head10"], note),
            metric_row("raw_seam_full_sign_match_head50", s["local_full_match_head50"], note),
            metric_row("numerical_refit_match_head50_linear", s["smooth_match_counts_head50"][0], note),
            metric_row("numerical_refit_match_head50_quadratic", s["smooth_match_counts_head50"][1], note),
            metric_row("numerical_refit_match_head50_cubic", s["smooth_match_counts_head50"][2], note),
            metric_row("numerical_refit_full_sign_match_head50", s["smooth_full_match_head50"], note),
            metric_row("projected_vs_empirical_match_head50_linear", p["projected_head50_term_matches"][0], note),
            metric_row("projected_vs_empirical_match_head50_quadratic", p["projected_head50_term_matches"][1], note),
            metric_row("projected_vs_empirical_match_head50_cubic", p["projected_head50_term_matches"][2], note),
            metric_row("projected_vs_empirical_full_match_head50", p["projected_head50_full_matches"], note),
            metric_row("projected_vs_empirical_full_match_count", p["projected_head50_full_matches"], note),
            metric_row("projected_vs_numerical_match_head50_linear", p["projected_vs_exact_smooth_term_matches"][0], note),
            metric_row("projected_vs_numerical_match_head50_quadratic", p["projected_vs_exact_smooth_term_matches"][1], note),
            metric_row("projected_vs_numerical_match_head50_cubic", p["projected_vs_exact_smooth_term_matches"][2], note),
            metric_row("projected_vs_numerical_full_match_head50", p["projected_vs_exact_smooth_full_matches"], note),
            metric_row("projected_vs_numerical_refit_full_match_count", p["projected_vs_exact_smooth_full_matches"], note),
            metric_row("second_order_match_head50_linear", q["second_order_head50_term_matches"][0], note),
            metric_row("second_order_match_head50_quadratic", q["second_order_head50_term_matches"][1], note),
            metric_row("second_order_match_head50_cubic", q["second_order_head50_term_matches"][2], note),
            metric_row("second_order_full_sign_match_head50", q["second_order_head50_full_matches"], note),
        ]
    )
    return rows


def build_consistency_rows(sign: dict[str, Any], proj: dict[str, Any], second: dict[str, Any]) -> list[dict[str, Any]]:
    s = sign["summary"]
    p = proj["summary"]
    q = second["summary"]
    checks = [
        ("local_head50_term_matches", s["local_match_counts_head50"], p["local_baseline_head50_term_matches"], q["local_baseline_head50_term_matches"]),
        ("local_head50_full_matches", s["local_full_match_head50"], p["local_baseline_head50_full_matches"], q["local_baseline_head50_full_matches"]),
        ("smooth_head50_term_matches", s["smooth_match_counts_head50"], p["smooth_exact_head50_term_matches"], q["exact_smooth_head50_term_matches"]),
        ("smooth_head50_full_matches", s["smooth_full_match_head50"], p["smooth_exact_head50_full_matches"], q["exact_smooth_head50_full_matches"]),
        ("projected_head50_term_matches", p["projected_head50_term_matches"], q["first_order_head50_term_matches"]),
        ("projected_head50_full_matches", p["projected_head50_full_matches"], q["first_order_head50_full_matches"]),
    ]
    rows = []
    for item in checks:
        name, *values = item
        rows.append(
            {
                "check_name": name,
                "values": as_json_cell(values),
                "consistent": all(value == values[0] for value in values[1:]),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    sign = load_json(SOURCE_BUNDLES["sign"])
    proj = load_json(SOURCE_BUNDLES["projection"])
    second = load_json(SOURCE_BUNDLES["second_order"])

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    per_corpus = build_per_corpus_rows(sign, proj, second)
    aggregates = build_aggregate_rows(sign, proj, second)
    consistency = build_consistency_rows(sign, proj, second)

    write_csv(OUTPUT_DIR / "seam_sign_checks.csv", per_corpus)
    write_csv(OUTPUT_DIR / "aggregate_statistics.csv", aggregates)
    write_csv(OUTPUT_DIR / "source_consistency_checks.csv", consistency)

    manifest = {
        "archive_id": "legacy_coupled_logistic",
        "parent_experiment": "4_analytic_seam_theory",
        "status": "first_class_legacy_prior_canon",
        "model": "coupled_logistic",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundles": {name: str(path.relative_to(REPO_ROOT)) for name, path in SOURCE_BUNDLES.items()},
        "outputs": {
            "seam_sign_checks.csv": {"rows": len(per_corpus)},
            "aggregate_statistics.csv": {"rows": len(aggregates)},
            "source_consistency_checks.csv": {"rows": len(consistency)},
        },
        "notes": "Archive means superseded by later canonical choice, not low-quality or invalid data.",
    }
    (OUTPUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Wrote {len(per_corpus)} per-corpus rows and {len(aggregates)} aggregate rows to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

