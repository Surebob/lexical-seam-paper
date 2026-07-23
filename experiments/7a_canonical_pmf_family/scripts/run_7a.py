"""Migrate canonical PMF-family outputs into experiment 7a.

The manuscript-facing Table 4 source is `zipf_v4_verification/table_a_fourway_pmf.csv`.
This is canonical-by-documentation rather than canonical-by-reproducible producer:
the table has no dedicated committed producer script, so this migration emits
provenance and source-diagnostic sidecars documenting the optimization-path
sensitivity behind the known three-row soft-k discrepancy.
"""

from __future__ import annotations

import csv
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any

import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))
from source_config import (  # noqa: E402
    CANONICAL_TABLE4_SOURCE,
    EXPERIMENT_DIR,
    EXPERIMENT_ID,
    OUTPUTS,
    SOFTK_PARAMETER_COUNT,
    SOURCE_BUNDLES,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None and rows:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames or [])
        writer.writeheader()
        writer.writerows(rows)


def f(value: Any) -> float:
    return float(value)


def metric(name: str, value: Any, display: str, notes: str) -> dict[str, Any]:
    return {"metric_name": name, "value": value, "display_format": display, "notes": notes}


def lambda_key(value: str) -> str:
    val = f(value)
    mapping = {
        1e-4: "lambda_count_1e-4",
        3e-4: "lambda_count_3e-4",
        1e-3: "lambda_count_1e-3",
        3e-3: "lambda_count_3e-3",
        1e-2: "lambda_count_1e-2",
        3e-2: "lambda_count_3e-2",
    }
    for k, name in mapping.items():
        if abs(val - k) < 1e-12:
            return name
    return f"lambda_count_{value}"


def build_table4(v4_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    out = []
    for row in v4_rows:
        out.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "token_count": int(row["token_count"]),
                "vocabulary_size": int(row["vocab_size"]),
                "zipf_test_avg_nll": f(row["zipf_test_avg_nll"]),
                "zm_test_avg_nll": f(row["zm_test_avg_nll"]),
                "moe_test_avg_nll": f(row["moe_test_avg_nll"]),
                "softk_test_avg_nll": f(row["softk_test_avg_nll"]),
                "softk_minus_moe": f(row["softk_minus_moe"]),
                "softk_minus_zm": f(row["softk_minus_zm"]),
                "softk_minus_zipf": f(row["softk_test_avg_nll"]) - f(row["zipf_test_avg_nll"]),
                "winner_family": row["winner"],
                "best_lambda_k": f(row["lambda_k"]),
                "canonical_source": "results/zipf_v4_verification/table_a_fourway_pmf.csv",
            }
        )
    return out


def build_provenance_rows(table4: list[dict[str, Any]], pmf: dict[str, dict[str, str]], split: dict[str, dict[str, str]], legacy: dict[str, dict[str, str]], softkw: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    rows = []
    for row in table4:
        slug = row["slug"]
        split_delta = row["softk_test_avg_nll"] - f(split[slug]["softk_test_avg_nll"])
        rows.append(
            {
                "slug": slug,
                "corpus": row["corpus"],
                "zipf_source": "zipf_seam_mandelbrot_pmf/seam_mandelbrot_table.csv",
                "zipf_matches_source": abs(row["zipf_test_avg_nll"] - f(pmf[slug]["zipf_test_avg_nll"])) < 1e-12,
                "zm_source": "zipf_seam_mandelbrot_pmf/seam_mandelbrot_table.csv",
                "zm_matches_source": abs(row["zm_test_avg_nll"] - f(pmf[slug]["zm_test_avg_nll"])) < 1e-12,
                "moe_source": "zipf_seam_mandelbrot_pmf/seam_mandelbrot_table.csv",
                "moe_matches_source": abs(row["moe_test_avg_nll"] - f(pmf[slug]["moe_test_avg_nll"])) < 1e-12,
                "softk_canonical_source": "zipf_v4_verification/table_a_fourway_pmf.csv",
                "softk_matches_legacy_softk": abs(row["softk_test_avg_nll"] - f(legacy[slug]["softk_test_avg_nll"])) < 1e-12,
                "softk_matches_softkw_repeated": abs(row["softk_test_avg_nll"] - f(softkw[slug]["softk_test_avg_nll"])) < 1e-12,
                "softk_matches_softk_splitfit": abs(split_delta) < 1e-12,
                "softk_splitfit_delta": split_delta,
                "provenance_note": "Canonical-by-documentation Table 4 row; soft-k local-optimum sensitivity documented in softk_source_diagnostic.csv.",
            }
        )
    return rows


def build_softk_diagnostic(table4: list[dict[str, Any]], split: dict[str, dict[str, str]], legacy: dict[str, dict[str, str]], softkw: dict[str, dict[str, str]], split_json: dict[str, dict], legacy_json: dict[str, dict], softkw_json: dict[str, dict]) -> list[dict[str, Any]]:
    rows = []
    for row in table4:
        slug = row["slug"]
        v4_val = row["softk_test_avg_nll"]
        split_val = f(split[slug]["softk_test_avg_nll"])
        legacy_val = f(legacy[slug]["softk_test_avg_nll"])
        softkw_repeated_val = f(softkw[slug]["softk_test_avg_nll"])
        rows.append(
            {
                "slug": slug,
                "corpus": row["corpus"],
                "v4_verification_softk": v4_val,
                "legacy_softk": legacy_val,
                "softkw_repeated_softk": softkw_repeated_val,
                "softk_splitfit": split_val,
                "delta_v4_minus_splitfit": v4_val - split_val,
                "lambda_k_v4": row["best_lambda_k"],
                "lambda_k_splitfit": f(split[slug]["best_lambda"]),
                "lambda_k_legacy": f(legacy[slug]["best_lambda"]),
                "legacy_k": legacy_json[slug]["best_lambda"]["params"]["k"],
                "legacy_w": legacy_json[slug]["best_lambda"]["params"]["w"],
                "splitfit_k": split_json[slug]["best_lambda"]["selection"]["params"]["k"],
                "splitfit_w": split_json[slug]["best_lambda"]["selection"]["params"]["w"],
                "softkw_inherits_legacy_softk": abs(softkw_repeated_val - legacy_val) < 1e-12,
                "is_known_local_optimum_discrepancy": abs(v4_val - split_val) > 1e-12,
                "diagnosis": "same_lambda_different_local_optimum" if abs(v4_val - split_val) > 1e-12 else "matches_splitfit",
            }
        )
    return rows


def build_variant_rows(table4: list[dict[str, Any]], pmf: dict[str, dict[str, str]], regularized: dict[str, dict[str, str]], split: dict[str, dict[str, str]], softkw: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    rows = []
    for row in table4:
        slug = row["slug"]
        p = pmf[slug]
        r = regularized[slug]
        s = split[slug]
        w = softkw[slug]
        rows.append(
            {
                "slug": slug,
                "corpus": row["corpus"],
                "free_test_avg_nll": f(p["seam_test_avg_nll"]),
                "fixedk_test_avg_nll": f(r["fixed_k_test_avg_nll"]),
                "fixedkw_test_avg_nll": f(r["fixed_kw_test_avg_nll"]),
                "softk_table4_test_avg_nll": row["softk_test_avg_nll"],
                "softk_splitfit_test_avg_nll": f(s["softk_test_avg_nll"]),
                "softkw_test_avg_nll": f(w["softkw_test_avg_nll"]),
                "softk_table4_minus_splitfit": row["softk_test_avg_nll"] - f(s["softk_test_avg_nll"]),
                "free_step2_gain": f(p["seam_step2_gain"]),
                "fixedk_step2_gain": f(r["fixed_k_step2_gain"]),
                "fixedkw_step2_gain": f(r["fixed_kw_step2_gain"]),
                "softk_splitfit_step2_gain": f(s["softk_step2_gain"]),
                "softkw_step2_gain": f(w["softkw_step2_gain"]),
                "best_lambda_k": row["best_lambda_k"],
                "best_lambda_w": f(w["best_lambda_w"]),
                "source_note": "Table4 soft-k is canonical-by-documentation; splitfit and softkw values are preserved for path-sensitivity audit.",
            }
        )
    return rows


def build_aggregate_rows(table4: list[dict[str, Any]], pmf_summary: dict, regularized_summary: dict, split_summary: dict) -> list[dict[str, Any]]:
    winners = Counter(row["winner_family"] for row in table4)
    lambda_counts = Counter(lambda_key(str(row["best_lambda_k"])) for row in table4)
    softk_minus_zipf = [row["softk_minus_zipf"] for row in table4]
    softk_minus_zm = [row["softk_minus_zm"] for row in table4]
    softk_minus_moe = [row["softk_minus_moe"] for row in table4]
    lambdas = [row["best_lambda_k"] for row in table4]

    rows = [
        metric("winner_count_zipf", winners.get("zipf", 0), "X/25 count", "Recomputed from canonical Table 4."),
        metric("winner_count_zm", winners.get("zm", 0), "X/25 count", "Recomputed from canonical Table 4."),
        metric("winner_count_moe", winners.get("moe", 0), "X/25 count", "Recomputed from canonical Table 4."),
        metric("winner_count_softk", winners.get("softk", 0), "X/25 count", "Recomputed from canonical Table 4."),
        metric("softk_beats_zipf_count", sum(row["softk_test_avg_nll"] < row["zipf_test_avg_nll"] for row in table4), "X/25 count", "Canonical Table 4 soft-k vs Zipf."),
        metric("softk_beats_zm_count", sum(row["softk_test_avg_nll"] < row["zm_test_avg_nll"] for row in table4), "X/25 count", "Canonical Table 4 soft-k vs ZM."),
        metric("softk_beats_moe_count", sum(row["softk_test_avg_nll"] < row["moe_test_avg_nll"] for row in table4), "X/25 count", "Canonical Table 4 soft-k vs MOEZipf."),
        metric("median_softk_minus_zipf", median(softk_minus_zipf), "signed decimal", "Median canonical Table 4 soft-k minus Zipf."),
        metric("median_softk_minus_zm", median(softk_minus_zm), "signed decimal", "Median canonical Table 4 soft-k minus ZM."),
        metric("median_softk_minus_moe", median(softk_minus_moe), "signed decimal", "Median canonical Table 4 soft-k minus MOEZipf."),
        metric("free_pmf_heldout_winner_count_zipf", pmf_summary["counts"]["test_winner_counts"].get("zipf", 0), "X/25 count", "Free PMF four-family held-out winner count."),
        metric("free_pmf_heldout_winner_count_zm", pmf_summary["counts"]["test_winner_counts"].get("zm", 0), "X/25 count", "Free PMF four-family held-out winner count."),
        metric("free_pmf_heldout_winner_count_moe", pmf_summary["counts"]["test_winner_counts"].get("moe", 0), "X/25 count", "Free PMF four-family held-out winner count."),
        metric("free_pmf_heldout_winner_count_seam", pmf_summary["counts"]["test_winner_counts"].get("seam", 0), "X/25 count", "Free PMF four-family held-out winner count."),
        metric("fixedk_beats_free_count", regularized_summary["counts"]["fixed_k_beats_free_test"], "X/25 count", "Fixed-k vs free held-out comparison."),
        metric("median_fixedk_minus_free", regularized_summary["medians"]["delta_fixed_k_vs_free_test"], "signed decimal", "Median fixed-k minus free held-out NLL."),
        metric("fixedk_step2_help_count", regularized_summary["counts"]["fixed_k_step2_help_count"], "X/25 count", "Full-corpus residual diagnostic."),
        metric("fixedk_step2_gain_median", regularized_summary["medians"]["fixed_k_step2_gain"], "signed decimal", "Median fixed-k step-2 gain."),
        metric("fixedkw_beats_free_count", regularized_summary["counts"]["fixed_kw_beats_free_test"], "X/25 count", "Fixed-k,w vs free held-out comparison."),
        metric("median_fixedkw_minus_free", regularized_summary["medians"]["delta_fixed_kw_vs_free_test"], "signed decimal", "Median fixed-k,w minus free held-out NLL."),
        metric("fixedkw_step2_help_count", regularized_summary["counts"]["fixed_kw_step2_help_count"], "X/25 count", "Full-corpus residual diagnostic."),
        metric("fixedkw_step2_gain_median", regularized_summary["medians"]["fixed_kw_step2_gain"], "signed decimal", "Median fixed-k,w step-2 gain."),
        metric("softk_beats_free_count", sum(row["softk_test_avg_nll"] < frow["free_test_avg_nll"] for row, frow in zip(table4, build_free_rows(table4, pmf_summary))), "X/25 count", "Canonical Table 4 soft-k vs free Seam-Mandelbrot."),
        metric("median_softk_minus_free", median([row["softk_test_avg_nll"] - frow["free_test_avg_nll"] for row, frow in zip(table4, build_free_rows(table4, pmf_summary))]), "signed decimal", "Median canonical Table 4 soft-k minus free Seam-Mandelbrot."),
        metric("softk_beats_fixedk_count", split_summary["counts"]["softk_beats_fixedk_test"], "X/25 count", "Named splitfit regularization comparison."),
        metric("median_softk_minus_fixedk", split_summary["medians"]["delta_best_vs_fixedk_test"], "signed decimal", "Named splitfit regularization comparison."),
        metric("softk_step2_help_count", split_summary["counts"]["softk_step2_help_count"], "X/25 count", "Full-corpus residual diagnostic from explicit splitfit schema."),
        metric("softk_step2_nonhelp_count", split_summary["counts"]["n_rows"] - split_summary["counts"]["softk_step2_help_count"], "X/25 count", "Complement of soft-k step-2 help count."),
        metric("softk_step2_gain_median", split_summary["medians"]["softk_step2_gain"], "signed decimal", "Median soft-k full-refit residual step-2 gain."),
        metric("free_pmf_step2_help_count", pmf_summary["counts"]["seam_step2_help_count"], "X/25 count", "Free Seam-Mandelbrot full-corpus residual step-2 help count."),
        metric("softk_parameter_count", SOFTK_PARAMETER_COUNT, "integer", "Protocol constant for soft-k Seam-Mandelbrot PMF."),
        metric("lambda_k_min", min(lambdas), "scientific notation", "Minimum selected lambda_k in canonical Table 4."),
        metric("lambda_k_max", max(lambdas), "scientific notation", "Maximum selected lambda_k in canonical Table 4."),
    ]
    for key in ["lambda_count_1e-4", "lambda_count_3e-4", "lambda_count_1e-3", "lambda_count_3e-3", "lambda_count_1e-2", "lambda_count_3e-2"]:
        rows.append(metric(key, lambda_counts.get(key, 0), "count", "Selected lambda_k distribution from canonical Table 4."))
    return rows


def build_free_rows(table4: list[dict[str, Any]], pmf_summary: dict) -> list[dict[str, Any]]:
    by_slug = {row["slug"]: row for row in pmf_summary["rows"]}
    return [{"slug": row["slug"], "free_test_avg_nll": by_slug[row["slug"]]["models"]["seam"]["test_avg_nll"]} for row in table4]


def build_fullrefit_rows(split_json: dict[str, dict], table4: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for row in table4:
        slug = row["slug"]
        full = split_json[slug]["best_lambda"]["full_refit"]
        rows.append(
            {
                "slug": slug,
                "corpus": row["corpus"],
                "softk_fullrefit_train_avg_nll": full["train_avg_nll"],
                "softk_fullrefit_bic": full["bic"],
                "softk_fullrefit_rmse": full["rmse"],
                "softk_fullrefit_step2_gain": full["step2_gain"],
                "softk_fullrefit_step2_helps": full["step2_helps"],
                "softk_fullrefit_k": full["params"]["k"],
                "softk_fullrefit_w": full["params"]["w"],
                "diagnostic_only": True,
            }
        )
    return rows


def main() -> None:
    pmf_rows = {row["slug"]: row for row in read_csv(SOURCE_BUNDLES["pmf_free"] / "seam_mandelbrot_table.csv")}
    reg_rows = {row["slug"]: row for row in read_csv(SOURCE_BUNDLES["regularized"] / "regularized_seam_table.csv")}
    split_rows = {row["slug"]: row for row in read_csv(SOURCE_BUNDLES["softk_splitfit"] / "softk_table.csv")}
    legacy_rows = {row["slug"]: row for row in read_csv(SOURCE_BUNDLES["softk_legacy"] / "softk_table.csv")}
    softkw_rows = {row["slug"]: row for row in read_csv(SOURCE_BUNDLES["softkw"] / "softkw_table.csv")}
    v4_rows = read_csv(CANONICAL_TABLE4_SOURCE)

    pmf_summary = load_json(SOURCE_BUNDLES["pmf_free"] / "summary.json")
    regularized_summary = load_json(SOURCE_BUNDLES["regularized"] / "summary.json")
    split_summary = load_json(SOURCE_BUNDLES["softk_splitfit"] / "summary.json")
    legacy_summary = load_json(SOURCE_BUNDLES["softk_legacy"] / "summary.json")
    softkw_summary = load_json(SOURCE_BUNDLES["softkw"] / "summary.json")

    split_json = {row["slug"]: row for row in split_summary["rows"]}
    legacy_json = {row["slug"]: row for row in legacy_summary["rows"]}
    softkw_json = {row["slug"]: row for row in softkw_summary["rows"]}

    table4 = build_table4(v4_rows)
    provenance = build_provenance_rows(table4, pmf_rows, split_rows, legacy_rows, softkw_rows)
    diagnostics = build_softk_diagnostic(table4, split_rows, legacy_rows, softkw_rows, split_json, legacy_json, softkw_json)
    variants = build_variant_rows(table4, pmf_rows, reg_rows, split_rows, softkw_rows)
    aggregates = build_aggregate_rows(table4, pmf_summary, regularized_summary, split_summary)
    fullrefit = build_fullrefit_rows(split_json, table4)

    write_csv(OUTPUTS["table4"], table4)
    write_csv(OUTPUTS["table4_provenance"], provenance)
    write_csv(OUTPUTS["softk_diagnostic"], diagnostics)
    write_csv(OUTPUTS["variant_per_corpus"], variants)
    write_csv(OUTPUTS["aggregate"], aggregates, ["metric_name", "value", "display_format", "notes"])
    write_csv(OUTPUTS["fullrefit_per_corpus"], fullrefit)
    fullrefit_aggregates = [
        metric("softk_fullrefit_step2_help_count", sum(bool(row["softk_fullrefit_step2_helps"]) for row in fullrefit), "X/25 count", "Diagnostic only; not held-out canonical."),
        metric("softk_fullrefit_step2_gain_median", median(row["softk_fullrefit_step2_gain"] for row in fullrefit), "signed decimal", "Diagnostic only; not held-out canonical."),
    ]
    write_csv(OUTPUTS["fullrefit_aggregate"], fullrefit_aggregates, ["metric_name", "value", "display_format", "notes"])

    archive = EXPERIMENT_DIR / "archive" / "provenance"
    archive.mkdir(parents=True, exist_ok=True)
    for bundle in SOURCE_BUNDLES.values():
        for src in bundle.glob("*"):
            if src.suffix in {".json", ".csv", ".md"}:
                shutil.copy2(src, archive / f"{bundle.name}_{src.name}")

    known_drift_count = sum(row["is_known_local_optimum_discrepancy"] for row in diagnostics)
    manifest = {
        "experiment_id": EXPERIMENT_ID,
        "status": "complete_canonical_by_documentation",
        "canonical_table4_source": str(CANONICAL_TABLE4_SOURCE.relative_to(EXPERIMENT_DIR.parents[1])),
        "canonical_reproducibility_status": "canonical-by-documentation, not canonical-by-reproducible-producer",
        "known_softk_local_optimum_discrepancy_rows": known_drift_count,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "outputs/splitfit/table4_fourway.csv": {"rows": len(table4)},
            "outputs/splitfit/table4_provenance.csv": {"rows": len(provenance)},
            "outputs/splitfit/softk_source_diagnostic.csv": {"rows": len(diagnostics)},
            "outputs/splitfit/pmf_variant_per_corpus.csv": {"rows": len(variants)},
            "outputs/splitfit/aggregate_statistics.csv": {"rows": len(aggregates)},
            "outputs/fullrefit/fourway_per_corpus.csv": {"rows": len(fullrefit)},
            "outputs/fullrefit/aggregate_statistics.csv": {"rows": len(fullrefit_aggregates)},
        },
    }
    OUTPUTS["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote 7a outputs: {len(table4)} table rows, {len(aggregates)} aggregate rows, {known_drift_count} soft-k source discrepancies")


if __name__ == "__main__":
    main()

