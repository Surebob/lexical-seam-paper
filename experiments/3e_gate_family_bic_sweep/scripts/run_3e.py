#!/usr/bin/env python3
"""Migrate S2 v3 decoupled gate-family sweep outputs into canonical experiment 3e."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
EXP_DIR = SCRIPT_DIR.parents[0]
ROOT = SCRIPT_DIR.parents[2]
sys.path.insert(0, str(EXP_DIR))

from source_config import INDEPENDENT_GATES, ORDERED_GATES, OUTPUTS, PROVENANCE_SCRIPTS, SOURCES  # noqa: E402


PARAM_FIELDS = ["a1", "b1", "c1", "a2", "b2", "c2"]


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def metric_map(rows: list[dict[str, str]]) -> dict[str, str]:
    return {row["metric_name"]: row["value"] for row in rows}


def boolish(value: str) -> bool:
    return str(value).strip().lower() == "true"


def copy_provenance() -> list[str]:
    copied = []
    archive = EXP_DIR / "archive" / "provenance"
    archive.mkdir(parents=True, exist_ok=True)
    for name, src in PROVENANCE_SCRIPTS.items():
        if src.exists():
            dst = archive / src.name
            shutil.copy2(src, dst)
            copied.append(str(dst.relative_to(ROOT)))
    return copied


def build_per_fit(output_dir: Path) -> tuple[list[dict[str, str]], list[str], dict[str, int]]:
    fit_rows = read_rows(SOURCES["real_gate_sweep_dir"] / "s2_v3_per_fit_results.csv")
    dispersion = {
        (row["slug"], row["gate"]): row
        for row in read_rows(SOURCES["real_gate_sweep_dir"] / "s2_v3_per_start_dispersion.csv")
    }
    logistic_params = {}
    if SOURCES["logistic_parameter_source"].exists():
        logistic_params = {
            row["slug"]: row for row in read_rows(SOURCES["logistic_parameter_source"])
        }

    canonical = []
    missing_param_cells = 0
    for row in fit_rows:
        gate = row["gate"]
        slug = row["slug"]
        drow = dispersion.get((slug, gate), {})
        out = {
            "slug": slug,
            "corpus": row["corpus"],
            "gate_family": gate,
            "seed": row.get("seed", ""),
            "a1": "",
            "b1": "",
            "c1": "",
            "a2": "",
            "b2": "",
            "c2": "",
            "k": row.get("k", ""),
            "w_gate": row.get("w_gate", ""),
            "w_tail": row.get("w_tail", ""),
            "BIC": row.get("bic", ""),
            "final_rmse": row.get("rmse", ""),
            "k_hit_lower_bound": row.get("k_hit_lower_bound", ""),
            "k_hit_upper_bound": row.get("k_hit_upper_bound", ""),
            "w_gate_hit_lower_bound": row.get("w_gate_hit_lower_bound", ""),
            "w_gate_hit_upper_bound": row.get("w_gate_hit_upper_bound", ""),
            "w_tail_hit_lower_bound": row.get("w_tail_hit_lower_bound", ""),
            "w_tail_hit_upper_bound": row.get("w_tail_hit_upper_bound", ""),
            "any_bounds_hit": str(
                any(
                    boolish(row.get(flag, "False"))
                    for flag in [
                        "k_hit_lower_bound",
                        "k_hit_upper_bound",
                        "w_gate_hit_lower_bound",
                        "w_gate_hit_upper_bound",
                        "w_tail_hit_lower_bound",
                        "w_tail_hit_upper_bound",
                    ]
                )
            ),
            "best_start_index": row.get("best_start_index", ""),
            "best_nfev": row.get("best_nfev", ""),
            "wall_clock_seconds": row.get("wall_clock_seconds", row.get("runtime_sec", "")),
            "peak_memory_mb": row.get("peak_memory_mb", ""),
            "starts_within_1_bic": row.get("starts_within_1_bic", ""),
            "optimizer_warm_start_info": "100 random starts; best objective selected; per-start dispersion summarized.",
            "dispersion_min_bic": drow.get("min_bic", ""),
            "dispersion_median_bic": drow.get("median_bic", ""),
            "dispersion_max_bic": drow.get("max_bic", ""),
            "dispersion_across_starts": (
                f"min={drow.get('min_bic','')};median={drow.get('median_bic','')};"
                f"max={drow.get('max_bic','')};within1={row.get('starts_within_1_bic','')}"
            ),
            "parameter_source": "s2_v3_per_fit_results.csv",
        }
        if gate == "logistic" and slug in logistic_params:
            lp = logistic_params[slug]
            for field in PARAM_FIELDS:
                out[field] = lp.get(f"logistic_{field}", "")
            out["parameter_source"] = "s2_v3_logistic_local_params.csv for a/b/c; s2_v3_per_fit_results.csv for k/w/BIC/RMSE"
        else:
            missing_param_cells += len(PARAM_FIELDS)
            out["parameter_source"] = "non-logistic a/b/c coefficients not saved in source bundle"
        canonical.append(out)

    fieldnames = [
        "slug",
        "corpus",
        "gate_family",
        "seed",
        *PARAM_FIELDS,
        "k",
        "w_gate",
        "w_tail",
        "BIC",
        "final_rmse",
        "k_hit_lower_bound",
        "k_hit_upper_bound",
        "w_gate_hit_lower_bound",
        "w_gate_hit_upper_bound",
        "w_tail_hit_lower_bound",
        "w_tail_hit_upper_bound",
        "any_bounds_hit",
        "best_start_index",
        "best_nfev",
        "wall_clock_seconds",
        "peak_memory_mb",
        "starts_within_1_bic",
        "optimizer_warm_start_info",
        "dispersion_min_bic",
        "dispersion_median_bic",
        "dispersion_max_bic",
        "dispersion_across_starts",
        "parameter_source",
    ]
    write_rows(output_dir / "per_fit_results.csv", canonical, fieldnames)
    return canonical, fieldnames, {"missing_nonlogistic_regime_parameter_cells": missing_param_cells}


def build_per_corpus(output_dir: Path) -> tuple[list[dict[str, str]], list[str]]:
    source = read_rows(SOURCES["real_gate_sweep_dir"] / "s2_v3_per_corpus_results.csv")
    canonical = []
    for row in source:
        winner = row["winner_gate"]
        out = {
            "slug": row["slug"],
            "corpus": row["corpus"],
            "seed": row.get("seed", ""),
            "vocab_size": row.get("vocab_size", ""),
            "winning_gate": winner,
            "winner_bic": row.get(f"{winner}_bic", ""),
            "independent_gate_spread": row.get("bic_spread", ""),
            "tanh_calibration_pass": row.get("tanh_calibration_pass", ""),
            "bic_tanh_minus_logistic": row.get("bic_tanh_minus_logistic", ""),
            "w_gate_tanh_over_2w_gate_logistic_ratio": row.get(
                "w_gate_tanh_over_2w_gate_logistic_ratio", ""
            ),
            "w_tail_tanh_over_w_tail_logistic_ratio": row.get(
                "w_tail_tanh_over_w_tail_logistic_ratio", ""
            ),
            "logistic_bic": row.get("logistic_bic", ""),
            "erf_bic": row.get("erf_bic", ""),
            "algebraic_bic": row.get("algebraic_bic", ""),
            "arctan_bic": row.get("arctan_bic", ""),
            "worst_gate": row.get("worst_gate", ""),
        }
        canonical.append(out)
    fieldnames = list(canonical[0].keys()) if canonical else []
    write_rows(output_dir / "per_corpus_results.csv", canonical, fieldnames)
    return canonical, fieldnames


def build_aggregate(output_dir: Path) -> tuple[list[dict[str, str]], list[str]]:
    src_rows = read_rows(SOURCES["real_gate_sweep_dir"] / "s2_v3_aggregate_statistics.csv")
    src = metric_map(src_rows)
    source_notes = {row["metric_name"]: row.get("notes", "") for row in src_rows}
    mapping = [
        ("erf_wins_count", "erf_bic_wins", "integer", "Erf wins among independent gates."),
        ("arctan_wins_count", "arctan_bic_wins", "integer", "Arctan wins among independent gates."),
        ("logistic_wins_count", "logistic_bic_wins", "integer", "Logistic wins among independent gates."),
        ("algebraic_wins_count", "algebraic_bic_wins", "integer", "Algebraic wins among independent gates."),
        ("tanh_calibration_pass_count", "tanh_calibration_pass_count", "integer", "Tanh/logistic calibration pass count."),
        (
            "median_independent_gate_bic_spread",
            "median_bic_spread",
            "decimal_2",
            "Median BIC spread across logistic, erf, algebraic, arctan.",
        ),
        (
            "mean_independent_gate_bic_spread",
            "mean_bic_spread",
            "decimal_2",
            "Mean BIC spread across logistic, erf, algebraic, arctan.",
        ),
        (
            "gates_indistinguishable_strict_count",
            "gates_indistinguishable_strict_count",
            "integer",
            "Count with independent-gate spread < 2.",
        ),
        (
            "gates_indistinguishable_positive_count",
            "gates_indistinguishable_positive_count",
            "integer",
            "Count with independent-gate spread < 6.",
        ),
        (
            "gates_indistinguishable_strong_count",
            "gates_indistinguishable_strong_count",
            "integer",
            "Count with independent-gate spread < 10.",
        ),
    ]
    rows = []
    for canonical_name, source_name, fmt, note in mapping:
        rows.append(
            {
                "metric_name": canonical_name,
                "value": src[source_name],
                "display_format": fmt,
                "source_metric_name": source_name,
                "notes": f"{note} Source note: {source_notes.get(source_name, '')}",
            }
        )
    fieldnames = ["metric_name", "value", "display_format", "source_metric_name", "notes"]
    write_rows(output_dir / "aggregate_statistics.csv", rows, fieldnames)
    return rows, fieldnames


def build_synthetic(output_dir: Path) -> tuple[list[dict[str, str]], list[str], list[dict[str, str]], list[str]]:
    source = read_rows(SOURCES["synthetic_recovery_dir"] / "synthetic_gate_recovery_per_corpus.csv")
    rows = []
    for row in source:
        match = row.get("winner_gate", "") == row.get("generator_gate", "logistic")
        erf_beats_true = float(row["erf_bic"]) < float(row["logistic_bic"])
        rows.append(
            {
                "synthetic_corpus_id": row["slug"],
                "corpus": row["corpus"],
                "true_gate": row.get("generator_gate", "logistic"),
                "recovered_gate_by_fitter": row.get("winner_gate", ""),
                "recovery_match": str(match),
                "erf_beats_true": str(erf_beats_true),
                "winner_bic": row.get(f"{row.get('winner_gate')}_bic", ""),
                "logistic_bic": row.get("logistic_bic", ""),
                "erf_bic": row.get("erf_bic", ""),
                "bic_erf_minus_logistic": row.get("bic_erf_minus_logistic", ""),
                "independent_gate_spread": row.get("bic_spread", ""),
            }
        )
    fieldnames = list(rows[0].keys()) if rows else []
    write_rows(output_dir / "synthetic_recovery.csv", rows, fieldnames)

    src_agg_rows = read_rows(SOURCES["synthetic_recovery_dir"] / "synthetic_gate_recovery_aggregate_statistics.csv")
    src = metric_map(src_agg_rows)
    agg = [
        {
            "metric_name": "logistic_recovery_count",
            "value": src["logistic_recovered_count"],
            "display_format": "X/Y",
            "source_metric_name": "logistic_recovered_count",
            "notes": "Count where logistic-generated synthetic data recovered logistic.",
        },
        {
            "metric_name": "erf_beats_true_count",
            "value": src["erf_beats_logistic_count"],
            "display_format": "X/Y",
            "source_metric_name": "erf_beats_logistic_count",
            "notes": "Count where erf beats the true logistic generator on logistic-generated data.",
        },
        {
            "metric_name": "synthetic_completed_corpus_count",
            "value": src["completed_corpus_count"],
            "display_format": "integer",
            "source_metric_name": "completed_corpus_count",
            "notes": "Number of synthetic corpora with all gate fits complete.",
        },
        {
            "metric_name": "synthetic_completed_fit_count",
            "value": src["completed_fit_count"],
            "display_format": "integer",
            "source_metric_name": "completed_fit_count",
            "notes": "Number of synthetic corpus/gate fits complete.",
        },
    ]
    agg_fields = ["metric_name", "value", "display_format", "source_metric_name", "notes"]
    write_rows(output_dir / "synthetic_recovery_aggregate.csv", agg, agg_fields)
    return rows, fieldnames, agg, agg_fields


def write_manifest(output_dir: Path, produced: dict[str, dict], provenance: list[str], limits: dict[str, int]) -> None:
    manifest = {
        "experiment_id": "3e_gate_family_bic_sweep",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve().relative_to(ROOT)),
        "source_bundles": {name: str(path.relative_to(ROOT)) for name, path in SOURCES.items()},
        "provenance_scripts_copied": provenance,
        "outputs": produced,
        "audit_limitations": [
            "Non-logistic a1,b1,c1,a2,b2,c2 coefficients were not saved in the source full-sweep CSVs.",
            "BIC/RMSE/winner claims are fully reconstructible from saved source bundles.",
            "Synthetic-recovery tanh calibration is not interpreted because near-zero RMSE synthetic fits magnify tiny tanh/logistic numerical differences in BIC.",
        ],
        "limits": limits,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS["aggregate"].parent)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    provenance = copy_provenance()
    per_fit, per_fit_fields, limits = build_per_fit(output_dir)
    per_corpus, per_corpus_fields = build_per_corpus(output_dir)
    aggregate, aggregate_fields = build_aggregate(output_dir)
    synthetic, synthetic_fields, synthetic_agg, synthetic_agg_fields = build_synthetic(output_dir)

    produced = {
        "per_fit_results.csv": {
            "rows": len(per_fit),
            "schema": per_fit_fields,
        },
        "per_corpus_results.csv": {
            "rows": len(per_corpus),
            "schema": per_corpus_fields,
        },
        "aggregate_statistics.csv": {
            "rows": len(aggregate),
            "schema": aggregate_fields,
        },
        "synthetic_recovery.csv": {
            "rows": len(synthetic),
            "schema": synthetic_fields,
        },
        "synthetic_recovery_aggregate.csv": {
            "rows": len(synthetic_agg),
            "schema": synthetic_agg_fields,
        },
    }
    write_manifest(output_dir, produced, provenance, limits)


if __name__ == "__main__":
    main()
