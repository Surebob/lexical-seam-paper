from __future__ import annotations

import csv
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

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


def load_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def f(value: object) -> float:
    return float(value)


def boolish(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def metric(name: str, value: object, display_format: str, notes: str) -> dict[str, object]:
    return {
        "metric_name": name,
        "value": value,
        "display_format": display_format,
        "notes": notes,
    }


def build_hybrid_headtail_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in summary["rows"]:
        full = item["heldout"]["full"]
        top100 = item["heldout"]["top100"]["avg_nll"]
        avg = full["avg_nll"]
        step2 = item["step2"]["full"]
        rows.append(
            {
                "slug": item["slug"],
                "corpus": item["name"],
                "token_count": item["token_count"],
                "vocabulary_size": item["vocab_size"],
                "selected_lambda": item["best_lambda"]["lambda"],
                "k": item["best_lambda"]["selection"]["params"]["k"],
                "w": item["best_lambda"]["selection"]["params"]["w"],
                "top100_hybrid": top100["hybrid"],
                "top100_moe": top100["moe"],
                "top100_softk": top100["softk"],
                "full_zipf": avg["zipf"],
                "full_zm": avg["zm"],
                "full_moe": avg["moe"],
                "full_softk": avg["softk"],
                "full_hybrid": avg["hybrid"],
                "hybrid_minus_moe": full["hybrid_minus_moe"],
                "hybrid_minus_zm": full["hybrid_minus_zm"],
                "hybrid_minus_softk": full["hybrid_minus_softk"],
                "hybrid_minus_zipf": full["hybrid_minus_zipf"],
                "canonical_hybrid_vs_softk_winner": "hybrid" if full["hybrid_minus_softk"] < 0 else "softk",
                "full_step2_help": step2["helps"],
                "full_step2_gain": step2["gain"],
                "full_step2_expression": step2["expr"],
                "canonical_source": "results/zipf_hybrid_headtail_splitfit/summary.json",
            }
        )
    return rows


def build_fully_ours_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in summary["rows"]:
        for variant in ["nested", "three_regime"]:
            source = item["variants"][variant]
            rows.append(
                {
                    "slug": item["slug"],
                    "corpus": item["name"],
                    "token_count": item["token_count"],
                    "vocabulary_size": item["vocab_size"],
                    "variant": variant,
                    "test_avg_nll": source["test_avg_nll"],
                    "bic": source["bic"],
                    "rmse": source["rmse"],
                    "step2_gain": source["step2_gain"],
                    "step2_help": source["step2_helps"],
                    "params_json": json.dumps(source["params"], sort_keys=True),
                    "canonical_source": "results/zipf_fully_ours_variants/summary.json",
                }
            )
    return rows


def build_mechpenalty_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in summary["corpora"]:
        for source in item["rows"]:
            fit = source["fit"]
            params = fit["params"]
            rows.append(
                {
                    "slug": item["slug"],
                    "corpus": item["name"],
                    "token_count": item["token_count"],
                    "vocabulary_size": item["vocab_size"],
                    "lambda_mech": source["lambda_mech"],
                    "test_avg_nll": source["test_avg_nll"],
                    "test_loglike": source["test_loglike"],
                    "train_step2_gain": fit["train_step2_gain"],
                    "train_step2_expression": fit["train_step2_expr"],
                    "heldout_step2_gain": source["heldout_step2_gain"],
                    "heldout_step2_help": source["heldout_step2_helps"],
                    "heldout_step2_expression": source["heldout_step2_expr"],
                    "k": params["k"],
                    "w": params["w"],
                    "alpha_tail": params["alpha_tail"],
                    "beta_tail": params["beta_tail"],
                    "transition_fraction": params["transition_fraction"],
                    "success": fit["success"],
                    "optimizer_message": fit["message"],
                    "canonical_source": "results/zipf_hybrid_mechpenalty/summary.json",
                }
            )
    return rows


def build_hybrid_vs_softk_diagnostic(
    analysis_rows: list[dict[str, str]],
    hybrid_rows: list[dict[str, object]],
    softk_source_diag: list[dict[str, str]],
) -> list[dict[str, object]]:
    canonical = {row["slug"]: row for row in hybrid_rows}
    source_diag = {row["slug"]: row for row in softk_source_diag}
    out: list[dict[str, object]] = []
    for row in analysis_rows:
        slug = row["slug"]
        c = canonical[slug]
        d = source_diag.get(slug, {})
        canonical_delta = f(c["hybrid_minus_softk"])
        analysis_delta = f(row["delta_hybrid_minus_softk"])
        out.append(
            {
                "slug": slug,
                "corpus": row["name"],
                "canonical_hybrid_full_nll": c["full_hybrid"],
                "canonical_softk_full_nll": c["full_softk"],
                "canonical_delta_hybrid_minus_softk": canonical_delta,
                "canonical_winner": "hybrid" if canonical_delta < 0 else "softk",
                "analysis_hybrid_full_nll": row["hybrid_full_nll"],
                "analysis_softk_full_nll": row["softk_full_nll"],
                "analysis_delta_hybrid_minus_softk": analysis_delta,
                "analysis_winner": row["winner"],
                "softk_comparator_delta_analysis_minus_canonical": f(row["softk_full_nll"]) - f(c["full_softk"]),
                "analysis_softk_matched_source_group": d.get("matched_source_group", "legacy_derivative_branch"),
                "analysis_source_verdict": d.get("verdict", "known_legacy_pre_splitfit_comparator"),
                "canonical_softk_generation": "splitfit_current",
                "analysis_softk_generation": "legacy_pre_splitfit_derivative_branch",
                "notes": "Different comparator generations explain the 22/25 canonical versus 16/25 legacy hybrid-win counts.",
            }
        )
    return out


def build_structure_summary(
    hybrid_summary: dict[str, object],
    fully_summary: dict[str, object],
    mech_summary: dict[str, object],
    analysis_summary: dict[str, object],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(component: str, metric_name: str, value: object, display: str, source: str, notes: str) -> None:
        rows.append(
            {
                "component": component,
                "metric_name": metric_name,
                "value": value,
                "display_format": display,
                "source": source,
                "notes": notes,
            }
        )

    full = hybrid_summary["cutoffs"]["full"]
    for key in [
        "hybrid_beats_moe",
        "hybrid_beats_zm",
        "hybrid_beats_softk",
        "hybrid_beats_zipf",
        "median_hybrid_minus_moe",
        "median_hybrid_minus_zm",
        "median_hybrid_minus_softk",
        "step2_help_count",
    ]:
        add("hybrid_headtail", key, full[key], "count or decimal", "zipf_hybrid_headtail_splitfit", "Canonical splitfit hybrid source.")

    for variant, stats in fully_summary["variants"].items():
        for key, value in stats.items():
            add(variant, key, json.dumps(value, sort_keys=True) if isinstance(value, dict) else value, "count, decimal, or JSON", "zipf_fully_ours_variants", "Fully-ours variant summary.")

    add(
        "hybrid_vs_softk_legacy_diagnostic",
        "legacy_hybrid_win_count",
        analysis_summary["by_winner"]["hybrid"]["count"],
        "X/25 count",
        "zipf_hybrid_vs_softk_analysis",
        "Legacy pre-splitfit comparator branch; preserved as diagnostic, not canonical headline.",
    )
    add(
        "hybrid_vs_softk_legacy_diagnostic",
        "legacy_softk_win_count",
        analysis_summary["by_winner"]["softk"]["count"],
        "X/25 count",
        "zipf_hybrid_vs_softk_analysis",
        "Legacy pre-splitfit comparator branch; preserved as diagnostic, not canonical headline.",
    )

    for lam, stats in mech_summary["lambdas"].items():
        for key, value in stats.items():
            add(
                "mechanism_penalty",
                f"lambda_{lam}_{key}",
                json.dumps(value, sort_keys=True) if isinstance(value, dict) else value,
                "count, decimal, or JSON",
                "zipf_hybrid_mechpenalty",
                "Mechanism-penalty sweep row.",
            )
    add(
        "mechanism_penalty",
        "frontier_decision",
        "no_lambda_satisfies_intended_rule",
        "categorical verdict",
        "zipf_hybrid_mechpenalty",
        "No penalty simultaneously keeps held-out step-2 help <=4/25 and beats MOE on >=18/25.",
    )

    for key, value in cfg.NESTED_PROTOCOL_CONSTANTS.items():
        add("nested_protocol_constants", key, value, "protocol constant", "claim_map_protocol_row", "Protocol constant preserved for Paper 2 provenance.")

    return rows


def build_aggregates(
    hybrid_summary: dict[str, object],
    fully_summary: dict[str, object],
    mech_summary: dict[str, object],
) -> list[dict[str, object]]:
    full = hybrid_summary["cutoffs"]["full"]
    nested = fully_summary["variants"]["nested"]
    three = fully_summary["variants"]["three_regime"]
    lambda_grid = [float(key) for key in mech_summary["lambdas"].keys()]
    lambda_grid = sorted(lambda_grid)
    all_help_counts = {
        str(key): value["heldout_step2_help_count"]
        for key, value in sorted(mech_summary["lambdas"].items(), key=lambda item: float(item[0]))
    }

    aggregates = [
        metric("hybrid_beats_moe_count", full["hybrid_beats_moe"], "X/25 count", "Canonical splitfit hybrid head-tail beats MOE on full held-out NLL."),
        metric("hybrid_beats_softk_count", full["hybrid_beats_softk"], "X/25 count", "Canonical splitfit hybrid head-tail beats current splitfit soft-k comparator."),
        metric("median_hybrid_minus_softk", full["median_hybrid_minus_softk"], "decimal NLL difference", "Median full held-out hybrid minus soft-k; negative favors hybrid."),
        metric("hybrid_step2_help_count", full["step2_help_count"], "X/25 count", "Count where residual step-2 search helps the hybrid model."),
        metric("nested_beats_moe_count", nested["beats_moe"], "X/25 count", "Nested seam architecture beats MOE on held-out NLL."),
        metric("nested_beats_softk_count", nested["beats_softk"], "X/25 count", "Nested seam architecture beats soft-k on held-out NLL."),
        metric("nested_beats_zm_count", nested["beats_zm"], "X/25 count", "Nested seam architecture beats ZM on held-out NLL."),
        metric("nested_fourway_winner_count", nested["winner_counts"]["nested"], "X/25 count", "Nested seam is four-way held-out winner among Zipf, ZM, MOE, and nested."),
        metric("nested_step2_help_count", nested["step2_help_count"], "X/25 count", "Count where residual step-2 search helps nested seam."),
        metric("three_regime_beats_softk_count", three["beats_softk"], "X/25 count", "Three-regime architecture beats soft-k on held-out NLL."),
        metric("three_regime_median_minus_softk", three["median_minus_softk"], "decimal NLL difference", "Median three-regime minus soft-k; negative favors three-regime."),
        metric("mechanism_penalty_grid", json.dumps(lambda_grid), "JSON numeric list", "Penalty strengths swept in the mechanism-penalty hybrid."),
        metric(
            "mechanism_penalty_step2_help_count_all_lambdas",
            json.dumps(all_help_counts, sort_keys=True),
            "JSON count map",
            "Held-out step-2 help count is 25/25 for every penalty strength.",
        ),
        metric(
            "mechanism_penalty_frontier_result",
            "no lambda satisfies help <=4/25 and beats MOE >=18/25",
            "verdict string",
            "Primary negative finding: penalty does not locate the intended likelihood-mechanism Pareto point.",
        ),
    ]
    for key, value in cfg.NESTED_PROTOCOL_CONSTANTS.items():
        aggregates.append(metric(key, value, "protocol constant", "Nested seam protocol constant from claim-map/manuscript method rows."))
    return aggregates


def archive_sources() -> None:
    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    source_bundles = {name: str(path.relative_to(cfg.ROOT)) for name, path in cfg.SOURCE_BUNDLES.items()}
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps({"source_bundles": source_bundles}, indent=2) + "\n",
        encoding="utf-8",
    )
    provenance = cfg.ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for name, bundle in cfg.SOURCE_BUNDLES.items():
        dest = provenance / name
        dest.mkdir(parents=True, exist_ok=True)
        for src in sorted(bundle.iterdir()):
            if src.is_file():
                shutil.copy2(src, dest / src.name)
    for src in [
        cfg.EXPERIMENT_DIR / "CONSOLIDATION_PLAN.md",
        cfg.EXPERIMENT_DIR / "SOFTK_COMPARATOR_DIAGNOSTIC.md",
        cfg.SOURCE_FILES["softk_comparator_source_diagnostic"],
    ]:
        if src.exists():
            shutil.copy2(src, provenance / src.name)


def write_manifest(outputs: dict[str, list[dict[str, object]]], schemas: dict[str, list[str]]) -> None:
    manifest = {
        "experiment_id": "9_likelihood_frontier",
        "status": "complete_from_historical_bundles",
        "migration_type": "multi_bundle_consolidation",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_bundles": {name: str(path.relative_to(cfg.ROOT)) for name, path in cfg.SOURCE_BUNDLES.items()},
        "canonical_decisions": {
            "hybrid_vs_softk": "zipf_hybrid_headtail_splitfit is canonical for manuscript-facing hybrid-vs-soft-k counts.",
            "legacy_diagnostic": "zipf_hybrid_vs_softk_analysis is preserved as a legitimate legacy pre-splitfit comparator sidecar.",
            "paper_scope": "PMF arc including experiment 9 is queued for Paper 2; Paper 1 will remove Pareto-frontier discussion.",
        },
        "outputs": {
            filename: {"rows": len(rows), "schema": schemas[filename]}
            for filename, rows in outputs.items()
        },
        "claim_map_rows_satisfied": [
            "hybrid-vs-MOE and hybrid-vs-soft-k held-out counts",
            "median hybrid-minus-soft-k held-out delta",
            "hybrid residual step-2 help count",
            "nested seam likelihood-frontier counts and protocol constants",
            "mechanism-penalty negative-frontier result",
        ],
    }
    cfg.OUTPUTS["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    cfg.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    hybrid_summary = load_json(cfg.SOURCE_FILES["hybrid_headtail_summary"])
    mech_summary = load_json(cfg.SOURCE_FILES["mechpenalty_summary"])
    analysis_summary = load_json(cfg.SOURCE_FILES["analysis_summary"])
    fully_summary = load_json(cfg.SOURCE_FILES["fully_ours_summary"])

    hybrid_rows = build_hybrid_headtail_rows(hybrid_summary)
    fully_rows = build_fully_ours_rows(fully_summary)
    mech_rows = build_mechpenalty_rows(mech_summary)
    analysis_rows = read_csv(cfg.SOURCE_FILES["analysis_table"])
    softk_source_diag = read_csv(cfg.SOURCE_FILES["softk_comparator_source_diagnostic"])
    diagnostic_rows = build_hybrid_vs_softk_diagnostic(analysis_rows, hybrid_rows, softk_source_diag)
    structure_rows = build_structure_summary(hybrid_summary, fully_summary, mech_summary, analysis_summary)
    aggregate_rows = build_aggregates(hybrid_summary, fully_summary, mech_summary)

    schemas = {
        "hybrid_headtail_per_corpus.csv": [
            "slug",
            "corpus",
            "token_count",
            "vocabulary_size",
            "selected_lambda",
            "k",
            "w",
            "top100_hybrid",
            "top100_moe",
            "top100_softk",
            "full_zipf",
            "full_zm",
            "full_moe",
            "full_softk",
            "full_hybrid",
            "hybrid_minus_moe",
            "hybrid_minus_zm",
            "hybrid_minus_softk",
            "hybrid_minus_zipf",
            "canonical_hybrid_vs_softk_winner",
            "full_step2_help",
            "full_step2_gain",
            "full_step2_expression",
            "canonical_source",
        ],
        "fully_ours_variants.csv": [
            "slug",
            "corpus",
            "token_count",
            "vocabulary_size",
            "variant",
            "test_avg_nll",
            "bic",
            "rmse",
            "step2_gain",
            "step2_help",
            "params_json",
            "canonical_source",
        ],
        "mechanism_penalty_sweep.csv": [
            "slug",
            "corpus",
            "token_count",
            "vocabulary_size",
            "lambda_mech",
            "test_avg_nll",
            "test_loglike",
            "train_step2_gain",
            "train_step2_expression",
            "heldout_step2_gain",
            "heldout_step2_help",
            "heldout_step2_expression",
            "k",
            "w",
            "alpha_tail",
            "beta_tail",
            "transition_fraction",
            "success",
            "optimizer_message",
            "canonical_source",
        ],
        "hybrid_vs_softk_diagnostic.csv": [
            "slug",
            "corpus",
            "canonical_hybrid_full_nll",
            "canonical_softk_full_nll",
            "canonical_delta_hybrid_minus_softk",
            "canonical_winner",
            "analysis_hybrid_full_nll",
            "analysis_softk_full_nll",
            "analysis_delta_hybrid_minus_softk",
            "analysis_winner",
            "softk_comparator_delta_analysis_minus_canonical",
            "analysis_softk_matched_source_group",
            "analysis_source_verdict",
            "canonical_softk_generation",
            "analysis_softk_generation",
            "notes",
        ],
        "hybrid_structure_summary.csv": [
            "component",
            "metric_name",
            "value",
            "display_format",
            "source",
            "notes",
        ],
        "aggregate_statistics.csv": ["metric_name", "value", "display_format", "notes"],
    }

    outputs = {
        "hybrid_headtail_per_corpus.csv": hybrid_rows,
        "fully_ours_variants.csv": fully_rows,
        "mechanism_penalty_sweep.csv": mech_rows,
        "hybrid_vs_softk_diagnostic.csv": diagnostic_rows,
        "hybrid_structure_summary.csv": structure_rows,
        "aggregate_statistics.csv": aggregate_rows,
    }

    write_csv(cfg.OUTPUTS["hybrid_headtail"], hybrid_rows, schemas["hybrid_headtail_per_corpus.csv"])
    write_csv(cfg.OUTPUTS["fully_ours"], fully_rows, schemas["fully_ours_variants.csv"])
    write_csv(cfg.OUTPUTS["mechpenalty"], mech_rows, schemas["mechanism_penalty_sweep.csv"])
    write_csv(cfg.OUTPUTS["hybrid_vs_softk_diagnostic"], diagnostic_rows, schemas["hybrid_vs_softk_diagnostic.csv"])
    write_csv(cfg.OUTPUTS["hybrid_structure_summary"], structure_rows, schemas["hybrid_structure_summary.csv"])
    write_csv(cfg.OUTPUTS["aggregate"], aggregate_rows, schemas["aggregate_statistics.csv"])
    archive_sources()
    write_manifest(outputs, schemas)

    print(
        "Wrote 9 outputs: "
        f"{len(hybrid_rows)} hybrid rows, "
        f"{len(fully_rows)} fully-ours rows, "
        f"{len(mech_rows)} penalty rows, "
        f"{len(diagnostic_rows)} comparator diagnostic rows, "
        f"{len(aggregate_rows)} aggregate rows"
    )


if __name__ == "__main__":
    main()
