from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
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


def build_hierk_summary(summary: dict[str, object]) -> list[dict[str, object]]:
    return [
        {
            "model": "hierarchical_k",
            "alpha": summary["alpha"],
            "sigma": summary["sigma"],
            "sigma_lower_bound": 0.05,
            "sigma_at_lower_bound": abs(summary["sigma"] - 0.05) < 1e-12,
            "iteration_count": len(summary["history"]),
            "source": "results/zipf_angle4_hierk/summary.json",
        }
    ]


def build_history(summary: dict[str, object]) -> list[dict[str, object]]:
    return [
        {
            "iteration": item["iteration"],
            "alpha": item["alpha"],
            "sigma": item["sigma"],
            "sigma_at_lower_bound": abs(item["sigma"] - 0.05) < 1e-12,
        }
        for item in summary["history"]
    ]


def build_head_window_rows(table_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in table_rows:
        rows.append(
            {
                "slug": row["slug"],
                "corpus": row["name"],
                "cutoff": row["cutoff"],
                "winner_family": row["winner"],
                "zipf_avg_nll": f(row["zipf"]),
                "zm_avg_nll": f(row["zm"]),
                "moe_avg_nll": f(row["moe"]),
                "hierk_avg_nll": f(row["hier"]),
                "hierk_minus_moe": f(row["hier_minus_moe"]),
                "hierk_minus_zm": f(row["hier_minus_zm"]),
                "step2_helps": bool(int(row["step2_help"])),
                "step2_gain": f(row["step2_gain"]),
                "canonical_source": "results/zipf_angle4_hierk/hierk_head_window_table.csv",
            }
        )
    return rows


def build_softk_comparison(head_rows: list[dict[str, object]], softk_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    softk = {row["slug"]: row for row in softk_rows}
    rows: list[dict[str, object]] = []
    for row in [item for item in head_rows if item["cutoff"] == "full"]:
        source = softk[row["slug"]]
        delta = f(row["hierk_avg_nll"]) - f(source["softk_test_avg_nll"])
        rows.append(
            {
                "slug": row["slug"],
                "corpus": row["corpus"],
                "hierk_full_avg_nll": row["hierk_avg_nll"],
                "softk_full_avg_nll": f(source["softk_test_avg_nll"]),
                "hierk_minus_softk": delta,
                "hierk_beats_softk": delta < 0,
                "softk_source": "experiments/7a_canonical_pmf_family/outputs/splitfit/table4_fourway.csv",
            }
        )
    return rows


def build_aggregates(summary: dict[str, object], softk_comparison: list[dict[str, object]]) -> list[dict[str, object]]:
    full = summary["cutoffs"]["full"]
    deltas = [row["hierk_minus_softk"] for row in softk_comparison]
    rows = [
        metric("alpha", summary["alpha"], "decimal exponent", "Fitted hierarchical prior exponent in log k_i ~ Normal(alpha log V_i, sigma^2)."),
        metric("sigma", summary["sigma"], "decimal", "Fitted hierarchical prior spread."),
        metric("sigma_lower_bound", 0.05, "decimal", "Lower bound used by historical hierarchical pooling fit."),
        metric("sigma_at_lower_bound", abs(summary["sigma"] - 0.05) < 1e-12, "boolean", "Whether fitted sigma is pinned at the lower bound."),
        metric("hierk_beats_moe_count", full["hier_beats_moe"], "X/25 count", "Full held-out NLL count where hierarchical-k beats MOEZipf."),
        metric("hierk_beats_zm_count", full["hier_beats_zm"], "X/25 count", "Full held-out NLL count where hierarchical-k beats ZM."),
        metric("hierk_beats_zipf_count", full["hier_beats_zipf"], "X/25 count", "Full held-out NLL count where hierarchical-k beats Zipf."),
        metric("median_hierk_minus_moe", full["median_hier_minus_moe"], "signed average NLL delta", "Median full held-out NLL delta, hierarchical-k minus MOEZipf."),
        metric("median_hierk_minus_zm", full["median_hier_minus_zm"], "signed average NLL delta", "Median full held-out NLL delta, hierarchical-k minus ZM."),
        metric("hierk_step2_help_count", full["step2_help_count"], "X/25 count", "Count where step-2 search helps on full hierarchical-k residuals."),
        metric("hierk_beats_softk_count", sum(row["hierk_beats_softk"] for row in softk_comparison), "X/25 count", "Cross-experiment comparison against 7a canonical v4 soft-k Table 4 values."),
        metric("median_hierk_minus_softk", statistics.median(deltas), "signed average NLL delta", "Cross-experiment median full held-out NLL delta, hierarchical-k minus 7a canonical soft-k."),
        metric(
            "hierk_vs_softk_summary",
            f"hierk beats 7a canonical soft-k on {sum(row['hierk_beats_softk'] for row in softk_comparison)}/25; median delta {statistics.median(deltas):.12g}",
            "text summary",
            "This data-level result conflicts with manuscript wording that hierarchical pooling does not improve over per-corpus soft-k; see README AUDIT.",
        ),
    ]
    return rows


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    table_rows = read_csv(cfg.SOURCE_TABLE)
    softk_rows = read_csv(cfg.SOFTK_CANONICAL_TABLE)
    summary = json.loads(cfg.SOURCE_SUMMARY.read_text(encoding="utf-8"))

    hierk_summary = build_hierk_summary(summary)
    history = build_history(summary)
    head_rows = build_head_window_rows(table_rows)
    softk_comparison = build_softk_comparison(head_rows, softk_rows)
    aggregates = build_aggregates(summary, softk_comparison)

    write_csv(output_dir / "hierk_summary.csv", hierk_summary, ["model", "alpha", "sigma", "sigma_lower_bound", "sigma_at_lower_bound", "iteration_count", "source"])
    write_csv(output_dir / "hierk_iteration_history.csv", history, ["iteration", "alpha", "sigma", "sigma_at_lower_bound"])
    write_csv(
        output_dir / "hierk_head_window_per_corpus.csv",
        head_rows,
        ["slug", "corpus", "cutoff", "winner_family", "zipf_avg_nll", "zm_avg_nll", "moe_avg_nll", "hierk_avg_nll", "hierk_minus_moe", "hierk_minus_zm", "step2_helps", "step2_gain", "canonical_source"],
    )
    write_csv(output_dir / "hierk_vs_softk_comparison.csv", softk_comparison, ["slug", "corpus", "hierk_full_avg_nll", "softk_full_avg_nll", "hierk_minus_softk", "hierk_beats_softk", "softk_source"])
    write_csv(output_dir / "aggregate_statistics.csv", aggregates, ["metric_name", "value", "display_format", "notes"])

    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps(
            {
                "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
                "cross_experiment_dependency": str(cfg.SOFTK_CANONICAL_TABLE.relative_to(cfg.ROOT)),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    provenance = cfg.ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for src in [cfg.SOURCE_TABLE, cfg.SOURCE_SUMMARY, cfg.SOURCE_REPORT]:
        shutil.copy2(src, provenance / src.name)

    manifest = {
        "experiment_id": "7e_hierarchical_k_pooling",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
        "cross_experiment_dependencies": [str(cfg.SOFTK_CANONICAL_TABLE.relative_to(cfg.ROOT))],
        "status": "complete_with_audit_flag",
        "outputs": {
            "hierk_summary.csv": {"rows": len(hierk_summary)},
            "hierk_iteration_history.csv": {"rows": len(history)},
            "hierk_head_window_per_corpus.csv": {"rows": len(head_rows)},
            "hierk_vs_softk_comparison.csv": {"rows": len(softk_comparison)},
            "aggregate_statistics.csv": {"rows": len(aggregates), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "audit_flags": [
            "Manuscript wording says hierarchical pooling does not improve over soft-k, but comparison to 7a canonical v4 soft-k shows hierk beats soft-k on 20/25 with median delta -0.000912303553.",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote 7e outputs: {len(head_rows)} head-window rows, {len(softk_comparison)} soft-k comparison rows, {len(aggregates)} aggregate rows")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()
