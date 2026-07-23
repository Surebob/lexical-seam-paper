from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from collections import Counter
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


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def as_float(value: object) -> float:
    return float(value)


def cutoff_rank(cutoff: str, vocab_size: int) -> int:
    if cutoff == "full":
        return vocab_size
    return int(cutoff.removeprefix("top"))


def metric(name: str, value: object, display_format: str, notes: str) -> dict[str, object]:
    return {
        "metric_name": name,
        "value": value,
        "display_format": display_format,
        "notes": notes,
    }


def modal_winner(rows: list[dict[str, object]]) -> str:
    counts = Counter(str(row["winner_family"]) for row in rows)
    return sorted(counts.items(), key=lambda item: (-item[1], cfg.MODEL_ORDER.index(item[0])))[0][0]


def build_per_corpus_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for item in summary["rows"]:
        slug = item["slug"]
        token_count = int(item["token_count"])
        vocab_size = int(item["vocab_size"])
        for cutoff in cfg.CUTOFF_ORDER:
            heldout = item["heldout"][cutoff]
            avg_nll = heldout["avg_nll"]
            step2 = item["step2"][cutoff]
            out.append(
                {
                    "slug": slug,
                    "corpus": item["name"],
                    "token_count": token_count,
                    "vocabulary_size": vocab_size,
                    "cutoff": cutoff,
                    "cutoff_rank": cutoff_rank(cutoff, vocab_size),
                    "zipf_avg_nll": avg_nll["zipf"],
                    "zm_avg_nll": avg_nll["zm"],
                    "moe_avg_nll": avg_nll["moe"],
                    "softk_avg_nll": avg_nll["softk"],
                    "winner_family": heldout["winner"],
                    "softk_minus_moe": heldout["softk_minus_moe"],
                    "softk_minus_zm": heldout["softk_minus_zm"],
                    "softk_minus_zipf": heldout["softk_minus_zipf"],
                    "step2_winner_expression": step2["expr"],
                    "step2_gain": step2["gain"],
                    "step2_helps": step2["helps"],
                    "canonical_source": "results/zipf_angle1_head_windows_splitfit/summary.json",
                }
            )
    return out


def build_aggregate_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    aggregates: list[dict[str, object]] = []
    by_cutoff = {cutoff: [row for row in rows if row["cutoff"] == cutoff] for cutoff in cfg.CUTOFF_ORDER}
    for cutoff in cfg.CUTOFF_ORDER:
        sub = by_cutoff[cutoff]
        winner_counts = Counter(str(row["winner_family"]) for row in sub)
        for model in cfg.MODEL_ORDER:
            aggregates.append(
                metric(
                    f"{cutoff}_winner_count_{model}",
                    winner_counts.get(model, 0),
                    "X/25 count",
                    f"Held-out NLL winner count for {cutoff}.",
                )
            )
        aggregates.append(
            metric(
                f"{cutoff}_winner_count_summary",
                json.dumps({model: winner_counts.get(model, 0) for model in cfg.MODEL_ORDER}, sort_keys=True),
                "JSON count map",
                f"Held-out NLL winner counts for {cutoff} in Zipf, ZM, MOE, soft-k order.",
            )
        )
        aggregates.append(
            metric(
                f"{cutoff}_modal_winner",
                modal_winner(sub),
                "categorical model label",
                f"Most frequent held-out winner for {cutoff}.",
            )
        )
        aggregates.append(
            metric(
                f"{cutoff}_softk_beats_moe_count",
                sum(as_float(row["softk_minus_moe"]) < 0 for row in sub),
                "X/25 count",
                f"Count where soft-k average held-out NLL is below MOEZipf for {cutoff}.",
            )
        )
        aggregates.append(
            metric(
                f"{cutoff}_softk_beats_zm_count",
                sum(as_float(row["softk_minus_zm"]) < 0 for row in sub),
                "X/25 count",
                f"Count where soft-k average held-out NLL is below ZM for {cutoff}.",
            )
        )
        aggregates.append(
            metric(
                f"{cutoff}_softk_step2_help_count",
                sum(as_bool(row["step2_helps"]) for row in sub),
                "X/25 count",
                f"Count where step-2 symbolic search improves the soft-k residual under {cutoff} scoring.",
            )
        )
    aggregates.append(
        metric(
            "max_head_window_rank",
            cfg.MAX_HEAD_WINDOW_RANK,
            "integer rank cutoff",
            "Largest finite top-K cutoff used before the full-vocabulary condition.",
        )
    )
    return aggregates


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(cfg.SOURCE_SUMMARY.read_text(encoding="utf-8"))
    rows = build_per_corpus_rows(summary)
    aggregates = build_aggregate_rows(rows)

    per_corpus_fields = [
        "slug",
        "corpus",
        "token_count",
        "vocabulary_size",
        "cutoff",
        "cutoff_rank",
        "zipf_avg_nll",
        "zm_avg_nll",
        "moe_avg_nll",
        "softk_avg_nll",
        "winner_family",
        "softk_minus_moe",
        "softk_minus_zm",
        "softk_minus_zipf",
        "step2_winner_expression",
        "step2_gain",
        "step2_helps",
        "canonical_source",
    ]
    write_csv(output_dir / "head_window_per_corpus.csv", rows, per_corpus_fields)
    write_csv(output_dir / "aggregate_statistics.csv", aggregates, ["metric_name", "value", "display_format", "notes"])

    cfg.ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.ARCHIVE_DIR / "source_bundles.json").write_text(
        json.dumps({"source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))]}, indent=2),
        encoding="utf-8",
    )
    provenance = cfg.ARCHIVE_DIR / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    for src in [cfg.SOURCE_TABLE, cfg.SOURCE_SUMMARY, cfg.SOURCE_REPORT]:
        shutil.copy2(src, provenance / src.name)

    manifest = {
        "experiment_id": "7b_pmf_head_window_evaluation",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [str(cfg.SOURCE_BUNDLE.relative_to(cfg.ROOT))],
        "status": "complete",
        "outputs": {
            "head_window_per_corpus.csv": {"rows": len(rows), "schema": per_corpus_fields},
            "aggregate_statistics.csv": {
                "rows": len(aggregates),
                "schema": ["metric_name", "value", "display_format", "notes"],
            },
        },
        "claim_map_rows_satisfied": [
            "line 409 head-window held-out winner counts and step-2 help counts",
            "line 504 modal winner summary by head window",
            "line 536 max head-window rank protocol constant",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote 7b outputs: {len(rows)} head-window rows, {len(aggregates)} aggregate rows")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=cfg.OUTPUT_DIR)
    args = parser.parse_args()
    migrate(args.output_dir)


if __name__ == "__main__":
    main()
