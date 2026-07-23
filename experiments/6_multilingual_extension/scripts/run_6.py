from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

import source_config as cfg


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def display_winner(expr: str) -> str:
    if expr == "sub[pow[x,x],sqrt[x]]":
        return "x^x - sqrt(x)"
    if expr == "eml[sub[x,1],eml[x,1]]":
        return "exp(x-1) - x"
    if expr == "sub[sub[x,1],log[x]]":
        return "(x-1) - log(x)"
    if expr == "mul[sub[1,x],sub[1,x]]":
        return "(1-x)^2"
    return expr


def is_bregman_winner(expr: str) -> bool:
    return expr == "eml[sub[x,1],eml[x,1]]"


def migrate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    romance = load_json(cfg.MULTILANG_ROMANCE_SUMMARY_JSON)
    verify = load_json(cfg.MULTILANG_VERIFY_SUMMARY_JSON)
    verify_by_slug = {row["slug"]: row for row in verify["corpora"]}
    rows: list[dict[str, object]] = []
    table_rows: list[dict[str, object]] = []
    for row in romance["rows"]:
        vrow = verify_by_slug[row["slug"]]
        transition_fraction = float(row["transition_fraction"])
        vocab = float(row["vocab_size"])
        k_over_sqrt_v = transition_fraction * math.sqrt(vocab)
        merged = {
            "slug": row["slug"],
            "language": row["language"],
            "corpus": row["corpus"],
            "token_count": row["token_count"],
            "raw_token_count": row["raw_token_count"],
            "vocabulary_size": row["vocab_size"],
            "zm_a": vrow["zm_a"],
            "zm_b": vrow["zm_b"],
            "zm_c": vrow["zm_c"],
            "single_zm_rmse": row["single_zm_rmse"],
            "single_zm_bic": row["single_zm_bic"],
            "moezipf_rmse": row["moezipf_rmse"],
            "moezipf_bic": row["moezipf_bic"],
            "smooth_7param_rmse": row["smooth_7param_rmse"],
            "smooth_7param_bic": row["smooth_7param_bic"],
            "smooth_8param_rmse": row["smooth_8param_rmse"],
            "smooth_8param_bic": row["smooth_8param_bic"],
            "transition_fraction": transition_fraction,
            "step2_winner": row["step2_winner"],
            "step2_winner_display": display_winner(row["step2_winner"]),
            "step2_rmse": row["step2_rmse"],
            "bregman_in_beam": str(row["bregman_in_beam"]),
            "artifact_note": row["artifact_note"],
            "high_bregman_rmse": vrow["high_bregman_rmse"],
            "low_bregman_rmse": vrow["low_bregman_rmse"],
            "euclidean_rmse": vrow["euclidean_rmse"],
            "xpow_rmse": vrow["xpow_rmse"],
            "gap_to_high_bregman": vrow["gap_to_high_bregman"],
            "gap_to_euclidean": vrow["gap_to_euclidean"],
            "k_over_sqrt_v": k_over_sqrt_v,
        }
        rows.append(merged)
        table_rows.append(
            {
                "language": row["language"],
                "corpus": row["corpus"],
                "vocabulary_size": row["vocab_size"],
                "zm_c": vrow["zm_c"],
                "single_zm_rmse": row["single_zm_rmse"],
                "smooth_rmse": row["smooth_8param_rmse"],
                "step2_winner_display": display_winner(row["step2_winner"]),
                "k_over_sqrt_v": k_over_sqrt_v,
            }
        )

    transition_values = [float(row["transition_fraction"]) for row in rows]
    smooth_beats = sum(float(row["smooth_8param_rmse"]) < float(row["single_zm_rmse"]) for row in rows)
    bregman_winners = sum(is_bregman_winner(str(row["step2_winner"])) for row in rows)
    low_named = {"latin_gallic_wars", "french_les_miserables", "spanish_don_quixote", "dutch_max_havelaar"}
    low_named_count = sum(float(row["zm_c"]) < cfg.LOW_C_THRESHOLD for row in rows if row["slug"] in low_named)
    very_low_count = sum(float(row["zm_c"]) < cfg.VERY_LOW_C_THRESHOLD for row in rows)
    near_half_count = sum(cfg.NEAR_HALF_LOW <= float(row["transition_fraction"]) <= cfg.NEAR_HALF_HIGH for row in rows)
    aggregate_rows = [
        ("non_english_corpus_count", len(rows), "integer", "Number of non-English corpora in the canonical multilingual extension."),
        ("smooth_beats_single_zm_count", smooth_beats, "integer", "Count where smooth 8-parameter historical fit has lower RMSE than single ZM."),
        ("transition_fraction_min", min(transition_values), "decimal_6", "Minimum transition fraction across the 7 non-English corpora."),
        ("transition_fraction_max", max(transition_values), "decimal_6", "Maximum transition fraction across the 7 non-English corpora."),
        ("transition_fraction_near_half_count", near_half_count, "integer", f"Count with transition fraction in [{cfg.NEAR_HALF_LOW}, {cfg.NEAR_HALF_HIGH}]."),
        ("multilang_bregman_winner_count", bregman_winners, "integer", "Count where the step-2 winner is exp(x-1)-x in the historical table."),
        ("latin_french_spanish_dutch_c_lt_5_count", low_named_count, "integer", "Count among Latin/French/Spanish/Dutch with fitted ZM c < 5."),
        ("very_low_c_multilang_count", very_low_count, "integer", f"Count with fitted ZM c < {cfg.VERY_LOW_C_THRESHOLD}."),
    ]
    aggregate = [
        {"metric_name": name, "value": repr(value), "display_format": fmt, "notes": notes}
        for name, value, fmt, notes in aggregate_rows
    ]

    per_fields = [
        "slug",
        "language",
        "corpus",
        "token_count",
        "raw_token_count",
        "vocabulary_size",
        "zm_a",
        "zm_b",
        "zm_c",
        "single_zm_rmse",
        "single_zm_bic",
        "moezipf_rmse",
        "moezipf_bic",
        "smooth_7param_rmse",
        "smooth_7param_bic",
        "smooth_8param_rmse",
        "smooth_8param_bic",
        "transition_fraction",
        "step2_winner",
        "step2_winner_display",
        "step2_rmse",
        "bregman_in_beam",
        "artifact_note",
        "high_bregman_rmse",
        "low_bregman_rmse",
        "euclidean_rmse",
        "xpow_rmse",
        "gap_to_high_bregman",
        "gap_to_euclidean",
        "k_over_sqrt_v",
    ]
    table_fields = [
        "language",
        "corpus",
        "vocabulary_size",
        "zm_c",
        "single_zm_rmse",
        "smooth_rmse",
        "step2_winner_display",
        "k_over_sqrt_v",
    ]
    write_csv(output_dir / "multilang_per_corpus.csv", rows, per_fields)
    write_csv(output_dir / "table3_multilang.csv", table_rows, table_fields)
    write_csv(output_dir / "aggregate_statistics.csv", aggregate, ["metric_name", "value", "display_format", "notes"])
    manifest = {
        "experiment_id": "6_multilingual_extension",
        "migration_type": "historical_bundle_consolidation",
        "source_bundles": [
            str(cfg.MULTILANG_ROMANCE_DIR.relative_to(cfg.ROOT)),
            str(cfg.MULTILANG_VERIFY_DIR.relative_to(cfg.ROOT)),
        ],
        "outputs": {
            "multilang_per_corpus.csv": {"rows": len(rows), "schema": per_fields},
            "table3_multilang.csv": {"rows": len(table_rows), "schema": table_fields},
            "aggregate_statistics.csv": {"rows": len(aggregate), "schema": ["metric_name", "value", "display_format", "notes"]},
        },
        "notes": [
            "Uses historical smooth_8param values from the multilingual bundle; decoupled-erf canon applies only to 3a/3c/4/5c per migration adjustment.",
            "T1 branch-shift result is not part of this Phase 2 migration output.",
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
