import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
EXPERIMENT_DIR = SCRIPT_DIR.parent
ROOT = EXPERIMENT_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

import eml_zipf_enriched_search as enriched
import source_config as cfg


TABLE1_FIELDS = [
    "corpus",
    "tokens",
    "vocabulary_size",
    "zm_a",
    "zm_b",
    "zm_c",
    "zm_rmse",
    "step2_winner_expression",
    "step2_winner_display",
    "step2_rmse",
    "rmse_improvement_absolute",
    "rmse_improvement_pct",
    "runner_up_expression",
    "runner_up_rmse",
    "winner_minus_runner_rmse",
    "winner_family",
    "contains_is_in_step2_beam",
    "contains_exp_in_step2_beam",
]

BEAM_FIELDS = [
    "corpus",
    "beam_rank",
    "expression",
    "expression_display",
    "rmse",
    "winner_family",
    "uses_new_operator",
]

DIFF_FIELDS = [
    "corpus",
    "field",
    "new_value",
    "historical_value",
    "abs_diff",
    "matches",
]

AGG_FIELDS = ["metric_name", "value", "display_format", "notes"]
SUMMARY_FIELDS = ["metric_name", "value", "notes"]

NUMERIC_COMPARE_FIELDS = [
    "tokens",
    "vocabulary_size",
    "zm_a",
    "zm_b",
    "zm_c",
    "zm_rmse",
    "step2_rmse",
    "rmse_improvement_absolute",
    "rmse_improvement_pct",
    "runner_up_rmse",
    "winner_minus_runner_rmse",
]

TEXT_COMPARE_FIELDS = [
    "step2_winner_expression",
    "step2_winner_display",
    "runner_up_expression",
    "winner_family",
    "contains_is_in_step2_beam",
    "contains_exp_in_step2_beam",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Run consolidated experiment 1a.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=cfg.OUTPUT_DIR,
        help="Directory where experiment outputs will be written.",
    )
    return parser.parse_args()


def step_summary_for(step_summaries, target_step):
    for item in step_summaries:
        if int(item["step"]) == int(target_step):
            return item
    raise ValueError(f"Missing step {target_step} in search summary")


def build_step2_rows(entry, bundle, step2_summary):
    top_candidates = step2_summary["top_candidates"]
    winner = top_candidates[0]
    runner = top_candidates[1]

    beam_exprs = {candidate["expr"] for candidate in top_candidates}
    zm_rmse = float(bundle["zm"]["rmse_full"])
    step2_rmse = float(winner["rmse"])
    improvement_abs = zm_rmse - step2_rmse
    improvement_pct = 100.0 * improvement_abs / zm_rmse if zm_rmse else float("nan")
    winner_runner_gap = float(runner["rmse"]) - step2_rmse

    per_corpus = {
        "corpus": entry["name"],
        "tokens": bundle["corpus"]["token_count"],
        "vocabulary_size": bundle["corpus"]["unique_words"],
        "zm_a": float(bundle["zm"]["a"]),
        "zm_b": float(bundle["zm"]["b"]),
        "zm_c": float(bundle["zm"]["c"]),
        "zm_rmse": zm_rmse,
        "step2_winner_expression": winner["expr"],
        "step2_winner_display": cfg.display_expression(winner["expr"], winner["math"]),
        "step2_rmse": step2_rmse,
        "rmse_improvement_absolute": improvement_abs,
        "rmse_improvement_pct": improvement_pct,
        "runner_up_expression": runner["expr"],
        "runner_up_rmse": float(runner["rmse"]),
        "winner_minus_runner_rmse": winner_runner_gap,
        "winner_family": cfg.winner_family(winner["expr"]),
        "contains_is_in_step2_beam": cfg.STEP2_OURS_EXPR in beam_exprs,
        "contains_exp_in_step2_beam": cfg.STEPL2_LOWC_EXPR in beam_exprs,
    }

    beam_rows = []
    for beam_rank, candidate in enumerate(top_candidates[:10], start=1):
        beam_rows.append(
            {
                "corpus": entry["name"],
                "beam_rank": beam_rank,
                "expression": candidate["expr"],
                "expression_display": cfg.display_expression(candidate["expr"], candidate["math"]),
                "rmse": float(candidate["rmse"]),
                "winner_family": cfg.winner_family(candidate["expr"]),
                "uses_new_operator": False,
            }
        )

    return per_corpus, beam_rows


def rerun_corpus(entry):
    corpus_path = cfg.corpus_path(entry)
    if not corpus_path.exists():
        raise FileNotFoundError(f"Missing corpus file: {corpus_path}")

    enriched.EXP_CLAMP = float(cfg.EXP_CLAMP)
    enriched.VALUE_ABS_LIMIT = float(cfg.VALUE_ABS_LIMIT)

    bundle = enriched.build_target_bundle(
        cfg.DUMMY_CORPUS_URL,
        corpus_path,
        cfg.SAMPLE_POINTS,
        cfg.X_LOW,
        cfg.X_HIGH,
    )
    zm_target = bundle["y"] - bundle["zm"]["prediction_sample"]
    zm_result = enriched.run_search(
        bundle["x"],
        zm_target,
        cfg.BEAM_WIDTH,
        cfg.REPLAY_MAX_STEPS,
        cfg.KEEP_ALL_UNTIL_STEP,
        cfg.DIVERSITY_WEIGHT,
        cfg.CONSTANT_VARIANCE_THRESHOLD,
    )
    step2_summary = step_summary_for(zm_result["steps"], 2)
    return build_step2_rows(entry, bundle, step2_summary)


def historical_rows(entry):
    summary_path = cfg.historical_summary_path(entry)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    step2_summary = step_summary_for(summary["zm_search"]["step_summary"], 2)

    top_candidates = step2_summary["top_candidates"]
    winner = top_candidates[0]
    runner = top_candidates[1]
    zm_rmse = float(summary["zm_baseline"]["rmse_full"])
    step2_rmse = float(winner["rmse"])
    improvement_abs = zm_rmse - step2_rmse
    improvement_pct = 100.0 * improvement_abs / zm_rmse if zm_rmse else float("nan")
    winner_runner_gap = float(runner["rmse"]) - step2_rmse
    beam_exprs = {candidate["expr"] for candidate in top_candidates}

    per_corpus = {
        "corpus": entry["name"],
        "tokens": summary["corpus"]["token_count"],
        "vocabulary_size": summary["corpus"]["unique_words"],
        "zm_a": float(summary["zm_baseline"]["a"]),
        "zm_b": float(summary["zm_baseline"]["b"]),
        "zm_c": float(summary["zm_baseline"]["c"]),
        "zm_rmse": zm_rmse,
        "step2_winner_expression": winner["expr"],
        "step2_winner_display": cfg.display_expression(winner["expr"], winner["math"]),
        "step2_rmse": step2_rmse,
        "rmse_improvement_absolute": improvement_abs,
        "rmse_improvement_pct": improvement_pct,
        "runner_up_expression": runner["expr"],
        "runner_up_rmse": float(runner["rmse"]),
        "winner_minus_runner_rmse": winner_runner_gap,
        "winner_family": cfg.winner_family(winner["expr"]),
        "contains_is_in_step2_beam": cfg.STEP2_OURS_EXPR in beam_exprs,
        "contains_exp_in_step2_beam": cfg.STEPL2_LOWC_EXPR in beam_exprs,
    }

    beam_rows = []
    for beam_rank, candidate in enumerate(top_candidates[:10], start=1):
        beam_rows.append(
            {
                "corpus": entry["name"],
                "beam_rank": beam_rank,
                "expression": candidate["expr"],
                "expression_display": cfg.display_expression(candidate["expr"], candidate["math"]),
                "rmse": float(candidate["rmse"]),
                "winner_family": cfg.winner_family(candidate["expr"]),
                "uses_new_operator": False,
            }
        )
    return per_corpus, beam_rows


def historical_terminal_metrics(entry):
    summary_path = cfg.historical_summary_path(entry)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    step_summaries = summary["zm_search"]["step_summary"]
    final_step = step_summaries[-1]
    winner = final_step["top_candidates"][0]
    runner = final_step["top_candidates"][1]
    zm_rmse = float(summary["zm_baseline"]["rmse_full"])
    terminal_rmse = float(summary["zm_search"]["best"]["composite_rmse"])
    terminal_improvement_abs = zm_rmse - terminal_rmse
    terminal_improvement_pct = 100.0 * terminal_improvement_abs / zm_rmse if zm_rmse else float("nan")
    final_gap = float(runner["rmse"]) - float(winner["rmse"])
    return {
        "corpus": entry["name"],
        "terminal_step_count": int(summary["zm_search"]["best"]["step"]),
        "terminal_improvement_absolute": terminal_improvement_abs,
        "terminal_improvement_pct": terminal_improvement_pct,
        "terminal_winner_runner_gap": final_gap,
    }


def compute_aggregates(per_corpus_rows, terminal_rows):
    winner_gaps = [float(row["winner_minus_runner_rmse"]) for row in per_corpus_rows]
    improvement_pcts = [float(row["rmse_improvement_pct"]) for row in per_corpus_rows]
    positive_step2_improvement_pcts = [pct for pct in improvement_pcts if pct > 0.0]
    terminal_improvement_pcts = [float(row["terminal_improvement_pct"]) for row in terminal_rows]
    terminal_winner_runner_gaps = [float(row["terminal_winner_runner_gap"]) for row in terminal_rows]
    c_by_family = {}
    for row in per_corpus_rows:
        c_by_family.setdefault(row["winner_family"], []).append(float(row["zm_c"]))

    exp_upper = max(c_by_family.get("exp", [float("nan")]))
    is_lower = min(c_by_family.get("is", [float("nan")]))
    transition_summary = (
        f"exp_dominant_below={exp_upper:.6f};"
        f"is_dominant_above={is_lower:.6f};"
        f"band_width={is_lower - exp_upper:.6f}"
    )

    both_present = sum(
        bool(row["contains_is_in_step2_beam"]) and bool(row["contains_exp_in_step2_beam"])
        for row in per_corpus_rows
    )

    metrics = [
        ("english_corpus_count", len(per_corpus_rows), "int", "Canonical English corpus count."),
        ("winner_family_is_count", sum(row["winner_family"] == "is" for row in per_corpus_rows), "int", "Count of IS step-2 winners."),
        ("winner_family_exp_count", sum(row["winner_family"] == "exp" for row in per_corpus_rows), "int", "Count of exponential step-2 winners."),
        ("winner_family_other_count", sum(row["winner_family"] not in {"is", "exp"} for row in per_corpus_rows), "int", "Count of non-Bregman winners."),
        ("winner_runner_gap_mean", statistics.fmean(winner_gaps), "float6", "Mean winner minus runner-up RMSE gap."),
        ("winner_runner_gap_median", statistics.median(winner_gaps), "float6", "Median winner minus runner-up RMSE gap."),
        ("winner_family_gap_below_0p01_count", sum(gap < 0.01 for gap in winner_gaps), "int", "Count of corpora with winner-runner gap below 0.01."),
        ("high_c_is_winner_count", sum(row["winner_family"] == "is" for row in per_corpus_rows), "int", "Count of IS winners in the high-c regime."),
        ("low_c_exp_winner_count", sum(row["winner_family"] == "exp" for row in per_corpus_rows), "int", "Count of exponential winners in the low-c regime."),
        ("c_transition_band_summary", transition_summary, "string", "Empirical c-band summary for the IS/exp split."),
        ("step2_improvement_pct_min", min(improvement_pcts), "pct1", "Minimum percentage RMSE improvement over single ZM."),
        ("step2_improvement_pct_max", max(improvement_pcts), "pct1", "Maximum percentage RMSE improvement over single ZM."),
        ("high_c_is_lower_band", is_lower, "float6", "Lower c bound of the IS-dominant region."),
        ("low_c_exp_upper_band", exp_upper, "float6", "Upper c bound of the exp-dominant region."),
        ("both_bregman_forms_present_in_step2_beam_count", both_present, "int", "Count of corpora where both Bregman generators appear in the canonical step-2 beam."),
        ("mean_step2_winner_runner_gap", statistics.fmean(winner_gaps), "float6", "Alias used by the manuscript claim map."),
        ("step2_all_improve_count", sum(pct > 0.0 for pct in improvement_pcts), "int", "Count of corpora where the canonical step-2 winner beats single ZM."),
        ("step2_improvement_pct_mean", statistics.fmean(improvement_pcts), "pct3", "Mean step-2 percentage improvement over single ZM across all 25 corpora."),
        ("step2_improvement_pct_median", statistics.median(improvement_pcts), "pct3", "Median step-2 percentage improvement over single ZM across all 25 corpora."),
        ("step2_positive_improve_count", len(positive_step2_improvement_pcts), "int", "Count of corpora where step-2 improvement is positive."),
        ("step2_positive_improvement_pct_min", min(positive_step2_improvement_pcts), "pct3", "Minimum positive step-2 percentage improvement over single ZM."),
        ("step2_positive_improvement_pct_max", max(positive_step2_improvement_pcts), "pct3", "Maximum positive step-2 percentage improvement over single ZM."),
        ("step2_positive_improvement_pct_mean", statistics.fmean(positive_step2_improvement_pcts), "pct3", "Mean positive-only step-2 percentage improvement over single ZM."),
        ("step2_positive_improvement_pct_median", statistics.median(positive_step2_improvement_pcts), "pct3", "Median positive-only step-2 percentage improvement over single ZM."),
        ("step2_winner_runner_gap_mean", statistics.fmean(winner_gaps), "float6", "Verified mean step-2 winner minus runner-up RMSE gap across all 25 corpora."),
        ("terminal_all_improve_count", sum(float(row["terminal_improvement_absolute"]) > 0.0 for row in terminal_rows), "int", "Count of corpora where the terminal enriched-search winner beats single ZM."),
        ("terminal_improvement_pct_min", min(terminal_improvement_pcts), "pct3", "Minimum terminal enriched-search percentage improvement over single ZM."),
        ("terminal_improvement_pct_max", max(terminal_improvement_pcts), "pct3", "Maximum terminal enriched-search percentage improvement over single ZM."),
        ("terminal_improvement_pct_mean", statistics.fmean(terminal_improvement_pcts), "pct3", "Mean terminal enriched-search percentage improvement over single ZM."),
        ("terminal_improvement_pct_median", statistics.median(terminal_improvement_pcts), "pct3", "Median terminal enriched-search percentage improvement over single ZM."),
        ("terminal_winner_runner_gap_mean", statistics.fmean(terminal_winner_runner_gaps), "float6", "Mean terminal-step winner minus runner-up RMSE gap across all 25 corpora."),
    ]
    return [
        {
            "metric_name": name,
            "value": value,
            "display_format": display_format,
            "notes": notes,
        }
        for name, value, display_format, notes in metrics
    ]


def compare_rows(new_rows, old_rows):
    diff_rows = []
    max_abs_numeric_diff = 0.0
    numeric_diff_field = ""
    mismatch_count = 0
    by_corpus_new = {row["corpus"]: row for row in new_rows}
    by_corpus_old = {row["corpus"]: row for row in old_rows}

    for corpus, new_row in by_corpus_new.items():
        old_row = by_corpus_old[corpus]
        for field in NUMERIC_COMPARE_FIELDS:
            new_value = float(new_row[field])
            old_value = float(old_row[field])
            abs_diff = abs(new_value - old_value)
            if abs_diff > max_abs_numeric_diff:
                max_abs_numeric_diff = abs_diff
                numeric_diff_field = f"{corpus}:{field}"
            diff_rows.append(
                {
                    "corpus": corpus,
                    "field": field,
                    "new_value": new_value,
                    "historical_value": old_value,
                    "abs_diff": abs_diff,
                    "matches": abs_diff <= 1e-12,
                }
            )
        for field in TEXT_COMPARE_FIELDS:
            new_value = new_row[field]
            old_value = old_row[field]
            matches = new_value == old_value
            mismatch_count += int(not matches)
            diff_rows.append(
                {
                    "corpus": corpus,
                    "field": field,
                    "new_value": new_value,
                    "historical_value": old_value,
                    "abs_diff": "",
                    "matches": matches,
                }
            )
    if max_abs_numeric_diff == 0.0 and not numeric_diff_field:
        numeric_diff_field = "all compared numeric per-corpus fields matched exactly"
    return diff_rows, max_abs_numeric_diff, numeric_diff_field, mismatch_count


def compare_top10(new_rows, old_rows):
    diff_rows = []
    max_abs_numeric_diff = 0.0
    numeric_diff_field = ""
    mismatch_count = 0
    by_key_new = {(row["corpus"], int(row["beam_rank"])): row for row in new_rows}
    by_key_old = {(row["corpus"], int(row["beam_rank"])): row for row in old_rows}

    for key, new_row in sorted(by_key_new.items()):
        old_row = by_key_old[key]
        rmse_diff = abs(float(new_row["rmse"]) - float(old_row["rmse"]))
        if rmse_diff > max_abs_numeric_diff:
            max_abs_numeric_diff = rmse_diff
            numeric_diff_field = f"{key[0]}:rank{key[1]}:rmse"
        diff_rows.append(
            {
                "corpus": key[0],
                "beam_rank": key[1],
                "field": "rmse",
                "new_value": new_row["rmse"],
                "historical_value": old_row["rmse"],
                "abs_diff": rmse_diff,
                "matches": rmse_diff <= 1e-12,
            }
        )
        for field in ("expression", "expression_display", "winner_family", "uses_new_operator"):
            matches = new_row[field] == old_row[field]
            mismatch_count += int(not matches)
            diff_rows.append(
                {
                    "corpus": key[0],
                    "beam_rank": key[1],
                    "field": field,
                    "new_value": new_row[field],
                    "historical_value": old_row[field],
                    "abs_diff": "",
                    "matches": matches,
                }
            )
    if max_abs_numeric_diff == 0.0 and not numeric_diff_field:
        numeric_diff_field = "all compared top-10 numeric fields matched exactly"
    return diff_rows, max_abs_numeric_diff, numeric_diff_field, mismatch_count


def write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    table_rows = []
    beam_rows = []
    historical_table_rows = []
    historical_beam_rows = []
    historical_terminal_rows = []

    for entry in cfg.SEARCHED_CORPORA:
        per_corpus_row, per_corpus_beam = rerun_corpus(entry)
        table_rows.append(per_corpus_row)
        beam_rows.extend(per_corpus_beam)

        historical_row, historical_beam = historical_rows(entry)
        historical_table_rows.append(historical_row)
        historical_beam_rows.extend(historical_beam)
        historical_terminal_rows.append(historical_terminal_metrics(entry))

    aggregate_rows = compute_aggregates(table_rows, historical_terminal_rows)
    per_corpus_diff_rows, table_max_diff, table_max_field, table_mismatch_count = compare_rows(
        table_rows, historical_table_rows
    )
    top10_diff_rows, top10_max_diff, top10_max_field, top10_mismatch_count = compare_top10(
        beam_rows, historical_beam_rows
    )
    diff_summary_rows = [
        {
            "metric_name": "table_row_count",
            "value": len(table_rows),
            "notes": "Expected 25 per-corpus rows.",
        },
        {
            "metric_name": "beam_row_count",
            "value": len(beam_rows),
            "notes": "Expected 250 top-10 beam rows.",
        },
        {
            "metric_name": "table_max_abs_numeric_diff",
            "value": table_max_diff,
            "notes": table_max_field or "no numeric fields compared",
        },
        {
            "metric_name": "table_text_mismatch_count",
            "value": table_mismatch_count,
            "notes": "Count of non-numeric per-corpus field mismatches versus historical bundles.",
        },
        {
            "metric_name": "top10_max_abs_numeric_diff",
            "value": top10_max_diff,
            "notes": top10_max_field or "no numeric fields compared",
        },
        {
            "metric_name": "top10_text_mismatch_count",
            "value": top10_mismatch_count,
            "notes": "Count of non-numeric top-10 beam field mismatches versus historical bundles.",
        },
    ]

    write_csv(output_dir / "table1_per_corpus.csv", TABLE1_FIELDS, table_rows)
    write_csv(output_dir / "table1_step2_beam_top10.csv", BEAM_FIELDS, beam_rows)
    write_csv(output_dir / "aggregate_statistics.csv", AGG_FIELDS, aggregate_rows)
    write_csv(output_dir / "historical_per_corpus_diff.csv", DIFF_FIELDS, per_corpus_diff_rows)
    write_csv(output_dir / "historical_top10_diff.csv", ["corpus", "beam_rank", "field", "new_value", "historical_value", "abs_diff", "matches"], top10_diff_rows)
    write_csv(output_dir / "historical_diff_summary.csv", SUMMARY_FIELDS, diff_summary_rows)

    manifest = {
        "experiment_id": "1a",
        "output_dir": str(output_dir),
        "corpus_count": len(table_rows),
        "beam_row_count": len(beam_rows),
        "historical_table_max_abs_numeric_diff": table_max_diff,
        "historical_table_text_mismatch_count": table_mismatch_count,
        "historical_top10_max_abs_numeric_diff": top10_max_diff,
        "historical_top10_text_mismatch_count": top10_mismatch_count,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote {output_dir / 'table1_per_corpus.csv'}")
    print(f"Wrote {output_dir / 'table1_step2_beam_top10.csv'}")
    print(f"Wrote {output_dir / 'aggregate_statistics.csv'}")
    print(f"Wrote {output_dir / 'historical_per_corpus_diff.csv'}")
    print(f"Wrote {output_dir / 'historical_top10_diff.csv'}")
    print(f"Wrote {output_dir / 'historical_diff_summary.csv'}")


if __name__ == "__main__":
    main()
