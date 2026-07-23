import csv
import importlib.util
import json
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

import source_config as cfg


COMMON_PATH = ROOT / "zipf_analysis_common.py"
CORRECT_MODEL_PATH = ROOT / "zipf_correct_model.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
RELAXED_PATH = ROOT / "zipf_reranked_model_all_corpora_relaxed.py"

SMOOTH_FIELDS = [
    "slug",
    "corpus",
    "token_count",
    "vocab_size",
    "single_zm_rmse",
    "piecewise_rmse",
    "step2_rmse",
    "step2_expression",
    "reranked_rmse",
    "beats_single_zm",
    "beats_piecewise",
    "beats_step2",
    "best_start_index",
    "best_nfev",
    "a1",
    "b1",
    "c1",
    "a2",
    "b2",
    "c2",
    "k",
    "w",
    "transition_fraction",
    "tried_count",
]

RELAXED_FIELDS = [
    "slug",
    "corpus",
    "original_rmse",
    "relaxed_rmse",
    "rmse_delta_vs_original",
    "original_transition_fraction",
    "relaxed_transition_fraction",
    "transition_fraction_delta",
    "original_k",
    "relaxed_k",
    "original_w",
    "relaxed_w",
    "best_start_index",
    "best_nfev",
    "tried_count",
    "relaxed_a1",
    "relaxed_b1",
    "relaxed_c1",
    "relaxed_a2",
    "relaxed_b2",
    "relaxed_c2",
]

AGG_FIELDS = ["metric_name", "value", "display_format", "notes"]
DIFF_FIELDS = ["slug", "field", "new_value", "historical_value", "abs_diff", "matches"]
SUMMARY_FIELDS = ["metric_name", "value", "notes"]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "experiment_3a_common")
correct_model = load_module(CORRECT_MODEL_PATH, "experiment_3a_correct_model")
reranked = load_module(RERANKED_PATH, "experiment_3a_reranked")
relaxed_mod = load_module(RELAXED_PATH, "experiment_3a_relaxed")


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Run consolidated experiment 3a.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=cfg.OUTPUT_DIR,
        help="Directory where experiment outputs will be written.",
    )
    return parser.parse_args()


def flatten_params(params: dict, prefix: str = ""):
    return {
        f"{prefix}a1": float(params["a1"]),
        f"{prefix}b1": float(params["b1"]),
        f"{prefix}c1": float(params["c1"]),
        f"{prefix}a2": float(params["a2"]),
        f"{prefix}b2": float(params["b2"]),
        f"{prefix}c2": float(params["c2"]),
        f"{prefix}k": float(params["k"]),
        f"{prefix}w": float(params["w"]),
        f"{prefix}transition_fraction": float(params["transition_fraction"]),
    }


def historical_rows(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {row["slug"]: row for row in payload["rows"]}


def normalize_historical_smooth_rows(rows_by_slug):
    normalized = {}
    for slug, row in rows_by_slug.items():
        params = row["params"]
        normalized[slug] = {
            "slug": slug,
            "corpus": row["name"],
            "token_count": int(row["token_count"]),
            "vocab_size": int(row["vocab_size"]),
            "single_zm_rmse": float(row["single_zm_rmse"]),
            "piecewise_rmse": float(row["piecewise_rmse"]),
            "step2_rmse": float(row["step2_rmse"]),
            "step2_expression": row["step2_expr"],
            "reranked_rmse": float(row["reranked_rmse"]),
            "beats_single_zm": bool(row["beats_single_zm"]),
            "beats_piecewise": bool(row["beats_piecewise"]),
            "beats_step2": bool(row["beats_step2"]),
            "best_start_index": int(row["best_start_index"]),
            "best_nfev": int(row["best_nfev"]),
            "tried_count": int(row["tried_count"]),
            "a1": float(params["a1"]),
            "b1": float(params["b1"]),
            "c1": float(params["c1"]),
            "a2": float(params["a2"]),
            "b2": float(params["b2"]),
            "c2": float(params["c2"]),
            "k": float(row["k"]),
            "w": float(row["w"]),
            "transition_fraction": float(row["transition_fraction"]),
        }
    return normalized


def normalize_historical_relaxed_rows(rows_by_slug):
    normalized = {}
    for slug, row in rows_by_slug.items():
        params = row["params"]
        normalized[slug] = {
            "slug": slug,
            "corpus": row["name"],
            "original_rmse": float(row["original_rmse"]),
            "relaxed_rmse": float(row["relaxed_rmse"]),
            "rmse_delta_vs_original": float(row["rmse_delta_vs_original"]),
            "original_transition_fraction": float(row["original_transition_fraction"]),
            "relaxed_transition_fraction": float(row["relaxed_transition_fraction"]),
            "transition_fraction_delta": float(row["transition_fraction_delta"]),
            "original_k": float(row["original_k"]),
            "relaxed_k": float(row["relaxed_k"]),
            "original_w": float(row["original_w"]),
            "relaxed_w": float(row["relaxed_w"]),
            "best_start_index": int(row["best_start_index"]),
            "best_nfev": int(row["best_nfev"]),
            "tried_count": int(row["tried_count"]),
            "relaxed_a1": float(params["a1"]),
            "relaxed_b1": float(params["b1"]),
            "relaxed_c1": float(params["c1"]),
            "relaxed_a2": float(params["a2"]),
            "relaxed_b2": float(params["b2"]),
            "relaxed_c2": float(params["c2"]),
        }
    return normalized


def build_smooth_row(spec: dict):
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    summary = common.load_enriched_summary(spec)
    piecewise = correct_model.fit_piecewise_zm(dataset, cfg.PIECEWISE_BREAKPOINT_RANK)
    best, tried = reranked.fit_reranked_model(dataset, n_starts=cfg.N_RANDOM_STARTS, max_nfev=cfg.MAX_NFEV)
    step2 = common.get_step2_candidate(summary)
    params = reranked.summarize_params(best["params"], dataset["unique_words"])
    row = {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "single_zm_rmse": float(summary["zm_baseline"]["rmse_full"]),
        "piecewise_rmse": float(piecewise["rmse"]),
        "step2_rmse": float(step2["rmse"]),
        "step2_expression": step2["expr"],
        "reranked_rmse": float(best["rmse"]),
        "beats_single_zm": bool(best["rmse"] < summary["zm_baseline"]["rmse_full"]),
        "beats_piecewise": bool(best["rmse"] < piecewise["rmse"]),
        "beats_step2": bool(best["rmse"] < step2["rmse"]),
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "tried_count": len(tried),
    }
    row.update(flatten_params(params))
    return row


def build_relaxed_row(spec: dict, idx: int, smooth_row: dict):
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    best, tried = relaxed_mod.fit_reranked_model_relaxed(
        dataset,
        n_starts=cfg.N_RANDOM_STARTS,
        max_nfev=cfg.MAX_NFEV,
        seed=cfg.RELAXED_SEED_BASE + idx,
    )
    params = reranked.summarize_params(best["params"], dataset["unique_words"])
    row = {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "original_rmse": float(smooth_row["reranked_rmse"]),
        "relaxed_rmse": float(best["rmse"]),
        "rmse_delta_vs_original": float(best["rmse"] - smooth_row["reranked_rmse"]),
        "original_transition_fraction": float(smooth_row["transition_fraction"]),
        "relaxed_transition_fraction": float(params["transition_fraction"]),
        "transition_fraction_delta": float(params["transition_fraction"] - smooth_row["transition_fraction"]),
        "original_k": float(smooth_row["k"]),
        "relaxed_k": float(params["k"]),
        "original_w": float(smooth_row["w"]),
        "relaxed_w": float(params["w"]),
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "tried_count": len(tried),
    }
    row.update(flatten_params(params, prefix="relaxed_"))
    return row


def compute_aggregates(smooth_rows, relaxed_rows):
    w_values = [float(row["w"]) for row in smooth_rows]
    k_values = [float(row["k"]) for row in smooth_rows]
    frac_values = [float(row["transition_fraction"]) for row in smooth_rows]
    sharp_rows = [row["corpus"] for row in smooth_rows if float(row["w"]) < 0.2]
    rmse_deltas = [float(row["rmse_delta_vs_original"]) for row in relaxed_rows]
    frac_deltas = [float(row["transition_fraction_delta"]) for row in relaxed_rows]
    metrics = [
        ("english_corpus_count", len(smooth_rows), "int", "Canonical English corpus count for the smooth-fit experiment."),
        ("smooth_beats_single_zm_count", sum(bool(row["beats_single_zm"]) for row in smooth_rows), "int", "Count of corpora where the smooth two-regime fit beats single ZM RMSE."),
        ("smooth_beats_piecewise_count", sum(bool(row["beats_piecewise"]) for row in smooth_rows), "int", "Count of corpora where the smooth two-regime fit beats hard piecewise ZM."),
        ("smooth_beats_step2_count", sum(bool(row["beats_step2"]) for row in smooth_rows), "int", "Count of corpora where the smooth two-regime fit beats ZM plus step-2."),
        ("mean_k", statistics.fmean(k_values), "float6", "Mean fitted transition centre k."),
        ("sd_k", statistics.pstdev(k_values), "float6", "Population standard deviation of fitted transition centre k."),
        ("min_k", min(k_values), "float6", "Minimum fitted transition centre k."),
        ("max_k", max(k_values), "float6", "Maximum fitted transition centre k."),
        ("mean_w", statistics.fmean(w_values), "float6", "Mean fitted smooth width w in log-rank units."),
        ("sd_w", statistics.pstdev(w_values), "float6", "Population standard deviation of fitted smooth width w."),
        ("min_w", min(w_values), "float6", "Minimum fitted smooth width w."),
        ("max_w", max(w_values), "float6", "Maximum fitted smooth width w."),
        ("mean_transition_fraction", statistics.fmean(frac_values), "float6", "Mean fitted transition fraction log(k)/log(V)."),
        ("sd_transition_fraction", statistics.pstdev(frac_values), "float6", "Population standard deviation of fitted transition fraction."),
        ("min_transition_fraction", min(frac_values), "float6", "Minimum fitted transition fraction."),
        ("max_transition_fraction", max(frac_values), "float6", "Maximum fitted transition fraction."),
        ("w_lt_0_2_count", len(sharp_rows), "int", "Count of corpora with fitted width w < 0.2."),
        ("w_lt_0_2_corpora", "; ".join(sharp_rows), "string", "Semicolon-delimited corpus list with w < 0.2."),
        ("relaxed_bounds_improves_rmse_count", sum(delta < 0.0 for delta in rmse_deltas), "int", "Count of corpora where relaxed bounds improve RMSE over the constrained 8-parameter fit."),
        ("relaxed_bounds_unchanged_count", sum(abs(delta) <= 1e-12 for delta in rmse_deltas), "int", "Count of corpora unchanged to numerical tolerance under relaxed bounds."),
        ("relaxed_bounds_degraded_count", sum(delta > 0.0 for delta in rmse_deltas), "int", "Count of corpora where relaxed bounds worsen RMSE."),
        ("relaxed_rmse_delta_mean", statistics.fmean(rmse_deltas), "float9", "Mean RMSE delta relaxed minus constrained."),
        ("relaxed_rmse_delta_median", statistics.median(rmse_deltas), "float9", "Median RMSE delta relaxed minus constrained."),
        ("relaxed_rmse_delta_min", min(rmse_deltas), "float9", "Minimum RMSE delta relaxed minus constrained."),
        ("relaxed_rmse_delta_max", max(rmse_deltas), "float9", "Maximum RMSE delta relaxed minus constrained."),
        ("relaxed_rmse_delta_sd", statistics.pstdev(rmse_deltas), "float9", "Population standard deviation of RMSE deltas relaxed minus constrained."),
        ("relaxed_transition_fraction_delta_mean", statistics.fmean(frac_deltas), "float9", "Mean transition-fraction delta relaxed minus constrained."),
        ("relaxed_transition_fraction_delta_median", statistics.median(frac_deltas), "float9", "Median transition-fraction delta relaxed minus constrained."),
    ]
    return [{"metric_name": name, "value": value, "display_format": fmt, "notes": notes} for name, value, fmt, notes in metrics]


def compare_rows(new_rows, old_rows, fields):
    diff_rows = []
    max_abs_numeric_diff = 0.0
    max_field = "all compared fields matched exactly"
    mismatch_count = 0
    for slug, new_row in {row["slug"]: row for row in new_rows}.items():
        old_row = old_rows[slug]
        for field in fields:
            new_value = new_row[field]
            old_value = old_row[field]
            if isinstance(new_value, bool):
                matches = bool(new_value) == bool(old_value)
                mismatch_count += int(not matches)
                diff_rows.append({"slug": slug, "field": field, "new_value": new_value, "historical_value": old_value, "abs_diff": "", "matches": matches})
            elif isinstance(new_value, (int, float)):
                old_float = float(old_value)
                abs_diff = abs(float(new_value) - old_float)
                if abs_diff > max_abs_numeric_diff:
                    max_abs_numeric_diff = abs_diff
                    max_field = f"{slug}:{field}"
                diff_rows.append({"slug": slug, "field": field, "new_value": new_value, "historical_value": old_value, "abs_diff": abs_diff, "matches": abs_diff <= 1e-12})
            else:
                matches = str(new_value) == str(old_value)
                mismatch_count += int(not matches)
                diff_rows.append({"slug": slug, "field": field, "new_value": new_value, "historical_value": old_value, "abs_diff": "", "matches": matches})
    return diff_rows, max_abs_numeric_diff, max_field, mismatch_count


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

    smooth_rows = []
    for spec in cfg.SEARCHED_CORPORA:
        smooth_rows.append(build_smooth_row(spec))

    smooth_by_slug = {row["slug"]: row for row in smooth_rows}
    relaxed_rows = []
    for idx, spec in enumerate(cfg.SEARCHED_CORPORA, start=1):
        relaxed_rows.append(build_relaxed_row(spec, idx, smooth_by_slug[spec["slug"]]))

    aggregate_rows = compute_aggregates(smooth_rows, relaxed_rows)

    historical_smooth = normalize_historical_smooth_rows(historical_rows(cfg.HISTORICAL_CONSTRAINED_SUMMARY))
    historical_relaxed = normalize_historical_relaxed_rows(historical_rows(cfg.HISTORICAL_RELAXED_SUMMARY))

    smooth_compare_fields = [
        "token_count",
        "vocab_size",
        "single_zm_rmse",
        "piecewise_rmse",
        "step2_rmse",
        "step2_expression",
        "reranked_rmse",
        "beats_single_zm",
        "beats_piecewise",
        "beats_step2",
        "best_start_index",
        "best_nfev",
        "a1",
        "b1",
        "c1",
        "a2",
        "b2",
        "c2",
        "k",
        "w",
        "transition_fraction",
        "tried_count",
    ]
    relaxed_compare_fields = [
        "original_rmse",
        "relaxed_rmse",
        "rmse_delta_vs_original",
        "original_transition_fraction",
        "relaxed_transition_fraction",
        "transition_fraction_delta",
        "original_k",
        "relaxed_k",
        "original_w",
        "relaxed_w",
        "best_start_index",
        "best_nfev",
        "relaxed_a1",
        "relaxed_b1",
        "relaxed_c1",
        "relaxed_a2",
        "relaxed_b2",
        "relaxed_c2",
        "tried_count",
    ]
    smooth_diff_rows, smooth_max_diff, smooth_max_field, smooth_text_mismatches = compare_rows(smooth_rows, historical_smooth, smooth_compare_fields)
    relaxed_diff_rows, relaxed_max_diff, relaxed_max_field, relaxed_text_mismatches = compare_rows(relaxed_rows, historical_relaxed, relaxed_compare_fields)

    diff_summary = [
        {"metric_name": "smooth_fit_row_count", "value": len(smooth_rows), "notes": "Expected 25 constrained smooth-fit rows."},
        {"metric_name": "bounds_robustness_row_count", "value": len(relaxed_rows), "notes": "Expected 25 relaxed-bounds rows."},
        {"metric_name": "smooth_max_abs_numeric_diff", "value": smooth_max_diff, "notes": smooth_max_field},
        {"metric_name": "smooth_text_mismatch_count", "value": smooth_text_mismatches, "notes": "Count of non-numeric mismatches against historical constrained summary."},
        {"metric_name": "relaxed_max_abs_numeric_diff", "value": relaxed_max_diff, "notes": relaxed_max_field},
        {"metric_name": "relaxed_text_mismatch_count", "value": relaxed_text_mismatches, "notes": "Count of non-numeric mismatches against historical relaxed summary."},
    ]

    write_csv(output_dir / "smooth_fit_per_corpus.csv", SMOOTH_FIELDS, smooth_rows)
    write_csv(output_dir / "bounds_robustness_per_corpus.csv", RELAXED_FIELDS, relaxed_rows)
    write_csv(output_dir / "aggregate_statistics.csv", AGG_FIELDS, aggregate_rows)
    write_csv(output_dir / "historical_smooth_diff.csv", DIFF_FIELDS, smooth_diff_rows)
    write_csv(output_dir / "historical_relaxed_diff.csv", DIFF_FIELDS, relaxed_diff_rows)
    write_csv(output_dir / "historical_diff_summary.csv", SUMMARY_FIELDS, diff_summary)

    manifest = {
        "experiment_id": "3a",
        "output_dir": str(output_dir),
        "smooth_row_count": len(smooth_rows),
        "bounds_row_count": len(relaxed_rows),
        "historical_smooth_max_abs_numeric_diff": smooth_max_diff,
        "historical_relaxed_max_abs_numeric_diff": relaxed_max_diff,
        "historical_smooth_text_mismatch_count": smooth_text_mismatches,
        "historical_relaxed_text_mismatch_count": relaxed_text_mismatches,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote {output_dir / 'smooth_fit_per_corpus.csv'}")
    print(f"Wrote {output_dir / 'bounds_robustness_per_corpus.csv'}")
    print(f"Wrote {output_dir / 'aggregate_statistics.csv'}")
    print(f"Wrote {output_dir / 'historical_smooth_diff.csv'}")
    print(f"Wrote {output_dir / 'historical_relaxed_diff.csv'}")
    print(f"Wrote {output_dir / 'historical_diff_summary.csv'}")


if __name__ == "__main__":
    main()
