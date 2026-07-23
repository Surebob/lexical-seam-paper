from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
BASE_PATH = ROOT / "zipf_seam_mandelbrot_pmf.py"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_angle1_head_windows"
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"

CUTOFFS = [50, 100, 200, 500, 1000, None]
STEP2_VARIANCE_THRESHOLD = 1e-10


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_angle1_common")
base = load_module(BASE_PATH, "zipf_angle1_base")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_angle1_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def parse_args():
    p = argparse.ArgumentParser(description="Head-window held-out NLL audit using canonical split-fit PMFs")
    p.add_argument("--softk-summary", type=Path, default=SOFTK_SUMMARY_PATH)
    p.add_argument("--base-summary", type=Path, default=BASE_SUMMARY_PATH)
    p.add_argument("--outdir", type=Path, default=OUTDIR)
    return p.parse_args()


def cutoff_label(cutoff: int | None) -> str:
    return "full" if cutoff is None else f"top{cutoff}"


def load_maps(base_summary_path: Path, softk_summary_path: Path):
    base_rows = json.loads(base_summary_path.read_text(encoding="utf-8"))["rows"]
    soft_rows = json.loads(softk_summary_path.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in base_rows}, {row["slug"]: row for row in soft_rows}


def logpmf_from_params(model_name: str, params: dict, ranks: np.ndarray) -> np.ndarray | None:
    if model_name == "zipf":
        return base.zipf_logpmf(float(params["alpha"]), ranks)
    if model_name == "zm":
        return base.zm_logpmf(float(params["b"]), float(params["c"]), ranks)
    if model_name == "moe":
        return base.moe_logpmf(float(params["alpha"]), float(params["beta"]), ranks)
    if model_name == "softk":
        arr = np.array(
            [
                float(params["b1"]),
                float(params["c1"]),
                float(params["b2"]),
                float(params["c2"]),
                float(params["k"]),
                float(params["w"]),
            ],
            dtype=np.float64,
        )
        return base.seam_logpmf(arr, ranks)
    raise KeyError(model_name)


def restricted_avg_nll(test_counts: np.ndarray, logpmf: np.ndarray | None, cutoff: int | None) -> float:
    if logpmf is None:
        return float("inf")
    end = len(test_counts) if cutoff is None else min(int(cutoff), len(test_counts))
    counts = np.asarray(test_counts[:end], dtype=np.float64)
    denom = float(np.sum(counts))
    if denom <= 0.0:
        return float("inf")
    return float(-np.dot(counts, logpmf[:end]) / denom)


def exact_step2_window(log_rank: np.ndarray, y: np.ndarray, prediction_log: np.ndarray) -> dict:
    x = common.normalize_x(log_rank)
    residual = y - prediction_log
    vocab0 = base.eml.initial_vocabulary(x, residual)
    step1 = base.eml.generate_candidates(vocab0, residual, 1)
    step1 = base.eml.dedupe_candidates(step1)
    step1 = base.eml.filter_candidates(step1, STEP2_VARIANCE_THRESHOLD)
    vocab = vocab0 + step1
    step2 = base.eml.generate_candidates(vocab, residual, 2)
    step2 = base.eml.dedupe_candidates(step2)
    step2 = base.eml.filter_candidates(step2, STEP2_VARIANCE_THRESHOLD)
    winner = min(step2, key=lambda row: (row["rmse"], row["expr"]))
    composite_rmse = base.eml.rmse(y, prediction_log + winner["values"])
    baseline_rmse = common.rmse(y, prediction_log)
    return {
        "expr": winner["expr"],
        "gain": float(baseline_rmse - composite_rmse),
        "helps": bool(composite_rmse + 1e-12 < baseline_rmse),
    }


def analyze_corpus(spec: dict, base_row: dict, soft_row: dict, corpus_index: int) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=base.SPLIT_SEED + corpus_index, train_fraction=base.TRAIN_FRACTION)
    train_ranks = np.arange(1, len(split["train_counts"]) + 1, dtype=np.float64)

    model_params = {
        "zipf": base_row["models"]["zipf"]["params"],
        "zm": base_row["models"]["zm"]["params"],
        "moe": base_row["models"]["moe"]["params"],
        "softk": protocol.split_fit_params(soft_row),
    }
    logpmfs = {name: logpmf_from_params(name, params, train_ranks) for name, params in model_params.items()}

    heldout = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        vals = {name: restricted_avg_nll(split["test_counts"], logpmf, cutoff) for name, logpmf in logpmfs.items()}
        winner = min(vals, key=vals.get)
        heldout[label] = {
            "avg_nll": vals,
            "winner": winner,
            "softk_minus_moe": float(vals["softk"] - vals["moe"]),
            "softk_minus_zm": float(vals["softk"] - vals["zm"]),
            "softk_minus_zipf": float(vals["softk"] - vals["zipf"]),
        }

    full_ranks = dataset["ranks"].astype(np.float64)
    soft_full_logpmf = logpmf_from_params("softk", protocol.full_refit_params(soft_row), full_ranks)
    soft_prediction_log = math.log(float(dataset["token_count"])) + soft_full_logpmf
    step2 = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        end = len(full_ranks) if cutoff is None else min(int(cutoff), len(full_ranks))
        step2[label] = exact_step2_window(
            dataset["log_rank"][:end],
            dataset["log_freq"][:end],
            soft_prediction_log[:end],
        )

    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "heldout": heldout,
        "step2": step2,
    }


def summarize(rows: list[dict]) -> dict:
    cutoff_summary = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        winner_counts = {"zipf": 0, "zm": 0, "moe": 0, "softk": 0}
        soft_vs_moe = 0
        soft_vs_zm = 0
        soft_vs_zipf = 0
        deltas_moe = []
        deltas_zm = []
        step2_help = 0
        for row in rows:
            info = row["heldout"][label]
            winner_counts[info["winner"]] += 1
            soft_vs_moe += info["avg_nll"]["softk"] < info["avg_nll"]["moe"]
            soft_vs_zm += info["avg_nll"]["softk"] < info["avg_nll"]["zm"]
            soft_vs_zipf += info["avg_nll"]["softk"] < info["avg_nll"]["zipf"]
            deltas_moe.append(info["softk_minus_moe"])
            deltas_zm.append(info["softk_minus_zm"])
            step2_help += int(row["step2"][label]["helps"])
        cutoff_summary[label] = {
            "winner_counts": winner_counts,
            "softk_beats_moe": int(soft_vs_moe),
            "softk_beats_zm": int(soft_vs_zm),
            "softk_beats_zipf": int(soft_vs_zipf),
            "median_softk_minus_moe": float(np.median(deltas_moe)),
            "median_softk_minus_zm": float(np.median(deltas_zm)),
            "step2_help_count": int(step2_help),
        }
    return {"rows": rows, "cutoffs": cutoff_summary}


def write_csv(rows: list[dict], path: Path):
    fieldnames = ["slug", "name", "cutoff", "winner", "zipf", "zm", "moe", "softk", "softk_minus_moe", "softk_minus_zm", "step2_help", "step2_gain"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            for cutoff in CUTOFFS:
                label = cutoff_label(cutoff)
                info = row["heldout"][label]
                step = row["step2"][label]
                writer.writerow(
                    {
                        "slug": row["slug"],
                        "name": row["name"],
                        "cutoff": label,
                        "winner": info["winner"],
                        "zipf": info["avg_nll"]["zipf"],
                        "zm": info["avg_nll"]["zm"],
                        "moe": info["avg_nll"]["moe"],
                        "softk": info["avg_nll"]["softk"],
                        "softk_minus_moe": info["softk_minus_moe"],
                        "softk_minus_zm": info["softk_minus_zm"],
                        "step2_help": int(step["helps"]),
                        "step2_gain": step["gain"],
                    }
                )


def build_report(summary: dict) -> str:
    lines = [
        "# Angle 1: Head-Window Held-out NLL",
        "",
        "- Same 80/20 binomial type-count split as the Seam-Mandelbrot PMF runs.",
        "- Held-out average NLL is restricted to train-rank windows `1-K` without refitting.",
        "- Step-2 help counts are recomputed on the soft-k residual restricted to the same head window.",
        "",
    ]
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        info = summary["cutoffs"][label]
        lines.extend(
            [
                f"## {label}",
                "",
                f"- soft-k beats MOE: `{info['softk_beats_moe']}` / 25",
                f"- soft-k beats ZM: `{info['softk_beats_zm']}` / 25",
                f"- soft-k beats Zipf: `{info['softk_beats_zipf']}` / 25",
                f"- median soft-k minus MOE held-out avg NLL: `{info['median_softk_minus_moe']:.12f}`",
                f"- median soft-k minus ZM held-out avg NLL: `{info['median_softk_minus_zm']:.12f}`",
                f"- winner counts: `{json.dumps(info['winner_counts'])}`",
                f"- soft-k step-2 help count: `{info['step2_help_count']}` / 25",
                "",
            ]
        )
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    base_map, soft_map = load_maps(args.base_summary, args.softk_summary)
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(analyze_corpus(spec, base_map[spec["slug"]], soft_map[spec["slug"]], corpus_index))
    summary = summarize(rows)
    summary["metadata"] = {
        "softk_summary": str(args.softk_summary),
        "base_summary": str(args.base_summary),
        "canonical_heldout_protocol": "soft-k held-out windows are evaluated with split-fit params only; full-refit params are used only for full-corpus residual step-2 diagnostics",
    }
    (args.outdir / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (args.outdir / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, args.outdir / "head_window_table.csv")


if __name__ == "__main__":
    main()
