from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
BASE_PATH = ROOT / "zipf_seam_mandelbrot_pmf.py"
ANGLE1_PATH = ROOT / "zipf_angle1_head_windows.py"
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_angle4_hierk"

CUTOFFS = [50, 100, 200, 500, 1000, None]
MAX_ITER = 4
TRAIN_STARTS = 6
FULL_STARTS = 4
SIGMA_MIN = 0.05
STEP2_VARIANCE_THRESHOLD = 1e-10


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_angle4_common")
base = load_module(BASE_PATH, "zipf_angle4_base")
angle1 = load_module(ANGLE1_PATH, "zipf_angle4_angle1")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_angle4_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def cutoff_label(cutoff: int | None) -> str:
    return "full" if cutoff is None else f"top{cutoff}"


def load_maps():
    base_rows = json.loads(BASE_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    soft_rows = json.loads(SOFTK_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    reranked_rows = json.loads(RERANKED_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return (
        {row["slug"]: row for row in base_rows},
        {row["slug"]: row for row in soft_rows},
        {row["slug"]: row for row in reranked_rows},
    )


def initial_alpha_sigma(soft_map: dict[str, dict]) -> tuple[float, float]:
    items = list(soft_map.values())
    logv = np.array([math.log(float(row["vocab_size"])) for row in items], dtype=np.float64)
    logk = np.array([math.log(float(protocol.split_fit_params(row)["k"])) for row in items], dtype=np.float64)
    alpha = float(np.dot(logk, logv) / np.dot(logv, logv))
    resid = logk - alpha * logv
    sigma = max(float(np.sqrt(np.mean(resid * resid))), SIGMA_MIN)
    return alpha, sigma


def hier_starts(vocab_size: int, base_row: dict, soft_row: dict, reranked_row: dict | None, prev_params: list[float] | None, seed: int, n_starts: int) -> list[np.ndarray]:
    lower = np.array(
        [base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.k_bounds(vocab_size)[0], base.SEAM_W_BOUNDS[0]],
        dtype=np.float64,
    )
    upper = np.array(
        [base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.k_bounds(vocab_size)[1], base.SEAM_W_BOUNDS[1]],
        dtype=np.float64,
    )
    starts = []
    if prev_params is not None:
        starts.append(np.asarray(prev_params, dtype=np.float64))
    soft_selection = protocol.split_fit_params(soft_row)
    starts.append(
        np.array(
            [
                float(soft_selection["b1"]),
                float(soft_selection["c1"]),
                float(soft_selection["b2"]),
                float(soft_selection["c2"]),
                float(soft_selection["k"]),
                float(soft_selection["w"]),
            ],
            dtype=np.float64,
        )
    )
    starts.append(
        np.array(
            [
                float(base_row["models"]["seam"]["params"]["b1"]),
                float(base_row["models"]["seam"]["params"]["c1"]),
                float(base_row["models"]["seam"]["params"]["b2"]),
                float(base_row["models"]["seam"]["params"]["c2"]),
                float(base_row["models"]["seam"]["params"]["k"]),
                float(base_row["models"]["seam"]["params"]["w"]),
            ],
            dtype=np.float64,
        )
    )
    if reranked_row is not None:
        params = reranked_row["params"]
        starts.append(
            np.array(
                [
                    float(params["b1"]),
                    float(params["c1"]),
                    float(params["b2"]),
                    float(params["c2"]),
                    float(min(max(params["k"], lower[4]), upper[4])),
                    float(params["w"]),
                ],
                dtype=np.float64,
            )
        )
    rng = np.random.default_rng(seed)
    while len(starts) < n_starts:
        starts.append(
            np.array(
                [
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    float(np.exp(rng.uniform(math.log(lower[4]), math.log(upper[4])))),
                    float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                ],
                dtype=np.float64,
            )
        )
    return [np.clip(start, lower, upper) for start in starts[:n_starts]]


def fit_hier_corpus(train_counts: np.ndarray, ranks: np.ndarray, vocab_size: int, alpha: float, sigma: float, base_row: dict, soft_row: dict, reranked_row: dict | None, seed: int, n_starts: int, prev_params: list[float] | None = None) -> dict:
    bounds = [
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.k_bounds(vocab_size),
        base.SEAM_W_BOUNDS,
    ]
    lower = np.array([b[0] for b in bounds], dtype=np.float64)
    upper = np.array([b[1] for b in bounds], dtype=np.float64)
    target = alpha * math.log(vocab_size)
    sigma2 = max(float(sigma) ** 2, SIGMA_MIN ** 2)
    n_tokens = max(float(np.sum(train_counts)), 1.0)
    best = None

    def objective(theta: np.ndarray) -> float:
        params = np.clip(theta, lower, upper)
        logpmf = base.seam_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        if not math.isfinite(ll):
            return math.inf
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens
        return float(avg_nll + (dev * dev) / (2.0 * sigma2))

    for start_index, start in enumerate(hier_starts(vocab_size, base_row, soft_row, reranked_row, prev_params, seed, n_starts), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = base.seam_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens if math.isfinite(ll) else math.inf
        penalized = avg_nll + (dev * dev) / (2.0 * sigma2) if math.isfinite(avg_nll) else math.inf
        record = {
            "params": {
                "b1": float(params[0]),
                "c1": float(params[1]),
                "b2": float(params[2]),
                "c2": float(params[3]),
                "k": float(params[4]),
                "w": float(params[5]),
                "transition_fraction": float(np.log(params[4]) / np.log(vocab_size)),
            },
            "params_vec": [float(v) for v in params],
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "train_avg_nll": float(avg_nll),
            "penalized_objective": float(penalized),
            "dev_logk": float(dev),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "start_index": int(start_index),
            "p": 6,
        }
        if best is None or record["penalized_objective"] < best["penalized_objective"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def update_alpha_sigma(fits: list[dict], vocab_sizes: list[int]) -> tuple[float, float]:
    logv = np.array([math.log(v) for v in vocab_sizes], dtype=np.float64)
    logk = np.array([math.log(float(fit["params"]["k"])) for fit in fits], dtype=np.float64)
    alpha = float(np.dot(logk, logv) / np.dot(logv, logv))
    resid = logk - alpha * logv
    sigma = max(float(np.sqrt(np.mean(resid * resid))), SIGMA_MIN)
    return alpha, sigma


def full_context(dataset: dict, alpha: float, sigma: float, base_row: dict, soft_row: dict, reranked_row: dict | None, prev_params: list[float] | None, seed: int) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    fit = fit_hier_corpus(counts, ranks, int(dataset["unique_words"]), alpha, sigma, base_row, soft_row, reranked_row, seed, FULL_STARTS, prev_params=prev_params)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {
        "fit": fit,
        "rmse": base.full_rank_rmse(counts, fit["logpmf"]),
        "prediction_log": pred_log,
    }


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


def restricted_avg_nll(test_counts: np.ndarray, logpmf: np.ndarray | None, cutoff: int | None) -> float:
    if logpmf is None:
        return float("inf")
    end = len(test_counts) if cutoff is None else min(int(cutoff), len(test_counts))
    counts = np.asarray(test_counts[:end], dtype=np.float64)
    denom = float(np.sum(counts))
    if denom <= 0.0:
        return float("inf")
    return float(-np.dot(counts, logpmf[:end]) / denom)


def analyze_all():
    base_map, soft_map, reranked_map = load_maps()
    datasets = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        dataset = common.build_zipf_dataset(common.corpus_path(spec))
        split = base.split_counts_by_train_rank(dataset, seed=base.SPLIT_SEED + corpus_index, train_fraction=base.TRAIN_FRACTION)
        datasets.append(
            {
                "corpus_index": corpus_index,
                "spec": spec,
                "dataset": dataset,
                "split": split,
                "base_row": base_map[spec["slug"]],
                "soft_row": soft_map[spec["slug"]],
                "reranked_row": reranked_map.get(spec["slug"]),
            }
        )

    alpha, sigma = initial_alpha_sigma(soft_map)
    prev_params = {item["spec"]["slug"]: None for item in datasets}
    fit_rows = None
    history = []
    for iteration in range(1, MAX_ITER + 1):
        fit_rows = []
        for item in datasets:
            ranks = np.arange(1, len(item["split"]["train_counts"]) + 1, dtype=np.float64)
            fit = fit_hier_corpus(
                item["split"]["train_counts"],
                ranks,
                item["split"]["vocab_size"],
                alpha,
                sigma,
                item["base_row"],
                item["soft_row"],
                item["reranked_row"],
                seed=base.SPLIT_SEED + 15000 + item["corpus_index"] * 10 + iteration,
                n_starts=TRAIN_STARTS,
                prev_params=prev_params[item["spec"]["slug"]],
            )
            prev_params[item["spec"]["slug"]] = fit["params_vec"]
            fit_rows.append(fit)
        alpha, sigma = update_alpha_sigma(fit_rows, [item["split"]["vocab_size"] for item in datasets])
        history.append({"iteration": iteration, "alpha": alpha, "sigma": sigma})

    rows = []
    for item, train_fit in zip(datasets, fit_rows):
        full = full_context(
            item["dataset"],
            alpha,
            sigma,
            item["base_row"],
            item["soft_row"],
            item["reranked_row"],
            prev_params[item["spec"]["slug"]],
            seed=base.SPLIT_SEED + 18000 + item["corpus_index"],
        )
        train_ranks = np.arange(1, len(item["split"]["train_counts"]) + 1, dtype=np.float64)
        logpmf_train = base.seam_logpmf(
            np.array(
                [
                    full["fit"]["params"]["b1"],
                    full["fit"]["params"]["c1"],
                    full["fit"]["params"]["b2"],
                    full["fit"]["params"]["c2"],
                    full["fit"]["params"]["k"],
                    full["fit"]["params"]["w"],
                ],
                dtype=np.float64,
            ),
            train_ranks,
        )
        heldout = {}
        for cutoff in CUTOFFS:
            label = cutoff_label(cutoff)
            vals = {
                "zipf": angle1.restricted_avg_nll(item["split"]["test_counts"], angle1.logpmf_from_params("zipf", item["base_row"]["models"]["zipf"]["params"], train_ranks), cutoff),
                "zm": angle1.restricted_avg_nll(item["split"]["test_counts"], angle1.logpmf_from_params("zm", item["base_row"]["models"]["zm"]["params"], train_ranks), cutoff),
                "moe": angle1.restricted_avg_nll(item["split"]["test_counts"], angle1.logpmf_from_params("moe", item["base_row"]["models"]["moe"]["params"], train_ranks), cutoff),
                "hier": angle1.restricted_avg_nll(item["split"]["test_counts"], logpmf_train, cutoff),
            }
            heldout[label] = {
                "avg_nll": vals,
                "winner": min(vals, key=vals.get),
                "hier_minus_moe": float(vals["hier"] - vals["moe"]),
                "hier_minus_zm": float(vals["hier"] - vals["zm"]),
                "hier_minus_zipf": float(vals["hier"] - vals["zipf"]),
            }

        step2 = {}
        for cutoff in CUTOFFS:
            label = cutoff_label(cutoff)
            end = len(item["dataset"]["ranks"]) if cutoff is None else min(int(cutoff), len(item["dataset"]["ranks"]))
            step2[label] = exact_step2_window(
                item["dataset"]["log_rank"][:end],
                item["dataset"]["log_freq"][:end],
                full["prediction_log"][:end],
            )

        rows.append(
            {
                "slug": item["spec"]["slug"],
                "name": item["spec"]["name"],
                "token_count": int(item["dataset"]["token_count"]),
                "vocab_size": int(item["dataset"]["unique_words"]),
                "params": full["fit"]["params"],
                "rmse": float(full["rmse"]),
                "heldout": heldout,
                "step2": step2,
            }
        )
    return rows, history, alpha, sigma


def summarize(rows: list[dict], history: list[dict], alpha: float, sigma: float) -> dict:
    cutoff_summary = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        winner_counts = {"zipf": 0, "zm": 0, "moe": 0, "hier": 0}
        hier_vs_moe = 0
        hier_vs_zm = 0
        hier_vs_zipf = 0
        deltas_moe = []
        deltas_zm = []
        step2_help = 0
        for row in rows:
            info = row["heldout"][label]
            winner_counts[info["winner"]] += 1
            hier_vs_moe += info["avg_nll"]["hier"] < info["avg_nll"]["moe"]
            hier_vs_zm += info["avg_nll"]["hier"] < info["avg_nll"]["zm"]
            hier_vs_zipf += info["avg_nll"]["hier"] < info["avg_nll"]["zipf"]
            deltas_moe.append(info["hier_minus_moe"])
            deltas_zm.append(info["hier_minus_zm"])
            step2_help += int(row["step2"][label]["helps"])
        cutoff_summary[label] = {
            "winner_counts": winner_counts,
            "hier_beats_moe": int(hier_vs_moe),
            "hier_beats_zm": int(hier_vs_zm),
            "hier_beats_zipf": int(hier_vs_zipf),
            "median_hier_minus_moe": float(np.median(deltas_moe)),
            "median_hier_minus_zm": float(np.median(deltas_zm)),
            "step2_help_count": int(step2_help),
        }
    return {"rows": rows, "history": history, "alpha": alpha, "sigma": sigma, "cutoffs": cutoff_summary}


def build_report(summary: dict) -> str:
    lines = [
        "# Angle 4: Empirical-Bayes Hierarchical k Pooling",
        "",
        "- Each corpus gets its own Seam-Mandelbrot PMF, but `log k_i` is pooled through a shared Gaussian prior.",
        "- The prior center and spread are fit iteratively across corpora.",
        f"- final alpha: `{summary['alpha']:.12f}`",
        f"- final sigma: `{summary['sigma']:.12f}`",
        "",
    ]
    for item in summary["history"]:
        lines.append(f"- iteration {item['iteration']}: alpha `{item['alpha']:.12f}`, sigma `{item['sigma']:.12f}`")
    lines.append("")
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        info = summary["cutoffs"][label]
        lines.extend(
            [
                f"## {label}",
                "",
                f"- hierarchical seam beats MOE: `{info['hier_beats_moe']}` / 25",
                f"- hierarchical seam beats ZM: `{info['hier_beats_zm']}` / 25",
                f"- hierarchical seam beats Zipf: `{info['hier_beats_zipf']}` / 25",
                f"- median hier minus MOE held-out avg NLL: `{info['median_hier_minus_moe']:.12f}`",
                f"- median hier minus ZM held-out avg NLL: `{info['median_hier_minus_zm']:.12f}`",
                f"- winner counts: `{json.dumps(info['winner_counts'])}`",
                f"- hierarchical seam step-2 help count: `{info['step2_help_count']}` / 25",
                "",
            ]
        )
    return "\n".join(lines) + "\n"


def write_csv(rows: list[dict], path: Path):
    fieldnames = ["slug", "name", "cutoff", "winner", "zipf", "zm", "moe", "hier", "hier_minus_moe", "hier_minus_zm", "step2_help", "step2_gain"]
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
                        "hier": info["avg_nll"]["hier"],
                        "hier_minus_moe": info["hier_minus_moe"],
                        "hier_minus_zm": info["hier_minus_zm"],
                        "step2_help": int(step["helps"]),
                        "step2_gain": step["gain"],
                    }
                )


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, history, alpha, sigma = analyze_all()
    summary = summarize(rows, history, alpha, sigma)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, OUTDIR / "hierk_head_window_table.csv")


if __name__ == "__main__":
    main()
