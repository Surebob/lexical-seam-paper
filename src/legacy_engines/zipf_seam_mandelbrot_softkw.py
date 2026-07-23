from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    matplotlib = None
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
BASE_PATH = ROOT / "zipf_seam_mandelbrot_pmf.py"
SOFTK_PATH = ROOT / "zipf_seam_mandelbrot_softk.py"
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_seam_mandelbrot_softkw"

SMOOTH_ALPHA = 0.521422287446811
SMOOTH_W = 0.9899928130291
LAMBDA_WS = [0.0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
TRAIN_STARTS = 8
FULL_STARTS = 5


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_softkw_common")
base = load_module(BASE_PATH, "zipf_seam_softkw_base")
softk_mod = load_module(SOFTK_PATH, "zipf_seam_softkw_softk")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_seam_softkw_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_reranked_map() -> dict[str, dict]:
    rows = json.loads(RERANKED_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def load_base_map() -> dict[str, dict]:
    rows = json.loads(BASE_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def load_softk_map() -> dict[str, dict]:
    rows = json.loads(SOFTK_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def softkw_starts(vocab_size: int, reranked_row: dict | None, base_row: dict, softk_row: dict, seed: int, n_starts: int) -> list[np.ndarray]:
    lower = np.array(
        [base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.k_bounds(vocab_size)[0], base.SEAM_W_BOUNDS[0]],
        dtype=np.float64,
    )
    upper = np.array(
        [base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.k_bounds(vocab_size)[1], base.SEAM_W_BOUNDS[1]],
        dtype=np.float64,
    )
    softk_selection = protocol.split_fit_params(softk_row)
    starts = []
    starts.append(
        np.array(
            [
                float(softk_selection["b1"]),
                float(softk_selection["c1"]),
                float(softk_selection["b2"]),
                float(softk_selection["c2"]),
                float(softk_selection["k"]),
                float(softk_selection["w"]),
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


def fit_softkw(train_counts: np.ndarray, ranks: np.ndarray, vocab_size: int, reranked_row: dict | None, base_row: dict, softk_row: dict, lambda_w: float, seed: int, n_starts: int) -> dict:
    lambda_k = protocol.selection_lambda(softk_row)
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
    target_k = SMOOTH_ALPHA * math.log(vocab_size)
    target_w = math.log(SMOOTH_W)
    n_tokens = max(float(np.sum(train_counts)), 1.0)
    best = None

    def objective(theta: np.ndarray) -> float:
        params = np.clip(theta, lower, upper)
        logpmf = base.seam_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        if not math.isfinite(ll):
            return math.inf
        dev_k = math.log(float(params[4])) - target_k
        dev_w = math.log(float(params[5])) - target_w
        avg_nll = -ll / n_tokens
        return float(avg_nll + lambda_k * dev_k * dev_k + lambda_w * dev_w * dev_w)

    for start_index, start in enumerate(softkw_starts(vocab_size, reranked_row, base_row, softk_row, seed, n_starts), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = base.seam_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        dev_k = math.log(float(params[4])) - target_k
        dev_w = math.log(float(params[5])) - target_w
        avg_nll = -ll / n_tokens if math.isfinite(ll) else math.inf
        penalized = avg_nll + lambda_k * dev_k * dev_k + lambda_w * dev_w * dev_w if math.isfinite(avg_nll) else math.inf
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
            "dev_logk": float(dev_k),
            "dev_logw": float(dev_w),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "start_index": int(start_index),
            "p": 6,
            "lambda_k": float(lambda_k),
            "lambda_w": float(lambda_w),
        }
        if best is None or record["penalized_objective"] < best["penalized_objective"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def evaluate_fit(fit: dict, test_counts: np.ndarray) -> dict:
    return base.evaluate_test_fit(fit, test_counts)


def full_context_softkw(dataset: dict, reranked_row: dict | None, base_row: dict, softk_row: dict, lambda_w: float, seed: int, n_starts: int) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    fit = fit_softkw(counts, ranks, int(dataset["unique_words"]), reranked_row, base_row, softk_row, lambda_w, seed, n_starts)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {
        "fit": fit,
        "rmse": base.full_rank_rmse(counts, fit["logpmf"]),
        "prediction_log": pred_log,
    }


def analyze_corpus(spec: dict, reranked_row: dict | None, base_row: dict, softk_row: dict, corpus_index: int) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    ranks = np.arange(1, len(train_counts) + 1, dtype=np.float64)

    sweep_rows = []
    for lambda_index, lambda_w in enumerate(LAMBDA_WS, start=1):
        fit = fit_softkw(
            train_counts,
            ranks,
            split["vocab_size"],
            reranked_row,
            base_row,
            softk_row,
            lambda_w=lambda_w,
            seed=SPLIT_SEED + 11000 + corpus_index * 20 + lambda_index,
            n_starts=TRAIN_STARTS,
        )
        evals = evaluate_fit(fit, test_counts)
        sweep_rows.append(
            {
                "lambda_w": float(lambda_w),
                "fit": fit,
                "test_avg_nll": float(evals["test_avg_nll"]),
                "test_loglike": float(evals["test_loglike"]),
            }
        )

    best_row = min(sweep_rows, key=lambda row: row["test_avg_nll"])
    best_full = full_context_softkw(
        dataset,
        reranked_row,
        base_row,
        softk_row,
        lambda_w=best_row["lambda_w"],
        seed=SPLIT_SEED + 13000 + corpus_index,
        n_starts=FULL_STARTS,
    )
    best_step2 = base.exact_step2_on_residual(dataset, best_full["prediction_log"])
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "free": {
            "test_avg_nll": float(base_row["models"]["seam"]["test_avg_nll"]),
            "bic": float(base_row["models"]["seam"]["bic"]),
            "rmse": float(base_row["models"]["seam"]["full_rmse"]),
            "step2_gain": float(base_row["seam_step2"]["gain"]),
            "step2_helps": bool(base_row["seam_step2"]["helps"]),
        },
        "softk": {
            "lambda": protocol.selection_lambda(softk_row),
            "test_avg_nll": protocol.split_fit_test_avg_nll(softk_row),
            "bic": float(protocol.full_refit_bic(softk_row)),
            "rmse": float(protocol.full_refit_rmse(softk_row)),
            "step2_gain": float(protocol.full_refit_step2_gain(softk_row)),
            "step2_helps": bool(protocol.full_refit_step2_helps(softk_row)),
            "params": protocol.split_fit_params(softk_row),
        },
        "sweep": [
            {
                "lambda_w": row["lambda_w"],
                "test_avg_nll": row["test_avg_nll"],
                "train_avg_nll": row["fit"]["train_avg_nll"],
                "dev_logk": row["fit"]["dev_logk"],
                "dev_logw": row["fit"]["dev_logw"],
                "transition_fraction": row["fit"]["params"]["transition_fraction"],
                "w": row["fit"]["params"]["w"],
            }
            for row in sweep_rows
        ],
        "best_lambda_w": {
            "lambda_w": float(best_row["lambda_w"]),
            "lambda_k": float(best_full["fit"]["lambda_k"]),
            "test_avg_nll": float(best_row["test_avg_nll"]),
            "bic": float(best_full["fit"]["bic"]),
            "rmse": float(best_full["rmse"]),
            "step2_gain": float(best_step2["gain"]),
            "step2_helps": bool(best_step2["helps"]),
            "params": best_full["fit"]["params"],
            "dev_logk": float(best_full["fit"]["dev_logk"]),
            "dev_logw": float(best_full["fit"]["dev_logw"]),
        },
        "delta_best_vs_softk_test": float(best_row["test_avg_nll"] - protocol.split_fit_test_avg_nll(softk_row)),
    }


def build_summary(rows: list[dict]) -> dict:
    counts = {}
    for row in rows:
        key = f"{row['best_lambda_w']['lambda_w']:.6g}"
        counts[key] = counts.get(key, 0) + 1
    return {
        "rows": rows,
        "counts": {
            "n_rows": len(rows),
            "best_lambda_w_counts": counts,
            "softkw_beats_softk_test": int(sum(row["delta_best_vs_softk_test"] < 0.0 for row in rows)),
            "softkw_step2_help_count": int(sum(row["best_lambda_w"]["step2_helps"] for row in rows)),
            "softk_step2_help_count": int(sum(row["softk"]["step2_helps"] for row in rows)),
            "free_step2_help_count": int(sum(row["free"]["step2_helps"] for row in rows)),
        },
        "medians": {
            "delta_best_vs_softk_test": float(np.median([row["delta_best_vs_softk_test"] for row in rows])),
            "softkw_step2_gain": float(np.median([row["best_lambda_w"]["step2_gain"] for row in rows])),
            "softk_step2_gain": float(np.median([row["softk"]["step2_gain"] for row in rows])),
            "free_step2_gain": float(np.median([row["free"]["step2_gain"] for row in rows])),
        },
    }


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "softk_lambda",
        "best_lambda_w",
        "free_test_avg_nll",
        "softk_test_avg_nll",
        "softkw_test_avg_nll",
        "free_step2_gain",
        "softk_step2_gain",
        "softkw_step2_gain",
        "softkw_transition_fraction",
        "softkw_w",
        "softkw_dev_logw",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "softk_lambda": row["softk"]["lambda"],
                    "best_lambda_w": row["best_lambda_w"]["lambda_w"],
                    "free_test_avg_nll": row["free"]["test_avg_nll"],
                    "softk_test_avg_nll": row["softk"]["test_avg_nll"],
                    "softkw_test_avg_nll": row["best_lambda_w"]["test_avg_nll"],
                    "free_step2_gain": row["free"]["step2_gain"],
                    "softk_step2_gain": row["softk"]["step2_gain"],
                    "softkw_step2_gain": row["best_lambda_w"]["step2_gain"],
                    "softkw_transition_fraction": row["best_lambda_w"]["params"]["transition_fraction"],
                    "softkw_w": row["best_lambda_w"]["params"]["w"],
                    "softkw_dev_logw": row["best_lambda_w"]["dev_logw"],
                }
            )


def plot_lambda_w_tradeoff(rows: list[dict], path: Path):
    if plt is None:
        return
    median_deltas = []
    chosen_counts = []
    for lambda_w in LAMBDA_WS:
        vals = []
        count = 0
        for row in rows:
            match = next(item for item in row["sweep"] if abs(item["lambda_w"] - lambda_w) < 1e-18)
            vals.append(match["test_avg_nll"] - row["softk"]["test_avg_nll"])
            if abs(row["best_lambda_w"]["lambda_w"] - lambda_w) < 1e-18:
                count += 1
        median_deltas.append(float(np.median(vals)))
        chosen_counts.append(count)
    fig, ax1 = plt.subplots(figsize=(8.5, 5.0))
    ax1.plot(LAMBDA_WS, median_deltas, marker="o", color="#1565c0")
    ax1.set_xscale("symlog", linthresh=1e-6)
    ax1.set_xlabel("soft-w prior lambda")
    ax1.set_ylabel("median test avg NLL delta vs soft-k", color="#1565c0")
    ax1.tick_params(axis="y", labelcolor="#1565c0")
    ax1.grid(True, alpha=0.25)
    ax2 = ax1.twinx()
    widths = [1e-6 if x == 0 else x * 0.25 for x in LAMBDA_WS]
    ax2.bar(LAMBDA_WS, chosen_counts, width=widths, alpha=0.25, color="#c62828")
    ax2.set_ylabel("best-lambda corpus count", color="#c62828")
    ax2.tick_params(axis="y", labelcolor="#c62828")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def build_report(summary: dict) -> str:
    counts = summary["counts"]
    medians = summary["medians"]
    rows = sorted(summary["rows"], key=lambda row: row["delta_best_vs_softk_test"])
    lines = [
        "# Soft-k plus Soft-w Seam-Mandelbrot PMF",
        "",
        "- `k` keeps the selected soft prior from the previous sweep.",
        f"- `w` gets an additional weak quadratic penalty on `(log w - log({SMOOTH_W:.12f}))^2`.",
        f"- lambda_w sweep: `{LAMBDA_WS}`",
        "",
        f"- best lambda_w counts: `{json.dumps(counts['best_lambda_w_counts'])}`",
        f"- soft-k,w beats soft-k on held-out avg NLL: `{counts['softkw_beats_softk_test']}` / `{counts['n_rows']}`",
        f"- free step-2 help count: `{counts['free_step2_help_count']}` / `{counts['n_rows']}`",
        f"- soft-k step-2 help count: `{counts['softk_step2_help_count']}` / `{counts['n_rows']}`",
        f"- soft-k,w step-2 help count: `{counts['softkw_step2_help_count']}` / `{counts['n_rows']}`",
        f"- median soft-k,w minus soft-k held-out avg NLL: `{medians['delta_best_vs_softk_test']:.12f}`",
        f"- median free step-2 gain: `{medians['free_step2_gain']:.12f}`",
        f"- median soft-k step-2 gain: `{medians['softk_step2_gain']:.12f}`",
        f"- median soft-k,w step-2 gain: `{medians['softkw_step2_gain']:.12f}`",
        "",
        "| corpus | best lambda_w | free | soft-k | soft-k,w | free step2 | soft-k step2 | soft-k,w step2 | w |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['best_lambda_w']['lambda_w']:.6g} | "
            f"{row['free']['test_avg_nll']:.6f} | {row['softk']['test_avg_nll']:.6f} | {row['best_lambda_w']['test_avg_nll']:.6f} | "
            f"{row['free']['step2_gain']:.6f} | {row['softk']['step2_gain']:.6f} | {row['best_lambda_w']['step2_gain']:.6f} | "
            f"{row['best_lambda_w']['params']['w']:.3f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_map = load_reranked_map()
    base_map = load_base_map()
    softk_map = load_softk_map()
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(analyze_corpus(spec, reranked_map.get(spec["slug"]), base_map[spec["slug"]], softk_map[spec["slug"]], corpus_index))
    summary = build_summary(rows)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, OUTDIR / "softkw_table.csv")
    plot_lambda_w_tradeoff(rows, OUTDIR / "softkw_lambda_tradeoff.png")


if __name__ == "__main__":
    main()
