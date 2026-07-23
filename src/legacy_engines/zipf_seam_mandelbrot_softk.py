from __future__ import annotations

import argparse
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
REG_PATH = ROOT / "zipf_seam_mandelbrot_regularized.py"
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
REG_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_regularized" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_seam_mandelbrot_softk"

SMOOTH_ALPHA = 0.521422287446811
LAMBDAS = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]
TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
TRAIN_STARTS = 10
FULL_STARTS = 6


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_softk_common")
base = load_module(BASE_PATH, "zipf_seam_softk_base")
reg = load_module(REG_PATH, "zipf_seam_softk_reg")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def parse_args():
    p = argparse.ArgumentParser(description="Soft-k Seam-Mandelbrot PMF with explicit split-fit/full-refit records")
    p.add_argument("--outdir", type=Path, default=OUTDIR)
    p.add_argument(
        "--reuse-best-lambdas-summary",
        type=Path,
        default=None,
        help="Optional legacy/new soft-k summary.json whose chosen best lambda per corpus should be reused instead of resweeping.",
    )
    return p.parse_args()


def load_reranked_map() -> dict[str, dict]:
    rows = json.loads(RERANKED_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def load_base_map() -> dict[str, dict]:
    rows = json.loads(BASE_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def load_reg_map() -> dict[str, dict]:
    rows = json.loads(REG_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def load_reused_lambda_map(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    out = {}
    for row in rows:
        best = row["best_lambda"]
        if isinstance(best, dict) and "selection" in best:
            out[row["slug"]] = float(best["lambda"])
        else:
            out[row["slug"]] = float(best["lambda"])
    return out


def k_target(vocab_size: int) -> float:
    return float(vocab_size ** SMOOTH_ALPHA)


def softk_starts(vocab_size: int, reranked_row: dict | None, base_row: dict, lambda_prev_params: list[float] | None, seed: int, n_starts: int) -> list[np.ndarray]:
    lower = np.array(
        [base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.k_bounds(vocab_size)[0], base.SEAM_W_BOUNDS[0]],
        dtype=np.float64,
    )
    upper = np.array(
        [base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.k_bounds(vocab_size)[1], base.SEAM_W_BOUNDS[1]],
        dtype=np.float64,
    )
    starts = []
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
    if lambda_prev_params is not None:
        starts.append(np.asarray(lambda_prev_params, dtype=np.float64))
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
    starts.append(np.array([1.2, 10.0, 1.8, 10.0, float(min(max(k_target(vocab_size), lower[4]), upper[4])), 1.2], dtype=np.float64))
    starts.append(np.array([1.0, 1.0, 1.4, 100.0, float(min(max(k_target(vocab_size), lower[4]), upper[4])), 0.8], dtype=np.float64))
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


def fit_softk(train_counts: np.ndarray, ranks: np.ndarray, vocab_size: int, reranked_row: dict | None, base_row: dict, lam: float, seed: int, n_starts: int, lambda_prev_params: list[float] | None = None) -> dict:
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
    log_v = math.log(vocab_size)
    target = SMOOTH_ALPHA * log_v
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
        return float(avg_nll + lam * dev * dev)

    for start_index, start in enumerate(softk_starts(vocab_size, reranked_row, base_row, lambda_prev_params, seed, n_starts), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = base.seam_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens if math.isfinite(ll) else math.inf
        penalized = avg_nll + lam * dev * dev if math.isfinite(avg_nll) else math.inf
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
            "lambda": float(lam),
        }
        if best is None or record["penalized_objective"] < best["penalized_objective"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def evaluate_fit(fit: dict, test_counts: np.ndarray) -> dict:
    return base.evaluate_test_fit(fit, test_counts)


def full_context_softk(dataset: dict, reranked_row: dict | None, base_row: dict, lam: float, seed: int, n_starts: int, lambda_prev_params: list[float] | None = None) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    fit = fit_softk(counts, ranks, int(dataset["unique_words"]), reranked_row, base_row, lam, seed=seed, n_starts=n_starts, lambda_prev_params=lambda_prev_params)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {
        "fit": fit,
        "rmse": base.full_rank_rmse(counts, fit["logpmf"]),
        "prediction_log": pred_log,
    }


def analyze_corpus(spec: dict, reranked_row: dict | None, base_row: dict, corpus_index: int, forced_lambda: float | None = None) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    ranks = np.arange(1, len(train_counts) + 1, dtype=np.float64)
    vocab_size = split["vocab_size"]

    if forced_lambda is None:
        lambda_schedule = list(enumerate(LAMBDAS, start=1))
    else:
        lambda_index = next(index for index, lam in enumerate(LAMBDAS, start=1) if abs(lam - forced_lambda) < 1e-18)
        lambda_schedule = [(lambda_index, forced_lambda)]
    lambda_rows = []
    prev_params = None
    for lambda_index, lam in lambda_schedule:
        fit = fit_softk(
            train_counts,
            ranks,
            vocab_size,
            reranked_row,
            base_row,
            lam,
            seed=SPLIT_SEED + 7000 + corpus_index * 50 + lambda_index,
            n_starts=TRAIN_STARTS,
            lambda_prev_params=prev_params,
        )
        prev_params = fit["params_vec"]
        evals = evaluate_fit(fit, test_counts)
        lambda_rows.append(
            {
                "lambda": float(lam),
                "fit": fit,
                "test_loglike": float(evals["test_loglike"]),
                "test_avg_nll": float(evals["test_avg_nll"]),
            }
        )

    best_lambda_row = min(lambda_rows, key=lambda row: row["test_avg_nll"])
    best_lambda = float(best_lambda_row["lambda"])
    best_full = full_context_softk(
        dataset,
        reranked_row,
        base_row,
        lam=best_lambda,
        seed=SPLIT_SEED + 9000 + corpus_index,
        n_starts=FULL_STARTS,
        lambda_prev_params=best_lambda_row["fit"]["params_vec"],
    )
    best_step2 = base.exact_step2_on_residual(dataset, best_full["prediction_log"])

    sweep_summary = [
        {
            "lambda": row["lambda"],
            "test_avg_nll": row["test_avg_nll"],
            "train_avg_nll": row["fit"]["train_avg_nll"],
            "dev_logk": row["fit"]["dev_logk"],
            "transition_fraction": row["fit"]["params"]["transition_fraction"],
            "w": row["fit"]["params"]["w"],
        }
        for row in lambda_rows
    ]
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
            "params": base_row["models"]["seam"]["params"],
        },
        "fixed_k": {
            "test_avg_nll": float(reg.load_reg_map()[spec["slug"]]["variants"]["fixed_k"]["test_avg_nll"]) if False else None,
        },
        "sweep": sweep_summary,
        "best_lambda": {
            "lambda": best_lambda,
            "selection": {
                "test_avg_nll": float(best_lambda_row["test_avg_nll"]),
                "test_loglike": float(best_lambda_row["test_loglike"]),
                "train_avg_nll": float(best_lambda_row["fit"]["train_avg_nll"]),
                "train_loglike": float(best_lambda_row["fit"]["train_loglike"]),
                "params": best_lambda_row["fit"]["params"],
                "params_vec": [float(v) for v in best_lambda_row["fit"]["params_vec"]],
                "dev_logk": float(best_lambda_row["fit"]["dev_logk"]),
            },
            "full_refit": {
                "train_avg_nll": float(best_full["fit"]["train_avg_nll"]),
                "train_loglike": float(best_full["fit"]["train_loglike"]),
                "bic": float(best_full["fit"]["bic"]),
                "rmse": float(best_full["rmse"]),
                "step2_gain": float(best_step2["gain"]),
                "step2_helps": bool(best_step2["helps"]),
                "params": best_full["fit"]["params"],
                "params_vec": [float(v) for v in best_full["fit"]["params_vec"]],
                "dev_logk": float(best_full["fit"]["dev_logk"]),
            },
            "protocol_note": "selection.* are split-fit quantities fit only on the train split and are canonical for held-out comparisons; full_refit.* are full-corpus diagnostics only.",
        },
        "delta_best_vs_free_test": float(best_lambda_row["test_avg_nll"] - base_row["models"]["seam"]["test_avg_nll"]),
    }


def attach_fixedk(rows: list[dict], reg_map: dict[str, dict]):
    for row in rows:
        fixed = reg_map[row["slug"]]["variants"]["fixed_k"]
        row["fixed_k"] = {
            "test_avg_nll": float(fixed["test_avg_nll"]),
            "bic": float(fixed["bic"]),
            "rmse": float(fixed["rmse"]),
            "step2_gain": float(fixed["step2_gain"]),
            "step2_helps": bool(fixed["step2_helps"]),
            "params": fixed["params"],
        }


def build_summary(rows: list[dict]) -> dict:
    lambda_counts = {}
    for row in rows:
        key = f"{row['best_lambda']['lambda']:.6g}"
        lambda_counts[key] = lambda_counts.get(key, 0) + 1
    return {
        "rows": rows,
        "counts": {
            "n_rows": len(rows),
            "best_lambda_counts": lambda_counts,
            "softk_beats_free_test": int(sum(row["delta_best_vs_free_test"] < 0.0 for row in rows)),
            "softk_beats_fixedk_test": int(sum(row["best_lambda"]["selection"]["test_avg_nll"] < row["fixed_k"]["test_avg_nll"] for row in rows)),
            "softk_step2_help_count": int(sum(row["best_lambda"]["full_refit"]["step2_helps"] for row in rows)),
            "free_step2_help_count": int(sum(row["free"]["step2_helps"] for row in rows)),
            "fixedk_step2_help_count": int(sum(row["fixed_k"]["step2_helps"] for row in rows)),
        },
        "medians": {
            "delta_best_vs_free_test": float(np.median([row["delta_best_vs_free_test"] for row in rows])),
            "delta_best_vs_fixedk_test": float(np.median([row["best_lambda"]["selection"]["test_avg_nll"] - row["fixed_k"]["test_avg_nll"] for row in rows])),
            "softk_step2_gain": float(np.median([row["best_lambda"]["full_refit"]["step2_gain"] for row in rows])),
            "free_step2_gain": float(np.median([row["free"]["step2_gain"] for row in rows])),
            "fixedk_step2_gain": float(np.median([row["fixed_k"]["step2_gain"] for row in rows])),
        },
    }


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "best_lambda",
        "free_test_avg_nll",
        "fixed_k_test_avg_nll",
        "softk_test_avg_nll",
        "free_step2_gain",
        "fixed_k_step2_gain",
        "softk_step2_gain",
        "softk_transition_fraction",
        "softk_w",
        "softk_dev_logk",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "best_lambda": row["best_lambda"]["lambda"],
                    "free_test_avg_nll": row["free"]["test_avg_nll"],
                    "fixed_k_test_avg_nll": row["fixed_k"]["test_avg_nll"],
                    "softk_test_avg_nll": row["best_lambda"]["selection"]["test_avg_nll"],
                    "free_step2_gain": row["free"]["step2_gain"],
                    "fixed_k_step2_gain": row["fixed_k"]["step2_gain"],
                    "softk_step2_gain": row["best_lambda"]["full_refit"]["step2_gain"],
                    "softk_transition_fraction": row["best_lambda"]["selection"]["params"]["transition_fraction"],
                    "softk_w": row["best_lambda"]["selection"]["params"]["w"],
                    "softk_dev_logk": row["best_lambda"]["selection"]["dev_logk"],
                }
            )


def plot_lambda_tradeoff(rows: list[dict], path: Path):
    if plt is None:
        return
    lam_values = LAMBDAS
    median_deltas = []
    chosen_counts = []
    for lam in lam_values:
        vals = []
        count = 0
        for row in rows:
            match = next((item for item in row["sweep"] if abs(item["lambda"] - lam) < 1e-15), None)
            if match is not None:
                vals.append(match["test_avg_nll"] - row["free"]["test_avg_nll"])
            if abs(row["best_lambda"]["lambda"] - lam) < 1e-15:
                count += 1
        median_deltas.append(float(np.median(vals)) if vals else float("nan"))
        chosen_counts.append(count)
    fig, ax1 = plt.subplots(figsize=(8.5, 5.0))
    ax1.plot(lam_values, median_deltas, marker="o", color="#1565c0")
    ax1.set_xscale("log")
    ax1.set_xlabel("soft-k prior lambda")
    ax1.set_ylabel("median test avg NLL delta vs free", color="#1565c0")
    ax1.tick_params(axis="y", labelcolor="#1565c0")
    ax1.grid(True, alpha=0.25)
    ax2 = ax1.twinx()
    ax2.bar(lam_values, chosen_counts, width=[lam * 0.25 for lam in lam_values], alpha=0.25, color="#c62828")
    ax2.set_ylabel("best-lambda corpus count", color="#c62828")
    ax2.tick_params(axis="y", labelcolor="#c62828")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def build_report(summary: dict) -> str:
    counts = summary["counts"]
    medians = summary["medians"]
    rows = sorted(summary["rows"], key=lambda row: row["delta_best_vs_free_test"])
    lines = [
        "# Soft-k Seam-Mandelbrot PMF",
        "",
        "- `k` stays free, but the train objective adds a quadratic penalty on `(log k - alpha log V)^2`.",
        f"- prior center: `alpha = {SMOOTH_ALPHA:.12f}`",
        f"- lambda sweep: `{LAMBDAS}`",
        "",
        f"- best-lambda counts: `{json.dumps(counts['best_lambda_counts'])}`",
        f"- soft-k beats free on held-out avg NLL: `{counts['softk_beats_free_test']}` / `{counts['n_rows']}`",
        f"- soft-k beats fixed-k on held-out avg NLL: `{counts['softk_beats_fixedk_test']}` / `{counts['n_rows']}`",
        f"- median soft-k minus free held-out avg NLL: `{medians['delta_best_vs_free_test']:.12f}`",
        f"- median soft-k minus fixed-k held-out avg NLL: `{medians['delta_best_vs_fixedk_test']:.12f}`",
        f"- free step-2 help count: `{counts['free_step2_help_count']}` / `{counts['n_rows']}`",
        f"- fixed-k step-2 help count: `{counts['fixedk_step2_help_count']}` / `{counts['n_rows']}`",
        f"- soft-k step-2 help count: `{counts['softk_step2_help_count']}` / `{counts['n_rows']}`",
        f"- median free step-2 gain: `{medians['free_step2_gain']:.12f}`",
        f"- median fixed-k step-2 gain: `{medians['fixedk_step2_gain']:.12f}`",
        f"- median soft-k step-2 gain: `{medians['softk_step2_gain']:.12f}`",
        "",
        "| corpus | best lambda | free | fixed-k | soft-k | free step2 | fixed-k step2 | soft-k step2 | frac | w |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['best_lambda']['lambda']:.6g} | "
            f"{row['free']['test_avg_nll']:.6f} | {row['fixed_k']['test_avg_nll']:.6f} | {row['best_lambda']['selection']['test_avg_nll']:.6f} | "
            f"{row['free']['step2_gain']:.6f} | {row['fixed_k']['step2_gain']:.6f} | {row['best_lambda']['full_refit']['step2_gain']:.6f} | "
            f"{row['best_lambda']['selection']['params']['transition_fraction']:.3f} | {row['best_lambda']['selection']['params']['w']:.3f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    reranked_map = load_reranked_map()
    base_map = load_base_map()
    reg_map = load_reg_map()
    reused_lambdas = load_reused_lambda_map(args.reuse_best_lambdas_summary)
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(
            analyze_corpus(
                spec,
                reranked_map.get(spec["slug"]),
                base_map[spec["slug"]],
                corpus_index,
                forced_lambda=reused_lambdas.get(spec["slug"]),
            )
        )
    attach_fixedk(rows, reg_map)
    summary = build_summary(rows)
    summary["metadata"] = {
        "outdir": str(args.outdir),
        "reuse_best_lambdas_summary": None if args.reuse_best_lambdas_summary is None else str(args.reuse_best_lambdas_summary),
        "canonical_heldout_protocol": "best_lambda.selection.* are split-fit quantities and are canonical for held-out comparisons",
    }
    (args.outdir / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (args.outdir / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, args.outdir / "softk_table.csv")
    plot_lambda_tradeoff(rows, args.outdir / "softk_lambda_tradeoff.png")


if __name__ == "__main__":
    main()
