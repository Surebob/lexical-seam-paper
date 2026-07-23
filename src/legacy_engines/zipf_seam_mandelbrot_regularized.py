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
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_seam_mandelbrot_regularized"

SMOOTH_ALPHA = 0.521422287446811
SMOOTH_W = 0.9899928130291
TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
STARTS = 18


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_regularized_common")
base = load_module(BASE_PATH, "zipf_seam_regularized_base")


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


def fixed_k(vocab_size: int) -> float:
    return float(vocab_size ** SMOOTH_ALPHA)


def fixed_k_bounds(vocab_size: int) -> tuple[float, float]:
    return base.k_bounds(vocab_size)


def starts_fixed_k(vocab_size: int, reranked_row: dict | None, seed: int) -> list[np.ndarray]:
    lower = np.array([base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_W_BOUNDS[0]], dtype=np.float64)
    upper = np.array([base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_W_BOUNDS[1]], dtype=np.float64)
    starts = []
    if reranked_row is not None:
        params = reranked_row["params"]
        starts.append(np.array([float(params["b1"]), float(params["c1"]), float(params["b2"]), float(params["c2"]), float(params["w"])], dtype=np.float64))
    starts.append(np.array([1.2, 10.0, 1.8, 10.0, 1.2], dtype=np.float64))
    starts.append(np.array([1.0, 1.0, 1.4, 100.0, 0.8], dtype=np.float64))
    rng = np.random.default_rng(seed)
    while len(starts) < STARTS:
        starts.append(
            np.array(
                [
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                ],
                dtype=np.float64,
            )
        )
    return [np.clip(start, lower, upper) for start in starts]


def starts_fixed_kw(reranked_row: dict | None, seed: int) -> list[np.ndarray]:
    lower = np.array([base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0]], dtype=np.float64)
    upper = np.array([base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1]], dtype=np.float64)
    starts = []
    if reranked_row is not None:
        params = reranked_row["params"]
        starts.append(np.array([float(params["b1"]), float(params["c1"]), float(params["b2"]), float(params["c2"])], dtype=np.float64))
    starts.append(np.array([1.2, 10.0, 1.8, 10.0], dtype=np.float64))
    starts.append(np.array([1.0, 1.0, 1.4, 100.0], dtype=np.float64))
    rng = np.random.default_rng(seed)
    while len(starts) < STARTS:
        starts.append(
            np.array(
                [
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                ],
                dtype=np.float64,
            )
        )
    return [np.clip(start, lower, upper) for start in starts]


def fixed_k_logpmf(theta: np.ndarray, ranks: np.ndarray, vocab_size: int):
    k = min(max(fixed_k(vocab_size), fixed_k_bounds(vocab_size)[0]), fixed_k_bounds(vocab_size)[1])
    params = np.array([float(theta[0]), float(theta[1]), float(theta[2]), float(theta[3]), k, float(theta[4])], dtype=np.float64)
    return base.seam_logpmf(params, ranks), params


def fixed_kw_logpmf(theta: np.ndarray, ranks: np.ndarray, vocab_size: int):
    k = min(max(fixed_k(vocab_size), fixed_k_bounds(vocab_size)[0]), fixed_k_bounds(vocab_size)[1])
    params = np.array([float(theta[0]), float(theta[1]), float(theta[2]), float(theta[3]), k, float(SMOOTH_W)], dtype=np.float64)
    return base.seam_logpmf(params, ranks), params


def fit_fixed_k(train_counts: np.ndarray, ranks: np.ndarray, vocab_size: int, reranked_row: dict | None, seed: int) -> dict:
    bounds = [
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.SEAM_W_BOUNDS,
    ]
    best = None

    def objective(theta: np.ndarray) -> float:
        logpmf, _ = fixed_k_logpmf(theta, ranks, vocab_size)
        ll = base.multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    for start_index, start in enumerate(starts_fixed_k(vocab_size, reranked_row, seed), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        logpmf, params = fixed_k_logpmf(result.x, ranks, vocab_size)
        ll = base.multinomial_loglike(train_counts, logpmf)
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
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "start_index": int(start_index),
            "p": 5,
        }
        if best is None or record["train_loglike"] > best["train_loglike"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def fit_fixed_kw(train_counts: np.ndarray, ranks: np.ndarray, vocab_size: int, reranked_row: dict | None, seed: int) -> dict:
    bounds = [
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
    ]
    best = None

    def objective(theta: np.ndarray) -> float:
        logpmf, _ = fixed_kw_logpmf(theta, ranks, vocab_size)
        ll = base.multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    for start_index, start in enumerate(starts_fixed_kw(reranked_row, seed), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        logpmf, params = fixed_kw_logpmf(result.x, ranks, vocab_size)
        ll = base.multinomial_loglike(train_counts, logpmf)
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
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "start_index": int(start_index),
            "p": 4,
        }
        if best is None or record["train_loglike"] > best["train_loglike"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def evaluate_fit(fit: dict, test_counts: np.ndarray) -> dict:
    return base.evaluate_test_fit(fit, test_counts)


def full_context_variant(dataset: dict, reranked_row: dict | None, variant: str, seed: int) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    vocab_size = int(dataset["unique_words"])
    if variant == "fixed_k":
        fit = fit_fixed_k(counts, ranks, vocab_size, reranked_row, seed)
    elif variant == "fixed_kw":
        fit = fit_fixed_kw(counts, ranks, vocab_size, reranked_row, seed)
    else:
        raise ValueError(variant)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {
        "params": fit["params"],
        "train_loglike": fit["train_loglike"],
        "aic": fit["aic"],
        "bic": fit["bic"],
        "logpmf": fit["logpmf"],
        "rmse": base.full_rank_rmse(counts, fit["logpmf"]),
        "prediction_log": pred_log,
    }


def analyze_corpus(spec: dict, reranked_row: dict | None, base_row: dict, corpus_index: int) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    ranks = np.arange(1, len(train_counts) + 1, dtype=np.float64)
    vocab_size = split["vocab_size"]

    fixed_k_fit = fit_fixed_k(train_counts, ranks, vocab_size, reranked_row, seed=SPLIT_SEED + 3000 + corpus_index)
    fixed_kw_fit = fit_fixed_kw(train_counts, ranks, vocab_size, reranked_row, seed=SPLIT_SEED + 4000 + corpus_index)
    fixed_k_eval = evaluate_fit(fixed_k_fit, test_counts)
    fixed_kw_eval = evaluate_fit(fixed_kw_fit, test_counts)

    fixed_k_full = full_context_variant(dataset, reranked_row, "fixed_k", seed=SPLIT_SEED + 5000 + corpus_index)
    fixed_kw_full = full_context_variant(dataset, reranked_row, "fixed_kw", seed=SPLIT_SEED + 6000 + corpus_index)
    fixed_k_step2 = base.exact_step2_on_residual(dataset, fixed_k_full["prediction_log"])
    fixed_kw_step2 = base.exact_step2_on_residual(dataset, fixed_kw_full["prediction_log"])

    variants = {
        "free": {
            "test_avg_nll": float(base_row["models"]["seam"]["test_avg_nll"]),
            "bic": float(base_row["models"]["seam"]["bic"]),
            "rmse": float(base_row["models"]["seam"]["full_rmse"]),
            "step2_gain": float(base_row["seam_step2"]["gain"]),
            "step2_helps": bool(base_row["seam_step2"]["helps"]),
            "params": base_row["models"]["seam"]["params"],
        },
        "fixed_k": {
            "test_avg_nll": float(fixed_k_eval["test_avg_nll"]),
            "bic": float(fixed_k_fit["bic"]),
            "rmse": float(fixed_k_full["rmse"]),
            "step2_gain": float(fixed_k_step2["gain"]),
            "step2_helps": bool(fixed_k_step2["helps"]),
            "params": fixed_k_fit["params"],
        },
        "fixed_kw": {
            "test_avg_nll": float(fixed_kw_eval["test_avg_nll"]),
            "bic": float(fixed_kw_fit["bic"]),
            "rmse": float(fixed_kw_full["rmse"]),
            "step2_gain": float(fixed_kw_step2["gain"]),
            "step2_helps": bool(fixed_kw_step2["helps"]),
            "params": fixed_kw_fit["params"],
        },
    }

    winner_test = min(variants, key=lambda key: variants[key]["test_avg_nll"])
    winner_bic = min(variants, key=lambda key: variants[key]["bic"])

    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "variants": variants,
        "winner_test": winner_test,
        "winner_bic": winner_bic,
        "delta_fixed_k_vs_free_test": float(variants["fixed_k"]["test_avg_nll"] - variants["free"]["test_avg_nll"]),
        "delta_fixed_kw_vs_free_test": float(variants["fixed_kw"]["test_avg_nll"] - variants["free"]["test_avg_nll"]),
    }


def build_summary(rows: list[dict]) -> dict:
    winner_test_counts = {}
    winner_bic_counts = {}
    for row in rows:
        winner_test_counts[row["winner_test"]] = winner_test_counts.get(row["winner_test"], 0) + 1
        winner_bic_counts[row["winner_bic"]] = winner_bic_counts.get(row["winner_bic"], 0) + 1
    return {
        "rows": rows,
        "counts": {
            "n_rows": len(rows),
            "winner_test_counts": winner_test_counts,
            "winner_bic_counts": winner_bic_counts,
            "fixed_k_beats_free_test": int(sum(row["delta_fixed_k_vs_free_test"] < 0.0 for row in rows)),
            "fixed_kw_beats_free_test": int(sum(row["delta_fixed_kw_vs_free_test"] < 0.0 for row in rows)),
            "free_step2_help_count": int(sum(row["variants"]["free"]["step2_helps"] for row in rows)),
            "fixed_k_step2_help_count": int(sum(row["variants"]["fixed_k"]["step2_helps"] for row in rows)),
            "fixed_kw_step2_help_count": int(sum(row["variants"]["fixed_kw"]["step2_helps"] for row in rows)),
        },
        "medians": {
            "delta_fixed_k_vs_free_test": float(np.median([row["delta_fixed_k_vs_free_test"] for row in rows])),
            "delta_fixed_kw_vs_free_test": float(np.median([row["delta_fixed_kw_vs_free_test"] for row in rows])),
            "free_step2_gain": float(np.median([row["variants"]["free"]["step2_gain"] for row in rows])),
            "fixed_k_step2_gain": float(np.median([row["variants"]["fixed_k"]["step2_gain"] for row in rows])),
            "fixed_kw_step2_gain": float(np.median([row["variants"]["fixed_kw"]["step2_gain"] for row in rows])),
        },
    }


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "winner_test",
        "winner_bic",
        "free_test_avg_nll",
        "fixed_k_test_avg_nll",
        "fixed_kw_test_avg_nll",
        "free_bic",
        "fixed_k_bic",
        "fixed_kw_bic",
        "free_rmse",
        "fixed_k_rmse",
        "fixed_kw_rmse",
        "free_step2_gain",
        "fixed_k_step2_gain",
        "fixed_kw_step2_gain",
        "fixed_k_transition_fraction",
        "fixed_kw_transition_fraction",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "winner_test": row["winner_test"],
                    "winner_bic": row["winner_bic"],
                    "free_test_avg_nll": row["variants"]["free"]["test_avg_nll"],
                    "fixed_k_test_avg_nll": row["variants"]["fixed_k"]["test_avg_nll"],
                    "fixed_kw_test_avg_nll": row["variants"]["fixed_kw"]["test_avg_nll"],
                    "free_bic": row["variants"]["free"]["bic"],
                    "fixed_k_bic": row["variants"]["fixed_k"]["bic"],
                    "fixed_kw_bic": row["variants"]["fixed_kw"]["bic"],
                    "free_rmse": row["variants"]["free"]["rmse"],
                    "fixed_k_rmse": row["variants"]["fixed_k"]["rmse"],
                    "fixed_kw_rmse": row["variants"]["fixed_kw"]["rmse"],
                    "free_step2_gain": row["variants"]["free"]["step2_gain"],
                    "fixed_k_step2_gain": row["variants"]["fixed_k"]["step2_gain"],
                    "fixed_kw_step2_gain": row["variants"]["fixed_kw"]["step2_gain"],
                    "fixed_k_transition_fraction": row["variants"]["fixed_k"]["params"]["transition_fraction"],
                    "fixed_kw_transition_fraction": row["variants"]["fixed_kw"]["params"]["transition_fraction"],
                }
            )


def plot_scatter(rows: list[dict], path: Path):
    if plt is None:
        return
    x = np.array([row["variants"]["free"]["test_avg_nll"] for row in rows], dtype=np.float64)
    y = np.array([row["variants"]["fixed_k"]["test_avg_nll"] for row in rows], dtype=np.float64)
    z = np.array([row["variants"]["fixed_kw"]["test_avg_nll"] for row in rows], dtype=np.float64)
    lo = min(float(np.min(x)), float(np.min(y)), float(np.min(z)))
    hi = max(float(np.max(x)), float(np.max(y)), float(np.max(z)))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), constrained_layout=True)
    for ax, yy, title in [(axes[0], y, "Fixed k vs Free"), (axes[1], z, "Fixed k,w vs Free")]:
        ax.scatter(x, yy, s=32, color="#1565c0")
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="#666666", linewidth=1.0)
        for row0, xx0, yy0 in zip(rows, x, yy):
            ax.annotate(row0["slug"], (xx0, yy0), fontsize=7, xytext=(3, 2), textcoords="offset points")
        ax.set_xlabel("Free Seam held-out avg NLL")
        ax.set_ylabel("Regularized held-out avg NLL")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def build_report(summary: dict) -> str:
    counts = summary["counts"]
    medians = summary["medians"]
    rows = sorted(summary["rows"], key=lambda row: row["delta_fixed_k_vs_free_test"])
    lines = [
        "# Regularized Seam-Mandelbrot PMF",
        "",
        "- Prior structure came from the smooth reranked family, not from the free PMF fit.",
        f"- shared transition law: `k = V^{SMOOTH_ALPHA:.12f}`",
        f"- shared width constant for the fixed-kw variant: `w = {SMOOTH_W:.12f}`",
        "",
        f"- held-out winner counts: `{json.dumps(counts['winner_test_counts'])}`",
        f"- train BIC winner counts: `{json.dumps(counts['winner_bic_counts'])}`",
        f"- fixed-k beats free on held-out avg NLL: `{counts['fixed_k_beats_free_test']}` / `{counts['n_rows']}`",
        f"- fixed-k,w beats free on held-out avg NLL: `{counts['fixed_kw_beats_free_test']}` / `{counts['n_rows']}`",
        f"- median fixed-k minus free held-out avg NLL: `{medians['delta_fixed_k_vs_free_test']:.12f}`",
        f"- median fixed-k,w minus free held-out avg NLL: `{medians['delta_fixed_kw_vs_free_test']:.12f}`",
        f"- free step-2 help count: `{counts['free_step2_help_count']}` / `{counts['n_rows']}`",
        f"- fixed-k step-2 help count: `{counts['fixed_k_step2_help_count']}` / `{counts['n_rows']}`",
        f"- fixed-k,w step-2 help count: `{counts['fixed_kw_step2_help_count']}` / `{counts['n_rows']}`",
        f"- median free step-2 gain: `{medians['free_step2_gain']:.12f}`",
        f"- median fixed-k step-2 gain: `{medians['fixed_k_step2_gain']:.12f}`",
        f"- median fixed-k,w step-2 gain: `{medians['fixed_kw_step2_gain']:.12f}`",
        "",
        "| corpus | winner (test) | free | fixed-k | fixed-k,w | free step2 | fixed-k step2 | fixed-k,w step2 | frac |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['winner_test']} | "
            f"{row['variants']['free']['test_avg_nll']:.6f} | "
            f"{row['variants']['fixed_k']['test_avg_nll']:.6f} | "
            f"{row['variants']['fixed_kw']['test_avg_nll']:.6f} | "
            f"{row['variants']['free']['step2_gain']:.6f} | "
            f"{row['variants']['fixed_k']['step2_gain']:.6f} | "
            f"{row['variants']['fixed_kw']['step2_gain']:.6f} | "
            f"{row['variants']['fixed_k']['params']['transition_fraction']:.3f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_map = load_reranked_map()
    base_map = load_base_map()
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(analyze_corpus(spec, reranked_map.get(spec["slug"]), base_map[spec["slug"]], corpus_index))
    summary = build_summary(rows)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, OUTDIR / "regularized_seam_table.csv")
    plot_scatter(rows, OUTDIR / "regularized_vs_free_test_nll.png")


if __name__ == "__main__":
    main()
