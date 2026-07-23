from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp, zeta as scipy_zeta

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    matplotlib = None
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
EML_PATH = ROOT / "eml_zipf_enriched_search.py"
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_seam_mandelbrot_pmf"

TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416

ZIPF_ALPHA_BOUNDS = (0.05, 5.0)
ZM_B_BOUNDS = (0.05, 5.0)
ZM_C_BOUNDS = (0.0, 5000.0)
MOE_ALPHA_BOUNDS = (1.000001, 6.0)
MOE_BETA_BOUNDS = (1e-3, 1e2)
SEAM_B_BOUNDS = (0.05, 5.0)
SEAM_C_BOUNDS = (0.0, 5000.0)
SEAM_W_BOUNDS = (0.05, 5.0)
SEAM_STARTS = 28
STEP2_VARIANCE_THRESHOLD = 1e-10


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_mandelbrot_common")
eml = load_module(EML_PATH, "zipf_seam_mandelbrot_eml")


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


def stable_train_order(words: list[str], full_counts: np.ndarray, train_counts: np.ndarray) -> np.ndarray:
    order = sorted(
        range(len(words)),
        key=lambda idx: (-int(train_counts[idx]), -int(full_counts[idx]), words[idx]),
    )
    return np.asarray(order, dtype=np.int64)


def split_counts_by_train_rank(dataset: dict, seed: int = SPLIT_SEED, train_fraction: float = TRAIN_FRACTION) -> dict:
    words = [word for word, _ in dataset["ranked"]]
    full_counts = dataset["freqs"].astype(np.int64)
    rng = np.random.default_rng(seed)
    train_counts = rng.binomial(full_counts, train_fraction).astype(np.int64)
    test_counts = full_counts - train_counts
    order = stable_train_order(words, full_counts, train_counts)
    ordered_words = [words[idx] for idx in order]
    return {
        "words": ordered_words,
        "full_counts": full_counts[order].astype(np.float64),
        "train_counts": train_counts[order].astype(np.float64),
        "test_counts": test_counts[order].astype(np.float64),
        "order": order,
        "token_count_full": int(full_counts.sum()),
        "token_count_train": int(train_counts.sum()),
        "token_count_test": int(test_counts.sum()),
        "vocab_size": int(len(full_counts)),
    }


def sigma_curve(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = np.clip((np.log(ranks) - math.log(k)) / w, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def smooth_tail_local_rank(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    scale = max(1.0, k * w)
    z = np.clip((ranks - k) / scale, -60.0, 60.0)
    return 1.0 + scale * np.log1p(np.exp(z))


def seam_log_scores(ranks: np.ndarray, params: np.ndarray) -> np.ndarray | None:
    b1, c1, b2, c2, k, w = [float(v) for v in params]
    if b1 <= 0.0 or b2 <= 0.0 or c1 < 0.0 or c2 < 0.0 or k <= 0.0 or w <= 0.0:
        return None
    sigma = sigma_curve(ranks, k, w)
    head = -b1 * np.log(ranks + c1)
    tail_rank = smooth_tail_local_rank(ranks, k, w)
    tail_at_k = float(smooth_tail_local_rank(np.array([k], dtype=np.float64), k, w)[0])
    head_at_k = -b1 * math.log(k + c1)
    d = head_at_k + b2 * math.log(tail_at_k + c2)
    tail = d - b2 * np.log(tail_rank + c2)
    scores = sigma * head + (1.0 - sigma) * tail
    return scores if np.all(np.isfinite(scores)) else None


def riemann_zeta(alpha: float) -> float:
    return float(scipy_zeta(alpha, 1.0))


def moe_base_logpmf(alpha: float, beta: float, x: np.ndarray) -> np.ndarray | None:
    z = riemann_zeta(alpha)
    if not math.isfinite(z) or z <= 0.0 or beta <= 0.0:
        return None
    hz_x = np.asarray(scipy_zeta(alpha, x), dtype=np.float64)
    hz_x1 = np.asarray(scipy_zeta(alpha, x + 1.0), dtype=np.float64)
    beta_bar = 1.0 - beta
    denom_x = z - beta_bar * hz_x
    denom_x1 = z - beta_bar * hz_x1
    if np.any(~np.isfinite(denom_x)) or np.any(~np.isfinite(denom_x1)) or np.any(denom_x <= 0.0) or np.any(denom_x1 <= 0.0):
        return None
    with np.errstate(all="ignore"):
        logp = -alpha * np.log(x) + math.log(beta) + math.log(z) - np.log(denom_x) - np.log(denom_x1)
    if np.any(~np.isfinite(logp)):
        return None
    return logp


def normalized_logpmf_from_scores(log_scores: np.ndarray | None) -> np.ndarray | None:
    if log_scores is None or np.any(~np.isfinite(log_scores)):
        return None
    return log_scores - logsumexp(log_scores)


def zipf_logpmf(alpha: float, ranks: np.ndarray) -> np.ndarray | None:
    if alpha <= 0.0:
        return None
    return normalized_logpmf_from_scores(-alpha * np.log(ranks))


def zm_logpmf(b: float, c: float, ranks: np.ndarray) -> np.ndarray | None:
    if b <= 0.0 or c < 0.0:
        return None
    return normalized_logpmf_from_scores(-b * np.log(ranks + c))


def moe_logpmf(alpha: float, beta: float, ranks: np.ndarray) -> np.ndarray | None:
    return normalized_logpmf_from_scores(moe_base_logpmf(alpha, beta, ranks))


def seam_logpmf(params: np.ndarray, ranks: np.ndarray) -> np.ndarray | None:
    return normalized_logpmf_from_scores(seam_log_scores(ranks, params))


def multinomial_loglike(counts: np.ndarray, logpmf: np.ndarray | None) -> float:
    if logpmf is None or np.any(~np.isfinite(logpmf)):
        return -math.inf
    return float(np.dot(counts.astype(np.float64), logpmf.astype(np.float64)))


def full_rank_rmse(counts: np.ndarray, logpmf: np.ndarray | None) -> float:
    if logpmf is None:
        return float("inf")
    log_counts = np.log(np.asarray(counts, dtype=np.float64))
    pred = math.log(float(np.sum(counts))) + logpmf
    return common.rmse(log_counts, pred)


def k_bounds(vocab_size: int) -> tuple[float, float]:
    return 5.0, float(min(max(vocab_size, 5), 5000))


def fit_zipf_rank(train_counts: np.ndarray, ranks: np.ndarray) -> dict:
    n = int(np.sum(train_counts))

    def objective(alpha: float) -> float:
        logpmf = zipf_logpmf(alpha, ranks)
        ll = multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    result = minimize_scalar(objective, bounds=ZIPF_ALPHA_BOUNDS, method="bounded", options={"xatol": 1e-6})
    alpha = float(result.x)
    logpmf = zipf_logpmf(alpha, ranks)
    ll = multinomial_loglike(train_counts, logpmf)
    return {
        "params": {"alpha": alpha},
        "logpmf": logpmf,
        "train_loglike": float(ll),
        "aic": float(2 - 2 * ll),
        "bic": float(math.log(n) * 1 - 2 * ll),
        "success": bool(result.success),
        "message": result.message,
        "p": 1,
    }


def fit_zm_rank(train_counts: np.ndarray, ranks: np.ndarray, dataset_hint: dict | None = None) -> dict:
    positive = train_counts > 0.0
    starts: list[np.ndarray] = []
    if dataset_hint is not None and positive.any():
        curve = common.fit_zipf_mandelbrot(ranks[positive], np.log(train_counts[positive]))
        starts.append(np.array([curve["b"], curve["c"]], dtype=np.float64))
    for b0 in [0.6, 0.9, 1.2, 1.6, 2.2]:
        for c0 in [0.0, 0.5, 2.0, 10.0, 50.0, 250.0]:
            starts.append(np.array([b0, c0], dtype=np.float64))

    best = None
    bounds = [ZM_B_BOUNDS, ZM_C_BOUNDS]

    def objective(theta: np.ndarray) -> float:
        logpmf = zm_logpmf(float(theta[0]), float(theta[1]), ranks)
        ll = multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    for start in starts:
        start = np.clip(start, [bounds[0][0], bounds[1][0]], [bounds[0][1], bounds[1][1]])
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        logpmf = zm_logpmf(float(result.x[0]), float(result.x[1]), ranks)
        ll = multinomial_loglike(train_counts, logpmf)
        record = {
            "params": {"b": float(result.x[0]), "c": float(result.x[1])},
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "p": 2,
        }
        if best is None or record["train_loglike"] > best["train_loglike"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def fit_moe_rank(train_counts: np.ndarray, ranks: np.ndarray, zipf_fit: dict) -> dict:
    starts = []
    for alpha0 in [zipf_fit["params"]["alpha"], 1.2, 1.5, 2.0, 3.0]:
        alpha0 = min(max(alpha0, MOE_ALPHA_BOUNDS[0]), MOE_ALPHA_BOUNDS[1])
        for beta0 in [0.3, 0.7, 1.0, 1.5, 2.5, 5.0]:
            starts.append(np.array([alpha0, beta0], dtype=np.float64))

    best = None
    bounds = [MOE_ALPHA_BOUNDS, MOE_BETA_BOUNDS]

    def objective(theta: np.ndarray) -> float:
        logpmf = moe_logpmf(float(theta[0]), float(theta[1]), ranks)
        ll = multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    for start in starts:
        start = np.clip(start, [bounds[0][0], bounds[1][0]], [bounds[0][1], bounds[1][1]])
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        logpmf = moe_logpmf(float(result.x[0]), float(result.x[1]), ranks)
        ll = multinomial_loglike(train_counts, logpmf)
        record = {
            "params": {"alpha": float(result.x[0]), "beta": float(result.x[1])},
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "p": 2,
        }
        if best is None or record["train_loglike"] > best["train_loglike"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def seam_starts(vocab_size: int, reranked_row: dict | None, seed: int) -> list[np.ndarray]:
    k_lo, k_hi = k_bounds(vocab_size)
    lower = np.array([SEAM_B_BOUNDS[0], SEAM_C_BOUNDS[0], SEAM_B_BOUNDS[0], SEAM_C_BOUNDS[0], k_lo, SEAM_W_BOUNDS[0]], dtype=np.float64)
    upper = np.array([SEAM_B_BOUNDS[1], SEAM_C_BOUNDS[1], SEAM_B_BOUNDS[1], SEAM_C_BOUNDS[1], k_hi, SEAM_W_BOUNDS[1]], dtype=np.float64)
    starts = []
    if reranked_row is not None:
        params = reranked_row["params"]
        starts.append(
            np.array(
                [
                    float(params["b1"]),
                    float(params["c1"]),
                    float(params["b2"]),
                    float(params["c2"]),
                    float(min(max(params["k"], k_lo), k_hi)),
                    float(params["w"]),
                ],
                dtype=np.float64,
            )
        )
    starts.append(np.array([1.2, 10.0, 1.8, 10.0, float(min(max(vocab_size ** 0.5, k_lo), k_hi)), 1.2], dtype=np.float64))
    starts.append(np.array([1.0, 1.0, 1.4, 100.0, float(min(max(vocab_size ** 0.52, k_lo), k_hi)), 0.8], dtype=np.float64))

    rng = np.random.default_rng(seed)
    while len(starts) < SEAM_STARTS:
        starts.append(
            np.array(
                [
                    rng.uniform(*SEAM_B_BOUNDS),
                    rng.uniform(*SEAM_C_BOUNDS),
                    rng.uniform(*SEAM_B_BOUNDS),
                    rng.uniform(*SEAM_C_BOUNDS),
                    float(np.exp(rng.uniform(math.log(k_lo), math.log(k_hi)))),
                    float(np.exp(rng.uniform(math.log(SEAM_W_BOUNDS[0]), math.log(SEAM_W_BOUNDS[1])))),
                ],
                dtype=np.float64,
            )
        )
    return [np.clip(start, lower, upper) for start in starts]


def fit_seam_rank(train_counts: np.ndarray, ranks: np.ndarray, vocab_size: int, reranked_row: dict | None, seed: int) -> dict:
    k_lo, k_hi = k_bounds(vocab_size)
    bounds = [
        SEAM_B_BOUNDS,
        SEAM_C_BOUNDS,
        SEAM_B_BOUNDS,
        SEAM_C_BOUNDS,
        (k_lo, k_hi),
        SEAM_W_BOUNDS,
    ]
    lower = np.array([b[0] for b in bounds], dtype=np.float64)
    upper = np.array([b[1] for b in bounds], dtype=np.float64)
    best = None

    def objective(theta: np.ndarray) -> float:
        logpmf = seam_logpmf(theta, ranks)
        ll = multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    for start_index, start in enumerate(seam_starts(vocab_size, reranked_row, seed=seed), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = seam_logpmf(params, ranks)
        ll = multinomial_loglike(train_counts, logpmf)
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
            "p": 6,
        }
        if best is None or record["train_loglike"] > best["train_loglike"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def evaluate_test_fit(fit: dict, test_counts: np.ndarray) -> dict:
    test_loglike = multinomial_loglike(test_counts, fit["logpmf"])
    token_count = max(int(np.sum(test_counts)), 1)
    return {
        "test_loglike": float(test_loglike),
        "test_avg_nll": float(-test_loglike / token_count),
    }


def fit_models_on_rank_counts(counts: np.ndarray, vocab_size: int, reranked_row: dict | None, seed: int, dataset_hint: dict | None = None) -> dict:
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    zipf_fit = fit_zipf_rank(counts, ranks)
    zm_fit = fit_zm_rank(counts, ranks, dataset_hint=dataset_hint)
    moe_fit = fit_moe_rank(counts, ranks, zipf_fit)
    seam_fit = fit_seam_rank(counts, ranks, vocab_size=vocab_size, reranked_row=reranked_row, seed=seed)
    return {
        "zipf": zipf_fit,
        "zm": zm_fit,
        "moe": moe_fit,
        "seam": seam_fit,
    }


def full_fit_context(dataset: dict, reranked_row: dict | None, seed: int) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    fits = fit_models_on_rank_counts(counts, int(dataset["unique_words"]), reranked_row=reranked_row, seed=seed, dataset_hint=dataset)
    contexts = {}
    for name, fit in fits.items():
        logpmf = fit["logpmf"]
        pred_log = math.log(float(np.sum(counts))) + logpmf
        contexts[name] = {
            "params": fit["params"],
            "train_loglike": fit["train_loglike"],
            "aic": fit["aic"],
            "bic": fit["bic"],
            "logpmf": logpmf,
            "rmse": full_rank_rmse(counts, logpmf),
            "prediction_log": pred_log,
        }
    return contexts


def exact_step2_on_residual(dataset: dict, prediction_log: np.ndarray) -> dict:
    x = common.normalize_x(dataset["log_rank"])
    y = dataset["log_freq"]
    residual = y - prediction_log
    vocab0 = eml.initial_vocabulary(x, residual)
    step1 = eml.generate_candidates(vocab0, residual, 1)
    step1 = eml.dedupe_candidates(step1)
    step1 = eml.filter_candidates(step1, STEP2_VARIANCE_THRESHOLD)
    vocab = vocab0 + step1
    step2 = eml.generate_candidates(vocab, residual, 2)
    step2 = eml.dedupe_candidates(step2)
    step2 = eml.filter_candidates(step2, STEP2_VARIANCE_THRESHOLD)
    winner = min(step2, key=lambda row: (row["rmse"], row["expr"]))
    composite_rmse = eml.rmse(y, prediction_log + winner["values"])
    baseline_rmse = common.rmse(y, prediction_log)
    return {
        "expr": winner["expr"],
        "math": winner["math"],
        "residual_rmse": float(winner["rmse"]),
        "composite_rmse": float(composite_rmse),
        "baseline_rmse": float(baseline_rmse),
        "gain": float(baseline_rmse - composite_rmse),
        "helps": bool(composite_rmse + 1e-12 < baseline_rmse),
    }


def analyze_corpus(spec: dict, reranked_row: dict | None, corpus_index: int) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]

    rank_fits = fit_models_on_rank_counts(
        train_counts,
        vocab_size=split["vocab_size"],
        reranked_row=reranked_row,
        seed=SPLIT_SEED + 1000 + corpus_index,
        dataset_hint=dataset,
    )
    evals = {name: evaluate_test_fit(fit, test_counts) for name, fit in rank_fits.items()}

    full_context = full_fit_context(dataset, reranked_row=reranked_row, seed=SPLIT_SEED + 2000 + corpus_index)
    seam_step2 = exact_step2_on_residual(dataset, full_context["seam"]["prediction_log"])

    train_bics = {name: fit["bic"] for name, fit in rank_fits.items()}
    test_nlls = {name: evals[name]["test_avg_nll"] for name in rank_fits}
    rmse_map = {name: full_context[name]["rmse"] for name in full_context}

    row = {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "split": {
            "train_tokens": split["token_count_train"],
            "test_tokens": split["token_count_test"],
        },
        "models": {
            name: {
                "params": fit["params"],
                "train_loglike": fit["train_loglike"],
                "test_loglike": evals[name]["test_loglike"],
                "test_avg_nll": evals[name]["test_avg_nll"],
                "aic": fit["aic"],
                "bic": fit["bic"],
                "full_rmse": full_context[name]["rmse"],
            }
            for name, fit in rank_fits.items()
        },
        "winners": {
            "test_avg_nll": min(test_nlls, key=test_nlls.get),
            "train_bic": min(train_bics, key=train_bics.get),
        },
        "seam_step2": seam_step2,
        "comparisons": {
            "seam_minus_moe_test_avg_nll": float(evals["seam"]["test_avg_nll"] - evals["moe"]["test_avg_nll"]),
            "seam_minus_zm_test_avg_nll": float(evals["seam"]["test_avg_nll"] - evals["zm"]["test_avg_nll"]),
            "seam_minus_moe_rmse": float(rmse_map["seam"] - rmse_map["moe"]),
        },
    }
    return row


def build_summary(rows: list[dict]) -> dict:
    test_winners = CounterLike([row["winners"]["test_avg_nll"] for row in rows])
    bic_winners = CounterLike([row["winners"]["train_bic"] for row in rows])
    seam_better_than_moe_test = [row["comparisons"]["seam_minus_moe_test_avg_nll"] < 0.0 for row in rows]
    seam_better_than_zm_test = [row["comparisons"]["seam_minus_zm_test_avg_nll"] < 0.0 for row in rows]
    seam_step2_gains = [row["seam_step2"]["gain"] for row in rows]
    return {
        "rows": rows,
        "counts": {
            "n_rows": len(rows),
            "test_winner_counts": dict(test_winners),
            "bic_winner_counts": dict(bic_winners),
            "seam_beats_moe_test": int(sum(seam_better_than_moe_test)),
            "seam_beats_zm_test": int(sum(seam_better_than_zm_test)),
            "seam_step2_help_count": int(sum(row["seam_step2"]["helps"] for row in rows)),
        },
        "medians": {
            "seam_minus_moe_test_avg_nll": float(np.median([row["comparisons"]["seam_minus_moe_test_avg_nll"] for row in rows])),
            "seam_minus_zm_test_avg_nll": float(np.median([row["comparisons"]["seam_minus_zm_test_avg_nll"] for row in rows])),
            "seam_step2_gain": float(np.median(seam_step2_gains)),
            "seam_transition_fraction": float(np.median([row["models"]["seam"]["params"]["transition_fraction"] for row in rows])),
        },
    }


class CounterLike(dict):
    def __init__(self, items):
        super().__init__()
        for item in items:
            self[item] = self.get(item, 0) + 1


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "token_count",
        "vocab_size",
        "winner_test_avg_nll",
        "winner_train_bic",
        "zipf_test_avg_nll",
        "zm_test_avg_nll",
        "moe_test_avg_nll",
        "seam_test_avg_nll",
        "zipf_bic",
        "zm_bic",
        "moe_bic",
        "seam_bic",
        "zipf_rmse",
        "zm_rmse",
        "moe_rmse",
        "seam_rmse",
        "seam_transition_fraction",
        "seam_step2_expr",
        "seam_step2_gain",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "token_count": row["token_count"],
                    "vocab_size": row["vocab_size"],
                    "winner_test_avg_nll": row["winners"]["test_avg_nll"],
                    "winner_train_bic": row["winners"]["train_bic"],
                    "zipf_test_avg_nll": row["models"]["zipf"]["test_avg_nll"],
                    "zm_test_avg_nll": row["models"]["zm"]["test_avg_nll"],
                    "moe_test_avg_nll": row["models"]["moe"]["test_avg_nll"],
                    "seam_test_avg_nll": row["models"]["seam"]["test_avg_nll"],
                    "zipf_bic": row["models"]["zipf"]["bic"],
                    "zm_bic": row["models"]["zm"]["bic"],
                    "moe_bic": row["models"]["moe"]["bic"],
                    "seam_bic": row["models"]["seam"]["bic"],
                    "zipf_rmse": row["models"]["zipf"]["full_rmse"],
                    "zm_rmse": row["models"]["zm"]["full_rmse"],
                    "moe_rmse": row["models"]["moe"]["full_rmse"],
                    "seam_rmse": row["models"]["seam"]["full_rmse"],
                    "seam_transition_fraction": row["models"]["seam"]["params"]["transition_fraction"],
                    "seam_step2_expr": row["seam_step2"]["expr"],
                    "seam_step2_gain": row["seam_step2"]["gain"],
                }
            )


def plot_test_scatter(rows: list[dict], path: Path):
    if plt is None:
        return
    x = np.array([row["models"]["moe"]["test_avg_nll"] for row in rows], dtype=np.float64)
    y = np.array([row["models"]["seam"]["test_avg_nll"] for row in rows], dtype=np.float64)
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    ax.scatter(x, y, color="#1565c0", s=36)
    lo = min(float(np.min(x)), float(np.min(y)))
    hi = max(float(np.max(x)), float(np.max(y)))
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="#666666", linewidth=1.0)
    for row, x0, y0 in zip(rows, x, y):
        ax.annotate(row["slug"], (x0, y0), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("MOE held-out avg NLL")
    ax.set_ylabel("Seam-Mandelbrot held-out avg NLL")
    ax.set_title("Held-out Likelihood: Seam-Mandelbrot vs MOE")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_step2_gains(rows: list[dict], path: Path):
    if plt is None:
        return
    gains = np.array([row["seam_step2"]["gain"] for row in rows], dtype=np.float64)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.hist(gains, bins=10, color="#c62828", edgecolor="white")
    ax.axvline(0.0, color="#333333", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Seam-PMF step-2 RMSE gain")
    ax.set_ylabel("count")
    ax.set_title("Residual Step-2 Gain After Seam-Mandelbrot PMF")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def build_report(summary: dict) -> str:
    counts = summary["counts"]
    medians = summary["medians"]
    rows = sorted(summary["rows"], key=lambda row: row["models"]["seam"]["test_avg_nll"] - row["models"]["moe"]["test_avg_nll"])
    lines = [
        "# Seam-Mandelbrot PMF Comparison",
        "",
        "- Finite-support rank PMFs were fit by train likelihood on an 80/20 token split induced by binomially splitting each type count.",
        "- Models: truncated Zipf, truncated Zipf-Mandelbrot, truncated MOEZipf, and the new finite-support Seam-Mandelbrot PMF.",
        "- The Seam-Mandelbrot PMF uses a smooth transition between two ZM-like slopes and is normalized directly over ranks `1..V`.",
        "",
        f"- test winner counts: `{json.dumps(counts['test_winner_counts'])}`",
        f"- train BIC winner counts: `{json.dumps(counts['bic_winner_counts'])}`",
        f"- Seam beats MOE on held-out avg NLL: `{counts['seam_beats_moe_test']}` / `{counts['n_rows']}`",
        f"- Seam beats ZM on held-out avg NLL: `{counts['seam_beats_zm_test']}` / `{counts['n_rows']}`",
        f"- Seam residual step-2 help count: `{counts['seam_step2_help_count']}` / `{counts['n_rows']}`",
        f"- median Seam minus MOE held-out avg NLL: `{medians['seam_minus_moe_test_avg_nll']:.12f}`",
        f"- median Seam minus ZM held-out avg NLL: `{medians['seam_minus_zm_test_avg_nll']:.12f}`",
        f"- median Seam step-2 gain: `{medians['seam_step2_gain']:.12f}`",
        f"- median Seam transition fraction: `{medians['seam_transition_fraction']:.12f}`",
        "",
        "| corpus | winner (test) | Zipf | ZM | MOE | Seam | winner (BIC) | Seam RMSE | MOE RMSE | Seam step-2 | gain | frac |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['winners']['test_avg_nll']} | "
            f"{row['models']['zipf']['test_avg_nll']:.6f} | {row['models']['zm']['test_avg_nll']:.6f} | "
            f"{row['models']['moe']['test_avg_nll']:.6f} | {row['models']['seam']['test_avg_nll']:.6f} | "
            f"{row['winners']['train_bic']} | {row['models']['seam']['full_rmse']:.6f} | {row['models']['moe']['full_rmse']:.6f} | "
            f"{row['seam_step2']['expr']} | {row['seam_step2']['gain']:.6f} | "
            f"{row['models']['seam']['params']['transition_fraction']:.3f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_map = load_reranked_map()
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(analyze_corpus(spec, reranked_row=reranked_map.get(spec["slug"]), corpus_index=corpus_index))
    summary = build_summary(rows)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, OUTDIR / "seam_mandelbrot_table.csv")
    plot_test_scatter(rows, OUTDIR / "seam_vs_moe_test_nll.png")
    plot_step2_gains(rows, OUTDIR / "seam_step2_gains.png")


if __name__ == "__main__":
    main()
