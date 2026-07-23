from __future__ import annotations

import math
import statistics

import numpy as np
import psutil
from scipy.optimize import least_squares

import fit_config as cfg
from gate_functions import GATE_FUNCS, gate_width_scale


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    diff = np.asarray(y_true, dtype=np.float64) - np.asarray(y_pred, dtype=np.float64)
    return float(math.sqrt(float(np.mean(diff * diff))))


def _solve_affine(z: np.ndarray, y: np.ndarray):
    design = np.column_stack([np.ones_like(z), z])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ coeffs
    mse = float(np.mean((pred - y) ** 2))
    return float(coeffs[0]), float(coeffs[1]), pred, mse


def fit_zipf_mandelbrot(ranks: np.ndarray, log_freq: np.ndarray) -> dict:
    max_rank = float(np.max(ranks))
    c_grid = np.concatenate(
        [np.array([0.0], dtype=np.float64), np.geomspace(1e-6, max_rank, 4096, dtype=np.float64)]
    )
    best = None
    for c in c_grid:
        z = np.log(ranks + c)
        intercept, slope, pred, mse = _solve_affine(z, log_freq)
        candidate = {
            "a": intercept,
            "b": -slope,
            "c": float(c),
            "mse": mse,
            "rmse": float(math.sqrt(max(mse, 0.0))),
            "prediction": pred,
        }
        if best is None or candidate["mse"] < best["mse"]:
            best = candidate
    return best


def fit_piecewise_zm(dataset: dict, k: int) -> dict:
    freqs = dataset["freqs"]
    log_freq = dataset["log_freq"]
    if k <= 0 or k >= len(freqs):
        raise ValueError(f"Invalid breakpoint K={k} for vocab size {len(freqs)}")

    top_log_freq = log_freq[:k]
    top_ranks = np.arange(1, k + 1, dtype=np.float64)
    top_fit = fit_zipf_mandelbrot(top_ranks, top_log_freq)

    tail_log_freq = log_freq[k:]
    tail_ranks = np.arange(1, len(tail_log_freq) + 1, dtype=np.float64)
    tail_fit = fit_zipf_mandelbrot(tail_ranks, tail_log_freq)

    composite_pred = np.concatenate([top_fit["prediction"], tail_fit["prediction"]])
    return {
        "k": int(k),
        "top": {
            "a": float(top_fit["a"]),
            "b": float(top_fit["b"]),
            "c": float(top_fit["c"]),
            "rmse": float(top_fit["rmse"]),
        },
        "tail": {
            "a": float(tail_fit["a"]),
            "b": float(tail_fit["b"]),
            "c": float(tail_fit["c"]),
            "rmse": float(tail_fit["rmse"]),
        },
        "prediction": composite_pred,
        "rmse": rmse(log_freq, composite_pred),
    }


def bic_from_rmse(rmse_value: float, p: int, n: int) -> float:
    mse = float(rmse_value) ** 2
    return p * math.log(n) + n * math.log(mse)


def decoupled_bounds() -> tuple[np.ndarray, np.ndarray]:
    lower = np.array(
        [-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, cfg.K_LOWER, cfg.W_GATE_LOWER, cfg.W_TAIL_LOWER],
        dtype=np.float64,
    )
    upper = np.array(
        [100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, cfg.K_UPPER, cfg.W_GATE_UPPER, cfg.W_TAIL_UPPER],
        dtype=np.float64,
    )
    return lower, upper


def smooth_tail_local_rank(ranks: np.ndarray, k: float, w_tail: float) -> np.ndarray:
    scale = max(1.0, k * w_tail)
    z = np.clip((ranks - k) / scale, -60.0, 60.0)
    return 1.0 + scale * np.log1p(np.exp(z))


def decoupled_prediction(ranks: np.ndarray, log_ranks: np.ndarray, params: np.ndarray, gate_name: str) -> np.ndarray:
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    sigma = GATE_FUNCS[gate_name](log_ranks, k, w_gate)
    head = a1 - b1 * np.log(ranks + max(c1, 0.0))
    tail_rank = smooth_tail_local_rank(ranks, k, w_tail)
    tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
    return sigma * head + (1.0 - sigma) * tail


def summarize_params(params: np.ndarray, vocab_size: int) -> dict:
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    return {
        "a1": a1,
        "b1": b1,
        "c1": c1,
        "a2": a2,
        "b2": b2,
        "c2": c2,
        "k": k,
        "w_gate": w_gate,
        "w_tail": w_tail,
        "transition_fraction": float(np.log(k) / np.log(vocab_size)),
    }


def gate_seed(gate_name: str, base_seed: int) -> int:
    return base_seed + {"logistic": 0, "tanh": 1, "erf": 2, "algebraic": 3, "arctan": 4}[gate_name]


def make_random_initializations(dataset: dict, gate_name: str, piece: dict, n_starts: int, seed: int) -> list[np.ndarray]:
    lower, upper = decoupled_bounds()
    rng = np.random.default_rng(seed)
    y = dataset["log_freq"]
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    gate_scale = gate_width_scale(gate_name)
    anchor_shape_width = 0.5
    anchor = np.array(
        [
            piece["top"]["a"],
            piece["top"]["b"],
            piece["top"]["c"],
            piece["tail"]["a"],
            piece["tail"]["b"],
            min(piece["tail"]["c"], upper[5]),
            500.0,
            gate_scale * anchor_shape_width,
            anchor_shape_width,
        ],
        dtype=np.float64,
    )
    anchor = np.clip(anchor, lower + 1e-8, upper - 1e-8)

    starts = []
    shape_gate_lower = max(cfg.W_GATE_LOWER / gate_scale, 1e-6)
    shape_gate_upper = cfg.W_GATE_UPPER / gate_scale
    for idx in range(n_starts):
        if idx < n_starts // 2:
            x0 = np.empty(9, dtype=np.float64)
            x0[:6] = anchor[:6] + rng.normal(
                0.0,
                np.array([8.0, 0.35, 120.0, 8.0, 0.35, 120.0], dtype=np.float64),
            )
            x0[6] = anchor[6] + rng.normal(0.0, 180.0)
            shape_gate = anchor_shape_width * math.exp(float(rng.normal(0.0, 0.45)))
            x0[7] = gate_scale * shape_gate
            x0[8] = anchor[8] * math.exp(float(rng.normal(0.0, 0.45)))
        else:
            x0 = np.empty(9, dtype=np.float64)
            x0[0] = rng.uniform(max(lower[0], y_min - 3.0), min(upper[0], y_max + 3.0))
            x0[1] = rng.uniform(lower[1], upper[1])
            x0[2] = rng.uniform(lower[2], upper[2])
            x0[3] = rng.uniform(max(lower[3], y_min - 3.0), min(upper[3], y_max + 3.0))
            x0[4] = rng.uniform(lower[4], upper[4])
            x0[5] = rng.uniform(lower[5], upper[5])
            x0[6] = float(np.exp(rng.uniform(np.log(lower[6]), np.log(upper[6]))))
            shape_gate = float(np.exp(rng.uniform(np.log(shape_gate_lower), np.log(shape_gate_upper))))
            x0[7] = gate_scale * shape_gate
            x0[8] = float(np.exp(rng.uniform(np.log(lower[8]), np.log(upper[8]))))
        x0 = np.clip(x0, lower + 1e-8, upper - 1e-8)
        starts.append(x0)
    return starts


def fit_with_gate(
    dataset: dict,
    gate_name: str,
    base_seed: int,
    n_starts: int = cfg.DEFAULT_N_STARTS,
    max_nfev: int = cfg.DEFAULT_MAX_NFEV,
) -> dict:
    ranks = dataset["ranks"]
    log_ranks = dataset["log_rank"]
    y = dataset["log_freq"]
    lower, upper = decoupled_bounds()
    best = None
    process = psutil.Process()
    peak_memory_mb = process.memory_info().rss / (1024.0 * 1024.0)
    start_bics = []

    piece = fit_piecewise_zm(dataset, 500)

    def residuals(params: np.ndarray) -> np.ndarray:
        return decoupled_prediction(ranks, log_ranks, params, gate_name) - y

    for idx, x0 in enumerate(
        make_random_initializations(dataset, gate_name, piece, n_starts, gate_seed(gate_name, base_seed)),
        start=1,
    ):
        result = least_squares(
            residuals,
            x0=x0,
            bounds=(lower, upper),
            method="trf",
            max_nfev=max_nfev,
        )
        pred = decoupled_prediction(ranks, log_ranks, result.x, gate_name)
        score = rmse(y, pred)
        bic = bic_from_rmse(score, cfg.PARAM_COUNT_DECOUPLED, int(dataset["unique_words"]))
        start_bics.append(float(bic))
        peak_memory_mb = max(peak_memory_mb, process.memory_info().rss / (1024.0 * 1024.0))
        if best is None or score < best["rmse"]:
            best = {
                "rmse": float(score),
                "params": [float(v) for v in result.x],
                "start_index": idx,
                "nfev": int(result.nfev),
            }

    params = summarize_params(np.asarray(best["params"], dtype=np.float64), int(dataset["unique_words"]))
    best_bic = float(bic_from_rmse(best["rmse"], cfg.PARAM_COUNT_DECOUPLED, int(dataset["unique_words"])))
    starts_within_1_bic = sum((bic - best_bic) < 1.0 for bic in start_bics)
    return {
        "rmse": float(best["rmse"]),
        "bic": best_bic,
        **params,
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "k_hit_lower_bound": bool(abs(float(params["k"]) - cfg.K_LOWER) <= cfg.BOUND_TOL),
        "k_hit_upper_bound": bool(abs(float(params["k"]) - cfg.K_UPPER) <= cfg.BOUND_TOL),
        "w_gate_hit_lower_bound": bool(abs(float(params["w_gate"]) - cfg.W_GATE_LOWER) <= cfg.BOUND_TOL),
        "w_gate_hit_upper_bound": bool(abs(float(params["w_gate"]) - cfg.W_GATE_UPPER) <= cfg.BOUND_TOL),
        "w_tail_hit_lower_bound": bool(abs(float(params["w_tail"]) - cfg.W_TAIL_LOWER) <= cfg.BOUND_TOL),
        "w_tail_hit_upper_bound": bool(abs(float(params["w_tail"]) - cfg.W_TAIL_UPPER) <= cfg.BOUND_TOL),
        "min_bic_across_starts": float(min(start_bics)),
        "median_bic_across_starts": float(statistics.median(start_bics)),
        "max_bic_across_starts": float(max(start_bics)),
        "starts_within_1_bic": int(starts_within_1_bic),
        "peak_memory_mb": float(peak_memory_mb),
    }

