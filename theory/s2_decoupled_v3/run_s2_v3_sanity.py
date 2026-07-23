from __future__ import annotations

import csv
import math
import statistics
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "s2_decoupled_v3"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import zipf_analysis_common as common
import zipf_correct_model as correct_model


N_RANDOM_STARTS = 100
MAX_NFEV = 12000
PARAM_COUNT_DECOUPLED = 9
PARAM_COUNT_EXISTING_3A = 8
K_LOWER = 20.0
K_UPPER = 1000.0
W_GATE_LOWER = 0.05
W_GATE_UPPER = 10.0
W_TAIL_LOWER = 0.05
W_TAIL_UPPER = 10.0
BOUND_TOL = 1e-6
TANH_BIC_THRESHOLD = 1.0
TANH_GATE_RATIO_THRESHOLD = 0.05


def bic_from_rmse(rmse: float, p: int, n: int) -> float:
    mse = float(rmse) ** 2
    return p * math.log(n) + n * math.log(mse)


def gate_width_scale(gate_name: str) -> float:
    return 2.0 if gate_name == "tanh" else 1.0


def logistic_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = np.clip((log_ranks - math.log(k)) / w_gate, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def tanh_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = (log_ranks - math.log(k)) / w_gate
    return 0.5 * (1.0 - np.tanh(z))


GATE_FUNCS = {
    "logistic": logistic_gate,
    "tanh": tanh_gate,
}


def decoupled_bounds() -> tuple[np.ndarray, np.ndarray]:
    lower = np.array(
        [-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, K_LOWER, W_GATE_LOWER, W_TAIL_LOWER],
        dtype=np.float64,
    )
    upper = np.array(
        [100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, K_UPPER, W_GATE_UPPER, W_TAIL_UPPER],
        dtype=np.float64,
    )
    return lower, upper


def summarize_params(params: np.ndarray, vocab_size: int) -> dict[str, float]:
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


def make_random_initializations(dataset: dict, gate_name: str, n_starts: int = N_RANDOM_STARTS, seed: int = 20260415) -> list[np.ndarray]:
    lower, upper = decoupled_bounds()
    rng = np.random.default_rng(seed)
    y = dataset["log_freq"]
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    piece = correct_model.fit_piecewise_zm(dataset, 500)
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

    starts: list[np.ndarray] = []
    shape_gate_lower = max(W_GATE_LOWER / gate_scale, 1e-6)
    shape_gate_upper = W_GATE_UPPER / gate_scale
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


def fit_with_gate(dataset: dict, gate_name: str, n_starts: int = N_RANDOM_STARTS, max_nfev: int = MAX_NFEV) -> dict[str, float | int | bool]:
    ranks = dataset["ranks"]
    log_ranks = dataset["log_rank"]
    y = dataset["log_freq"]
    lower, upper = decoupled_bounds()
    best = None
    tried = []

    def residuals(params: np.ndarray) -> np.ndarray:
        return decoupled_prediction(ranks, log_ranks, params, gate_name) - y

    for idx, x0 in enumerate(make_random_initializations(dataset, gate_name, n_starts=n_starts), start=1):
        result = least_squares(
            residuals,
            x0=x0,
            bounds=(lower, upper),
            method="trf",
            max_nfev=max_nfev,
        )
        pred = decoupled_prediction(ranks, log_ranks, result.x, gate_name)
        score = correct_model.rmse(y, pred)
        record = {
            "start_index": idx,
            "rmse": float(score),
            "success": bool(result.success),
            "nfev": int(result.nfev),
            "status": int(result.status),
            "message": result.message,
            "params": [float(v) for v in result.x],
        }
        tried.append(record)
        if best is None or score < best["rmse"]:
            best = {
                "rmse": float(score),
                "prediction": pred,
                "params": [float(v) for v in result.x],
                "start_index": idx,
                "nfev": int(result.nfev),
                "success": bool(result.success),
                "status": int(result.status),
                "message": result.message,
            }

    params = summarize_params(np.asarray(best["params"], dtype=np.float64), int(dataset["unique_words"]))
    return {
        "rmse": float(best["rmse"]),
        "bic": float(bic_from_rmse(best["rmse"], PARAM_COUNT_DECOUPLED, int(dataset["unique_words"]))),
        "k": float(params["k"]),
        "w_gate": float(params["w_gate"]),
        "w_tail": float(params["w_tail"]),
        "transition_fraction": float(params["transition_fraction"]),
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "tried_count": len(tried),
        "k_hit_lower_bound": bool(abs(float(params["k"]) - K_LOWER) <= BOUND_TOL),
        "k_hit_upper_bound": bool(abs(float(params["k"]) - K_UPPER) <= BOUND_TOL),
        "w_gate_hit_lower_bound": bool(abs(float(params["w_gate"]) - W_GATE_LOWER) <= BOUND_TOL),
        "w_gate_hit_upper_bound": bool(abs(float(params["w_gate"]) - W_GATE_UPPER) <= BOUND_TOL),
        "w_tail_hit_lower_bound": bool(abs(float(params["w_tail"]) - W_TAIL_LOWER) <= BOUND_TOL),
        "w_tail_hit_upper_bound": bool(abs(float(params["w_tail"]) - W_TAIL_UPPER) <= BOUND_TOL),
    }


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_existing_3a_bics() -> dict[str, dict[str, float]]:
    path = ROOT / "experiments" / "3a_smooth_two_regime_fits" / "outputs" / "smooth_fit_per_corpus.csv"
    out: dict[str, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            vocab_size = int(row["vocab_size"])
            reranked_rmse = float(row["reranked_rmse"])
            out[row["slug"]] = {
                "bic": float(bic_from_rmse(reranked_rmse, PARAM_COUNT_EXISTING_3A, vocab_size)),
                "rmse": reranked_rmse,
                "vocab_size": vocab_size,
            }
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    existing_3a = load_existing_3a_bics()
    corpus_specs = common.SEARCHED_CORPORA

    sanity_rows: list[dict[str, str]] = []
    comparison_rows: list[dict[str, str]] = []
    logistic_vs_3a_deltas: list[float] = []
    tanh_bic_diffs: list[float] = []
    tanh_gate_ratio_devs: list[float] = []
    tanh_tail_ratio_devs: list[float] = []
    tanh_failures: list[str] = []

    total = len(corpus_specs)
    for idx, spec in enumerate(corpus_specs, start=1):
        print(f"[{idx}/{total}] decoupled sanity fit {spec['slug']} ...", flush=True)
        dataset = common.build_zipf_dataset(common.corpus_path(spec))
        fits = {gate: fit_with_gate(dataset, gate) for gate in GATE_FUNCS}

        tanh_bic_diff = abs(float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"]))
        tanh_gate_ratio = float(fits["tanh"]["w_gate"]) / (2.0 * float(fits["logistic"]["w_gate"]))
        tanh_gate_ratio_dev = abs(tanh_gate_ratio - 1.0)
        tanh_tail_ratio = float(fits["tanh"]["w_tail"]) / float(fits["logistic"]["w_tail"])
        tanh_tail_ratio_dev = abs(tanh_tail_ratio - 1.0)
        tanh_pass = tanh_bic_diff < TANH_BIC_THRESHOLD and tanh_gate_ratio_dev < TANH_GATE_RATIO_THRESHOLD
        if not tanh_pass:
            tanh_failures.append(spec["name"])

        existing = existing_3a[spec["slug"]]
        delta_vs_existing = float(fits["logistic"]["bic"]) - float(existing["bic"])
        logistic_vs_3a_deltas.append(delta_vs_existing)
        tanh_bic_diffs.append(tanh_bic_diff)
        tanh_gate_ratio_devs.append(tanh_gate_ratio_dev)
        tanh_tail_ratio_devs.append(tanh_tail_ratio_dev)

        sanity_rows.append(
            {
                "corpus": spec["name"],
                "slug": spec["slug"],
                "vocab_size": str(int(dataset["unique_words"])),
                "logistic_bic": repr(float(fits["logistic"]["bic"])),
                "logistic_rmse": repr(float(fits["logistic"]["rmse"])),
                "logistic_k": repr(float(fits["logistic"]["k"])),
                "logistic_w_gate": repr(float(fits["logistic"]["w_gate"])),
                "logistic_w_tail": repr(float(fits["logistic"]["w_tail"])),
                "logistic_k_hit_lower_bound": str(bool(fits["logistic"]["k_hit_lower_bound"])),
                "logistic_k_hit_upper_bound": str(bool(fits["logistic"]["k_hit_upper_bound"])),
                "logistic_w_gate_hit_lower_bound": str(bool(fits["logistic"]["w_gate_hit_lower_bound"])),
                "logistic_w_gate_hit_upper_bound": str(bool(fits["logistic"]["w_gate_hit_upper_bound"])),
                "logistic_w_tail_hit_lower_bound": str(bool(fits["logistic"]["w_tail_hit_lower_bound"])),
                "logistic_w_tail_hit_upper_bound": str(bool(fits["logistic"]["w_tail_hit_upper_bound"])),
                "tanh_bic": repr(float(fits["tanh"]["bic"])),
                "tanh_rmse": repr(float(fits["tanh"]["rmse"])),
                "tanh_k": repr(float(fits["tanh"]["k"])),
                "tanh_w_gate": repr(float(fits["tanh"]["w_gate"])),
                "tanh_w_tail": repr(float(fits["tanh"]["w_tail"])),
                "tanh_k_hit_lower_bound": str(bool(fits["tanh"]["k_hit_lower_bound"])),
                "tanh_k_hit_upper_bound": str(bool(fits["tanh"]["k_hit_upper_bound"])),
                "tanh_w_gate_hit_lower_bound": str(bool(fits["tanh"]["w_gate_hit_lower_bound"])),
                "tanh_w_gate_hit_upper_bound": str(bool(fits["tanh"]["w_gate_hit_upper_bound"])),
                "tanh_w_tail_hit_lower_bound": str(bool(fits["tanh"]["w_tail_hit_lower_bound"])),
                "tanh_w_tail_hit_upper_bound": str(bool(fits["tanh"]["w_tail_hit_upper_bound"])),
                "bic_tanh_minus_logistic": repr(float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"])),
                "w_gate_tanh_over_2w_gate_logistic_ratio": repr(float(tanh_gate_ratio)),
                "w_tail_tanh_over_w_tail_logistic_ratio": repr(float(tanh_tail_ratio)),
                "tanh_calibration_pass": str(tanh_pass),
            }
        )
        comparison_rows.append(
            {
                "corpus": spec["name"],
                "slug": spec["slug"],
                "existing_3a_coupled_bic": repr(float(existing["bic"])),
                "existing_3a_coupled_rmse": repr(float(existing["rmse"])),
                "decoupled_logistic_bic": repr(float(fits["logistic"]["bic"])),
                "decoupled_logistic_rmse": repr(float(fits["logistic"]["rmse"])),
                "delta_decoupled_minus_existing": repr(float(delta_vs_existing)),
            }
        )

    aggregate_rows = [
        {
            "metric_name": "tanh_calibration_pass_count",
            "value": str(len(corpus_specs) - len(tanh_failures)),
            "notes": "Count of corpora satisfying both decoupled tanh calibration checks.",
        },
        {
            "metric_name": "max_tanh_logistic_bic_diff",
            "value": repr(float(max(tanh_bic_diffs))),
            "notes": "Maximum absolute BIC difference between tanh and logistic under decoupled parameterization.",
        },
        {
            "metric_name": "max_tanh_logistic_w_gate_ratio_deviation",
            "value": repr(float(max(tanh_gate_ratio_devs))),
            "notes": "Maximum absolute deviation of w_gate_tanh / (2*w_gate_logistic) from 1.0.",
        },
        {
            "metric_name": "max_tanh_logistic_w_tail_ratio_deviation",
            "value": repr(float(max(tanh_tail_ratio_devs))),
            "notes": "Maximum absolute deviation of w_tail_tanh / w_tail_logistic from 1.0.",
        },
        {
            "metric_name": "median_delta_decoupled_minus_existing_3a_bic",
            "value": repr(float(statistics.median(logistic_vs_3a_deltas))),
            "notes": "Median BIC delta comparing the decoupled logistic model against the coupled 3a logistic model; negative favors decoupling.",
        },
        {
            "metric_name": "mean_delta_decoupled_minus_existing_3a_bic",
            "value": repr(float(statistics.fmean(logistic_vs_3a_deltas))),
            "notes": "Mean BIC delta comparing the decoupled logistic model against the coupled 3a logistic model; negative favors decoupling.",
        },
        {
            "metric_name": "decoupled_logistic_beats_existing_3a_count",
            "value": str(sum(1 for delta in logistic_vs_3a_deltas if delta < 0.0)),
            "notes": "Count of corpora where the decoupled logistic model has lower BIC than the existing coupled 3a logistic model.",
        },
    ]

    write_csv(
        OUTDIR / "s2_v3_sanity_per_corpus.csv",
        sanity_rows,
        [
            "corpus",
            "slug",
            "vocab_size",
            "logistic_bic",
            "logistic_rmse",
            "logistic_k",
            "logistic_w_gate",
            "logistic_w_tail",
            "logistic_k_hit_lower_bound",
            "logistic_k_hit_upper_bound",
            "logistic_w_gate_hit_lower_bound",
            "logistic_w_gate_hit_upper_bound",
            "logistic_w_tail_hit_lower_bound",
            "logistic_w_tail_hit_upper_bound",
            "tanh_bic",
            "tanh_rmse",
            "tanh_k",
            "tanh_w_gate",
            "tanh_w_tail",
            "tanh_k_hit_lower_bound",
            "tanh_k_hit_upper_bound",
            "tanh_w_gate_hit_lower_bound",
            "tanh_w_gate_hit_upper_bound",
            "tanh_w_tail_hit_lower_bound",
            "tanh_w_tail_hit_upper_bound",
            "bic_tanh_minus_logistic",
            "w_gate_tanh_over_2w_gate_logistic_ratio",
            "w_tail_tanh_over_w_tail_logistic_ratio",
            "tanh_calibration_pass",
        ],
    )
    write_csv(
        OUTDIR / "s2_v3_vs_existing_3a_per_corpus.csv",
        comparison_rows,
        [
            "corpus",
            "slug",
            "existing_3a_coupled_bic",
            "existing_3a_coupled_rmse",
            "decoupled_logistic_bic",
            "decoupled_logistic_rmse",
            "delta_decoupled_minus_existing",
        ],
    )
    write_csv(
        OUTDIR / "s2_v3_aggregate_statistics.csv",
        aggregate_rows,
        ["metric_name", "value", "notes"],
    )

    lines = [
        "# S2 v3 Decoupled Sanity Check",
        "",
        "- This stage reruns only `logistic` and `tanh` with decoupled `w_gate` and `w_tail`.",
        "- Pass criterion: `|BIC_tanh - BIC_logistic| < 1` and `|w_gate_tanh / (2 * w_gate_logistic) - 1| < 0.05` per corpus.",
        "",
        f"- tanh calibration pass count: `{len(corpus_specs) - len(tanh_failures)}/{len(corpus_specs)}`",
        f"- max absolute tanh/logistic BIC diff: `{max(tanh_bic_diffs):.12f}`",
        f"- max absolute tanh/logistic gate-width ratio deviation: `{max(tanh_gate_ratio_devs):.12f}`",
        f"- max absolute tanh/logistic tail-width ratio deviation: `{max(tanh_tail_ratio_devs):.12f}`",
        f"- median decoupled-logistic minus existing-3a BIC delta: `{statistics.median(logistic_vs_3a_deltas):.12f}`",
        f"- mean decoupled-logistic minus existing-3a BIC delta: `{statistics.fmean(logistic_vs_3a_deltas):.12f}`",
        "",
    ]
    if tanh_failures:
        lines.extend(
            [
                "## Calibration failures",
                "",
                *[f"- {name}" for name in tanh_failures],
                "",
                "- Interpretation: decoupling did not fully restore logistic/tanh equivalence under the current optimizer configuration.",
            ]
        )
    else:
        lines.extend(
            [
                "## Calibration outcome",
                "",
                "- Interpretation: decoupling restores the expected logistic/tanh equivalence, so the gate-family sweep can proceed.",
            ]
        )

    (OUTDIR / "s2_v3_summary_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
