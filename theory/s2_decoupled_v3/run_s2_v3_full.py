from __future__ import annotations

import csv
import math
import re
import statistics
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import special
from scipy.optimize import least_squares


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "s2_v3_decoupled_five_gate_2026-04-18"
DATA_DIR = ROOT / "data" / "zipf"
THREEA_PATH = ROOT / "experiments" / "3a_smooth_two_regime_fits" / "outputs" / "smooth_fit_per_corpus.csv"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = [
    "*** START OF THE PROJECT GUTENBERG EBOOK",
    "*** START OF THIS PROJECT GUTENBERG EBOOK",
]
END_MARKERS = [
    "*** END OF THE PROJECT GUTENBERG EBOOK",
    "*** END OF THIS PROJECT GUTENBERG EBOOK",
]

SEARCHED_CORPORA = [
    {"slug": "shakespeare", "name": "Complete Works of Shakespeare", "filename": "pg100.txt"},
    {"slug": "war_and_peace", "name": "War and Peace", "filename": "pg2600.txt"},
    {"slug": "moby_dick", "name": "Moby Dick", "filename": "pg2701.txt"},
    {"slug": "king_james_bible", "name": "King James Bible", "filename": "pg10.txt"},
    {"slug": "federalist_papers", "name": "Federalist Papers", "filename": "pg1404.txt"},
    {"slug": "grimms_fairy_tales", "name": "Grimm's Fairy Tales", "filename": "pg2591.txt"},
    {"slug": "don_quixote", "name": "Don Quixote", "filename": "pg996.txt"},
    {"slug": "pride_and_prejudice", "name": "Pride and Prejudice", "filename": "pg1342.txt"},
    {"slug": "canterbury_tales", "name": "Canterbury Tales", "filename": "pg2383.txt"},
    {"slug": "arabian_nights_vol1", "name": "Arabian Nights (Vol 1)", "filename": "pg3435.txt"},
    {"slug": "aesops_fables", "name": "Aesop's Fables", "filename": "pg11339.txt"},
    {"slug": "complete_sherlock_holmes", "name": "Complete Sherlock Holmes", "filename": "pg1661.txt"},
    {"slug": "jane_eyre", "name": "Jane Eyre", "filename": "pg1260.txt"},
    {"slug": "dubliners", "name": "Dubliners", "filename": "pg2814.txt"},
    {"slug": "the_iliad", "name": "The Iliad", "filename": "pg6130.txt"},
    {"slug": "democracy_in_america", "name": "Democracy in America", "filename": "pg815.txt"},
    {"slug": "origin_of_species", "name": "Origin of Species", "filename": "pg1228.txt"},
    {"slug": "wealth_of_nations", "name": "Wealth of Nations", "filename": "pg3300.txt"},
    {"slug": "les_miserables", "name": "Les Miserables", "filename": "pg135.txt"},
    {"slug": "decline_and_fall_vol1", "name": "Decline and Fall Vol 1", "filename": "pg731.txt"},
    {"slug": "emile", "name": "Emile", "filename": "pg5427.txt"},
    {"slug": "ulysses", "name": "Ulysses", "filename": "pg4300.txt"},
    {"slug": "collected_poe", "name": "Collected Poe", "filename": "pg2147.txt"},
    {"slug": "principia_ethica", "name": "Principia Ethica", "filename": "pg53430.txt"},
    {"slug": "critique_of_pure_reason", "name": "Critique of Pure Reason", "filename": "pg4280.txt"},
]
SEARCHED_SLUGS = [spec["slug"] for spec in SEARCHED_CORPORA]
CORPUS_NAME_BY_SLUG = {spec["slug"]: spec["name"] for spec in SEARCHED_CORPORA}

N_RANDOM_STARTS = 100
MAX_NFEV = 12000
BASE_SEED = 20260415

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

LOW_C_CLUSTER = {"grimms_fairy_tales", "aesops_fables", "dubliners", "critique_of_pure_reason"}
SHARED_PARAMS = ["a1", "b1", "c1", "a2", "b2", "c2", "k", "w_gate"]
HEAD_PARAMS = ["a1", "b1", "c1", "k"]
TAIL_PARAMS = ["a2", "b2", "c2"]


def strip_gutenberg_boilerplate(text: str) -> str:
    start = 0
    end = len(text)
    for marker in START_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            line_end = text.find("\n", idx)
            start = line_end + 1 if line_end != -1 else idx
            break
    for marker in END_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            end = idx
            break
    return text[start:end]


def tokenize_text(text: str) -> list[str]:
    return TOKEN_RE.findall(text.lower())


def build_zipf_dataset(corpus_path: Path) -> dict:
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = strip_gutenberg_boilerplate(raw_text)
    counts = Counter(tokenize_text(clean_text))
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    return {
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": np.log(ranks),
        "log_freq": np.log(freqs),
        "token_count": int(sum(counts.values())),
        "unique_words": len(freqs),
    }


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


def gate_width_scale(gate_name: str) -> float:
    return 2.0 if gate_name == "tanh" else 1.0


def logistic_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = np.clip((log_ranks - math.log(k)) / w_gate, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def tanh_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = (log_ranks - math.log(k)) / w_gate
    return 0.5 * (1.0 - np.tanh(z))


def erf_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = (log_ranks - math.log(k)) / w_gate
    return 0.5 * (1.0 - special.erf(z))


def algebraic_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = log_ranks - math.log(k)
    return 0.5 * (1.0 - z / np.sqrt((w_gate * w_gate) + (z * z)))


def arctan_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = log_ranks - math.log(k)
    return 0.5 - np.arctan(z / w_gate) / math.pi


GATE_FUNCS = {
    "logistic": logistic_gate,
    "tanh": tanh_gate,
    "erf": erf_gate,
    "algebraic": algebraic_gate,
    "arctan": arctan_gate,
}
INDEPENDENT_GATES = ["logistic", "erf", "algebraic", "arctan"]


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


def gate_seed(gate_name: str) -> int:
    return BASE_SEED + {"logistic": 0, "tanh": 1, "erf": 2, "algebraic": 3, "arctan": 4}[gate_name]


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


def fit_with_gate(dataset: dict, gate_name: str, n_starts: int = N_RANDOM_STARTS, max_nfev: int = MAX_NFEV) -> dict:
    ranks = dataset["ranks"]
    log_ranks = dataset["log_rank"]
    y = dataset["log_freq"]
    lower, upper = decoupled_bounds()
    best = None

    piece = fit_piecewise_zm(dataset, 500)

    def residuals(params: np.ndarray) -> np.ndarray:
        return decoupled_prediction(ranks, log_ranks, params, gate_name) - y

    for idx, x0 in enumerate(
        make_random_initializations(dataset, gate_name, piece, n_starts, gate_seed(gate_name)), start=1
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
        if best is None or score < best["rmse"]:
            best = {
                "rmse": float(score),
                "params": [float(v) for v in result.x],
                "start_index": idx,
                "nfev": int(result.nfev),
            }

    params = summarize_params(np.asarray(best["params"], dtype=np.float64), int(dataset["unique_words"]))
    return {
        "rmse": float(best["rmse"]),
        "bic": float(bic_from_rmse(best["rmse"], PARAM_COUNT_DECOUPLED, int(dataset["unique_words"]))),
        **params,
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "k_hit_lower_bound": bool(abs(float(params["k"]) - K_LOWER) <= BOUND_TOL),
        "k_hit_upper_bound": bool(abs(float(params["k"]) - K_UPPER) <= BOUND_TOL),
        "w_gate_hit_lower_bound": bool(abs(float(params["w_gate"]) - W_GATE_LOWER) <= BOUND_TOL),
        "w_gate_hit_upper_bound": bool(abs(float(params["w_gate"]) - W_GATE_UPPER) <= BOUND_TOL),
        "w_tail_hit_lower_bound": bool(abs(float(params["w_tail"]) - W_TAIL_LOWER) <= BOUND_TOL),
        "w_tail_hit_upper_bound": bool(abs(float(params["w_tail"]) - W_TAIL_UPPER) <= BOUND_TOL),
    }


def ordered_slug_rows(rows_by_slug: dict[str, dict]) -> list[dict]:
    return [rows_by_slug[slug] for slug in SEARCHED_SLUGS if slug in rows_by_slug]


def load_csv_rows(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["slug"]: row for row in csv.DictReader(handle)}


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_threea_rows() -> dict[str, dict]:
    with THREEA_PATH.open(newline="", encoding="utf-8") as handle:
        return {row["slug"]: row for row in csv.DictReader(handle)}


def percent_change(new_value: float, old_value: float) -> float:
    if old_value == 0.0:
        return float("nan")
    return 100.0 * (new_value - old_value) / old_value


def build_parameter_shift_rows(params_by_slug: dict[str, dict], threea_by_slug: dict[str, dict]) -> list[dict]:
    corpus_rows = []
    shifts_by_param = {param: [] for param in SHARED_PARAMS}

    for slug in SEARCHED_SLUGS:
        if slug not in params_by_slug:
            continue
        p = params_by_slug[slug]
        threea = threea_by_slug[slug]
        row = {
            "row_type": "corpus",
            "label": "",
            "slug": slug,
            "corpus": CORPUS_NAME_BY_SLUG[slug],
        }
        for param in SHARED_PARAMS:
            old_value = float(threea["w"] if param == "w_gate" else threea[param])
            new_value = float(p[f"logistic_{param}"])
            pct = percent_change(new_value, old_value)
            shifts_by_param[param].append(pct)
            row[f"{param}_pct_change"] = repr(float(pct))
        head_abs = statistics.fmean(abs(float(row[f"{param}_pct_change"])) for param in HEAD_PARAMS)
        tail_abs = statistics.fmean(abs(float(row[f"{param}_pct_change"])) for param in TAIL_PARAMS)
        row["head_mean_abs_pct_change"] = repr(float(head_abs))
        row["tail_mean_abs_pct_change"] = repr(float(tail_abs))
        row["notes"] = ""
        corpus_rows.append(row)

    if not corpus_rows:
        return []

    aggregate_rows = []
    for label, reducer in [
        ("mean_percent_change", statistics.fmean),
        ("median_percent_change", statistics.median),
        ("max_abs_percent_change", lambda vals: max(abs(v) for v in vals)),
    ]:
        row = {
            "row_type": "aggregate",
            "label": label,
            "slug": "",
            "corpus": "",
        }
        for param in SHARED_PARAMS:
            row[f"{param}_pct_change"] = repr(float(reducer(shifts_by_param[param])))
        row["head_mean_abs_pct_change"] = ""
        row["tail_mean_abs_pct_change"] = ""
        row["notes"] = "Across completed corpora."
        aggregate_rows.append(row)

    head_abs_per_corpus = [float(row["head_mean_abs_pct_change"]) for row in corpus_rows]
    tail_abs_per_corpus = [float(row["tail_mean_abs_pct_change"]) for row in corpus_rows]
    aggregate_rows.append(
        {
            "row_type": "aggregate",
            "label": "mean_abs_head_vs_tail",
            "slug": "",
            "corpus": "",
            **{f"{param}_pct_change": "" for param in SHARED_PARAMS},
            "head_mean_abs_pct_change": repr(float(statistics.fmean(head_abs_per_corpus))),
            "tail_mean_abs_pct_change": repr(float(statistics.fmean(tail_abs_per_corpus))),
            "notes": "Head group = a1, b1, c1, k. Tail group = a2, b2, c2.",
        }
    )
    return corpus_rows + aggregate_rows


def build_compare_rows(results_by_slug: dict[str, dict], threea_by_slug: dict[str, dict]) -> list[dict]:
    rows = []
    for slug in SEARCHED_SLUGS:
        if slug not in results_by_slug:
            continue
        result = results_by_slug[slug]
        threea = threea_by_slug[slug]
        vocab_size = int(threea["vocab_size"])
        existing_bic = bic_from_rmse(float(threea["reranked_rmse"]), PARAM_COUNT_EXISTING_3A, vocab_size)
        rows.append(
            {
                "slug": slug,
                "corpus": CORPUS_NAME_BY_SLUG[slug],
                "existing_3a_coupled_bic": repr(float(existing_bic)),
                "existing_3a_coupled_rmse": threea["reranked_rmse"],
                "decoupled_logistic_bic": result["logistic_bic"],
                "decoupled_logistic_rmse": result["logistic_rmse"],
                "delta_decoupled_minus_existing": repr(float(result["logistic_bic"]) - float(existing_bic)),
            }
        )
    return rows


def build_aggregate_rows(results_by_slug: dict[str, dict], params_by_slug: dict[str, dict], compare_rows: list[dict]) -> list[dict]:
    completed = ordered_slug_rows(results_by_slug)
    if not completed:
        return []

    independent_winner_counts = {gate: 0 for gate in INDEPENDENT_GATES}
    bic_spreads = []
    tanh_bic_diffs = []
    tanh_gate_ratio_devs = []
    tanh_tail_ratio_devs = []
    strict_count = 0
    positive_count = 0
    strong_count = 0

    gate_param_values = {
        gate: {"k": [], "w_gate": [], "w_tail": []}
        for gate in GATE_FUNCS
    }
    gate_bound_counts = {
        gate: {
            "k_hit_lower_bound_count": 0,
            "k_hit_upper_bound_count": 0,
            "w_gate_hit_lower_bound_count": 0,
            "w_gate_hit_upper_bound_count": 0,
            "w_tail_hit_lower_bound_count": 0,
            "w_tail_hit_upper_bound_count": 0,
        }
        for gate in GATE_FUNCS
    }

    for row in completed:
        tanh_bic_diffs.append(abs(float(row["bic_tanh_minus_logistic"])))
        tanh_gate_ratio_devs.append(abs(float(row["w_gate_tanh_over_2w_gate_logistic_ratio"]) - 1.0))
        tanh_tail_ratio_devs.append(abs(float(row["w_tail_tanh_over_w_tail_logistic_ratio"]) - 1.0))
        independent_bics = {gate: float(row[f"{gate}_bic"]) for gate in INDEPENDENT_GATES}
        winner_gate = min(independent_bics, key=independent_bics.get)
        independent_winner_counts[winner_gate] += 1
        bic_spread = max(independent_bics.values()) - min(independent_bics.values())
        bic_spreads.append(bic_spread)
        if bic_spread < 2.0:
            strict_count += 1
        if bic_spread < 6.0:
            positive_count += 1
        if bic_spread < 10.0:
            strong_count += 1

    for slug in SEARCHED_SLUGS:
        if slug not in params_by_slug:
            continue
        p = params_by_slug[slug]
        for gate in GATE_FUNCS:
            for name in ["k", "w_gate", "w_tail"]:
                gate_param_values[gate][name].append(float(p[f"{gate}_{name}"]))
            for name in gate_bound_counts[gate]:
                gate_bound_counts[gate][name] += int(p[f"{gate}_{name.replace('_count', '')}"] == "True")

    compare_deltas = [float(row["delta_decoupled_minus_existing"]) for row in compare_rows]

    rows = [
        {
            "metric_name": "tanh_calibration_pass_count",
            "value": str(sum(row["tanh_calibration_pass"] == "True" for row in completed)),
            "notes": "Count of corpora satisfying the decoupled tanh calibration checks.",
        },
        {
            "metric_name": "max_tanh_logistic_bic_diff",
            "value": repr(float(max(tanh_bic_diffs))),
            "notes": "Maximum absolute BIC difference between tanh and logistic.",
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
    ]
    for gate in INDEPENDENT_GATES:
        rows.append(
            {
                "metric_name": f"{gate}_bic_wins",
                "value": str(independent_winner_counts[gate]),
                "notes": f"Count of corpora where {gate} has the lowest BIC among the four independent gates.",
            }
        )
    rows.extend(
        [
            {
                "metric_name": "median_bic_spread",
                "value": repr(float(statistics.median(bic_spreads))),
                "notes": "Median of max_bic - min_bic across logistic, erf, algebraic, and arctan.",
            },
            {
                "metric_name": "mean_bic_spread",
                "value": repr(float(statistics.fmean(bic_spreads))),
                "notes": "Mean of max_bic - min_bic across logistic, erf, algebraic, and arctan.",
            },
            {
                "metric_name": "gates_indistinguishable_strong_count",
                "value": str(strong_count),
                "notes": "Count of corpora with independent-gate BIC spread < 10.",
            },
            {
                "metric_name": "gates_indistinguishable_positive_count",
                "value": str(positive_count),
                "notes": "Count of corpora with independent-gate BIC spread < 6.",
            },
            {
                "metric_name": "gates_indistinguishable_strict_count",
                "value": str(strict_count),
                "notes": "Count of corpora with independent-gate BIC spread < 2.",
            },
        ]
    )
    if compare_deltas:
        rows.extend(
            [
                {
                    "metric_name": "median_delta_decoupled_minus_existing_3a_bic",
                    "value": repr(float(statistics.median(compare_deltas))),
                    "notes": "Median BIC delta comparing decoupled logistic to coupled 3a logistic; negative favors decoupling.",
                },
                {
                    "metric_name": "mean_delta_decoupled_minus_existing_3a_bic",
                    "value": repr(float(statistics.fmean(compare_deltas))),
                    "notes": "Mean BIC delta comparing decoupled logistic to coupled 3a logistic; negative favors decoupling.",
                },
                {
                    "metric_name": "decoupled_logistic_beats_existing_3a_count",
                    "value": str(sum(delta < 0.0 for delta in compare_deltas)),
                    "notes": "Count of corpora where decoupled logistic has lower BIC than coupled 3a logistic.",
                },
            ]
        )
    for gate in GATE_FUNCS:
        for param in ["k", "w_gate", "w_tail"]:
            vals = gate_param_values[gate][param]
            rows.extend(
                [
                    {
                        "metric_name": f"{gate}_{param}_min",
                        "value": repr(float(min(vals))),
                        "notes": f"Minimum fitted {param} for {gate}.",
                    },
                    {
                        "metric_name": f"{gate}_{param}_max",
                        "value": repr(float(max(vals))),
                        "notes": f"Maximum fitted {param} for {gate}.",
                    },
                    {
                        "metric_name": f"{gate}_{param}_median",
                        "value": repr(float(statistics.median(vals))),
                        "notes": f"Median fitted {param} for {gate}.",
                    },
                ]
            )
        for name, count in gate_bound_counts[gate].items():
            rows.append(
                {
                    "metric_name": f"{gate}_{name}",
                    "value": str(count),
                    "notes": f"Count of corpora where {gate} hit the {name.replace('_count', '')} flag.",
                }
            )
    return rows


def build_summary_report(results_by_slug: dict[str, dict], aggregate_rows: list[dict], parameter_shift_rows: list[dict]) -> str:
    completed = ordered_slug_rows(results_by_slug)
    if not completed:
        return "# S2 v3 Decoupled Five-Gate Sweep\n\nNo completed corpora yet.\n"

    agg = {row["metric_name"]: row["value"] for row in aggregate_rows}
    coupled_winners = [
        row["corpus"]
        for row in parameter_shift_rows
        if row["row_type"] == "corpus" and float(results_by_slug[row["slug"]]["delta_decoupled_minus_existing"]) >= 0.0
    ]
    lines = [
        "# S2 v3 Decoupled Five-Gate Sweep",
        "",
        f"- completed corpora: `{len(completed)}/25`",
        f"- tanh calibration pass count: `{agg.get('tanh_calibration_pass_count', 'n/a')}`",
        f"- logistic wins: `{agg.get('logistic_bic_wins', 'n/a')}`",
        f"- erf wins: `{agg.get('erf_bic_wins', 'n/a')}`",
        f"- algebraic wins: `{agg.get('algebraic_bic_wins', 'n/a')}`",
        f"- arctan wins: `{agg.get('arctan_bic_wins', 'n/a')}`",
        f"- median independent-gate BIC spread: `{agg.get('median_bic_spread', 'n/a')}`",
        f"- decoupled logistic beats existing 3a count: `{agg.get('decoupled_logistic_beats_existing_3a_count', 'n/a')}`",
    ]
    if coupled_winners:
        lines.extend(["", "## Coupled 3a still wins", "", *[f"- {name}" for name in coupled_winners]])
    return "\n".join(lines) + "\n"


def load_complete_rows(existing_rows: dict[str, dict]) -> dict[str, dict]:
    out = {}
    required = ["logistic_bic", "tanh_bic", "erf_bic", "algebraic_bic", "arctan_bic"]
    for slug, row in existing_rows.items():
        if all(row.get(key) not in ("", None) for key in required):
            out[slug] = row
    return out


def rewrite_outputs(
    results_by_slug: dict[str, dict],
    params_by_slug: dict[str, dict],
    threea_by_slug: dict[str, dict],
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    compare_rows = build_compare_rows(results_by_slug, threea_by_slug)
    parameter_shift_rows = build_parameter_shift_rows(params_by_slug, threea_by_slug)
    aggregate_rows = build_aggregate_rows(results_by_slug, params_by_slug, compare_rows)

    per_corpus_rows = ordered_slug_rows(results_by_slug)
    param_rows = ordered_slug_rows(params_by_slug)

    per_corpus_fieldnames = [
        "slug", "corpus", "vocab_size",
    ]
    for gate in GATE_FUNCS:
        per_corpus_fieldnames.extend(
            [
                f"{gate}_bic",
                f"{gate}_rmse",
                f"{gate}_k",
                f"{gate}_w_gate",
                f"{gate}_w_tail",
                f"{gate}_k_hit_lower_bound",
                f"{gate}_k_hit_upper_bound",
                f"{gate}_w_gate_hit_lower_bound",
                f"{gate}_w_gate_hit_upper_bound",
                f"{gate}_w_tail_hit_lower_bound",
                f"{gate}_w_tail_hit_upper_bound",
            ]
        )
    per_corpus_fieldnames.extend(
        [
            "bic_tanh_minus_logistic",
            "w_gate_tanh_over_2w_gate_logistic_ratio",
            "w_tail_tanh_over_w_tail_logistic_ratio",
            "tanh_calibration_pass",
            "winner_gate",
            "worst_gate",
            "bic_spread",
            "delta_decoupled_minus_existing",
        ]
    )

    gate_param_fieldnames = ["slug", "corpus"]
    for gate in GATE_FUNCS:
        gate_param_fieldnames.extend(
            [
                f"{gate}_a1",
                f"{gate}_b1",
                f"{gate}_c1",
                f"{gate}_a2",
                f"{gate}_b2",
                f"{gate}_c2",
                f"{gate}_k",
                f"{gate}_w_gate",
                f"{gate}_w_tail",
                f"{gate}_transition_fraction",
                f"{gate}_best_start_index",
                f"{gate}_best_nfev",
                f"{gate}_k_hit_lower_bound",
                f"{gate}_k_hit_upper_bound",
                f"{gate}_w_gate_hit_lower_bound",
                f"{gate}_w_gate_hit_upper_bound",
                f"{gate}_w_tail_hit_lower_bound",
                f"{gate}_w_tail_hit_upper_bound",
            ]
        )

    parameter_shift_fieldnames = [
        "row_type",
        "label",
        "slug",
        "corpus",
        *[f"{param}_pct_change" for param in SHARED_PARAMS],
        "head_mean_abs_pct_change",
        "tail_mean_abs_pct_change",
        "notes",
    ]

    write_csv(OUTDIR / "s2_v3_per_corpus_results.csv", per_corpus_rows, per_corpus_fieldnames)
    write_csv(OUTDIR / "s2_v3_gate_params_per_corpus.csv", param_rows, gate_param_fieldnames)
    write_csv(
        OUTDIR / "s2_v3_vs_existing_3a_per_corpus.csv",
        compare_rows,
        [
            "slug",
            "corpus",
            "existing_3a_coupled_bic",
            "existing_3a_coupled_rmse",
            "decoupled_logistic_bic",
            "decoupled_logistic_rmse",
            "delta_decoupled_minus_existing",
        ],
    )
    write_csv(OUTDIR / "s2_v3_parameter_shifts.csv", parameter_shift_rows, parameter_shift_fieldnames)
    write_csv(OUTDIR / "s2_v3_aggregate_statistics.csv", aggregate_rows, ["metric_name", "value", "notes"])
    (OUTDIR / "s2_v3_summary_report.md").write_text(
        build_summary_report(results_by_slug, aggregate_rows, parameter_shift_rows), encoding="utf-8"
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    progress_path = OUTDIR / "progress.log"
    if not progress_path.exists():
        progress_path.write_text("", encoding="utf-8")

    results_by_slug = load_complete_rows(load_csv_rows(OUTDIR / "s2_v3_per_corpus_results.csv"))
    params_by_slug = load_csv_rows(OUTDIR / "s2_v3_gate_params_per_corpus.csv")
    threea_by_slug = load_threea_rows()

    completed = set(results_by_slug) & set(params_by_slug)
    total = len(SEARCHED_CORPORA)

    for idx, spec in enumerate(SEARCHED_CORPORA, start=1):
        slug = spec["slug"]
        if slug in completed:
            continue
        print(f"[{idx}/{total}] decoupled five-gate fit {slug} ...", flush=True)
        dataset = build_zipf_dataset(DATA_DIR / spec["filename"])
        fits = {gate: fit_with_gate(dataset, gate) for gate in GATE_FUNCS}

        tanh_bic_diff = abs(float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"]))
        tanh_gate_ratio = float(fits["tanh"]["w_gate"]) / (2.0 * float(fits["logistic"]["w_gate"]))
        tanh_tail_ratio = float(fits["tanh"]["w_tail"]) / float(fits["logistic"]["w_tail"])
        tanh_pass = tanh_bic_diff < TANH_BIC_THRESHOLD and abs(tanh_gate_ratio - 1.0) < TANH_GATE_RATIO_THRESHOLD

        independent_bics = {gate: float(fits[gate]["bic"]) for gate in INDEPENDENT_GATES}
        winner_gate = min(independent_bics, key=independent_bics.get)
        worst_gate = max(independent_bics, key=independent_bics.get)
        bic_spread = max(independent_bics.values()) - min(independent_bics.values())

        vocab_size = int(dataset["unique_words"])
        existing_bic = bic_from_rmse(
            float(threea_by_slug[slug]["reranked_rmse"]), PARAM_COUNT_EXISTING_3A, int(threea_by_slug[slug]["vocab_size"])
        )
        results_row = {
            "slug": slug,
            "corpus": spec["name"],
            "vocab_size": str(vocab_size),
        }
        params_row = {"slug": slug, "corpus": spec["name"]}
        for gate in GATE_FUNCS:
            for key in [
                "bic", "rmse", "k", "w_gate", "w_tail",
                "k_hit_lower_bound", "k_hit_upper_bound",
                "w_gate_hit_lower_bound", "w_gate_hit_upper_bound",
                "w_tail_hit_lower_bound", "w_tail_hit_upper_bound",
            ]:
                results_row[f"{gate}_{key}"] = repr(float(fits[gate][key])) if isinstance(fits[gate][key], float) else str(fits[gate][key])
            for key in [
                "a1", "b1", "c1", "a2", "b2", "c2", "k", "w_gate", "w_tail",
                "transition_fraction", "best_start_index", "best_nfev",
                "k_hit_lower_bound", "k_hit_upper_bound",
                "w_gate_hit_lower_bound", "w_gate_hit_upper_bound",
                "w_tail_hit_lower_bound", "w_tail_hit_upper_bound",
            ]:
                val = fits[gate][key]
                params_row[f"{gate}_{key}"] = repr(float(val)) if isinstance(val, float) else str(val)
        results_row["bic_tanh_minus_logistic"] = repr(float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"]))
        results_row["w_gate_tanh_over_2w_gate_logistic_ratio"] = repr(float(tanh_gate_ratio))
        results_row["w_tail_tanh_over_w_tail_logistic_ratio"] = repr(float(tanh_tail_ratio))
        results_row["tanh_calibration_pass"] = str(tanh_pass)
        results_row["winner_gate"] = winner_gate
        results_row["worst_gate"] = worst_gate
        results_row["bic_spread"] = repr(float(bic_spread))
        results_row["delta_decoupled_minus_existing"] = repr(float(fits["logistic"]["bic"]) - float(existing_bic))

        results_by_slug[slug] = results_row
        params_by_slug[slug] = params_row
        rewrite_outputs(results_by_slug, params_by_slug, threea_by_slug)
        with progress_path.open("a", encoding="utf-8") as handle:
            handle.write(f"completed {idx}/{total}: {slug}\n")


if __name__ == "__main__":
    main()
