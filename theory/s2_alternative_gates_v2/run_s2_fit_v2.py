from __future__ import annotations

import csv
import importlib.util
import json
import math
import statistics
import sys
from pathlib import Path

import numpy as np
from scipy import special


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "s2_alternative_gates_v2"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import zipf_analysis_common as common


RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
N_RANDOM_STARTS = 100
MAX_NFEV = 12000
PARAM_COUNT = 8
W_LOWER = 0.05
W_UPPER = 10.0
BOUND_TOL = 1e-6


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reranked = load_module(RERANKED_PATH, "s2_alt_gate_v2_reranked")


def bic_from_rmse(rmse: float, p: int, n: int) -> float:
    mse = float(rmse) ** 2
    return p * math.log(n) + n * math.log(mse)


def widened_bounds():
    lower = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, W_LOWER], dtype=np.float64)
    upper = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 1000.0, W_UPPER], dtype=np.float64)
    return lower, upper


def logistic_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = np.clip((np.log(ranks) - math.log(k)) / w, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def tanh_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = (np.log(ranks) - math.log(k)) / w
    return 0.5 * (1.0 - np.tanh(z))


def erf_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = (np.log(ranks) - math.log(k)) / w
    return 0.5 * (1.0 - special.erf(z))


def algebraic_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = np.log(ranks) - math.log(k)
    return 0.5 * (1.0 - z / np.sqrt((w * w) + (z * z)))


def arctan_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = np.log(ranks) - math.log(k)
    return 0.5 - np.arctan(z / w) / math.pi


GATE_FUNCS = {
    "logistic": logistic_gate,
    "tanh": tanh_gate,
    "erf": erf_gate,
    "algebraic": algebraic_gate,
    "arctan": arctan_gate,
}


def fit_with_gate(dataset: dict, gate_name: str) -> dict[str, float | int | bool]:
    original_sigma = reranked.sigma_curve
    original_bounds = reranked.constrained_bounds
    reranked.sigma_curve = GATE_FUNCS[gate_name]
    reranked.constrained_bounds = widened_bounds
    try:
        best, tried = reranked.fit_reranked_model(dataset, n_starts=N_RANDOM_STARTS, max_nfev=MAX_NFEV)
    finally:
        reranked.sigma_curve = original_sigma
        reranked.constrained_bounds = original_bounds

    params = reranked.summarize_params(best["params"], dataset["unique_words"])
    w_value = float(params["w"])
    return {
        "rmse": float(best["rmse"]),
        "bic": float(bic_from_rmse(best["rmse"], PARAM_COUNT, int(dataset["unique_words"]))),
        "k": float(params["k"]),
        "w": w_value,
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "tried_count": len(tried),
        "hit_lower_bound": bool(abs(w_value - W_LOWER) <= BOUND_TOL),
        "hit_upper_bound": bool(abs(w_value - W_UPPER) <= BOUND_TOL),
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    corpus_specs = common.SEARCHED_CORPORA

    result_rows: list[dict[str, str]] = []
    param_rows: list[dict[str, str]] = []

    independent_winner_counts = {gate: 0 for gate in ["logistic", "erf", "algebraic", "arctan"]}
    bic_spreads: list[float] = []
    tanh_bic_diffs: list[float] = []
    tanh_ratio_devs: list[float] = []
    tanh_failures: list[str] = []
    gate_w_values = {gate: [] for gate in GATE_FUNCS}
    gate_upper_hits = {gate: 0 for gate in GATE_FUNCS}
    gate_lower_hits = {gate: 0 for gate in GATE_FUNCS}
    strict_count = 0
    positive_count = 0
    strong_count = 0

    total = len(corpus_specs)
    for idx, spec in enumerate(corpus_specs, start=1):
        print(f"[{idx}/{total}] fitting {spec['slug']} ...", flush=True)
        dataset = common.build_zipf_dataset(common.corpus_path(spec))
        fits = {gate: fit_with_gate(dataset, gate) for gate in GATE_FUNCS}

        for gate, fit in fits.items():
            gate_w_values[gate].append(float(fit["w"]))
            gate_upper_hits[gate] += int(bool(fit["hit_upper_bound"]))
            gate_lower_hits[gate] += int(bool(fit["hit_lower_bound"]))

        tanh_bic_diff = abs(float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"]))
        tanh_ratio = float(fits["tanh"]["w"]) / (2.0 * float(fits["logistic"]["w"]))
        tanh_ratio_dev = abs(tanh_ratio - 1.0)
        tanh_pass = tanh_bic_diff < 1.0 and tanh_ratio_dev < 0.05
        if not tanh_pass:
            tanh_failures.append(spec["name"])

        independent_bics = {gate: float(fits[gate]["bic"]) for gate in ["logistic", "erf", "algebraic", "arctan"]}
        winner_gate = min(independent_bics, key=independent_bics.get)
        worst_gate = max(independent_bics, key=independent_bics.get)
        bic_spread = max(independent_bics.values()) - min(independent_bics.values())
        independent_winner_counts[winner_gate] += 1
        bic_spreads.append(bic_spread)
        tanh_bic_diffs.append(tanh_bic_diff)
        tanh_ratio_devs.append(tanh_ratio_dev)
        if bic_spread < 2.0:
            strict_count += 1
        if bic_spread < 6.0:
            positive_count += 1
        if bic_spread < 10.0:
            strong_count += 1

        result_rows.append(
            {
                "corpus": spec["name"],
                "logistic_bic": repr(float(fits["logistic"]["bic"])),
                "logistic_rmse": repr(float(fits["logistic"]["rmse"])),
                "logistic_k": repr(float(fits["logistic"]["k"])),
                "logistic_w": repr(float(fits["logistic"]["w"])),
                "tanh_bic": repr(float(fits["tanh"]["bic"])),
                "tanh_rmse": repr(float(fits["tanh"]["rmse"])),
                "tanh_k": repr(float(fits["tanh"]["k"])),
                "tanh_w": repr(float(fits["tanh"]["w"])),
                "erf_bic": repr(float(fits["erf"]["bic"])),
                "erf_rmse": repr(float(fits["erf"]["rmse"])),
                "erf_k": repr(float(fits["erf"]["k"])),
                "erf_w": repr(float(fits["erf"]["w"])),
                "algebraic_bic": repr(float(fits["algebraic"]["bic"])),
                "algebraic_rmse": repr(float(fits["algebraic"]["rmse"])),
                "algebraic_k": repr(float(fits["algebraic"]["k"])),
                "algebraic_w": repr(float(fits["algebraic"]["w"])),
                "arctan_bic": repr(float(fits["arctan"]["bic"])),
                "arctan_rmse": repr(float(fits["arctan"]["rmse"])),
                "arctan_k": repr(float(fits["arctan"]["k"])),
                "arctan_w": repr(float(fits["arctan"]["w"])),
                "bic_tanh_minus_logistic": repr(float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"])),
                "w_tanh_over_2w_logistic_ratio": repr(tanh_ratio),
                "tanh_calibration_pass": str(tanh_pass),
                "winner_gate": winner_gate,
                "worst_gate": worst_gate,
                "bic_spread": repr(float(bic_spread)),
            }
        )
        param_rows.append(
            {
                "corpus": spec["name"],
                "logistic_k": repr(float(fits["logistic"]["k"])),
                "logistic_w": repr(float(fits["logistic"]["w"])),
                "logistic_hit_upper_bound": str(bool(fits["logistic"]["hit_upper_bound"])),
                "logistic_hit_lower_bound": str(bool(fits["logistic"]["hit_lower_bound"])),
                "tanh_k": repr(float(fits["tanh"]["k"])),
                "tanh_w": repr(float(fits["tanh"]["w"])),
                "tanh_hit_upper_bound": str(bool(fits["tanh"]["hit_upper_bound"])),
                "tanh_hit_lower_bound": str(bool(fits["tanh"]["hit_lower_bound"])),
                "erf_k": repr(float(fits["erf"]["k"])),
                "erf_w": repr(float(fits["erf"]["w"])),
                "erf_hit_upper_bound": str(bool(fits["erf"]["hit_upper_bound"])),
                "erf_hit_lower_bound": str(bool(fits["erf"]["hit_lower_bound"])),
                "algebraic_k": repr(float(fits["algebraic"]["k"])),
                "algebraic_w": repr(float(fits["algebraic"]["w"])),
                "algebraic_hit_upper_bound": str(bool(fits["algebraic"]["hit_upper_bound"])),
                "algebraic_hit_lower_bound": str(bool(fits["algebraic"]["hit_lower_bound"])),
                "arctan_k": repr(float(fits["arctan"]["k"])),
                "arctan_w": repr(float(fits["arctan"]["w"])),
                "arctan_hit_upper_bound": str(bool(fits["arctan"]["hit_upper_bound"])),
                "arctan_hit_lower_bound": str(bool(fits["arctan"]["hit_lower_bound"])),
            }
        )

    aggregate_rows = [
        {
            "metric_name": "logistic_bic_wins",
            "value": str(independent_winner_counts["logistic"]),
            "notes": "Count of corpora where logistic has the lowest BIC among the four independent gates.",
        },
        {
            "metric_name": "erf_bic_wins",
            "value": str(independent_winner_counts["erf"]),
            "notes": "Count of corpora where erf has the lowest BIC among the four independent gates.",
        },
        {
            "metric_name": "algebraic_bic_wins",
            "value": str(independent_winner_counts["algebraic"]),
            "notes": "Count of corpora where algebraic has the lowest BIC among the four independent gates.",
        },
        {
            "metric_name": "arctan_bic_wins",
            "value": str(independent_winner_counts["arctan"]),
            "notes": "Count of corpora where arctan has the lowest BIC among the four independent gates.",
        },
        {
            "metric_name": "tanh_calibration_pass_count",
            "value": str(len(corpus_specs) - len(tanh_failures)),
            "notes": "Count of corpora satisfying both tanh calibration checks.",
        },
        {
            "metric_name": "max_tanh_logistic_bic_diff",
            "value": repr(float(max(tanh_bic_diffs))),
            "notes": "Maximum absolute BIC difference between tanh and logistic across corpora.",
        },
        {
            "metric_name": "max_tanh_logistic_ratio_deviation",
            "value": repr(float(max(tanh_ratio_devs))),
            "notes": "Maximum absolute deviation of w_tanh / (2*w_logistic) from 1.0 across corpora.",
        },
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

    for gate in GATE_FUNCS:
        values = gate_w_values[gate]
        aggregate_rows.extend(
            [
                {
                    "metric_name": f"{gate}_w_min",
                    "value": repr(float(min(values))),
                    "notes": f"Minimum fitted w for {gate}.",
                },
                {
                    "metric_name": f"{gate}_w_max",
                    "value": repr(float(max(values))),
                    "notes": f"Maximum fitted w for {gate}.",
                },
                {
                    "metric_name": f"{gate}_w_median",
                    "value": repr(float(statistics.median(values))),
                    "notes": f"Median fitted w for {gate}.",
                },
                {
                    "metric_name": f"{gate}_hit_upper_bound_count",
                    "value": str(gate_upper_hits[gate]),
                    "notes": f"Count of corpora where {gate} hit the widened upper w bound.",
                },
                {
                    "metric_name": f"{gate}_hit_lower_bound_count",
                    "value": str(gate_lower_hits[gate]),
                    "notes": f"Count of corpora where {gate} hit the widened lower w bound.",
                },
            ]
        )

    lines = [
        "# S2 v2 Alternative Gate Fit Sweep",
        "",
        "- Sanity check summary remains in `s2_pre_fit_sanity_check_v2.csv`.",
        "- Full fit sweep uses five gates: logistic, tanh, erf, algebraic, arctan.",
        "- Bounds widened only in `w`: `[0.05, 10.0]`; all other constraints and optimizer settings match the historical reranked smooth-model fit.",
        "",
        f"- tanh calibration pass count: `{len(corpus_specs) - len(tanh_failures)}` / `{len(corpus_specs)}`",
        f"- max |BIC_tanh - BIC_logistic|: `{max(tanh_bic_diffs):.12f}`",
        f"- max |w_tanh/(2*w_logistic) - 1|: `{max(tanh_ratio_devs):.12f}`",
        "",
    ]

    if tanh_failures:
        lines.append("## Tanh calibration failures")
        lines.append("")
        for name in tanh_failures:
            lines.append(f"- `{name}`")
        lines.append("")

    if len(tanh_failures) > 2:
        lines.append("- Interpretation: more than 2 corpora fail tanh calibration, so the optimizer-correctness check does not pass and the gate-comparison results should not yet be interpreted.")
    else:
        lines.extend(
            [
                "## Independent-gate winner counts",
                "",
                f"- logistic: `{independent_winner_counts['logistic']}`",
                f"- erf: `{independent_winner_counts['erf']}`",
                f"- algebraic: `{independent_winner_counts['algebraic']}`",
                f"- arctan: `{independent_winner_counts['arctan']}`",
                "",
                f"- median BIC spread (excluding tanh): `{statistics.median(bic_spreads):.12f}`",
                f"- mean BIC spread (excluding tanh): `{statistics.fmean(bic_spreads):.12f}`",
                f"- strict indistinguishable count (<2): `{strict_count}`",
                f"- positive indistinguishable count (<6): `{positive_count}`",
                f"- strong indistinguishable count (<10): `{strong_count}`",
            ]
        )

    write_csv(
        OUTDIR / "s2_per_corpus_results.csv",
        result_rows,
        [
            "corpus",
            "logistic_bic",
            "logistic_rmse",
            "logistic_k",
            "logistic_w",
            "tanh_bic",
            "tanh_rmse",
            "tanh_k",
            "tanh_w",
            "erf_bic",
            "erf_rmse",
            "erf_k",
            "erf_w",
            "algebraic_bic",
            "algebraic_rmse",
            "algebraic_k",
            "algebraic_w",
            "arctan_bic",
            "arctan_rmse",
            "arctan_k",
            "arctan_w",
            "bic_tanh_minus_logistic",
            "w_tanh_over_2w_logistic_ratio",
            "tanh_calibration_pass",
            "winner_gate",
            "worst_gate",
            "bic_spread",
        ],
    )
    write_csv(
        OUTDIR / "s2_gate_params_per_corpus.csv",
        param_rows,
        [
            "corpus",
            "logistic_k",
            "logistic_w",
            "logistic_hit_upper_bound",
            "logistic_hit_lower_bound",
            "tanh_k",
            "tanh_w",
            "tanh_hit_upper_bound",
            "tanh_hit_lower_bound",
            "erf_k",
            "erf_w",
            "erf_hit_upper_bound",
            "erf_hit_lower_bound",
            "algebraic_k",
            "algebraic_w",
            "algebraic_hit_upper_bound",
            "algebraic_hit_lower_bound",
            "arctan_k",
            "arctan_w",
            "arctan_hit_upper_bound",
            "arctan_hit_lower_bound",
        ],
    )
    write_csv(
        OUTDIR / "s2_aggregate_statistics.csv",
        aggregate_rows,
        ["metric_name", "value", "notes"],
    )
    (OUTDIR / "s2_summary_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
