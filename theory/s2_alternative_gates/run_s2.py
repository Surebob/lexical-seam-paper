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
OUTDIR = ROOT / "phase2_addon" / "s2_alternative_gates"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import zipf_analysis_common as common


RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
SMOOTH_FIT_CSV = ROOT / "experiments" / "3a_smooth_two_regime_fits" / "outputs" / "smooth_fit_per_corpus.csv"

N_RANDOM_STARTS = 100
MAX_NFEV = 12000
PARAM_COUNT = 8


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reranked = load_module(RERANKED_PATH, "s2_alt_gate_reranked")


def bic_from_rmse(rmse: float, p: int, n: int) -> float:
    mse = float(rmse) ** 2
    return p * math.log(n) + n * math.log(mse)


def tanh_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = (np.log(ranks) - math.log(k)) / w
    return 0.5 * (1.0 - np.tanh(z))


def erf_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = (np.log(ranks) - math.log(k)) / w
    return 0.5 * (1.0 - special.erf(z))


def algebraic_gate(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    z = np.log(ranks) - math.log(k)
    return 0.5 * (1.0 - z / np.sqrt((w * w) + (z * z)))


GATE_FUNCS = {
    "tanh": tanh_gate,
    "erf": erf_gate,
    "algebraic": algebraic_gate,
}


def load_logistic_rows() -> list[dict[str, str]]:
    with SMOOTH_FIT_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fit_with_gate(dataset: dict, gate_name: str) -> dict[str, float | int]:
    original_sigma = reranked.sigma_curve
    reranked.sigma_curve = GATE_FUNCS[gate_name]
    try:
        best, tried = reranked.fit_reranked_model(dataset, n_starts=N_RANDOM_STARTS, max_nfev=MAX_NFEV)
    finally:
        reranked.sigma_curve = original_sigma
    params = reranked.summarize_params(best["params"], dataset["unique_words"])
    return {
        "rmse": float(best["rmse"]),
        "bic": float(bic_from_rmse(best["rmse"], PARAM_COUNT, int(dataset["unique_words"]))),
        "k": float(params["k"]),
        "w": float(params["w"]),
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "tried_count": len(tried),
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    logistic_rows = load_logistic_rows()

    result_rows: list[dict[str, str]] = []
    param_rows: list[dict[str, str]] = []
    bic_spreads: list[float] = []
    winner_counts = {gate: 0 for gate in ["logistic", "tanh", "erf", "algebraic"]}
    indistinguishable_count = 0

    for logistic_row in logistic_rows:
        slug = logistic_row["slug"]
        spec = common.get_corpus_spec(slug)
        dataset = common.build_zipf_dataset(common.corpus_path(spec))
        vocab_size = int(logistic_row["vocab_size"])

        logistic_rmse = float(logistic_row["reranked_rmse"])
        logistic_bic = bic_from_rmse(logistic_rmse, PARAM_COUNT, vocab_size)
        logistic_k = float(logistic_row["k"])
        logistic_w = float(logistic_row["w"])

        alt = {name: fit_with_gate(dataset, name) for name in GATE_FUNCS}
        bics = {
            "logistic": logistic_bic,
            "tanh": float(alt["tanh"]["bic"]),
            "erf": float(alt["erf"]["bic"]),
            "algebraic": float(alt["algebraic"]["bic"]),
        }
        winner_gate = min(bics, key=bics.get)
        worst_gate = max(bics, key=bics.get)
        bic_spread = max(bics.values()) - min(bics.values())
        bic_spreads.append(bic_spread)
        winner_counts[winner_gate] += 1
        if bic_spread < 10.0:
            indistinguishable_count += 1

        result_rows.append(
            {
                "corpus": logistic_row["corpus"],
                "logistic_bic": repr(logistic_bic),
                "logistic_rmse": repr(logistic_rmse),
                "tanh_bic": repr(float(alt["tanh"]["bic"])),
                "tanh_rmse": repr(float(alt["tanh"]["rmse"])),
                "erf_bic": repr(float(alt["erf"]["bic"])),
                "erf_rmse": repr(float(alt["erf"]["rmse"])),
                "algebraic_bic": repr(float(alt["algebraic"]["bic"])),
                "algebraic_rmse": repr(float(alt["algebraic"]["rmse"])),
                "winner_gate": winner_gate,
                "worst_gate": worst_gate,
                "bic_spread": repr(float(bic_spread)),
            }
        )
        param_rows.append(
            {
                "corpus": logistic_row["corpus"],
                "logistic_k": repr(logistic_k),
                "logistic_w": repr(logistic_w),
                "tanh_k": repr(float(alt["tanh"]["k"])),
                "tanh_w": repr(float(alt["tanh"]["w"])),
                "erf_k": repr(float(alt["erf"]["k"])),
                "erf_w": repr(float(alt["erf"]["w"])),
                "algebraic_k": repr(float(alt["algebraic"]["k"])),
                "algebraic_w": repr(float(alt["algebraic"]["w"])),
            }
        )

    aggregate_rows = [
        {
            "metric_name": "logistic_bic_wins",
            "value": str(winner_counts["logistic"]),
            "notes": "Count of corpora where logistic gate has the lowest BIC.",
        },
        {
            "metric_name": "tanh_bic_wins",
            "value": str(winner_counts["tanh"]),
            "notes": "Count of corpora where tanh gate has the lowest BIC.",
        },
        {
            "metric_name": "erf_bic_wins",
            "value": str(winner_counts["erf"]),
            "notes": "Count of corpora where erf gate has the lowest BIC.",
        },
        {
            "metric_name": "algebraic_bic_wins",
            "value": str(winner_counts["algebraic"]),
            "notes": "Count of corpora where algebraic gate has the lowest BIC.",
        },
        {
            "metric_name": "median_bic_spread",
            "value": repr(float(statistics.median(bic_spreads))),
            "notes": "Median of max_bic - min_bic across the four gates per corpus.",
        },
        {
            "metric_name": "gates_indistinguishable_count",
            "value": str(indistinguishable_count),
            "notes": "Count of corpora where max_bic - min_bic < 10.",
        },
        {
            "metric_name": "mean_bic_spread",
            "value": repr(float(statistics.fmean(bic_spreads))),
            "notes": "Mean of max_bic - min_bic across the four gates per corpus.",
        },
        {
            "metric_name": "max_bic_spread",
            "value": repr(float(max(bic_spreads))),
            "notes": "Maximum of max_bic - min_bic across the four gates per corpus.",
        },
        {
            "metric_name": "min_bic_spread",
            "value": repr(float(min(bic_spreads))),
            "notes": "Minimum of max_bic - min_bic across the four gates per corpus.",
        },
    ]

    lines = [
        "# S2 Alternative Smooth Gate Families",
        "",
        "- Alternative gates rerun with the historical reranked smooth-fit engine and optimizer settings.",
        "- Logistic baseline BIC is computed from canonical 3a RMSE with the historical BIC formula `p * ln(n) + n * ln(MSE)` using `p = 8` and `n = vocab_size`.",
        "",
        "## Winner counts",
        "",
        f"- logistic: `{winner_counts['logistic']}`",
        f"- tanh: `{winner_counts['tanh']}`",
        f"- erf: `{winner_counts['erf']}`",
        f"- algebraic: `{winner_counts['algebraic']}`",
        "",
        f"- median BIC spread: `{statistics.median(bic_spreads):.12f}`",
        f"- gates indistinguishable count (`spread < 10`): `{indistinguishable_count}`",
        "",
    ]
    if winner_counts["logistic"] >= 20:
        lines.append("- Interpretation: logistic wins on 20+ corpora, so the specific functional form is empirically supported.")
    elif indistinguishable_count > 12:
        lines.append("- Interpretation: most corpora fall below the BIC-noise threshold, so the gate family is underdetermined by the current data.")
    elif max(winner_counts["tanh"], winner_counts["erf"], winner_counts["algebraic"]) > winner_counts["logistic"]:
        lines.append("- Interpretation: a non-logistic gate wins on a plurality/majority of corpora; this is unexpected and weakens a logistic-specific physics analogy.")
    else:
        lines.append("- Interpretation: no gate family dominates strongly enough for a clean functional-form claim.")

    write_csv(
        OUTDIR / "s2_per_corpus_results.csv",
        result_rows,
        [
            "corpus",
            "logistic_bic",
            "logistic_rmse",
            "tanh_bic",
            "tanh_rmse",
            "erf_bic",
            "erf_rmse",
            "algebraic_bic",
            "algebraic_rmse",
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
            "tanh_k",
            "tanh_w",
            "erf_k",
            "erf_w",
            "algebraic_k",
            "algebraic_w",
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
