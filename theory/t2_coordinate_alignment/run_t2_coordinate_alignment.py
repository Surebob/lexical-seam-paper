from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import poisson


ROOT = Path("/Volumes/External2TB/emlexperiment")
T2_SCRIPT = ROOT / "phase2_addon" / "t2_redesigned" / "run_t2_redesigned.py"
OUTDIR = ROOT / "phase2_addon" / "t2_coordinate_alignment"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from zipf_analysis_common import (  # noqa: E402
    build_zipf_dataset,
    corpus_path,
    get_corpus_spec,
    load_enriched_summary,
    normalize_x,
    zm_prediction,
)


def load_t2_module():
    spec = importlib.util.spec_from_file_location("t2_redesigned_for_coordinate_alignment", T2_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


t2 = load_t2_module()
search = t2.search


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0 or not math.isfinite(denom):
        return float("nan")
    return float(np.dot(a, b) / denom)


def residual_from_probs(q: np.ndarray, p: np.ndarray) -> np.ndarray:
    eps = 1e-300
    return np.log(np.maximum(q, eps)) - np.log(np.maximum(p, eps))


def poisson_mixture_rates_5_20():
    lambda_a = 5.0
    lambda_b = 20.0
    fitted_lambda = 0.5 * lambda_a + 0.5 * lambda_b
    variance = 0.5 * lambda_a + 0.5 * lambda_b + 0.25 * (lambda_b - lambda_a) ** 2
    support_hi = int(math.ceil(fitted_lambda + 12.0 * math.sqrt(variance)))
    support = np.arange(0, support_hi + 1, dtype=np.int64)
    q = 0.5 * poisson.pmf(support, lambda_a) + 0.5 * poisson.pmf(support, lambda_b)
    p = poisson.pmf(support, fitted_lambda)
    q = q / float(np.sum(q))
    p = p / float(np.sum(p))
    return support, q, p, fitted_lambda


def eval_node(node, x: np.ndarray):
    if isinstance(node, str):
        if node == "x":
            return x
        if node == "1":
            return np.ones_like(x)
        raise ValueError(f"Unsupported terminal {node!r}")
    op = node[0]
    if op in search.UNARY_FUNCS:
        child = eval_node(node[1], x)
        if child is None:
            return None
        return search.UNARY_FUNCS[op](child)
    if op in search.BINARY_FUNCS:
        left = eval_node(node[1], x)
        right = eval_node(node[2], x)
        if left is None or right is None:
            return None
        return search.BINARY_FUNCS[op](left, right)
    raise ValueError(f"Unsupported op {op!r}")


def bregman_boundary_diagnostics(candidate: dict) -> dict:
    x_dense = np.linspace(0.05, 1.0, 5000, dtype=np.float64)
    values = eval_node(candidate["node"], x_dense)
    if values is None or not np.all(np.isfinite(values)):
        return {
            "f_at_1": float("nan"),
            "fprime_at_1": float("nan"),
            "fsecond_min": float("nan"),
            "bregman_conditions_pass": False,
        }
    first = np.gradient(values, x_dense, edge_order=2)
    second = np.gradient(first, x_dense, edge_order=2)
    f_at_1 = float(values[-1])
    fprime_at_1 = float(first[-1])
    fsecond_min = float(np.min(second))
    return {
        "f_at_1": f_at_1,
        "fprime_at_1": fprime_at_1,
        "fsecond_min": fsecond_min,
        "bregman_conditions_pass": bool(abs(f_at_1) <= 1e-5 and abs(fprime_at_1) <= 1e-4 and fsecond_min > -1e-5),
    }


def beam_top20_rows(test_name: str, kind: str, steps: list[dict], x: np.ndarray):
    predicted = t2.generator_values(kind, x)
    rows = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"][:20], start=1):
            values = np.asarray(candidate["values"], dtype=np.float64)
            diag = bregman_boundary_diagnostics(candidate)
            rows.append(
                {
                    "test": test_name,
                    "step": int(step_payload["step"]),
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    "is_predicted_generator": bool(t2.is_predicted(kind, values, x)),
                    "cosine_vs_predicted": f"{cosine(values, predicted):.17g}",
                    "cosine_vs_negative_predicted": f"{cosine(values, -predicted):.17g}",
                    "f_at_1": f"{diag['f_at_1']:.17g}",
                    "fprime_at_1": f"{diag['fprime_at_1']:.17g}",
                    "fsecond_min": f"{diag['fsecond_min']:.17g}",
                    "bregman_conditions_pass": diag["bregman_conditions_pass"],
                }
            )
    return rows


def search_summary(test_name: str, kind: str, x: np.ndarray, residual: np.ndarray, max_steps: int):
    steps, global_best, occurrences = t2.run_search_detailed(kind, x, residual, max_steps)
    predicted = t2.generator_values(kind, x)
    winner = np.asarray(global_best["values"], dtype=np.float64)
    top20 = beam_top20_rows(test_name, kind, steps, x)
    bregman_count = sum(1 for row in top20 if row["bregman_conditions_pass"])
    return {
        "steps": steps,
        "top20_rows": top20,
        "summary": {
            "test": test_name,
            "kind": kind,
            "support_size": int(len(x)),
            "residual_rmse": float(np.sqrt(np.mean(np.asarray(residual) ** 2))),
            "winner_expression": global_best["expr"],
            "winner_math": global_best["math"],
            "winner_rmse": float(global_best["rmse"]),
            "winner_cosine_vs_predicted": cosine(winner, predicted),
            "winner_cosine_vs_negative_predicted": cosine(winner, -predicted),
            "winner_max_abs_cosine_vs_predicted": max(abs(cosine(winner, predicted)), abs(cosine(winner, -predicted))),
            "predicted_occurrences": occurrences,
            "top20_bregman_count": int(bregman_count),
        },
    }


def run_poisson_sorted_rank():
    support, q, p, fitted_lambda = poisson_mixture_rates_5_20()
    residual = residual_from_probs(q, p)
    order = sorted(range(len(support)), key=lambda idx: (-float(p[idx]), int(support[idx])))
    order = np.asarray(order, dtype=np.int64)
    support_sorted = support[order]
    q_sorted = q[order]
    p_sorted = p[order]
    residual_sorted = residual[order]
    x_sorted_rank = t2.normalized_sr_coordinate(len(support_sorted))
    count_over_lambda = support_sorted.astype(np.float64) / fitted_lambda
    with np.errstate(divide="ignore", invalid="ignore"):
        deviance = np.where(
            support_sorted > 0,
            support_sorted.astype(np.float64) * np.log(support_sorted.astype(np.float64) / fitted_lambda) - support_sorted + fitted_lambda,
            fitted_lambda,
        )
    rows = []
    for idx, count in enumerate(support_sorted):
        rows.append(
            {
                "sorted_position": idx + 1,
                "count": int(count),
                "empirical_prob_mixture": f"{q_sorted[idx]:.17g}",
                "fitted_prob_single_poisson": f"{p_sorted[idx]:.17g}",
                "residual_log_q_minus_log_p": f"{residual_sorted[idx]:.17g}",
                "sr_x_sorted_rank": f"{x_sorted_rank[idx]:.17g}",
                "count_over_lambda": f"{count_over_lambda[idx]:.17g}",
                "poisson_deviance_count_lambda": f"{deviance[idx]:.17g}",
            }
        )
    write_csv(
        OUTDIR / "sorted_rank_poisson_residual.csv",
        rows,
        [
            "sorted_position",
            "count",
            "empirical_prob_mixture",
            "fitted_prob_single_poisson",
            "residual_log_q_minus_log_p",
            "sr_x_sorted_rank",
            "count_over_lambda",
            "poisson_deviance_count_lambda",
        ],
    )
    result = search_summary("poisson_sorted_rank", "poisson_kl", x_sorted_rank, residual_sorted, 3)
    write_csv(
        OUTDIR / "poisson_sorted_rank_beam_top20.csv",
        result["top20_rows"],
        [
            "test",
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "is_predicted_generator",
            "cosine_vs_predicted",
            "cosine_vs_negative_predicted",
            "f_at_1",
            "fprime_at_1",
            "fsecond_min",
            "bregman_conditions_pass",
        ],
    )
    result["summary"].update(
        {
            "fitted_lambda": fitted_lambda,
            "sort_rule": "descending fitted single-Poisson probability, then ascending count",
            "residual_vs_deviance_cosine_in_sorted_order": cosine(residual_sorted, deviance),
            "residual_vs_kl_sr_coordinate_cosine": cosine(residual_sorted, t2.generator_values("poisson_kl", x_sorted_rank)),
        }
    )
    return result["summary"]


def run_zm_natural_parameter():
    spec = get_corpus_spec("shakespeare")
    dataset = build_zipf_dataset(corpus_path(spec))
    summary = load_enriched_summary(spec)
    ranks = dataset["ranks"]
    log_freq = dataset["log_freq"]
    pred = zm_prediction(summary, ranks)
    residual = log_freq - pred
    c = float(summary["zm_baseline"]["c"])
    natural_raw = c / (ranks + c)
    natural_x = normalize_x(natural_raw, 0.05, 1.0)
    log_rank_x = normalize_x(np.log(ranks), 0.05, 1.0)
    log_rplusc_x = normalize_x(np.log(ranks + c), 0.05, 1.0)
    rows = []
    for idx, (word, freq) in enumerate(dataset["ranked"]):
        rows.append(
            {
                "rank": int(ranks[idx]),
                "word": word,
                "frequency": int(freq),
                "log_frequency": f"{log_freq[idx]:.17g}",
                "zm_log_prediction": f"{pred[idx]:.17g}",
                "residual_log_freq_minus_zm": f"{residual[idx]:.17g}",
                "sr_x_log_rank": f"{log_rank_x[idx]:.17g}",
                "natural_raw_c_over_r_plus_c": f"{natural_raw[idx]:.17g}",
                "natural_param_x_c_over_r_plus_c": f"{natural_x[idx]:.17g}",
                "reference_x_log_r_plus_c": f"{log_rplusc_x[idx]:.17g}",
            }
        )
    write_csv(
        OUTDIR / "natural_param_zm_residual.csv",
        rows,
        [
            "rank",
            "word",
            "frequency",
            "log_frequency",
            "zm_log_prediction",
            "residual_log_freq_minus_zm",
            "sr_x_log_rank",
            "natural_raw_c_over_r_plus_c",
            "natural_param_x_c_over_r_plus_c",
            "reference_x_log_r_plus_c",
        ],
    )
    # Expressions are evaluated on each coordinate's values paired with the original residuals.
    result = search_summary("zm_natural_param_c_over_r_plus_c", "gamma_is", natural_x, residual, 2)
    log_rank_control = search_summary("zm_log_rank_control", "gamma_is", log_rank_x, residual, 2)
    log_rplusc_reference = search_summary("zm_log_rplusc_reference", "gamma_is", log_rplusc_x, residual, 2)
    write_csv(
        OUTDIR / "zm_natural_param_beam_top20.csv",
        result["top20_rows"],
        [
            "test",
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "is_predicted_generator",
            "cosine_vs_predicted",
            "cosine_vs_negative_predicted",
            "f_at_1",
            "fprime_at_1",
            "fsecond_min",
            "bregman_conditions_pass",
        ],
    )
    for filename, payload in [
        ("zm_log_rank_control_beam_top20.csv", log_rank_control),
        ("zm_log_rplusc_reference_beam_top20.csv", log_rplusc_reference),
    ]:
        write_csv(
            OUTDIR / filename,
            payload["top20_rows"],
            [
                "test",
                "step",
                "rank",
                "expression",
                "math",
                "rmse",
                "is_predicted_generator",
                "cosine_vs_predicted",
                "cosine_vs_negative_predicted",
                "f_at_1",
                "fprime_at_1",
                "fsecond_min",
                "bregman_conditions_pass",
            ],
        )
    result["summary"].update(
        {
            "corpus": spec["name"],
            "vocabulary_size": int(dataset["unique_words"]),
            "zm_c": c,
            "coordinate_primary": "linear-normalized c/(r+c)",
            "reference_log_rank_is_cosine": cosine(residual, t2.generator_values("gamma_is", log_rank_x)),
            "natural_param_is_cosine": cosine(residual, t2.generator_values("gamma_is", natural_x)),
            "reference_log_rplusc_is_cosine": cosine(residual, t2.generator_values("gamma_is", log_rplusc_x)),
            "log_rank_control": log_rank_control["summary"],
            "log_rplusc_reference": log_rplusc_reference["summary"],
        }
    )
    return result["summary"]


def write_report(poisson_summary: dict, zm_summary: dict):
    def occurrence_text(summary: dict) -> str:
        if not summary["predicted_occurrences"]:
            return "not found"
        first = summary["predicted_occurrences"][0]
        return f"step {first['step']}, RMSE rank {first['rmse_rank']}, RMSE {first['rmse']:.6g}"

    if poisson_summary["predicted_occurrences"] and zm_summary["predicted_occurrences"]:
        matrix = "Both predicted generators appeared under the tested coordinates."
    elif poisson_summary["predicted_occurrences"] and not zm_summary["predicted_occurrences"]:
        matrix = "KL appears in sorted-rank Poisson, but IS does not appear in the tested ZM natural-coordinate frame."
    elif (not poisson_summary["predicted_occurrences"]) and zm_summary["predicted_occurrences"]:
        matrix = "KL does not appear in sorted-rank Poisson, while IS appears in the tested ZM natural-coordinate frame."
    else:
        matrix = "Neither predicted generator appears under the tested coordinate transforms."

    lines = [
        "# T2 Coordinate-Alignment Verification",
        "",
        "This addendum tests whether the T2 Poisson failure is primarily a coordinate-alignment issue. It compares a sorted-rank projection of the Poisson residual with a ZM residual expressed in a Mandelbrot-shift natural-coordinate proxy.",
        "",
        "## Test 1: Poisson sorted-rank coordinate",
        "",
        "The `(5, 20)` equal-weight Poisson mixture was fit by a single Poisson with lambda `12.5`. Counts were sorted by descending fitted single-Poisson probability and mapped to the manuscript normalized log-rank coordinate.",
        "",
        f"- winner: `{poisson_summary['winner_expression']}` ({poisson_summary['winner_math']})",
        f"- residual RMSE: `{poisson_summary['residual_rmse']:.6g}`",
        f"- KL occurrence: `{occurrence_text(poisson_summary)}`",
        f"- winner max |cosine vs ±KL|: `{poisson_summary['winner_max_abs_cosine_vs_predicted']:.6g}`",
        f"- residual vs Poisson deviance cosine after sorting: `{poisson_summary['residual_vs_deviance_cosine_in_sorted_order']:.6g}`",
        f"- residual vs KL on sorted-rank SR coordinate cosine: `{poisson_summary['residual_vs_kl_sr_coordinate_cosine']:.6g}`",
        f"- top-20 Bregman-condition count: `{poisson_summary['top20_bregman_count']}`",
        "",
        "## Test 2: ZM residual in `c/(r+c)` coordinate",
        "",
        "The Shakespeare single-ZM residual was paired with a linear-normalized `c/(r+c)` coordinate. This is a local coordinate induced by the Mandelbrot shift parameter; the output CSV also includes log-rank and log(r+c) reference coordinates.",
        "",
        f"- winner: `{zm_summary['winner_expression']}` ({zm_summary['winner_math']})",
        f"- residual RMSE: `{zm_summary['residual_rmse']:.6g}`",
        f"- IS occurrence: `{occurrence_text(zm_summary)}`",
        f"- winner max |cosine vs ±IS|: `{zm_summary['winner_max_abs_cosine_vs_predicted']:.6g}`",
        f"- residual vs IS on log-rank x cosine: `{zm_summary['reference_log_rank_is_cosine']:.6g}`",
        f"- residual vs IS on `c/(r+c)` x cosine: `{zm_summary['natural_param_is_cosine']:.6g}`",
        f"- residual vs IS on log(r+c) x cosine: `{zm_summary['reference_log_rplusc_is_cosine']:.6g}`",
        f"- top-20 Bregman-condition count: `{zm_summary['top20_bregman_count']}`",
        "",
        "### ZM coordinate controls",
        "",
        f"- log-rank control winner: `{zm_summary['log_rank_control']['winner_expression']}`; IS occurrence: `{occurrence_text(zm_summary['log_rank_control'])}`",
        f"- log(r+c) reference winner: `{zm_summary['log_rplusc_reference']['winner_expression']}`; IS occurrence: `{occurrence_text(zm_summary['log_rplusc_reference'])}`",
        "",
        "## Interpretation Matrix",
        "",
        matrix,
        "",
        "Under this implementation, the Bregman generator recovery is coordinate-dependent rather than invariant. The Poisson deviance remains present in natural count/lambda space, but the sorted-rank projection does not make KL a competitive SR expression. Conversely, the Shakespeare ZM residual does not recover IS under the tested `c/(r+c)` coordinate, even though the manuscript log-rank coordinate is the one where IS is known to be recovered.",
    ]
    (OUTDIR / "coordinate_alignment_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    poisson_summary = run_poisson_sorted_rank()
    zm_summary = run_zm_natural_parameter()
    payload = {
        "poisson_sorted_rank": poisson_summary,
        "zm_natural_parameter": zm_summary,
    }
    write_json(OUTDIR / "coordinate_alignment_summary.json", payload)
    write_report(poisson_summary, zm_summary)
    print(json.dumps({"outdir": str(OUTDIR), "poisson_kl_found": bool(poisson_summary["predicted_occurrences"]), "zm_is_found": bool(zm_summary["predicted_occurrences"])}, indent=2))


if __name__ == "__main__":
    main()
