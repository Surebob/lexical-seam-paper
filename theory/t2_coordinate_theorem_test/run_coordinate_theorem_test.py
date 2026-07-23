from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import poisson


ROOT = Path("/Volumes/External2TB/emlexperiment")
T2_SCRIPT = ROOT / "phase2_addon" / "t2_redesigned" / "run_t2_redesigned.py"
OUTDIR = ROOT / "phase2_addon" / "t2_coordinate_theorem_test"


def load_t2_module():
    spec = importlib.util.spec_from_file_location("t2_redesigned_for_coordinate_theorem", T2_SCRIPT)
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
            "fsecond_positive_on_domain": False,
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
        "fsecond_positive_on_domain": bool(fsecond_min > -1e-5),
        "bregman_conditions_pass": bool(abs(f_at_1) <= 1e-5 and abs(fprime_at_1) <= 1e-4 and fsecond_min > -1e-5),
    }


def poisson_mixture_residual_99pct():
    lambda_a = 5.0
    lambda_b = 20.0
    weight = 0.5
    fitted_lambda = weight * lambda_a + weight * lambda_b

    # Choose the smallest k_max whose prefix support captures at least 99% of mixture mass.
    k = np.arange(0, 200, dtype=np.int64)
    mixture = weight * poisson.pmf(k, lambda_a) + weight * poisson.pmf(k, lambda_b)
    cumulative = np.cumsum(mixture)
    k_max = int(np.searchsorted(cumulative, 0.99, side="left"))
    support = np.arange(0, k_max + 1, dtype=np.int64)
    q = weight * poisson.pmf(support, lambda_a) + weight * poisson.pmf(support, lambda_b)
    p = poisson.pmf(support, fitted_lambda)
    residual = np.log(q) - np.log(p)
    x_linear = 0.05 + 0.95 * support.astype(np.float64) / float(k_max)
    rows = []
    for idx, count in enumerate(support):
        rows.append(
            {
                "k": int(count),
                "x_linear": f"{x_linear[idx]:.17g}",
                "empirical_mixture_prob": f"{q[idx]:.17g}",
                "fitted_poisson_prob": f"{p[idx]:.17g}",
                "residual": f"{residual[idx]:.17g}",
            }
        )
    return {
        "lambda_a": lambda_a,
        "lambda_b": lambda_b,
        "weight_a": weight,
        "weight_b": weight,
        "fitted_lambda": fitted_lambda,
        "k_max": k_max,
        "captured_mass": float(np.sum(q)),
        "support": support,
        "x_linear": x_linear,
        "residual": residual,
        "rows": rows,
    }


def top20_rows(steps: list[dict], x: np.ndarray):
    predicted = t2.generator_values("poisson_kl", x)
    rows = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"][:20], start=1):
            values = np.asarray(candidate["values"], dtype=np.float64)
            diag = bregman_boundary_diagnostics(candidate)
            rows.append(
                {
                    "step": int(step_payload["step"]),
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    "is_shifted_kl_generator": bool(t2.is_predicted("poisson_kl", values, x)),
                    "cosine_vs_shifted_kl": f"{cosine(values, predicted):.17g}",
                    "cosine_vs_negative_shifted_kl": f"{cosine(values, -predicted):.17g}",
                    "f_at_1": f"{diag['f_at_1']:.17g}",
                    "fprime_at_1": f"{diag['fprime_at_1']:.17g}",
                    "fsecond_min": f"{diag['fsecond_min']:.17g}",
                    "fsecond_positive_on_domain": diag["fsecond_positive_on_domain"],
                    "bregman_conditions_pass": diag["bregman_conditions_pass"],
                }
            )
    return rows


def target_occurrences(steps: list[dict], target: np.ndarray):
    hits = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"], start=1):
            if np.max(np.abs(np.asarray(candidate["values"], dtype=np.float64) - target)) <= 1e-8:
                hits.append(
                    {
                        "step": int(step_payload["step"]),
                        "rmse_rank": rank,
                        "expr": candidate["expr"],
                        "math": candidate["math"],
                        "rmse": float(candidate["rmse"]),
                        "in_top5_by_rmse": bool(rank <= 5),
                        "in_top20_by_rmse": bool(rank <= 20),
                    }
                )
                break
    return hits


def best_cosine_candidate(steps: list[dict], target: np.ndarray):
    best = None
    for rank, candidate in enumerate(steps[-1]["generated"], start=1):
        values = np.asarray(candidate["values"], dtype=np.float64)
        cos = cosine(values, target)
        if best is None or abs(cos) > abs(best["cosine"]):
            best = {
                "step": int(steps[-1]["step"]),
                "rmse_rank": rank,
                "expr": candidate["expr"],
                "math": candidate["math"],
                "rmse": float(candidate["rmse"]),
                "cosine": cos,
            }
    return best


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    data = poisson_mixture_residual_99pct()
    write_csv(
        OUTDIR / "poisson_count_coordinate_residual.csv",
        data["rows"],
        ["k", "x_linear", "empirical_mixture_prob", "fitted_poisson_prob", "residual"],
    )

    x = data["x_linear"]
    residual = data["residual"]
    steps, global_best, occurrences = t2.run_search_detailed("poisson_kl", x, residual, 3)
    beam_rows = top20_rows(steps, x)
    write_csv(
        OUTDIR / "poisson_count_coordinate_top20_beam.csv",
        beam_rows,
        [
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "is_shifted_kl_generator",
            "cosine_vs_shifted_kl",
            "cosine_vs_negative_shifted_kl",
            "f_at_1",
            "fprime_at_1",
            "fsecond_min",
            "fsecond_positive_on_domain",
            "bregman_conditions_pass",
        ],
    )

    shifted_kl = t2.generator_values("poisson_kl", x)
    unshifted_negative_entropy = x * np.log(x) - x
    winner_values = np.asarray(global_best["values"], dtype=np.float64)
    winner_diag = bregman_boundary_diagnostics(global_best)
    unshifted_occurrences = target_occurrences(steps, unshifted_negative_entropy)
    predicted_in_top5 = any(item["step"] == 3 and item["rmse_rank"] <= 5 for item in occurrences + unshifted_occurrences)
    predicted_in_top20 = any(item["step"] == 3 and item["rmse_rank"] <= 20 for item in occurrences + unshifted_occurrences)
    top5_step3 = [
        {
            "rank": idx,
            "expression": item["expr"],
            "math": item["math"],
            "rmse": float(item["rmse"]),
            "is_shifted_kl_generator": bool(t2.is_predicted("poisson_kl", item["values"], x)),
            "cosine_vs_shifted_kl": cosine(np.asarray(item["values"], dtype=np.float64), shifted_kl),
            "cosine_vs_unshifted_negative_entropy": cosine(np.asarray(item["values"], dtype=np.float64), unshifted_negative_entropy),
        }
        for idx, item in enumerate(steps[-1]["generated"][:5], start=1)
    ]
    verdict = {
        "test": "poisson_count_coordinate_theorem_prediction",
        "rates": [data["lambda_a"], data["lambda_b"]],
        "weights": [data["weight_a"], data["weight_b"]],
        "fitted_lambda": data["fitted_lambda"],
        "k_max_99pct_mass": data["k_max"],
        "captured_mass": data["captured_mass"],
        "support_size": int(len(x)),
        "residual_rmse": float(np.sqrt(np.mean(residual**2))),
        "predicted_shifted_kl_generator": "x*log(x)-x+1",
        "predicted_unshifted_negative_entropy": "x*log(x)-x",
        "shifted_kl_occurrences": occurrences,
        "unshifted_negative_entropy_occurrences": unshifted_occurrences,
        "predicted_in_top5_at_depth3": bool(predicted_in_top5),
        "predicted_in_top20_at_depth3": bool(predicted_in_top20),
        "winner": {
            "expression": global_best["expr"],
            "math": global_best["math"],
            "rmse": float(global_best["rmse"]),
            "cosine_vs_shifted_kl": cosine(winner_values, shifted_kl),
            "cosine_vs_negative_shifted_kl": cosine(winner_values, -shifted_kl),
            "max_abs_cosine_vs_shifted_kl": max(abs(cosine(winner_values, shifted_kl)), abs(cosine(winner_values, -shifted_kl))),
            "cosine_vs_unshifted_negative_entropy": cosine(winner_values, unshifted_negative_entropy),
            **winner_diag,
        },
        "top5_depth3": top5_step3,
        "best_cosine_candidate_vs_shifted_kl": best_cosine_candidate(steps, shifted_kl),
        "best_cosine_candidate_vs_unshifted_negative_entropy": best_cosine_candidate(steps, unshifted_negative_entropy),
        "top20_bregman_count": sum(1 for row in beam_rows if row["bregman_conditions_pass"]),
        "primary_verdict": "CONFIRMED" if predicted_in_top5 else "FAILED_PRIMARY_TOP5_CRITERION",
    }
    write_json(OUTDIR / "coordinate_theorem_verdict.json", verdict)

    occurrence_text = "not found"
    if occurrences:
        first = occurrences[0]
        occurrence_text = f"step {first['step']}, RMSE rank {first['rmse_rank']}, RMSE {first['rmse']:.6g}"
    unshifted_occurrence_text = "not found"
    if unshifted_occurrences:
        first = unshifted_occurrences[0]
        unshifted_occurrence_text = f"step {first['step']}, RMSE rank {first['rmse_rank']}, RMSE {first['rmse']:.6g}"
    lines = [
        "# T2 Coordinate-Theorem Prediction Test",
        "",
        "This primary test expresses the Poisson two-component mixture residual in affine-normalized count coordinate, `x = 0.05 + 0.95*k/k_max`, rather than log-rank or sorted-rank coordinate.",
        "",
        f"The `(5, 20)` equal-weight mixture is fit by a single Poisson with lambda `{data['fitted_lambda']}`. The support `k=0..{data['k_max']}` captures `{data['captured_mass']:.6f}` of mixture probability mass.",
        "",
        f"- residual RMSE: `{verdict['residual_rmse']:.6g}`",
        f"- shifted KL occurrence `x log(x)-x+1`: `{occurrence_text}`",
        f"- unshifted negative-entropy occurrence `x log(x)-x`: `{unshifted_occurrence_text}`",
        f"- either predicted form in top 5 at depth 3: `{predicted_in_top5}`",
        f"- either predicted form in top 20 at depth 3: `{predicted_in_top20}`",
        f"- winner: `{global_best['expr']}` ({global_best['math']})",
        f"- winner RMSE: `{float(global_best['rmse']):.6g}`",
        f"- winner max |cosine vs ± shifted KL|: `{verdict['winner']['max_abs_cosine_vs_shifted_kl']:.6g}`",
        f"- best cosine candidate vs shifted KL: rank `{verdict['best_cosine_candidate_vs_shifted_kl']['rmse_rank']}`, cosine `{verdict['best_cosine_candidate_vs_shifted_kl']['cosine']:.6g}`, RMSE `{verdict['best_cosine_candidate_vs_shifted_kl']['rmse']:.6g}`",
        f"- winner satisfies Bregman boundary conditions: `{winner_diag['bregman_conditions_pass']}`",
        "",
        "## Top 5 at depth 3",
        "",
    ]
    for row in top5_step3:
        lines.append(
            f"- rank `{row['rank']}`: `{row['expression']}` RMSE `{row['rmse']:.6g}`, cosine vs shifted KL `{row['cosine_vs_shifted_kl']:.6g}`"
        )
    lines.extend(
        [
            "",
            "## Verdict",
            "",
            (
            "The primary theorem prediction is confirmed under this affine count-coordinate test: an allowed KL/negative-entropy form appears in the top 5 at depth 3."
            if predicted_in_top5
            else "The primary theorem prediction fails under this affine count-coordinate test: neither the shifted KL generator nor the unshifted negative-entropy form appears in the top 5 at depth 3."
        ),
        ]
    )
    (OUTDIR / "coordinate_theorem_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"outdir": str(OUTDIR), "primary_verdict": verdict["primary_verdict"]}, indent=2))


if __name__ == "__main__":
    main()
