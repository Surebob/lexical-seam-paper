from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import nbinom, poisson


ROOT = Path("/Volumes/External2TB/emlexperiment")
T2_SCRIPT = ROOT / "phase2_addon" / "t2_redesigned" / "run_t2_redesigned.py"
OUTDIR = ROOT / "phase2_addon" / "t2_poisson_exhaustive"


def load_t2_module():
    spec = importlib.util.spec_from_file_location("t2_redesigned_for_poisson_exhaustive", T2_SCRIPT)
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


def support_for_poisson_like(mean: float, variance: float | None = None, tail_eps: float = 1e-12) -> np.ndarray:
    variance = float(mean if variance is None else variance)
    hi = int(max(30, math.ceil(mean + 12.0 * math.sqrt(max(variance, 1.0)))))
    return np.arange(0, hi + 1, dtype=np.int64)


def normalize_probs(probs: np.ndarray) -> np.ndarray:
    probs = np.asarray(probs, dtype=np.float64)
    total = float(np.sum(probs))
    if total <= 0.0 or not math.isfinite(total):
        raise RuntimeError("Probability vector has invalid total mass.")
    return probs / total


def poisson_mixture_distribution(lambda_a: float, lambda_b: float):
    mean = 0.5 * lambda_a + 0.5 * lambda_b
    variance = 0.5 * lambda_a + 0.5 * lambda_b + 0.25 * (lambda_a - lambda_b) ** 2
    support = support_for_poisson_like(mean, variance)
    q = 0.5 * poisson.pmf(support, lambda_a) + 0.5 * poisson.pmf(support, lambda_b)
    p = poisson.pmf(support, mean)
    return support, normalize_probs(q), normalize_probs(p), mean, variance


def negbinom_distribution(mean: float = 10.0, dispersion: float = 2.0):
    variance = mean * (1.0 + mean / dispersion)
    p_success = dispersion / (dispersion + mean)
    support = support_for_poisson_like(mean, variance)
    q = nbinom.pmf(support, dispersion, p_success)
    p = poisson.pmf(support, mean)
    return support, normalize_probs(q), normalize_probs(p), mean, variance, p_success


def residual_from_probs(q: np.ndarray, p: np.ndarray) -> np.ndarray:
    eps = 1e-300
    return np.log(np.maximum(q, eps)) - np.log(np.maximum(p, eps))


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


def summarize_search(label: str, residual: np.ndarray):
    x = t2.normalized_sr_coordinate(len(residual))
    steps, global_best, occurrences = t2.run_search_detailed("poisson_kl", x, residual, 3)
    kl = t2.generator_values("poisson_kl", x)
    winner = np.asarray(global_best["values"], dtype=np.float64)
    cos_pos = cosine(winner, kl)
    cos_neg = cosine(winner, -kl)
    occurrence = occurrences[0] if occurrences else None
    occurrence_diag = None
    if occurrence is not None:
        for candidate in steps[occurrence["step"] - 1]["generated"]:
            if candidate["expr"] == occurrence["expr"]:
                occurrence_diag = bregman_boundary_diagnostics(candidate)
                break
    top5 = [
        {
            "step": int(steps[-1]["step"]),
            "rank": idx,
            "expression": item["expr"],
            "math": item["math"],
            "rmse": float(item["rmse"]),
            "is_kl": bool(t2.is_predicted("poisson_kl", item["values"], x)),
        }
        for idx, item in enumerate(steps[-1]["generated"][:5], start=1)
    ]
    top20_bregman = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"][:20], start=1):
            diag = bregman_boundary_diagnostics(candidate)
            if diag["bregman_conditions_pass"]:
                top20_bregman.append(
                    {
                        "step": int(step_payload["step"]),
                        "rank": rank,
                        "expression": candidate["expr"],
                        "math": candidate["math"],
                        "rmse": float(candidate["rmse"]),
                        **diag,
                    }
                )
    return {
        "label": label,
        "support_size": int(len(residual)),
        "residual_rmse": float(np.sqrt(np.mean(residual**2))),
        "winner_expression": global_best["expr"],
        "winner_math": global_best["math"],
        "winner_rmse": float(global_best["rmse"]),
        "winner_cosine_kl": cos_pos,
        "winner_cosine_neg_kl": cos_neg,
        "winner_max_abs_cosine_kl": max(abs(cos_pos), abs(cos_neg)),
        "top5_step3": top5,
        "kl_occurrence": occurrence,
        "kl_occurrence_bregman_diag": occurrence_diag,
        "kl_in_top20_with_bregman_conditions": bool(
            occurrence is not None
            and occurrence["step"] == 3
            and occurrence["rmse_rank"] <= 20
            and occurrence_diag is not None
            and occurrence_diag["bregman_conditions_pass"]
        ),
        "top20_bregman_count": len(top20_bregman),
        "top20_bregman": top20_bregman,
    }


def test_a_rate_ratio_sweep():
    rows = []
    summaries = []
    for ratio in [1.2, 1.5, 2.0, 3.0, 4.0]:
        lambda_a = 10.0
        lambda_b = 10.0 * ratio
        support, q, p, fitted_lambda, variance = poisson_mixture_distribution(lambda_a, lambda_b)
        residual = residual_from_probs(q, p)
        summary = summarize_search(f"ratio_{ratio:g}", residual)
        occurrence = summary["kl_occurrence"] or {}
        summaries.append(
            {
                "ratio": ratio,
                "lambda_a": lambda_a,
                "lambda_b": lambda_b,
                "support": support.tolist(),
                "empirical_prob": q.tolist(),
                "fitted_prob": p.tolist(),
                **summary,
            }
        )
        rows.append(
            {
                "rate_ratio": f"{ratio:.17g}",
                "lambda_a": f"{lambda_a:.17g}",
                "lambda_b": f"{lambda_b:.17g}",
                "fitted_lambda": f"{fitted_lambda:.17g}",
                "support_size": summary["support_size"],
                "residual_rmse": f"{summary['residual_rmse']:.17g}",
                "winner_expression": summary["winner_expression"],
                "winner_math": summary["winner_math"],
                "winner_rmse": f"{summary['winner_rmse']:.17g}",
                "winner_cosine_kl": f"{summary['winner_cosine_kl']:.17g}",
                "winner_cosine_neg_kl": f"{summary['winner_cosine_neg_kl']:.17g}",
                "winner_max_abs_cosine_kl": f"{summary['winner_max_abs_cosine_kl']:.17g}",
                "kl_step": occurrence.get("step", ""),
                "kl_rank": occurrence.get("rmse_rank", ""),
                "kl_rmse": f"{occurrence.get('rmse', float('nan')):.17g}" if occurrence else "",
                "kl_in_top20_with_bregman_conditions": summary["kl_in_top20_with_bregman_conditions"],
                "top20_bregman_count": summary["top20_bregman_count"],
                "top5_step3_json": json.dumps(summary["top5_step3"], separators=(",", ":")),
            }
        )
    write_csv(
        OUTDIR / "t2_poisson_rate_ratio_sweep.csv",
        rows,
        [
            "rate_ratio",
            "lambda_a",
            "lambda_b",
            "fitted_lambda",
            "support_size",
            "residual_rmse",
            "winner_expression",
            "winner_math",
            "winner_rmse",
            "winner_cosine_kl",
            "winner_cosine_neg_kl",
            "winner_max_abs_cosine_kl",
            "kl_step",
            "kl_rank",
            "kl_rmse",
            "kl_in_top20_with_bregman_conditions",
            "top20_bregman_count",
            "top5_step3_json",
        ],
    )
    return summaries


def test_b_truncation():
    support, q, p, fitted_lambda, variance = poisson_mixture_distribution(5.0, 20.0)
    residual = residual_from_probs(q, p)
    by_closeness = np.argsort(np.abs(support.astype(float) - fitted_lambda))
    cumulative = 0.0
    selected = []
    for idx in by_closeness:
        selected.append(idx)
        cumulative += float(q[idx])
        if cumulative >= 0.80:
            break
    selected = np.asarray(sorted(selected), dtype=np.int64)
    truncated_residual = residual[selected]
    summary = summarize_search("original_5_20_truncated_80pct_mass_closest_to_fit", truncated_residual)
    payload = {
        "test": "B_residual_truncation",
        "original_rates": [5.0, 20.0],
        "fitted_lambda": fitted_lambda,
        "selected_mass": float(np.sum(q[selected])),
        "selected_support_min": int(support[selected].min()),
        "selected_support_max": int(support[selected].max()),
        "selected_support_values": support[selected].astype(int).tolist(),
        "selected_support_size": int(len(selected)),
        "full_residual_rmse": float(np.sqrt(np.mean(residual**2))),
        **summary,
    }
    write_json(OUTDIR / "t2_poisson_truncation_test.json", payload)
    return payload


def test_c_negbinom():
    support, q, p, fitted_lambda, variance, p_success = negbinom_distribution(mean=10.0, dispersion=2.0)
    residual = residual_from_probs(q, p)
    summary = summarize_search("negative_binomial_mean10_dispersion2_fit_poisson", residual)
    payload = {
        "test": "C_negative_binomial_fit_by_poisson",
        "negative_binomial_mean": 10.0,
        "dispersion": 2.0,
        "negative_binomial_variance": variance,
        "negative_binomial_success_probability": p_success,
        "fitted_poisson_lambda": fitted_lambda,
        "support_min": int(support.min()),
        "support_max": int(support.max()),
        **summary,
    }
    write_json(OUTDIR / "t2_poisson_negbinom_fit.json", payload)
    return payload


def test_d_theoretical_expansion():
    support, q, p, fitted_lambda, variance = poisson_mixture_distribution(5.0, 20.0)
    residual = residual_from_probs(q, p)
    sr_x = t2.normalized_sr_coordinate(len(support))
    ratio = support.astype(np.float64) / fitted_lambda
    with np.errstate(divide="ignore", invalid="ignore"):
        poisson_deviance = np.where(support > 0, support * np.log(support / fitted_lambda) - support + fitted_lambda, fitted_lambda)
        normalized_deviance = poisson_deviance / fitted_lambda
    kl_sr = t2.generator_values("poisson_kl", sr_x)
    pointwise_kl_contribution = q * residual
    centered_residual = residual - float(np.mean(residual))
    centered_deviance = poisson_deviance - float(np.mean(poisson_deviance))
    centered_kl_sr = kl_sr - float(np.mean(kl_sr))
    rows = []
    for idx, k in enumerate(support):
        rows.append(
            {
                "k": int(k),
                "empirical_prob": f"{q[idx]:.17g}",
                "fitted_prob": f"{p[idx]:.17g}",
                "log_density_residual": f"{residual[idx]:.17g}",
                "pointwise_kl_contribution_q_log_q_over_p": f"{pointwise_kl_contribution[idx]:.17g}",
                "count_over_lambda": f"{ratio[idx]:.17g}",
                "poisson_deviance_D_k_lambda": f"{poisson_deviance[idx]:.17g}",
                "normalized_poisson_deviance": f"{normalized_deviance[idx]:.17g}",
                "sr_x": f"{sr_x[idx]:.17g}",
                "kl_generator_on_sr_x": f"{kl_sr[idx]:.17g}",
            }
        )
    write_csv(
        OUTDIR / "t2_poisson_theoretical_expansion.csv",
        rows,
        [
            "k",
            "empirical_prob",
            "fitted_prob",
            "log_density_residual",
            "pointwise_kl_contribution_q_log_q_over_p",
            "count_over_lambda",
            "poisson_deviance_D_k_lambda",
            "normalized_poisson_deviance",
            "sr_x",
            "kl_generator_on_sr_x",
        ],
    )
    payload = {
        "test": "D_bregman_divergence_expansion_reference",
        "rates": [5.0, 20.0],
        "fitted_lambda": fitted_lambda,
        "derivation_summary": [
            "For the Poisson exponential family, A(theta)=exp(theta) and A*(mu)=mu*log(mu)-mu, so the mean-parameter Bregman divergence is D(mu||lambda)=mu*log(mu/lambda)-mu+lambda.",
            "For a support count k this gives the Poisson deviance D(k||lambda).",
            "The SR residual used here is the pointwise log-density ratio log q(k)-log p_lambda(k), not D(k||lambda). The full distribution KL is sum_k q(k)*[log q(k)-log p_lambda(k)]. Therefore the KL generator does not directly equal the pointwise residual; this is a methodology/theory-object mismatch to test numerically.",
        ],
        "raw_cosine_residual_vs_poisson_deviance": cosine(residual, poisson_deviance),
        "centered_cosine_residual_vs_poisson_deviance": cosine(centered_residual, centered_deviance),
        "raw_cosine_residual_vs_kl_on_sr_x": cosine(residual, kl_sr),
        "centered_cosine_residual_vs_kl_on_sr_x": cosine(centered_residual, centered_kl_sr),
        "raw_cosine_residual_vs_negative_poisson_deviance": cosine(residual, -poisson_deviance),
        "total_distribution_kl_sum_q_log_q_over_p": float(np.sum(pointwise_kl_contribution)),
        "poisson_deviance_exists_in_mean_parameter_space": True,
        "poisson_deviance_directly_equals_pointwise_log_density_residual": False,
    }
    return payload


def report(test_a: list[dict], test_b: dict, test_c: dict, test_d: dict):
    clean_kl = []
    for item in test_a:
        if item["kl_in_top20_with_bregman_conditions"] or item["winner_max_abs_cosine_kl"] >= 0.95:
            clean_kl.append(f"rate ratio {item['ratio']}")
    if test_b["kl_in_top20_with_bregman_conditions"] or test_b["winner_max_abs_cosine_kl"] >= 0.95:
        clean_kl.append("80% mass truncation")
    if test_c["kl_in_top20_with_bregman_conditions"] or test_c["winner_max_abs_cosine_kl"] >= 0.95:
        clean_kl.append("negative binomial misspecification")

    deviance_match = (
        abs(test_d["raw_cosine_residual_vs_poisson_deviance"]) >= 0.95
        or abs(test_d["centered_cosine_residual_vs_poisson_deviance"]) >= 0.95
    )
    sr_coordinate_miss = abs(test_d["raw_cosine_residual_vs_kl_on_sr_x"]) < 0.50
    if clean_kl:
        verdict = "support_under_restricted_condition"
    elif deviance_match and sr_coordinate_miss:
        verdict = "methodology_scope_theory_matches_deviance_coordinate_not_sr_coordinate"
    elif not test_d["poisson_deviance_directly_equals_pointwise_log_density_residual"]:
        verdict = "ambiguous_methodology_theory_object_mismatch"
    else:
        verdict = "problematic_for_subclaim_1"

    lines = [
        "# T2 Poisson Exhaustive Follow-Up",
        "",
        "This follow-up tests whether the redesigned T2 Poisson failure is caused by severity, tail truncation, misspecification choice, or a mismatch between the Poisson KL generator and the pointwise log-density residual given to symbolic regression.",
        "",
        "## Test A: rate-ratio sweep",
        "",
    ]
    for item in test_a:
        occurrence = item["kl_occurrence"] or {}
        lines.append(
            f"- ratio `{item['ratio']}`: residual RMSE `{item['residual_rmse']:.6g}`, winner `{item['winner_expression']}`, max |cosine vs ±KL| `{item['winner_max_abs_cosine_kl']:.4f}`, KL rank `{occurrence.get('rmse_rank', 'not found')}`."
        )
    lines.extend(
        [
            "",
            "## Test B: truncation",
            "",
            f"The original `(5, 20)` mixture restricted to the 80% empirical-mass support closest to the fitted rate used `{test_b['selected_support_size']}` support points from `{test_b['selected_support_min']}` to `{test_b['selected_support_max']}`. Residual RMSE dropped from `{test_b['full_residual_rmse']:.6g}` to `{test_b['residual_rmse']:.6g}`. Winner: `{test_b['winner_expression']}`; max |cosine vs ±KL| `{test_b['winner_max_abs_cosine_kl']:.4f}`; KL occurrence: `{test_b['kl_occurrence']}`.",
            "",
            "## Test C: negative-binomial misspecification",
            "",
            f"NB(mean=10, dispersion=2) fitted by Poisson(lambda=10) gave residual RMSE `{test_c['residual_rmse']:.6g}`. Winner: `{test_c['winner_expression']}`; max |cosine vs ±KL| `{test_c['winner_max_abs_cosine_kl']:.4f}`; KL occurrence: `{test_c['kl_occurrence']}`.",
            "",
            "## Test D: theoretical expansion reference",
            "",
            "For Poisson, the KL/Bregman generator is well-defined in mean-parameter/deviance space: `D(mu||lambda)=mu log(mu/lambda)-mu+lambda`. The SR residual in T2 is instead the pointwise log-density ratio `log q(k)-log p_lambda(k)`. The full distribution KL is the weighted sum of that log-ratio, not the same function as the per-count deviance. This means the predicted KL generator exists, but it is not algebraically identical to the object passed to SR.",
            "",
            f"Raw cosine between pointwise residual and Poisson deviance: `{test_d['raw_cosine_residual_vs_poisson_deviance']:.6g}`. Centered cosine: `{test_d['centered_cosine_residual_vs_poisson_deviance']:.6g}`. Raw cosine against KL on the SR coordinate: `{test_d['raw_cosine_residual_vs_kl_on_sr_x']:.6g}`; centered cosine `{test_d['centered_cosine_residual_vs_kl_on_sr_x']:.6g}`.",
            "",
            "## Verdict",
            "",
        ]
    )
    if clean_kl:
        lines.append(f"Poisson shows KL support under: {', '.join(clean_kl)}.")
    elif deviance_match and sr_coordinate_miss:
        lines.append(
            "No SR follow-up produced a clean KL winner, but Test D shows the empirical pointwise residual strongly agrees with the Poisson deviance in count/lambda space while strongly disagreeing with KL evaluated on the manuscript-style SR coordinate. This is a methodology-scope finding: the Poisson Bregman prediction is visible in the natural deviance coordinate, but the current density-frame SR encoding does not expose that coordinate to the grammar."
        )
    else:
        lines.append(
            "No follow-up produced a clean KL winner or a winner sign-equivalent to KL. The strongest finding is methodological: the Poisson KL prediction lives in mean-parameter/deviance space, while the current SR target is a pointwise log-density residual dominated by tail/log-ratio structure."
        )
    lines.append(f"Final interpretation: `{verdict}`.")
    (OUTDIR / "t2_poisson_exhaustive_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    verdict_payload = {
        "final_interpretation": verdict,
        "clean_kl_support_conditions": clean_kl,
        "test_a_any_kl_support": any(
            item["kl_in_top20_with_bregman_conditions"] or item["winner_max_abs_cosine_kl"] >= 0.95 for item in test_a
        ),
        "test_b_kl_support": bool(test_b["kl_in_top20_with_bregman_conditions"] or test_b["winner_max_abs_cosine_kl"] >= 0.95),
        "test_c_kl_support": bool(test_c["kl_in_top20_with_bregman_conditions"] or test_c["winner_max_abs_cosine_kl"] >= 0.95),
        "test_d": test_d,
        "test_d_deviance_match": deviance_match,
        "test_d_sr_coordinate_miss": sr_coordinate_miss,
    }
    write_json(OUTDIR / "t2_poisson_exhaustive_verdict.json", verdict_payload)
    return verdict_payload


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    test_a = test_a_rate_ratio_sweep()
    test_b = test_b_truncation()
    test_c = test_c_negbinom()
    test_d = test_d_theoretical_expansion()
    verdict = report(test_a, test_b, test_c, test_d)
    print(json.dumps({"outdir": str(OUTDIR), "final_interpretation": verdict["final_interpretation"]}, indent=2))


if __name__ == "__main__":
    main()
