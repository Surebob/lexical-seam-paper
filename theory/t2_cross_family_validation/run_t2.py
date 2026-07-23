from __future__ import annotations

import csv
import importlib.util
import json
import math
import shutil
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.special import gammaln


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "t2_cross_family_validation"
SEARCH_PATH = ROOT / "eml_zipf_enriched_search.py"

BEAM_WIDTH = 50
MAX_STEPS = 4
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
CONSTANT_VARIANCE_THRESHOLD = 1e-10
X_LOW = 0.05
X_HIGH = 1.0

PREDICTED_EXPRESSIONS = {
    "gaussian_euclidean": {
        "label": "Euclidean",
        "expected_first_step": 2,
        "canonical": "(1-x)^2",
    },
    "poisson_kl": {
        "label": "generalized KL",
        "expected_first_step": 3,
        "canonical": "x*log(x)-x+1",
    },
    "gamma_is": {
        "label": "Itakura-Saito",
        "expected_first_step": 2,
        "canonical": "(x-1)-log(x)",
    },
}


def load_search_module():
    spec = importlib.util.spec_from_file_location("t2_enriched_search", SEARCH_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


search = load_search_module()


def normalize_rank_axis(vocab_size: int) -> np.ndarray:
    ranks = np.arange(1, vocab_size + 1, dtype=np.float64)
    return X_LOW + (X_HIGH - X_LOW) * np.log(ranks) / np.log(float(vocab_size))


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def safe_log(values: np.ndarray) -> np.ndarray:
    return np.log(np.clip(np.asarray(values, dtype=np.float64), 1e-300, None))


def sort_desc(values: np.ndarray):
    order = np.argsort(np.asarray(values))[::-1]
    return np.asarray(values, dtype=np.float64)[order], order


def gaussian_log_pdf(t: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    return -0.5 * ((t - mu) / sigma) ** 2 - np.log(sigma) - 0.5 * np.log(2.0 * np.pi)


def poisson_log_pmf(k: np.ndarray, lam: float) -> np.ndarray:
    return k * np.log(lam) - lam - gammaln(k + 1.0)


def gamma_log_pdf(t: np.ndarray, shape: float, rate: float) -> np.ndarray:
    return shape * np.log(rate) - gammaln(shape) + (shape - 1.0) * np.log(t) - rate * t


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def candidate_values(expr: dict) -> np.ndarray:
    return np.asarray(expr["values"], dtype=np.float64)


def predicted_values(kind: str, x: np.ndarray) -> np.ndarray:
    if kind == "gaussian_euclidean":
        return (1.0 - x) ** 2
    if kind == "poisson_kl":
        return x * np.log(x) - x + 1.0
    if kind == "gamma_is":
        return (x - 1.0) - np.log(x)
    raise ValueError(kind)


def is_predicted_candidate(kind: str, values: np.ndarray, x: np.ndarray) -> bool:
    target = predicted_values(kind, x)
    return bool(np.max(np.abs(np.asarray(values) - target)) <= 1e-8)


def classify_expression(expr: str, values: np.ndarray, x: np.ndarray) -> str:
    for kind, info in PREDICTED_EXPRESSIONS.items():
        if is_predicted_candidate(kind, values, x):
            return info["label"]
    if expr in {"1", "x"}:
        return "trivial"
    return "other"


def run_search_with_values(x: np.ndarray, target: np.ndarray):
    current_vocabulary = search.initial_vocabulary(x, target)
    steps = []
    global_best = None

    for step in range(1, MAX_STEPS + 1):
        generated = search.generate_candidates(current_vocabulary, target, step)
        generated = search.dedupe_candidates(generated)
        generated = search.filter_candidates(generated, CONSTANT_VARIANCE_THRESHOLD)
        generated = sorted(generated, key=lambda item: (item["rmse"], item["expr"]))
        if not generated:
            break

        step_best = generated[0]
        if global_best is None or (step_best["rmse"], step_best["expr"]) < (global_best["rmse"], global_best["expr"]):
            global_best = step_best

        steps.append({"step": step, "candidates": generated})

        if step <= KEEP_ALL_UNTIL_STEP:
            if step < KEEP_ALL_UNTIL_STEP:
                current_vocabulary = current_vocabulary + generated
            else:
                current_vocabulary = search.select_diverse_beam(generated, BEAM_WIDTH, DIVERSITY_WEIGHT)
        else:
            current_vocabulary = search.select_diverse_beam(generated, BEAM_WIDTH, DIVERSITY_WEIGHT)

    return steps, global_best


def predicted_rank_by_step(kind: str, steps: list[dict], x: np.ndarray):
    found = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["candidates"], start=1):
            if is_predicted_candidate(kind, candidate_values(candidate), x):
                found.append(
                    {
                        "step": step_payload["step"],
                        "beam_rank": rank,
                        "expr": candidate["expr"],
                        "math": candidate["math"],
                        "rmse": float(candidate["rmse"]),
                    }
                )
                break
    return found


def curvature_diagnostic(log_frequency: np.ndarray) -> dict:
    second = np.diff(np.asarray(log_frequency, dtype=np.float64), n=2)
    return {
        "max_abs_second_difference": float(np.max(np.abs(second))),
        "mean_abs_second_difference": float(np.mean(np.abs(second))),
        "sign_changes_second_difference": int(np.sum(np.signbit(second[:-1]) != np.signbit(second[1:]))),
    }


def check_not_degenerate(name: str, component_a: np.ndarray, component_b: np.ndarray, log_frequency: np.ndarray):
    corr = float(np.corrcoef(component_a, component_b)[0, 1])
    rel_l2 = float(np.linalg.norm(component_a - component_b) / max(np.linalg.norm(component_a), 1e-300))
    curv = curvature_diagnostic(log_frequency)
    if rel_l2 < 0.05 or abs(corr) > 0.9999:
        raise RuntimeError(f"{name}: synthetic components are degenerate (corr={corr}, rel_l2={rel_l2}).")
    if curv["max_abs_second_difference"] < 1e-8:
        raise RuntimeError(f"{name}: rank-frequency curve lacks visible curvature.")
    return {"component_corr": corr, "component_relative_l2": rel_l2, **curv}


def build_gaussian_case():
    t = np.linspace(-5.0, 5.0, 1400)
    params = {
        "component_1": {"mu": -1.20, "sigma": 0.65, "weight": 0.5},
        "component_2": {"mu": 1.35, "sigma": 1.15, "weight": 0.5},
    }
    comp1 = params["component_1"]["weight"] * np.exp(gaussian_log_pdf(t, params["component_1"]["mu"], params["component_1"]["sigma"]))
    comp2 = params["component_2"]["weight"] * np.exp(gaussian_log_pdf(t, params["component_2"]["mu"], params["component_2"]["sigma"]))
    freqs_unsorted = comp1 + comp2
    freqs, order = sort_desc(freqs_unsorted)
    t_sorted = t[order]
    y = safe_log(freqs)

    def residual(theta):
        log_amp, mu, log_sigma = theta
        sigma = np.exp(log_sigma)
        model_unsorted = log_amp + gaussian_log_pdf(t, mu, sigma)
        return np.sort(model_unsorted)[::-1] - y

    result = least_squares(
        residual,
        x0=np.array([0.0, 0.0, np.log(1.2)]),
        bounds=([-20.0, -8.0, np.log(0.15)], [20.0, 8.0, np.log(5.0)]),
        max_nfev=20000,
    )
    if not result.success:
        raise RuntimeError(f"T2a Gaussian fit failed: {result.message}")
    log_amp, mu, log_sigma = result.x
    sigma = float(np.exp(log_sigma))
    pred_unsorted = log_amp + gaussian_log_pdf(t, mu, sigma)
    fitted = np.sort(pred_unsorted)[::-1]
    diagnostic = check_not_degenerate("T2a Gaussian", comp1[order], comp2[order], y)
    return {
        "kind": "gaussian_euclidean",
        "case_id": "t2a_gaussian_euclidean",
        "title": "T2a Gaussian -> Euclidean",
        "support": t_sorted,
        "frequency": freqs,
        "log_frequency": y,
        "fitted_log_frequency": fitted,
        "fit": {
            "method": "nonlinear_least_squares_on_sorted_log_density",
            "success": bool(result.success),
            "message": result.message,
            "cost": float(result.cost),
            "params": {"log_amp": float(log_amp), "mu": float(mu), "sigma": sigma},
            "at_bounds": {
                "mu": bool(abs(mu - (-8.0)) < 1e-6 or abs(mu - 8.0) < 1e-6),
                "sigma": bool(abs(sigma - 0.15) < 1e-6 or abs(sigma - 5.0) < 1e-6),
            },
        },
        "mixture_params": params,
        "diagnostic": diagnostic,
    }


def build_poisson_case():
    k = np.arange(0, 90, dtype=np.float64)
    params = {
        "component_1": {"lambda": 6.0, "weight": 0.5},
        "component_2": {"lambda": 24.0, "weight": 0.5},
    }
    comp1 = params["component_1"]["weight"] * np.exp(poisson_log_pmf(k, params["component_1"]["lambda"]))
    comp2 = params["component_2"]["weight"] * np.exp(poisson_log_pmf(k, params["component_2"]["lambda"]))
    freqs_unsorted = comp1 + comp2
    freqs, order = sort_desc(freqs_unsorted)
    k_sorted = k[order]
    y = safe_log(freqs)

    def residual(theta):
        log_amp, log_lam = theta
        lam = np.exp(log_lam)
        model_unsorted = log_amp + poisson_log_pmf(k, lam)
        return np.sort(model_unsorted)[::-1] - y

    result = least_squares(
        residual,
        x0=np.array([0.0, np.log(15.0)]),
        bounds=([-20.0, np.log(0.2)], [20.0, np.log(80.0)]),
        max_nfev=20000,
    )
    if not result.success:
        raise RuntimeError(f"T2b Poisson fit failed: {result.message}")
    log_amp, log_lam = result.x
    lam = float(np.exp(log_lam))
    pred_unsorted = log_amp + poisson_log_pmf(k, lam)
    fitted = np.sort(pred_unsorted)[::-1]
    diagnostic = check_not_degenerate("T2b Poisson", comp1[order], comp2[order], y)
    mle_lambda = params["component_1"]["weight"] * params["component_1"]["lambda"] + params["component_2"]["weight"] * params["component_2"]["lambda"]
    return {
        "kind": "poisson_kl",
        "case_id": "t2b_poisson_kl",
        "title": "T2b Poisson -> generalized KL",
        "support": k_sorted,
        "frequency": freqs,
        "log_frequency": y,
        "fitted_log_frequency": fitted,
        "fit": {
            "method": "nonlinear_least_squares_on_sorted_log_pmf",
            "success": bool(result.success),
            "message": result.message,
            "cost": float(result.cost),
            "params": {"log_amp": float(log_amp), "lambda": lam, "mixture_mean_mle_lambda": float(mle_lambda)},
            "at_bounds": {"lambda": bool(abs(lam - 0.2) < 1e-6 or abs(lam - 80.0) < 1e-6)},
        },
        "mixture_params": params,
        "diagnostic": diagnostic,
    }


def build_gamma_case():
    t = np.linspace(0.02, 18.0, 1400)
    params = {
        "component_1": {"shape": 0.70, "rate": 1.05, "weight": 0.5},
        "component_2": {"shape": 3.60, "rate": 1.05, "weight": 0.5},
    }
    comp1 = params["component_1"]["weight"] * np.exp(gamma_log_pdf(t, params["component_1"]["shape"], params["component_1"]["rate"]))
    comp2 = params["component_2"]["weight"] * np.exp(gamma_log_pdf(t, params["component_2"]["shape"], params["component_2"]["rate"]))
    freqs_unsorted = comp1 + comp2
    weights = freqs_unsorted / np.sum(freqs_unsorted)
    freqs, order = sort_desc(freqs_unsorted)
    t_sorted = t[order]
    y = safe_log(freqs)

    def nll(theta):
        log_shape, log_rate = theta
        shape = np.exp(log_shape)
        rate = np.exp(log_rate)
        log_pdf = gamma_log_pdf(t, shape, rate)
        return -float(np.sum(weights * log_pdf))

    result = minimize(
        nll,
        x0=np.array([np.log(1.6), np.log(1.0)]),
        bounds=[(np.log(0.15), np.log(10.0)), (np.log(0.05), np.log(10.0))],
        method="L-BFGS-B",
        options={"maxiter": 20000},
    )
    if not result.success:
        raise RuntimeError(f"T2c Gamma fit failed: {result.message}")
    log_shape, log_rate = result.x
    shape = float(np.exp(log_shape))
    rate = float(np.exp(log_rate))
    pred_unsorted_no_amp = gamma_log_pdf(t, shape, rate)
    pred_sorted_no_amp = np.sort(pred_unsorted_no_amp)[::-1]
    log_amp = float(np.mean(y - pred_sorted_no_amp))
    fitted = log_amp + pred_sorted_no_amp
    diagnostic = check_not_degenerate("T2c Gamma", comp1[order], comp2[order], y)
    return {
        "kind": "gamma_is",
        "case_id": "t2c_gamma_is",
        "title": "T2c Gamma -> Itakura-Saito",
        "support": t_sorted,
        "frequency": freqs,
        "log_frequency": y,
        "fitted_log_frequency": fitted,
        "fit": {
            "method": "expected_mle_on_unsorted_gamma_grid_plus_log_amplitude_alignment",
            "success": bool(result.success),
            "message": str(result.message),
            "objective": float(result.fun),
            "params": {"log_amp": log_amp, "shape": shape, "rate": rate},
            "at_bounds": {
                "shape": bool(abs(shape - 0.15) < 1e-6 or abs(shape - 10.0) < 1e-6),
                "rate": bool(abs(rate - 0.05) < 1e-6 or abs(rate - 10.0) < 1e-6),
            },
        },
        "mixture_params": params,
        "diagnostic": diagnostic,
    }


def top5_rows(case: dict, steps: list[dict], x: np.ndarray):
    rows = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["candidates"][:5], start=1):
            rows.append(
                {
                    "step": step_payload["step"],
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    "classification": classify_expression(candidate["expr"], candidate_values(candidate), x),
                }
            )
    return rows


def write_case_outputs(case: dict):
    case_dir = OUTDIR / case["case_id"]
    case_dir.mkdir(parents=True, exist_ok=True)
    x = normalize_rank_axis(len(case["frequency"]))
    residual = case["log_frequency"] - case["fitted_log_frequency"]
    residual_rmse = float(np.sqrt(np.mean(residual * residual)))

    steps, global_best = run_search_with_values(x, residual)
    predicted_found = predicted_rank_by_step(case["kind"], steps, x)
    predicted_best = min(predicted_found, key=lambda item: (item["rmse"], item["step"], item["beam_rank"])) if predicted_found else None
    predicted_retained_top50 = [item for item in predicted_found if int(item["beam_rank"]) <= BEAM_WIDTH]
    predicted_top5 = [item for item in predicted_found if int(item["beam_rank"]) <= 5]
    predicted_wins = bool(
        predicted_best is not None
        and global_best is not None
        and abs(float(predicted_best["rmse"]) - float(global_best["rmse"])) <= 1e-10
    )

    synthetic_rows = []
    for rank, (xi, support_value, freq, y, fitted, resid) in enumerate(
        zip(x, case["support"], case["frequency"], case["log_frequency"], case["fitted_log_frequency"], residual),
        start=1,
    ):
        synthetic_rows.append(
            {
                "rank": rank,
                "x": f"{float(xi):.17g}",
                "support_value_sorted_by_frequency": f"{float(support_value):.17g}",
                "frequency": f"{float(freq):.17g}",
                "log_frequency": f"{float(y):.17g}",
                "fitted_log_frequency": f"{float(fitted):.17g}",
                "residual": f"{float(resid):.17g}",
            }
        )
    write_csv(
        case_dir / "synthetic_data.csv",
        synthetic_rows,
        ["rank", "x", "support_value_sorted_by_frequency", "frequency", "log_frequency", "fitted_log_frequency", "residual"],
    )

    beam_rows = top5_rows(case, steps, x)
    write_csv(case_dir / "per_step_beam.csv", beam_rows, ["step", "rank", "expression", "math", "rmse", "classification"])

    fit_payload = {
        "case_id": case["case_id"],
        "title": case["title"],
        "mixture_params": case["mixture_params"],
        "fit": case["fit"],
        "residual_rmse": residual_rmse,
        "diagnostic": case["diagnostic"],
        "search_config": {
            "beam_width": BEAM_WIDTH,
            "max_steps": MAX_STEPS,
            "keep_all_until_step": KEEP_ALL_UNTIL_STEP,
            "diversity_weight": DIVERSITY_WEIGHT,
            "constant_variance_threshold": CONSTANT_VARIANCE_THRESHOLD,
            "exp_clamp": search.EXP_CLAMP,
            "value_abs_limit": search.VALUE_ABS_LIMIT,
            "semantic_hash_rounding_digits": search.SIG_ROUND,
            "unary_ops": [op["name"] for op in search.UNARY_OPS],
            "binary_ops": [op["name"] for op in search.BINARY_OPS],
        },
    }
    (case_dir / "single_fit_parameters.json").write_text(json.dumps(fit_payload, indent=2) + "\n", encoding="utf-8")

    global_best_payload = None
    if global_best is not None:
        global_best_payload = {
            "step": int(global_best["step"]),
            "expr": global_best["expr"],
            "math": global_best["math"],
            "rmse": float(global_best["rmse"]),
            "classification": classify_expression(global_best["expr"], candidate_values(global_best), x),
        }

    pass_fail = {
        "predicted_generator": PREDICTED_EXPRESSIONS[case["kind"]]["label"],
        "predicted_canonical": PREDICTED_EXPRESSIONS[case["kind"]]["canonical"],
        "expected_first_reachable_step": PREDICTED_EXPRESSIONS[case["kind"]]["expected_first_step"],
        "predicted_generated": bool(predicted_found),
        "predicted_first_generated": predicted_found[0] if predicted_found else None,
        "predicted_best_generated": predicted_best,
        "predicted_retained_in_top50_rmse_beam": bool(predicted_retained_top50),
        "predicted_first_retained_top50": predicted_retained_top50[0] if predicted_retained_top50 else None,
        "predicted_in_reported_top5": bool(predicted_top5),
        "predicted_first_top5": predicted_top5[0] if predicted_top5 else None,
        "global_best": global_best_payload,
        "predicted_wins_by_global_rmse": predicted_wins,
        "rmse_gap_global_minus_predicted": None
        if predicted_best is None or global_best is None
        else float(global_best["rmse"] - predicted_best["rmse"]),
    }

    if predicted_wins:
        verdict = "PASS"
    elif predicted_retained_top50:
        verdict = "STOP_CONDITION_PREDICTED_APPEARS_BUT_DOES_NOT_WIN"
    elif predicted_found:
        verdict = "FAIL_PREDICTED_GENERATED_BUT_RANKED_OUTSIDE_TOP50"
    else:
        verdict = "FAIL_PREDICTED_NOT_IN_BEAM"
    pass_fail["verdict"] = verdict

    report_lines = [
        f"# {case['title']}",
        "",
        f"- predicted generator: `{pass_fail['predicted_generator']}` (`{pass_fail['predicted_canonical']}`)",
        f"- residual RMSE after single-family fit: `{residual_rmse:.12g}`",
        f"- mixture diagnostic max abs second difference: `{case['diagnostic']['max_abs_second_difference']:.12g}`",
        f"- predicted generated by depth 4: `{pass_fail['predicted_generated']}`",
        f"- predicted retained in top-50 RMSE beam: `{pass_fail['predicted_retained_in_top50_rmse_beam']}`",
        f"- predicted appears in reported top-5 rows: `{pass_fail['predicted_in_reported_top5']}`",
        f"- predicted wins by global RMSE through depth 4: `{pass_fail['predicted_wins_by_global_rmse']}`",
        f"- verdict: `{verdict}`",
        "",
    ]
    if predicted_found:
        first = predicted_found[0]
        report_lines.extend(
            [
                f"The predicted generator is first generated at step `{first['step']}` with RMSE rank `{first['beam_rank']}` and RMSE `{first['rmse']:.12g}`.",
                "",
            ]
        )
    if global_best_payload:
        report_lines.extend(
            [
                f"The global best expression through depth 4 is `{global_best_payload['expr']}` at step `{global_best_payload['step']}` with RMSE `{global_best_payload['rmse']:.12g}`.",
                "",
            ]
        )
    if verdict.startswith("STOP"):
        report_lines.append("Stop condition triggered: the predicted generator appears in the beam but is not the lowest-RMSE expression through depth 4.")
    elif verdict == "FAIL_PREDICTED_GENERATED_BUT_RANKED_OUTSIDE_TOP50":
        report_lines.append("The predicted generator is generated, but it ranks outside the top-50 RMSE beam and is not competitive under this residual.")
    elif verdict == "FAIL_PREDICTED_NOT_IN_BEAM":
        report_lines.append("The predicted generator did not appear in the retained beam through depth 4 under the live search protocol.")
    else:
        report_lines.append("The predicted generator is the lowest-RMSE expression observed through depth 4.")
    (case_dir / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    return {
        "case_id": case["case_id"],
        "title": case["title"],
        "pass_fail": pass_fail,
        "residual_rmse": residual_rmse,
        "diagnostic": case["diagnostic"],
        "output_dir": str(case_dir.relative_to(ROOT)),
    }


def write_methodology_note():
    note = [
        "# T2 Methodology Note",
        "",
        "T2 uses the live Section 2.4 enumerative search implementation from `eml_zipf_enriched_search.py` with `beam_width=50`, `max_steps=4`, `keep_all_until_step=2`, semantic hashing rounded to ten decimals, `EXP_CLAMP=30`, and `VALUE_ABS_LIMIT=1e6`.",
        "",
        "The depth-4 revision is required because the Poisson generalized-KL generator `x*log(x)-x+1` is not reachable at step 2. It first becomes operationally reachable at step 3 as the equivalent expression `sub[add[1,neg[x]],mul[neg[x],log[x]]]` when the target residual is exactly KL-shaped. Euclidean and Itakura-Saito are reachable at step 2.",
        "",
        "For each synthetic family, the generated frequencies are sorted descending to form a rank-frequency curve, the corresponding single-family curve is fitted, and the residual `log_frequency - fitted_log_frequency` is passed to the parameter-free search on the normalized rank axis `x = 0.05 + 0.95*log(r)/log(V)`.",
    ]
    (OUTDIR / "t2_methodology_note.md").write_text("\n".join(note) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for old in ["t2a_gaussian_euclidean", "t2b_poisson_kl", "t2c_gamma_is"]:
        path = OUTDIR / old
        if path.exists():
            shutil.rmtree(path)

    write_methodology_note()
    cases = [build_gaussian_case(), build_poisson_case(), build_gamma_case()]
    summaries = []
    for case in cases:
        summaries.append(write_case_outputs(case))

    pass_count = sum(1 for item in summaries if item["pass_fail"]["verdict"] == "PASS")
    aggregate = {
        "status": "complete",
        "max_steps": MAX_STEPS,
        "beam_width": BEAM_WIDTH,
        "pass_count": pass_count,
        "case_count": len(summaries),
        "overall_verdict": "strong_support" if pass_count == 3 else "partial_or_no_support",
        "cases": summaries,
    }
    (OUTDIR / "t2_aggregate_summary.json").write_text(json.dumps(aggregate, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# T2 Cross-Family Synthetic Validation",
        "",
        f"Depth-4 search completed for `{len(summaries)}` synthetic families. Predicted generator wins: `{pass_count}/{len(summaries)}`.",
        "",
    ]
    for item in summaries:
        pf = item["pass_fail"]
        lines.extend(
            [
                f"## {item['title']}",
                "",
                f"- predicted: `{pf['predicted_generator']}`",
                f"- observed global best: `{pf['global_best']['expr'] if pf['global_best'] else None}`",
                f"- predicted first generated: `{pf['predicted_first_generated']}`",
                f"- predicted first retained top-50: `{pf['predicted_first_retained_top50']}`",
                f"- verdict: `{pf['verdict']}`",
                "",
            ]
        )
    if pass_count == 3:
        lines.append("All three synthetic families selected the predicted generator once reachable, supporting the cross-family generator-identity claim under this depth-4 protocol.")
    elif pass_count == 2:
        lines.append("Two of three synthetic families selected the predicted generator. This is partial support; the failing family should be inspected before the cross-family claim is strengthened.")
    else:
        lines.append("Zero or one synthetic family selected the predicted generator. This does not support the cross-family generalization as stated under this protocol.")
    (OUTDIR / "t2_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({"outdir": str(OUTDIR), "pass_count": pass_count, "case_count": len(summaries)}, indent=2))


if __name__ == "__main__":
    main()
