from __future__ import annotations

import csv
import importlib.util
import json
import math
import shutil
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import gaussian_kde


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "t2_redesigned"
SEARCH_PATH = ROOT / "eml_zipf_enriched_search.py"

SEED = 20260419
N_SAMPLES = 10_000
X_LOW = 0.05
X_HIGH = 1.0
BEAM_WIDTH = 50
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
CONSTANT_VARIANCE_THRESHOLD = 1e-10


def load_search_module():
    spec = importlib.util.spec_from_file_location("t2_redesigned_search", SEARCH_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


search = load_search_module()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalized_sr_coordinate(n: int) -> np.ndarray:
    index = np.arange(1, n + 1, dtype=np.float64)
    if n <= 1:
        return np.full(n, X_LOW, dtype=np.float64)
    return X_LOW + (X_HIGH - X_LOW) * np.log(index) / np.log(float(n))


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def fit_noise_diagnostics(residual: np.ndarray) -> dict:
    residual = np.asarray(residual, dtype=np.float64)
    second = np.diff(residual, n=2)
    sign_changes = int(np.sum(np.signbit(second[:-1]) != np.signbit(second[1:]))) if len(second) > 1 else 0
    sign_change_fraction = float(sign_changes / max(len(second) - 1, 1))
    centered = residual - float(np.mean(residual))
    if len(centered) > 1 and np.std(centered[:-1]) > 0 and np.std(centered[1:]) > 0:
        lag1 = float(np.corrcoef(centered[:-1], centered[1:])[0, 1])
    else:
        lag1 = float("nan")
    return {
        "residual_rmse": float(np.sqrt(np.mean(residual * residual))),
        "residual_std": float(np.std(residual)),
        "second_difference_sign_changes": sign_changes,
        "second_difference_sign_change_fraction": sign_change_fraction,
        "lag1_autocorrelation": lag1,
        "noise_dominated_flag": bool(sign_change_fraction > 0.40 and (not math.isfinite(lag1) or lag1 < 0.70)),
    }


def generator_values(kind: str, x: np.ndarray) -> np.ndarray:
    if kind == "gaussian_euclidean":
        return (1.0 - x) ** 2
    if kind == "poisson_kl":
        return x * np.log(x) - x + 1.0
    if kind == "gamma_is":
        return (x - 1.0) - np.log(x)
    raise ValueError(kind)


def is_predicted(kind: str, values: np.ndarray, x: np.ndarray) -> bool:
    return bool(np.max(np.abs(np.asarray(values) - generator_values(kind, x))) <= 1e-8)


PREDICTIONS = {
    "gaussian_euclidean": {
        "expected_generator": "Euclidean",
        "canonical": "(1-x)^2",
        "expected_first_reachable_step": 2,
        "max_steps": 2,
    },
    "poisson_kl": {
        "expected_generator": "generalized KL",
        "canonical": "x*log(x)-x+1",
        "expected_first_reachable_step": 3,
        "max_steps": 3,
    },
    "gamma_is": {
        "expected_generator": "Itakura-Saito",
        "canonical": "(x-1)-log(x)",
        "expected_first_reachable_step": 2,
        "max_steps": 2,
    },
}


def run_search_detailed(kind: str, x: np.ndarray, target: np.ndarray, max_steps: int):
    current_vocabulary = search.initial_vocabulary(x, target)
    steps = []
    global_best = None

    for step in range(1, max_steps + 1):
        generated = search.generate_candidates(current_vocabulary, target, step)
        generated = search.dedupe_candidates(generated)
        generated = search.filter_candidates(generated, CONSTANT_VARIANCE_THRESHOLD)
        generated = sorted(generated, key=lambda item: (item["rmse"], item["expr"]))
        if not generated:
            break

        diverse = search.select_diverse_beam(generated, BEAM_WIDTH, DIVERSITY_WEIGHT)
        step_best = generated[0]
        if global_best is None or (step_best["rmse"], step_best["expr"]) < (global_best["rmse"], global_best["expr"]):
            global_best = step_best

        steps.append(
            {
                "step": step,
                "generated": generated,
                "diverse": diverse,
            }
        )

        if step <= KEEP_ALL_UNTIL_STEP:
            if step < KEEP_ALL_UNTIL_STEP:
                current_vocabulary = current_vocabulary + generated
            else:
                current_vocabulary = diverse
        else:
            current_vocabulary = diverse

    occurrences = []
    for step_payload in steps:
        diverse_exprs = {item["expr"] for item in step_payload["diverse"]}
        for rank, candidate in enumerate(step_payload["generated"], start=1):
            if is_predicted(kind, candidate["values"], x):
                occurrences.append(
                    {
                        "step": int(step_payload["step"]),
                        "rmse_rank": rank,
                        "expr": candidate["expr"],
                        "math": candidate["math"],
                        "rmse": float(candidate["rmse"]),
                        "in_top5_by_rmse": bool(rank <= 5),
                        "in_diversity_beam": bool(candidate["expr"] in diverse_exprs),
                        "diversity_rank": next(
                            (idx for idx, item in enumerate(step_payload["diverse"], start=1) if item["expr"] == candidate["expr"]),
                            None,
                        ),
                    }
                )
                break
    return steps, global_best, occurrences


def top_rows(kind: str, steps: list[dict], x: np.ndarray, mode: str):
    rows = []
    for step_payload in steps:
        candidates = step_payload["generated"][:5] if mode == "rmse" else step_payload["diverse"][:5]
        for rank, candidate in enumerate(candidates, start=1):
            rows.append(
                {
                    "step": int(step_payload["step"]),
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    "is_predicted_generator": bool(is_predicted(kind, candidate["values"], x)),
                }
            )
    return rows


def gaussian_case(rng: np.random.Generator) -> dict:
    components = [
        {"mean": -1.25, "sigma": 1.0, "weight": 0.5},
        {"mean": 1.25, "sigma": 1.0, "weight": 0.5},
    ]
    labels = rng.choice(2, size=N_SAMPLES, p=[0.5, 0.5])
    samples = np.empty(N_SAMPLES, dtype=np.float64)
    for idx, comp in enumerate(components):
        mask = labels == idx
        samples[mask] = rng.normal(comp["mean"], comp["sigma"], size=int(np.sum(mask)))

    mu = float(np.mean(samples))
    sigma = float(np.sqrt(np.mean((samples - mu) ** 2)))
    if not math.isfinite(sigma) or sigma <= 0:
        raise RuntimeError("Gaussian MLE failed: non-positive sigma.")

    lo, hi = np.quantile(samples, [0.005, 0.995])
    grid = np.linspace(float(lo), float(hi), 600)
    kde = gaussian_kde(samples)
    log_emp = np.log(np.clip(kde(grid), 1e-300, None))
    log_fit = -0.5 * ((grid - mu) / sigma) ** 2 - np.log(sigma) - 0.5 * np.log(2.0 * np.pi)
    return {
        "case_id": "t2a_gaussian_euclidean",
        "kind": "gaussian_euclidean",
        "family": "Gaussian",
        "samples": samples,
        "sample_labels": labels,
        "grid": grid,
        "log_empirical": log_emp,
        "log_fitted": log_fit,
        "mixture_params": components,
        "fit": {
            "method": "closed_form_gaussian_mle",
            "converged": True,
            "params": {"mean": mu, "sigma": sigma},
            "at_bounds": False,
        },
        "distinguishability": {
            "mean_gap": abs(components[1]["mean"] - components[0]["mean"]),
            "pooled_sigma": 1.0,
            "separation_in_sigma": abs(components[1]["mean"] - components[0]["mean"]) / 1.0,
        },
    }


def poisson_case(rng: np.random.Generator) -> dict:
    components = [
        {"lambda": 5.0, "weight": 0.5},
        {"lambda": 20.0, "weight": 0.5},
    ]
    labels = rng.choice(2, size=N_SAMPLES, p=[0.5, 0.5])
    samples = np.empty(N_SAMPLES, dtype=np.int64)
    for idx, comp in enumerate(components):
        mask = labels == idx
        samples[mask] = rng.poisson(comp["lambda"], size=int(np.sum(mask)))

    lam = float(np.mean(samples))
    if not math.isfinite(lam) or lam <= 0:
        raise RuntimeError("Poisson MLE failed: non-positive lambda.")

    max_k = max(int(np.max(samples)), int(math.ceil(max(c["lambda"] for c in components) + 8 * math.sqrt(max(c["lambda"] for c in components)))))
    support = np.arange(0, max_k + 1, dtype=np.float64)
    counts = np.bincount(samples, minlength=max_k + 1).astype(np.float64)
    smoothed = ndimage.gaussian_filter1d(counts + 0.5, sigma=1.25, mode="nearest")
    p_emp = smoothed / np.sum(smoothed)
    log_emp = np.log(np.clip(p_emp, 1e-300, None))
    log_fit = support * np.log(lam) - lam - gammaln(support + 1.0)
    return {
        "case_id": "t2b_poisson_kl",
        "kind": "poisson_kl",
        "family": "Poisson",
        "samples": samples.astype(np.float64),
        "sample_labels": labels,
        "grid": support,
        "log_empirical": log_emp,
        "log_fitted": log_fit,
        "mixture_params": components,
        "fit": {
            "method": "closed_form_poisson_mle_plus_smoothed_histogram_density_estimate",
            "converged": True,
            "params": {"lambda": lam},
            "at_bounds": False,
        },
        "distinguishability": {
            "lambda_gap": abs(components[1]["lambda"] - components[0]["lambda"]),
            "sqrt_largest_lambda": math.sqrt(max(c["lambda"] for c in components)),
            "separation_in_sqrt_largest_lambda": abs(components[1]["lambda"] - components[0]["lambda"]) / math.sqrt(max(c["lambda"] for c in components)),
        },
    }


def gamma_neg_loglik(theta: np.ndarray, samples: np.ndarray) -> float:
    log_shape, log_rate = theta
    shape = math.exp(float(log_shape))
    rate = math.exp(float(log_rate))
    return -float(np.sum(shape * np.log(rate) - gammaln(shape) + (shape - 1.0) * np.log(samples) - rate * samples))


def gamma_case(rng: np.random.Generator) -> dict:
    components = [
        {"shape": 0.85, "rate": 1.0, "weight": 0.5},
        {"shape": 4.00, "rate": 1.0, "weight": 0.5},
    ]
    labels = rng.choice(2, size=N_SAMPLES, p=[0.5, 0.5])
    samples = np.empty(N_SAMPLES, dtype=np.float64)
    for idx, comp in enumerate(components):
        mask = labels == idx
        samples[mask] = rng.gamma(shape=comp["shape"], scale=1.0 / comp["rate"], size=int(np.sum(mask)))
    samples = np.clip(samples, 1e-12, None)

    mean = float(np.mean(samples))
    var = float(np.var(samples))
    shape0 = max(mean * mean / max(var, 1e-12), 0.1)
    rate0 = max(mean / max(var, 1e-12), 0.1)
    result = minimize(
        gamma_neg_loglik,
        x0=np.log([shape0, rate0]),
        args=(samples,),
        method="L-BFGS-B",
        bounds=[(np.log(0.05), np.log(20.0)), (np.log(0.05), np.log(20.0))],
        options={"maxiter": 20000},
    )
    if not result.success:
        raise RuntimeError(f"Gamma MLE failed: {result.message}")
    shape, rate = np.exp(result.x)
    shape = float(shape)
    rate = float(rate)
    at_bounds = bool(abs(shape - 0.05) < 1e-6 or abs(shape - 20.0) < 1e-6 or abs(rate - 0.05) < 1e-6 or abs(rate - 20.0) < 1e-6)
    if at_bounds:
        raise RuntimeError("Gamma MLE hit parameter bounds.")

    lo, hi = np.quantile(samples, [0.005, 0.995])
    grid = np.linspace(max(float(lo), 1e-8), float(hi), 600)
    kde = gaussian_kde(samples)
    log_emp = np.log(np.clip(kde(grid), 1e-300, None))
    log_fit = shape * np.log(rate) - gammaln(shape) + (shape - 1.0) * np.log(grid) - rate * grid
    return {
        "case_id": "t2c_gamma_is",
        "kind": "gamma_is",
        "family": "Gamma",
        "samples": samples,
        "sample_labels": labels,
        "grid": grid,
        "log_empirical": log_emp,
        "log_fitted": log_fit,
        "mixture_params": components,
        "fit": {
            "method": "gamma_mle_shape_rate_loc_fixed_zero",
            "converged": True,
            "message": str(result.message),
            "params": {"shape": shape, "rate": rate},
            "at_bounds": at_bounds,
        },
        "distinguishability": {
            "shape_gap": abs(components[1]["shape"] - components[0]["shape"]),
            "shared_rate": 1.0,
        },
    }


def reachability_preflight():
    x = normalized_sr_coordinate(600)
    rows = []
    first_seen = {}
    for kind in PREDICTIONS:
        target = generator_values(kind, x)
        current = search.initial_vocabulary(x, target)
        first_seen[kind] = None
        for step in range(1, 6):
            generated = search.generate_candidates(current, target, step)
            generated = search.dedupe_candidates(generated)
            generated = search.filter_candidates(generated, CONSTANT_VARIANCE_THRESHOLD)
            generated = sorted(generated, key=lambda item: (item["rmse"], item["expr"]))
            for rank, candidate in enumerate(generated, start=1):
                if is_predicted(kind, candidate["values"], x):
                    first_seen[kind] = {
                        "first_reachable_step": step,
                        "rmse_rank_on_generator_target": rank,
                        "expr": candidate["expr"],
                        "math": candidate["math"],
                    }
                    break
            if first_seen[kind] is not None:
                break
            if step < KEEP_ALL_UNTIL_STEP:
                current = current + generated
            else:
                current = search.select_diverse_beam(generated, BEAM_WIDTH, DIVERSITY_WEIGHT)
    for kind, item in first_seen.items():
        rows.append({"kind": kind, **(item or {})})
    unreachable = [kind for kind, item in first_seen.items() if item is None]
    if unreachable:
        raise RuntimeError(f"Predicted generators not reachable by depth 5: {unreachable}")
    return first_seen


def run_case(case: dict):
    case_dir = OUTDIR / case["case_id"]
    case_dir.mkdir(parents=True, exist_ok=True)
    x_sr = normalized_sr_coordinate(len(case["grid"]))
    residual = case["log_empirical"] - case["log_fitted"]
    noise = fit_noise_diagnostics(residual)
    if noise["residual_rmse"] < 1e-4:
        raise RuntimeError(f"{case['case_id']}: residual is trivially small.")
    if noise["noise_dominated_flag"]:
        raise RuntimeError(f"{case['case_id']}: residual appears noise dominated: {noise}.")

    max_steps = PREDICTIONS[case["kind"]]["max_steps"]
    steps, global_best, occurrences = run_search_detailed(case["kind"], x_sr, residual, max_steps)

    sample_rows = [
        {
            "sample_index": idx,
            "value": f"{float(value):.17g}",
            "component": int(label),
        }
        for idx, (value, label) in enumerate(zip(case["samples"], case["sample_labels"]), start=1)
    ]
    write_csv(case_dir / "synthetic_samples.csv", sample_rows, ["sample_index", "value", "component"])

    density_rows = [
        {
            "grid_index": idx,
            "domain_value": f"{float(v):.17g}",
            "sr_x": f"{float(x):.17g}",
            "log_p_emp": f"{float(emp):.17g}",
            "log_p_fit": f"{float(fit):.17g}",
            "residual": f"{float(res):.17g}",
        }
        for idx, (v, x, emp, fit, res) in enumerate(
            zip(case["grid"], x_sr, case["log_empirical"], case["log_fitted"], residual),
            start=1,
        )
    ]
    write_csv(case_dir / "kde_vs_fitted.csv", density_rows, ["grid_index", "domain_value", "sr_x", "log_p_emp", "log_p_fit", "residual"])

    write_csv(case_dir / "per_step_top5_by_rmse.csv", top_rows(case["kind"], steps, x_sr, "rmse"), ["step", "rank", "expression", "math", "rmse", "is_predicted_generator"])
    write_csv(case_dir / "per_step_top5_by_diversity.csv", top_rows(case["kind"], steps, x_sr, "diversity"), ["step", "rank", "expression", "math", "rmse", "is_predicted_generator"])

    predicted_first_step = PREDICTIONS[case["kind"]]["expected_first_reachable_step"]
    predicted_at_first = [item for item in occurrences if int(item["step"]) == predicted_first_step]
    pass_top5 = any(item["in_top5_by_rmse"] or (item["diversity_rank"] is not None and item["diversity_rank"] <= 5) for item in predicted_at_first)
    global_payload = None
    if global_best is not None:
        global_payload = {
            "step": int(global_best["step"]),
            "expr": global_best["expr"],
            "math": global_best["math"],
            "rmse": float(global_best["rmse"]),
            "is_predicted_generator": bool(is_predicted(case["kind"], global_best["values"], x_sr)),
        }

    fit_payload = {
        "case_id": case["case_id"],
        "family": case["family"],
        "prediction": PREDICTIONS[case["kind"]],
        "mixture_params": case["mixture_params"],
        "distinguishability": case["distinguishability"],
        "fit": case["fit"],
        "density_estimation": {
            "continuous": case["family"] in {"Gaussian", "Gamma"},
            "method": "gaussian_kde_on_iid_samples" if case["family"] in {"Gaussian", "Gamma"} else "gaussian_smoothed_histogram_with_0.5_pseudocount",
            "grid_points": int(len(case["grid"])),
        },
        "residual_diagnostics": noise,
        "search_config": {
            "max_steps": max_steps,
            "beam_width": BEAM_WIDTH,
            "keep_all_until_step": KEEP_ALL_UNTIL_STEP,
            "diversity_weight": DIVERSITY_WEIGHT,
            "semantic_hash_rounding_digits": search.SIG_ROUND,
            "exp_clamp": search.EXP_CLAMP,
            "value_abs_limit": search.VALUE_ABS_LIMIT,
            "unary_ops": [op["name"] for op in search.UNARY_OPS],
            "binary_ops": [op["name"] for op in search.BINARY_OPS],
        },
        "sr_result": {
            "predicted_occurrences": occurrences,
            "predicted_at_first_reachable_step": predicted_at_first,
            "pass_top5_rmse_or_diversity_at_first_step": bool(pass_top5),
            "global_best": global_payload,
        },
    }
    (case_dir / "fitted_parameters.json").write_text(json.dumps(fit_payload, indent=2) + "\n", encoding="utf-8")

    if pass_top5:
        verdict = "PASS"
    elif predicted_at_first:
        verdict = "FAIL_PREDICTED_PRESENT_BUT_NOT_TOP5"
    elif occurrences:
        verdict = "FAIL_PREDICTED_APPEARS_AFTER_FIRST_STEP_ONLY"
    else:
        verdict = "FAIL_PREDICTED_NOT_OBSERVED"

    report = [
        f"# {case['case_id']}",
        "",
        f"- family: `{case['family']}`",
        f"- predicted generator: `{PREDICTIONS[case['kind']]['expected_generator']}` (`{PREDICTIONS[case['kind']]['canonical']}`)",
        f"- search depth: `{max_steps}`",
        f"- residual RMSE: `{noise['residual_rmse']:.12g}`",
        f"- residual lag-1 autocorrelation: `{noise['lag1_autocorrelation']:.6g}`",
        f"- residual second-difference sign-change fraction: `{noise['second_difference_sign_change_fraction']:.6g}`",
        f"- predicted occurrences: `{occurrences}`",
        f"- global best: `{global_payload}`",
        f"- verdict: `{verdict}`",
        "",
    ]
    if verdict != "PASS":
        report.append("The predicted generator did not appear in the top-5 by RMSE or top-5 by diversity at its first reachable step under this density-estimation protocol.")
    else:
        report.append("The predicted generator appears in the requested top-5 criterion at its first reachable step.")
    (case_dir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    return {
        "case_id": case["case_id"],
        "family": case["family"],
        "predicted_generator": PREDICTIONS[case["kind"]]["expected_generator"],
        "max_steps": max_steps,
        "verdict": verdict,
        "pass": bool(pass_top5),
        "predicted_occurrences": occurrences,
        "global_best": global_payload,
        "residual_diagnostics": noise,
    }


def write_methodology_note(reachability: dict):
    lines = [
        "# T2 Redesigned Methodology Note",
        "",
        "This T2 restart uses the density-estimation frame rather than the rank-frequency frame. Each case generates `N=10000` iid samples from a two-component mixture, fits a single-family density by maximum likelihood, estimates the empirical log density on a grid, and runs the Section 2.4 deterministic search on `log p_empirical - log p_fitted`.",
        "",
        "The SR coordinate is a monotone log-index coordinate over the density-evaluation grid: `x = 0.05 + 0.95*log(i)/log(M)` for grid index `i=1..M`. This preserves the manuscript search domain convention while acknowledging that density-estimation support is not word-rank.",
        "",
        "The live Section 2.4 implementation includes unary operators `neg, inv, sqr, sqrt, exp, log` and binary operators `eml, add, sub, mul, div, pow`. Beam width is 50, semantic hashing rounds output vectors to ten decimals, `EXP_CLAMP=30`, and `VALUE_ABS_LIMIT=1e6`.",
        "",
        "## Reachability",
        "",
    ]
    for kind, item in reachability.items():
        lines.append(f"- `{kind}` first reachable at step `{item['first_reachable_step']}` as `{item['expr']}` ({item['math']}).")
    (OUTDIR / "t2_redesigned_methodology_note.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for child in ["t2a_gaussian_euclidean", "t2b_poisson_kl", "t2c_gamma_is"]:
        path = OUTDIR / child
        if path.exists():
            shutil.rmtree(path)

    rng = np.random.default_rng(SEED)
    reachability = reachability_preflight()
    write_methodology_note(reachability)
    (OUTDIR / "t2_redesigned_reachability.json").write_text(json.dumps(reachability, indent=2) + "\n", encoding="utf-8")

    cases = [gaussian_case(rng), poisson_case(rng), gamma_case(rng)]
    summaries = [run_case(case) for case in cases]
    pass_count = sum(1 for item in summaries if item["pass"])
    aggregate = {
        "status": "complete",
        "seed": SEED,
        "sample_count_per_experiment": N_SAMPLES,
        "pass_count": pass_count,
        "case_count": len(summaries),
        "overall_verdict": "strong_support" if pass_count == 3 else "partial_support" if pass_count == 2 else "no_support_under_redesigned_protocol",
        "cases": summaries,
    }
    (OUTDIR / "t2_redesigned_aggregate.json").write_text(json.dumps(aggregate, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# T2 Redesigned Report",
        "",
        f"Density-frame T2 completed for `{len(summaries)}` families. Predicted generators passed the top-5 criterion in `{pass_count}/{len(summaries)}` cases.",
        "",
    ]
    for item in summaries:
        lines.extend(
            [
                f"## {item['case_id']}",
                "",
                f"- family: `{item['family']}`",
                f"- predicted generator: `{item['predicted_generator']}`",
                f"- verdict: `{item['verdict']}`",
                f"- global best: `{item['global_best']}`",
                f"- predicted occurrences: `{item['predicted_occurrences']}`",
                "",
            ]
        )
    if pass_count == 3:
        lines.append("All three families pass the redesigned density-frame criterion.")
    elif pass_count == 2:
        lines.append("Two of three families pass; the result is partial support and the failed family should be inspected.")
    else:
        lines.append("The redesigned density-frame test does not support the cross-family generator-identity claim as stated.")
    (OUTDIR / "t2_redesigned_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({"outdir": str(OUTDIR), "pass_count": pass_count, "case_count": len(summaries)}, indent=2))


if __name__ == "__main__":
    main()
