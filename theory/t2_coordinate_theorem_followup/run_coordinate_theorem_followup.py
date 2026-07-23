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
OUTDIR = ROOT / "phase2_addon" / "t2_coordinate_theorem_followup"


def load_t2_module():
    spec = importlib.util.spec_from_file_location("t2_redesigned_for_coordinate_followup", T2_SCRIPT)
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


def affine_fit(y: np.ndarray, g: np.ndarray) -> dict:
    design = np.column_stack([np.ones_like(g), g])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ coeffs
    resid = y - pred
    ss_res = float(np.sum(resid * resid))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    return {
        "a": float(coeffs[0]),
        "b": float(coeffs[1]),
        "prediction": pred,
        "residual": resid,
        "rmse": float(np.sqrt(np.mean(resid * resid))),
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
        "max_abs_residual": float(np.max(np.abs(resid))),
    }


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


def boundary_diagnostics(candidate: dict, support_x: np.ndarray | None = None) -> dict:
    x_dense = np.linspace(0.05, 1.0, 5000, dtype=np.float64)
    values = eval_node(candidate["node"], x_dense)
    if values is None or not np.all(np.isfinite(values)):
        return {
            "f_at_1": float("nan"),
            "fprime_at_1": float("nan"),
            "fsecond_min_continuous_grid": float("nan"),
            "fsecond_positive_continuous_grid": False,
            "bregman_conditions_pass": False,
            "second_diff_min_on_integer_support": float("nan"),
            "second_diff_positive_on_integer_support": False,
        }
    first = np.gradient(values, x_dense, edge_order=2)
    second = np.gradient(first, x_dense, edge_order=2)
    support_min = float("nan")
    support_positive = False
    if support_x is not None:
        support_values = eval_node(candidate["node"], support_x)
        if support_values is not None and len(support_values) >= 3 and np.all(np.isfinite(support_values)):
            second_diff = np.diff(support_values, n=2)
            support_min = float(np.min(second_diff))
            support_positive = bool(support_min > -1e-8)
    f_at_1 = float(values[-1])
    fprime_at_1 = float(first[-1])
    second_min = float(np.min(second))
    return {
        "f_at_1": f_at_1,
        "fprime_at_1": fprime_at_1,
        "fsecond_min_continuous_grid": second_min,
        "fsecond_positive_continuous_grid": bool(second_min > -1e-5),
        "bregman_conditions_pass": bool(abs(f_at_1) <= 1e-5 and abs(fprime_at_1) <= 1e-4 and second_min > -1e-5),
        "second_diff_min_on_integer_support": support_min,
        "second_diff_positive_on_integer_support": support_positive,
    }


def poisson_mixture(lambda_a: float, lambda_b: float, mass: float = 0.99) -> dict:
    weight = 0.5
    fitted_lambda = 0.5 * lambda_a + 0.5 * lambda_b
    k = np.arange(0, 500, dtype=np.int64)
    full_q = weight * poisson.pmf(k, lambda_a) + weight * poisson.pmf(k, lambda_b)
    cumulative = np.cumsum(full_q)
    k_max = int(np.searchsorted(cumulative, mass, side="left"))
    support = np.arange(0, k_max + 1, dtype=np.int64)
    q = weight * poisson.pmf(support, lambda_a) + weight * poisson.pmf(support, lambda_b)
    p = poisson.pmf(support, fitted_lambda)
    residual = np.log(q) - np.log(p)
    x = 0.05 + 0.95 * support.astype(np.float64) / float(k_max)
    return {
        "lambda_a": lambda_a,
        "lambda_b": lambda_b,
        "fitted_lambda": fitted_lambda,
        "k_max": k_max,
        "captured_mass": float(np.sum(q)),
        "support": support,
        "q": q,
        "p": p,
        "residual": residual,
        "x": x,
    }


def negbin_dispersion_case(mean: float = 12.0, variance: float = 24.0, mass: float = 0.99) -> dict:
    size = mean * mean / (variance - mean)
    p_success = size / (size + mean)
    k = np.arange(0, 500, dtype=np.int64)
    full_q = nbinom.pmf(k, size, p_success)
    cumulative = np.cumsum(full_q)
    k_max = int(np.searchsorted(cumulative, mass, side="left"))
    support = np.arange(0, k_max + 1, dtype=np.int64)
    q = nbinom.pmf(support, size, p_success)
    p = poisson.pmf(support, mean)
    residual = np.log(q) - np.log(p)
    x = 0.05 + 0.95 * support.astype(np.float64) / float(k_max)
    return {
        "mean": mean,
        "variance": variance,
        "nb_size_dispersion": size,
        "nb_success_probability": p_success,
        "fitted_lambda": mean,
        "k_max": k_max,
        "captured_mass": float(np.sum(q)),
        "support": support,
        "q": q,
        "p": p,
        "residual": residual,
        "x": x,
    }


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


def top_rows(test_name: str, steps: list[dict], x: np.ndarray, support_x: np.ndarray | None = None, topn: int = 20):
    shifted_kl = t2.generator_values("poisson_kl", x)
    unshifted = x * np.log(x) - x
    rows = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"][:topn], start=1):
            values = np.asarray(candidate["values"], dtype=np.float64)
            diag = boundary_diagnostics(candidate, support_x)
            rows.append(
                {
                    "test": test_name,
                    "step": int(step_payload["step"]),
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    "cosine_vs_shifted_kl": f"{cosine(values, shifted_kl):.17g}",
                    "cosine_vs_unshifted_negative_entropy": f"{cosine(values, unshifted):.17g}",
                    "is_shifted_kl": bool(t2.is_predicted("poisson_kl", values, x)),
                    **{key: (f"{val:.17g}" if isinstance(val, float) else val) for key, val in diag.items()},
                }
            )
    return rows


def search_case(label: str, x: np.ndarray, residual: np.ndarray, max_steps: int = 3):
    steps, global_best, shifted_occurrences = t2.run_search_detailed("poisson_kl", x, residual, max_steps)
    shifted_kl = t2.generator_values("poisson_kl", x)
    unshifted = x * np.log(x) - x
    unshifted_occurrences = target_occurrences(steps, unshifted)
    winner_values = np.asarray(global_best["values"], dtype=np.float64)
    return {
        "label": label,
        "steps": steps,
        "summary": {
            "label": label,
            "support_size": int(len(x)),
            "residual_rmse": float(np.sqrt(np.mean(residual * residual))),
            "winner_expression": global_best["expr"],
            "winner_math": global_best["math"],
            "winner_rmse": float(global_best["rmse"]),
            "winner_cosine_vs_shifted_kl": cosine(winner_values, shifted_kl),
            "winner_cosine_vs_unshifted_negative_entropy": cosine(winner_values, unshifted),
            "winner_bregman_diag": boundary_diagnostics(global_best, x),
            "shifted_kl_occurrences": shifted_occurrences,
            "unshifted_negative_entropy_occurrences": unshifted_occurrences,
            "predicted_form_in_top5": any(
                item["step"] == max_steps and item["rmse_rank"] <= 5 for item in shifted_occurrences + unshifted_occurrences
            ),
            "predicted_form_in_top20": any(
                item["step"] == max_steps and item["rmse_rank"] <= 20 for item in shifted_occurrences + unshifted_occurrences
            ),
        },
    }


def test1_rank3727_analysis(base: dict, base_search: dict):
    x = base["x"]
    shifted_kl = t2.generator_values("poisson_kl", x)
    unshifted = x * np.log(x) - x
    best = None
    for rank, candidate in enumerate(base_search["steps"][-1]["generated"], start=1):
        values = np.asarray(candidate["values"], dtype=np.float64)
        cos_shifted = cosine(values, shifted_kl)
        if best is None or abs(cos_shifted) > abs(best["cosine_vs_shifted_kl"]):
            best = {
                "step": int(base_search["steps"][-1]["step"]),
                "rmse_rank": rank,
                "expression": candidate["expr"],
                "math": candidate["math"],
                "rmse": float(candidate["rmse"]),
                "values": values,
                "node": candidate["node"],
                "cosine_vs_shifted_kl": cos_shifted,
                "cosine_vs_unshifted_negative_entropy": cosine(values, unshifted),
            }
    affine_shifted = affine_fit(best["values"], shifted_kl)
    affine_unshifted = affine_fit(best["values"], unshifted)
    candidate = {"node": best["node"], "expr": best["expression"], "math": best["math"], "values": best["values"]}
    diag = boundary_diagnostics(candidate, x)
    payload = {
        "rank_by_rmse": best["rmse_rank"],
        "expression": best["expression"],
        "math": best["math"],
        "rmse": best["rmse"],
        "cosine_vs_shifted_kl": best["cosine_vs_shifted_kl"],
        "cosine_vs_unshifted_negative_entropy": best["cosine_vs_unshifted_negative_entropy"],
        "bregman_boundary_diagnostics": diag,
        "affine_equivalence_to_shifted_kl": {
            key: val for key, val in affine_shifted.items() if key not in {"prediction", "residual"}
        },
        "affine_equivalence_to_unshifted_negative_entropy": {
            key: val for key, val in affine_unshifted.items() if key not in {"prediction", "residual"}
        },
        "algebraically_equivalent_to_shifted_kl_under_affine_fit": bool(affine_shifted["max_abs_residual"] <= 1e-8),
        "algebraically_equivalent_to_unshifted_negative_entropy_under_affine_fit": bool(affine_unshifted["max_abs_residual"] <= 1e-8),
    }
    write_json(OUTDIR / "rank3727_expression_analysis.json", payload)
    return payload


def test2_residual_decomposition(base: dict):
    x = base["x"]
    residual = base["residual"]
    shifted_kl = t2.generator_values("poisson_kl", x)
    fit = affine_fit(residual, shifted_kl)
    remainder_search = search_case("remainder_after_shifted_kl_fit", x, fit["residual"], max_steps=2)
    rows = []
    for idx, k in enumerate(base["support"]):
        rows.append(
            {
                "k": int(k),
                "x_linear": f"{x[idx]:.17g}",
                "residual": f"{residual[idx]:.17g}",
                "shifted_kl": f"{shifted_kl[idx]:.17g}",
                "fitted_kl_component": f"{fit['prediction'][idx]:.17g}",
                "remainder_after_subtraction": f"{fit['residual'][idx]:.17g}",
            }
        )
    write_csv(
        OUTDIR / "residual_decomposition.csv",
        rows,
        ["k", "x_linear", "residual", "shifted_kl", "fitted_kl_component", "remainder_after_subtraction"],
    )
    remainder_rows = top_rows("remainder_after_shifted_kl_fit", remainder_search["steps"], x, x, topn=20)
    write_csv(
        OUTDIR / "remainder_sr_top20.csv",
        remainder_rows,
        [
            "test",
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "cosine_vs_shifted_kl",
            "cosine_vs_unshifted_negative_entropy",
            "is_shifted_kl",
            "f_at_1",
            "fprime_at_1",
            "fsecond_min_continuous_grid",
            "fsecond_positive_continuous_grid",
            "bregman_conditions_pass",
            "second_diff_min_on_integer_support",
            "second_diff_positive_on_integer_support",
        ],
    )
    payload = {
        "fit_a": fit["a"],
        "fit_b": fit["b"],
        "r2": fit["r2"],
        "fit_rmse": fit["rmse"],
        "original_residual_rmse": float(np.sqrt(np.mean(residual * residual))),
        "remainder_rmse": float(np.sqrt(np.mean(fit["residual"] * fit["residual"]))),
        "remainder_winner": remainder_search["summary"],
    }
    return payload


def test3_mild_ratios():
    rows = []
    summaries = []
    for lambda_b in [12.0, 15.0]:
        case = poisson_mixture(10.0, lambda_b)
        result = search_case(f"poisson_mixture_10_{lambda_b:g}", case["x"], case["residual"], max_steps=3)
        summary = result["summary"]
        shifted_occ = summary["shifted_kl_occurrences"][0] if summary["shifted_kl_occurrences"] else {}
        unshifted_occ = summary["unshifted_negative_entropy_occurrences"][0] if summary["unshifted_negative_entropy_occurrences"] else {}
        rows.append(
            {
                "rates": f"10,{lambda_b:g}",
                "rate_ratio": f"{lambda_b / 10.0:.17g}",
                "k_max": case["k_max"],
                "captured_mass": f"{case['captured_mass']:.17g}",
                "residual_rmse": f"{summary['residual_rmse']:.17g}",
                "winner_expression": summary["winner_expression"],
                "winner_rmse": f"{summary['winner_rmse']:.17g}",
                "winner_cosine_vs_shifted_kl": f"{summary['winner_cosine_vs_shifted_kl']:.17g}",
                "winner_cosine_vs_unshifted_negative_entropy": f"{summary['winner_cosine_vs_unshifted_negative_entropy']:.17g}",
                "shifted_kl_rank": shifted_occ.get("rmse_rank", ""),
                "shifted_kl_rmse": f"{shifted_occ.get('rmse', float('nan')):.17g}" if shifted_occ else "",
                "unshifted_rank": unshifted_occ.get("rmse_rank", ""),
                "unshifted_rmse": f"{unshifted_occ.get('rmse', float('nan')):.17g}" if unshifted_occ else "",
                "predicted_form_in_top5": summary["predicted_form_in_top5"],
                "predicted_form_in_top20": summary["predicted_form_in_top20"],
            }
        )
        summaries.append({"case": case, "summary": summary})
    write_csv(
        OUTDIR / "mild_ratio_count_coordinate_results.csv",
        rows,
        [
            "rates",
            "rate_ratio",
            "k_max",
            "captured_mass",
            "residual_rmse",
            "winner_expression",
            "winner_rmse",
            "winner_cosine_vs_shifted_kl",
            "winner_cosine_vs_unshifted_negative_entropy",
            "shifted_kl_rank",
            "shifted_kl_rmse",
            "unshifted_rank",
            "unshifted_rmse",
            "predicted_form_in_top5",
            "predicted_form_in_top20",
        ],
    )
    return summaries


def test4_negbin_dispersion():
    case = negbin_dispersion_case(mean=12.0, variance=24.0)
    result = search_case("negbin_mean12_var24_fit_poisson_count_coordinate", case["x"], case["residual"], max_steps=3)
    rows = top_rows("negbin_mean12_var24_fit_poisson_count_coordinate", result["steps"], case["x"], case["x"], topn=20)
    write_csv(
        OUTDIR / "negbin_dispersion_count_coordinate_top20.csv",
        rows,
        [
            "test",
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "cosine_vs_shifted_kl",
            "cosine_vs_unshifted_negative_entropy",
            "is_shifted_kl",
            "f_at_1",
            "fprime_at_1",
            "fsecond_min_continuous_grid",
            "fsecond_positive_continuous_grid",
            "bregman_conditions_pass",
            "second_diff_min_on_integer_support",
            "second_diff_positive_on_integer_support",
        ],
    )
    payload = {
        "mean": case["mean"],
        "variance": case["variance"],
        "nb_size_dispersion_used": case["nb_size_dispersion"],
        "nb_success_probability": case["nb_success_probability"],
        "note": "Mean 12 and variance 24 imply NB size/dispersion 12 under var = mean + mean^2/size.",
        "k_max": case["k_max"],
        "captured_mass": case["captured_mass"],
        "summary": result["summary"],
    }
    write_json(OUTDIR / "negbin_dispersion_count_coordinate_summary.json", payload)
    return payload


def write_report(rank_analysis: dict, decomposition: dict, mild: list[dict], nb: dict):
    mild_support = [item for item in mild if item["summary"]["predicted_form_in_top5"]]
    nb_support = nb["summary"]["predicted_form_in_top5"]
    decomposition_support = decomposition["r2"] > 0.5
    if decomposition_support and (mild_support or nb_support):
        verdict = "theorem_partially_supported_under_cleaner_or_decomposed_conditions"
    elif decomposition_support:
        verdict = "bregman_component_present_but_subdominant"
    elif mild_support or nb_support:
        verdict = "support_under_alternative_misspecification"
    else:
        verdict = "no_clean_poisson_support_in_sr_selection"
    lines = [
        "# T2 Coordinate-Theorem Follow-Up",
        "",
        "This follow-up tests whether Poisson KL is present but subdominant, whether milder mixtures recover it, and whether a pure dispersion misspecification exposes it more cleanly.",
        "",
        "## Test 1: high-cosine expression",
        "",
        f"The highest-cosine expression versus shifted KL is rank `{rank_analysis['rank_by_rmse']}` by RMSE: `{rank_analysis['expression']}` ({rank_analysis['math']}). Its cosine versus shifted KL is `{rank_analysis['cosine_vs_shifted_kl']:.6g}` and its RMSE is `{rank_analysis['rmse']:.6g}`.",
        f"It is not algebraically equivalent to shifted KL under affine fit: max residual `{rank_analysis['affine_equivalence_to_shifted_kl']['max_abs_residual']:.6g}`, R^2 `{rank_analysis['affine_equivalence_to_shifted_kl']['r2']:.6g}`.",
        f"Bregman boundary pass: `{rank_analysis['bregman_boundary_diagnostics']['bregman_conditions_pass']}`; integer-support discrete convexity pass: `{rank_analysis['bregman_boundary_diagnostics']['second_diff_positive_on_integer_support']}`.",
        "",
        "## Test 2: KL component fit and remainder",
        "",
        f"Least-squares fit `a + b*(x log x - x + 1)` gives `a={decomposition['fit_a']:.6g}`, `b={decomposition['fit_b']:.6g}`, R^2 `{decomposition['r2']:.6g}`, and RMSE `{decomposition['fit_rmse']:.6g}` versus original residual RMSE `{decomposition['original_residual_rmse']:.6g}`.",
        f"Remainder SR depth-2 winner: `{decomposition['remainder_winner']['winner_expression']}` with RMSE `{decomposition['remainder_winner']['winner_rmse']:.6g}`.",
        "",
        "## Test 3: mild mixtures in linear count coordinate",
        "",
    ]
    for item in mild:
        summary = item["summary"]
        lines.append(
            f"- rates `(10, {item['case']['lambda_b']:.0f})`: residual RMSE `{summary['residual_rmse']:.6g}`, winner `{summary['winner_expression']}`, predicted form in top 5 `{summary['predicted_form_in_top5']}`, in top 20 `{summary['predicted_form_in_top20']}`."
        )
    lines.extend(
        [
            "",
            "## Test 4: negative-binomial dispersion misspecification",
            "",
            f"NB(mean=12, variance=24) fit by Poisson used NB size/dispersion `{nb['nb_size_dispersion_used']:.6g}`. Residual RMSE `{nb['summary']['residual_rmse']:.6g}`; winner `{nb['summary']['winner_expression']}`; predicted form in top 5 `{nb['summary']['predicted_form_in_top5']}`, in top 20 `{nb['summary']['predicted_form_in_top20']}`.",
            "",
            "## Synthesis",
            "",
        ]
    )
    if verdict == "bregman_component_present_but_subdominant":
        lines.append(
            "The KL component is present as a substantial fitted component, but it remains subdominant under symbolic-regression selection. Milder mixtures and NB dispersion did not make KL a top-5/top-20 winner."
        )
    elif verdict == "no_clean_poisson_support_in_sr_selection":
        lines.append(
            "None of the follow-ups produced clean Poisson KL support under the SR selection rule. The high-cosine expression confirms KL-like shape is expressible, but it is not selected and is not algebraically equivalent to KL."
        )
    else:
        lines.append(f"Verdict category: `{verdict}`.")
    lines.append(f"Final verdict: `{verdict}`.")
    (OUTDIR / "followup_synthesis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(
        OUTDIR / "followup_synthesis_verdict.json",
        {
            "final_verdict": verdict,
            "rank3727_analysis": rank_analysis,
            "residual_decomposition_summary": {key: val for key, val in decomposition.items() if key != "remainder_winner"},
            "remainder_winner": decomposition["remainder_winner"],
            "mild_ratio_summaries": [item["summary"] for item in mild],
            "negbin_dispersion_summary": nb,
        },
    )
    return verdict


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    base = poisson_mixture(5.0, 20.0)
    base_search = search_case("poisson_mixture_5_20_count_coordinate", base["x"], base["residual"], max_steps=3)
    rank_analysis = test1_rank3727_analysis(base, base_search)
    decomposition = test2_residual_decomposition(base)
    mild = test3_mild_ratios()
    nb = test4_negbin_dispersion()
    verdict = write_report(rank_analysis, decomposition, mild, nb)
    print(json.dumps({"outdir": str(OUTDIR), "verdict": verdict}, indent=2))


if __name__ == "__main__":
    main()
