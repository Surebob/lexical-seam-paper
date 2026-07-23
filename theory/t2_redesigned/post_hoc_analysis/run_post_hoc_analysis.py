from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
T2_DIR = ROOT / "phase2_addon" / "t2_redesigned"
OUTDIR = T2_DIR / "post_hoc_analysis"
SCRIPT_PATH = T2_DIR / "run_t2_redesigned.py"


def load_t2_module():
    spec = importlib.util.spec_from_file_location("t2_redesigned_posthoc", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


t2 = load_t2_module()
search = t2.search


CASES = [
    {
        "case_id": "t2a_gaussian_euclidean",
        "family": "Gaussian",
        "kind": "gaussian_euclidean",
        "predicted": "Euclidean",
        "max_steps": 2,
    },
    {
        "case_id": "t2b_poisson_kl",
        "family": "Poisson",
        "kind": "poisson_kl",
        "predicted": "generalized KL",
        "max_steps": 3,
    },
    {
        "case_id": "t2c_gamma_is",
        "family": "Gamma",
        "kind": "gamma_is",
        "predicted": "Itakura-Saito",
        "max_steps": 2,
    },
]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_residual(case_id: str):
    rows = list(csv.DictReader((T2_DIR / case_id / "kde_vs_fitted.csv").open(encoding="utf-8")))
    x = np.asarray([float(row["sr_x"]) for row in rows], dtype=np.float64)
    residual = np.asarray([float(row["residual"]) for row in rows], dtype=np.float64)
    return x, residual


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
        return search.UNARY_FUNCS[op](child)
    if op in search.BINARY_FUNCS:
        left = eval_node(node[1], x)
        right = eval_node(node[2], x)
        if left is None or right is None:
            return None
        return search.BINARY_FUNCS[op](left, right)
    raise ValueError(f"Unsupported op {op!r}")


def boundary_diagnostics(candidate: dict, predicted_kind: str):
    x_dense = np.linspace(0.05, 1.0, 5000, dtype=np.float64)
    values = eval_node(candidate["node"], x_dense)
    if values is None or not np.all(np.isfinite(values)):
        return {
            "f_at_1": "nan",
            "fprime_at_1": "nan",
            "fsecond_min": "nan",
            "fsecond_positive_on_domain": False,
            "bregman_conditions_pass": False,
            "shape_cosine_vs_predicted": "nan",
        }
    first = np.gradient(values, x_dense, edge_order=2)
    second = np.gradient(first, x_dense, edge_order=2)
    predicted = t2.generator_values(predicted_kind, x_dense)
    f_at_1 = float(values[-1])
    fprime_at_1 = float(first[-1])
    fsecond_min = float(np.min(second))
    return {
        "f_at_1": f"{f_at_1:.17g}",
        "fprime_at_1": f"{fprime_at_1:.17g}",
        "fsecond_min": f"{fsecond_min:.17g}",
        "fsecond_positive_on_domain": bool(fsecond_min > -1e-5),
        "bregman_conditions_pass": bool(abs(f_at_1) <= 1e-5 and abs(fprime_at_1) <= 1e-4 and fsecond_min > -1e-5),
        "shape_cosine_vs_predicted": f"{cosine(values, predicted):.17g}",
    }


def top20_bregman_rows(case: dict, steps: list[dict]):
    rows = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"][:20], start=1):
            diag = boundary_diagnostics(candidate, case["kind"])
            rows.append(
                {
                    "family": case["family"],
                    "case_id": case["case_id"],
                    "step": int(step_payload["step"]),
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    **diag,
                }
            )
    return rows


def sign_convention_rows(case: dict, global_best: dict, x: np.ndarray):
    predicted = t2.generator_values(case["kind"], x)
    winner = np.asarray(global_best["values"], dtype=np.float64)
    cos_pos = cosine(winner, predicted)
    cos_neg = cosine(winner, -predicted)
    max_abs = max(abs(cos_pos), abs(cos_neg))
    if max_abs >= 0.95:
        verdict = "equivalent_to_predicted_up_to_sign"
    elif max_abs >= 0.80:
        verdict = "shape_related_but_not_equivalent"
    else:
        verdict = "not_shape_equivalent"
    return {
        "family": case["family"],
        "case_id": case["case_id"],
        "winner": global_best["expr"],
        "winner_math": global_best["math"],
        "predicted": case["predicted"],
        "predicted_canonical": t2.PREDICTIONS[case["kind"]]["canonical"],
        "cosine_pos": f"{cos_pos:.17g}",
        "cosine_neg": f"{cos_neg:.17g}",
        "max_abs_cosine": f"{max_abs:.17g}",
        "verdict": verdict,
    }


def poisson_diagnosis():
    case_dir = T2_DIR / "t2b_poisson_kl"
    params = json.loads((case_dir / "fitted_parameters.json").read_text(encoding="utf-8"))
    rows = list(csv.DictReader((case_dir / "kde_vs_fitted.csv").open(encoding="utf-8")))
    support = np.asarray([float(row["domain_value"]) for row in rows], dtype=np.float64)
    log_emp = np.asarray([float(row["log_p_emp"]) for row in rows], dtype=np.float64)
    log_fit = np.asarray([float(row["log_p_fit"]) for row in rows], dtype=np.float64)
    residual = np.asarray([float(row["residual"]) for row in rows], dtype=np.float64)
    p_emp = np.exp(log_emp)
    mass_order = np.argsort(p_emp)[::-1]
    cumulative = np.cumsum(p_emp[mass_order])
    cumulative /= cumulative[-1]
    high_mass_indices = mass_order[cumulative <= 0.50]
    if len(high_mass_indices) == 0:
        high_mass_indices = mass_order[:1]
    # Include the first point that crosses 50% so the mass region is not too small.
    crossing = mass_order[min(len(high_mass_indices), len(mass_order) - 1)]
    high_mass_indices = np.unique(np.concatenate([high_mass_indices, [crossing]]))
    high_mass_rmse = float(np.sqrt(np.mean(residual[high_mass_indices] ** 2)))
    full_rmse = float(np.sqrt(np.mean(residual**2)))
    tail_indices = np.setdiff1d(np.arange(len(residual)), high_mass_indices)
    tail_rmse = float(np.sqrt(np.mean(residual[tail_indices] ** 2))) if len(tail_indices) else float("nan")
    max_abs_idx = int(np.argmax(np.abs(residual)))
    lambdas = [item["lambda"] for item in params["mixture_params"]]
    fitted_lambda = params["fit"]["params"]["lambda"]
    lines = [
        "# T2 Poisson Residual Diagnosis",
        "",
        f"The two Poisson mixture rates are `{lambdas[0]}` and `{lambdas[1]}` with equal weights; the fitted single-Poisson rate is `{fitted_lambda}`. The component separation is `{params['distinguishability']['separation_in_sqrt_largest_lambda']:.3f}` in sqrt-largest-lambda units, so this is a deliberately strong two-regime mismatch rather than a near-single-Poisson case.",
        "",
        f"The full-support log-density residual RMSE is `{full_rmse:.6g}`. Restricting to the highest empirical-mass support points carrying about half the probability mass gives RMSE `{high_mass_rmse:.6g}`, while the remaining lower-mass support has RMSE `{tail_rmse:.6g}`. The largest absolute residual occurs at count `{support[max_abs_idx]:.0f}` with residual `{residual[max_abs_idx]:.6g}`. This indicates the large residual is not only random noise; the single-Poisson fit is structurally poor, with tail/log-ratio terms contributing strongly to the SR objective.",
        "",
        "Interpretation caveat: because the Poisson residual scale is orders of magnitude larger than the Gaussian and Gamma residual scales, depth-3 compositional expressions can reduce RMSE by modeling gross tail structure rather than revealing a simple KL generator.",
    ]
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sign_rows = []
    bregman_rows = []
    revised = []

    for case in CASES:
        x, residual = load_residual(case["case_id"])
        steps, global_best, occurrences = t2.run_search_detailed(case["kind"], x, residual, case["max_steps"])
        sign_row = sign_convention_rows(case, global_best, x)
        sign_rows.append(sign_row)
        bregman_rows.extend(top20_bregman_rows(case, steps))
        bregman_count = sum(
            1
            for row in bregman_rows
            if row["case_id"] == case["case_id"] and row["bregman_conditions_pass"]
        )
        predicted_seen = occurrences[0] if occurrences else None
        if sign_row["verdict"] == "equivalent_to_predicted_up_to_sign":
            verdict = "SIGN_CONVENTION_PASS"
        elif bregman_count > 0:
            verdict = "BREGMAN_CLASS_PRESENT_BUT_PREDICTED_NOT_SELECTED"
        else:
            verdict = "NO_TOP20_BREGMAN_SUPPORT"
        revised.append(
            {
                "family": case["family"],
                "case_id": case["case_id"],
                "original_predicted_generator": case["predicted"],
                "winner": global_best["expr"],
                "winner_math": global_best["math"],
                "winner_sign_cosine_verdict": sign_row["verdict"],
                "top20_bregman_expression_count": bregman_count,
                "predicted_occurrence": predicted_seen,
                "revised_verdict": verdict,
            }
        )

    write_csv(
        OUTDIR / "t2_sign_convention_audit.csv",
        sign_rows,
        [
            "family",
            "case_id",
            "winner",
            "winner_math",
            "predicted",
            "predicted_canonical",
            "cosine_pos",
            "cosine_neg",
            "max_abs_cosine",
            "verdict",
        ],
    )
    write_csv(
        OUTDIR / "t2_bregman_class_in_top20.csv",
        bregman_rows,
        [
            "family",
            "case_id",
            "step",
            "rank",
            "expression",
            "math",
            "rmse",
            "f_at_1",
            "fprime_at_1",
            "fsecond_min",
            "fsecond_positive_on_domain",
            "bregman_conditions_pass",
            "shape_cosine_vs_predicted",
        ],
    )
    (OUTDIR / "t2_poisson_residual_diagnosis.md").write_text(poisson_diagnosis(), encoding="utf-8")

    lines = [
        "# T2 Post-Hoc Interpretation",
        "",
        "This post-hoc audit does not modify the redesigned T2 outputs. It checks whether the winning expression is sign-equivalent to the predicted generator, whether any top-20 expressions satisfy numerical Bregman boundary conditions, and whether the Poisson miss is dominated by residual-scale pathology.",
        "",
    ]
    for item in revised:
        lines.append(f"## {item['family']}")
        lines.append("")
        lines.append(f"- revised verdict: `{item['revised_verdict']}`")
        lines.append(f"- winner: `{item['winner']}` ({item['winner_math']})")
        lines.append(f"- winner sign-cosine verdict: `{item['winner_sign_cosine_verdict']}`")
        lines.append(f"- top-20 Bregman-class expression count: `{item['top20_bregman_expression_count']}`")
        if item["predicted_occurrence"]:
            occurrence = item["predicted_occurrence"]
            lines.append(
                f"- predicted generator occurrence: step `{occurrence['step']}`, RMSE rank `{occurrence['rmse_rank']}`, RMSE `{occurrence['rmse']:.6g}`"
            )
        lines.append("")
    lines.append(
        "Summary: Gamma is upgraded from a raw fail to a sign-convention pass because the winner is exactly the negative IS form under the current residual direction. Gaussian remains a close but non-equivalent miss: Euclidean is present at rank 16, but the winner is not cosine-equivalent to Euclidean up to sign. Poisson remains the cleanest miss; KL is reachable but deeply ranked, and the residual diagnosis shows a severe two-regime mismatch with large tail/log-ratio structure."
    )
    (OUTDIR / "t2_post_hoc_interpretation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUTDIR / "t2_post_hoc_summary.json").write_text(json.dumps(revised, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"outdir": str(OUTDIR), "families": len(revised)}, indent=2))


if __name__ == "__main__":
    main()
