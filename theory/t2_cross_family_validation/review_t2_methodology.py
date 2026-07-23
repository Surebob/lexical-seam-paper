from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
T2_DIR = ROOT / "phase2_addon" / "t2_cross_family_validation"
RUN_T2_PATH = T2_DIR / "run_t2.py"


def load_run_t2():
    spec = importlib.util.spec_from_file_location("run_t2_review_module", RUN_T2_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


t2 = load_run_t2()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_search_detailed(x: np.ndarray, target: np.ndarray, max_steps: int):
    current_vocabulary = t2.search.initial_vocabulary(x, target)
    steps = []
    global_best = None

    for step in range(1, max_steps + 1):
        generated = t2.search.generate_candidates(current_vocabulary, target, step)
        generated = t2.search.dedupe_candidates(generated)
        generated = t2.search.filter_candidates(generated, t2.CONSTANT_VARIANCE_THRESHOLD)
        generated = sorted(generated, key=lambda item: (item["rmse"], item["expr"]))
        if not generated:
            break

        step_best = generated[0]
        if global_best is None or (step_best["rmse"], step_best["expr"]) < (global_best["rmse"], global_best["expr"]):
            global_best = step_best

        if step <= t2.KEEP_ALL_UNTIL_STEP:
            if step < t2.KEEP_ALL_UNTIL_STEP:
                selected_for_next = current_vocabulary + generated
                selected_generated = generated
            else:
                selected_generated = t2.search.select_diverse_beam(generated, t2.BEAM_WIDTH, t2.DIVERSITY_WEIGHT)
                selected_for_next = selected_generated
        else:
            selected_generated = t2.search.select_diverse_beam(generated, t2.BEAM_WIDTH, t2.DIVERSITY_WEIGHT)
            selected_for_next = selected_generated

        steps.append(
            {
                "step": step,
                "generated": generated,
                "selected_generated": selected_generated,
                "selected_exprs": {item["expr"] for item in selected_generated},
            }
        )
        current_vocabulary = selected_for_next

    return steps, global_best


def predicted_occurrences(kind: str, steps: list[dict], x: np.ndarray):
    rows = []
    for step_payload in steps:
        selected_exprs = step_payload["selected_exprs"]
        for rank, candidate in enumerate(step_payload["generated"], start=1):
            if t2.is_predicted_candidate(kind, t2.candidate_values(candidate), x):
                rows.append(
                    {
                        "step": step_payload["step"],
                        "rmse_rank": rank,
                        "expr": candidate["expr"],
                        "math": candidate["math"],
                        "rmse": float(candidate["rmse"]),
                        "selected_by_diversity_for_next_step": candidate["expr"] in selected_exprs,
                    }
                )
                break
    return rows


def top_rows(case: dict, steps: list[dict], x: np.ndarray, top_n: int = 5):
    rows = []
    for step_payload in steps:
        for rank, candidate in enumerate(step_payload["generated"][:top_n], start=1):
            rows.append(
                {
                    "step": step_payload["step"],
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": f"{float(candidate['rmse']):.17g}",
                    "classification": t2.classify_expression(candidate["expr"], t2.candidate_values(candidate), x),
                    "selected_by_diversity_for_next_step": candidate["expr"] in step_payload["selected_exprs"],
                }
            )
    return rows


def fit_candidate_decomposition(x: np.ndarray, candidate_values: np.ndarray, predicted_values: np.ndarray, degree: int = 5):
    powers = [np.ones_like(x)]
    for d in range(1, degree + 1):
        powers.append(x**d)
    poly = np.column_stack(powers)
    pred_plus_poly = np.column_stack([predicted_values, poly])

    def metrics(design: np.ndarray):
        coef, *_ = np.linalg.lstsq(design, candidate_values, rcond=None)
        fitted = design @ coef
        ss_res = float(np.sum((candidate_values - fitted) ** 2))
        ss_tot = float(np.sum((candidate_values - np.mean(candidate_values)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        err = float(np.sqrt(np.mean((candidate_values - fitted) ** 2)))
        return coef, fitted, r2, err

    coef_poly, _, r2_poly, rmse_poly = metrics(poly)
    coef_pred_poly, _, r2_pred_poly, rmse_pred_poly = metrics(pred_plus_poly)
    return {
        "poly_degree": degree,
        "poly_only_r2": r2_poly,
        "poly_only_rmse": rmse_poly,
        "predicted_plus_poly_r2": r2_pred_poly,
        "predicted_plus_poly_rmse": rmse_pred_poly,
        "predicted_generator_coefficient": float(coef_pred_poly[0]),
        "poly_coefficients": [float(v) for v in coef_poly],
        "predicted_plus_poly_coefficients": [float(v) for v in coef_pred_poly],
    }


def build_case_runtime(case: dict, max_steps: int):
    x = t2.normalize_rank_axis(len(case["frequency"]))
    residual = case["log_frequency"] - case["fitted_log_frequency"]
    steps, global_best = run_search_detailed(x, residual, max_steps=max_steps)
    occurrences = predicted_occurrences(case["kind"], steps, x)
    return x, residual, steps, global_best, occurrences


def main():
    outdir = T2_DIR
    depth2_dir = outdir / "t2_depth2_review"
    decomp_dir = outdir / "t2_depth4_decomposition"
    depth2_dir.mkdir(parents=True, exist_ok=True)
    decomp_dir.mkdir(parents=True, exist_ok=True)

    cases = [t2.build_gaussian_case(), t2.build_poisson_case(), t2.build_gamma_case()]

    depth2_summary = []
    decomp_rows = []
    protocol_lines = [
        "# T2 Protocol Review",
        "",
        "This review separates three questions that the first depth-4 aggregate can conflate: data/fitting space, canonical depth-2 behavior, and whether depth-4 winners are interpretable or polynomial-like overfits.",
        "",
        "## Live Search Semantics",
        "",
        "- Generated candidates are sorted by RMSE for reporting.",
        "- The diversity beam is the selected vocabulary carried to the next step; it is not the same as the top-50-by-RMSE list.",
        "- The manuscript's step-2 winner language corresponds to the lowest-RMSE step-2 candidate, while diversity selection only controls what expressions remain available for deeper search.",
        "",
        "## Synthetic Data/Fitting Space",
        "",
    ]

    for case in cases:
        x2, residual2, steps2, global2, occ2 = build_case_runtime(case, max_steps=2)
        case_depth2_dir = depth2_dir / case["case_id"]
        case_depth2_dir.mkdir(parents=True, exist_ok=True)
        write_csv(
            case_depth2_dir / "per_step_beam_depth2.csv",
            top_rows(case, steps2, x2, top_n=5),
            ["step", "rank", "expression", "math", "rmse", "classification", "selected_by_diversity_for_next_step"],
        )
        predicted_wins_depth2 = bool(
            occ2
            and global2 is not None
            and abs(float(occ2[0]["rmse"]) - float(global2["rmse"])) <= 1e-10
        )
        depth2_item = {
            "case_id": case["case_id"],
            "title": case["title"],
            "predicted_generator": t2.PREDICTED_EXPRESSIONS[case["kind"]]["label"],
            "predicted_generated_depth2": bool(occ2),
            "predicted_occurrences_depth2": occ2,
            "global_best_depth2": None
            if global2 is None
            else {
                "step": int(global2["step"]),
                "expr": global2["expr"],
                "math": global2["math"],
                "rmse": float(global2["rmse"]),
                "classification": t2.classify_expression(global2["expr"], t2.candidate_values(global2), x2),
            },
            "predicted_wins_depth2": predicted_wins_depth2,
        }
        depth2_summary.append(depth2_item)
        (case_depth2_dir / "report_depth2.md").write_text(
            "\n".join(
                [
                    f"# Depth-2 Review: {case['title']}",
                    "",
                    f"- predicted generator: `{depth2_item['predicted_generator']}`",
                    f"- predicted generated by depth 2: `{depth2_item['predicted_generated_depth2']}`",
                    f"- predicted wins by depth 2: `{depth2_item['predicted_wins_depth2']}`",
                    f"- predicted occurrences: `{occ2}`",
                    f"- global best depth 2: `{depth2_item['global_best_depth2']}`",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        # Depth-4 decomposition of the operational global best expression.
        x4, residual4, steps4, global4, occ4 = build_case_runtime(case, max_steps=4)
        if global4 is not None:
            pred_vals = t2.predicted_values(case["kind"], x4)
            metrics = fit_candidate_decomposition(x4, t2.candidate_values(global4), pred_vals, degree=5)
            decomp_rows.append(
                {
                    "case_id": case["case_id"],
                    "global_best_expr": global4["expr"],
                    "global_best_rmse_to_residual": f"{float(global4['rmse']):.17g}",
                    "predicted_generator": t2.PREDICTED_EXPRESSIONS[case["kind"]]["label"],
                    "predicted_generated_depth4": bool(occ4),
                    "predicted_first_occurrence": json.dumps(occ4[0]) if occ4 else "",
                    "poly_only_r2": f"{metrics['poly_only_r2']:.17g}",
                    "poly_only_rmse": f"{metrics['poly_only_rmse']:.17g}",
                    "predicted_plus_poly_r2": f"{metrics['predicted_plus_poly_r2']:.17g}",
                    "predicted_plus_poly_rmse": f"{metrics['predicted_plus_poly_rmse']:.17g}",
                    "predicted_generator_coefficient": f"{metrics['predicted_generator_coefficient']:.17g}",
                }
            )

        protocol_lines.extend(
            [
                f"### {case['title']}",
                "",
                f"- generated data space: sorted rank-frequency curve formed from a two-component `{case['title'].split(' -> ')[0].split(' ', 1)[1]}` density/PMF evaluated on a fixed support grid.",
                f"- fitted model space: `{case['fit']['method']}`.",
                "- SR residual: `log_frequency - fitted_log_frequency` on the sorted rank-frequency curve.",
                f"- fitted parameters: `{case['fit']['params']}`.",
                f"- bounds hit: `{case['fit']['at_bounds']}`.",
                f"- mixture curvature diagnostic: `{case['diagnostic']}`.",
                "",
            ]
        )

    depth2_payload = {
        "max_steps": 2,
        "note": "Depth-2 review uses the same live search implementation and reports both RMSE ranking and diversity-beam retention.",
        "cases": depth2_summary,
        "predicted_win_count_depth2": sum(1 for item in depth2_summary if item["predicted_wins_depth2"]),
    }
    (depth2_dir / "t2_depth2_summary.json").write_text(json.dumps(depth2_payload, indent=2) + "\n", encoding="utf-8")

    write_csv(
        decomp_dir / "t2_depth4_decomposition_summary.csv",
        decomp_rows,
        [
            "case_id",
            "global_best_expr",
            "global_best_rmse_to_residual",
            "predicted_generator",
            "predicted_generated_depth4",
            "predicted_first_occurrence",
            "poly_only_r2",
            "poly_only_rmse",
            "predicted_plus_poly_r2",
            "predicted_plus_poly_rmse",
            "predicted_generator_coefficient",
        ],
    )
    (outdir / "t2_protocol_review.md").write_text("\n".join(protocol_lines) + "\n", encoding="utf-8")

    report_lines = [
        "# T2 Methodology Review",
        "",
        f"Depth-2 predicted wins: `{depth2_payload['predicted_win_count_depth2']}/3`.",
        "",
        "Depth-2 did not rescue the generator-identity prediction under this synthetic setup: Gaussian and Gamma predicted generators were generated but ranked outside the RMSE winners, and Poisson KL was not generated by depth 2.",
        "",
        "The data/fitting protocol is not an exact analogue of the ZM rank-frequency construction. Each family was evaluated as a density/PMF on its own support grid, sorted into a rank-frequency curve, fit by the corresponding single-family curve in sorted log-density/log-PMF space, then searched on the sorted residual.",
        "",
        "Depth-4 decomposition diagnostics are saved in `t2_depth4_decomposition/t2_depth4_decomposition_summary.csv`.",
    ]
    (outdir / "t2_methodology_review.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "depth2_predicted_win_count": depth2_payload["predicted_win_count_depth2"],
                "depth2_summary": depth2_summary,
                "decomposition_rows": decomp_rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
