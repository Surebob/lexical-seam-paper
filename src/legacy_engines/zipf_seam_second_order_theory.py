from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
LOCAL_THEORY_SUMMARY_PATH = ROOT / "results" / "zipf_seam_sign_theory" / "summary.json"
PROJECTION_THEORY_SUMMARY_PATH = ROOT / "results" / "zipf_seam_projection_theory" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_seam_second_order_theory"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_second_order_common")
reranked = load_module(RERANKED_PATH, "zipf_seam_second_order_reranked")


def load_rows():
    return {row["slug"]: row for row in json.loads(RERANKED_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]}


def centered_poly_signs(values: np.ndarray, head_k: int = 50) -> tuple[list[float], list[int]]:
    head_k = min(head_k, len(values))
    h = np.log(np.arange(1, head_k + 1, dtype=np.float64)) / math.log(head_k)
    y = values[:head_k] - values[0]
    X = np.column_stack([h, h * h, h * h * h])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    signs = np.sign(coeffs).astype(int).tolist()
    return coeffs.tolist(), signs


def head_branch(ranks: np.ndarray, params: dict) -> np.ndarray:
    return params["a1"] - params["b1"] * np.log(ranks + max(params["c1"], 0.0))


def exact_smooth_residual(ranks: np.ndarray, smooth_prediction: np.ndarray) -> np.ndarray:
    fit = common.fit_zipf_mandelbrot(ranks, smooth_prediction)
    return smooth_prediction - fit["prediction"]


def seam_curve(ranks: np.ndarray, params: dict) -> np.ndarray:
    sigma = reranked.sigma_curve(ranks, params["k"], params["w"])
    tail_rank = reranked.smooth_tail_local_rank(ranks, params["k"], params["w"])
    head = head_branch(ranks, params)
    tail = params["a2"] - params["b2"] * np.log(tail_rank + max(params["c2"], 0.0))
    return (1.0 - sigma) * (tail - head)


def tangent_matrix(ranks: np.ndarray, params: dict) -> np.ndarray:
    c1 = max(params["c1"], 0.0)
    b1 = params["b1"]
    denom = ranks + c1
    return np.column_stack(
        [
            np.ones_like(ranks),
            -np.log(denom),
            -b1 / denom,
        ]
    )


def quadratic_term(ranks: np.ndarray, params: dict, delta: np.ndarray) -> np.ndarray:
    c1 = max(params["c1"], 0.0)
    b1 = params["b1"]
    db = float(delta[1])
    dc = float(delta[2])
    denom = ranks + c1
    return (-2.0 * db * dc / denom) + (b1 * dc * dc / (denom * denom))


def delta_hessian_columns(ranks: np.ndarray, params: dict, delta: np.ndarray) -> np.ndarray:
    c1 = max(params["c1"], 0.0)
    b1 = params["b1"]
    db = float(delta[1])
    dc = float(delta[2])
    denom = ranks + c1
    return np.column_stack(
        [
            np.zeros_like(ranks),
            -dc / denom,
            -db / denom + (b1 * dc / (denom * denom)),
        ]
    )


def first_order_projection(ranks: np.ndarray, params: dict, seam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    J = tangent_matrix(ranks, params)
    delta1, _, _, _ = np.linalg.lstsq(J, seam, rcond=None)
    residual = seam - J @ delta1
    return delta1, residual


def second_order_projection(ranks: np.ndarray, params: dict, seam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    J = tangent_matrix(ranks, params)
    jt_j = J.T @ J
    delta1, _, _, _ = np.linalg.lstsq(J, seam, rcond=None)
    r1 = seam - J @ delta1
    b11 = quadratic_term(ranks, params, delta1)
    hcols = delta_hessian_columns(ranks, params, delta1)
    rhs = hcols.T @ r1 - 0.5 * (J.T @ b11)
    delta2 = np.linalg.solve(jt_j, rhs)
    delta = delta1 + delta2
    residual = seam - J @ delta - 0.5 * b11
    return delta, residual


def quadratic_surrogate_projection(ranks: np.ndarray, params: dict, seam: np.ndarray, x0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    J = tangent_matrix(ranks, params)

    def residuals(delta: np.ndarray):
        return seam - J @ delta - 0.5 * quadratic_term(ranks, params, delta)

    result = least_squares(residuals, x0=x0, method="lm", max_nfev=2000)
    delta = result.x
    return delta, residuals(delta)


def count_matches(rows: list[dict], key: str) -> list[int]:
    return [int(sum(row[key][i] for row in rows)) for i in range(3)]


def full_match_count(rows: list[dict], key: str) -> int:
    return int(sum(all(row[key]) for row in rows))


def analyze_corpus(spec: dict, reranked_row: dict, local_row: dict) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    summary = common.load_enriched_summary(spec)
    params = reranked_row["params"]
    param_vec = np.array(
        [params["a1"], params["b1"], params["c1"], params["a2"], params["b2"], params["c2"], params["k"], params["w"]],
        dtype=np.float64,
    )

    empirical_residual = dataset["log_freq"] - common.zm_prediction(summary, dataset["ranks"])
    _, empirical_signs = centered_poly_signs(empirical_residual, 50)

    smooth_prediction = reranked.reranked_prediction(dataset["ranks"], param_vec)
    smooth_residual = exact_smooth_residual(dataset["ranks"], smooth_prediction)
    _, smooth_signs = centered_poly_signs(smooth_residual, 50)

    seam = seam_curve(dataset["ranks"], params)
    _, proj_residual = first_order_projection(dataset["ranks"], params, seam)
    _, proj_signs = centered_poly_signs(proj_residual, 50)

    delta2, second_residual = second_order_projection(dataset["ranks"], params, seam)
    _, second_signs = centered_poly_signs(second_residual, 50)

    _, quad_residual = quadratic_surrogate_projection(dataset["ranks"], params, seam, delta2)
    _, quad_signs = centered_poly_signs(quad_residual, 50)

    family = "high_family" if common.get_step2_expr(summary) == common.STEP2_OURS_EXPR else "low_family"
    local_signs = local_row["local_signs"]
    return {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "family": family,
        "c": float(summary["zm_baseline"]["c"]),
        "local_signs": local_signs,
        "projected_signs": proj_signs,
        "second_order_signs": second_signs,
        "quadratic_surrogate_signs": quad_signs,
        "smooth_signs": smooth_signs,
        "empirical_signs": empirical_signs,
        "match_proj_vs_emp": [proj_signs[i] == empirical_signs[i] for i in range(3)],
        "match_second_vs_emp": [second_signs[i] == empirical_signs[i] for i in range(3)],
        "match_quad_vs_emp": [quad_signs[i] == empirical_signs[i] for i in range(3)],
        "match_smooth_vs_emp": [smooth_signs[i] == empirical_signs[i] for i in range(3)],
        "match_second_vs_smooth": [second_signs[i] == smooth_signs[i] for i in range(3)],
        "match_quad_vs_smooth": [quad_signs[i] == smooth_signs[i] for i in range(3)],
    }


def summarize(rows: list[dict], local_summary: dict, projection_summary: dict) -> dict:
    return {
        "local_baseline_head50_term_matches": local_summary["summary"]["local_match_counts_head50"],
        "local_baseline_head50_full_matches": local_summary["summary"]["local_full_match_head50"],
        "first_order_head50_term_matches": projection_summary["summary"]["projected_head50_term_matches"],
        "first_order_head50_full_matches": projection_summary["summary"]["projected_head50_full_matches"],
        "exact_smooth_head50_term_matches": local_summary["summary"]["smooth_match_counts_head50"],
        "exact_smooth_head50_full_matches": local_summary["summary"]["smooth_full_match_head50"],
        "second_order_head50_term_matches": count_matches(rows, "match_second_vs_emp"),
        "second_order_head50_full_matches": full_match_count(rows, "match_second_vs_emp"),
        "quadratic_surrogate_head50_term_matches": count_matches(rows, "match_quad_vs_emp"),
        "quadratic_surrogate_head50_full_matches": full_match_count(rows, "match_quad_vs_emp"),
        "second_order_vs_smooth_term_matches": count_matches(rows, "match_second_vs_smooth"),
        "second_order_vs_smooth_full_matches": full_match_count(rows, "match_second_vs_smooth"),
        "quadratic_surrogate_vs_smooth_term_matches": count_matches(rows, "match_quad_vs_smooth"),
        "quadratic_surrogate_vs_smooth_full_matches": full_match_count(rows, "match_quad_vs_smooth"),
    }


def build_report(summary: dict, rows: list[dict]) -> str:
    lines = [
        "# Seam Second-Order Theory",
        "",
        "## Quadratic ZM Refit Expansion",
        "",
        "Start from the head-branch ZM `Y0(r) = a1 - b1 log(r + c1)` and let the smooth seam perturbation be `s(r)`.",
        "",
        "For a single-ZM refit with parameter shift `delta = (da, db, dc)`, the rankwise expansion is:",
        "",
        "`Y(theta0 + delta) = Y0 + J delta + 1/2 H[delta, delta] + O(||delta||^3)`",
        "",
        "with ZM tangent basis",
        "",
        "- `ga(r) = 1`",
        "- `gb(r) = -log(r + c1)`",
        "- `gc(r) = -b1 / (r + c1)`",
        "",
        "and only two nonzero second derivatives:",
        "",
        "- `d^2Y / db dc = -1 / (r + c1)`",
        "- `d^2Y / dc^2 = b1 / (r + c1)^2`",
        "",
        "So the quadratic curvature term along `delta` is",
        "",
        "`H[delta,delta] = -2 db dc / (r + c1) + b1 dc^2 / (r + c1)^2`.",
        "",
        "The first-order projected theorem subtracts only the tangent-space component:",
        "",
        "`delta1 = (J^T J)^(-1) J^T s`",
        "",
        "`r_proj = s - J delta1`",
        "",
        "The second-order correction adds the quadratic normal-equation term:",
        "",
        "`delta2 = (J^T J)^(-1) [ K(delta1)^T r1 - 1/2 J^T H[delta1,delta1] ]`",
        "",
        "where `r1 = s - J delta1` and the columns of `K(delta1)` are `H[delta1, e_j]`.",
        "",
        "This gives the second-order residual predictor",
        "",
        "`r_second = s - J(delta1 + delta2) - 1/2 H[delta1,delta1]`.",
        "",
        "I also solved the full quadratic surrogate exactly by least squares as a ceiling check:",
        "",
        "`min_delta || s - J delta - 1/2 H[delta,delta] ||_2`.",
        "",
        "## Head-50 Sign Match Counts Against Empirical Single-ZM Residual",
        "",
        f"- raw local seam theorem: `{summary['local_baseline_head50_term_matches']}` term matches, `{summary['local_baseline_head50_full_matches']}` full matches",
        f"- first-order tangent projection: `{summary['first_order_head50_term_matches']}` term matches, `{summary['first_order_head50_full_matches']}` full matches",
        f"- second-order projected theorem: `{summary['second_order_head50_term_matches']}` term matches, `{summary['second_order_head50_full_matches']}` full matches",
        f"- exact quadratic surrogate: `{summary['quadratic_surrogate_head50_term_matches']}` term matches, `{summary['quadratic_surrogate_head50_full_matches']}` full matches",
        f"- exact smooth-model refit residual: `{summary['exact_smooth_head50_term_matches']}` term matches, `{summary['exact_smooth_head50_full_matches']}` full matches",
        "",
        "## Match Counts Against Exact Smooth-Model Refit Residual",
        "",
        f"- second-order projected theorem vs exact smooth residual: `{summary['second_order_vs_smooth_term_matches']}` term matches, `{summary['second_order_vs_smooth_full_matches']}` full matches",
        f"- quadratic surrogate vs exact smooth residual: `{summary['quadratic_surrogate_vs_smooth_term_matches']}` term matches, `{summary['quadratic_surrogate_vs_smooth_full_matches']}` full matches",
        "",
        "## Per-Corpus Head-50 Sign Vectors",
        "",
        "| corpus | family | c | first order | second order | quadratic surrogate | exact smooth residual | empirical residual |",
        "| --- | --- | ---: | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['corpus']} | {row['family']} | {row['c']:.3f} | `{row['projected_signs']}` | `{row['second_order_signs']}` | `{row['quadratic_surrogate_signs']}` | `{row['smooth_signs']}` | `{row['empirical_signs']}` |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_rows = load_rows()
    local_summary = json.loads(LOCAL_THEORY_SUMMARY_PATH.read_text(encoding="utf-8"))
    projection_summary = json.loads(PROJECTION_THEORY_SUMMARY_PATH.read_text(encoding="utf-8"))
    local_rows = {row["slug"]: row for row in local_summary["rows"]}
    rows = [analyze_corpus(spec, reranked_rows[spec["slug"]], local_rows[spec["slug"]]) for spec in common.SEARCHED_CORPORA]
    summary = summarize(rows, local_summary, projection_summary)
    payload = {"summary": summary, "rows": rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary, rows), encoding="utf-8")


if __name__ == "__main__":
    main()
