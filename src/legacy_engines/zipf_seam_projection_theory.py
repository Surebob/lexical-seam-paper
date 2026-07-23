from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
LOCAL_THEORY_SUMMARY_PATH = ROOT / "results" / "zipf_seam_sign_theory" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_seam_projection_theory"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_projection_common")
reranked = load_module(RERANKED_PATH, "zipf_seam_projection_reranked")


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


def fast_fit_zm_prediction(ranks: np.ndarray, y: np.ndarray) -> np.ndarray:
    max_rank = float(np.max(ranks))
    c_grid = np.concatenate([np.array([0.0], dtype=np.float64), np.geomspace(1e-6, max_rank, 512, dtype=np.float64)])
    best_pred = None
    best_mse = None
    for c in c_grid:
        z = np.log(ranks + c)
        design = np.column_stack([np.ones_like(z), z])
        coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        pred = design @ coeffs
        mse = float(np.mean((pred - y) ** 2))
        if best_mse is None or mse < best_mse:
            best_mse = mse
            best_pred = pred
    assert best_pred is not None
    return best_pred


def seam_curve(ranks: np.ndarray, params: dict) -> np.ndarray:
    sigma = reranked.sigma_curve(ranks, params["k"], params["w"])
    tail_rank = reranked.smooth_tail_local_rank(ranks, params["k"], params["w"])
    head = params["a1"] - params["b1"] * np.log(ranks + max(params["c1"], 0.0))
    tail = params["a2"] - params["b2"] * np.log(tail_rank + max(params["c2"], 0.0))
    return (1.0 - sigma) * (tail - head)


def tangent_matrix(ranks: np.ndarray, params: dict) -> np.ndarray:
    c1 = max(params["c1"], 0.0)
    b1 = params["b1"]
    return np.column_stack(
        [
            np.ones_like(ranks),
            -np.log(ranks + c1),
            -b1 / (ranks + c1),
        ]
    )


def projected_residual(ranks: np.ndarray, params: dict) -> np.ndarray:
    seam = seam_curve(ranks, params)
    G = tangent_matrix(ranks, params)
    coeffs, _, _, _ = np.linalg.lstsq(G, seam, rcond=None)
    return seam - G @ coeffs


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
    smooth_zm_prediction = fast_fit_zm_prediction(dataset["ranks"], smooth_prediction)
    smooth_residual = smooth_prediction - smooth_zm_prediction
    _, smooth_signs = centered_poly_signs(smooth_residual, 50)

    proj_residual = projected_residual(dataset["ranks"], params)
    _, proj_signs = centered_poly_signs(proj_residual, 50)

    local_signs = local_row["local_signs"]
    family = "high_family" if common.get_step2_expr(summary) == common.STEP2_OURS_EXPR else "low_family"
    return {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "family": family,
        "c": float(summary["zm_baseline"]["c"]),
        "local_signs": local_signs,
        "projected_signs": proj_signs,
        "smooth_signs": smooth_signs,
        "empirical_signs": empirical_signs,
        "match_proj_vs_emp": [proj_signs[i] == empirical_signs[i] for i in range(3)],
        "match_smooth_vs_emp": [smooth_signs[i] == empirical_signs[i] for i in range(3)],
        "match_proj_vs_smooth": [proj_signs[i] == smooth_signs[i] for i in range(3)],
    }


def summarize(rows: list[dict], local_summary: dict) -> dict:
    def count_matches(key: str):
        return [int(sum(row[key][i] for row in rows)) for i in range(3)]

    proj_full_match = int(sum(all(row["match_proj_vs_emp"]) for row in rows))
    smooth_full_match = int(sum(all(row["match_smooth_vs_emp"]) for row in rows))
    proj_smooth_full_match = int(sum(all(row["match_proj_vs_smooth"]) for row in rows))
    return {
        "local_baseline_head50_term_matches": local_summary["summary"]["local_match_counts_head50"],
        "local_baseline_head50_full_matches": local_summary["summary"]["local_full_match_head50"],
        "smooth_exact_head50_term_matches": local_summary["summary"]["smooth_match_counts_head50"],
        "smooth_exact_head50_full_matches": local_summary["summary"]["smooth_full_match_head50"],
        "projected_head50_term_matches": count_matches("match_proj_vs_emp"),
        "projected_head50_full_matches": proj_full_match,
        "projected_vs_exact_smooth_term_matches": count_matches("match_proj_vs_smooth"),
        "projected_vs_exact_smooth_full_matches": proj_smooth_full_match,
        "smooth_repeated_head50_full_matches": smooth_full_match,
    }


def build_report(summary: dict, rows: list[dict]) -> str:
    lines = [
        "# Seam Projection Theory",
        "",
        "## First-order Projection Formula",
        "",
        "Let the head-branch ZM be `Y_0(r) = a_1 - b_1 log(r + c_1)` and the smooth seam perturbation be `s(r)`.",
        "",
        "The ZM tangent basis at the head branch is:",
        "",
        "- `g_a(r) = 1`",
        "- `g_b(r) = -log(r + c_1)`",
        "- `g_c(r) = -b_1 / (r + c_1)`",
        "",
        "To first order, re-fitting a single ZM subtracts the tangent-space projection of the seam:",
        "",
        "`s_proj = s - G (G^T G)^(-1) G^T s`",
        "",
        "where `G = [g_a, g_b, g_c]` evaluated over the rank grid.",
        "",
        "## Match Counts At Head-50",
        "",
        f"- raw local seam theorem vs empirical residual: `{summary['local_baseline_head50_term_matches']}` term matches, `{summary['local_baseline_head50_full_matches']}` full matches",
        f"- exact smooth-model misspecification residual vs empirical residual: `{summary['smooth_exact_head50_term_matches']}` term matches, `{summary['smooth_exact_head50_full_matches']}` full matches",
        f"- tangent-projected seam theorem vs empirical residual: `{summary['projected_head50_term_matches']}` term matches, `{summary['projected_head50_full_matches']}` full matches",
        f"- tangent-projected seam vs exact smooth residual: `{summary['projected_vs_exact_smooth_term_matches']}` term matches, `{summary['projected_vs_exact_smooth_full_matches']}` full matches",
        "",
        "## Per-Corpus Head-50 Sign Vectors",
        "",
        "| corpus | family | c | local theorem | projected theorem | exact smooth residual | empirical residual |",
        "| --- | --- | ---: | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['corpus']} | {row['family']} | {row['c']:.3f} | `{row['local_signs']}` | `{row['projected_signs']}` | `{row['smooth_signs']}` | `{row['empirical_signs']}` |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_rows = load_rows()
    local_summary = json.loads(LOCAL_THEORY_SUMMARY_PATH.read_text(encoding="utf-8"))
    local_rows = {row["slug"]: row for row in local_summary["rows"]}
    rows = [analyze_corpus(spec, reranked_rows[spec["slug"]], local_rows[spec["slug"]]) for spec in common.SEARCHED_CORPORA]
    summary = summarize(rows, local_summary)
    payload = {"summary": summary, "rows": rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary, rows), encoding="utf-8")


if __name__ == "__main__":
    main()
