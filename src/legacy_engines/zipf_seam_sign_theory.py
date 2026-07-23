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
OUTDIR = ROOT / "results" / "zipf_seam_sign_theory"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_seam_sign_common")
reranked = load_module(RERANKED_PATH, "zipf_seam_sign_reranked")


def load_rows():
    return {row["slug"]: row for row in json.loads(RERANKED_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]}


def head_branch(ranks: np.ndarray, params: dict) -> np.ndarray:
    return params["a1"] - params["b1"] * np.log(ranks + max(params["c1"], 0.0))


def centered_poly_signs(values: np.ndarray, head_k: int) -> tuple[list[float], list[int]]:
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


def tau_derivatives_at_head(k: float, w: float) -> dict[str, float]:
    tau0 = 1.0 / (1.0 + math.exp(math.log(k) / w))
    tau1 = tau0 * (1.0 - tau0) / w
    tau2 = tau0 * (1.0 - tau0) * (1.0 - 2.0 * tau0) / (w * w)
    tau3 = tau0 * (1.0 - tau0) * (1.0 - 6.0 * tau0 + 6.0 * tau0 * tau0) / (w**3)
    return {"tau0": tau0, "tau1": tau1, "tau2": tau2, "tau3": tau3}


def delta_value(z: float, params: dict) -> float:
    r = math.exp(z)
    head = params["a1"] - params["b1"] * math.log(r + max(params["c1"], 0.0))
    tail_rank = reranked.smooth_tail_local_rank(np.array([r], dtype=np.float64), params["k"], params["w"])[0]
    tail = params["a2"] - params["b2"] * math.log(tail_rank + max(params["c2"], 0.0))
    return tail - head


def delta_derivatives_at_head(params: dict) -> dict[str, float]:
    eps = 1e-5
    d0 = delta_value(0.0, params)
    d1 = (delta_value(eps, params) - delta_value(-eps, params)) / (2.0 * eps)
    d2 = (delta_value(eps, params) - 2.0 * d0 + delta_value(-eps, params)) / (eps**2)
    d3 = (
        delta_value(2.0 * eps, params)
        - 2.0 * delta_value(eps, params)
        + 2.0 * delta_value(-eps, params)
        - delta_value(-2.0 * eps, params)
    ) / (2.0 * eps**3)
    return {"delta0": d0, "delta1": d1, "delta2": d2, "delta3": d3}


def local_seam_coefficients(params: dict) -> dict[str, float]:
    tau = tau_derivatives_at_head(params["k"], params["w"])
    delta = delta_derivatives_at_head(params)
    a1 = tau["tau1"] * delta["delta0"] + tau["tau0"] * delta["delta1"]
    a2 = 0.5 * (tau["tau2"] * delta["delta0"] + 2.0 * tau["tau1"] * delta["delta1"] + tau["tau0"] * delta["delta2"])
    a3 = (
        tau["tau3"] * delta["delta0"]
        + 3.0 * tau["tau2"] * delta["delta1"]
        + 3.0 * tau["tau1"] * delta["delta2"]
        + tau["tau0"] * delta["delta3"]
    ) / 6.0
    return {**tau, **delta, "a1_local": a1, "a2_local": a2, "a3_local": a3}


def analyze_corpus(spec: dict, reranked_row: dict) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    summary = common.load_enriched_summary(spec)
    params = reranked_row["params"]
    param_vec = np.array(
        [params["a1"], params["b1"], params["c1"], params["a2"], params["b2"], params["c2"], params["k"], params["w"]],
        dtype=np.float64,
    )

    empirical_residual = dataset["log_freq"] - common.zm_prediction(summary, dataset["ranks"])
    empirical_10, empirical_signs_10 = centered_poly_signs(empirical_residual, 10)
    empirical_50, empirical_signs_50 = centered_poly_signs(empirical_residual, 50)

    smooth_prediction = reranked.reranked_prediction(dataset["ranks"], param_vec)
    smooth_zm_prediction = fast_fit_zm_prediction(dataset["ranks"], smooth_prediction)
    smooth_misspec_residual = smooth_prediction - smooth_zm_prediction
    smooth_50, smooth_signs_50 = centered_poly_signs(smooth_misspec_residual, 50)

    seam_curve = smooth_prediction - head_branch(dataset["ranks"], params)
    seam_50, seam_signs_50 = centered_poly_signs(seam_curve, 50)

    local = local_seam_coefficients(params)
    local_signs = [int(np.sign(local["a1_local"])), int(np.sign(local["a2_local"])), int(np.sign(local["a3_local"]))]
    family = "high_family" if common.get_step2_expr(summary) == common.STEP2_OURS_EXPR else "low_family"

    return {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "family": family,
        "c": float(summary["zm_baseline"]["c"]),
        "tau0": local["tau0"],
        "tau1": local["tau1"],
        "tau2": local["tau2"],
        "tau3": local["tau3"],
        "delta0": local["delta0"],
        "delta1": local["delta1"],
        "delta2": local["delta2"],
        "delta3": local["delta3"],
        "local_signs": local_signs,
        "empirical_coeffs_10": empirical_10,
        "empirical_signs_10": empirical_signs_10,
        "empirical_coeffs_50": empirical_50,
        "empirical_signs_50": empirical_signs_50,
        "smooth_coeffs_50": smooth_50,
        "smooth_signs_50": smooth_signs_50,
        "seam_coeffs_50": seam_50,
        "seam_signs_50": seam_signs_50,
        "match_local_vs_emp10": [local_signs[i] == empirical_signs_10[i] for i in range(3)],
        "match_local_vs_emp50": [local_signs[i] == empirical_signs_50[i] for i in range(3)],
        "match_smooth_vs_emp50": [smooth_signs_50[i] == empirical_signs_50[i] for i in range(3)],
    }


def summarize(rows: list[dict]) -> dict:
    def count_matches(key: str):
        return [int(sum(row[key][i] for row in rows)) for i in range(3)]

    all_negative_local = int(sum(row["local_signs"] == [-1, -1, -1] for row in rows))
    smooth_full_match = int(sum(all(row["match_smooth_vs_emp50"]) for row in rows))
    local_full_match_10 = int(sum(all(row["match_local_vs_emp10"]) for row in rows))
    local_full_match_50 = int(sum(all(row["match_local_vs_emp50"]) for row in rows))
    return {
        "local_match_counts_head10": count_matches("match_local_vs_emp10"),
        "local_match_counts_head50": count_matches("match_local_vs_emp50"),
        "smooth_match_counts_head50": count_matches("match_smooth_vs_emp50"),
        "local_full_match_head10": local_full_match_10,
        "local_full_match_head50": local_full_match_50,
        "smooth_full_match_head50": smooth_full_match,
        "local_all_negative_count": all_negative_local,
    }


def build_report(summary: dict, rows: list[dict]) -> str:
    lines = [
        "# Seam Sign Theory",
        "",
        "## Analytic Setup",
        "",
        "Write the smooth model in head log-rank coordinate `z = log(rank)` as:",
        "",
        "`Y(z) = Y_head(z) + tau(z) * Delta(z)`",
        "",
        "where `tau(z)` is the logistic tail activation and `Delta(z) = Y_tail(z) - Y_head(z)`.",
        "",
        "At the head (`z = 0`), the first three centered Taylor coefficients of the seam term are:",
        "",
        "- `a1 = tau'(0) Delta(0) + tau(0) Delta'(0)`",
        "- `a2 = 1/2 [ tau''(0) Delta(0) + 2 tau'(0) Delta'(0) + tau(0) Delta''(0) ]`",
        "- `a3 = 1/6 [ tau'''(0) Delta(0) + 3 tau''(0) Delta'(0) + 3 tau'(0) Delta''(0) + tau(0) Delta'''(0) ]`",
        "",
        "For the logistic gate, if `tau(0)` is small enough then `tau'(0), tau''(0), tau'''(0)` are all positive. In that gate-dominated regime, if `Delta(0) < 0` and the derivatives of `Delta` are not too large, the local sign law is simply:",
        "",
        "- `a1 < 0`",
        "- `a2 < 0`",
        "- `a3 < 0`",
        "",
        "## What The Theory Predicts",
        "",
        f"- corpora with local seam sign exactly `(-,-,-)`: `{summary['local_all_negative_count']}` / `{len(rows)}`",
        "",
        "## Sign Matches Against Empirical Single-ZM Residual",
        "",
        f"- local theorem vs empirical head-10 signs: `{summary['local_match_counts_head10']}` / 25 for `(a1,a2,a3)`",
        f"- local theorem vs empirical head-50 signs: `{summary['local_match_counts_head50']}` / 25 for `(a1,a2,a3)`",
        f"- full 3-sign match count at head-10: `{summary['local_full_match_head10']}` / 25",
        f"- full 3-sign match count at head-50: `{summary['local_full_match_head50']}` / 25",
        "",
        "## Sign Matches After Re-fitting ZM To The Smooth Model",
        "",
        f"- smooth-model-predicted residual vs empirical head-50 signs: `{summary['smooth_match_counts_head50']}` / 25 for `(a1,a2,a3)`",
        f"- full 3-sign match count at head-50: `{summary['smooth_full_match_head50']}` / 25",
        "",
        "## Per-Corpus Head-50 Sign Vectors",
        "",
        "| corpus | family | c | local theorem | smooth-model residual | empirical residual |",
        "| --- | --- | ---: | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['corpus']} | {row['family']} | {row['c']:.3f} | `{row['local_signs']}` | `{row['smooth_signs_50']}` | `{row['empirical_signs_50']}` |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_rows = load_rows()
    rows = [analyze_corpus(spec, reranked_rows[spec["slug"]]) for spec in common.SEARCHED_CORPORA]
    summary = summarize(rows)
    payload = {"summary": summary, "rows": rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary, rows), encoding="utf-8")


if __name__ == "__main__":
    main()
