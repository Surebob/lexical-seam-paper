from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
MOE_SUMMARY_PATH = ROOT / "results" / "zipf_moezipf_comparison" / "summary.json"
BIC_SUMMARY_PATH = ROOT / "results" / "zipf_bic_comparison" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_continuous_piecewise"

B_LO = 0.5
B_HI = 5.0
C_LO = 0.0
C_HI = 5000.0


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_continuous_piecewise_common")


def bic_from_rmse(rmse: float, p: int, n: int) -> float:
    mse = max(float(rmse) ** 2, 1e-300)
    return float(p * math.log(n) + n * math.log(mse))


def load_reference_rows():
    moe_rows = {row["slug"]: row for row in json.loads(MOE_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]}
    bic_rows = {row["slug"]: row for row in json.loads(BIC_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]}
    return moe_rows, bic_rows


def valid_k_bounds(n: int) -> tuple[int, int]:
    lo = 5
    hi = max(lo, n - 5)
    return lo, hi


def coarse_k_grid(n: int) -> list[int]:
    lo, hi = valid_k_bounds(n)
    if hi <= lo:
        return [lo]
    candidates = {
        lo,
        hi,
        int(round(math.sqrt(n))),
        10,
        20,
        50,
        100,
        200,
        500,
        1000,
    }
    logs = np.exp(np.linspace(math.log(lo), math.log(hi), 20))
    candidates.update(int(round(v)) for v in logs)
    return sorted(k for k in candidates if lo <= k <= hi)


def refine_k_grid(n: int, seeds: list[int]) -> list[int]:
    lo, hi = valid_k_bounds(n)
    candidates = set()
    for k in seeds:
        for delta in range(-15, 16):
            kk = k + delta
            if lo <= kk <= hi:
                candidates.add(kk)
        local_lo = max(lo, int(round(0.8 * k)))
        local_hi = min(hi, int(round(1.2 * k)))
        if local_hi > local_lo:
            linear = np.linspace(local_lo, local_hi, 15)
            candidates.update(int(round(v)) for v in linear)
    return sorted(candidates)


def offsets_for_theta(ranks: np.ndarray, k: int, theta: np.ndarray) -> np.ndarray | None:
    b1, c1, b2, c2 = [float(v) for v in theta]
    if not (B_LO <= b1 <= B_HI and B_LO <= b2 <= B_HI and C_LO <= c1 <= C_HI and C_LO <= c2 <= C_HI):
        return None
    kf = float(k)
    with np.errstate(all="ignore"):
        left = -b1 * np.log(ranks + c1)
        continuity_shift = -b1 * math.log(kf + c1) + b2 * math.log(kf + c2)
        right = continuity_shift - b2 * np.log(ranks + c2)
    offsets = np.where(ranks <= kf, left, right)
    if not np.all(np.isfinite(offsets)):
        return None
    return offsets


def prediction_from_theta(ranks: np.ndarray, y: np.ndarray, k: int, theta: np.ndarray):
    offsets = offsets_for_theta(ranks, k, theta)
    if offsets is None:
        return None
    a1 = float(np.mean(y - offsets))
    pred = a1 + offsets
    if not np.all(np.isfinite(pred)):
        return None
    return a1, pred


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    diff = np.asarray(y_true, dtype=np.float64) - np.asarray(y_pred, dtype=np.float64)
    return float(math.sqrt(float(np.mean(diff * diff))))


def fit_zm_segment(ranks: np.ndarray, y: np.ndarray):
    max_rank = min(float(np.max(ranks)), C_HI)

    def objective(c: float) -> float:
        z = np.log(ranks + c)
        _, _, _, mse = common._solve_affine(z, y)
        return mse

    best_c = 0.0
    best_mse = objective(0.0)

    if max_rank > 0.0:
        opt = minimize_scalar(objective, bounds=(0.0, max_rank), method="bounded", options={"xatol": 1e-3})
        if opt.success and float(opt.fun) < best_mse:
            best_c = float(opt.x)
            best_mse = float(opt.fun)

    z = np.log(ranks + best_c)
    intercept, slope, pred, mse = common._solve_affine(z, y)
    return {
        "a": float(intercept),
        "b": float(-slope),
        "c": float(best_c),
        "mse": float(mse),
        "rmse": float(math.sqrt(max(float(mse), 0.0))),
        "prediction": pred,
    }


def fit_for_fixed_k(dataset: dict, k: int, full_fit: dict):
    ranks = dataset["ranks"]
    y = dataset["log_freq"]
    lo_idx = k
    top_fit = fit_zm_segment(ranks[:lo_idx], y[:lo_idx])
    tail_fit = fit_zm_segment(ranks[lo_idx:], y[lo_idx:])

    starts = [
        np.array([top_fit["b"], min(top_fit["c"], C_HI), tail_fit["b"], min(tail_fit["c"], C_HI)], dtype=np.float64),
        np.array([full_fit["b"], min(full_fit["c"], C_HI), full_fit["b"], min(full_fit["c"], C_HI)], dtype=np.float64),
        np.array([1.5, 0.0, 1.5, 0.0], dtype=np.float64),
    ]
    bounds = (
        np.array([B_LO, C_LO, B_LO, C_LO], dtype=np.float64),
        np.array([B_HI, C_HI, B_HI, C_HI], dtype=np.float64),
    )

    best = None

    def residuals(theta: np.ndarray):
        result = prediction_from_theta(ranks, y, k, theta)
        if result is None:
            return np.full_like(y, 1e6, dtype=np.float64)
        _, pred = result
        return pred - y

    for idx, x0 in enumerate(starts, start=1):
        x0 = np.clip(x0, bounds[0] + 1e-8, bounds[1] - 1e-8)
        result = least_squares(residuals, x0=x0, bounds=bounds, method="trf", max_nfev=4000)
        pred_result = prediction_from_theta(ranks, y, k, result.x)
        if pred_result is None:
            continue
        a1, pred = pred_result
        score = rmse(y, pred)
        candidate = {
            "k": int(k),
            "a1": float(a1),
            "b1": float(result.x[0]),
            "c1": float(result.x[1]),
            "b2": float(result.x[2]),
            "c2": float(result.x[3]),
            "rmse": float(score),
            "prediction": pred,
            "start_index": idx,
            "nfev": int(result.nfev),
            "success": bool(result.success),
        }
        if best is None or candidate["rmse"] < best["rmse"]:
            best = candidate
    return best


def fit_continuous_piecewise(dataset: dict):
    n = dataset["unique_words"]
    full_fit = fit_zm_segment(dataset["ranks"], dataset["log_freq"])
    coarse = coarse_k_grid(n)
    coarse_fits = []
    for k in coarse:
        fit = fit_for_fixed_k(dataset, k, full_fit)
        if fit is not None:
            coarse_fits.append(fit)
    if not coarse_fits:
        raise RuntimeError("No valid coarse continuous-piecewise fit found")
    coarse_fits.sort(key=lambda row: row["rmse"])

    seed_ks = [row["k"] for row in coarse_fits[:2]]
    refined = refine_k_grid(n, seed_ks)
    best = coarse_fits[0]
    for k in refined:
        fit = fit_for_fixed_k(dataset, k, full_fit)
        if fit is not None and fit["rmse"] < best["rmse"]:
            best = fit
    return best


def winner_counts(rows: list[dict]) -> dict:
    counts = {
        "zipf": 0,
        "moezipf": 0,
        "piecewise_k500": 0,
        "continuous_piecewise": 0,
        "reranked_7param_sqrtv": 0,
        "reranked_8param": 0,
    }
    for row in rows:
        counts[row["best_model"]] += 1
    return counts


def write_partial_outputs(rows: list[dict]) -> None:
    counts = winner_counts(rows)
    summary = {
        "rows": rows,
        "winner_counts": counts,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(rows, counts), encoding="utf-8")
    write_csv(rows, OUTDIR / "continuous_piecewise_table.csv")


def build_rows():
    moe_map, bic_map = load_reference_rows()
    rows = []
    for spec in common.SEARCHED_CORPORA:
        print(f"[continuous-piecewise] fitting {spec['slug']}")
        dataset = common.build_zipf_dataset(common.corpus_path(spec))
        fit = fit_continuous_piecewise(dataset)
        vocab = int(dataset["unique_words"])
        cont_bic = bic_from_rmse(fit["rmse"], 6, vocab)
        ref_moe = moe_map[spec["slug"]]
        ref_bic = bic_map[spec["slug"]]
        all_bics = {
            "zipf": float(ref_moe["zipf_rank_bic"]),
            "moezipf": float(ref_moe["moe_rank_bic"]),
            "piecewise_k500": float(ref_bic["piecewise_k500_bic"]),
            "continuous_piecewise": float(cont_bic),
            "reranked_7param_sqrtv": float(ref_bic["reranked_7param_sqrtv_bic"]),
            "reranked_8param": float(ref_bic["reranked_8param_bic"]),
        }
        best_model = min(all_bics, key=all_bics.get)
        rows.append(
            {
                "slug": spec["slug"],
                "name": spec["name"],
                "vocab_size": vocab,
                "continuous_piecewise_rmse": float(fit["rmse"]),
                "continuous_piecewise_k": int(fit["k"]),
                "continuous_piecewise_params": {
                    "a1": float(fit["a1"]),
                    "b1": float(fit["b1"]),
                    "c1": float(fit["c1"]),
                    "b2": float(fit["b2"]),
                    "c2": float(fit["c2"]),
                },
                "zipf_bic": all_bics["zipf"],
                "moezipf_bic": all_bics["moezipf"],
                "piecewise_k500_bic": all_bics["piecewise_k500"],
                "continuous_piecewise_bic": all_bics["continuous_piecewise"],
                "reranked_7param_sqrtv_bic": all_bics["reranked_7param_sqrtv"],
                "reranked_8param_bic": all_bics["reranked_8param"],
                "best_model": best_model,
            }
        )
        rows = sorted(rows, key=lambda row: row["slug"])
        write_partial_outputs(rows)
    return sorted(rows, key=lambda row: row["slug"])


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "slug",
                "name",
                "vocab_size",
                "continuous_piecewise_rmse",
                "continuous_piecewise_k",
                "zipf_bic",
                "moezipf_bic",
                "piecewise_k500_bic",
                "continuous_piecewise_bic",
                "reranked_7param_sqrtv_bic",
                "reranked_8param_bic",
                "best_model",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_report(rows: list[dict], counts: dict) -> str:
    labels = {
        "zipf": "Single Zipf",
        "moezipf": "MOEZipf",
        "piecewise_k500": "Hard piecewise ZM (K=500)",
        "continuous_piecewise": "Continuous piecewise ZM",
        "reranked_7param_sqrtv": "Smooth 7-param k=sqrt(V)",
        "reranked_8param": "Smooth 8-param",
    }
    lines = [
        "# Continuous Piecewise ZM Comparison",
        "",
        "- Continuous piecewise model: hard breakpoint with C0 continuity; `a2` is determined by continuity at `K`.",
        "- BIC uses the same RMSE-based formula as the earlier model-comparison runs: `p * ln(n) + n * ln(MSE)`.",
        "",
        "## Winner Counts",
        "",
    ]
    for key in ["zipf", "moezipf", "piecewise_k500", "continuous_piecewise", "reranked_7param_sqrtv", "reranked_8param"]:
        lines.append(f"- {labels[key]}: `{counts[key]}`")
    lines.extend(
        [
            "",
            "| Corpus | Zipf | MOEZipf | Hard piecewise | Continuous piecewise | Smooth 7-param | Smooth 8-param | Best BIC |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['zipf_bic']:.3f} | {row['moezipf_bic']:.3f} | {row['piecewise_k500_bic']:.3f} | {row['continuous_piecewise_bic']:.3f} | {row['reranked_7param_sqrtv_bic']:.3f} | {row['reranked_8param_bic']:.3f} | {labels[row['best_model']]} |"
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    write_partial_outputs(rows)


if __name__ == "__main__":
    main()
