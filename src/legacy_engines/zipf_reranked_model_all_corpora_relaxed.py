import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
CORRECT_MODEL_PATH = ROOT / "zipf_correct_model.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
ORIGINAL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_reranked_model_all_corpora_relaxed"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_reranked_relaxed_common")
correct_model = load_module(CORRECT_MODEL_PATH, "zipf_reranked_relaxed_correct_model")
reranked = load_module(RERANKED_PATH, "zipf_reranked_relaxed_reranked")


def relaxed_bounds():
    lower = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.1], dtype=np.float64)
    upper = np.array([100.0, 5.0, 5000.0, 100.0, 5.0, 5000.0, 2000.0, 5.0], dtype=np.float64)
    return lower, upper


def make_random_initializations(dataset: dict, n_starts: int, seed: int):
    lower, upper = relaxed_bounds()
    rng = np.random.default_rng(seed)
    y = dataset["log_freq"]
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    piece = correct_model.fit_piecewise_zm(dataset, 500)
    anchor = np.array(
        [
            piece["top"]["a"],
            piece["top"]["b"],
            piece["top"]["c"],
            piece["tail"]["a"],
            piece["tail"]["b"],
            min(piece["tail"]["c"], upper[5]),
            500.0,
            0.5,
        ],
        dtype=np.float64,
    )
    anchor = np.clip(anchor, lower + 1e-8, upper - 1e-8)
    local_scale = np.array([8.0, 0.5, 220.0, 8.0, 0.5, 220.0, 260.0, 0.55], dtype=np.float64)

    starts = []
    for idx in range(n_starts):
        if idx < n_starts // 2:
            x0 = anchor + rng.normal(0.0, local_scale)
        else:
            x0 = np.empty(8, dtype=np.float64)
            x0[0] = rng.uniform(max(lower[0], y_min - 3.0), min(upper[0], y_max + 3.0))
            x0[1] = rng.uniform(lower[1], upper[1])
            x0[2] = rng.uniform(lower[2], upper[2])
            x0[3] = rng.uniform(max(lower[3], y_min - 3.0), min(upper[3], y_max + 3.0))
            x0[4] = rng.uniform(lower[4], upper[4])
            x0[5] = rng.uniform(lower[5], upper[5])
            x0[6] = float(np.exp(rng.uniform(np.log(lower[6]), np.log(upper[6]))))
            x0[7] = float(np.exp(rng.uniform(np.log(lower[7]), np.log(upper[7]))))
        starts.append(np.clip(x0, lower + 1e-8, upper - 1e-8))
    return starts


def fit_reranked_model_relaxed(dataset: dict, n_starts: int, max_nfev: int, seed: int):
    ranks = dataset["ranks"]
    y = dataset["log_freq"]
    lower, upper = relaxed_bounds()
    best = None
    tried = []

    def residuals(params):
        return reranked.reranked_prediction(ranks, params) - y

    for idx, x0 in enumerate(make_random_initializations(dataset, n_starts=n_starts, seed=seed), start=1):
        result = least_squares(
            residuals,
            x0=x0,
            bounds=(lower, upper),
            method="trf",
            max_nfev=max_nfev,
        )
        pred = reranked.reranked_prediction(ranks, result.x)
        score = correct_model.rmse(y, pred)
        record = {
            "start_index": idx,
            "rmse": float(score),
            "success": bool(result.success),
            "nfev": int(result.nfev),
            "status": int(result.status),
            "message": result.message,
            "params": [float(v) for v in result.x],
        }
        tried.append(record)
        if best is None or score < best["rmse"]:
            best = {
                "rmse": float(score),
                "prediction": pred,
                "params": [float(v) for v in result.x],
                "start_index": idx,
                "nfev": int(result.nfev),
                "success": bool(result.success),
                "status": int(result.status),
                "message": result.message,
            }

    return best, tried


def load_original_reference():
    payload = json.loads(ORIGINAL_SUMMARY_PATH.read_text(encoding="utf-8"))
    return {row["slug"]: row for row in payload["rows"]}


def bound_contact_counts(rows: list[dict], lower: np.ndarray, upper: np.ndarray):
    labels = ["a1", "b1", "c1", "a2", "b2", "c2", "k", "w"]
    summary = {}
    for idx, label in enumerate(labels):
        lo_tol = max(1e-6, 1e-3 * max(1.0, abs(lower[idx])))
        hi_tol = max(1e-6, 1e-3 * max(1.0, abs(upper[idx])))
        lower_hits = sum(abs(row["params"][label] - float(lower[idx])) <= lo_tol for row in rows)
        upper_hits = sum(abs(row["params"][label] - float(upper[idx])) <= hi_tol for row in rows)
        summary[label] = {"lower": int(lower_hits), "upper": int(upper_hits)}
    return summary


def build_summary(rows: list[dict], original_rows: list[dict]):
    rmse_deltas = [row["rmse_delta_vs_original"] for row in rows]
    frac_deltas = [row["transition_fraction_delta"] for row in rows]
    lower, upper = relaxed_bounds()
    improved = [row for row in rows if row["rmse_delta_vs_original"] < 0.0]
    unchanged = [row for row in rows if abs(row["rmse_delta_vs_original"]) <= 1e-9]
    degraded = [row for row in rows if row["rmse_delta_vs_original"] > 0.0]
    return {
        "rows": rows,
        "counts": {
            "n_rows": len(rows),
            "rmse_improved": len(improved),
            "rmse_unchanged": len(unchanged),
            "rmse_degraded": len(degraded),
        },
        "rmse_delta_vs_original": {
            "mean": float(np.mean(rmse_deltas)),
            "median": float(np.median(rmse_deltas)),
            "min": float(np.min(rmse_deltas)),
            "max": float(np.max(rmse_deltas)),
            "std": float(np.std(rmse_deltas)),
        },
        "transition_fraction_delta_vs_original": {
            "mean": float(np.mean(frac_deltas)),
            "median": float(np.median(frac_deltas)),
            "min": float(np.min(frac_deltas)),
            "max": float(np.max(frac_deltas)),
            "std": float(np.std(frac_deltas)),
        },
        "bound_contacts_relaxed": bound_contact_counts(rows, lower, upper),
        "bound_contacts_original_box": bound_contact_counts(original_rows, np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.1], dtype=np.float64), np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 1000.0, 3.0], dtype=np.float64)),
    }


def build_report(summary: dict):
    counts = summary["counts"]
    rmse = summary["rmse_delta_vs_original"]
    frac = summary["transition_fraction_delta_vs_original"]
    bound_relaxed = summary["bound_contacts_relaxed"]
    bound_orig = summary["bound_contacts_original_box"]
    rows = sorted(summary["rows"], key=lambda row: row["relaxed_rmse"])
    lines = [
        "# Reranked Smooth Model on All 25 Corpora with Relaxed Bounds",
        "",
        "- Relaxed bounds:",
        "  - `b1, b2 in [0.5, 5.0]`",
        "  - `c1, c2 in [0, 5000]`",
        "  - `k in [20, 2000]`",
        "  - `w in [0.1, 5.0]`",
        "",
        f"- corpora with improved RMSE vs original tighter-bounds run: `{counts['rmse_improved']}` / `{counts['n_rows']}`",
        f"- corpora unchanged to numerical tolerance: `{counts['rmse_unchanged']}` / `{counts['n_rows']}`",
        f"- corpora with worse RMSE: `{counts['rmse_degraded']}` / `{counts['n_rows']}`",
        f"- mean RMSE delta vs original: `{rmse['mean']:.12f}`",
        f"- median RMSE delta vs original: `{rmse['median']:.12f}`",
        f"- min/max RMSE delta vs original: `{rmse['min']:.12f}` / `{rmse['max']:.12f}`",
        f"- mean transition-fraction delta vs original: `{frac['mean']:.12f}`",
        f"- median transition-fraction delta vs original: `{frac['median']:.12f}`",
        f"- min/max transition-fraction delta vs original: `{frac['min']:.12f}` / `{frac['max']:.12f}`",
        "",
        "## Bound-contact comparison",
        "",
        "| param | original lower | original upper | relaxed lower | relaxed upper |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for label in ["b1", "b2", "c1", "c2", "k", "w"]:
        lines.append(
            f"| {label} | {bound_orig[label]['lower']} | {bound_orig[label]['upper']} | {bound_relaxed[label]['lower']} | {bound_relaxed[label]['upper']} |"
        )
    lines.extend(
        [
            "",
            "| corpus | original RMSE | relaxed RMSE | delta | original frac | relaxed frac | delta frac | original k | relaxed k | original w | relaxed w |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['original_rmse']:.12f} | {row['relaxed_rmse']:.12f} | {row['rmse_delta_vs_original']:.12f} | {row['original_transition_fraction']:.6f} | {row['relaxed_transition_fraction']:.6f} | {row['transition_fraction_delta']:.6f} | {row['original_k']:.3f} | {row['relaxed_k']:.3f} | {row['original_w']:.3f} | {row['relaxed_w']:.3f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    original_ref = load_original_reference()
    rows = []
    original_rows = []

    for idx, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        dataset = common.build_zipf_dataset(common.corpus_path(spec))
        original = original_ref[spec["slug"]]
        best, tried = fit_reranked_model_relaxed(dataset, n_starts=100, max_nfev=12000, seed=20260800 + idx)
        params = reranked.summarize_params(best["params"], dataset["unique_words"])

        row = {
            "slug": spec["slug"],
            "name": spec["name"],
            "original_rmse": float(original["reranked_rmse"]),
            "relaxed_rmse": float(best["rmse"]),
            "rmse_delta_vs_original": float(best["rmse"] - original["reranked_rmse"]),
            "original_transition_fraction": float(original["transition_fraction"]),
            "relaxed_transition_fraction": float(params["transition_fraction"]),
            "transition_fraction_delta": float(params["transition_fraction"] - original["transition_fraction"]),
            "original_k": float(original["k"]),
            "relaxed_k": float(params["k"]),
            "original_w": float(original["w"]),
            "relaxed_w": float(params["w"]),
            "best_start_index": int(best["start_index"]),
            "best_nfev": int(best["nfev"]),
            "params": params,
            "tried_count": len(tried),
        }
        rows.append(row)
        original_rows.append({"params": original["params"]})

    summary = build_summary(rows, original_rows)
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")


if __name__ == "__main__":
    main()
