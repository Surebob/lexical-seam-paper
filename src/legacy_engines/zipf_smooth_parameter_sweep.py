from __future__ import annotations

import importlib.util
import json
import math
import statistics
from collections import Counter
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
except ModuleNotFoundError:
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
PHASE_COORDINATE_SUMMARY_PATH = ROOT / "results" / "zipf_phase_coordinate" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_smooth_parameter_sweep"

FRAC_GRID = np.array([0.45, 0.50, 0.55, 0.60, 0.65, 0.70], dtype=np.float64)
W_GRID = np.array([0.20, 0.50, 0.80, 1.10, 1.40, 1.70], dtype=np.float64)
B2_GRID = np.array([1.60, 1.80, 2.00, 2.20, 2.40, 2.60], dtype=np.float64)
LAMBDA_GRID = np.linspace(0.0, 1.0, 21)
HEAD_TOP_K = 100
ANGLE_HEAD_K = 200
FAST_C_GRID_SIZE = 256


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_smooth_sweep_common")
reranked = load_module(RERANKED_PATH, "zipf_smooth_sweep_reranked")


def load_reranked_rows():
    rows = json.loads(RERANKED_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return rows


def load_phase_low_rows():
    return json.loads(PHASE_COORDINATE_SUMMARY_PATH.read_text(encoding="utf-8"))["low_rows"]


def exp_bregman(x: np.ndarray) -> np.ndarray:
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def xpow_minus_sqrt(x: np.ndarray) -> np.ndarray:
    return np.power(x, x) - np.sqrt(x)


def is_bregman(x: np.ndarray) -> np.ndarray:
    return (x - 1.0) - np.log(x)


FORMULAS = {
    "is": is_bregman,
    "exp": exp_bregman,
    "xpow": xpow_minus_sqrt,
}
WINNER_CODE = {"is": 0, "exp": 1, "xpow": 2}
WINNER_COLOR = ListedColormap(["#2563eb", "#059669", "#dc2626"]) if plt is not None else None


def representative_lowc_template() -> dict:
    rows = load_reranked_rows()
    by_slug = {row["slug"]: row for row in rows}
    low_specs = []
    for spec in common.SEARCHED_CORPORA:
        summary = common.load_enriched_summary(spec)
        if common.get_step2_expr(summary) == common.STEPL2_LOWC_EXPR:
            low_specs.append(by_slug[spec["slug"]])
    params = {key: statistics.median([row["params"][key] for row in low_specs]) for key in ["a1", "b1", "c1", "a2", "b2", "c2", "w", "transition_fraction"]}
    vocab_size = int(round(statistics.median([row["vocab_size"] for row in low_specs])))
    return {"vocab_size": vocab_size, **params}


def fast_fit_zm_prediction(ranks: np.ndarray, y: np.ndarray) -> dict:
    max_rank = float(np.max(ranks))
    c_grid = np.concatenate(
        [np.array([0.0], dtype=np.float64), np.geomspace(1e-6, max_rank, FAST_C_GRID_SIZE - 1, dtype=np.float64)]
    )
    best = None
    for c in c_grid:
        z = np.log(ranks + c)
        design = np.column_stack([np.ones_like(z), z])
        coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        pred = design @ coeffs
        mse = float(np.mean((pred - y) ** 2))
        if best is None or mse < best["mse"]:
            best = {"a": float(coeffs[0]), "b": float(-coeffs[1]), "c": float(c), "prediction": pred, "mse": mse}
    assert best is not None
    return best


def weighted_rmse(target: np.ndarray, pred: np.ndarray, weights: np.ndarray) -> float:
    diff = target - pred
    return float(np.sqrt(np.sum(weights * diff * diff) / np.sum(weights)))


def angle_fit(target: np.ndarray, f_exp: np.ndarray, f_xpow: np.ndarray, head_k: int = ANGLE_HEAD_K) -> dict:
    k = min(head_k, len(target))
    y = target[:k] - np.mean(target[:k])
    e = f_exp[:k] - np.mean(f_exp[:k])
    xw = f_xpow[:k] - np.mean(f_xpow[:k])
    e /= np.linalg.norm(e)
    xw /= np.linalg.norm(xw)
    X = np.column_stack([e, xw])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    sse = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else 1.0
    return {
        "theta_deg": float(np.degrees(np.arctan2(beta[1], beta[0]))),
        "beta_exp": float(beta[0]),
        "beta_xpow": float(beta[1]),
        "r2": float(r2),
    }


def winner_sweep(target: np.ndarray, curves: dict[str, np.ndarray], top_k: int = HEAD_TOP_K) -> dict:
    n = len(target)
    top_k = min(top_k, n)
    head = np.zeros(n, dtype=np.float64)
    head[:top_k] = 1.0
    rows = []
    flip_lambda = None
    for lam in LAMBDA_GRID:
        weights = (1.0 - lam) * np.ones(n, dtype=np.float64) + lam * head
        scores = {name: weighted_rmse(target, curve, weights) for name, curve in curves.items()}
        winner = min(scores, key=lambda key: (scores[key], key))
        if flip_lambda is None and scores["xpow"] < scores["exp"]:
            flip_lambda = float(lam)
        rows.append({"lambda": float(lam), "winner": winner, "scores": scores})
    return {"rows": rows, "flip_lambda_exp_to_xpow": flip_lambda}


def standardize_regression(rows: list[dict], y_key: str) -> dict:
    X = np.array([[row["frac"], row["w"], row["b2"]] for row in rows], dtype=np.float64)
    y = np.array([row[y_key] for row in rows], dtype=np.float64)
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    y = (y - y.mean()) / y.std()
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    return {"frac": float(beta[0]), "w": float(beta[1]), "b2": float(beta[2])}


def correlation(xs: list[float], ys: list[float]) -> float:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    if len(x) < 2 or float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def build_smooth_curve(template: dict, frac: float, w: float, b2: float) -> dict:
    V = int(template["vocab_size"])
    ranks = np.arange(1, V + 1, dtype=np.float64)
    a1 = float(template["a1"])
    b1 = float(template["b1"])
    c1 = float(template["c1"])
    c2 = float(template["c2"])
    k = float(V**frac)
    sigma = reranked.sigma_curve(ranks, k, w)
    head = a1 - b1 * np.log(ranks + c1)
    tail_rank = reranked.smooth_tail_local_rank(ranks, k, w)
    tail_rank_k = reranked.smooth_tail_local_rank(np.array([k], dtype=np.float64), k, w)[0]
    head_k = a1 - b1 * math.log(k + c1)
    a2 = float(head_k + b2 * math.log(tail_rank_k + c2))
    tail = a2 - b2 * np.log(tail_rank + c2)
    y = sigma * head + (1.0 - sigma) * tail
    return {
        "ranks": ranks,
        "log_rank": np.log(ranks),
        "log_freq": y,
        "params": {"a1": a1, "b1": b1, "c1": c1, "a2": a2, "b2": b2, "c2": c2, "k": k, "w": w},
    }


def analyze_curve(curve: dict, frac: float, w: float, b2: float) -> dict:
    zm_fit = fast_fit_zm_prediction(curve["ranks"], curve["log_freq"])
    target = curve["log_freq"] - zm_fit["prediction"]
    x = common.normalize_x(curve["log_rank"], 0.05, 1.0)
    formula_curves = {name: fn(x) for name, fn in FORMULAS.items()}
    sweep = winner_sweep(target, formula_curves)
    angle = angle_fit(target, formula_curves["exp"], formula_curves["xpow"])
    return {
        "frac": float(frac),
        "w": float(w),
        "b2": float(b2),
        "delta_b": float(curve["params"]["b1"] - b2),
        "k": float(curve["params"]["k"]),
        "zm_c": float(zm_fit["c"]),
        "theta_deg": angle["theta_deg"],
        "angle_r2": angle["r2"],
        "winner_full": sweep["rows"][0]["winner"],
        "winner_top100": sweep["rows"][-1]["winner"],
        "flip_lambda": 1.05 if sweep["flip_lambda_exp_to_xpow"] is None else float(sweep["flip_lambda_exp_to_xpow"]),
    }


def plot_theta_scatter(rows: list[dict], outpath: Path):
    if plt is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    sc = axes[0].scatter([row["delta_b"] for row in rows], [row["theta_deg"] for row in rows], c=[row["frac"] for row in rows], cmap="viridis", alpha=0.85)
    axes[0].set_xlabel("delta_b = b1 - b2")
    axes[0].set_ylabel("theta (deg)")
    axes[0].set_title("Synthetic theta vs slope contrast")
    axes[0].grid(True, alpha=0.25)
    fig.colorbar(sc, ax=axes[0], label="transition fraction")

    sc2 = axes[1].scatter([row["w"] for row in rows], [row["flip_lambda"] for row in rows], c=[row["b2"] for row in rows], cmap="plasma", alpha=0.85)
    axes[1].set_xlabel("w")
    axes[1].set_ylabel("lambda flip exp->xpow")
    axes[1].set_title("Flip lambda vs seam width")
    axes[1].grid(True, alpha=0.25)
    fig.colorbar(sc2, ax=axes[1], label="b2")
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_top100_heatmaps(rows: list[dict], outpath: Path):
    if plt is None:
        return
    b2_values = sorted(set(row["b2"] for row in rows))
    frac_values = sorted(set(row["frac"] for row in rows))
    w_values = sorted(set(row["w"] for row in rows))
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for ax, b2 in zip(axes.flat, b2_values):
        mat = np.zeros((len(w_values), len(frac_values)), dtype=np.int64)
        for i, w in enumerate(w_values):
            for j, frac in enumerate(frac_values):
                row = next(item for item in rows if item["b2"] == b2 and item["w"] == w and item["frac"] == frac)
                mat[i, j] = WINNER_CODE[row["winner_top100"]]
        ax.imshow(mat, aspect="auto", cmap=WINNER_COLOR, interpolation="nearest", origin="lower")
        ax.set_xticks(range(len(frac_values)))
        ax.set_xticklabels([f"{v:.2f}" for v in frac_values], rotation=90)
        ax.set_yticks(range(len(w_values)))
        ax.set_yticklabels([f"{v:.2f}" for v in w_values])
        ax.set_title(f"b2 = {b2:.2f}")
        ax.set_xlabel("transition fraction")
        ax.set_ylabel("w")
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(color="#2563eb", label="IS"),
        Patch(color="#059669", label="exp"),
        Patch(color="#dc2626", label="xpow"),
    ]
    fig.legend(handles=legend_handles, loc="upper right")
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def build_report(template: dict, rows: list[dict], empirical_low_rows: list[dict]) -> str:
    full_counts = Counter(row["winner_full"] for row in rows)
    top_counts = Counter(row["winner_top100"] for row in rows)
    theta_reg = standardize_regression(rows, "theta_deg")
    flip_reg = standardize_regression(rows, "flip_lambda")
    lines = [
        "# Smooth Parameter Sweep",
        "",
        "Synthetic low-c smooth curves were generated from a continuity-normalized template using the median low-c English parameters, while varying:",
        "",
        "- transition fraction `log(k)/log(V)`",
        "- sigmoid width `w`",
        "- tail slope `b2` (with head slope `b1` fixed at the low-c median)",
        "",
        f"- template vocab size: `{template['vocab_size']}`",
        f"- template medians: `a1={template['a1']:.6f}, b1={template['b1']:.6f}, c1={template['c1']:.6f}, c2={template['c2']:.6f}`",
        "",
        "## Winner Counts Across The Synthetic Grid",
        "",
        f"- full-RMSE winners: `{dict(full_counts)}`",
        f"- top-100 winners: `{dict(top_counts)}`",
        "",
        "## Parameter Influence Rankings",
        "",
        f"- corr(theta, frac): `{correlation([row['frac'] for row in rows], [row['theta_deg'] for row in rows]):.6f}`",
        f"- corr(theta, w): `{correlation([row['w'] for row in rows], [row['theta_deg'] for row in rows]):.6f}`",
        f"- corr(theta, b2): `{correlation([row['b2'] for row in rows], [row['theta_deg'] for row in rows]):.6f}`",
        f"- standardized theta regression weights: `{theta_reg}`",
        "",
        f"- corr(flip_lambda, frac): `{correlation([row['frac'] for row in rows], [row['flip_lambda'] for row in rows]):.6f}`",
        f"- corr(flip_lambda, w): `{correlation([row['w'] for row in rows], [row['flip_lambda'] for row in rows]):.6f}`",
        f"- corr(flip_lambda, b2): `{correlation([row['b2'] for row in rows], [row['flip_lambda'] for row in rows]):.6f}`",
        f"- standardized flip-lambda regression weights: `{flip_reg}`",
        "",
        "## Comparison To Empirical Low-c Corpora",
        "",
        f"- empirical low-c theta median: `{statistics.median(row['angle']['theta_deg'] for row in empirical_low_rows):.6f}`",
        f"- synthetic theta range: `[{min(row['theta_deg'] for row in rows):.6f}, {max(row['theta_deg'] for row in rows):.6f}]`",
        f"- empirical low-c flip-lambda median: `{statistics.median((1.05 if row['flip_lambda_exp_to_xpow'] is None else row['flip_lambda_exp_to_xpow']) for row in empirical_low_rows):.6f}`",
        f"- synthetic flip-lambda range: `[{min(row['flip_lambda'] for row in rows):.6f}, {max(row['flip_lambda'] for row in rows):.6f}]`",
        "",
        "Interpretation: if one parameter dominates both theta and flip-lambda, it is the main low-c manifold control knob.",
        "",
        "## Sample Rows",
        "",
        "| frac | w | b2 | delta_b | zm_c | theta | winner full | winner top100 | flip lambda |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: |",
    ]
    sample_rows = sorted(rows, key=lambda row: (row["frac"], row["w"], row["b2"]))[:: max(1, len(rows) // 24)]
    for row in sample_rows[:24]:
        lines.append(
            f"| {row['frac']:.2f} | {row['w']:.2f} | {row['b2']:.2f} | {row['delta_b']:.2f} | {row['zm_c']:.3f} | {row['theta_deg']:.3f} | `{row['winner_full']}` | `{row['winner_top100']}` | {row['flip_lambda']:.2f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    template = representative_lowc_template()
    phase_summary = load_phase_low_rows()
    rows = []
    for frac in FRAC_GRID:
        for w in W_GRID:
            for b2 in B2_GRID:
                curve = build_smooth_curve(template, frac, w, b2)
                rows.append(analyze_curve(curve, frac, w, b2))
    summary = {
        "template": template,
        "n_rows": len(rows),
        "full_winner_counts": dict(Counter(row["winner_full"] for row in rows)),
        "top100_winner_counts": dict(Counter(row["winner_top100"] for row in rows)),
        "corr_theta_frac": correlation([row["frac"] for row in rows], [row["theta_deg"] for row in rows]),
        "corr_theta_w": correlation([row["w"] for row in rows], [row["theta_deg"] for row in rows]),
        "corr_theta_b2": correlation([row["b2"] for row in rows], [row["theta_deg"] for row in rows]),
        "corr_flip_frac": correlation([row["frac"] for row in rows], [row["flip_lambda"] for row in rows]),
        "corr_flip_w": correlation([row["w"] for row in rows], [row["flip_lambda"] for row in rows]),
        "corr_flip_b2": correlation([row["b2"] for row in rows], [row["flip_lambda"] for row in rows]),
        "theta_regression": standardize_regression(rows, "theta_deg"),
        "flip_regression": standardize_regression(rows, "flip_lambda"),
        "theta_range": [float(min(row["theta_deg"] for row in rows)), float(max(row["theta_deg"] for row in rows))],
        "flip_range": [float(min(row["flip_lambda"] for row in rows)), float(max(row["flip_lambda"] for row in rows))],
    }
    payload = {"summary": summary, "rows": rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(template, rows, phase_summary), encoding="utf-8")
    plot_theta_scatter(rows, OUTDIR / "theta_flip_scatter.png")
    plot_top100_heatmaps(rows, OUTDIR / "top100_phase_heatmaps.png")


if __name__ == "__main__":
    main()
