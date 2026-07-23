from __future__ import annotations

import importlib.util
import json
import math
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
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
SIM_RECOVERY_SUMMARY_PATH = ROOT / "results" / "zipf_simulation_recovery" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_phase_coordinate"

LAMBDA_GRID = np.linspace(0.0, 1.0, 21)
HEAD_TOP_K = 100
ANGLE_HEAD_K = 200


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_phase_coordinate_common")


def load_reranked_rows():
    rows = json.loads(RERANKED_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def exp_bregman(x: np.ndarray) -> np.ndarray:
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def xpow_minus_sqrt(x: np.ndarray) -> np.ndarray:
    return np.power(x, x) - np.sqrt(x)


def is_bregman(x: np.ndarray) -> np.ndarray:
    return (x - 1.0) - np.log(x)


FORMULAS = {
    "is": ("sub[sub[x,1],log[x]]", is_bregman),
    "exp": ("eml[sub[x,1],eml[x,1]]", exp_bregman),
    "xpow": ("sub[pow[x,x],sqrt[x]]", xpow_minus_sqrt),
}

WINNER_CODE = {"is": 0, "exp": 1, "xpow": 2}
WINNER_COLOR = ListedColormap(["#2563eb", "#059669", "#dc2626"]) if plt is not None else None


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
    theta_deg = float(np.degrees(np.arctan2(beta[1], beta[0])))
    amplitude = float(np.linalg.norm(beta))
    return {
        "theta_deg": theta_deg,
        "beta_exp": float(beta[0]),
        "beta_xpow": float(beta[1]),
        "amplitude": amplitude,
        "r2": float(r2),
    }


def winner_sweep(target: np.ndarray, curves: dict[str, np.ndarray], top_k: int = HEAD_TOP_K) -> dict:
    n = len(target)
    top_k = min(top_k, n)
    head = np.zeros(n, dtype=np.float64)
    head[:top_k] = 1.0
    rows = []
    flip_lambda = None
    exp_curve = curves["exp"]
    xpow_curve = curves["xpow"]
    for lam in LAMBDA_GRID:
        weights = (1.0 - lam) * np.ones(n, dtype=np.float64) + lam * head
        scores = {name: weighted_rmse(target, curve, weights) for name, curve in curves.items()}
        winner = min(scores, key=lambda key: (scores[key], key))
        if flip_lambda is None and scores["xpow"] < scores["exp"]:
            flip_lambda = float(lam)
        rows.append({"lambda": float(lam), "winner": winner, "scores": scores})
    return {"rows": rows, "flip_lambda_exp_to_xpow": flip_lambda}


def analyze_corpus(spec: dict, reranked_row: dict) -> dict:
    summary = common.load_enriched_summary(spec)
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    target = dataset["log_freq"] - common.zm_prediction(summary, dataset["ranks"])
    curves = {name: fn(x) for name, (_, fn) in FORMULAS.items()}
    sweep = winner_sweep(target, curves)
    angle = angle_fit(target, curves["exp"], curves["xpow"])
    empirical_expr = common.get_step2_expr(summary)
    empirical_winner = next(name for name, (expr, _) in FORMULAS.items() if expr == empirical_expr)
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "c": float(summary["zm_baseline"]["c"]),
        "transition_fraction": float(reranked_row["transition_fraction"]),
        "empirical_winner": empirical_winner,
        "angle": angle,
        "sweep": sweep,
        "winner_full": sweep["rows"][0]["winner"],
        "winner_top100": sweep["rows"][-1]["winner"],
        "flip_lambda_exp_to_xpow": sweep["flip_lambda_exp_to_xpow"],
    }


def correlation(xs: list[float], ys: list[float]) -> float:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    if len(x) < 2 or float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def plot_lowc_theta(rows: list[dict], outpath: Path):
    if plt is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    xs = [row["c"] for row in rows]
    ys = [row["angle"]["theta_deg"] for row in rows]
    flips = [row["flip_lambda_exp_to_xpow"] if row["flip_lambda_exp_to_xpow"] is not None else 1.05 for row in rows]
    axes[0].scatter(xs, ys, color="#2563eb", alpha=0.9)
    for row in rows:
        axes[0].annotate(row["slug"], (row["c"], row["angle"]["theta_deg"]), fontsize=7, alpha=0.7)
    axes[0].set_xlabel("ZM c")
    axes[0].set_ylabel("theta (deg)")
    axes[0].set_title("Low-c angle vs c")
    axes[0].grid(True, alpha=0.25)

    axes[1].scatter(flips, ys, color="#dc2626", alpha=0.9)
    for row, flip in zip(rows, flips):
        axes[1].annotate(row["slug"], (flip, row["angle"]["theta_deg"]), fontsize=7, alpha=0.7)
    axes[1].set_xlabel("lambda flip (exp -> xpow)")
    axes[1].set_ylabel("theta (deg)")
    axes[1].set_title("Low-c angle vs flip lambda")
    axes[1].grid(True, alpha=0.25)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_winner_heatmap(rows: list[dict], outpath: Path):
    if plt is None:
        return
    ordered = sorted(rows, key=lambda row: (row["empirical_winner"], row["c"]))
    matrix = np.array([[WINNER_CODE[cell["winner"]] for cell in row["sweep"]["rows"]] for row in ordered], dtype=np.int64)
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    im = ax.imshow(matrix, aspect="auto", cmap=WINNER_COLOR, interpolation="nearest")
    ax.set_xticks(range(len(LAMBDA_GRID)))
    ax.set_xticklabels([f"{lam:.2f}" for lam in LAMBDA_GRID], rotation=90)
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels([row["slug"] for row in ordered])
    ax.set_xlabel("Head-weight lambda (0 = full RMSE, 1 = top-100 RMSE)")
    ax.set_title("Winner sweep across lambda")
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(color="#2563eb", label="IS"),
        Patch(color="#059669", label="exp"),
        Patch(color="#dc2626", label="xpow"),
    ]
    ax.legend(handles=legend_handles, loc="upper right")
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def build_report(rows: list[dict], low_rows: list[dict], sim_summary: dict) -> str:
    low_counter_full = Counter(row["winner_full"] for row in low_rows)
    low_counter_top = Counter(row["winner_top100"] for row in low_rows)
    high_rows = [row for row in rows if row["empirical_winner"] == "is"]
    high_counter_full = Counter(row["winner_full"] for row in high_rows)
    high_counter_top = Counter(row["winner_top100"] for row in high_rows)
    sim_low_rows = [row for row in sim_summary["rows"] if row["empirical_winner"] == common.STEPL2_LOWC_EXPR]
    expr_to_name = {expr: name for name, (expr, _) in FORMULAS.items()}
    expr_to_name["sub[sqrt[x],pow[x,x]]"] = "xpow"
    smooth_match_full = float(
        np.mean([expr_to_name.get(row["smooth_modal_winner"], row["smooth_modal_winner"]) == "exp" for row in sim_low_rows])
    )
    smooth_match_top = float(
        np.mean(
            [
                expr_to_name.get(row["smooth_modal_winner"], row["smooth_modal_winner"])
                == next(low["winner_top100"] for low in low_rows if low["slug"] == row["slug"])
                for row in sim_low_rows
            ]
        )
    )
    zm_match_top = float(
        np.mean(
            [
                expr_to_name.get(row["single_zm_modal_winner"], row["single_zm_modal_winner"])
                == next(low["winner_top100"] for low in low_rows if low["slug"] == row["slug"])
                for row in sim_low_rows
            ]
        )
    )
    lines = [
        "# Phase Coordinate Analysis",
        "",
        "A weighted-loss sweep interpolates between full-curve RMSE (`lambda = 0`) and top-100-only RMSE (`lambda = 1`).",
        "",
        "## High-vs-Low Winner Stability",
        "",
        f"- high-c/IS block full winners: `{dict(high_counter_full)}`",
        f"- high-c/IS block top-100 winners: `{dict(high_counter_top)}`",
        f"- low-c/exp block full winners: `{dict(low_counter_full)}`",
        f"- low-c/exp block top-100 winners: `{dict(low_counter_top)}`",
        "",
        "## Low-c Phase Coordinate",
        "",
        f"- median theta over low-c corpora: `{np.median([row['angle']['theta_deg'] for row in low_rows]):.12f}` deg",
        f"- median head-200 fit R^2 in span{{exp,xpow}}: `{np.median([row['angle']['r2'] for row in low_rows]):.12f}`",
        f"- corr(theta, c): `{correlation([row['c'] for row in low_rows], [row['angle']['theta_deg'] for row in low_rows]):.6f}`",
        f"- corr(theta, transition_fraction): `{correlation([row['transition_fraction'] for row in low_rows], [row['angle']['theta_deg'] for row in low_rows]):.6f}`",
        f"- corr(theta, lambda_flip): `{correlation([row['flip_lambda_exp_to_xpow'] if row['flip_lambda_exp_to_xpow'] is not None else 1.05 for row in low_rows], [row['angle']['theta_deg'] for row in low_rows]):.6f}`",
        "",
        "## Link To Smooth Synthetic Recovery",
        "",
        f"- on the low-c family, smooth synthetic modal winner matches empirical full-RMSE winner only `{smooth_match_full:.6f}` of the time",
        f"- but it matches the empirical top-100 winner `{smooth_match_top:.6f}` of the time",
        f"- single-ZM control matches the empirical top-100 winner `{zm_match_top:.6f}` of the time",
        "",
        "This is the key asymmetry: the high-c side is a stable IS phase, while the low-c side is a rotating manifold where the named winner depends strongly on the scoring functional.",
        "",
        "## Per-Corpus Low-c Table",
        "",
        "| corpus | c | transition frac | theta (deg) | span R^2 | full winner | top-100 winner | lambda flip exp->xpow |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- | ---: |",
    ]
    for row in low_rows:
        flip = row["flip_lambda_exp_to_xpow"]
        flip_text = "never" if flip is None else f"{flip:.2f}"
        lines.append(
            f"| {row['name']} | {row['c']:.3f} | {row['transition_fraction']:.6f} | {row['angle']['theta_deg']:.3f} | {row['angle']['r2']:.6f} | `{row['winner_full']}` | `{row['winner_top100']}` | {flip_text} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_rows = load_reranked_rows()
    rows = [analyze_corpus(spec, reranked_rows[spec["slug"]]) for spec in common.SEARCHED_CORPORA]
    low_rows = [row for row in rows if row["empirical_winner"] == "exp"]
    sim_summary = json.loads(SIM_RECOVERY_SUMMARY_PATH.read_text(encoding="utf-8"))
    summary = {
        "n_corpora": len(rows),
        "n_lowc_corpora": len(low_rows),
        "high_full_counts": dict(Counter(row["winner_full"] for row in rows if row["empirical_winner"] == "is")),
        "high_top100_counts": dict(Counter(row["winner_top100"] for row in rows if row["empirical_winner"] == "is")),
        "low_full_counts": dict(Counter(row["winner_full"] for row in low_rows)),
        "low_top100_counts": dict(Counter(row["winner_top100"] for row in low_rows)),
        "corr_theta_c": correlation([row["c"] for row in low_rows], [row["angle"]["theta_deg"] for row in low_rows]),
        "corr_theta_transition_fraction": correlation(
            [row["transition_fraction"] for row in low_rows],
            [row["angle"]["theta_deg"] for row in low_rows],
        ),
        "corr_theta_lambda_flip": correlation(
            [row["flip_lambda_exp_to_xpow"] if row["flip_lambda_exp_to_xpow"] is not None else 1.05 for row in low_rows],
            [row["angle"]["theta_deg"] for row in low_rows],
        ),
    }
    payload = {"summary": summary, "rows": rows, "low_rows": low_rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(rows, low_rows, sim_summary), encoding="utf-8")
    plot_lowc_theta(low_rows, OUTDIR / "lowc_theta_phase.png")
    plot_winner_heatmap(rows, OUTDIR / "winner_sweep_heatmap.png")


if __name__ == "__main__":
    main()
