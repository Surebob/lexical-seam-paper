from __future__ import annotations

import importlib.util
import json
import statistics
from collections import Counter
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
SIM_RECOVERY_SUMMARY_PATH = ROOT / "results" / "zipf_simulation_recovery" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_lowc_manifold_analysis"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_lowc_manifold_common")


def exp_bregman(x: np.ndarray) -> np.ndarray:
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def xpow_minus_sqrt(x: np.ndarray) -> np.ndarray:
    return np.power(x, x) - np.sqrt(x)


def is_bregman(x: np.ndarray) -> np.ndarray:
    return (x - 1.0) - np.log(x)


def rmse(y: np.ndarray, p: np.ndarray) -> float:
    return common.rmse(y, p)


def topk_rmse(y: np.ndarray, p: np.ndarray, k: int) -> float:
    k = min(k, len(y))
    return rmse(y[:k], p[:k])


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    a0 = a - np.mean(a)
    b0 = b - np.mean(b)
    denom = float(np.linalg.norm(a0) * np.linalg.norm(b0))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(a0, b0) / denom)


def span_r2(y: np.ndarray, cols: list[np.ndarray]) -> tuple[float, list[float]]:
    X = np.column_stack(cols)
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coeffs
    sse = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else 1.0
    return float(r2), [float(v) for v in coeffs]


def winner_name(scores: dict[str, float]) -> str:
    return min(scores, key=lambda key: (scores[key], key))


def analyze_corpus(spec: dict) -> dict:
    summary = common.load_enriched_summary(spec)
    if common.get_step2_expr(summary) != common.STEPL2_LOWC_EXPR:
        raise ValueError(f"{spec['slug']} is not a low-c exponential-winner corpus")
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    target = dataset["log_freq"] - common.zm_prediction(summary, dataset["ranks"])
    f_exp = exp_bregman(x)
    f_xpow = xpow_minus_sqrt(x)
    f_is = is_bregman(x)

    full_scores = {
        "exp": rmse(target, f_exp),
        "xpow": rmse(target, f_xpow),
        "is": rmse(target, f_is),
    }
    top50_scores = {
        "exp": topk_rmse(target, f_exp, 50),
        "xpow": topk_rmse(target, f_xpow, 50),
        "is": topk_rmse(target, f_is, 50),
    }
    top100_scores = {
        "exp": topk_rmse(target, f_exp, 100),
        "xpow": topk_rmse(target, f_xpow, 100),
        "is": topk_rmse(target, f_is, 100),
    }
    top200_scores = {
        "exp": topk_rmse(target, f_exp, 200),
        "xpow": topk_rmse(target, f_xpow, 200),
        "is": topk_rmse(target, f_is, 200),
    }
    head_k = min(200, len(x))
    cos_exp_xpow = cosine_similarity(f_exp[:head_k], f_xpow[:head_k])
    r2_expxpow, beta_expxpow = span_r2(target[:head_k], [f_exp[:head_k], f_xpow[:head_k]])
    r2_expis, beta_expis = span_r2(target[:head_k], [f_exp[:head_k], f_is[:head_k]])
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "c": float(summary["zm_baseline"]["c"]),
        "full_winner": winner_name(full_scores),
        "top50_winner": winner_name(top50_scores),
        "top100_winner": winner_name(top100_scores),
        "top200_winner": winner_name(top200_scores),
        "full_scores": full_scores,
        "top50_scores": top50_scores,
        "top100_scores": top100_scores,
        "top200_scores": top200_scores,
        "cos_exp_xpow_top200": cos_exp_xpow,
        "r2_span_expxpow_top200": r2_expxpow,
        "beta_span_expxpow_top200": beta_expxpow,
        "r2_span_expis_top200": r2_expis,
        "beta_span_expis_top200": beta_expis,
        "xpow_minus_exp_full": float(full_scores["xpow"] - full_scores["exp"]),
        "xpow_minus_exp_top100": float(top100_scores["xpow"] - top100_scores["exp"]),
    }


def plot_rmse_flip(rows: list[dict], outpath: Path):
    if plt is None:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    full_exp = [row["full_scores"]["exp"] for row in rows]
    full_xpow = [row["full_scores"]["xpow"] for row in rows]
    top_exp = [row["top100_scores"]["exp"] for row in rows]
    top_xpow = [row["top100_scores"]["xpow"] for row in rows]
    for ax, ex, xp, title in [
        (axes[0], full_exp, full_xpow, "Full RMSE"),
        (axes[1], top_exp, top_xpow, "Top-100 RMSE"),
    ]:
        lo = min(ex + xp)
        hi = max(ex + xp)
        ax.scatter(ex, xp, color="#2563eb", alpha=0.85)
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="#666666", linewidth=1.0)
        ax.set_xlabel("exp(x-1)-x RMSE")
        ax.set_ylabel("x^x - sqrt(x) RMSE")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def build_report(rows: list[dict], sim_summary: dict) -> str:
    full_counts = Counter(row["full_winner"] for row in rows)
    top50_counts = Counter(row["top50_winner"] for row in rows)
    top100_counts = Counter(row["top100_winner"] for row in rows)
    top200_counts = Counter(row["top200_winner"] for row in rows)
    lines = [
        "# Low-c Manifold Analysis",
        "",
        "These are the 14 English corpora whose empirical step-2 winner is the exponential Bregman generator `exp(x-1)-x`.",
        "",
        f"- median cosine between `exp(x-1)-x` and `x^x-sqrt(x)` over top-200 head coordinates: `{statistics.median(row['cos_exp_xpow_top200'] for row in rows):.12f}`",
        f"- median head-200 R^2 of span{{exp, xpow}}: `{statistics.median(row['r2_span_expxpow_top200'] for row in rows):.12f}`",
        f"- median head-200 R^2 of span{{exp, IS}}: `{statistics.median(row['r2_span_expis_top200'] for row in rows):.12f}`",
        f"- median `(xpow - exp)` RMSE on full curve: `{statistics.median(row['xpow_minus_exp_full'] for row in rows):.12f}`",
        f"- median `(xpow - exp)` RMSE on top-100: `{statistics.median(row['xpow_minus_exp_top100'] for row in rows):.12f}`",
        "",
        "## Winner Counts On The Same 14 Corpora",
        "",
        f"- full RMSE winners: `{dict(full_counts)}`",
        f"- top-50 RMSE winners: `{dict(top50_counts)}`",
        f"- top-100 RMSE winners: `{dict(top100_counts)}`",
        f"- top-200 RMSE winners: `{dict(top200_counts)}`",
        "",
        "## Link To Simulation Recovery",
        "",
        f"- smooth synthetic exact winner match rate over all corpora: `{sim_summary['summary']['smooth']['replicate_exact_winner_match_rate']:.6f}`",
        f"- smooth synthetic exact winner match rate on the empirical low-c family is only `{sum(row['smooth_exact_match_rate'] for row in sim_summary['rows'] if row['empirical_winner']==common.STEPL2_LOWC_EXPR)/len(rows):.6f}`",
        f"- smooth synthetic modal winners on the empirical low-c family: `{dict(Counter(row['smooth_modal_winner'] for row in sim_summary['rows'] if row['empirical_winner']==common.STEPL2_LOWC_EXPR))}`",
        "",
        "This says the low-c side is not a clean exp-vs-not-exp phase split. It is a near-degenerate head manifold where full-curve scoring picks `exp`, but head-focused scoring usually picks `x^x-sqrt(x)`.",
        "",
        "## Per-Corpus Table",
        "",
        "| corpus | c | cosine(exp,xpow) | full winner | top-100 winner | xpow-exp full | xpow-exp top100 | R^2 span{exp,xpow} | R^2 span{exp,IS} |",
        "| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['c']:.3f} | {row['cos_exp_xpow_top200']:.6f} | `{row['full_winner']}` | `{row['top100_winner']}` | {row['xpow_minus_exp_full']:.6f} | {row['xpow_minus_exp_top100']:.6f} | {row['r2_span_expxpow_top200']:.6f} | {row['r2_span_expis_top200']:.6f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = [analyze_corpus(spec) for spec in common.SEARCHED_CORPORA if common.get_step2_expr(common.load_enriched_summary(spec)) == common.STEPL2_LOWC_EXPR]
    sim_summary = json.loads(SIM_RECOVERY_SUMMARY_PATH.read_text(encoding="utf-8"))
    summary = {
        "n_lowc_corpora": len(rows),
        "median_cos_exp_xpow_top200": float(statistics.median(row["cos_exp_xpow_top200"] for row in rows)),
        "median_r2_span_expxpow_top200": float(statistics.median(row["r2_span_expxpow_top200"] for row in rows)),
        "median_r2_span_expis_top200": float(statistics.median(row["r2_span_expis_top200"] for row in rows)),
        "median_xpow_minus_exp_full": float(statistics.median(row["xpow_minus_exp_full"] for row in rows)),
        "median_xpow_minus_exp_top100": float(statistics.median(row["xpow_minus_exp_top100"] for row in rows)),
        "full_winner_counts": dict(Counter(row["full_winner"] for row in rows)),
        "top50_winner_counts": dict(Counter(row["top50_winner"] for row in rows)),
        "top100_winner_counts": dict(Counter(row["top100_winner"] for row in rows)),
        "top200_winner_counts": dict(Counter(row["top200_winner"] for row in rows)),
    }
    payload = {"summary": summary, "rows": rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(rows, sim_summary), encoding="utf-8")
    plot_rmse_flip(rows, OUTDIR / "lowc_full_vs_head_flip.png")


if __name__ == "__main__":
    main()
