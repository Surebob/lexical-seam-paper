import importlib.util
import json
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    matplotlib = None
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
CORRECT_MODEL_PATH = ROOT / "zipf_correct_model.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
OUTDIR = ROOT / "results" / "zipf_correct_model_reranked_all"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_correct_model_reranked_all_common")
correct_model = load_module(CORRECT_MODEL_PATH, "zipf_correct_model_reranked_all_base")
reranked = load_module(RERANKED_PATH, "zipf_correct_model_reranked_impl")


def summarize_row(spec: dict, n_starts: int = 100):
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    summary = common.load_enriched_summary(spec)
    piecewise = correct_model.fit_piecewise_zm(dataset, 500)
    best, tried = reranked.fit_reranked_model(dataset, n_starts=n_starts, max_nfev=12000)
    step2 = common.get_step2_candidate(summary)
    params = reranked.summarize_params(best["params"], dataset["unique_words"])
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "single_zm_rmse": float(summary["zm_baseline"]["rmse_full"]),
        "piecewise_rmse": float(piecewise["rmse"]),
        "step2_rmse": float(step2["rmse"]),
        "step2_expr": step2["expr"],
        "reranked_rmse": float(best["rmse"]),
        "beats_single_zm": bool(best["rmse"] < summary["zm_baseline"]["rmse_full"]),
        "beats_piecewise": bool(best["rmse"] < piecewise["rmse"]),
        "beats_step2": bool(best["rmse"] < step2["rmse"]),
        "best_start_index": int(best["start_index"]),
        "best_nfev": int(best["nfev"]),
        "params": params,
        "k": float(params["k"]),
        "w": float(params["w"]),
        "transition_fraction": float(params["transition_fraction"]),
        "tried_count": len(tried),
    }


def build_summary(rows: list[dict]):
    k_values = [row["k"] for row in rows]
    w_values = [row["w"] for row in rows]
    frac_values = [row["transition_fraction"] for row in rows]
    return {
        "rows": rows,
        "parameter_consistency": {
            "k_mean": float(np.mean(k_values)),
            "k_std": float(np.std(k_values)),
            "k_min": float(np.min(k_values)),
            "k_max": float(np.max(k_values)),
            "w_mean": float(np.mean(w_values)),
            "w_std": float(np.std(w_values)),
            "w_min": float(np.min(w_values)),
            "w_max": float(np.max(w_values)),
            "transition_fraction_mean": float(np.mean(frac_values)),
            "transition_fraction_std": float(np.std(frac_values)),
            "transition_fraction_min": float(np.min(frac_values)),
            "transition_fraction_max": float(np.max(frac_values)),
        },
        "counts": {
            "beats_single_zm": int(sum(row["beats_single_zm"] for row in rows)),
            "beats_piecewise": int(sum(row["beats_piecewise"] for row in rows)),
            "beats_step2": int(sum(row["beats_step2"] for row in rows)),
            "n_rows": len(rows),
        },
    }


def build_report(summary: dict):
    counts = summary["counts"]
    stats = summary["parameter_consistency"]
    rows = sorted(summary["rows"], key=lambda row: row["k"])
    lines = [
        "# Reranked Smooth Model on All 25 Corpora",
        "",
        f"- corpora beating single ZM: `{counts['beats_single_zm']}` / `{counts['n_rows']}`",
        f"- corpora beating hard two-component ZM (`K=500`): `{counts['beats_piecewise']}` / `{counts['n_rows']}`",
        f"- corpora beating step-2 correction: `{counts['beats_step2']}` / `{counts['n_rows']}`",
        f"- k mean/std: `{stats['k_mean']:.6f}` / `{stats['k_std']:.6f}`",
        f"- k min/max: `{stats['k_min']:.6f}` / `{stats['k_max']:.6f}`",
        f"- w mean/std: `{stats['w_mean']:.6f}` / `{stats['w_std']:.6f}`",
        f"- transition fraction mean/std: `{stats['transition_fraction_mean']:.6f}` / `{stats['transition_fraction_std']:.6f}`",
        "",
        "| corpus | single ZM | piecewise | step-2 | reranked | k | w | frac | step-2 winner |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['single_zm_rmse']:.12f} | {row['piecewise_rmse']:.12f} | {row['step2_rmse']:.12f} | {row['reranked_rmse']:.12f} | {row['k']:.3f} | {row['w']:.3f} | {row['transition_fraction']:.3f} | {row['step2_expr']} |"
        )
    return "\n".join(lines) + "\n"


def plot_k_vs_vocab(rows: list[dict], outpath: Path):
    if plt is None:
        return
    rows = sorted(rows, key=lambda row: row["vocab_size"])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = [row["vocab_size"] for row in rows]
    y = [row["k"] for row in rows]
    ax.scatter(x, y, color="#1565c0", s=28)
    for row in rows:
        ax.annotate(row["slug"], (row["vocab_size"], row["k"]), fontsize=7, alpha=0.8, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("vocabulary size")
    ax.set_ylabel("fitted transition k")
    ax.set_title("Reranked Smooth Model Transition Center vs Vocabulary Size")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_fraction_hist(rows: list[dict], outpath: Path):
    if plt is None:
        return
    fractions = [row["transition_fraction"] for row in rows]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.hist(fractions, bins=10, color="#c62828", edgecolor="white")
    ax.set_xlabel("transition fraction")
    ax.set_ylabel("count")
    ax.set_title("Reranked Smooth Model Transition Fraction Distribution")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in common.SEARCHED_CORPORA:
        rows.append(summarize_row(spec, n_starts=100))
    summary = build_summary(rows)
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    plot_k_vs_vocab(rows, OUTDIR / "k_vs_vocab.png")
    plot_fraction_hist(rows, OUTDIR / "transition_fraction_hist.png")


if __name__ == "__main__":
    main()
