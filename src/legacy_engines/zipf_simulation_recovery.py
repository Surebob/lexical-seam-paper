from __future__ import annotations

import csv
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
except ModuleNotFoundError:
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
ENRICHED_PATH = ROOT / "eml_zipf_enriched_search.py"
RERANKED_PATH = ROOT / "zipf_correct_model_reranked.py"
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_simulation_recovery"

N_REPLICATES = 10
HEAD_K = 200
TOP_K_GAP = 100


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_sim_recovery_common")
enriched = load_module(ENRICHED_PATH, "zipf_sim_recovery_enriched")
reranked = load_module(RERANKED_PATH, "zipf_sim_recovery_reranked")


def load_reranked_rows():
    return {row["slug"]: row for row in json.loads(RERANKED_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]}


def rank_dataset_from_freqs(freqs: np.ndarray) -> dict:
    freqs = np.asarray(freqs, dtype=np.float64)
    freqs = freqs[freqs > 0]
    freqs = np.sort(freqs)[::-1]
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    return {
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": np.log(ranks),
        "log_freq": np.log(freqs),
        "token_count": int(np.sum(freqs)),
        "unique_words": int(len(freqs)),
    }


def sample_rank_dataset(log_pred: np.ndarray, token_count: int, seed: int) -> dict:
    shifted = np.asarray(log_pred, dtype=np.float64) - float(np.max(log_pred))
    probs = np.exp(shifted)
    probs /= np.sum(probs)
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(token_count, probs)
    return rank_dataset_from_freqs(counts)


def fit_head_basis(x: np.ndarray, target: np.ndarray, head_k: int = HEAD_K) -> dict:
    k = min(head_k, len(x))
    u = 1.0 - x[:k]
    y = target[:k]
    X = np.column_stack([u, u * u, u * u * u])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coeffs
    sse = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else 1.0
    return {
        "linear_u": float(coeffs[0]),
        "quadratic_u2": float(coeffs[1]),
        "cubic_u3": float(coeffs[2]),
        "head_rmse": common.rmse(y, pred),
        "head_r2": float(r2),
        "head_k": int(k),
    }


def is_bregman(x: np.ndarray) -> np.ndarray:
    return (x - 1.0) - np.log(x)


def exp_bregman(x: np.ndarray) -> np.ndarray:
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def euclidean(x: np.ndarray) -> np.ndarray:
    return (1.0 - x) ** 2


def xpow_minus_sqrt(x: np.ndarray) -> np.ndarray:
    return np.power(x, x) - np.sqrt(x)


FORMULAS = {
    "is_bregman": ("sub[sub[x,1],log[x]]", is_bregman),
    "exp_bregman": ("eml[sub[x,1],eml[x,1]]", exp_bregman),
    "euclidean": ("mul[sub[1,x],sub[1,x]]", euclidean),
    "xpow_sqrt": ("sub[pow[x,x],sqrt[x]]", xpow_minus_sqrt),
}


def exact_step2_search(x: np.ndarray, target: np.ndarray) -> dict:
    current = enriched.initial_vocabulary(x, target)
    step1 = enriched.generate_candidates(current, target, 1)
    step1 = enriched.dedupe_candidates(step1)
    step1 = enriched.filter_candidates(step1, 1e-10)
    current = current + step1
    step2 = enriched.generate_candidates(current, target, 2)
    step2 = enriched.dedupe_candidates(step2)
    step2 = enriched.filter_candidates(step2, 1e-10)
    if not step2:
        raise RuntimeError("No valid step-2 candidates generated")
    best = min(step2, key=lambda item: (item["rmse"], item["expr"]))
    return best


def formula_gap_metrics(x: np.ndarray, target: np.ndarray, winner_values: np.ndarray, winner_rmse: float) -> dict:
    euclid_vals = euclidean(x)
    euclid_full = common.rmse(target, euclid_vals)
    k = min(TOP_K_GAP, len(x))
    winner_top = common.rmse(target[:k], winner_values[:k])
    euclid_top = common.rmse(target[:k], euclid_vals[:k])
    return {
        "winner_vs_euclidean_gap_full": float(euclid_full - winner_rmse),
        "winner_vs_euclidean_gap_top100": float(euclid_top - winner_top),
    }


def analyze_dataset(dataset: dict) -> dict:
    zm_fit = common.fit_zipf_mandelbrot(dataset["ranks"], dataset["log_freq"])
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    target = dataset["log_freq"] - zm_fit["prediction"]
    winner = exact_step2_search(x, target)
    basis = fit_head_basis(x, target)
    gaps = formula_gap_metrics(x, target, winner["values"], winner["rmse"])
    composite_rmse = common.rmse(dataset["log_freq"], zm_fit["prediction"] + winner["values"])
    simple_scores = {}
    for key, (_, fn) in FORMULAS.items():
        vals = fn(x)
        simple_scores[key] = {
            "full_rmse": common.rmse(target, vals),
            "top100_rmse": common.rmse(target[: min(TOP_K_GAP, len(x))], vals[: min(TOP_K_GAP, len(x))]),
        }
    return {
        "zm_rmse": float(zm_fit["rmse"]),
        "zm_c": float(zm_fit["c"]),
        "winner_expr": winner["expr"],
        "winner_rmse": float(winner["rmse"]),
        "step2_helpful": bool(composite_rmse < zm_fit["rmse"]),
        "step2_improvement": float(zm_fit["rmse"] - composite_rmse),
        "basis": basis,
        "gaps": gaps,
        "simple_scores": simple_scores,
    }


def reranked_prediction_for_dataset(dataset: dict, reranked_row: dict) -> np.ndarray:
    p = reranked_row["params"]
    vec = np.array([p["a1"], p["b1"], p["c1"], p["a2"], p["b2"], p["c2"], p["k"], p["w"]], dtype=np.float64)
    return reranked.reranked_prediction(dataset["ranks"], vec)


def empirical_summary_for_corpus(spec: dict) -> tuple[dict, dict]:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    analysis = analyze_dataset(dataset)
    return dataset, analysis


def correlation(xs: list[float], ys: list[float]) -> float:
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    if len(x) < 2 or float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def modal_winner(rows: list[dict]) -> tuple[str, float]:
    counts = Counter(row["winner_expr"] for row in rows)
    expr, count = counts.most_common(1)[0]
    return expr, float(count / len(rows))


def median_basis(rows: list[dict]) -> dict:
    return {
        "linear_u": float(np.median([row["basis"]["linear_u"] for row in rows])),
        "quadratic_u2": float(np.median([row["basis"]["quadratic_u2"] for row in rows])),
        "cubic_u3": float(np.median([row["basis"]["cubic_u3"] for row in rows])),
        "head_r2": float(np.median([row["basis"]["head_r2"] for row in rows])),
    }


def summarize_model_group(corpus_rows: list[dict], empirical_rows: list[dict], key: str) -> dict:
    replicate_rows = [rep for row in corpus_rows for rep in row[key]]
    exact_match_rate = float(np.mean([rep["winner_expr"] == row["empirical_winner"] for row in corpus_rows for rep in row[key]]))
    help_rate = float(np.mean([rep["step2_helpful"] for rep in replicate_rows]))
    winner_gap_full = float(np.median([rep["gaps"]["winner_vs_euclidean_gap_full"] for rep in replicate_rows]))
    winner_gap_top100 = float(np.median([rep["gaps"]["winner_vs_euclidean_gap_top100"] for rep in replicate_rows]))
    majority_match_count = int(sum(row[f"{key}_modal_winner"] == row["empirical_winner"] for row in corpus_rows))
    coeff_names = ["linear_u", "quadratic_u2", "cubic_u3"]
    basis_corrs = {
        name: correlation(
            [row["empirical_basis"][name] for row in corpus_rows],
            [row[f"{key}_median_basis"][name] for row in corpus_rows],
        )
        for name in coeff_names
    }
    return {
        "replicate_exact_winner_match_rate": exact_match_rate,
        "replicate_step2_help_rate": help_rate,
        "replicate_median_winner_vs_euclidean_gap_full": winner_gap_full,
        "replicate_median_winner_vs_euclidean_gap_top100": winner_gap_top100,
        "majority_winner_match_count": majority_match_count,
        "basis_correlations": basis_corrs,
    }


def plot_basis_recovery(corpus_rows: list[dict], outpath: Path):
    if plt is None:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    names = [("linear_u", "Linear"), ("quadratic_u2", "Quadratic"), ("cubic_u3", "Cubic")]
    for ax, (name, title) in zip(axes, names):
        empirical = [row["empirical_basis"][name] for row in corpus_rows]
        smooth = [row["smooth_median_basis"][name] for row in corpus_rows]
        control = [row["single_zm_median_basis"][name] for row in corpus_rows]
        ax.scatter(empirical, smooth, label="Smooth synthetic", color="#2563eb", alpha=0.85)
        ax.scatter(empirical, control, label="Single-ZM synthetic", color="#dc2626", alpha=0.75)
        lo = min(empirical + smooth + control)
        hi = max(empirical + smooth + control)
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="#666666", linewidth=1.0)
        ax.set_xlabel("Empirical")
        ax.set_ylabel("Synthetic median")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    axes[0].legend(loc="best")
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def write_csv(corpus_rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "empirical_winner",
        "empirical_zm_c",
        "empirical_linear",
        "empirical_quadratic",
        "empirical_cubic",
        "smooth_modal_winner",
        "smooth_modal_share",
        "smooth_exact_match_rate",
        "smooth_help_rate",
        "smooth_median_gap_full",
        "smooth_median_gap_top100",
        "single_zm_modal_winner",
        "single_zm_modal_share",
        "single_zm_exact_match_rate",
        "single_zm_help_rate",
        "single_zm_median_gap_full",
        "single_zm_median_gap_top100",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in corpus_rows:
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "empirical_winner": row["empirical_winner"],
                    "empirical_zm_c": row["empirical_zm_c"],
                    "empirical_linear": row["empirical_basis"]["linear_u"],
                    "empirical_quadratic": row["empirical_basis"]["quadratic_u2"],
                    "empirical_cubic": row["empirical_basis"]["cubic_u3"],
                    "smooth_modal_winner": row["smooth_modal_winner"],
                    "smooth_modal_share": row["smooth_modal_share"],
                    "smooth_exact_match_rate": row["smooth_exact_match_rate"],
                    "smooth_help_rate": row["smooth_help_rate"],
                    "smooth_median_gap_full": row["smooth_median_gap_full"],
                    "smooth_median_gap_top100": row["smooth_median_gap_top100"],
                    "single_zm_modal_winner": row["single_zm_modal_winner"],
                    "single_zm_modal_share": row["single_zm_modal_share"],
                    "single_zm_exact_match_rate": row["single_zm_exact_match_rate"],
                    "single_zm_help_rate": row["single_zm_help_rate"],
                    "single_zm_median_gap_full": row["single_zm_median_gap_full"],
                    "single_zm_median_gap_top100": row["single_zm_median_gap_top100"],
                }
            )


def build_report(summary: dict, corpus_rows: list[dict]) -> str:
    lines = [
        "# Simulation Recovery Test",
        "",
        "Synthetic corpora were sampled from fitted rank-probability curves, then re-analyzed with single ZM plus exact step-2 search.",
        "",
        f"- replicates per corpus per generator: `{summary['n_replicates']}`",
        f"- corpora analyzed: `{summary['n_corpora']}`",
        "",
        "## Overall Recovery",
        "",
        f"- empirical median winner-vs-Euclidean gap (full): `{summary['empirical_median_gap_full']:.12f}`",
        f"- empirical median winner-vs-Euclidean gap (top-100): `{summary['empirical_median_gap_top100']:.12f}`",
        "",
        "### Smooth-model synthetic corpora",
        "",
        f"- replicate exact winner match rate: `{summary['smooth']['replicate_exact_winner_match_rate']:.6f}`",
        f"- majority winner match count: `{summary['smooth']['majority_winner_match_count']}` / `{summary['n_corpora']}`",
        f"- replicate step-2 help rate: `{summary['smooth']['replicate_step2_help_rate']:.6f}`",
        f"- median winner-vs-Euclidean gap (full): `{summary['smooth']['replicate_median_winner_vs_euclidean_gap_full']:.12f}`",
        f"- median winner-vs-Euclidean gap (top-100): `{summary['smooth']['replicate_median_winner_vs_euclidean_gap_top100']:.12f}`",
        f"- basis corr linear/quadratic/cubic: `{summary['smooth']['basis_correlations']}`",
        "",
        "### Single-ZM synthetic controls",
        "",
        f"- replicate exact winner match rate: `{summary['single_zm']['replicate_exact_winner_match_rate']:.6f}`",
        f"- majority winner match count: `{summary['single_zm']['majority_winner_match_count']}` / `{summary['n_corpora']}`",
        f"- replicate step-2 help rate: `{summary['single_zm']['replicate_step2_help_rate']:.6f}`",
        f"- median winner-vs-Euclidean gap (full): `{summary['single_zm']['replicate_median_winner_vs_euclidean_gap_full']:.12f}`",
        f"- median winner-vs-Euclidean gap (top-100): `{summary['single_zm']['replicate_median_winner_vs_euclidean_gap_top100']:.12f}`",
        f"- basis corr linear/quadratic/cubic: `{summary['single_zm']['basis_correlations']}`",
        "",
        "## Per-Corpus Winner Recovery",
        "",
        "| corpus | empirical winner | smooth modal (share) | smooth exact-match rate | smooth help rate | ZM modal (share) | ZM exact-match rate | ZM help rate |",
        "| --- | --- | --- | ---: | ---: | --- | ---: | ---: |",
    ]
    for row in corpus_rows:
        lines.append(
            f"| {row['name']} | `{row['empirical_winner']}` | `{row['smooth_modal_winner']}` ({row['smooth_modal_share']:.2f}) | {row['smooth_exact_match_rate']:.2f} | {row['smooth_help_rate']:.2f} | `{row['single_zm_modal_winner']}` ({row['single_zm_modal_share']:.2f}) | {row['single_zm_exact_match_rate']:.2f} | {row['single_zm_help_rate']:.2f} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    reranked_rows = load_reranked_rows()
    corpus_rows = []
    all_empirical = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA):
        dataset, empirical = empirical_summary_for_corpus(spec)
        reranked_row = reranked_rows[spec["slug"]]
        smooth_log_pred = reranked_prediction_for_dataset(dataset, reranked_row)
        zm_log_pred = common.zm_prediction(common.load_enriched_summary(spec), dataset["ranks"])
        smooth_reps = []
        zm_reps = []
        for rep in range(N_REPLICATES):
            smooth_ds = sample_rank_dataset(smooth_log_pred, dataset["token_count"], seed=100000 + corpus_index * 100 + rep)
            zm_ds = sample_rank_dataset(zm_log_pred, dataset["token_count"], seed=200000 + corpus_index * 100 + rep)
            smooth_reps.append(analyze_dataset(smooth_ds))
            zm_reps.append(analyze_dataset(zm_ds))
        smooth_modal_winner, smooth_modal_share = modal_winner(smooth_reps)
        zm_modal_winner, zm_modal_share = modal_winner(zm_reps)
        row = {
            "slug": spec["slug"],
            "name": spec["name"],
            "empirical_winner": empirical["winner_expr"],
            "empirical_zm_c": empirical["zm_c"],
            "empirical_basis": empirical["basis"],
            "smooth": smooth_reps,
            "single_zm": zm_reps,
            "smooth_modal_winner": smooth_modal_winner,
            "smooth_modal_share": smooth_modal_share,
            "smooth_median_basis": median_basis(smooth_reps),
            "single_zm_modal_winner": zm_modal_winner,
            "single_zm_modal_share": zm_modal_share,
            "single_zm_median_basis": median_basis(zm_reps),
            "smooth_exact_match_rate": float(np.mean([rep["winner_expr"] == empirical["winner_expr"] for rep in smooth_reps])),
            "smooth_help_rate": float(np.mean([rep["step2_helpful"] for rep in smooth_reps])),
            "smooth_median_gap_full": float(np.median([rep["gaps"]["winner_vs_euclidean_gap_full"] for rep in smooth_reps])),
            "smooth_median_gap_top100": float(np.median([rep["gaps"]["winner_vs_euclidean_gap_top100"] for rep in smooth_reps])),
            "single_zm_exact_match_rate": float(np.mean([rep["winner_expr"] == empirical["winner_expr"] for rep in zm_reps])),
            "single_zm_help_rate": float(np.mean([rep["step2_helpful"] for rep in zm_reps])),
            "single_zm_median_gap_full": float(np.median([rep["gaps"]["winner_vs_euclidean_gap_full"] for rep in zm_reps])),
            "single_zm_median_gap_top100": float(np.median([rep["gaps"]["winner_vs_euclidean_gap_top100"] for rep in zm_reps])),
        }
        corpus_rows.append(row)
        all_empirical.append(empirical)

    summary = {
        "n_replicates": N_REPLICATES,
        "n_corpora": len(corpus_rows),
        "empirical_median_gap_full": float(np.median([row["gaps"]["winner_vs_euclidean_gap_full"] for row in all_empirical])),
        "empirical_median_gap_top100": float(np.median([row["gaps"]["winner_vs_euclidean_gap_top100"] for row in all_empirical])),
        "smooth": summarize_model_group(corpus_rows, all_empirical, "smooth"),
        "single_zm": summarize_model_group(corpus_rows, all_empirical, "single_zm"),
    }
    payload = {"summary": summary, "rows": corpus_rows}
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary, corpus_rows), encoding="utf-8")
    write_csv(corpus_rows, OUTDIR / "simulation_recovery_table.csv")
    plot_basis_recovery(corpus_rows, OUTDIR / "basis_recovery_scatter.png")


if __name__ == "__main__":
    main()
