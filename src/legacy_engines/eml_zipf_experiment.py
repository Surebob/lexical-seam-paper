import argparse
import importlib.util
import json
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np
try:
    import torch
except ModuleNotFoundError:
    torch = None

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    matplotlib = None
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
PORT_PATH = ROOT / "eml_tree_port.py"
DEFAULT_CORPUS_URL = "https://www.gutenberg.org/cache/epub/100/pg100.txt"
DEFAULT_CORPUS_PATH = ROOT / "data" / "zipf" / "pg100.txt"
DEFAULT_OUTDIR = ROOT / "results" / "zipf_shakespeare"
TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if torch is not None:
    port = load_module(PORT_PATH, "eml_tree_port_module_zipf")
else:
    port = None


def sanitize_jsonable(obj):
    if isinstance(obj, dict):
        return {k: sanitize_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj


def parse_args():
    p = argparse.ArgumentParser(description="EML Zipf residual experiment")
    p.add_argument("--corpus-url", type=str, default=DEFAULT_CORPUS_URL)
    p.add_argument("--corpus-path", type=Path, default=DEFAULT_CORPUS_PATH)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--depths", type=int, nargs="*", default=[2, 3, 4])
    p.add_argument("--modes", type=str, nargs="*", default=["direct", "residual"], choices=["direct", "residual"])
    p.add_argument("--sample-points", type=int, default=200)
    p.add_argument("--seed0", type=int, default=137)
    p.add_argument("--seeds", type=int, default=16)
    p.add_argument("--search-iters", type=int, default=6000)
    p.add_argument("--hardening-iters", type=int, default=2000)
    p.add_argument("--tau-search", type=float, default=2.5)
    p.add_argument("--tau-hard", type=float, default=0.01)
    p.add_argument("--eval-every", type=int, default=200)
    p.add_argument("--tail-eval-every", type=int, default=50)
    p.add_argument("--x-low", type=float, default=0.05)
    p.add_argument("--x-high", type=float, default=1.0)
    p.add_argument("--skip-plot", action="store_true")
    return p.parse_args()


def ensure_corpus(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return dest
    with urllib.request.urlopen(url) as response:
        data = response.read()
    dest.write_bytes(data)
    return dest


def strip_gutenberg_boilerplate(text: str):
    start_markers = [
        "*** START OF THE PROJECT GUTENBERG EBOOK",
        "*** START OF THIS PROJECT GUTENBERG EBOOK",
    ]
    end_markers = [
        "*** END OF THE PROJECT GUTENBERG EBOOK",
        "*** END OF THIS PROJECT GUTENBERG EBOOK",
    ]
    start = 0
    end = len(text)
    for marker in start_markers:
        idx = text.find(marker)
        if idx != -1:
            line_end = text.find("\n", idx)
            start = line_end + 1 if line_end != -1 else idx
            break
    for marker in end_markers:
        idx = text.find(marker)
        if idx != -1:
            end = idx
            break
    return text[start:end]


def tokenize_text(text: str):
    return TOKEN_RE.findall(text.lower())


def build_zipf_dataset(corpus_path: Path):
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = strip_gutenberg_boilerplate(raw_text)
    tokens = tokenize_text(clean_text)
    counts = Counter(tokens)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)
    return {
        "tokens": tokens,
        "counts": counts,
        "ranked": ranked,
        "ranks": ranks,
        "freqs": freqs,
        "log_rank": log_rank,
        "log_freq": log_freq,
        "unique_words": len(freqs),
        "token_count": len(tokens),
        "top_words": ranked[:20],
    }


def fit_linear_zipf(log_rank: np.ndarray, log_freq: np.ndarray):
    slope, intercept = np.polyfit(log_rank, log_freq, 1)
    pred = intercept + slope * log_rank
    mse = float(np.mean((pred - log_freq) ** 2))
    return {
        "intercept": float(intercept),
        "slope": float(slope),
        "mse": mse,
        "rmse": float(math.sqrt(max(mse, 0.0))),
        "prediction": pred,
    }


def _solve_affine(z: np.ndarray, y: np.ndarray):
    design = np.column_stack([np.ones_like(z), z])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ coeffs
    mse = float(np.mean((pred - y) ** 2))
    return float(coeffs[0]), float(coeffs[1]), pred, mse


def fit_zipf_mandelbrot(ranks: np.ndarray, log_freq: np.ndarray):
    max_rank = float(np.max(ranks))
    c_grid = np.concatenate(
        [
            np.array([0.0], dtype=np.float64),
            np.geomspace(1e-6, max_rank, 2048, dtype=np.float64),
        ]
    )
    best = None
    for c in c_grid:
        z = np.log(ranks + c)
        intercept, slope, pred, mse = _solve_affine(z, log_freq)
        if best is None or mse < best["mse"]:
            best = {
                "a": intercept,
                "b": -slope,
                "c": float(c),
                "mse": mse,
                "rmse": float(math.sqrt(max(mse, 0.0))),
                "prediction": pred,
            }
    return best


def select_log_spaced_indices(n_points: int, sample_points: int):
    if sample_points >= n_points:
        return np.arange(n_points, dtype=np.int64)

    target = min(sample_points, n_points)
    oversample = target
    unique = np.array([], dtype=np.int64)
    while len(unique) < target:
        oversample = max(oversample * 2, target + 8)
        grid = np.exp(np.linspace(0.0, np.log(n_points), oversample)) - 1.0
        idx = np.unique(np.clip(np.round(grid).astype(np.int64), 0, n_points - 1))
        unique = idx
        if oversample > n_points * 16:
            break
    if len(unique) > target:
        pick = np.linspace(0, len(unique) - 1, target)
        unique = unique[np.round(pick).astype(np.int64)]
    return unique


def normalize_x(values: np.ndarray, low: float, high: float):
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if vmax <= vmin:
        return np.full_like(values, (low + high) * 0.5)
    scaled = (values - vmin) / (vmax - vmin)
    return low + (high - low) * scaled


def make_args(cli, depth: int):
    return SimpleNamespace(
        target_fn="",
        depth=depth,
        init_scale=1.0,
        init_strategy="all",
        init_expr="",
        init_blend="",
        init_leaves="",
        init_k=32.0,
        init_noise=0.0,
        seed0=cli.seed0,
        seeds=cli.seeds,
        data_mode="custom-univariate",
        data_lo=1.0,
        data_hi=3.0,
        data_step=0.1,
        gen_lo=0.5,
        gen_hi=5.0,
        generalization_points=4000,
        harmonic_json="",
        harmonic_index=1,
        harmonic_fit_time_min=0.0,
        search_iters=cli.search_iters,
        hardening_iters=cli.hardening_iters,
        lr=0.01,
        tau_search=cli.tau_search,
        tau_hard=cli.tau_hard,
        hardening_tau_power=2.0,
        hardening_lr_floor=0.01,
        patience=4200,
        patience_threshold=1e-2,
        plateau_rtol=1e-3,
        lam_ent_hard=2e-2,
        lam_bin_hard=2e-2,
        lam_inter=1e-4,
        inter_threshold=50.0,
        eml_clamp=1e300,
        eval_every=cli.eval_every,
        tail_eval_every=cli.tail_eval_every,
        tail_eval_tau=0.2,
        early_stop_count=10,
        hard_trigger_mse=1e-20,
        hard_trigger_count=3,
        nan_restart_patience=50,
        max_nan_restarts=100,
        fit_success_thr=1e-6,
        success_thr=1e-20,
        snap_threshold=0.01,
        max_uncertain_success=0,
        lbfgs_steps=0,
        lbfgs_lr=0.6,
        save_prefix="",
        export_m="",
        skip_plot=True,
        loss_y_min=1e-16,
        loss_y_max=1e1,
        plot_dpi=300,
        plot_title_fontsize=13.0,
        plot_label_fontsize=15.0,
        plot_tick_fontsize=12.0,
        plot_legend_fontsize=13.0,
        plot_title="",
    )


def render_snapped_formula(tree):
    leaf_choices = torch.argmax(tree.leaf_logits, dim=1).cpu().numpy().tolist()
    gate_choices = (tree.blend_logits >= 0).to(torch.int64).cpu().numpy()
    leaf_names = tree.leaf_choice_names

    def build_node(level_from_top: int, pos_in_level: int):
        if level_from_top == tree.depth:
            return leaf_names[leaf_choices[pos_in_level]]

        flat_idx = port.flat_node_idx(tree.depth, level_from_top, pos_in_level)
        left_gate = int(gate_choices[flat_idx, 0])
        right_gate = int(gate_choices[flat_idx, 1])
        left_expr = "1" if left_gate else build_node(level_from_top + 1, 2 * pos_in_level)
        right_expr = "1" if right_gate else build_node(level_from_top + 1, 2 * pos_in_level + 1)
        return f"EML[{left_expr},{right_expr}]"

    return build_node(0, 0)


def predict_snapped(snapped_tree, x_values: np.ndarray):
    x_tensor = torch.tensor(x_values, dtype=port.REAL_DTYPE)
    with torch.no_grad():
        pred, _, _, _ = snapped_tree(x_tensor, None, tau_leaf=0.01, tau_gate=0.01)
    return pred.real.detach().cpu().numpy()


def compute_rmse(y_true: np.ndarray, y_pred: np.ndarray):
    mse = float(np.mean((y_true - y_pred) ** 2))
    return mse, float(math.sqrt(max(mse, 0.0)))


def run_depth(dataset_name: str, depth: int, cli, x_values: np.ndarray, y_values: np.ndarray):
    args = make_args(cli, depth)
    x_tensor = torch.tensor(x_values, dtype=port.REAL_DTYPE)
    t_tensor = torch.tensor(y_values, dtype=port.DTYPE)
    run_results = []
    best = None

    for seed in range(cli.seed0, cli.seed0 + cli.seeds):
        for strategy in port.INIT_STRATEGIES_ALL:
            print(f"[{dataset_name}] depth={depth} seed={seed} strategy={strategy}")
            tree, snapped_tree, hist, summary = port.train_one_seed(
                seed=seed,
                strategy=strategy,
                args=args,
                x_train=x_tensor,
                y_train=None,
                t_train=t_tensor,
                manual_init_fn=None,
            )
            prediction = predict_snapped(snapped_tree, x_values)
            mse, rmse = compute_rmse(y_values, prediction)
            result = {
                "seed": seed,
                "strategy": strategy,
                "formula": render_snapped_formula(snapped_tree),
                "summary": summary,
                "best_soft_rmse": float(min(hist["best_soft_rmse"])) if hist["best_soft_rmse"] else None,
                "prediction": prediction.tolist(),
                "rmse": rmse,
                "mse": mse,
            }
            run_results.append(result)
            if best is None or mse < best["mse"]:
                best = result

    return {
        "dataset_name": dataset_name,
        "depth": depth,
        "n_runs": len(run_results),
        "fit_success_count": sum(1 for item in run_results if item["summary"]["fit_success"]),
        "symbol_success_count": sum(1 for item in run_results if item["summary"]["symbol_success"]),
        "stable_symbol_success_count": sum(1 for item in run_results if item["summary"]["stable_symbol_success"]),
        "best_seed": best["seed"],
        "best_strategy": best["strategy"],
        "best_formula": best["formula"],
        "best_mse": best["mse"],
        "best_rmse": best["rmse"],
        "best_prediction": best["prediction"],
        "best_soft_rmse": best["best_soft_rmse"],
        "per_run": run_results,
    }


def run_generalization_check(depth_result, cli, x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, y_test: np.ndarray):
    args = make_args(cli, depth_result["depth"])
    x_train_tensor = torch.tensor(x_train, dtype=port.REAL_DTYPE)
    t_train_tensor = torch.tensor(y_train, dtype=port.DTYPE)
    x_test_tensor = torch.tensor(x_test, dtype=port.REAL_DTYPE)
    t_test_tensor = torch.tensor(y_test, dtype=port.DTYPE)
    tree, snapped_tree, _, summary = port.train_one_seed(
        seed=depth_result["best_seed"],
        strategy=depth_result["best_strategy"],
        args=args,
        x_train=x_train_tensor,
        y_train=None,
        t_train=t_train_tensor,
        manual_init_fn=None,
    )
    train_pred = predict_snapped(snapped_tree, x_train)
    test_pred = predict_snapped(snapped_tree, x_test)
    train_mse, train_rmse = compute_rmse(y_train, train_pred)
    test_mse, test_rmse = compute_rmse(y_test, test_pred)
    return {
        "formula": render_snapped_formula(snapped_tree),
        "summary": summary,
        "train_mse": train_mse,
        "train_rmse": train_rmse,
        "test_mse": test_mse,
        "test_rmse": test_rmse,
        "test_prediction": test_pred.tolist(),
    }


def make_summary_plots(outdir: Path, full_dataset, sampled_dataset, linear_fit, zm_fit, results_by_mode):
    fig, ax = plt.subplots(figsize=(9.2, 5.6), constrained_layout=True)
    ax.scatter(full_dataset["log_rank"], full_dataset["log_freq"], s=5, alpha=0.08, color="#111111", label="full corpus")
    ax.scatter(sampled_dataset["log_rank"], sampled_dataset["log_freq"], s=18, alpha=0.9, color="#1d4ed8", label="sampled points")
    order = np.argsort(sampled_dataset["log_rank"])
    ax.plot(sampled_dataset["log_rank"][order], linear_fit["prediction_sampled"][order], color="#059669", linewidth=1.6, label="linear Zipf")
    ax.plot(sampled_dataset["log_rank"][order], zm_fit["prediction_sampled"][order], color="#b45309", linewidth=1.6, linestyle="--", label="Zipf-Mandelbrot")
    direct_best = results_by_mode["direct"]["best_overall"]
    ax.plot(sampled_dataset["log_rank"][order], np.asarray(direct_best["best_prediction"])[order], color="#dc2626", linewidth=1.4, label=f"best EML direct d{direct_best['depth']}")
    ax.set_title("Zipf Log-Log Curve")
    ax.set_xlabel("ln(rank)")
    ax.set_ylabel("ln(frequency)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(outdir / "zipf_loglog.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.2, 5.6), constrained_layout=True)
    ax.scatter(full_dataset["log_rank"], full_dataset["residual_linear"], s=5, alpha=0.08, color="#111111", label="full residual")
    ax.scatter(sampled_dataset["log_rank"], sampled_dataset["residual_linear"], s=18, alpha=0.9, color="#1d4ed8", label="sampled residual")
    residual_best = results_by_mode["residual"]["best_overall"]
    ax.plot(sampled_dataset["log_rank"][order], np.asarray(residual_best["best_prediction"])[order], color="#7c3aed", linewidth=1.4, label=f"best EML residual d{residual_best['depth']}")
    ax.axhline(0.0, color="#6b7280", linestyle="--", linewidth=1.0)
    ax.set_title("Residual After Linear Zipf Fit")
    ax.set_xlabel("ln(rank)")
    ax.set_ylabel("Residual ln(frequency)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(outdir / "zipf_residual.png", dpi=180)
    plt.close(fig)


def write_report(outdir: Path, corpus_info, sampled_dataset, linear_fit, zm_fit, results_by_mode, generalization):
    lines = [
        "# Zipf EML Experiment",
        "",
        f"- Corpus path: `{corpus_info['corpus_path']}`",
        f"- Tokens: `{corpus_info['token_count']}`",
        f"- Unique words: `{corpus_info['unique_words']}`",
        f"- Sampled points for EML: `{len(sampled_dataset['sample_indices'])}`",
        f"- Linear Zipf fit: `ln(f) = {linear_fit['intercept']:.6f} + ({linear_fit['slope']:.6f}) * ln(rank)`",
        f"- Linear Zipf RMSE: `{linear_fit['rmse_sampled']:.6e}` (sampled), `{linear_fit['rmse_full']:.6e}` (full)",
        f"- Zipf-Mandelbrot fit: `ln(f) = {zm_fit['a']:.6f} - {zm_fit['b']:.6f} * ln(rank + {zm_fit['c']:.6f})`",
        f"- Zipf-Mandelbrot RMSE: `{zm_fit['rmse_sampled']:.6e}` (sampled), `{zm_fit['rmse_full']:.6e}` (full)",
        "",
        "## Top words",
        "",
        f"`{corpus_info['top_words']}`",
        "",
    ]

    for mode_name in ("direct", "residual"):
        if mode_name not in results_by_mode:
            continue
        lines.append(f"## Mode: {mode_name}")
        lines.append("")
        for depth_result in results_by_mode[mode_name]["depth_results"]:
            lines.extend(
                [
                    f"### Depth {depth_result['depth']}",
                    "",
                    f"- runs: `{depth_result['n_runs']}`",
                    f"- fit_success: `{depth_result['fit_success_count']}/{depth_result['n_runs']}`",
                    f"- symbol_success: `{depth_result['symbol_success_count']}/{depth_result['n_runs']}`",
                    f"- stable_symbol_success: `{depth_result['stable_symbol_success_count']}/{depth_result['n_runs']}`",
                    f"- best RMSE: `{depth_result['best_rmse']:.6e}`",
                    f"- best seed/strategy: `{depth_result['best_seed']} / {depth_result['best_strategy']}`",
                    f"- best formula: `{depth_result['best_formula']}`",
                    "",
                ]
            )
        overall = results_by_mode[mode_name]["best_overall"]
        lines.extend(
            [
                f"- Best overall: depth `{overall['depth']}` with RMSE `{overall['best_rmse']:.6e}`",
                "",
            ]
        )

    if generalization:
        lines.append("## Generalization")
        lines.append("")
        for key, value in generalization.items():
            lines.extend(
                [
                    f"### {key}",
                    "",
                    f"- formula: `{value['formula']}`",
                    f"- train RMSE: `{value['train_rmse']:.6e}`",
                    f"- test RMSE: `{value['test_rmse']:.6e}`",
                    "",
                ]
            )

    (outdir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def summarize_mode(depth_results):
    best_overall = min(depth_results, key=lambda item: item["best_mse"])
    return {
        "depth_results": depth_results,
        "best_overall": best_overall,
    }


def main():
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)

    corpus_path = ensure_corpus(cli.corpus_url, cli.corpus_path)
    corpus = build_zipf_dataset(corpus_path)
    sample_indices = select_log_spaced_indices(corpus["unique_words"], cli.sample_points)

    log_rank_sampled = corpus["log_rank"][sample_indices]
    log_freq_sampled = corpus["log_freq"][sample_indices]

    linear_fit = fit_linear_zipf(corpus["log_rank"], corpus["log_freq"])
    linear_fit["prediction_full"] = linear_fit["prediction"]
    linear_fit["prediction_sampled"] = linear_fit["intercept"] + linear_fit["slope"] * log_rank_sampled
    linear_fit["residual_full"] = corpus["log_freq"] - linear_fit["prediction_full"]
    linear_fit["residual_sampled"] = log_freq_sampled - linear_fit["prediction_sampled"]
    linear_fit["mse_full"] = linear_fit["mse"]
    linear_fit["rmse_full"] = linear_fit["rmse"]
    linear_fit["mse_sampled"], linear_fit["rmse_sampled"] = compute_rmse(log_freq_sampled, linear_fit["prediction_sampled"])

    zm_fit = fit_zipf_mandelbrot(corpus["ranks"], corpus["log_freq"])
    zm_fit["prediction_full"] = zm_fit["prediction"]
    zm_fit["prediction_sampled"] = zm_fit["a"] - zm_fit["b"] * np.log(corpus["ranks"][sample_indices] + zm_fit["c"])
    zm_fit["mse_full"] = zm_fit["mse"]
    zm_fit["rmse_full"] = zm_fit["rmse"]
    zm_fit["mse_sampled"], zm_fit["rmse_sampled"] = compute_rmse(log_freq_sampled, zm_fit["prediction_sampled"])

    x_sampled = normalize_x(log_rank_sampled, cli.x_low, cli.x_high)
    sampled_dataset = {
        "sample_indices": sample_indices.tolist(),
        "ranks": corpus["ranks"][sample_indices].tolist(),
        "freqs": corpus["freqs"][sample_indices].tolist(),
        "log_rank": log_rank_sampled,
        "log_freq": log_freq_sampled,
        "x_eml": x_sampled,
        "residual_linear": linear_fit["residual_sampled"],
    }
    full_dataset = {
        "log_rank": corpus["log_rank"],
        "log_freq": corpus["log_freq"],
        "residual_linear": linear_fit["residual_full"],
    }

    results_direct = []
    results_residual = []
    for depth in cli.depths:
        if "direct" in cli.modes:
            results_direct.append(run_depth("direct", depth, cli, x_sampled, log_freq_sampled))
        if "residual" in cli.modes:
            results_residual.append(run_depth("residual", depth, cli, x_sampled, linear_fit["residual_sampled"]))

    results_by_mode = {}
    if results_direct:
        results_by_mode["direct"] = summarize_mode(results_direct)
    if results_residual:
        results_by_mode["residual"] = summarize_mode(results_residual)

    generalization = {}
    mid = len(sample_indices) // 2
    x_train = x_sampled[:mid]
    x_test = x_sampled[mid:]
    y_direct_train = log_freq_sampled[:mid]
    y_direct_test = log_freq_sampled[mid:]
    y_resid_train = linear_fit["residual_sampled"][:mid]
    y_resid_test = linear_fit["residual_sampled"][mid:]
    for mode_name, y_train, y_test in (
        ("direct", y_direct_train, y_direct_test),
        ("residual", y_resid_train, y_resid_test),
    ):
        if mode_name not in results_by_mode:
            continue
        if any(item["stable_symbol_success_count"] > 0 for item in results_by_mode[mode_name]["depth_results"]):
            generalization[mode_name] = run_generalization_check(
                results_by_mode[mode_name]["best_overall"],
                cli,
                x_train,
                y_train,
                x_test,
                y_test,
            )

    corpus_info = {
        "corpus_path": str(corpus_path),
        "token_count": corpus["token_count"],
        "unique_words": corpus["unique_words"],
        "top_words": corpus["top_words"],
    }

    payload = {
        "corpus": corpus_info,
        "sampled_dataset": {
            "sample_indices": sampled_dataset["sample_indices"],
            "ranks": sampled_dataset["ranks"],
            "freqs": sampled_dataset["freqs"],
            "log_rank": sampled_dataset["log_rank"].tolist(),
            "log_freq": sampled_dataset["log_freq"].tolist(),
            "x_eml": sampled_dataset["x_eml"].tolist(),
            "residual_linear": sampled_dataset["residual_linear"].tolist(),
        },
        "linear_zipf": {
            "intercept": linear_fit["intercept"],
            "slope": linear_fit["slope"],
            "mse_full": linear_fit["mse_full"],
            "rmse_full": linear_fit["rmse_full"],
            "mse_sampled": linear_fit["mse_sampled"],
            "rmse_sampled": linear_fit["rmse_sampled"],
            "prediction_sampled": linear_fit["prediction_sampled"].tolist(),
        },
        "zipf_mandelbrot": {
            "a": zm_fit["a"],
            "b": zm_fit["b"],
            "c": zm_fit["c"],
            "mse_full": zm_fit["mse_full"],
            "rmse_full": zm_fit["rmse_full"],
            "mse_sampled": zm_fit["mse_sampled"],
            "rmse_sampled": zm_fit["rmse_sampled"],
            "prediction_sampled": zm_fit["prediction_sampled"].tolist(),
        },
        "results_by_mode": results_by_mode,
        "generalization": generalization,
        "args": sanitize_jsonable(port._sanitize_for_json(vars(cli))),
    }
    summary_path = cli.outdir / "summary.json"
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    write_report(cli.outdir, corpus_info, sampled_dataset, linear_fit, zm_fit, results_by_mode, generalization)
    if not cli.skip_plot and "direct" in results_by_mode and "residual" in results_by_mode:
        make_summary_plots(cli.outdir, full_dataset, sampled_dataset, linear_fit, zm_fit, results_by_mode)

    print(f"Saved {summary_path}")
    print(f"Saved {cli.outdir / 'report.md'}")
    if not cli.skip_plot and "direct" in results_by_mode and "residual" in results_by_mode:
        print(f"Saved {cli.outdir / 'zipf_loglog.png'}")
        print(f"Saved {cli.outdir / 'zipf_residual.png'}")


if __name__ == "__main__":
    main()
