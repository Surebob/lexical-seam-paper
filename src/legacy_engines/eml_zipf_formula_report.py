import argparse
import json
import math
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
DEFAULT_SUMMARY = ROOT / "results" / "zipf_shakespeare_full" / "summary.json"


def parse_args():
    p = argparse.ArgumentParser(description="Extract and compare discovered Zipf EML formulas")
    p.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    p.add_argument("--outdir", type=Path, default=None)
    p.add_argument("--escape-rmse-threshold", type=float, default=2.0)
    p.add_argument("--depth4-rmse-threshold", type=float, default=1.0)
    return p.parse_args()


def load_summary(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def eml_to_math(expr: str):
    simplified = {
        "EML[EML[1,EML[x,1]],EML[EML[x,1],EML[1,x]]]": "exp(e - x) - log(exp(exp(x)) - log(e - log(x)))",
        "EML[1,EML[EML[1,x],EML[x,1]]]": "e - log(e^e / x - x)",
        "EML[EML[1,EML[EML[1,x],EML[1,1]]],EML[EML[1,EML[1,1]],EML[EML[1,1],EML[x,1]]]]": "exp(e - log(e^e / x - 1)) - log(exp(e - 1) - log(e^e - x))",
        "EML[EML[1,EML[EML[1,x],EML[1,1]]],EML[1,1]]": "exp(e - log(e^e / x - 1)) - 1",
        "EML[EML[1,EML[EML[1,x],EML[1,x]]],EML[1,1]]": "exp(e - log(exp(e - log(x)) - log(e - log(x)))) - 1",
    }
    if expr in simplified:
        return simplified[expr]

    def parse_node(s: str):
        s = s.strip()
        if s in {"1", "x", "y"}:
            return s
        if not s.startswith("EML[") or not s.endswith("]"):
            return s
        inner = s[4:-1]
        depth = 0
        split = None
        for idx, ch in enumerate(inner):
            if ch == "[":
                depth += 1
            elif ch == "]":
                depth -= 1
            elif ch == "," and depth == 0:
                split = idx
                break
        if split is None:
            return s
        left = parse_node(inner[:split])
        right = parse_node(inner[split + 1 :])
        return ("EML", left, right)

    def render(node):
        if isinstance(node, str):
            return node
        _, left, right = node
        return f"(exp({render(left)}) - log({render(right)}))"

    return render(parse_node(expr))


def rmse(y_true, y_pred):
    return float(math.sqrt(max(float(np.mean((y_true - y_pred) ** 2)), 0.0)))


def parse_expr(expr: str):
    expr = expr.strip()
    if expr in {"1", "x", "y"}:
        return expr
    if not expr.startswith("EML[") or not expr.endswith("]"):
        return expr
    inner = expr[4:-1]
    depth = 0
    split = None
    for idx, ch in enumerate(inner):
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        elif ch == "," and depth == 0:
            split = idx
            break
    if split is None:
        return expr
    return ("EML", parse_expr(inner[:split]), parse_expr(inner[split + 1 :]))


def eval_expr(node, x_values):
    if node == "1":
        return np.ones_like(x_values, dtype=np.complex128)
    if node == "x":
        return x_values.astype(np.complex128)
    _, left, right = node
    left_val = eval_expr(left, x_values)
    right_val = eval_expr(right, x_values)
    return np.exp(left_val) - np.log(right_val)


def build_full_dataset(summary):
    corpus_path = Path(summary["corpus"]["corpus_path"])
    text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    for marker in [
        "*** START OF THE PROJECT GUTENBERG EBOOK",
        "*** START OF THIS PROJECT GUTENBERG EBOOK",
    ]:
        idx = text.find(marker)
        if idx != -1:
            line_end = text.find("\n", idx)
            text = text[line_end + 1 :] if line_end != -1 else text[idx:]
            break
    for marker in [
        "*** END OF THE PROJECT GUTENBERG EBOOK",
        "*** END OF THIS PROJECT GUTENBERG EBOOK",
    ]:
        idx = text.find(marker)
        if idx != -1:
            text = text[:idx]
            break

    import re
    from collections import Counter

    words = re.findall(r"[a-z]+(?:'[a-z]+)?", text.lower())
    counts = Counter(words)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)

    sampled = summary["sampled_dataset"]
    sample_log_rank = np.asarray(sampled["log_rank"], dtype=np.float64)
    sample_x = np.asarray(sampled["x_eml"], dtype=np.float64)
    lo = float(sample_log_rank.min())
    hi = float(sample_log_rank.max())
    x_lo = float(sample_x.min())
    x_hi = float(sample_x.max())
    x_full = x_lo + (x_hi - x_lo) * (log_rank - lo) / (hi - lo)

    return {
        "ranks": ranks,
        "log_rank": log_rank,
        "log_freq": log_freq,
        "x_full": x_full,
    }


def evaluate_formula(summary, run, full_dataset):
    sampled = summary["sampled_dataset"]
    sample_x = np.asarray(sampled["x_eml"], dtype=np.float64)
    sample_log_freq = np.asarray(sampled["log_freq"], dtype=np.float64)
    linear_sample = np.asarray(summary["linear_zipf"]["prediction_sampled"], dtype=np.float64)
    formula_node = parse_expr(run["formula"])

    pred_sample = eval_expr(formula_node, sample_x).real
    pred_full = eval_expr(formula_node, full_dataset["x_full"]).real
    if run["mode"] == "residual":
        pred_sample_total = linear_sample + pred_sample
        linear_full = summary["linear_zipf"]["intercept"] + summary["linear_zipf"]["slope"] * full_dataset["log_rank"]
        pred_full_total = linear_full + pred_full
    else:
        pred_sample_total = pred_sample
        pred_full_total = pred_full

    mid = len(sample_log_freq) // 2
    return {
        "sampled_total_rmse": rmse(sample_log_freq, pred_sample_total),
        "sampled_head_rmse": rmse(sample_log_freq[:mid], pred_sample_total[:mid]),
        "sampled_tail_rmse": rmse(sample_log_freq[mid:], pred_sample_total[mid:]),
        "full_total_rmse": rmse(full_dataset["log_freq"], pred_full_total),
        "prediction_sample_total": pred_sample_total.tolist(),
    }


def summarize_strategy(depth_result):
    summary = {}
    for run in depth_result["per_run"]:
        key = run["strategy"]
        summary.setdefault(
            key,
            {
                "runs": 0,
                "fit_successes": 0,
                "symbol_successes": 0,
                "stable_symbol_successes": 0,
                "escapes": [],
            },
        )
        bucket = summary[key]
        bucket["runs"] += 1
        bucket["fit_successes"] += int(bool(run["summary"]["fit_success"]))
        bucket["symbol_successes"] += int(bool(run["summary"]["symbol_success"]))
        bucket["stable_symbol_successes"] += int(bool(run["summary"]["stable_symbol_success"]))
    return summary


def add_escape_runs(strategy_summary, depth_result, threshold):
    for run in depth_result["per_run"]:
        if float(run["rmse"]) <= threshold:
            strategy_summary[run["strategy"]]["escapes"].append(
                {
                    "seed": run["seed"],
                    "rmse": run["rmse"],
                    "formula_eml": run["formula"],
                    "formula_math": eml_to_math(run["formula"]),
                }
            )


def make_overlay_plot(outdir: Path, summary, selected_runs):
    sampled = summary["sampled_dataset"]
    x = np.asarray(sampled["log_rank"], dtype=np.float64)
    y = np.asarray(sampled["log_freq"], dtype=np.float64)
    order = np.argsort(x)
    linear_pred = np.asarray(summary["linear_zipf"]["prediction_sampled"], dtype=np.float64)
    zm_pred = np.asarray(summary["zipf_mandelbrot"]["prediction_sampled"], dtype=np.float64)

    fig, ax = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)
    ax.scatter(x, y, s=18, color="#111111", alpha=0.8, label="sampled Zipf data")
    ax.plot(x[order], linear_pred[order], color="#059669", linewidth=1.8, label="linear Zipf")
    ax.plot(x[order], zm_pred[order], color="#b45309", linewidth=1.8, linestyle="--", label="Zipf-Mandelbrot")

    palette = ["#dc2626", "#2563eb", "#7c3aed", "#0891b2"]
    for idx, run in enumerate(selected_runs):
        pred = np.asarray(run["prediction_sample_total"], dtype=np.float64)
        ax.plot(
            x[order],
            pred[order],
            linewidth=1.5,
            color=palette[idx % len(palette)],
            label=f"{run['label']} ({run['sampled_total_rmse']:.3f})",
        )

    ax.set_title("Zipf Data With Analytic Baselines and Selected EML Fits")
    ax.set_xlabel("ln(rank)")
    ax.set_ylabel("ln(frequency)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(outdir / "formula_overlay.png", dpi=180)
    plt.close(fig)


def write_report(outdir: Path, summary, direct_d3, selected_runs, strategy_summary_d3):
    linear_rmse = float(summary["linear_zipf"]["rmse_sampled"])
    zm_rmse = float(summary["zipf_mandelbrot"]["rmse_sampled"])
    lines = [
        "# Zipf Formula Extraction",
        "",
        f"- Source summary: `{summary['args'].get('outdir', '')}`",
        f"- Linear Zipf sampled RMSE: `{linear_rmse:.6e}`",
        f"- Zipf-Mandelbrot sampled RMSE: `{zm_rmse:.6e}`",
        "",
        "## Depth-3 Direct Escape Statistics",
        "",
    ]

    for strategy, bucket in strategy_summary_d3.items():
        seeds = [item["seed"] for item in bucket["escapes"]]
        lines.extend(
            [
                f"### {strategy}",
                "",
                f"- runs: `{bucket['runs']}`",
                f"- fit_successes: `{bucket['fit_successes']}`",
                f"- symbol_successes: `{bucket['symbol_successes']}`",
                f"- stable_symbol_successes: `{bucket['stable_symbol_successes']}`",
                f"- escapes (RMSE <= threshold): `{seeds}`",
                "",
            ]
        )
        for item in bucket["escapes"]:
            lines.extend(
                [
                f"- seed `{item['seed']}` RMSE `{item['rmse']:.6e}`",
                f"  EML: `{item['formula_eml']}`",
                f"  Math: `{item['formula_math']}`",
            ]
        )
        if bucket["escapes"]:
            lines.append("")

    lines.extend(
        [
            "## Selected Formulas",
            "",
        ]
    )
    for run in selected_runs:
        lines.extend(
            [
                f"### {run['label']}",
                "",
                f"- seeds: `{run['seeds']}`",
                f"- strategies: `{run['strategies']}`",
                f"- mode/depth: `{run['mode']} / {run['depth']}`",
                f"- sampled total RMSE: `{run['sampled_total_rmse']:.6e}`",
                f"- full total RMSE: `{run['full_total_rmse']:.6e}`",
                f"- sampled head/tail RMSE: `{run['sampled_head_rmse']:.6e}` / `{run['sampled_tail_rmse']:.6e}`",
                f"- EML: `{run['formula']}`",
                f"- Math: `{eml_to_math(run['formula'])}`",
                "",
            ]
        )

    lines.extend(
        [
            "## RMSE Comparison",
            "",
            "| Model | Sampled RMSE |",
            "| --- | ---: |",
            f"| Linear Zipf (sampled) | {linear_rmse:.6f} |",
            f"| Zipf-Mandelbrot (sampled) | {zm_rmse:.6f} |",
        ]
    )
    for run in selected_runs:
        lines.append(f"| {run['label']} | {run['sampled_total_rmse']:.6f} |")

    lines.extend(
        [
            "",
            "| Model | Full-corpus RMSE |",
            "| --- | ---: |",
            f"| Linear Zipf | {summary['linear_zipf_full_rmse']:.6f} |",
            f"| Zipf-Mandelbrot | {summary['zipf_mandelbrot_full_rmse']:.6f} |",
        ]
    )
    for run in selected_runs:
        lines.append(f"| {run['label']} | {run['full_total_rmse']:.6f} |")

    (outdir / "formula_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    cli = parse_args()
    summary = load_summary(cli.summary)
    outdir = cli.outdir or cli.summary.parent
    outdir.mkdir(parents=True, exist_ok=True)
    full_dataset = build_full_dataset(summary)

    summary["linear_zipf_full_rmse"] = rmse(
        full_dataset["log_freq"],
        summary["linear_zipf"]["intercept"] + summary["linear_zipf"]["slope"] * full_dataset["log_rank"],
    )
    summary["zipf_mandelbrot_full_rmse"] = rmse(
        full_dataset["log_freq"],
        summary["zipf_mandelbrot"]["a"] - summary["zipf_mandelbrot"]["b"] * np.log(full_dataset["ranks"] + summary["zipf_mandelbrot"]["c"]),
    )

    direct_depths = {item["depth"]: item for item in summary["results_by_mode"]["direct"]["depth_results"]}
    residual_depths = {item["depth"]: item for item in summary["results_by_mode"]["residual"]["depth_results"]}

    direct_d3 = direct_depths.get(3)
    if direct_d3 is None:
        raise ValueError("Depth-3 direct results not found in summary.")

    strategy_summary_d3 = summarize_strategy(direct_d3)
    add_escape_runs(strategy_summary_d3, direct_d3, cli.escape_rmse_threshold)

    selected_runs = []
    for run in direct_d3["per_run"]:
        if float(run["rmse"]) <= cli.escape_rmse_threshold:
            enriched = {"label": f"direct d3 seed {run['seed']} {run['strategy']}", "mode": "direct", "depth": 3, **run}
            enriched.update(evaluate_formula(summary, enriched, full_dataset))
            selected_runs.append(enriched)

    direct_d4 = direct_depths.get(4)
    if direct_d4 is not None:
        for run in direct_d4["per_run"]:
            if float(run["rmse"]) <= cli.depth4_rmse_threshold:
                enriched = {"label": f"direct d4 seed {run['seed']} {run['strategy']}", "mode": "direct", "depth": 4, **run}
                enriched.update(evaluate_formula(summary, enriched, full_dataset))
                selected_runs.append(enriched)

    residual_d3 = residual_depths.get(3)
    if residual_d3 is not None:
        for run in residual_d3["per_run"]:
            if float(run["rmse"]) <= cli.depth4_rmse_threshold:
                enriched = {"label": f"residual d3 seed {run['seed']} {run['strategy']}", "mode": "residual", "depth": 3, **run}
                enriched.update(evaluate_formula(summary, enriched, full_dataset))
                selected_runs.append(enriched)

    residual_d4 = residual_depths.get(4)
    if residual_d4 is not None:
        for run in residual_d4["per_run"]:
            if float(run["rmse"]) <= cli.depth4_rmse_threshold:
                enriched = {"label": f"residual d4 seed {run['seed']} {run['strategy']}", "mode": "residual", "depth": 4, **run}
                enriched.update(evaluate_formula(summary, enriched, full_dataset))
                selected_runs.append(enriched)

    selected_runs.sort(key=lambda item: item["sampled_total_rmse"])
    grouped = {}
    for item in selected_runs:
        key = (item["mode"], item["depth"], item["formula"])
        bucket = grouped.setdefault(
            key,
            {
                "label": f"{item['mode']} d{item['depth']}",
                "mode": item["mode"],
                "depth": item["depth"],
                "formula": item["formula"],
                "sampled_total_rmse": item["sampled_total_rmse"],
                "full_total_rmse": item["full_total_rmse"],
                "sampled_head_rmse": item["sampled_head_rmse"],
                "sampled_tail_rmse": item["sampled_tail_rmse"],
                "prediction_sample_total": item["prediction_sample_total"],
                "seeds": [],
                "strategies": [],
            },
        )
        bucket["seeds"].append(item["seed"])
        if item["strategy"] not in bucket["strategies"]:
            bucket["strategies"].append(item["strategy"])

    deduped_runs = []
    for bucket in grouped.values():
        bucket["seeds"] = sorted(set(bucket["seeds"]))
        bucket["label"] = f"{bucket['mode']} d{bucket['depth']} ({len(bucket['seeds'])} runs)"
        deduped_runs.append(bucket)

    deduped_runs.sort(key=lambda item: item["sampled_total_rmse"])
    make_overlay_plot(outdir, summary, deduped_runs[:4])
    write_report(outdir, summary, direct_d3, deduped_runs, strategy_summary_d3)

    payload = {
        "direct_depth3_strategy_summary": strategy_summary_d3,
        "selected_runs": deduped_runs,
    }
    (outdir / "formula_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Saved {outdir / 'formula_report.md'}")
    print(f"Saved {outdir / 'formula_report.json'}")
    print(f"Saved {outdir / 'formula_overlay.png'}")


if __name__ == "__main__":
    main()
