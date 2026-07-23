import argparse
import functools
import hashlib
import importlib.util
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
ZIPF_EXPERIMENT_PATH = ROOT / "eml_zipf_experiment.py"
FORMULA_REPORT_PATH = ROOT / "eml_zipf_formula_report.py"
DEFAULT_CORPUS_URL = "https://www.gutenberg.org/cache/epub/100/pg100.txt"
DEFAULT_CORPUS_PATH = ROOT / "data" / "zipf" / "pg100.txt"
DEFAULT_OUTDIR = ROOT / "results" / "zipf_discrete_search"
EXP_CLAMP = 60.0
LOG_EPS = 1e-12
LOG_DEN_EPS = 1e-9
SIG_ROUND = 10
TERNARY_SAFE_X_HIGH = 0.95


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


zipf = load_module(ZIPF_EXPERIMENT_PATH, "zipf_discrete_search_zipf")
formula_tools = load_module(FORMULA_REPORT_PATH, "zipf_discrete_search_formula")


def parse_args():
    p = argparse.ArgumentParser(description="Discrete greedy/beam EML search on Zipf residuals")
    p.add_argument("--corpus-url", type=str, default=DEFAULT_CORPUS_URL)
    p.add_argument("--corpus-path", type=Path, default=DEFAULT_CORPUS_PATH)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--sample-points", type=int, default=2000)
    p.add_argument("--beam-width", type=int, default=5)
    p.add_argument("--max-steps", type=int, default=20)
    p.add_argument("--keep-all-until-step", type=int, default=2)
    p.add_argument("--search-strategy", type=str, default="library_pairwise", choices=["library_pairwise", "primitive_beam"])
    p.add_argument("--diversity-weight", type=float, default=0.35)
    p.add_argument("--x-low", type=float, default=0.05)
    p.add_argument("--x-high", type=float, default=1.0)
    p.add_argument("--mode", type=str, default="residual", choices=["residual", "direct"])
    p.add_argument("--operator", type=str, default="eml", choices=["eml", "ternary"])
    p.add_argument("--constant-variance-threshold", type=float, default=1e-10)
    p.add_argument("--improvement-epsilon", type=float, default=1e-10)
    p.add_argument("--skip-plot", action="store_true")
    p.add_argument("--residual-baseline", type=str, default="linear", choices=["linear", "zipf_mandelbrot"])
    return p.parse_args()


def build_bundle(corpus_url: str, corpus_path: Path, sample_points: int, x_low: float, x_high: float):
    corpus_path = zipf.ensure_corpus(corpus_url, corpus_path)
    corpus = zipf.build_zipf_dataset(corpus_path)

    if sample_points <= 0 or sample_points >= corpus["unique_words"]:
        sample_indices = np.arange(corpus["unique_words"], dtype=np.int64)
    else:
        sample_indices = zipf.select_log_spaced_indices(corpus["unique_words"], sample_points)

    log_rank_sampled = corpus["log_rank"][sample_indices]
    log_freq_sampled = corpus["log_freq"][sample_indices]

    linear_fit = zipf.fit_linear_zipf(corpus["log_rank"], corpus["log_freq"])
    linear_fit["prediction_full"] = linear_fit["prediction"]
    linear_fit["prediction_sampled"] = linear_fit["intercept"] + linear_fit["slope"] * log_rank_sampled
    linear_fit["residual_full"] = corpus["log_freq"] - linear_fit["prediction_full"]
    linear_fit["residual_sampled"] = log_freq_sampled - linear_fit["prediction_sampled"]
    linear_fit["mse_full"] = linear_fit["mse"]
    linear_fit["rmse_full"] = linear_fit["rmse"]
    linear_fit["mse_sampled"], linear_fit["rmse_sampled"] = zipf.compute_rmse(
        log_freq_sampled, linear_fit["prediction_sampled"]
    )

    zm_fit = zipf.fit_zipf_mandelbrot(corpus["ranks"], corpus["log_freq"])
    zm_fit["prediction_full"] = zm_fit["prediction"]
    zm_fit["prediction_sampled"] = zm_fit["a"] - zm_fit["b"] * np.log(corpus["ranks"][sample_indices] + zm_fit["c"])
    zm_fit["mse_full"] = zm_fit["mse"]
    zm_fit["rmse_full"] = zm_fit["rmse"]
    zm_fit["mse_sampled"], zm_fit["rmse_sampled"] = zipf.compute_rmse(log_freq_sampled, zm_fit["prediction_sampled"])

    x_sampled = zipf.normalize_x(log_rank_sampled, x_low, x_high)

    lo = float(log_rank_sampled.min())
    hi = float(log_rank_sampled.max())
    x_lo = float(x_sampled.min())
    x_hi = float(x_sampled.max())
    if hi <= lo:
        x_full = np.full_like(corpus["log_rank"], (x_low + x_high) * 0.5)
    else:
        x_full = x_lo + (x_hi - x_lo) * (corpus["log_rank"] - lo) / (hi - lo)

    return {
        "corpus_path": corpus_path,
        "corpus": corpus,
        "sample_indices": sample_indices,
        "log_rank_sampled": log_rank_sampled,
        "log_freq_sampled": log_freq_sampled,
        "x_sampled": x_sampled,
        "x_full": x_full,
        "linear_fit": linear_fit,
        "zm_fit": zm_fit,
    }


def get_residual_bundle(bundle, residual_baseline: str):
    if residual_baseline == "zipf_mandelbrot":
        baseline_sample = bundle["zm_fit"]["prediction_sampled"]
        baseline_full = bundle["zm_fit"]["prediction_full"]
        baseline_name = "zipf_mandelbrot"
    else:
        baseline_sample = bundle["linear_fit"]["prediction_sampled"]
        baseline_full = bundle["linear_fit"]["prediction_full"]
        baseline_name = "linear"
    return {
        "name": baseline_name,
        "prediction_sampled": np.asarray(baseline_sample, dtype=np.float64),
        "prediction_full": np.asarray(baseline_full, dtype=np.float64),
        "residual_sampled": bundle["log_freq_sampled"] - np.asarray(baseline_sample, dtype=np.float64),
        "residual_full": bundle["corpus"]["log_freq"] - np.asarray(baseline_full, dtype=np.float64),
        "rmse_sampled": formula_tools.rmse(bundle["log_freq_sampled"], baseline_sample),
        "rmse_full": formula_tools.rmse(bundle["corpus"]["log_freq"], baseline_full),
    }


def safe_eml(left: np.ndarray, right: np.ndarray):
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        return None
    if np.any(right <= LOG_EPS):
        return None
    out = np.exp(np.clip(left, -EXP_CLAMP, EXP_CLAMP)) - np.log(right)
    if not np.all(np.isfinite(out)):
        return None
    return out


def safe_ternary(left: np.ndarray, middle: np.ndarray, right: np.ndarray):
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(middle)) or not np.all(np.isfinite(right)):
        return None
    if np.any(left <= LOG_EPS) or np.any(right <= LOG_EPS):
        return None
    log_left = np.log(left)
    if np.any(~np.isfinite(log_left)) or np.any(np.abs(log_left) <= LOG_DEN_EPS):
        return None
    log_right = np.log(right)
    if np.any(~np.isfinite(log_right)):
        return None
    out = np.exp(np.clip(left - middle, -EXP_CLAMP, EXP_CLAMP)) * (log_right / log_left)
    if not np.all(np.isfinite(out)):
        return None
    return out


def parse_expr(expr: str):
    def parse_at(index: int):
        if index >= len(expr):
            raise ValueError(f"Unexpected end of expression: {expr}")

        if expr[index].isalpha():
            end = index
            while end < len(expr) and expr[end].isalpha():
                end += 1
            token = expr[index:end]
            if end < len(expr) and expr[end] == "[":
                end += 1
                children = []
                while True:
                    child, end = parse_at(end)
                    children.append(child)
                    if end >= len(expr):
                        raise ValueError(f"Unterminated expression: {expr}")
                    if expr[end] == ",":
                        end += 1
                        continue
                    if expr[end] == "]":
                        end += 1
                        break
                    raise ValueError(f"Unexpected token in expression {expr!r} at offset {end}")
                return (token, *children), end
            return token, end

        if expr[index].isdigit():
            end = index
            while end < len(expr) and expr[end].isdigit():
                end += 1
            return expr[index:end], end

        raise ValueError(f"Unexpected token in expression {expr!r} at offset {index}")

    node, offset = parse_at(0)
    if offset != len(expr):
        raise ValueError(f"Trailing content in expression {expr!r} at offset {offset}")
    return node


def render_expr(node):
    if isinstance(node, str):
        return node
    op = node[0]
    children = ",".join(render_expr(child) for child in node[1:])
    return f"{op}[{children}]"


def expr_to_math(node):
    if isinstance(node, str):
        return node
    op = node[0]
    if op == "EML":
        left = expr_to_math(node[1])
        right = expr_to_math(node[2])
        return f"(exp({left}) - log({right}))"
    if op == "T":
        left = expr_to_math(node[1])
        middle = expr_to_math(node[2])
        right = expr_to_math(node[3])
        return f"((exp({left}) / log({left})) * (log({right}) / exp({middle})))"
    raise ValueError(f"Unsupported operator {op!r}")


def format_expr(operator: str, *children: str):
    if operator == "eml":
        return f"EML[{children[0]},{children[1]}]"
    if operator == "ternary":
        return f"T[{children[0]},{children[1]},{children[2]}]"
    raise ValueError(f"Unsupported operator {operator!r}")


def apply_operator(operator: str, *values: np.ndarray):
    if operator == "eml":
        return safe_eml(values[0], values[1])
    if operator == "ternary":
        return safe_ternary(values[0], values[1], values[2])
    raise ValueError(f"Unsupported operator {operator!r}")


def value_signature(values: np.ndarray):
    rounded = np.round(np.asarray(values, dtype=np.float64), SIG_ROUND)
    return hashlib.sha1(rounded.tobytes()).hexdigest()


@functools.lru_cache(maxsize=None)
def expr_feature_set(expr: str):
    node = parse_expr(expr)
    features = set()

    def walk(cur):
        if isinstance(cur, str):
            features.add(cur)
            return
        expr_cur = render_expr(cur)
        features.add(expr_cur)
        for child in cur[1:]:
            walk(child)

    walk(node)
    return features


def expr_similarity(expr_a: str, expr_b: str):
    fa = expr_feature_set(expr_a)
    fb = expr_feature_set(expr_b)
    union = len(fa | fb)
    if union == 0:
        return 0.0
    return len(fa & fb) / union


def make_candidate(expr: str, sample_values: np.ndarray, full_values: np.ndarray, bundle, mode: str, step: int):
    if mode == "residual":
        sample_total = bundle["residual_fit"]["prediction_sampled"] + sample_values
        full_total = bundle["residual_fit"]["prediction_full"] + full_values
        target_sample = bundle["log_freq_sampled"]
        target_full = bundle["corpus"]["log_freq"]
    else:
        sample_total = sample_values
        full_total = full_values
        target_sample = bundle["log_freq_sampled"]
        target_full = bundle["corpus"]["log_freq"]

    sampled_rmse = formula_tools.rmse(target_sample, sample_total)
    full_rmse = formula_tools.rmse(target_full, full_total)
    return {
        "expr": expr,
        "math": expr_to_math(parse_expr(expr)),
        "step": step,
        "sample_values": sample_values,
        "full_values": full_values,
        "sampled_rmse": sampled_rmse,
        "full_rmse": full_rmse,
        "sample_variance": float(np.var(sample_values)),
        "full_variance": float(np.var(full_values)),
        "sample_total": sample_total,
        "full_total": full_total,
        "signature": value_signature(sample_values),
    }


def dedupe_candidates(candidates):
    best_by_signature = {}
    for candidate in candidates:
        incumbent = best_by_signature.get(candidate["signature"])
        if incumbent is None or (
            candidate["sampled_rmse"],
            candidate["full_rmse"],
            candidate["expr"],
        ) < (
            incumbent["sampled_rmse"],
            incumbent["full_rmse"],
            incumbent["expr"],
        ):
            best_by_signature[candidate["signature"]] = candidate
    return sorted(best_by_signature.values(), key=lambda item: (item["sampled_rmse"], item["full_rmse"], item["expr"]))


def is_effectively_constant(candidate, threshold: float):
    return candidate["full_variance"] < threshold


def filter_nonconstant(candidates, threshold: float):
    return [candidate for candidate in candidates if not is_effectively_constant(candidate, threshold)]


def initial_library(bundle, mode: str, operator: str):
    ones_sample = np.ones_like(bundle["x_sampled"], dtype=np.float64)
    ones_full = np.ones_like(bundle["x_full"], dtype=np.float64)
    x_sample = np.asarray(bundle["x_sampled"], dtype=np.float64)
    x_full = np.asarray(bundle["x_full"], dtype=np.float64)
    if operator == "ternary":
        return [make_candidate("x", x_sample, x_full, bundle, mode, step=0)]
    return [
        make_candidate("1", ones_sample, ones_full, bundle, mode, step=0),
        make_candidate("x", x_sample, x_full, bundle, mode, step=0),
    ]


def iterate_library_combinations(library, operator: str):
    if operator == "ternary":
        for first in library:
            for second in library:
                for third in library:
                    yield (first, second, third)
        return
    for left in library:
        for right in library:
            yield (left, right)


def iterate_beam_extensions(current, library, operator: str):
    if operator == "ternary":
        for first in library:
            for second in library:
                combos = {
                    (current["expr"], first["expr"], second["expr"]): (current, first, second),
                    (first["expr"], current["expr"], second["expr"]): (first, current, second),
                    (first["expr"], second["expr"], current["expr"]): (first, second, current),
                }
                for combo in combos.values():
                    yield combo
        return

    for primitive in library:
        yield (current, primitive)
        yield (primitive, current)


def bootstrap_library(bundle, mode: str, keep_all_until_step: int, operator: str, constant_variance_threshold: float):
    temp_library = initial_library(bundle, mode, operator)
    seen_expr = {item["expr"] for item in temp_library}
    seen_sig = {item["signature"] for item in temp_library}
    raw_additions = []

    for step in range(1, keep_all_until_step + 1):
        generated = []
        for combo in iterate_library_combinations(temp_library, operator):
                expr = format_expr(operator, *(item["expr"] for item in combo))
                if expr in seen_expr:
                    continue
                sample_values = apply_operator(operator, *(item["sample_values"] for item in combo))
                if sample_values is None:
                    continue
                full_values = apply_operator(operator, *(item["full_values"] for item in combo))
                if full_values is None:
                    continue
                sig = value_signature(sample_values)
                if sig in seen_sig:
                    continue
                candidate = make_candidate(expr, sample_values, full_values, bundle, mode, step=step)
                generated.append(candidate)

        if not generated:
            break

        generated = dedupe_candidates(generated)
        chosen = []
        local_seen = set()
        for candidate in generated:
            if candidate["signature"] in local_seen:
                continue
            local_seen.add(candidate["signature"])
            chosen.append(candidate)

        if not chosen:
            break

        for candidate in chosen:
            temp_library.append(candidate)
            raw_additions.append(candidate)
            seen_expr.add(candidate["expr"])
            seen_sig.add(candidate["signature"])

    filtered_library = [item for item in temp_library if item["expr"] == "x" or not is_effectively_constant(item, constant_variance_threshold)]
    filtered_additions = filter_nonconstant(raw_additions, constant_variance_threshold)

    return {
        "library": filtered_library,
        "additions": filtered_additions,
        "seed_candidates": filtered_additions,
    }


def select_diverse_beam(candidates, beam_width: int, diversity_weight: float):
    if not candidates:
        return []
    if len(candidates) <= beam_width:
        return candidates

    rmses = np.asarray([item["sampled_rmse"] for item in candidates], dtype=np.float64)
    lo = float(rmses.min())
    hi = float(rmses.max())
    if hi <= lo:
        quality_scores = {item["expr"]: 1.0 for item in candidates}
    else:
        quality_scores = {item["expr"]: (hi - item["sampled_rmse"]) / (hi - lo) for item in candidates}

    selected = [min(candidates, key=lambda item: (item["sampled_rmse"], item["full_rmse"], item["expr"]))]
    remaining = [item for item in candidates if item["expr"] != selected[0]["expr"]]
    while remaining and len(selected) < beam_width:
        best_item = None
        best_score = None
        for candidate in remaining:
            novelty = min(1.0 - expr_similarity(candidate["expr"], picked["expr"]) for picked in selected)
            score = (1.0 - diversity_weight) * quality_scores[candidate["expr"]] + diversity_weight * novelty
            if best_score is None or score > best_score or (
                math.isclose(score, best_score) and (candidate["sampled_rmse"], candidate["full_rmse"], candidate["expr"]) < (best_item["sampled_rmse"], best_item["full_rmse"], best_item["expr"])
            ):
                best_score = score
                best_item = candidate
        selected.append(best_item)
        remaining = [item for item in remaining if item["expr"] != best_item["expr"]]
    return selected


def primitive_beam_search(bundle, mode: str, max_steps: int, beam_width: int, keep_all_until_step: int, diversity_weight: float, operator: str, constant_variance_threshold: float, improvement_epsilon: float):
    boot = bootstrap_library(bundle, mode, keep_all_until_step, operator, constant_variance_threshold)
    library = boot["library"]
    additions = list(boot["additions"])
    best_by_sample = None
    best_by_full = None
    beam = list(boot["seed_candidates"])
    if not beam:
        beam = [item for item in library if item["expr"] != "x"]
    beam = select_diverse_beam(sorted(beam, key=lambda item: (item["sampled_rmse"], item["full_rmse"], item["expr"])), beam_width, diversity_weight)

    seen_expr = {item["expr"] for item in library}
    seen_sig = {item["signature"] for item in library}
    best_threshold = bundle["residual_fit"]["rmse_sampled"] if mode == "residual" else float("inf")

    for step in range(keep_all_until_step + 1, max_steps + 1):
        generated = []
        for current in beam:
            for combo in iterate_beam_extensions(current, library, operator):
                    expr = format_expr(operator, *(item["expr"] for item in combo))
                    if expr in seen_expr:
                        continue
                    sample_values = apply_operator(operator, *(item["sample_values"] for item in combo))
                    if sample_values is None:
                        continue
                    full_values = apply_operator(operator, *(item["full_values"] for item in combo))
                    if full_values is None:
                        continue
                    sig = value_signature(sample_values)
                    if sig in seen_sig:
                        continue
                    candidate = make_candidate(expr, sample_values, full_values, bundle, mode, step=step)
                    generated.append(candidate)

        if not generated:
            break

        generated = dedupe_candidates(generated)
        generated = filter_nonconstant(generated, constant_variance_threshold)
        generated = [
            candidate
            for candidate in generated
            if candidate["sampled_rmse"] < best_threshold - improvement_epsilon
        ]
        if not generated:
            break
        beam = select_diverse_beam(generated, beam_width, diversity_weight)

        for candidate in beam:
            additions.append(candidate)
            seen_expr.add(candidate["expr"])
            seen_sig.add(candidate["signature"])
            if best_by_sample is None or candidate["sampled_rmse"] < best_by_sample["sampled_rmse"]:
                best_by_sample = candidate
            if best_by_full is None or candidate["full_rmse"] < best_by_full["full_rmse"]:
                best_by_full = candidate
        best_threshold = min(best_threshold, min(candidate["sampled_rmse"] for candidate in beam))

    return {
        "library": library,
        "additions": additions,
        "best_by_sample": best_by_sample,
        "best_by_full": best_by_full,
    }


def discrete_search(bundle, mode: str, max_steps: int, beam_width: int, keep_all_until_step: int, search_strategy: str, diversity_weight: float, operator: str, constant_variance_threshold: float, improvement_epsilon: float):
    if search_strategy == "primitive_beam":
        return primitive_beam_search(bundle, mode, max_steps, beam_width, keep_all_until_step, diversity_weight, operator, constant_variance_threshold, improvement_epsilon)
    return bootstrap_library(bundle, mode, max_steps if max_steps < keep_all_until_step else keep_all_until_step, operator, constant_variance_threshold) if max_steps <= keep_all_until_step else _library_pairwise_search(bundle, mode, max_steps, beam_width, keep_all_until_step, operator, constant_variance_threshold)


def _library_pairwise_search(bundle, mode: str, max_steps: int, beam_width: int, keep_all_until_step: int, operator: str, constant_variance_threshold: float):
    library = initial_library(bundle, mode, operator)
    seen_expr = {item["expr"] for item in library}
    seen_sig = {item["signature"] for item in library}
    additions = []
    best_by_sample = None
    best_by_full = None

    for step in range(1, max_steps + 1):
        generated = []
        for combo in iterate_library_combinations(library, operator):
                expr = format_expr(operator, *(item["expr"] for item in combo))
                if expr in seen_expr:
                    continue
                sample_values = apply_operator(operator, *(item["sample_values"] for item in combo))
                if sample_values is None:
                    continue
                full_values = apply_operator(operator, *(item["full_values"] for item in combo))
                if full_values is None:
                    continue
                sig = value_signature(sample_values)
                if sig in seen_sig:
                    continue
                candidate = make_candidate(expr, sample_values, full_values, bundle, mode, step=step)
                generated.append(candidate)

        if not generated:
            break

        generated = dedupe_candidates(generated)
        generated = filter_nonconstant(generated, constant_variance_threshold)
        chosen = []
        local_seen = set()
        limit = None if step <= keep_all_until_step else beam_width
        for candidate in generated:
            if candidate["signature"] in local_seen:
                continue
            local_seen.add(candidate["signature"])
            chosen.append(candidate)
            if limit is not None and len(chosen) >= limit:
                break

        if not chosen:
            break

        for candidate in chosen:
            library.append(candidate)
            additions.append(candidate)
            seen_expr.add(candidate["expr"])
            seen_sig.add(candidate["signature"])
            if best_by_sample is None or candidate["sampled_rmse"] < best_by_sample["sampled_rmse"]:
                best_by_sample = candidate
            if best_by_full is None or candidate["full_rmse"] < best_by_full["full_rmse"]:
                best_by_full = candidate

    return {
        "library": library,
        "additions": additions,
        "best_by_sample": best_by_sample,
        "best_by_full": best_by_full,
    }


def summarize_steps(additions):
    per_step = {}
    for item in additions:
        per_step.setdefault(item["step"], []).append(item)
    summary = []
    for step in sorted(per_step):
        rows = sorted(per_step[step], key=lambda item: (item["sampled_rmse"], item["full_rmse"], item["expr"]))
        summary.append(
            {
                "step": step,
                "top_candidates": [
                    {
                        "expr": row["expr"],
                        "math": row["math"],
                        "sampled_rmse": row["sampled_rmse"],
                        "full_rmse": row["full_rmse"],
                    }
                    for row in rows
                ],
                "best_sampled_rmse": rows[0]["sampled_rmse"],
                "best_full_rmse": min(row["full_rmse"] for row in rows),
            }
        )
    return summary


def make_plot(outdir: Path, step_summary, baseline_full_rmse):
    if plt is None:
        return
    steps = [item["step"] for item in step_summary]
    sampled = [item["best_sampled_rmse"] for item in step_summary]
    full = [min(row["full_rmse"] for row in item["top_candidates"]) for item in step_summary]

    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    ax.plot(steps, sampled, marker="o", linewidth=1.5, color="#2563eb", label="best sampled RMSE")
    ax.plot(steps, full, marker="s", linewidth=1.5, color="#dc2626", label="best full RMSE")
    ax.axhline(baseline_full_rmse, color="#b45309", linestyle="--", linewidth=1.3, label="Zipf-Mandelbrot full")
    ax.set_xlabel("Discrete search step")
    ax.set_ylabel("RMSE")
    ax.set_title("Discrete EML Search on Zipf Data")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(outdir / "rmse_vs_step.png", dpi=180)
    plt.close(fig)


def sanitize_candidate(item):
    if item is None:
        return None
    return {
        "expr": item["expr"],
        "math": item["math"],
        "step": item["step"],
        "sampled_rmse": item["sampled_rmse"],
        "full_rmse": item["full_rmse"],
    }


def write_report(outdir: Path, args, bundle, search_result, step_summary):
    best_sample = search_result["best_by_sample"]
    best_full = search_result["best_by_full"]
    lines = [
        "# Discrete Operator Search",
        "",
        f"- mode: `{args.mode}`",
        f"- operator: `{args.operator}`",
        f"- sample points: `{len(bundle['sample_indices'])}`",
        f"- beam width: `{args.beam_width}`",
        f"- max steps: `{args.max_steps}`",
        f"- search strategy: `{args.search_strategy}`",
        f"- keep all until step: `{args.keep_all_until_step}`",
        f"- diversity weight: `{args.diversity_weight}`",
        f"- constant variance threshold: `{args.constant_variance_threshold}`",
        f"- improvement epsilon: `{args.improvement_epsilon}`",
        f"- residual baseline: `{args.residual_baseline}`",
        f"- linear full RMSE: `{bundle['linear_fit']['rmse_full']:.6e}`",
        f"- Zipf-Mandelbrot full RMSE: `{bundle['zm_fit']['rmse_full']:.6e}`",
        f"- selected residual baseline full RMSE: `{bundle['residual_fit']['rmse_full']:.6e}`",
        "",
        "## Best By Sampled RMSE",
        "",
    ]

    if best_sample is None:
        lines.extend(
            [
                "- no accepted candidate beat the current best by the improvement threshold.",
                "",
                "## Best By Full RMSE",
                "",
                "- no accepted candidate beat the current best by the improvement threshold.",
                "",
                "## Step Summary",
                "",
            ]
        )
    else:
        lines.extend(
            [
                f"- step: `{best_sample['step']}`",
                f"- formula: `{best_sample['expr']}`",
                f"- math: `{best_sample['math']}`",
                f"- sampled RMSE: `{best_sample['sampled_rmse']:.6e}`",
                f"- full RMSE: `{best_sample['full_rmse']:.6e}`",
                "",
                "## Best By Full RMSE",
                "",
                f"- step: `{best_full['step']}`",
                f"- formula: `{best_full['expr']}`",
                f"- math: `{best_full['math']}`",
                f"- sampled RMSE: `{best_full['sampled_rmse']:.6e}`",
                f"- full RMSE: `{best_full['full_rmse']:.6e}`",
                "",
                "## Step Summary",
                "",
            ]
        )

    for item in step_summary:
        lines.extend(
            [
                f"### Step {item['step']}",
                "",
                f"- best sampled RMSE: `{item['best_sampled_rmse']:.6e}`",
                f"- best full RMSE in step: `{item['best_full_rmse']:.6e}`",
            ]
        )
        for row in item["top_candidates"]:
            lines.extend(
                [
                    f"- `{row['expr']}`",
                    f"  math: `{row['math']}`",
                    f"  sampled/full RMSE: `{row['sampled_rmse']:.6e}` / `{row['full_rmse']:.6e}`",
                ]
            )
        lines.append("")

    (outdir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def write_checkpoint(outdir: Path, payload):
    (outdir / "checkpoint.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main():
    args = parse_args()
    if args.operator == "ternary" and args.x_high >= 1.0:
        args.x_high = TERNARY_SAFE_X_HIGH
    args.outdir.mkdir(parents=True, exist_ok=True)

    bundle = build_bundle(args.corpus_url, args.corpus_path, args.sample_points, args.x_low, args.x_high)
    bundle["residual_fit"] = get_residual_bundle(bundle, args.residual_baseline)
    search_result = discrete_search(
        bundle,
        args.mode,
        args.max_steps,
        args.beam_width,
        args.keep_all_until_step,
        args.search_strategy,
        args.diversity_weight,
        args.operator,
        args.constant_variance_threshold,
        args.improvement_epsilon,
    )
    step_summary = summarize_steps(search_result["additions"])

    payload = {
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "corpus": {
            "corpus_path": str(bundle["corpus_path"]),
            "token_count": bundle["corpus"]["token_count"],
            "unique_words": bundle["corpus"]["unique_words"],
        },
        "linear_zipf": {
            "rmse_full": bundle["linear_fit"]["rmse_full"],
            "rmse_sampled": bundle["linear_fit"]["rmse_sampled"],
        },
        "zipf_mandelbrot": {
            "rmse_full": bundle["zm_fit"]["rmse_full"],
            "rmse_sampled": bundle["zm_fit"]["rmse_sampled"],
        },
        "residual_baseline": {
            "name": bundle["residual_fit"]["name"],
            "rmse_full": bundle["residual_fit"]["rmse_full"],
            "rmse_sampled": bundle["residual_fit"]["rmse_sampled"],
        },
        "best_by_sample": sanitize_candidate(search_result["best_by_sample"]),
        "best_by_full": sanitize_candidate(search_result["best_by_full"]),
        "step_summary": step_summary,
    }
    write_checkpoint(args.outdir, payload)
    (args.outdir / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(args.outdir, args, bundle, search_result, step_summary)
    if not args.skip_plot:
        make_plot(args.outdir, step_summary, bundle["zm_fit"]["rmse_full"])

    print(f"Saved {args.outdir / 'summary.json'}")
    print(f"Saved {args.outdir / 'report.md'}")
    if not args.skip_plot:
        print(f"Saved {args.outdir / 'rmse_vs_step.png'}")


if __name__ == "__main__":
    main()
