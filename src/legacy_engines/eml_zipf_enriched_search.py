import argparse
import functools
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
ZIPF_EXPERIMENT_PATH = ROOT / "eml_zipf_experiment.py"
DEFAULT_OUTDIR = ROOT / "results" / "zipf_enriched_search"
DEFAULT_CORPUS_URL = "https://www.gutenberg.org/cache/epub/100/pg100.txt"
DEFAULT_CORPUS_PATH = ROOT / "data" / "zipf" / "pg100.txt"
VALUE_EPS = 1e-12
EXP_CLAMP = 30.0
SIG_ROUND = 10
VALUE_ABS_LIMIT = 1.0e6


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


zipf = load_module(ZIPF_EXPERIMENT_PATH, "zipf_enriched_search_zipf")


UNARY_OPS = [
    {"name": "neg", "arity": 1, "math": lambda a: f"(-{a})"},
    {"name": "inv", "arity": 1, "math": lambda a: f"(1/({a}))"},
    {"name": "sqr", "arity": 1, "math": lambda a: f"(({a})^2)"},
    {"name": "sqrt", "arity": 1, "math": lambda a: f"sqrt({a})"},
    {"name": "exp", "arity": 1, "math": lambda a: f"exp({a})"},
    {"name": "log", "arity": 1, "math": lambda a: f"log({a})"},
]

BINARY_OPS = [
    {"name": "eml", "arity": 2, "commutative": False, "math": lambda a, b: f"EML({a},{b})"},
    {"name": "add", "arity": 2, "commutative": True, "math": lambda a, b: f"({a}+{b})"},
    {"name": "sub", "arity": 2, "commutative": False, "math": lambda a, b: f"({a}-{b})"},
    {"name": "mul", "arity": 2, "commutative": True, "math": lambda a, b: f"({a}*{b})"},
    {"name": "div", "arity": 2, "commutative": False, "math": lambda a, b: f"({a}/{b})"},
    {"name": "pow", "arity": 2, "commutative": False, "math": lambda a, b: f"pow({a},{b})"},
]


def parse_args():
    parser = argparse.ArgumentParser(description="Enriched macro search on Zipf residuals")
    parser.add_argument("--corpus-url", type=str, default=DEFAULT_CORPUS_URL)
    parser.add_argument("--corpus-path", type=Path, default=DEFAULT_CORPUS_PATH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--beam-width", type=int, default=50)
    parser.add_argument("--max-steps", type=int, default=10)
    parser.add_argument("--keep-all-until-step", type=int, default=2)
    parser.add_argument("--diversity-weight", type=float, default=0.35)
    parser.add_argument("--x-low", type=float, default=0.05)
    parser.add_argument("--x-high", type=float, default=1.0)
    parser.add_argument("--constant-variance-threshold", type=float, default=1e-10)
    parser.add_argument("--sample-points", type=int, default=0)
    parser.add_argument("--exp-clamp", type=float, default=EXP_CLAMP)
    parser.add_argument("--value-abs-limit", type=float, default=VALUE_ABS_LIMIT)
    return parser.parse_args()


def value_signature(values: np.ndarray):
    rounded = np.round(np.asarray(values, dtype=np.float64), SIG_ROUND)
    return rounded.tobytes()


def rmse(target: np.ndarray, prediction: np.ndarray):
    with np.errstate(all="ignore"):
        mse = float(np.mean((target - prediction) ** 2))
    if not math.isfinite(mse):
        return float("inf")
    return float(math.sqrt(max(mse, 0.0)))


def finite_and_bounded(values: np.ndarray):
    return np.all(np.isfinite(values)) and np.all(np.abs(values) <= VALUE_ABS_LIMIT)


def render_expr(node):
    if isinstance(node, str):
        return node
    op = node[0]
    children = ",".join(render_expr(child) for child in node[1:])
    return f"{op}[{children}]"


def render_math(node):
    if isinstance(node, str):
        return node
    op = node[0]
    if op in {"neg", "inv", "sqr", "sqrt", "exp", "log"}:
        child = render_math(node[1])
        return next(item["math"](child) for item in UNARY_OPS if item["name"] == op)
    left = render_math(node[1])
    right = render_math(node[2])
    return next(item["math"](left, right) for item in BINARY_OPS if item["name"] == op)


@functools.lru_cache(maxsize=None)
def expr_feature_set(expr: str):
    features = set()

    def parse_at(index: int):
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
                    if expr[end] == ",":
                        end += 1
                        continue
                    if expr[end] == "]":
                        end += 1
                        break
                return (token, *children), end
            return token, end
        if expr[index].isdigit():
            end = index
            while end < len(expr) and expr[end].isdigit():
                end += 1
            return expr[index:end], end
        raise ValueError(f"Bad expression token near {expr[index:index+10]!r}")

    def walk(node):
        if isinstance(node, str):
            features.add(node)
            return
        features.add(render_expr(node))
        for child in node[1:]:
            walk(child)

    node, offset = parse_at(0)
    if offset != len(expr):
        raise ValueError(f"Trailing expression content in {expr!r}")
    walk(node)
    return features


def expr_similarity(expr_a: str, expr_b: str):
    fa = expr_feature_set(expr_a)
    fb = expr_feature_set(expr_b)
    union = len(fa | fb)
    if union == 0:
        return 0.0
    return len(fa & fb) / union


def safe_eml(a: np.ndarray, b: np.ndarray):
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)) or np.any(b <= VALUE_EPS):
        return None
    with np.errstate(all="ignore"):
        out = np.exp(np.clip(a, -EXP_CLAMP, EXP_CLAMP)) - np.log(b)
    if not finite_and_bounded(out):
        return None
    return out


def safe_add(a: np.ndarray, b: np.ndarray):
    with np.errstate(all="ignore"):
        out = a + b
    return out if finite_and_bounded(out) else None


def safe_sub(a: np.ndarray, b: np.ndarray):
    with np.errstate(all="ignore"):
        out = a - b
    return out if finite_and_bounded(out) else None


def safe_mul(a: np.ndarray, b: np.ndarray):
    with np.errstate(all="ignore"):
        out = a * b
    return out if finite_and_bounded(out) else None


def safe_div(a: np.ndarray, b: np.ndarray):
    if np.any(np.abs(b) <= VALUE_EPS):
        return None
    with np.errstate(all="ignore"):
        out = a / b
    return out if finite_and_bounded(out) else None


def safe_neg(a: np.ndarray):
    with np.errstate(all="ignore"):
        out = -a
    return out if finite_and_bounded(out) else None


def safe_inv(a: np.ndarray):
    if np.any(np.abs(a) <= VALUE_EPS):
        return None
    with np.errstate(all="ignore"):
        out = 1.0 / a
    return out if finite_and_bounded(out) else None


def safe_sqr(a: np.ndarray):
    with np.errstate(all="ignore"):
        out = a * a
    return out if finite_and_bounded(out) else None


def safe_sqrt(a: np.ndarray):
    if np.any(a < 0.0):
        return None
    with np.errstate(all="ignore"):
        out = np.sqrt(a)
    return out if finite_and_bounded(out) else None


def safe_pow(a: np.ndarray, b: np.ndarray):
    if np.any(a <= VALUE_EPS):
        return None
    with np.errstate(all="ignore"):
        out = np.power(a, b)
    return out if finite_and_bounded(out) else None


def safe_exp(a: np.ndarray):
    with np.errstate(all="ignore"):
        out = np.exp(np.clip(a, -EXP_CLAMP, EXP_CLAMP))
    return out if finite_and_bounded(out) else None


def safe_log(a: np.ndarray):
    if np.any(a <= VALUE_EPS):
        return None
    with np.errstate(all="ignore"):
        out = np.log(a)
    return out if finite_and_bounded(out) else None


UNARY_FUNCS = {
    "neg": safe_neg,
    "inv": safe_inv,
    "sqr": safe_sqr,
    "sqrt": safe_sqrt,
    "exp": safe_exp,
    "log": safe_log,
}

BINARY_FUNCS = {
    "eml": safe_eml,
    "add": safe_add,
    "sub": safe_sub,
    "mul": safe_mul,
    "div": safe_div,
    "pow": safe_pow,
}


def make_candidate(node, values: np.ndarray, target: np.ndarray, step: int):
    if not finite_and_bounded(values):
        return None
    expr = render_expr(node)
    candidate_rmse = rmse(target, values)
    if not math.isfinite(candidate_rmse):
        return None
    return {
        "node": node,
        "expr": expr,
        "math": render_math(node),
        "values": values,
        "signature": value_signature(values),
        "variance": float(np.var(values)),
        "rmse": candidate_rmse,
        "step": step,
    }


def dedupe_candidates(candidates):
    best_by_signature = {}
    for candidate in candidates:
        incumbent = best_by_signature.get(candidate["signature"])
        if incumbent is None or (candidate["rmse"], candidate["expr"]) < (incumbent["rmse"], incumbent["expr"]):
            best_by_signature[candidate["signature"]] = candidate
    return sorted(best_by_signature.values(), key=lambda item: (item["rmse"], item["expr"]))


def filter_candidates(candidates, constant_variance_threshold: float):
    filtered = []
    for candidate in candidates:
        if candidate["variance"] < constant_variance_threshold:
            continue
        filtered.append(candidate)
    return filtered


def select_diverse_beam(candidates, beam_width: int, diversity_weight: float):
    if len(candidates) <= beam_width:
        return candidates

    rmses = np.asarray([item["rmse"] for item in candidates], dtype=np.float64)
    lo = float(rmses.min())
    hi = float(rmses.max())
    if hi <= lo:
        quality_scores = {item["expr"]: 1.0 for item in candidates}
    else:
        quality_scores = {item["expr"]: (hi - item["rmse"]) / (hi - lo) for item in candidates}

    selected = [min(candidates, key=lambda item: (item["rmse"], item["expr"]))]
    remaining = [item for item in candidates if item["expr"] != selected[0]["expr"]]
    while remaining and len(selected) < beam_width:
        best_item = None
        best_score = None
        for candidate in remaining:
            novelty = min(1.0 - expr_similarity(candidate["expr"], chosen["expr"]) for chosen in selected)
            score = (1.0 - diversity_weight) * quality_scores[candidate["expr"]] + diversity_weight * novelty
            if best_score is None or score > best_score or (
                math.isclose(score, best_score) and (candidate["rmse"], candidate["expr"]) < (best_item["rmse"], best_item["expr"])
            ):
                best_score = score
                best_item = candidate
        selected.append(best_item)
        remaining = [item for item in remaining if item["expr"] != best_item["expr"]]
    return selected


def build_target_bundle(corpus_url: str, corpus_path: Path, sample_points: int, x_low: float, x_high: float):
    corpus_path = zipf.ensure_corpus(corpus_url, corpus_path)
    corpus = zipf.build_zipf_dataset(corpus_path)

    if sample_points <= 0 or sample_points >= corpus["unique_words"]:
        indices = np.arange(corpus["unique_words"], dtype=np.int64)
    else:
        indices = zipf.select_log_spaced_indices(corpus["unique_words"], sample_points)

    x_full = zipf.normalize_x(corpus["log_rank"], x_low, x_high)
    x = x_full[indices]
    y = corpus["log_freq"][indices]

    linear = zipf.fit_linear_zipf(corpus["log_rank"], corpus["log_freq"])
    zm = zipf.fit_zipf_mandelbrot(corpus["ranks"], corpus["log_freq"])

    linear_sample_pred = linear["intercept"] + linear["slope"] * corpus["log_rank"][indices]
    zm_sample_pred = zm["a"] - zm["b"] * np.log(corpus["ranks"][indices] + zm["c"])

    return {
        "corpus": corpus,
        "indices": indices,
        "x_full": x_full,
        "x": x,
        "y": y,
        "linear": {
            "intercept": float(linear["intercept"]),
            "slope": float(linear["slope"]),
            "prediction_full": linear["prediction"],
            "prediction_sample": linear_sample_pred,
            "rmse_full": float(linear["rmse"]),
            "rmse_sample": rmse(y, linear_sample_pred),
        },
        "zm": {
            "a": float(zm["a"]),
            "b": float(zm["b"]),
            "c": float(zm["c"]),
            "prediction_full": zm["prediction"],
            "prediction_sample": zm_sample_pred,
            "rmse_full": float(zm["rmse"]),
            "rmse_sample": rmse(y, zm_sample_pred),
        },
    }


def initial_vocabulary(x_values: np.ndarray, target: np.ndarray):
    ones = np.ones_like(x_values, dtype=np.float64)
    return [
        make_candidate("1", ones, target, step=0),
        make_candidate("x", np.asarray(x_values, dtype=np.float64), target, step=0),
    ]


def generate_candidates(vocabulary, target: np.ndarray, step: int):
    generated = []
    for item in vocabulary:
        for op in UNARY_OPS:
            values = UNARY_FUNCS[op["name"]](item["values"])
            if values is None:
                continue
            node = (op["name"], item["node"])
            candidate = make_candidate(node, values, target, step)
            if candidate is not None:
                generated.append(candidate)

    for op in BINARY_OPS:
        for left_index, left in enumerate(vocabulary):
            right_range = range(left_index, len(vocabulary)) if op["commutative"] else range(len(vocabulary))
            for right_index in right_range:
                right = vocabulary[right_index]
                values = BINARY_FUNCS[op["name"]](left["values"], right["values"])
                if values is None:
                    continue
                node = (op["name"], left["node"], right["node"])
                candidate = make_candidate(node, values, target, step)
                if candidate is not None:
                    generated.append(candidate)
                if not op["commutative"] and left_index != right_index:
                    values_rev = BINARY_FUNCS[op["name"]](right["values"], left["values"])
                    if values_rev is None:
                        continue
                    node_rev = (op["name"], right["node"], left["node"])
                    candidate_rev = make_candidate(node_rev, values_rev, target, step)
                    if candidate_rev is not None:
                        generated.append(candidate_rev)

    return generated


def run_search(x_values: np.ndarray, target: np.ndarray, beam_width: int, max_steps: int, keep_all_until_step: int, diversity_weight: float, constant_variance_threshold: float):
    current_vocabulary = initial_vocabulary(x_values, target)
    additions = []
    best = None

    for step in range(1, max_steps + 1):
        generated = generate_candidates(current_vocabulary, target, step)
        generated = dedupe_candidates(generated)
        generated = filter_candidates(generated, constant_variance_threshold)
        if not generated:
            break

        additions.extend(generated)
        current_best = min(generated, key=lambda item: (item["rmse"], item["expr"]))
        if best is None or current_best["rmse"] < best["rmse"]:
            best = current_best

        if step <= keep_all_until_step:
            if step < keep_all_until_step:
                current_vocabulary = current_vocabulary + generated
            else:
                current_vocabulary = select_diverse_beam(generated, beam_width, diversity_weight)
        else:
            current_vocabulary = select_diverse_beam(generated, beam_width, diversity_weight)

    return {
        "best": best,
        "steps": summarize_steps(additions, beam_width),
    }


def summarize_steps(additions, beam_width: int):
    grouped = {}
    for item in additions:
        grouped.setdefault(item["step"], []).append(item)

    summary = []
    for step in sorted(grouped):
        rows = sorted(grouped[step], key=lambda item: (item["rmse"], item["expr"]))
        summary.append(
            {
                "step": step,
                "best_rmse": rows[0]["rmse"],
                "top_candidates": [
                    {
                        "expr": row["expr"],
                        "math": row["math"],
                        "rmse": row["rmse"],
                    }
                    for row in rows[: min(len(rows), beam_width)]
                ],
            }
        )
    return summary


def sanitize_candidate(candidate):
    if candidate is None:
        return None
    payload = {
        "expr": candidate["expr"],
        "math": candidate["math"],
        "residual_rmse": candidate["rmse"],
        "step": candidate["step"],
    }
    if "composite_rmse" in candidate:
        payload["composite_rmse"] = candidate["composite_rmse"]
    return payload


def write_report(outdir: Path, bundle, linear_result, zm_result):
    lines = [
        "# Enriched Zipf Search",
        "",
        f"- corpus tokens: `{bundle['corpus']['token_count']}`",
        f"- unique words: `{bundle['corpus']['unique_words']}`",
        f"- sample points used: `{len(bundle['indices'])}`",
        f"- linear baseline full RMSE: `{bundle['linear']['rmse_full']:.6e}`",
        f"- Zipf-Mandelbrot full RMSE: `{bundle['zm']['rmse_full']:.6e}`",
        f"- binary EML best full RMSE: `2.012251e-01`",
        "",
        "## Linear Residual",
        "",
    ]

    if linear_result["best"] is None:
        lines.append("- no valid enriched candidate found.")
    else:
        lines.extend(
            [
                f"- best step: `{linear_result['best']['step']}`",
                f"- best formula: `{linear_result['best']['expr']}`",
                f"- math: `{linear_result['best']['math']}`",
                f"- composite RMSE: `{linear_result['best']['composite_rmse']:.6e}`",
            ]
        )

    lines.extend(
        [
            "",
            "## Zipf-Mandelbrot Residual",
            "",
        ]
    )

    if zm_result["best"] is None:
        lines.append("- no valid enriched candidate found.")
    else:
        lines.extend(
            [
                f"- best step: `{zm_result['best']['step']}`",
                f"- best formula: `{zm_result['best']['expr']}`",
                f"- math: `{zm_result['best']['math']}`",
                f"- composite RMSE: `{zm_result['best']['composite_rmse']:.6e}`",
            ]
        )

    (outdir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    globals()["EXP_CLAMP"] = float(args.exp_clamp)
    globals()["VALUE_ABS_LIMIT"] = float(args.value_abs_limit)
    args.outdir.mkdir(parents=True, exist_ok=True)

    bundle = build_target_bundle(args.corpus_url, args.corpus_path, args.sample_points, args.x_low, args.x_high)

    linear_target = bundle["y"] - bundle["linear"]["prediction_sample"]
    linear_result = run_search(
        bundle["x"],
        linear_target,
        args.beam_width,
        args.max_steps,
        args.keep_all_until_step,
        args.diversity_weight,
        args.constant_variance_threshold,
    )
    if linear_result["best"] is not None:
        linear_result["best"]["composite_rmse"] = rmse(
            bundle["y"],
            bundle["linear"]["prediction_sample"] + linear_result["best"]["values"],
        )

    zm_target = bundle["y"] - bundle["zm"]["prediction_sample"]
    zm_result = run_search(
        bundle["x"],
        zm_target,
        args.beam_width,
        args.max_steps,
        args.keep_all_until_step,
        args.diversity_weight,
        args.constant_variance_threshold,
    )
    if zm_result["best"] is not None:
        zm_result["best"]["composite_rmse"] = rmse(
            bundle["y"],
            bundle["zm"]["prediction_sample"] + zm_result["best"]["values"],
        )

    payload = {
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "corpus": {
            "token_count": bundle["corpus"]["token_count"],
            "unique_words": bundle["corpus"]["unique_words"],
        },
        "linear_baseline": {
            "intercept": bundle["linear"]["intercept"],
            "slope": bundle["linear"]["slope"],
            "rmse_full": bundle["linear"]["rmse_full"],
            "rmse_sample": bundle["linear"]["rmse_sample"],
        },
        "zm_baseline": {
            "a": bundle["zm"]["a"],
            "b": bundle["zm"]["b"],
            "c": bundle["zm"]["c"],
            "rmse_full": bundle["zm"]["rmse_full"],
            "rmse_sample": bundle["zm"]["rmse_sample"],
        },
        "linear_search": {
            "best": sanitize_candidate(linear_result["best"]),
            "step_summary": linear_result["steps"],
        },
        "zm_search": {
            "best": sanitize_candidate(zm_result["best"]),
            "step_summary": zm_result["steps"],
        },
    }
    (args.outdir / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(args.outdir, bundle, linear_result, zm_result)
    print(f"Saved {args.outdir / 'summary.json'}")
    print(f"Saved {args.outdir / 'report.md'}")


if __name__ == "__main__":
    main()
