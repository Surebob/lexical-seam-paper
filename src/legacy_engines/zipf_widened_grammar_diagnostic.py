import json
import math
from pathlib import Path

import numpy as np
from scipy import special

import eml_zipf_enriched_search as enriched
import zipf_analysis_common as common


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_widened_grammar_diagnostic"

BEAM_WIDTH = 50
MAX_STEPS = 10
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
CONSTANT_VARIANCE_THRESHOLD = 1e-10

NEW_UNARY_OP_NAMES = {
    "sin",
    "cos",
    "tan",
    "sinh",
    "cosh",
    "tanh",
    "erf",
    "gamma",
    "besselj",
}

CORPORA = [
    ("shakespeare", "zm_residual"),
    ("shakespeare", "post_is_residual"),
    ("war_and_peace", "zm_residual"),
    ("war_and_peace", "post_is_residual"),
    ("pride_and_prejudice", "zm_residual"),
    ("pride_and_prejudice", "post_is_residual"),
]


def _unary_math_lookup():
    return {
        "neg": lambda a: f"(-{a})",
        "inv": lambda a: f"(1/({a}))",
        "sqr": lambda a: f"(({a})^2)",
        "sqrt": lambda a: f"sqrt({a})",
        "exp": lambda a: f"exp({a})",
        "log": lambda a: f"log({a})",
        "sin": lambda a: f"sin({a})",
        "cos": lambda a: f"cos({a})",
        "tan": lambda a: f"tan({a})",
        "sinh": lambda a: f"sinh({a})",
        "cosh": lambda a: f"cosh({a})",
        "tanh": lambda a: f"tanh({a})",
        "erf": lambda a: f"erf({a})",
        "gamma": lambda a: f"gamma({a})",
        "besselj": lambda a: f"J0({a})",
    }


def _binary_math_lookup():
    return {
        "eml": lambda a, b: f"EML({a},{b})",
        "add": lambda a, b: f"({a}+{b})",
        "sub": lambda a, b: f"({a}-{b})",
        "mul": lambda a, b: f"({a}*{b})",
        "div": lambda a, b: f"({a}/{b})",
        "pow": lambda a, b: f"pow({a},{b})",
    }


def make_safe_unary_numpy(fn):
    def wrapped(a: np.ndarray):
        with np.errstate(all="ignore"):
            out = fn(a)
        return out if enriched.finite_and_bounded(out) else None

    return wrapped


def safe_tan(a: np.ndarray):
    with np.errstate(all="ignore"):
        out = np.tan(a)
    return out if enriched.finite_and_bounded(out) else None


def safe_gamma(a: np.ndarray):
    if np.any(a <= 0.0):
        return None
    with np.errstate(all="ignore"):
        out = special.gamma(a)
    return out if enriched.finite_and_bounded(out) else None


def build_extended_operators():
    unary_math = _unary_math_lookup()
    binary_math = _binary_math_lookup()
    unary_ops = [
        {"name": "neg", "arity": 1, "math": unary_math["neg"]},
        {"name": "inv", "arity": 1, "math": unary_math["inv"]},
        {"name": "sqr", "arity": 1, "math": unary_math["sqr"]},
        {"name": "sqrt", "arity": 1, "math": unary_math["sqrt"]},
        {"name": "exp", "arity": 1, "math": unary_math["exp"]},
        {"name": "log", "arity": 1, "math": unary_math["log"]},
        {"name": "sin", "arity": 1, "math": unary_math["sin"]},
        {"name": "cos", "arity": 1, "math": unary_math["cos"]},
        {"name": "tan", "arity": 1, "math": unary_math["tan"]},
        {"name": "sinh", "arity": 1, "math": unary_math["sinh"]},
        {"name": "cosh", "arity": 1, "math": unary_math["cosh"]},
        {"name": "tanh", "arity": 1, "math": unary_math["tanh"]},
        {"name": "erf", "arity": 1, "math": unary_math["erf"]},
        {"name": "gamma", "arity": 1, "math": unary_math["gamma"]},
        {"name": "besselj", "arity": 1, "math": unary_math["besselj"]},
    ]
    binary_ops = [
        {"name": "eml", "arity": 2, "commutative": False, "math": binary_math["eml"]},
        {"name": "add", "arity": 2, "commutative": True, "math": binary_math["add"]},
        {"name": "sub", "arity": 2, "commutative": False, "math": binary_math["sub"]},
        {"name": "mul", "arity": 2, "commutative": True, "math": binary_math["mul"]},
        {"name": "div", "arity": 2, "commutative": False, "math": binary_math["div"]},
        {"name": "pow", "arity": 2, "commutative": False, "math": binary_math["pow"]},
    ]
    unary_funcs = {
        "neg": enriched.safe_neg,
        "inv": enriched.safe_inv,
        "sqr": enriched.safe_sqr,
        "sqrt": enriched.safe_sqrt,
        "exp": enriched.safe_exp,
        "log": enriched.safe_log,
        "sin": make_safe_unary_numpy(np.sin),
        "cos": make_safe_unary_numpy(np.cos),
        "tan": safe_tan,
        "sinh": make_safe_unary_numpy(np.sinh),
        "cosh": make_safe_unary_numpy(np.cosh),
        "tanh": make_safe_unary_numpy(np.tanh),
        "erf": make_safe_unary_numpy(special.erf),
        "gamma": safe_gamma,
        "besselj": make_safe_unary_numpy(special.j0),
    }
    binary_funcs = {
        "eml": enriched.safe_eml,
        "add": enriched.safe_add,
        "sub": enriched.safe_sub,
        "mul": enriched.safe_mul,
        "div": enriched.safe_div,
        "pow": enriched.safe_pow,
    }
    return unary_ops, binary_ops, unary_funcs, binary_funcs, unary_math, binary_math


def install_extended_grammar():
    unary_ops, binary_ops, unary_funcs, binary_funcs, unary_math, binary_math = build_extended_operators()

    def render_math(node):
        if isinstance(node, str):
            return node
        op = node[0]
        if op in unary_math:
            return unary_math[op](render_math(node[1]))
        if op in binary_math:
            return binary_math[op](render_math(node[1]), render_math(node[2]))
        raise KeyError(f"Unsupported op in render_math: {op}")

    enriched.UNARY_OPS = unary_ops
    enriched.BINARY_OPS = binary_ops
    enriched.UNARY_FUNCS = unary_funcs
    enriched.BINARY_FUNCS = binary_funcs
    enriched.render_math = render_math
    enriched.expr_feature_set.cache_clear()


def uses_new_operator(expr: str) -> bool:
    return any(f"{name}[" in expr for name in NEW_UNARY_OP_NAMES)


def corpus_bundle(slug: str):
    spec = common.get_corpus_spec(slug)
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    summary = common.load_enriched_summary(spec)
    x = common.normalize_x(dataset["log_rank"])
    zm_pred = common.zm_prediction(summary, dataset["ranks"])
    step2 = common.step2_formula(x)
    return {
        "spec": spec,
        "dataset": dataset,
        "summary": summary,
        "x": x,
        "zm_pred": zm_pred,
        "step2_is": step2,
        "orig_step2": common.get_step2_candidate(summary),
    }


def target_for_mode(bundle: dict, mode: str):
    if mode == "zm_residual":
        return bundle["dataset"]["log_freq"] - bundle["zm_pred"]
    if mode == "post_is_residual":
        return bundle["dataset"]["log_freq"] - bundle["zm_pred"] - bundle["step2_is"]
    raise KeyError(mode)


def top10_for_step(step_summary, step: int):
    for row in step_summary:
        if row["step"] == step:
            return row["top_candidates"][:10]
    return []


def analyze_case(slug: str, mode: str):
    bundle = corpus_bundle(slug)
    target = target_for_mode(bundle, mode)
    search_result = enriched.run_search(
        bundle["x"],
        target,
        beam_width=BEAM_WIDTH,
        max_steps=MAX_STEPS,
        keep_all_until_step=KEEP_ALL_UNTIL_STEP,
        diversity_weight=DIVERSITY_WEIGHT,
        constant_variance_threshold=CONSTANT_VARIANCE_THRESHOLD,
    )
    step_rows = []
    any_new_winner = False
    any_new_in_top10 = False
    for step in range(2, MAX_STEPS + 1):
        top10 = top10_for_step(search_result["steps"], step)
        if not top10:
            continue
        winner = top10[0]
        winner_uses_new = uses_new_operator(winner["expr"])
        top10_has_new = any(uses_new_operator(item["expr"]) for item in top10)
        any_new_winner = any_new_winner or winner_uses_new
        any_new_in_top10 = any_new_in_top10 or top10_has_new
        step_rows.append(
            {
                "step": step,
                "winner_expr": winner["expr"],
                "winner_math": winner["math"],
                "winner_rmse": float(winner["rmse"]),
                "winner_uses_new_operator": winner_uses_new,
                "top10_has_new_operator": top10_has_new,
                "top10": [
                    {
                        "rank": idx + 1,
                        "expr": item["expr"],
                        "math": item["math"],
                        "rmse": float(item["rmse"]),
                        "uses_new_operator": uses_new_operator(item["expr"]),
                    }
                    for idx, item in enumerate(top10)
                ],
            }
        )
    return {
        "slug": slug,
        "name": bundle["spec"]["name"],
        "mode": mode,
        "zm_c": float(bundle["summary"]["zm_baseline"]["c"]),
        "original_step2_expr": bundle["orig_step2"]["expr"],
        "original_step2_math": bundle["orig_step2"]["math"],
        "original_step2_rmse": float(bundle["orig_step2"]["rmse"]),
        "widened_step2_expr": step_rows[0]["winner_expr"] if step_rows else None,
        "widened_step2_math": step_rows[0]["winner_math"] if step_rows else None,
        "widened_step2_rmse": step_rows[0]["winner_rmse"] if step_rows else None,
        "widened_step2_matches_original": bool(step_rows and step_rows[0]["winner_expr"] == bundle["orig_step2"]["expr"]),
        "widened_step2_is_bregman": bool(
            step_rows and step_rows[0]["winner_expr"] in {common.STEP2_OURS_EXPR, common.STEPL2_LOWC_EXPR}
        ),
        "any_winner_uses_new_operator": any_new_winner,
        "any_top10_uses_new_operator": any_new_in_top10,
        "step_rows": step_rows,
    }


def write_report(results):
    lines = [
        "# Widened Grammar Diagnostic",
        "",
        "- Deterministic enumerative beam search with beam `50`, steps `10`, `keep_all_until_step=2`, seeds `{x, 1}`, and no continuous coefficient fitting.",
        "- Added unary operators: `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `erf`, `gamma`, `J0`.",
        "- Residual modes: `zm_residual` and `post_is_residual = log(f) - ZM - ((x-1)-log(x))`.",
        "",
    ]
    for result in results:
        lines.extend(
            [
                f"## {result['name']} / {result['mode']}",
                "",
                f"- single-ZM `c`: `{result['zm_c']:.12f}`",
                f"- original step-2 winner: `{result['original_step2_expr']}`",
                f"- original step-2 math: `{result['original_step2_math']}`",
                f"- widened step-2 winner: `{result['widened_step2_expr']}`",
                f"- widened step-2 math: `{result['widened_step2_math']}`",
                f"- widened step-2 matches original: `{result['widened_step2_matches_original']}`",
                f"- widened step-2 is prior Bregman winner: `{result['widened_step2_is_bregman']}`",
                f"- any widened winner uses new operator: `{result['any_winner_uses_new_operator']}`",
                f"- any widened top-10 uses new operator: `{result['any_top10_uses_new_operator']}`",
                "",
            ]
        )
        for row in result["step_rows"]:
            lines.extend(
                [
                    f"### Step {row['step']}",
                    "",
                    f"- winner: `{row['winner_expr']}`",
                    f"- winner math: `{row['winner_math']}`",
                    f"- winner RMSE: `{row['winner_rmse']:.12f}`",
                    f"- winner uses new operator: `{row['winner_uses_new_operator']}`",
                    f"- top-10 has new operator: `{row['top10_has_new_operator']}`",
                    "",
                    "| rank | expr | rmse | uses new op |",
                    "| --- | --- | ---: | --- |",
                ]
            )
            for item in row["top10"]:
                lines.append(
                    f"| {item['rank']} | `{item['expr']}` | {item['rmse']:.12f} | {item['uses_new_operator']} |"
                )
            lines.append("")
    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    install_extended_grammar()
    results = [analyze_case(slug, mode) for slug, mode in CORPORA]
    payload = {
        "beam_width": BEAM_WIDTH,
        "max_steps": MAX_STEPS,
        "keep_all_until_step": KEEP_ALL_UNTIL_STEP,
        "new_unary_ops": sorted(NEW_UNARY_OP_NAMES),
        "results": results,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(results)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
