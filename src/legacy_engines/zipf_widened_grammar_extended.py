from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_widened_grammar_extended"

COMMON_PATH = ROOT / "zipf_analysis_common.py"
WIDENED_PATH = ROOT / "zipf_widened_grammar_diagnostic.py"
MULTILANG_VERIFY_PATH = ROOT / "zipf_multilang_verify.py"

HEAD_K = 200
STEP_LIMIT = 2


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_widened_extended_common")
widen = load_module(WIDENED_PATH, "zipf_widened_extended_widen")
multilang_verify = load_module(MULTILANG_VERIFY_PATH, "zipf_widened_extended_multilang_verify")


def exp_bregman(x: np.ndarray) -> np.ndarray:
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def xpow_minus_sqrt(x: np.ndarray) -> np.ndarray:
    return np.power(x, x) - np.sqrt(x)


def is_bregman(x: np.ndarray) -> np.ndarray:
    return (x - 1.0) - np.log(x)


def weighted_center(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    weights = np.asarray(weights, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    mean = float(np.sum(weights * values) / np.sum(weights))
    return values - mean


def weighted_centered_cosine(a: np.ndarray, b: np.ndarray, weights: np.ndarray) -> float:
    a0 = weighted_center(a, weights)
    b0 = weighted_center(b, weights)
    denom = float(
        math.sqrt(float(np.sum(weights * a0 * a0))) * math.sqrt(float(np.sum(weights * b0 * b0)))
    )
    if denom <= 0.0:
        return float("nan")
    return float(np.sum(weights * a0 * b0) / denom)


def weighted_span_r2(y: np.ndarray, cols: list[np.ndarray], weights: np.ndarray) -> tuple[float, list[float]]:
    y0 = weighted_center(y, weights)
    centered_cols = [weighted_center(col, weights) for col in cols]
    X = np.column_stack(centered_cols)
    sqrt_w = np.sqrt(weights)
    coeffs, _, _, _ = np.linalg.lstsq(sqrt_w[:, None] * X, sqrt_w * y0, rcond=None)
    pred = X @ coeffs
    sse = float(np.sum(weights * (y0 - pred) ** 2))
    sst = float(np.sum(weights * y0 * y0))
    r2 = 1.0 if sst <= 0.0 else 1.0 - sse / sst
    return float(r2), [float(v) for v in coeffs]


def search_with_raw_steps(x_values: np.ndarray, target: np.ndarray) -> dict:
    enriched = widen.enriched
    current_vocabulary = enriched.initial_vocabulary(x_values, target)
    additions = []
    best = None

    for step in range(1, STEP_LIMIT + 1):
        generated = enriched.generate_candidates(current_vocabulary, target, step)
        generated = enriched.dedupe_candidates(generated)
        generated = enriched.filter_candidates(generated, widen.CONSTANT_VARIANCE_THRESHOLD)
        if not generated:
            break

        additions.extend(generated)
        current_best = min(generated, key=lambda item: (item["rmse"], item["expr"]))
        if best is None or current_best["rmse"] < best["rmse"]:
            best = current_best

        if step <= widen.KEEP_ALL_UNTIL_STEP:
            if step < widen.KEEP_ALL_UNTIL_STEP:
                current_vocabulary = current_vocabulary + generated
            else:
                current_vocabulary = enriched.select_diverse_beam(
                    generated,
                    widen.BEAM_WIDTH,
                    widen.DIVERSITY_WEIGHT,
                )
        else:
            current_vocabulary = enriched.select_diverse_beam(
                generated,
                widen.BEAM_WIDTH,
                widen.DIVERSITY_WEIGHT,
            )

    by_step: dict[int, list[dict]] = {}
    for item in additions:
        by_step.setdefault(int(item["step"]), []).append(item)
    for step in list(by_step):
        by_step[step] = sorted(by_step[step], key=lambda item: (item["rmse"], item["expr"]))
    return {"best": best, "by_step": by_step}


def english_lowc_specs() -> list[dict]:
    specs = []
    for spec in common.SEARCHED_CORPORA:
        summary = common.load_enriched_summary(spec)
        if common.get_step2_expr(summary) == common.STEPL2_LOWC_EXPR:
            specs.append(spec)
    return specs


def english_bundle(spec: dict) -> dict:
    summary = common.load_enriched_summary(spec)
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    zm_pred = common.zm_prediction(summary, dataset["ranks"])
    target = dataset["log_freq"] - zm_pred
    original = common.get_step2_candidate(summary)
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "language": "English",
        "dataset": dataset,
        "x": x,
        "target": target,
        "zm_c": float(summary["zm_baseline"]["c"]),
        "original_step2_expr": original["expr"],
        "original_step2_math": original["math"],
        "selection_note": "empirical English low-c family (original step-2 winner is exp-Bregman)",
    }


def multilingual_bundle(slug: str) -> dict:
    spec, dataset, zm_fit, x, target, raw_token_count = multilang_verify.load_dataset_for_slug(slug)
    row = next(item for item in multilang_verify.load_rows() if item["slug"] == slug)
    return {
        "slug": slug,
        "name": row["corpus"],
        "language": row["language"],
        "dataset": dataset,
        "x": x,
        "target": target,
        "zm_c": float(zm_fit["c"]),
        "original_step2_expr": row["step2_winner"],
        "original_step2_math": row["step2_winner"],
        "selection_note": "multilingual c≈0 test case",
        "raw_token_count": int(raw_token_count),
    }


def manifold_metrics(winner_values: np.ndarray, x: np.ndarray, weights: np.ndarray) -> dict:
    exp_vals = exp_bregman(x)
    xpow_vals = xpow_minus_sqrt(x)
    is_vals = is_bregman(x)
    cos_exp = weighted_centered_cosine(winner_values, exp_vals, weights)
    cos_xpow = weighted_centered_cosine(winner_values, xpow_vals, weights)
    cos_is = weighted_centered_cosine(winner_values, is_vals, weights)
    span_r2, coeffs = weighted_span_r2(winner_values, [exp_vals, xpow_vals], weights)
    if span_r2 > 0.9:
        verdict = "yes"
    elif span_r2 > 0.7:
        verdict = "partial"
    else:
        verdict = "no"
    return {
        "cos_vs_exp": float(cos_exp),
        "cos_vs_xx": float(cos_xpow),
        "cos_vs_is": float(cos_is),
        "span_r2": float(span_r2),
        "span_coeffs": coeffs,
        "verdict": verdict,
    }


def analyze_bundle(bundle: dict) -> dict:
    result = search_with_raw_steps(bundle["x"], bundle["target"])
    step2_rows = result["by_step"].get(2, [])
    if not step2_rows:
        raise RuntimeError(f"No step-2 candidates for {bundle['slug']}")
    top10 = step2_rows[:10]
    winner = top10[0]
    head_k = min(HEAD_K, len(bundle["x"]))
    weights = np.asarray(bundle["dataset"]["freqs"][:head_k], dtype=np.float64)
    metrics = manifold_metrics(
        np.asarray(winner["values"][:head_k], dtype=np.float64),
        np.asarray(bundle["x"][:head_k], dtype=np.float64),
        weights,
    )
    widened_is_bregman = winner["expr"] in {common.STEP2_OURS_EXPR, common.STEPL2_LOWC_EXPR}
    return {
        "slug": bundle["slug"],
        "name": bundle["name"],
        "language": bundle["language"],
        "zm_c": float(bundle["zm_c"]),
        "selection_note": bundle["selection_note"],
        "original_step2_expr": bundle["original_step2_expr"],
        "original_step2_math": bundle["original_step2_math"],
        "widened_step2_expr": winner["expr"],
        "widened_step2_math": winner["math"],
        "widened_step2_rmse": float(winner["rmse"]),
        "widened_matches_original": bool(winner["expr"] == bundle["original_step2_expr"]),
        "widened_is_bregman": bool(widened_is_bregman),
        "top10": [
            {
                "rank": index + 1,
                "expr": row["expr"],
                "math": row["math"],
                "rmse": float(row["rmse"]),
                "uses_new_operator": bool(widen.uses_new_operator(row["expr"])),
            }
            for index, row in enumerate(top10)
        ],
        "manifold": metrics,
    }


def build_report(rows: list[dict], selection_note: str) -> str:
    yes_count = sum(row["manifold"]["verdict"] == "yes" for row in rows)
    partial_count = sum(row["manifold"]["verdict"] == "partial" for row in rows)
    no_count = sum(row["manifold"]["verdict"] == "no" for row in rows)
    lines = [
        "# Extended Widened-Grammar Low-c Manifold Check",
        "",
        f"- Selection note: {selection_note}",
        "- Protocol matches the original widened diagnostic through step 2: deterministic enumerative beam search, beam 50, keep-all-until-step-2, seeds `{x, 1}`, no coefficient fitting, extended grammar with trig/hyperbolic/erf/gamma/J0.",
        "- Implementation note: this extension materializes only steps `1..2`, because later steps do not affect the step-2 beam under `keep-all-until-step-2` and the requested outputs are step-2-only.",
        "- This extension uses raw ZM residuals only.",
        "",
        f"- manifold holds across `{yes_count}/{len(rows)}` low-c corpora under the strict `R^2 > 0.9` rule.",
        f"- verdict counts: `yes={yes_count}`, `partial={partial_count}`, `no={no_count}`",
        "",
        "| corpus | ZM c | widened step-2 winner | cos_vs_exp | cos_vs_xx | 2D_span_R^2 | manifold_verdict |",
        "| --- | ---: | --- | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['zm_c']:.6f} | `{row['widened_step2_expr']}` | {row['manifold']['cos_vs_exp']:.6f} | {row['manifold']['cos_vs_xx']:.6f} | {row['manifold']['span_r2']:.6f} | `{row['manifold']['verdict']}` |"
        )

    for row in rows:
        lines.extend(
            [
                "",
                f"## {row['name']}",
                "",
                f"- language: `{row['language']}`",
                f"- ZM c: `{row['zm_c']:.12f}`",
                f"- original-grammar step-2 winner: `{row['original_step2_expr']}`",
                f"- widened step-2 winner: `{row['widened_step2_expr']}`",
                f"- widened step-2 math: `{row['widened_step2_math']}`",
                f"- widened step-2 RMSE: `{row['widened_step2_rmse']:.12f}`",
                f"- widened matches original grammar winner: `{row['widened_matches_original']}`",
                f"- widened winner is Bregman (`IS` or `exp`): `{row['widened_is_bregman']}`",
                f"- weighted centered cosine vs exp-Bregman: `{row['manifold']['cos_vs_exp']:.12f}`",
                f"- weighted centered cosine vs x^x-sqrt(x): `{row['manifold']['cos_vs_xx']:.12f}`",
                f"- weighted centered cosine vs IS-Bregman: `{row['manifold']['cos_vs_is']:.12f}`",
                f"- weighted centered span R^2 in span{{exp, x^x-sqrt(x)}}: `{row['manifold']['span_r2']:.12f}`",
                f"- manifold verdict: `{row['manifold']['verdict']}`",
                "",
                "| rank | expr | math | RMSE | new-op? |",
                "| ---: | --- | --- | ---: | --- |",
            ]
        )
        for cand in row["top10"]:
            lines.append(
                f"| {cand['rank']} | `{cand['expr']}` | `{cand['math']}` | {cand['rmse']:.12f} | {cand['uses_new_operator']} |"
            )
    return "\n".join(lines) + "\n"


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "language",
        "zm_c",
        "original_step2_expr",
        "widened_step2_expr",
        "widened_matches_original",
        "widened_is_bregman",
        "cos_vs_exp",
        "cos_vs_xx",
        "cos_vs_is",
        "span_r2",
        "manifold_verdict",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "language": row["language"],
                    "zm_c": row["zm_c"],
                    "original_step2_expr": row["original_step2_expr"],
                    "widened_step2_expr": row["widened_step2_expr"],
                    "widened_matches_original": int(row["widened_matches_original"]),
                    "widened_is_bregman": int(row["widened_is_bregman"]),
                    "cos_vs_exp": row["manifold"]["cos_vs_exp"],
                    "cos_vs_xx": row["manifold"]["cos_vs_xx"],
                    "cos_vs_is": row["manifold"]["cos_vs_is"],
                    "span_r2": row["manifold"]["span_r2"],
                    "manifold_verdict": row["manifold"]["verdict"],
                }
            )


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    widen.install_extended_grammar()

    english_specs = english_lowc_specs()
    multilingual_slugs = [
        "russian_war_and_peace",
        "mandarin_three_kingdoms",
        "arabic_1001_nights",
    ]

    bundles = [english_bundle(spec) for spec in english_specs]
    bundles.extend(multilingual_bundle(slug) for slug in multilingual_slugs)

    rows = [analyze_bundle(bundle) for bundle in bundles]
    rows = sorted(rows, key=lambda row: (row["language"], row["name"]))

    selection_note = (
        "English side uses the previously established low-c family (original step-2 exp-Bregman winner) "
        "because the literal saved single-ZM c<2.5 threshold only selects Ulysses and does not match the named low-c set."
    )
    summary = {
        "selection_note": selection_note,
        "n_rows": len(rows),
        "yes_count": int(sum(row["manifold"]["verdict"] == "yes" for row in rows)),
        "partial_count": int(sum(row["manifold"]["verdict"] == "partial" for row in rows)),
        "no_count": int(sum(row["manifold"]["verdict"] == "no" for row in rows)),
        "rows": rows,
    }

    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(rows, OUTDIR / "widened_lowc_summary.csv")
    (OUTDIR / "report.md").write_text(build_report(rows, selection_note), encoding="utf-8")


if __name__ == "__main__":
    main()
