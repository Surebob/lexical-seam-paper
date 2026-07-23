from __future__ import annotations

import json
import math
from pathlib import Path
import importlib.util

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_multilang_verify"
MULTILANG_PATH = ROOT / "zipf_multilang.py"
ROMANCE_PATH = ROOT / "zipf_multilang_romance.py"
MULTILANG_SUMMARY_PATH = ROOT / "results" / "zipf_multilang_romance" / "summary.json"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


base = load_module(MULTILANG_PATH, "zipf_multilang_verify_base")
romance = load_module(ROMANCE_PATH, "zipf_multilang_verify_romance")


def combined_specs():
    prev_specs = {spec["slug"]: spec for spec in base.CORPORA}
    prev_dirs = {spec["slug"]: base.DATA_DIR / spec["slug"] for spec in base.CORPORA}
    new_specs = {spec["slug"]: spec for spec in romance.CORPORA}
    new_dirs = {spec["slug"]: romance.DATA_DIR / spec["slug"] for spec in romance.CORPORA}
    specs = {}
    specs.update(prev_specs)
    specs.update(new_specs)
    data_dirs = {}
    data_dirs.update(prev_dirs)
    data_dirs.update(new_dirs)
    return specs, data_dirs


def load_rows():
    return json.loads(MULTILANG_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]


def eml_identity_values(x: np.ndarray):
    inner = np.exp(x) - np.log(np.ones_like(x))
    return np.exp(x - 1.0) - np.log(inner)


def low_bregman(x: np.ndarray):
    return np.exp(x - 1.0) - x


def high_bregman(x: np.ndarray):
    return (x - 1.0) - np.log(x)


def euclidean(x: np.ndarray):
    return (1.0 - x) ** 2


def xpow_minus_sqrt(x: np.ndarray):
    return np.power(x, x) - np.sqrt(x)


def formula_rmse(expr: str, x: np.ndarray, target: np.ndarray) -> float:
    table = {
        "eml[sub[x,1],eml[x,1]]": low_bregman(x),
        "sub[sub[x,1],log[x]]": high_bregman(x),
        "mul[sub[1,x],sub[1,x]]": euclidean(x),
        "sub[pow[x,x],sqrt[x]]": xpow_minus_sqrt(x),
    }
    if expr not in table:
        raise KeyError(f"Unsupported expr for direct verification: {expr}")
    return float(base.correct_model.rmse(target, table[expr]))


def load_dataset_for_slug(slug: str):
    specs, data_dirs = combined_specs()
    spec = specs[slug]
    source_dir = data_dirs[slug]
    text = (source_dir / "combined_clean.txt").read_text(encoding="utf-8", errors="ignore")
    raw_tokens = base.tokenize_text(spec, text)
    max_tokens = spec.get("max_tokens")
    tokens = raw_tokens[:max_tokens] if max_tokens else raw_tokens
    dataset = base.build_dataset(tokens)
    zm_fit = base.common.fit_zipf_mandelbrot(dataset["ranks"], dataset["log_freq"])
    x = base.common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    target = dataset["log_freq"] - zm_fit["prediction"]
    return spec, dataset, zm_fit, x, target, len(raw_tokens)


def analyze_corpora(rows):
    analyzed = []
    for row in rows:
        spec, dataset, zm_fit, x, target, raw_token_count = load_dataset_for_slug(row["slug"])
        winner_expr = row["step2_winner"]
        winner_rmse = formula_rmse(winner_expr, x, target)
        high_rmse = formula_rmse("sub[sub[x,1],log[x]]", x, target)
        euclid_rmse = formula_rmse("mul[sub[1,x],sub[1,x]]", x, target)
        low_rmse = formula_rmse("eml[sub[x,1],eml[x,1]]", x, target)
        xpow_rmse = formula_rmse("sub[pow[x,x],sqrt[x]]", x, target)
        analyzed.append(
            {
                "slug": row["slug"],
                "language": row["language"],
                "corpus": row["corpus"],
                "token_count": dataset["token_count"],
                "raw_token_count": raw_token_count,
                "vocab_size": dataset["unique_words"],
                "zm_a": float(zm_fit["a"]),
                "zm_b": float(zm_fit["b"]),
                "zm_c": float(zm_fit["c"]),
                "step2_winner": winner_expr,
                "winner_rmse": winner_rmse,
                "high_bregman_rmse": high_rmse,
                "low_bregman_rmse": low_rmse,
                "euclidean_rmse": euclid_rmse,
                "xpow_rmse": xpow_rmse,
                "gap_to_high_bregman": float(high_rmse - winner_rmse),
                "gap_to_euclidean": float(euclid_rmse - winner_rmse),
            }
        )
    return analyzed


def build_identity_section():
    xs = np.linspace(0.05, 1.0, 10)
    eml_vals = eml_identity_values(xs)
    closed_vals = low_bregman(xs)
    diffs = np.abs(eml_vals - closed_vals)
    rows = []
    for x_val, eml_val, closed_val, diff in zip(xs, eml_vals, closed_vals, diffs):
        rows.append(
            {
                "x": float(x_val),
                "eml_expr": float(eml_val),
                "exp_minus_x": float(closed_val),
                "abs_diff": float(diff),
            }
        )
    return {
        "rows": rows,
        "max_abs_diff": float(np.max(diffs)),
    }


def build_boundary_section():
    grid = np.linspace(0.05, 1.0, 10000)
    second_grid = np.power(grid, grid) * ((np.log(grid) + 1.0) ** 2 + 1.0 / grid) + 1.0 / (4.0 * np.power(grid, 1.5))
    return {
        "analytic": {
            "f_1": 0.0,
            "fprime_1": 0.5,
            "fdoubleprime_1": 2.25,
            "fprime_expr": "x^x * (log(x) + 1) - 1 / (2*sqrt(x))",
            "fdoubleprime_expr": "x^x * ((log(x)+1)^2 + 1/x) + 1/(4*x^(3/2))",
        },
        "numeric": {
            "min_fdoubleprime_on_[0.05,1]": float(np.min(second_grid)),
            "all_positive_on_[0.05,1]": bool(np.all(second_grid > 0.0)),
        },
        "is_bregman_generator_at_x1": False,
        "reason": "Fails the boundary condition f'(1)=0.",
    }


def build_report(analyzed, identity, boundary):
    english_low_cut = 66.0
    romance_slugs = {"latin_gallic_wars", "french_les_miserables", "spanish_don_quixote", "dutch_max_havelaar"}
    romance_rows = [row for row in analyzed if row["slug"] in romance_slugs]
    low_c_count = sum(row["zm_c"] < english_low_cut for row in romance_rows)
    intermediate_count = len(romance_rows) - low_c_count

    lines = [
        "# Multilingual Verification",
        "",
        "## 1. ZM Parameters",
        "",
        "| Language | Corpus | a | b | c | English-range note |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ]
    for row in analyzed:
        note = "within English 5-245 range" if 5.0 <= row["zm_c"] <= 245.0 else "outside English 5-245 range"
        lines.append(
            f"| {row['language']} | {row['corpus']} | {row['zm_a']:.6f} | {row['zm_b']:.6f} | {row['zm_c']:.6f} | {note} |"
        )

    lines.extend(
        [
            "",
            "## 2. Symbolic Identity Check",
            "",
            "`eml[sub[x,1],eml[x,1]] = exp(x-1) - x` because `eml(a,b)=exp(a)-log(b)` and `eml[x,1]=exp(x)`.",
            "",
            "| x | eml[sub[x,1],eml[x,1]] | exp(x-1)-x | abs diff |",
            "| ---: | ---: | ---: | ---: |",
        ]
    )
    for row in identity["rows"]:
        lines.append(
            f"| {row['x']:.6f} | {row['eml_expr']:.12f} | {row['exp_minus_x']:.12f} | {row['abs_diff']:.3e} |"
        )
    lines.extend(["", f"- max absolute difference: `{identity['max_abs_diff']:.3e}`"])

    lines.extend(
        [
            "",
            "## 3. Bregman Boundary Check for x^x - sqrt(x)",
            "",
            "- `f(1) = 0`",
            "- `f'(x) = x^x (log(x)+1) - 1/(2*sqrt(x))`, so `f'(1) = 0.5`",
            "- `f''(x) = x^x ((log(x)+1)^2 + 1/x) + 1/(4*x^(3/2))`, so `f''(1) = 2.25`",
            f"- min `f''(x)` on `[0.05,1]`: `{boundary['numeric']['min_fdoubleprime_on_[0.05,1]']:.6f}`",
            f"- all positive on `[0.05,1]`: `{boundary['numeric']['all_positive_on_[0.05,1]']}`",
            f"- Bregman-generator verdict at `x=1`: `{boundary['is_bregman_generator_at_x1']}` ({boundary['reason']})",
            "",
            "## 4. Intermediate-zone Check for Latin/French/Spanish/Dutch",
            "",
            f"- low-c count (`c < 66`): `{low_c_count}` / `{len(romance_rows)}`",
            f"- not-low count (`c >= 66`): `{intermediate_count}` / `{len(romance_rows)}`",
            "",
            "| Language | Corpus | c | Low-c by English rule? |",
            "| --- | --- | ---: | --- |",
        ]
    )
    for row in romance_rows:
        lines.append(
            f"| {row['language']} | {row['corpus']} | {row['zm_c']:.6f} | {row['zm_c'] < english_low_cut} |"
        )

    lines.extend(
        [
            "",
            "## 5. RMSE Gaps to IS-Bregman and Euclidean",
            "",
            "| Language | Winner | winner RMSE | gap to (x-1)-log(x) | gap to (1-x)^2 |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in analyzed:
        lines.append(
            f"| {row['language']} | {row['step2_winner']} | {row['winner_rmse']:.12f} | {row['gap_to_high_bregman']:.12f} | {row['gap_to_euclidean']:.12f} |"
        )

    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    analyzed = analyze_corpora(rows)
    identity = build_identity_section()
    boundary = build_boundary_section()

    low_c_cut = 66.0
    summary = {
        "corpora": analyzed,
        "identity_check": identity,
        "xpow_minus_sqrt_boundary_check": boundary,
        "latin_french_spanish_dutch_lowc_count": sum(
            row["zm_c"] < low_c_cut
            for row in analyzed
            if row["slug"] in {"latin_gallic_wars", "french_les_miserables", "spanish_don_quixote", "dutch_max_havelaar"}
        ),
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(analyzed, identity, boundary), encoding="utf-8")


if __name__ == "__main__":
    main()
