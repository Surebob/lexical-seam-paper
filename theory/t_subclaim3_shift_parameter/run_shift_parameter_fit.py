from __future__ import annotations

import csv
import importlib.util
import json
import math
import re
from collections import Counter
from pathlib import Path

import jieba
import matplotlib
import numpy as np
import regex

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "t_subclaim3_shift_parameter"
TABLE_1 = ROOT / "experiments" / "1a_per_corpus_enriched_search" / "outputs" / "table1_per_corpus.csv"
COMMON_PATH = ROOT / "zipf_analysis_common.py"
MULTILANG_DATA = ROOT / "data" / "zipf_multilang"
MULTILANG_SUMMARY = ROOT / "results" / "zipf_multilang_romance" / "summary.json"

LOWC_WINNER = "eml[sub[x,1],eml[x,1]]"
HEAD_N = 200
LAMBDA_DERIV_TOL = 1e-6

UNICODE_WORD_RE = regex.compile(r"[\p{L}\p{M}]+(?:['’\-][\p{L}\p{M}]+)*", flags=regex.VERSION1)
HAN_RE = regex.compile(r"\p{Script=Han}", flags=regex.VERSION1)

MULTILANG_TARGETS = {
    "russian_war_and_peace": {
        "language": "Russian",
        "language_family": "Slavic",
        "corpus": "War and Peace (Russian, Wikisource)",
        "tokenizer": "unicode_words",
        "max_tokens": 150000,
    },
    "mandarin_three_kingdoms": {
        "language": "Mandarin",
        "language_family": "Sinitic",
        "corpus": "Romance of the Three Kingdoms (Chinese, Gutenberg 23950)",
        "tokenizer": "jieba",
        "max_tokens": 80000,
    },
    "arabic_1001_nights": {
        "language": "Arabic",
        "language_family": "Semitic",
        "corpus": "One Thousand and One Nights (Arabic, Wikisource)",
        "tokenizer": "unicode_words",
        "max_tokens": 150000,
    },
}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "subclaim3_common")


def safe_slug(name: str) -> str:
    text = name.lower().replace("&", "and")
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def read_table1_lowc_rows() -> list[dict]:
    with TABLE_1.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    lowc = []
    for row in rows:
        if float(row["zm_c"]) < 66.0 and row["step2_winner_expression"] == LOWC_WINNER:
            lowc.append(row)
    return lowc


def tokenize_unicode_words(text: str) -> list[str]:
    return [tok.casefold() for tok in UNICODE_WORD_RE.findall(text)]


def tokenize_jieba_words(text: str) -> list[str]:
    tokens = []
    for token in jieba.cut(text, cut_all=False):
        token = token.strip()
        if not token:
            continue
        if not HAN_RE.search(token):
            continue
        token = regex.sub(r"[^\p{Script=Han}]+", "", token)
        if token:
            tokens.append(token)
    return tokens


def build_dataset_from_tokens(tokens: list[str]) -> dict:
    counts = Counter(tokens)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    return {
        "ranked": ranked,
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": np.log(ranks),
        "log_freq": np.log(freqs),
        "token_count": len(tokens),
        "unique_words": len(freqs),
    }


def load_multilang_dataset(slug: str, spec: dict) -> dict:
    path = MULTILANG_DATA / slug / "combined_clean.txt"
    if not path.exists():
        raise FileNotFoundError(path)
    text = path.read_text(encoding="utf-8", errors="ignore")
    if spec["tokenizer"] == "jieba":
        raw_tokens = tokenize_jieba_words(text)
    else:
        raw_tokens = tokenize_unicode_words(text)
    max_tokens = spec.get("max_tokens")
    tokens = raw_tokens[:max_tokens] if max_tokens else raw_tokens
    return build_dataset_from_tokens(tokens)


def english_dataset_for_name(corpus_name: str) -> dict:
    for spec in common.SEARCHED_CORPORA:
        if spec["name"] == corpus_name:
            return common.build_zipf_dataset(common.corpus_path(spec))
    raise KeyError(f"No SEARCHED_CORPORA match for {corpus_name!r}")


def shifted_xx_sqrt(x: np.ndarray, lam: float) -> np.ndarray:
    return np.power(x, x) - np.sqrt(x) - lam * (x - 1.0)


def fit_lambda(x: np.ndarray, residual: np.ndarray) -> dict:
    base = np.power(x, x) - np.sqrt(x)
    direction = x - 1.0
    denom = float(np.dot(direction, direction))
    lam = float(np.dot(direction, base - residual) / denom)
    fitted = shifted_xx_sqrt(x, lam)
    err = residual - fitted
    sse = float(np.dot(err, err))
    centered = residual - float(np.mean(residual))
    sst = float(np.dot(centered, centered))
    r2 = float(1.0 - sse / sst) if sst > 0 else float("nan")
    rmse = float(math.sqrt(sse / len(residual)))
    return {"lambda": lam, "fitted": fitted, "rmse": rmse, "r2": r2}


def fsecond_min() -> float:
    x = np.linspace(0.05, 1.0, 1000, dtype=np.float64)
    d2 = np.power(x, x) * (np.log(x) + 1.0) ** 2 + np.power(x, x) / x + 0.25 * np.power(x, -1.5)
    return float(np.min(d2))


def bregman_status(lam: float, min_second: float) -> tuple[bool, float, float]:
    f1 = 0.0
    fp1 = 0.5 - lam
    ok = abs(f1) <= 1e-12 and abs(fp1) <= LAMBDA_DERIV_TOL and min_second > 0.0
    return ok, f1, fp1


def evaluate_corpus(
    *,
    slug: str,
    corpus: str,
    group: str,
    language: str,
    language_family: str,
    dataset: dict,
    zm_a: float,
    zm_b: float,
    zm_c: float,
    source: str,
    min_second: float,
) -> dict:
    ranks = dataset["ranks"]
    log_freq = dataset["log_freq"]
    pred = zm_a - zm_b * np.log(ranks + zm_c)
    residual = log_freq - pred
    x = common.normalize_x(dataset["log_rank"], 0.05, 1.0)
    n = min(HEAD_N, len(x))
    fit = fit_lambda(x[:n], residual[:n])
    lam = fit["lambda"]
    ok, f1, fp1 = bregman_status(lam, min_second)
    if abs(lam) > 10:
        raise ValueError(f"{corpus} fitted lambda outside reasonable range: {lam}")
    return {
        "slug": slug,
        "corpus": corpus,
        "group": group,
        "language": language,
        "language_family": language_family,
        "V": int(dataset["unique_words"]),
        "token_count": int(dataset["token_count"]),
        "single_ZM_a": float(zm_a),
        "single_ZM_b": float(zm_b),
        "single_ZM_c": float(zm_c),
        "fitted_λ": lam,
        "lambda_minus_0p5": lam - 0.5,
        "R²": fit["r2"],
        "RMSE_after_shift": fit["rmse"],
        "f_at_1": f1,
        "fprime_at_1": fp1,
        "fsecond_min_on_0p05_1": min_second,
        "bregman_conditions_satisfied": ok,
        "head_points": n,
        "source": source,
    }


def write_csv(rows: list[dict], path: Path) -> None:
    fieldnames = [
        "slug",
        "corpus",
        "group",
        "language",
        "language_family",
        "V",
        "token_count",
        "single_ZM_a",
        "single_ZM_b",
        "single_ZM_c",
        "fitted_λ",
        "lambda_minus_0p5",
        "R²",
        "RMSE_after_shift",
        "f_at_1",
        "fprime_at_1",
        "fsecond_min_on_0p05_1",
        "bregman_conditions_satisfied",
        "head_points",
        "source",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return float("nan")
    x = np.asarray(xs, dtype=np.float64)
    y = np.asarray(ys, dtype=np.float64)
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def make_plots(rows: list[dict]) -> None:
    lambdas = [float(r["fitted_λ"]) for r in rows]
    cs = [float(r["single_ZM_c"]) for r in rows]
    vs = [float(r["V"]) for r in rows]
    colors = {"English": "#3f6fb5", "Russian": "#b5533f", "Mandarin": "#4f9a64", "Arabic": "#8c5fbf"}
    row_colors = [colors.get(r["language"], "#666666") for r in rows]

    plt.figure(figsize=(8, 5))
    plt.hist(lambdas, bins=min(10, max(4, len(rows) // 2)), color="#526a83", edgecolor="white")
    plt.axvline(0.5, color="#b23a48", linestyle="--", label="λ = 0.5")
    plt.xlabel("fitted λ")
    plt.ylabel("corpus count")
    plt.title("Shift-parameter distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTDIR / "shift_lambda_distribution.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.scatter(cs, lambdas, c=row_colors, s=52)
    plt.axhline(0.5, color="#b23a48", linestyle="--", linewidth=1)
    plt.xlabel("single-ZM c")
    plt.ylabel("fitted λ")
    plt.title("λ vs single-ZM c")
    plt.tight_layout()
    plt.savefig(OUTDIR / "lambda_vs_c_scatter.png", dpi=180)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.scatter(vs, lambdas, c=row_colors, s=52)
    plt.axhline(0.5, color="#b23a48", linestyle="--", linewidth=1)
    plt.xlabel("vocabulary size V")
    plt.ylabel("fitted λ")
    plt.title("λ vs vocabulary size")
    plt.tight_layout()
    plt.savefig(OUTDIR / "lambda_vs_V_scatter.png", dpi=180)
    plt.close()

    families = sorted(set(r["language_family"] for r in rows))
    positions = {fam: idx for idx, fam in enumerate(families)}
    plt.figure(figsize=(8, 5))
    xs = [positions[r["language_family"]] for r in rows]
    plt.scatter(xs, lambdas, c=row_colors, s=52)
    plt.axhline(0.5, color="#b23a48", linestyle="--", linewidth=1)
    plt.xticks(range(len(families)), families, rotation=25, ha="right")
    plt.ylabel("fitted λ")
    plt.title("λ by language family")
    plt.tight_layout()
    plt.savefig(OUTDIR / "lambda_by_language_family.png", dpi=180)
    plt.close()


def build_report(rows: list[dict], skipped: list[dict]) -> str:
    lambdas = [float(r["fitted_λ"]) for r in rows]
    english = [r for r in rows if r["group"] == "english_low_c"]
    multilang = [r for r in rows if r["group"] == "multilingual_c0"]
    r_c = pearson([float(r["fitted_λ"]) for r in rows], [float(r["single_ZM_c"]) for r in rows])
    r_v = pearson([float(r["fitted_λ"]) for r in rows], [float(r["V"]) for r in rows])
    bregman_pass = sum(1 for r in rows if str(r["bregman_conditions_satisfied"]) == "True")
    std = float(np.std(lambdas, ddof=1)) if len(lambdas) > 1 else 0.0
    mean = float(np.mean(lambdas)) if lambdas else float("nan")
    median = float(np.median(lambdas)) if lambdas else float("nan")
    min_row = min(rows, key=lambda r: float(r["fitted_λ"]))
    max_row = max(rows, key=lambda r: float(r["fitted_λ"]))

    def group_stats(group_rows: list[dict]) -> dict:
        vals = [float(r["fitted_λ"]) for r in group_rows]
        return {
            "n": len(vals),
            "mean": float(np.mean(vals)) if vals else float("nan"),
            "median": float(np.median(vals)) if vals else float("nan"),
            "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
            "min": float(np.min(vals)) if vals else float("nan"),
            "max": float(np.max(vals)) if vals else float("nan"),
            "r_c": pearson([float(r["fitted_λ"]) for r in group_rows], [float(r["single_ZM_c"]) for r in group_rows]),
            "r_v": pearson([float(r["fitted_λ"]) for r in group_rows], [float(r["V"]) for r in group_rows]),
        }

    english_stats = group_stats(english)
    multilang_stats = group_stats(multilang)

    lines = [
        "# Subclaim 3 Shift-Parameter Fit",
        "",
        "## Protocol",
        "",
        f"Fit `f(x) = x^x - sqrt(x) - λ(x - 1)` to the top-{HEAD_N} single-ZM log-frequency residual with λ as the only free parameter.",
        "The normalized coordinate is the paper's log-rank coordinate over the full vocabulary, `x = 0.05 + 0.95 log(r) / log(V)`.",
        "English targets are the canonical Section 3.1 low-c corpora with `c < 66` and step-2 winner `exp(x - 1) - x`; multilingual targets are Russian, Mandarin, and Arabic from the cached multilingual corpus run.",
        "",
        "## Pre-flight",
        "",
        f"- English low-c targets loaded: `{len(english)}`.",
        f"- Multilingual c≈0 targets loaded: `{len(multilang)}`.",
        f"- Skipped targets: `{len(skipped)}`.",
        f"- Minimum analytic/numerical `f''` on `[0.05, 1]`: `{rows[0]['fsecond_min_on_0p05_1']:.12f}`.",
        "",
        "## Aggregate Results",
        "",
        f"- N fitted corpora: `{len(rows)}`.",
        f"- λ mean / median / std: `{mean:.6f}` / `{median:.6f}` / `{std:.6f}`.",
        f"- λ range: `{float(min_row['fitted_λ']):.6f}` ({min_row['corpus']}) to `{float(max_row['fitted_λ']):.6f}` ({max_row['corpus']}).",
        f"- Pearson r(λ, single-ZM c): `{r_c:.6f}`.",
        f"- Pearson r(λ, V): `{r_v:.6f}`.",
        f"- Bregman boundary conditions satisfied at fitted λ: `{bregman_pass}` / `{len(rows)}`.",
        "",
        "## Group Breakdown",
        "",
        f"- English low-c λ mean / median / std: `{english_stats['mean']:.6f}` / `{english_stats['median']:.6f}` / `{english_stats['std']:.6f}`; range `{english_stats['min']:.6f}` to `{english_stats['max']:.6f}`.",
        f"- English low-c Pearson r(λ, c) / r(λ, V): `{english_stats['r_c']:.6f}` / `{english_stats['r_v']:.6f}`.",
        f"- Multilingual c≈0 λ mean / median / std: `{multilang_stats['mean']:.6f}` / `{multilang_stats['median']:.6f}` / `{multilang_stats['std']:.6f}`; range `{multilang_stats['min']:.6f}` to `{multilang_stats['max']:.6f}`.",
        f"- Multilingual c≈0 Pearson r(λ, V): `{multilang_stats['r_v']:.6f}`; r(λ, c) is undefined because all three fitted c values are zero.",
        "",
        "## Interpretation Against Requested Framework",
        "",
    ]
    if std < 0.1:
        lines.append("The fitted λ values cluster tightly by the requested `std < 0.1` criterion.")
    elif abs(r_c) > 0.5 or abs(r_v) > 0.5:
        lines.append("The fitted λ values vary systematically by the requested correlation criterion (`|r| > 0.5` for c or V).")
    else:
        lines.append("The fitted λ values do not cluster tightly around 0.5 and do not show a strong linear correlation with c or V under the requested thresholds.")
    lines.extend(
        [
            "Because `f'(1) = 0.5 - λ`, exact Bregman recentering requires λ = 0.5; the free residual fits should therefore be read as a data-fit diagnostic, not as automatically Bregman-valid generators.",
            "",
            "## Per-Corpus Summary",
            "",
            "| Corpus | Group | V | c | λ | R² | RMSE | Bregman conditions? |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for r in rows:
        lines.append(
            f"| {r['corpus']} | {r['group']} | {r['V']} | {float(r['single_ZM_c']):.6f} | {float(r['fitted_λ']):.6f} | {float(r['R²']):.6f} | {float(r['RMSE_after_shift']):.6f} | {r['bregman_conditions_satisfied']} |"
        )
    if skipped:
        lines.extend(["", "## Skipped Targets", ""])
        for item in skipped:
            lines.append(f"- `{item['slug']}`: {item['reason']}")
    return "\n".join(lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    min_second = fsecond_min()
    rows: list[dict] = []
    skipped: list[dict] = []

    for table_row in read_table1_lowc_rows():
        corpus = table_row["corpus"]
        dataset = english_dataset_for_name(corpus)
        rows.append(
            evaluate_corpus(
                slug=safe_slug(corpus),
                corpus=corpus,
                group="english_low_c",
                language="English",
                language_family="Germanic/Romance/Other English corpus",
                dataset=dataset,
                zm_a=float(table_row["zm_a"]),
                zm_b=float(table_row["zm_b"]),
                zm_c=float(table_row["zm_c"]),
                source="experiments/1a_per_corpus_enriched_search/outputs/table1_per_corpus.csv",
                min_second=min_second,
            )
        )

    if not MULTILANG_SUMMARY.exists():
        for slug in MULTILANG_TARGETS:
            skipped.append({"slug": slug, "reason": f"missing {MULTILANG_SUMMARY}"})
    else:
        summary = json.loads(MULTILANG_SUMMARY.read_text(encoding="utf-8"))
        summary_rows = {row["slug"]: row for row in summary.get("rows", [])}
        for slug, spec in MULTILANG_TARGETS.items():
            if slug not in summary_rows:
                skipped.append({"slug": slug, "reason": "missing from zipf_multilang_romance summary"})
                continue
            try:
                dataset = load_multilang_dataset(slug, spec)
                zm = common.fit_zipf_mandelbrot(dataset["ranks"], dataset["log_freq"])
                rows.append(
                    evaluate_corpus(
                        slug=slug,
                        corpus=spec["corpus"],
                        group="multilingual_c0",
                        language=spec["language"],
                        language_family=spec["language_family"],
                        dataset=dataset,
                        zm_a=float(zm["a"]),
                        zm_b=float(zm["b"]),
                        zm_c=float(zm["c"]),
                        source=f"{MULTILANG_DATA / slug / 'combined_clean.txt'}; fit recomputed with zipf_analysis_common.fit_zipf_mandelbrot",
                        min_second=min_second,
                    )
                )
            except Exception as exc:
                skipped.append({"slug": slug, "reason": f"{type(exc).__name__}: {exc}"})

    if not rows:
        raise RuntimeError("No target corpora could be fitted.")

    write_csv(rows, OUTDIR / "per_corpus_shift_fit.csv")
    make_plots(rows)
    (OUTDIR / "subclaim3_shift_parameter_report.md").write_text(build_report(rows, skipped), encoding="utf-8")
    print(f"wrote {len(rows)} rows to {OUTDIR / 'per_corpus_shift_fit.csv'}")
    print(f"skipped {len(skipped)} targets")


if __name__ == "__main__":
    main()
