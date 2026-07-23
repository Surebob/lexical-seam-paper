from __future__ import annotations

import csv
import importlib.util
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import spacy
from scipy.optimize import curve_fit
from scipy.stats import t, ttest_1samp


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
MANUAL_SUMMARY_PATH = ROOT / "results" / "zipf_pos_manual_v2" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_pos_all_corpora"
TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?$")
TOP_N = 500

CLOSED_TAGS = {"ADP", "AUX", "CCONJ", "DET", "PART", "PRON", "SCONJ"}
OPEN_TAGS = {"ADJ", "ADV", "NOUN", "PROPN", "VERB"}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_pos_all_corpora_common")


def load_spacy_model():
    nlp = spacy.load("en_core_web_sm", disable=["parser", "ner", "lemmatizer"])
    nlp.max_length = max(nlp.max_length, 2_000_000)
    return nlp


def iter_text_chunks(text: str, chunk_size: int = 250_000):
    start = 0
    length = len(text)
    while start < length:
        stop = min(start + chunk_size, length)
        if stop < length:
            split = text.rfind("\n", start, stop)
            if split > start:
                stop = split + 1
        yield text[start:stop]
        start = stop


def collect_pos_counts(nlp, text: str, target_words: set[str]) -> dict[str, Counter]:
    counts: dict[str, Counter] = defaultdict(Counter)
    for doc in nlp.pipe(iter_text_chunks(text), batch_size=4):
        for token in doc:
            word = token.text.lower()
            if word in target_words and TOKEN_RE.fullmatch(word):
                counts[word][token.pos_] += 1
    return counts


def classify_word(pos_counter: Counter) -> tuple[str, str]:
    majority_tag = pos_counter.most_common(1)[0][0] if pos_counter else "UNK"
    family = "closed" if majority_tag in CLOSED_TAGS else "open"
    return family, majority_tag


def fit_forced_alpha(rows: list[dict]) -> dict:
    log_v = np.array([row["log_V"] for row in rows], dtype=float)
    log_k = np.array([row["log_k"] for row in rows], dtype=float)

    def forced_scaling_model(x, alpha):
        return alpha * x

    popt, pcov = curve_fit(forced_scaling_model, log_v, log_k, p0=(0.5,), maxfev=10000)
    alpha = float(popt[0])
    se = float(math.sqrt(pcov[0, 0]))
    df = len(rows) - 1
    tcrit = float(t.ppf(0.975, df))
    return {
        "alpha": alpha,
        "alpha_se": se,
        "alpha_ci_95": [alpha - tcrit * se, alpha + tcrit * se],
        "df": df,
    }


def summarize_alpha_distribution(rows: list[dict]) -> dict:
    alphas = np.array([row["alpha_per_corpus"] for row in rows], dtype=float)
    mean_alpha = float(np.mean(alphas))
    median_alpha = float(np.median(alphas))
    se = float(np.std(alphas, ddof=1) / math.sqrt(len(alphas)))
    tcrit = float(t.ppf(0.975, len(alphas) - 1))
    mean_ci = [mean_alpha - tcrit * se, mean_alpha + tcrit * se]
    ttest_result = ttest_1samp(alphas, popmean=0.5)
    outliers = [
        row["name"]
        for row in rows
        if row["alpha_per_corpus"] > 0.60 or row["alpha_per_corpus"] < 0.45
    ]
    return {
        "mean_alpha": mean_alpha,
        "median_alpha": median_alpha,
        "mean_alpha_ci_95": mean_ci,
        "t_statistic": float(ttest_result.statistic),
        "p_value": float(ttest_result.pvalue),
        "df": int(len(alphas) - 1),
        "outliers": outliers,
    }


def manual_reference() -> dict:
    summary = json.loads(MANUAL_SUMMARY_PATH.read_text(encoding="utf-8"))
    forced = summary["forced_scaling_fit"]
    per_corpus = summary["per_corpus_alpha_test"]
    return {
        "forced_alpha": float(forced["alpha"]),
        "forced_alpha_ci_95": [float(v) for v in forced["alpha_ci_95"]],
        "mean_alpha": float(per_corpus["mean_alpha"]),
        "median_alpha": None,
        "p_value": float(per_corpus["p_value"]),
    }


def analyze_corpus(nlp, spec: dict) -> dict:
    corpus_path = common.corpus_path(spec)
    dataset = common.build_zipf_dataset(corpus_path)
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = common.strip_gutenberg_boilerplate(raw_text)
    ranked = dataset["ranked"][:TOP_N]
    target_words = {word for word, _ in ranked}
    pos_counts = collect_pos_counts(nlp, clean_text, target_words)

    rows = []
    closed_seen = 0
    crossover_rank = None
    sqrt_v = math.sqrt(dataset["unique_words"])
    for rank, (word, freq) in enumerate(ranked, start=1):
        family, majority_tag = classify_word(pos_counts[word])
        if family == "closed":
            closed_seen += 1
        running_fraction = closed_seen / rank
        if crossover_rank is None and running_fraction <= 0.5:
            crossover_rank = rank
        rows.append(
            {
                "rank": rank,
                "word": word,
                "freq": int(freq),
                "majority_pos": majority_tag,
                "classification": family,
                "running_closed_fraction": running_fraction,
            }
        )

    censored = False
    if crossover_rank is None:
        crossover_rank = TOP_N
        censored = True

    vocab = int(dataset["unique_words"])
    alpha = math.log(crossover_rank) / math.log(vocab)
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "V": vocab,
        "sqrt_V": float(sqrt_v),
        "k_crossover": int(crossover_rank),
        "alpha_per_corpus": float(alpha),
        "log_V": float(math.log(vocab)),
        "log_k": float(math.log(crossover_rank)),
        "crossover_deviation_from_sqrt_V": float((crossover_rank - sqrt_v) / sqrt_v),
        "censored_at_top_n": censored,
        "rows": rows,
    }


def write_top_words_csv(result: dict) -> None:
    path = OUTDIR / f"{result['slug']}_top500_pos.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "rank",
                "word",
                "freq",
                "majority_pos",
                "classification",
                "running_closed_fraction",
            ],
        )
        writer.writeheader()
        writer.writerows(result["rows"])


def write_points_csv(rows: list[dict]) -> None:
    path = OUTDIR / "pos_all_corpora_points.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "slug",
                "name",
                "V",
                "sqrt_V",
                "k_crossover",
                "alpha_per_corpus",
                "log_V",
                "log_k",
                "crossover_deviation_from_sqrt_V",
                "censored_at_top_n",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_scaling(rows: list[dict], forced_fit: dict) -> None:
    x = np.array([row["log_V"] for row in rows], dtype=float)
    y = np.array([row["log_k"] for row in rows], dtype=float)
    outlier_mask = np.array(
        [row["alpha_per_corpus"] > 0.60 or row["alpha_per_corpus"] < 0.45 for row in rows],
        dtype=bool,
    )
    xline = np.linspace(float(x.min()) - 0.05, float(x.max()) + 0.05, 200)

    plt.figure(figsize=(9.2, 6.6))
    plt.scatter(x[~outlier_mask], y[~outlier_mask], color="#1f77b4", s=42, label="corpora")
    if np.any(outlier_mask):
        plt.scatter(x[outlier_mask], y[outlier_mask], color="#d62728", s=50, label="outliers")
        for row in [row for row in rows if row["alpha_per_corpus"] > 0.60 or row["alpha_per_corpus"] < 0.45]:
            plt.annotate(row["name"], (row["log_V"], row["log_k"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    plt.plot(xline, forced_fit["alpha"] * xline, color="#174a7e", linewidth=2.2, label=f"forced fit: alpha={forced_fit['alpha']:.3f}")
    plt.plot(xline, 0.5 * xline, color="#6a4c93", linewidth=1.8, linestyle=":", label="sqrt(V) reference: alpha=0.5")
    plt.xlabel("log(V)")
    plt.ylabel("log(k crossover)")
    plt.title("All-Corpus POS Crossover Scaling")
    plt.legend(loc="best", fontsize=9)
    plt.tight_layout()
    plt.savefig(OUTDIR / "logk_vs_logv_all_corpora.png", dpi=220)
    plt.close()


def build_report(rows: list[dict], forced_fit: dict, alpha_stats: dict, manual_ref: dict) -> str:
    lines = [
        "# POS Validation Across All 25 Corpora",
        "",
        "| Corpus | V | sqrt(V) | POS crossover rank | alpha = log(k)/log(V) | deviation from sqrt(V) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['V']} | {row['sqrt_V']:.3f} | {row['k_crossover']} | {row['alpha_per_corpus']:.6f} | {row['crossover_deviation_from_sqrt_V']:.6%} |"
        )

    lines.extend(
        [
            "",
            "## All-25 Scaling Fit",
            "",
            f"- Forced `k = V^alpha` fit: `{forced_fit['alpha']:.12f}`",
            f"- Forced-fit 95% CI: `[{forced_fit['alpha_ci_95'][0]:.12f}, {forced_fit['alpha_ci_95'][1]:.12f}]`",
            f"- Mean per-corpus alpha: `{alpha_stats['mean_alpha']:.12f}`",
            f"- Median per-corpus alpha: `{alpha_stats['median_alpha']:.12f}`",
            f"- Mean-alpha 95% CI: `[{alpha_stats['mean_alpha_ci_95'][0]:.12f}, {alpha_stats['mean_alpha_ci_95'][1]:.12f}]`",
            f"- One-sample t-test vs `0.50`: `t={alpha_stats['t_statistic']:.12f}`, `p={alpha_stats['p_value']:.12f}`",
            f"- Forced-fit CI {'includes' if forced_fit['alpha_ci_95'][0] <= 0.5 <= forced_fit['alpha_ci_95'][1] else 'excludes'} `0.5`.",
            "",
            "## Outliers",
            "",
            f"- Outliers (`alpha > 0.60` or `alpha < 0.45`): `{alpha_stats['outliers']}`",
            "",
            "## Comparison To Earlier 4-Corpus Manual Result",
            "",
            f"- Earlier manual forced alpha: `{manual_ref['forced_alpha']:.12f}`",
            f"- Earlier manual 95% CI: `[{manual_ref['forced_alpha_ci_95'][0]:.12f}, {manual_ref['forced_alpha_ci_95'][1]:.12f}]`",
            f"- All-25 minus manual forced alpha: `{forced_fit['alpha'] - manual_ref['forced_alpha']:.12f}`",
            f"- Earlier manual mean alpha: `{manual_ref['mean_alpha']:.12f}`",
            f"- Earlier manual p-value vs `0.50`: `{manual_ref['p_value']:.12f}`",
            "",
            "Per-corpus top-500 POS tables:",
        ]
    )
    for row in rows:
        lines.append(f"- `{row['slug']}_top500_pos.csv`")
    return "\n".join(lines) + "\n"


def write_partial(rows: list[dict], forced_fit: dict | None = None, alpha_stats: dict | None = None) -> None:
    write_points_csv(rows)
    summary = {"rows": rows}
    if forced_fit is not None:
        summary["forced_fit"] = forced_fit
    if alpha_stats is not None:
        summary["alpha_stats"] = alpha_stats
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    nlp = load_spacy_model()
    rows = []
    for spec in common.SEARCHED_CORPORA:
        print(f"[pos-all] analyzing {spec['slug']}")
        result = analyze_corpus(nlp, spec)
        write_top_words_csv(result)
        rows.append({key: value for key, value in result.items() if key != "rows"})
        rows = sorted(rows, key=lambda row: row["slug"])
        write_partial(rows)

    forced_fit = fit_forced_alpha(rows)
    alpha_stats = summarize_alpha_distribution(rows)
    manual_ref = manual_reference()
    plot_scaling(rows, forced_fit)
    summary = {
        "rows": rows,
        "forced_fit": forced_fit,
        "alpha_stats": alpha_stats,
        "manual_reference": manual_ref,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(rows, forced_fit, alpha_stats, manual_ref), encoding="utf-8")


if __name__ == "__main__":
    main()
