from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t, ttest_1samp


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_pos_manual_v2"
SPACY_ALPHA_PATH = ROOT / "results" / "zipf_alpha_fit" / "summary.json"

TARGET_SLUGS = [
    "moby_dick",
    "federalist_papers",
    "wealth_of_nations",
]

SHAKESPEARE_MANUAL = {
    "name": "Shakespeare",
    "slug": "shakespeare",
    "V": 24458.0,
    "k_crossover": 160.0,
    "source": "manual POS v2",
}


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


common = load_module(ROOT / "zipf_analysis_common.py", "zipf_pos_manual_v2_common")


def normalize_token(word: str) -> str:
    return word.lower().replace("'", "")


BASE_CLOSED_WORDS = {
    # Articles / determiners
    "the", "a", "an", "this", "that", "these", "those", "my", "your", "his", "her", "its", "our",
    "their", "which", "what", "some", "any", "no", "every", "all", "each", "both", "th", "whose",
    "another", "anothers", "such", "same", "certain",
    # Pronouns
    "i", "you", "he", "she", "it", "we", "they", "me", "him", "us", "them", "thou", "thee", "ye",
    "myself", "himself", "herself", "itself", "ourselves", "themselves", "yourself", "mine", "thine",
    "none",
    # Prepositions
    "of", "in", "to", "for", "with", "at", "by", "from", "about", "into", "through", "between", "among",
    "before", "after", "during", "under", "over", "above", "below", "upon", "unto", "against", "within",
    "without", "till", "until", "since", "behind", "beside", "besides", "beyond", "across", "toward",
    "towards", "per",
    # Conjunctions
    "and", "but", "or", "nor", "if", "when", "while", "because", "although", "though", "as", "than",
    "whether", "whereas", "therefore", "so", "yet", "ere",
    # Auxiliary / modal verbs
    "is", "am", "are", "was", "were", "be", "been", "being", "have", "has", "had", "do", "does", "did",
    "will", "would", "shall", "should", "can", "could", "may", "might", "must", "ought", "hath", "doth",
    "art", "wilt", "shalt", "dost", "didst", "hast", "canst",
    # Particles / negation / contractions / vocatives
    "not", "s", "d", "ll", "t", "o", "st", "er", "tis", "twas", "nay", "ay", "aye",
}

# Manual extension: "on" is an obvious preposition omitted from the provided base list.
CORPUS_EXTENSIONS = {
    "moby_dick": {"on"},
    "federalist_papers": {"on"},
    "wealth_of_nations": {"on"},
}

# Special context-sensitive override requested by the prompt.
CORPUS_OPEN_OVERRIDES = {
    "moby_dick": {"one"},
    "federalist_papers": {"one"},
    "wealth_of_nations": {"one"},
}


def classify_word(slug: str, word: str) -> tuple[str, str]:
    normalized = normalize_token(word)
    if normalized in CORPUS_OPEN_OVERRIDES.get(slug, set()):
        return "open", "manual_open_override"
    if normalized in CORPUS_EXTENSIONS.get(slug, set()):
        return "closed", "manual_extension"
    if normalized in BASE_CLOSED_WORDS:
        return "closed", "base_closed"
    return "open", "lexical_open"


def analyze_corpus(spec: dict) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    rows = []
    closed_seen = 0
    crossover_rank = None
    sqrt_v = math.sqrt(dataset["unique_words"])

    for rank, (word, freq) in enumerate(dataset["ranked"][:300], start=1):
        family, source = classify_word(spec["slug"], word)
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
                "classification": family,
                "classification_source": source,
                "running_closed_fraction": running_fraction,
            }
        )

    if crossover_rank is None:
        crossover_rank = 300

    alpha = math.log(crossover_rank) / math.log(dataset["unique_words"])
    deviation = (crossover_rank - sqrt_v) / sqrt_v
    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "V": float(dataset["unique_words"]),
        "sqrt_v": float(sqrt_v),
        "k_crossover": float(crossover_rank),
        "alpha": float(alpha),
        "deviation_from_sqrt_v": float(deviation),
        "extensions_used": sorted(CORPUS_EXTENSIONS.get(spec["slug"], set())),
        "open_overrides_used": sorted(CORPUS_OPEN_OVERRIDES.get(spec["slug"], set())),
        "rows": rows,
    }


def write_corpus_csv(result: dict) -> None:
    path = OUTDIR / f"{result['slug']}_top300_manual_pos.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["rank", "word", "freq", "classification", "classification_source", "running_closed_fraction"],
        )
        writer.writeheader()
        writer.writerows(result["rows"])


def forced_scaling_model(v, alpha):
    return alpha * v


def fit_forced_alpha(points: list[dict]) -> dict:
    log_v = np.array([math.log(point["V"]) for point in points], dtype=float)
    log_k = np.array([math.log(point["k_crossover"]) for point in points], dtype=float)
    popt, pcov = curve_fit(forced_scaling_model, log_v, log_k, p0=(0.5,), maxfev=10000)
    alpha = float(popt[0])
    se = float(math.sqrt(pcov[0, 0]))
    df = len(log_v) - 1
    tcrit = float(t.ppf(0.975, df))
    return {
        "alpha": alpha,
        "alpha_se": se,
        "alpha_ci_95": [alpha - tcrit * se, alpha + tcrit * se],
        "df": df,
    }


def build_points(new_results: list[dict]) -> list[dict]:
    return [
        SHAKESPEARE_MANUAL,
        *[
            {
                "slug": result["slug"],
                "name": result["name"],
                "V": result["V"],
                "k_crossover": result["k_crossover"],
                "source": "manual POS explicit list",
            }
            for result in new_results
        ],
    ]


def write_points_csv(points: list[dict]) -> None:
    path = OUTDIR / "manual_alpha_points.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["slug", "name", "source", "V", "k_crossover"])
        writer.writeheader()
        writer.writerows(points)


def build_report(new_results: list[dict], points: list[dict], forced_fit: dict, ttest_result, spaCy_alpha: float) -> str:
    per_alphas = [math.log(point["k_crossover"]) / math.log(point["V"]) for point in points]
    mean_alpha = float(np.mean(per_alphas))
    se_mean = float(np.std(per_alphas, ddof=1) / math.sqrt(len(per_alphas)))
    tcrit = float(t.ppf(0.975, len(per_alphas) - 1))
    mean_ci = [mean_alpha - tcrit * se_mean, mean_alpha + tcrit * se_mean]
    shift_vs_spacy = forced_fit["alpha"] - spaCy_alpha

    lines = [
        "# Manual POS Validation v2 on Three Additional Corpora",
        "",
        "## Corpus Results",
        "",
        "| Corpus | V | sqrt(V) | manual crossover rank | alpha = log(k)/log(V) | deviation from sqrt(V) | extensions | open overrides |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for result in new_results:
        lines.append(
            f"| {result['name']} | {int(result['V'])} | {result['sqrt_v']:.3f} | {int(result['k_crossover'])} | {result['alpha']:.12f} | {result['deviation_from_sqrt_v']:.6%} | {', '.join(result['extensions_used']) or '-'} | {', '.join(result['open_overrides_used']) or '-'} |"
        )
    lines.extend(
        [
            "",
            "Per-corpus top-300 manual classification tables:",
        ]
    )
    for result in new_results:
        lines.append(f"- `{result['slug']}_top300_manual_pos.csv`")
    lines.extend(
        [
            "",
            "## Four-Corpus Manual Alpha Fit",
            "",
            "| Corpus | V | k_crossover | source | alpha_per_corpus |",
            "| --- | ---: | ---: | --- | ---: |",
        ]
    )
    for point in points:
        lines.append(
            f"| {point['name']} | {int(point['V'])} | {int(point['k_crossover'])} | {point['source']} | {math.log(point['k_crossover']) / math.log(point['V']):.12f} |"
        )
    lines.extend(
        [
            "",
            "### Forced scaling law",
            "",
            "- Model: `k = V^alpha`",
            f"- alpha: `{forced_fit['alpha']:.12f}`",
            f"- alpha 95% CI: `[{forced_fit['alpha_ci_95'][0]:.12f}, {forced_fit['alpha_ci_95'][1]:.12f}]`",
            "",
            "### One-sample t-test on per-corpus alpha values vs 0.50",
            "",
            f"- mean per-corpus alpha: `{mean_alpha:.12f}`",
            f"- mean alpha 95% CI: `[{mean_ci[0]:.12f}, {mean_ci[1]:.12f}]`",
            f"- t statistic: `{float(ttest_result.statistic):.12f}`",
            f"- p value: `{float(ttest_result.pvalue):.12f}`",
            "",
            "### Comparison to previous spaCy-based 4-corpus fit",
            "",
            f"- previous spaCy forced-law alpha: `{spaCy_alpha:.12f}`",
            f"- manual forced-law alpha minus spaCy alpha: `{shift_vs_spacy:.12f}`",
            f"- 95% CI {'includes' if forced_fit['alpha_ci_95'][0] <= 0.5 <= forced_fit['alpha_ci_95'][1] else 'excludes'} `0.5`.",
            f"- Manual classification shifts alpha {'toward' if abs(forced_fit['alpha'] - 0.5) < abs(spaCy_alpha - 0.5) else 'away from'} `0.50` relative to the spaCy result.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results = []
    for slug in TARGET_SLUGS:
        spec = common.get_corpus_spec(slug)
        result = analyze_corpus(spec)
        results.append(result)
        write_corpus_csv(result)

    points = build_points(results)
    forced_fit = fit_forced_alpha(points)
    per_alphas = np.array([math.log(point["k_crossover"]) / math.log(point["V"]) for point in points], dtype=float)
    ttest_result = ttest_1samp(per_alphas, popmean=0.5)
    spaCy_alpha = float(json.loads(SPACY_ALPHA_PATH.read_text(encoding="utf-8"))["forced_scaling_fit"]["alpha"])

    summary = {
        "corpus_results": [
            {key: value for key, value in result.items() if key != "rows"}
            for result in results
        ],
        "manual_alpha_points": points,
        "forced_scaling_fit": forced_fit,
        "per_corpus_alpha_test": {
            "mean_alpha": float(np.mean(per_alphas)),
            "t_statistic": float(ttest_result.statistic),
            "p_value": float(ttest_result.pvalue),
            "df": int(len(per_alphas) - 1),
        },
        "spacy_forced_alpha_reference": spaCy_alpha,
    }

    write_points_csv(points)
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(results, points, forced_fit, ttest_result, spaCy_alpha), encoding="utf-8")


if __name__ == "__main__":
    main()
