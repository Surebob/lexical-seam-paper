import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_wls_test"
TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = [
    "*** START OF THE PROJECT GUTENBERG EBOOK",
    "*** START OF THIS PROJECT GUTENBERG EBOOK",
]
END_MARKERS = [
    "*** END OF THE PROJECT GUTENBERG EBOOK",
    "*** END OF THIS PROJECT GUTENBERG EBOOK",
]
CORPORA = [
    ("shakespeare", "Shakespeare", ROOT / "data" / "zipf" / "pg100.txt"),
    ("war_and_peace", "War and Peace", ROOT / "data" / "zipf" / "pg2600.txt"),
    ("bible", "King James Bible", ROOT / "data" / "zipf" / "pg10.txt"),
    ("moby_dick", "Moby Dick", ROOT / "data" / "zipf" / "pg2701.txt"),
]
METHODS = [
    ("ols", "OLS"),
    ("freq_weighted", "Frequency-Weighted"),
    ("rank_weighted", "Rank-Weighted"),
]


def strip_gutenberg_boilerplate(text: str) -> str:
    start = 0
    end = len(text)
    for marker in START_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            line_end = text.find("\n", idx)
            start = line_end + 1 if line_end != -1 else idx
            break
    for marker in END_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            end = idx
            break
    return text[start:end]


def build_zipf_dataset(corpus_path: Path):
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = strip_gutenberg_boilerplate(raw_text)
    tokens = TOKEN_RE.findall(clean_text.lower())
    counts = Counter(tokens)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)
    return {
        "token_count": len(tokens),
        "unique_words": len(freqs),
        "ranks": ranks,
        "freqs": freqs,
        "log_rank": log_rank,
        "log_freq": log_freq,
    }


def normalize_x(values: np.ndarray, low: float = 0.05, high: float = 1.0):
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    scaled = (values - vmin) / (vmax - vmin)
    return low + (high - low) * scaled


def step2_formula(x: np.ndarray):
    return (x - 1.0) - np.log(x)


def rmse(y_true: np.ndarray, y_pred: np.ndarray):
    diff = np.asarray(y_true, dtype=np.float64) - np.asarray(y_pred, dtype=np.float64)
    return float(math.sqrt(float(np.mean(diff * diff))))


def weighted_affine_fit(z: np.ndarray, y: np.ndarray, weights: np.ndarray):
    sqrt_w = np.sqrt(weights)
    design = np.column_stack([np.ones_like(z), z])
    weighted_design = design * sqrt_w[:, None]
    weighted_y = y * sqrt_w
    coeffs, _, _, _ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
    pred = design @ coeffs
    weighted_mse = float(np.sum(weights * (pred - y) ** 2) / np.sum(weights))
    return float(coeffs[0]), float(coeffs[1]), pred, weighted_mse


def fit_zipf_mandelbrot_weighted(ranks: np.ndarray, log_freq: np.ndarray, weights: np.ndarray):
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
        intercept, slope, pred, weighted_mse = weighted_affine_fit(z, log_freq, weights)
        if best is None or weighted_mse < best["weighted_mse"]:
            best = {
                "a": intercept,
                "b": -slope,
                "c": float(c),
                "weighted_mse": weighted_mse,
                "prediction": pred,
            }
    best["rmse_unweighted"] = rmse(log_freq, best["prediction"])
    return best


def get_weights(method_slug: str, dataset):
    if method_slug == "ols":
        weights = np.ones_like(dataset["ranks"], dtype=np.float64)
    elif method_slug == "freq_weighted":
        weights = dataset["freqs"].astype(np.float64)
    elif method_slug == "rank_weighted":
        weights = 1.0 / dataset["ranks"].astype(np.float64)
    else:
        raise ValueError(f"Unknown method {method_slug}")
    return weights / float(np.mean(weights))


def run_corpus(slug: str, name: str, corpus_path: Path):
    dataset = build_zipf_dataset(corpus_path)
    x = normalize_x(dataset["log_rank"])
    step2 = step2_formula(x)
    corpus_result = {
        "slug": slug,
        "name": name,
        "corpus_path": str(corpus_path),
        "token_count": dataset["token_count"],
        "unique_words": dataset["unique_words"],
        "methods": [],
    }
    for method_slug, method_name in METHODS:
        weights = get_weights(method_slug, dataset)
        fit = fit_zipf_mandelbrot_weighted(dataset["ranks"], dataset["log_freq"], weights)
        corrected_pred = fit["prediction"] + step2
        corrected_rmse = rmse(dataset["log_freq"], corrected_pred)
        corpus_result["methods"].append(
            {
                "slug": method_slug,
                "name": method_name,
                "a": fit["a"],
                "b": fit["b"],
                "c": fit["c"],
                "weighted_mse": fit["weighted_mse"],
                "zm_rmse_unweighted": fit["rmse_unweighted"],
                "zm_plus_step2_rmse_unweighted": corrected_rmse,
                "step2_delta": corrected_rmse - fit["rmse_unweighted"],
            }
        )
    return corpus_result


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results = [run_corpus(slug, name, path) for slug, name, path in CORPORA]
    (OUTDIR / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

    lines = [
        "# Zipf Weighted-Least-Squares Test",
        "",
        "All fits use the same log-frequency target and are evaluated with the same unweighted RMSE metric for fair comparison.",
        "",
    ]
    csv_lines = [
        "corpus,method,tokens,vocab,a,b,c,zm_rmse_unweighted,zm_plus_step2_rmse_unweighted,step2_delta"
    ]
    for corpus in results:
        lines.extend(
            [
                f"## {corpus['name']}",
                "",
                f"- tokens: `{corpus['token_count']}`",
                f"- unique words: `{corpus['unique_words']}`",
                "",
                "| Method | a | b | c | ZM RMSE | ZM + step2 RMSE | step2 delta |",
                "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ]
        )
        for method in corpus["methods"]:
            lines.append(
                f"| {method['name']} | {method['a']:.12f} | {method['b']:.12f} | {method['c']:.12f} | "
                f"{method['zm_rmse_unweighted']:.12f} | {method['zm_plus_step2_rmse_unweighted']:.12f} | "
                f"{method['step2_delta']:.12f} |"
            )
            csv_lines.append(
                ",".join(
                    [
                        corpus["slug"],
                        method["slug"],
                        str(corpus["token_count"]),
                        str(corpus["unique_words"]),
                        f"{method['a']:.12f}",
                        f"{method['b']:.12f}",
                        f"{method['c']:.12f}",
                        f"{method['zm_rmse_unweighted']:.12f}",
                        f"{method['zm_plus_step2_rmse_unweighted']:.12f}",
                        f"{method['step2_delta']:.12f}",
                    ]
                )
            )
        lines.append("")

    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")
    (OUTDIR / "summary.csv").write_text("\n".join(csv_lines) + "\n", encoding="utf-8")
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")
    print(f"Saved {OUTDIR / 'summary.csv'}")


if __name__ == "__main__":
    main()
