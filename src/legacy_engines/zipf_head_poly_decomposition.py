import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_head_poly_decomposition"
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
    {
        "name": "Shakespeare",
        "slug": "shakespeare",
        "corpus_path": ROOT / "data" / "zipf" / "pg100.txt",
        "summary_path": ROOT / "results" / "zipf_enriched_search_full" / "summary.json",
    },
    {
        "name": "War and Peace",
        "slug": "war_and_peace",
        "corpus_path": ROOT / "data" / "zipf" / "pg2600.txt",
        "summary_path": ROOT / "results" / "zipf_enriched_war_and_peace_full_seq" / "summary.json",
    },
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


def build_dataset(corpus_path: Path):
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = strip_gutenberg_boilerplate(raw_text)
    tokens = TOKEN_RE.findall(clean_text.lower())
    counts = Counter(tokens)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)
    return log_rank, log_freq, ranks


def normalize_x(values: np.ndarray, low: float = 0.05, high: float = 1.0):
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    scaled = (values - vmin) / (vmax - vmin)
    return low + (high - low) * scaled


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    diff = np.asarray(y_true, dtype=np.float64) - np.asarray(y_pred, dtype=np.float64)
    return float(math.sqrt(float(np.mean(diff * diff))))


def step2_formula(x: np.ndarray):
    return (x - 1.0) - np.log(x)


def run_one(entry):
    summary = json.loads(entry["summary_path"].read_text())
    log_rank, log_freq, ranks = build_dataset(entry["corpus_path"])
    x = normalize_x(log_rank)

    zm = summary["zm_baseline"]
    zm_pred = zm["a"] - zm["b"] * np.log(ranks + zm["c"])
    step2 = step2_formula(x)
    base_residual = log_freq - zm_pred - step2

    result = {
        "name": entry["name"],
        "slug": entry["slug"],
        "zm_rmse": rmse(log_freq, zm_pred),
        "step2_rmse": rmse(log_freq, zm_pred + step2),
        "step10_rmse": summary["zm_search"]["best"]["composite_rmse"],
        "poly_models": [],
    }

    for degree in (3, 4, 5):
        coeffs = np.polyfit(x, base_residual, degree)
        poly = np.polyval(coeffs, x)
        composite = zm_pred + step2 + poly
        result["poly_models"].append(
            {
                "degree": degree,
                "rmse": rmse(log_freq, composite),
                "coefficients": [float(c) for c in coeffs],
            }
        )

    return result


def write_report(results):
    lines = [
        "# Zipf Head Polynomial Decomposition",
        "",
        "- Model tested: `ZM + ((x - 1) - log(x)) + poly_d(x)` fit on the full corpus.",
        "- `poly_d(x)` is least-squares fit on the remaining full-corpus residual after ZM plus the fixed step-2 term.",
        "",
    ]

    for result in results:
        lines.extend(
            [
                f"## {result['name']}",
                "",
                f"- ZM alone RMSE: `{result['zm_rmse']:.12f}`",
                f"- Step-2 only RMSE: `{result['step2_rmse']:.12f}`",
                f"- Step-10 monster RMSE: `{result['step10_rmse']:.12f}`",
            ]
        )
        for model in result["poly_models"]:
            coeffs = ", ".join(f"{c:.12f}" for c in model["coefficients"])
            lines.append(f"- Step-2 + poly degree {model['degree']} RMSE: `{model['rmse']:.12f}`")
            lines.append(f"  coefficients: `{coeffs}`")
        lines.append("")

    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results = [run_one(entry) for entry in CORPORA]
    (OUTDIR / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    write_report(results)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
