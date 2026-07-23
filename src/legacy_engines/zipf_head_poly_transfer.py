import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "results" / "zipf_head_poly_transfer"
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
DECOMP_SUMMARY = ROOT / "results" / "zipf_head_poly_decomposition" / "summary.json"


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


def load_corpus_records():
    decomp = {item["slug"]: item for item in json.loads(DECOMP_SUMMARY.read_text())}
    records = {}
    for entry in CORPORA:
        summary = json.loads(entry["summary_path"].read_text())
        log_rank, log_freq, ranks = build_dataset(entry["corpus_path"])
        x = normalize_x(log_rank)
        zm = summary["zm_baseline"]
        zm_pred = zm["a"] - zm["b"] * np.log(ranks + zm["c"])
        step2 = step2_formula(x)
        poly5 = next(model for model in decomp[entry["slug"]]["poly_models"] if model["degree"] == 5)
        records[entry["slug"]] = {
            "name": entry["name"],
            "slug": entry["slug"],
            "log_freq": log_freq,
            "x": x,
            "zm_pred": zm_pred,
            "zm_rmse": rmse(log_freq, zm_pred),
            "step2_rmse": rmse(log_freq, zm_pred + step2),
            "step10_rmse": summary["zm_search"]["best"]["composite_rmse"],
            "poly5_rmse": poly5["rmse"],
            "poly5_coefficients": poly5["coefficients"],
        }
    return records


def evaluate_transfer(source, target):
    coeffs = np.asarray(source["poly5_coefficients"], dtype=np.float64)
    target_poly = np.polyval(coeffs, target["x"])
    prediction = target["zm_pred"] + step2_formula(target["x"]) + target_poly
    return {
        "source": source["name"],
        "target": target["name"],
        "source_poly5_coefficients": [float(c) for c in coeffs],
        "target_zm_rmse": target["zm_rmse"],
        "target_step2_rmse": target["step2_rmse"],
        "target_in_domain_poly5_rmse": target["poly5_rmse"],
        "target_step10_rmse": target["step10_rmse"],
        "zero_shot_rmse": rmse(target["log_freq"], prediction),
    }


def write_report(records, transfers):
    lines = [
        "# Zipf Head Polynomial Transfer",
        "",
        "- Zero-shot test: fit the degree-5 head polynomial on one corpus, keep those coefficients fixed, and evaluate on the other corpus with that corpus's own ZM baseline plus the universal step-2 term.",
        "",
        "## In-Domain Reference",
        "",
    ]
    for record in records.values():
        lines.extend(
            [
                f"### {record['name']}",
                f"- ZM RMSE: `{record['zm_rmse']:.12f}`",
                f"- Step-2 RMSE: `{record['step2_rmse']:.12f}`",
                f"- Step-2 + in-domain poly5 RMSE: `{record['poly5_rmse']:.12f}`",
                f"- Step-10 monster RMSE: `{record['step10_rmse']:.12f}`",
                "",
            ]
        )

    lines.extend(["## Zero-Shot Transfer", ""])
    for transfer in transfers:
        lines.extend(
            [
                f"### {transfer['source']} poly5 -> {transfer['target']}",
                f"- Target ZM RMSE: `{transfer['target_zm_rmse']:.12f}`",
                f"- Target Step-2 RMSE: `{transfer['target_step2_rmse']:.12f}`",
                f"- Target in-domain poly5 RMSE: `{transfer['target_in_domain_poly5_rmse']:.12f}`",
                f"- Target step-10 RMSE: `{transfer['target_step10_rmse']:.12f}`",
                f"- Zero-shot transferred poly5 RMSE: `{transfer['zero_shot_rmse']:.12f}`",
                "",
            ]
        )

    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")


def sanitize_record(record):
    return {
        "name": record["name"],
        "slug": record["slug"],
        "zm_rmse": record["zm_rmse"],
        "step2_rmse": record["step2_rmse"],
        "step10_rmse": record["step10_rmse"],
        "poly5_rmse": record["poly5_rmse"],
        "poly5_coefficients": record["poly5_coefficients"],
    }


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    records = load_corpus_records()
    transfers = [
        evaluate_transfer(records["shakespeare"], records["war_and_peace"]),
        evaluate_transfer(records["war_and_peace"], records["shakespeare"]),
    ]
    payload = {
        "records": {slug: sanitize_record(record) for slug, record in records.items()},
        "transfers": transfers,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(records, transfers)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
