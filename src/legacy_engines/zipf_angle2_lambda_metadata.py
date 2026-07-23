from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
ANGLE6_SUMMARY_PATH = ROOT / "results" / "zipf_angle6_bible_books" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_angle2_lambda_metadata"


# Coarse structural metadata estimates for the 25 English corpora.
# These are intentionally approximate and order-of-magnitude oriented:
# the goal is to test whether lambda_k varies with corpus heterogeneity at all.
METADATA = {
    "shakespeare": {"structure": "single_author_collection", "unit_count": 37, "author_count": 1, "era_span_years": 24},
    "war_and_peace": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 4},
    "moby_dick": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 1},
    "king_james_bible": {"structure": "multi_author_composite", "unit_count": 66, "author_count": 30, "era_span_years": 1000},
    "federalist_papers": {"structure": "multi_author_composite", "unit_count": 85, "author_count": 3, "era_span_years": 1},
    "grimms_fairy_tales": {"structure": "multi_source_collection", "unit_count": 200, "author_count": 50, "era_span_years": 45},
    "don_quixote": {"structure": "single_work", "unit_count": 2, "author_count": 1, "era_span_years": 10},
    "pride_and_prejudice": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 2},
    "canterbury_tales": {"structure": "framed_multi_unit", "unit_count": 24, "author_count": 1, "era_span_years": 13},
    "arabian_nights_vol1": {"structure": "multi_source_collection", "unit_count": 30, "author_count": 10, "era_span_years": 500},
    "aesops_fables": {"structure": "multi_source_collection", "unit_count": 140, "author_count": 10, "era_span_years": 600},
    "complete_sherlock_holmes": {"structure": "single_author_collection", "unit_count": 60, "author_count": 1, "era_span_years": 40},
    "jane_eyre": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 1},
    "dubliners": {"structure": "single_author_collection", "unit_count": 15, "author_count": 1, "era_span_years": 10},
    "the_iliad": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 100},
    "democracy_in_america": {"structure": "single_work", "unit_count": 2, "author_count": 1, "era_span_years": 5},
    "origin_of_species": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 1},
    "wealth_of_nations": {"structure": "single_work", "unit_count": 5, "author_count": 1, "era_span_years": 1},
    "les_miserables": {"structure": "single_work", "unit_count": 5, "author_count": 1, "era_span_years": 17},
    "decline_and_fall_vol1": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 1},
    "emile": {"structure": "single_work", "unit_count": 5, "author_count": 1, "era_span_years": 1},
    "ulysses": {"structure": "single_work", "unit_count": 18, "author_count": 1, "era_span_years": 8},
    "collected_poe": {"structure": "single_author_collection", "unit_count": 70, "author_count": 1, "era_span_years": 20},
    "principia_ethica": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 1},
    "critique_of_pure_reason": {"structure": "single_work", "unit_count": 1, "author_count": 1, "era_span_years": 6},
}

STRUCTURE_ORDER = {
    "single_work": 0,
    "framed_multi_unit": 1,
    "single_author_collection": 2,
    "multi_source_collection": 3,
    "multi_author_composite": 4,
}


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def rankdata(values: list[float]) -> np.ndarray:
    order = np.argsort(values)
    ranks = np.empty(len(values), dtype=np.float64)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = rank
        i = j + 1
    return ranks


def pearsonr(x: list[float], y: list[float]) -> float:
    x_arr = np.asarray(x, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    x_center = x_arr - x_arr.mean()
    y_center = y_arr - y_arr.mean()
    denom = math.sqrt(float(np.dot(x_center, x_center) * np.dot(y_center, y_center)))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(x_center, y_center) / denom)


def spearmanr(x: list[float], y: list[float]) -> float:
    return pearsonr(rankdata(x).tolist(), rankdata(y).tolist())


def load_softk_rows():
    return json.loads(SOFTK_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]


def build_rows():
    rows = []
    for row in load_softk_rows():
        slug = row["slug"]
        meta = METADATA[slug]
        lam = float(row["best_lambda"]["lambda"])
        rows.append(
            {
                "slug": slug,
                "name": row["name"],
                "lambda_k": lam,
                "log10_lambda_k": float(math.log10(lam)),
                "token_count": int(row["token_count"]),
                "vocab_size": int(row["vocab_size"]),
                "structure": meta["structure"],
                "structure_code": STRUCTURE_ORDER[meta["structure"]],
                "unit_count": int(meta["unit_count"]),
                "author_count": int(meta["author_count"]),
                "era_span_years": int(meta["era_span_years"]),
                "log_unit_count": float(math.log1p(meta["unit_count"])),
                "log_author_count": float(math.log1p(meta["author_count"])),
                "log_era_span_years": float(math.log1p(meta["era_span_years"])),
                "heterogeneity_score": float(
                    math.log1p(meta["unit_count"]) + math.log1p(meta["author_count"]) + 0.5 * math.log1p(meta["era_span_years"])
                ),
            }
        )
    return rows


def summarize(rows: list[dict]) -> dict:
    y = [row["log10_lambda_k"] for row in rows]
    metrics = {}
    for key in ["structure_code", "log_unit_count", "log_author_count", "log_era_span_years", "heterogeneity_score", "token_count", "vocab_size"]:
        x = [row[key] for row in rows]
        metrics[key] = {
            "pearson": pearsonr(x, y),
            "spearman": spearmanr(x, y),
        }

    by_structure = {}
    for structure in sorted({row["structure"] for row in rows}, key=lambda name: STRUCTURE_ORDER[name]):
        subset = [row for row in rows if row["structure"] == structure]
        by_structure[structure] = {
            "count": len(subset),
            "median_lambda_k": float(np.median([row["lambda_k"] for row in subset])),
            "median_log10_lambda_k": float(np.median([row["log10_lambda_k"] for row in subset])),
        }

    highest = sorted(rows, key=lambda row: row["lambda_k"], reverse=True)[:5]
    lowest = sorted(rows, key=lambda row: row["lambda_k"])[:5]
    return {
        "rows": rows,
        "correlations": metrics,
        "by_structure": by_structure,
        "highest_lambda": highest,
        "lowest_lambda": lowest,
    }


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "lambda_k",
        "structure",
        "unit_count",
        "author_count",
        "era_span_years",
        "heterogeneity_score",
        "token_count",
        "vocab_size",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fieldnames})


def build_report(summary: dict) -> str:
    lines = [
        "# Angle 2: Structural Predictors of lambda_k",
        "",
        "- Best soft-k `lambda_k` values were taken from the existing soft-k sweep.",
        "- Structural metadata are coarse hand-coded estimates meant to test whether lambda varies with corpus heterogeneity at all.",
        "",
        "## Correlations With log10(lambda_k)",
        "",
    ]
    for key, label in [
        ("structure_code", "structure code"),
        ("log_unit_count", "log unit count"),
        ("log_author_count", "log author count"),
        ("log_era_span_years", "log era span"),
        ("heterogeneity_score", "heterogeneity score"),
        ("token_count", "token count"),
        ("vocab_size", "vocabulary size"),
    ]:
        info = summary["correlations"][key]
        lines.append(f"- {label}: Pearson `{info['pearson']:.6f}`, Spearman `{info['spearman']:.6f}`")
    lines.extend(
        [
            "",
            "## By Structure",
            "",
        ]
    )
    for structure, info in summary["by_structure"].items():
        lines.append(f"- {structure}: n=`{info['count']}`, median lambda_k=`{info['median_lambda_k']:.6g}`")
    lines.extend(
        [
            "",
            "## Highest lambda_k",
            "",
        ]
    )
    for row in summary["highest_lambda"]:
        lines.append(
            f"- {row['name']}: lambda_k=`{row['lambda_k']:.6g}`, structure=`{row['structure']}`, units=`{row['unit_count']}`, authors=`{row['author_count']}`, span=`{row['era_span_years']}`"
        )
    lines.extend(
        [
            "",
            "## Lowest lambda_k",
            "",
        ]
    )
    for row in summary["lowest_lambda"]:
        lines.append(
            f"- {row['name']}: lambda_k=`{row['lambda_k']:.6g}`, structure=`{row['structure']}`, units=`{row['unit_count']}`, authors=`{row['author_count']}`, span=`{row['era_span_years']}`"
        )
    if ANGLE6_SUMMARY_PATH.exists():
        bible = json.loads(ANGLE6_SUMMARY_PATH.read_text(encoding="utf-8"))
        lines.extend(
            [
                "",
                "## Bible Decomposition Context",
                "",
                f"- whole-Bible failure is structurally real: per-book soft-k beats MOE on `{bible['counts']['softk_beats_moe']}` / `{bible['counts']['n_books']}` books, and median per-book soft-k minus MOE is `{bible['medians']['softk_minus_moe']:.12f}`.",
                f"- aggregate per-book held-out avg NLL is `{bible['aggregate']['bookwise_softk_test_avg_nll']:.12f}` vs whole-Bible single-fit `{bible['aggregate']['whole_bible_softk_test_avg_nll']:.12f}`.",
            ]
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    summary = summarize(rows)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, OUTDIR / "lambda_metadata_table.csv")


if __name__ == "__main__":
    main()
