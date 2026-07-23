import json
import math
from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
SQRT_V_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_sqrt_v_all_corpora" / "summary.json"
SIX_PARAM_SUMMARY_PATH = ROOT / "results" / "zipf_6param_test" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_bic_comparison"


MODEL_PARAM_COUNTS = {
    "single_zm": 3,
    "piecewise_k500": 6,
    "reranked_8param": 8,
    "reranked_7param_sqrtv": 7,
    "reranked_6param_sqrtv_w12": 6,
}


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def bic_from_rmse(rmse: float, p: int, n: int) -> float:
    mse = float(rmse) ** 2
    return p * math.log(n) + n * math.log(mse)


def build_rows():
    reranked_all = load_json(RERANKED_ALL_SUMMARY_PATH)["rows"]
    sqrt_v_map = {row["slug"]: row for row in load_json(SQRT_V_ALL_SUMMARY_PATH)["rows"]}
    six_param_map = {row["slug"]: row for row in load_json(SIX_PARAM_SUMMARY_PATH)["rows"]}

    rows = []
    winner_counts = {name: 0 for name in MODEL_PARAM_COUNTS}

    for row in reranked_all:
        slug = row["slug"]
        vocab_size = int(row["vocab_size"])
        model_rmses = {
            "single_zm": float(row["single_zm_rmse"]),
            "piecewise_k500": float(row["piecewise_rmse"]),
            "reranked_8param": float(row["reranked_rmse"]),
            "reranked_7param_sqrtv": float(sqrt_v_map[slug]["sqrt_v_rmse"]),
            "reranked_6param_sqrtv_w12": float(six_param_map[slug]["rmse_6param"]),
        }
        model_bics = {
            model_name: bic_from_rmse(rmse, MODEL_PARAM_COUNTS[model_name], vocab_size)
            for model_name, rmse in model_rmses.items()
        }
        best_model = min(model_bics, key=model_bics.get)
        winner_counts[best_model] += 1

        rows.append(
            {
                "slug": slug,
                "name": row["name"],
                "vocab_size": vocab_size,
                "single_zm_bic": model_bics["single_zm"],
                "piecewise_k500_bic": model_bics["piecewise_k500"],
                "reranked_8param_bic": model_bics["reranked_8param"],
                "reranked_7param_sqrtv_bic": model_bics["reranked_7param_sqrtv"],
                "reranked_6param_sqrtv_w12_bic": model_bics["reranked_6param_sqrtv_w12"],
                "best_model": best_model,
            }
        )
    return rows, winner_counts


def write_csv(rows):
    path = OUTDIR / "bic_table.csv"
    headers = [
        "slug",
        "name",
        "vocab_size",
        "single_zm_bic",
        "piecewise_k500_bic",
        "reranked_8param_bic",
        "reranked_7param_sqrtv_bic",
        "reranked_6param_sqrtv_w12_bic",
        "best_model",
    ]
    with path.open("w", encoding="utf-8") as f:
        f.write(",".join(headers) + "\n")
        for row in rows:
            f.write(
                ",".join(
                    [
                        row["slug"],
                        json.dumps(row["name"]),
                        str(row["vocab_size"]),
                        f"{row['single_zm_bic']:.12f}",
                        f"{row['piecewise_k500_bic']:.12f}",
                        f"{row['reranked_8param_bic']:.12f}",
                        f"{row['reranked_7param_sqrtv_bic']:.12f}",
                        f"{row['reranked_6param_sqrtv_w12_bic']:.12f}",
                        row["best_model"],
                    ]
                )
                + "\n"
            )


def build_report(rows, winner_counts):
    model_labels = {
        "single_zm": "Single ZM",
        "piecewise_k500": "Hard two-component K=500",
        "reranked_8param": "Reranked smooth 8-param",
        "reranked_7param_sqrtv": "Reranked smooth k=sqrt(V)",
        "reranked_6param_sqrtv_w12": "Reranked smooth k=sqrt(V), w=1.2",
    }
    lines = [
        "# BIC Comparison Across 25 Corpora",
        "",
        "- BIC formula: `p * ln(n) + n * ln(MSE)`",
        "- `n` is vocabulary size, `MSE = RMSE^2` on the log-frequency fit.",
        "",
        "## Winner Counts",
        "",
    ]
    for model_name, count in winner_counts.items():
        lines.append(f"- {model_labels[model_name]}: `{count}`")
    lines.extend(
        [
            "",
            "| corpus | Single ZM | Hard two-component K=500 | Reranked smooth 8-param | Reranked smooth k=sqrt(V) | Reranked smooth k=sqrt(V), w=1.2 | best BIC |",
            "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['single_zm_bic']:.3f} | {row['piecewise_k500_bic']:.3f} | {row['reranked_8param_bic']:.3f} | {row['reranked_7param_sqrtv_bic']:.3f} | {row['reranked_6param_sqrtv_w12_bic']:.3f} | {model_labels[row['best_model']]} |"
        )
    return "\n".join(lines) + "\n"


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, winner_counts = build_rows()
    rows.sort(key=lambda row: row["slug"])
    summary = {
        "rows": rows,
        "winner_counts": winner_counts,
        "param_counts": MODEL_PARAM_COUNTS,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(rows, winner_counts), encoding="utf-8")
    write_csv(rows)


if __name__ == "__main__":
    main()
