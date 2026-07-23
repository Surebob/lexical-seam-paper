import argparse
import importlib.util
import json
import math
import zipfile
from pathlib import Path
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


ROOT = Path("/Volumes/External2TB/emlexperiment")
ENRICHED_PATH = ROOT / "eml_zipf_enriched_search.py"
DEFAULT_DATA_URL = "https://download.geonames.org/export/dump/cities15000.zip"
DEFAULT_ZIP_PATH = ROOT / "data" / "zipf" / "cities15000.zip"
DEFAULT_TXT_PATH = ROOT / "data" / "zipf" / "cities15000.txt"
DEFAULT_OUTDIR = ROOT / "results" / "zipf_enriched_city_populations"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


enriched = load_module(ENRICHED_PATH, "city_population_enriched_search")
zipf = enriched.zipf


def parse_args():
    parser = argparse.ArgumentParser(description="Run enriched search on ranked city populations")
    parser.add_argument("--dataset-url", type=str, default=DEFAULT_DATA_URL)
    parser.add_argument("--zip-path", type=Path, default=DEFAULT_ZIP_PATH)
    parser.add_argument("--txt-path", type=Path, default=DEFAULT_TXT_PATH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--beam-width", type=int, default=50)
    parser.add_argument("--max-steps", type=int, default=10)
    parser.add_argument("--keep-all-until-step", type=int, default=2)
    parser.add_argument("--diversity-weight", type=float, default=0.35)
    parser.add_argument("--x-low", type=float, default=0.05)
    parser.add_argument("--x-high", type=float, default=1.0)
    parser.add_argument("--constant-variance-threshold", type=float, default=1e-10)
    parser.add_argument("--sample-points", type=int, default=0)
    parser.add_argument("--exp-clamp", type=float, default=enriched.EXP_CLAMP)
    parser.add_argument("--value-abs-limit", type=float, default=enriched.VALUE_ABS_LIMIT)
    return parser.parse_args()


def ensure_city_dataset(dataset_url: str, zip_path: Path, txt_path: Path):
    if txt_path.exists():
        return txt_path
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    if not zip_path.exists():
        with urlopen(dataset_url, timeout=60) as response:
            zip_path.write_bytes(response.read())
    with zipfile.ZipFile(zip_path) as zf:
        names = [name for name in zf.namelist() if name.endswith(".txt")]
        if not names:
            raise RuntimeError(f"no txt payload found in {zip_path}")
        payload = names[0]
        txt_path.write_bytes(zf.read(payload))
    return txt_path


def load_city_dataset(txt_path: Path):
    cities = []
    for line in txt_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.split("\t")
        if len(parts) < 15:
            continue
        try:
            population = int(parts[14])
        except ValueError:
            continue
        if population <= 0:
            continue
        cities.append(
            {
                "name": parts[1],
                "country_code": parts[8],
                "population": population,
            }
        )
    cities.sort(key=lambda row: (-row["population"], row["name"], row["country_code"]))
    populations = np.array([row["population"] for row in cities], dtype=np.float64)
    ranks = np.arange(1, len(populations) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_pop = np.log(populations)
    probabilities = populations / float(np.sum(populations))
    shannon_entropy = float(-np.sum(probabilities * np.log(probabilities)))
    return {
        "rows": cities,
        "populations": populations,
        "ranks": ranks,
        "log_rank": log_rank,
        "log_pop": log_pop,
        "city_count": len(cities),
        "total_population": int(np.sum(populations)),
        "shannon_entropy": shannon_entropy,
        "top_cities": cities[:20],
    }


def fit_zipf_mandelbrot_curve(ranks: np.ndarray, log_values: np.ndarray):
    coarse = zipf.fit_zipf_mandelbrot(ranks, log_values)

    def model(rank, a, b, c):
        return a - b * np.log(rank + c)

    try:
        params, _ = curve_fit(
            model,
            ranks,
            log_values,
            p0=(coarse["a"], max(coarse["b"], 1e-8), max(coarse["c"], 1e-8)),
            bounds=([-np.inf, 0.0, 0.0], [np.inf, np.inf, ranks[-1] * 100.0]),
            maxfev=200000,
        )
        prediction = model(ranks, *params)
        mse = float(np.mean((prediction - log_values) ** 2))
        return {
            "a": float(params[0]),
            "b": float(params[1]),
            "c": float(params[2]),
            "mse": mse,
            "rmse": float(math.sqrt(max(mse, 0.0))),
            "prediction": prediction,
        }
    except Exception:
        return coarse


def build_city_bundle(dataset_url: str, zip_path: Path, txt_path: Path, sample_points: int, x_low: float, x_high: float):
    txt_path = ensure_city_dataset(dataset_url, zip_path, txt_path)
    dataset = load_city_dataset(txt_path)

    if sample_points <= 0 or sample_points >= dataset["city_count"]:
        indices = np.arange(dataset["city_count"], dtype=np.int64)
    else:
        indices = zipf.select_log_spaced_indices(dataset["city_count"], sample_points)

    x_full = zipf.normalize_x(dataset["log_rank"], x_low, x_high)
    x = x_full[indices]
    y = dataset["log_pop"][indices]

    linear = zipf.fit_linear_zipf(dataset["log_rank"], dataset["log_pop"])
    zm = fit_zipf_mandelbrot_curve(dataset["ranks"], dataset["log_pop"])

    linear_sample_pred = linear["intercept"] + linear["slope"] * dataset["log_rank"][indices]
    zm_sample_pred = zm["a"] - zm["b"] * np.log(dataset["ranks"][indices] + zm["c"])

    return {
        "dataset_path": str(txt_path),
        "corpus": {
            "token_count": dataset["total_population"],
            "unique_words": dataset["city_count"],
            "entropy": dataset["shannon_entropy"],
            "top_cities": dataset["top_cities"],
        },
        "rows": dataset["rows"],
        "ranks": dataset["ranks"],
        "log_rank": dataset["log_rank"],
        "log_pop": dataset["log_pop"],
        "indices": indices,
        "x_full": x_full,
        "x": x,
        "y": y,
        "linear": {
            "intercept": float(linear["intercept"]),
            "slope": float(linear["slope"]),
            "prediction_full": linear["prediction"],
            "prediction_sample": linear_sample_pred,
            "rmse_full": float(linear["rmse"]),
            "rmse_sample": enriched.rmse(y, linear_sample_pred),
        },
        "zm": {
            "a": float(zm["a"]),
            "b": float(zm["b"]),
            "c": float(zm["c"]),
            "prediction_full": zm["prediction"],
            "prediction_sample": zm_sample_pred,
            "rmse_full": float(zm["rmse"]),
            "rmse_sample": enriched.rmse(y, zm_sample_pred),
        },
    }


def plot_city_distribution(outdir: Path, bundle):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(bundle["log_rank"], bundle["log_pop"], s=6, alpha=0.25, color="#264653", label="cities")
    ax.plot(bundle["log_rank"], bundle["linear"]["prediction_full"], color="#e76f51", linewidth=2.0, label="linear")
    ax.plot(bundle["log_rank"], bundle["zm"]["prediction_full"], color="#2a9d8f", linewidth=2.0, label="Zipf-Mandelbrot")
    ax.set_xlabel("log(rank)")
    ax.set_ylabel("log(population)")
    ax.set_title("City populations in log-log space")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "city_loglog.png", dpi=220)
    plt.close(fig)


def write_report(outdir: Path, bundle, linear_result, zm_result):
    step2_candidates = zm_result["steps"][1]["top_candidates"][:5] if len(zm_result["steps"]) > 1 else []
    lines = [
        "# Enriched Search on City Populations",
        "",
        f"- dataset path: `{bundle['dataset_path']}`",
        f"- total population (token equivalent): `{bundle['corpus']['token_count']}`",
        f"- city count (vocabulary equivalent): `{bundle['corpus']['unique_words']}`",
        f"- Shannon entropy of population shares: `{bundle['corpus']['entropy']:.6f}`",
        f"- sample points used: `{len(bundle['indices'])}`",
        f"- linear baseline full RMSE: `{bundle['linear']['rmse_full']:.6e}`",
        f"- Zipf-Mandelbrot full RMSE: `{bundle['zm']['rmse_full']:.6e}`",
        f"- fitted ZM parameters: `a={bundle['zm']['a']:.6f}`, `b={bundle['zm']['b']:.6f}`, `c={bundle['zm']['c']:.6f}`",
        "",
        "## Step-2 ZM Candidates",
        "",
    ]
    if not step2_candidates:
        lines.append("- no valid step-2 candidates found.")
    else:
        for idx, row in enumerate(step2_candidates, 1):
            lines.append(f"{idx}. `{row['expr']}` | RMSE `{row['rmse']:.12f}`")

    lines.extend(["", "## Linear Residual", ""])
    if linear_result["best"] is None:
        lines.append("- no valid enriched candidate found.")
    else:
        lines.extend(
            [
                f"- best step: `{linear_result['best']['step']}`",
                f"- best formula: `{linear_result['best']['expr']}`",
                f"- math: `{linear_result['best']['math']}`",
                f"- composite RMSE: `{linear_result['best']['composite_rmse']:.6e}`",
            ]
        )

    lines.extend(["", "## Zipf-Mandelbrot Residual", ""])
    if zm_result["best"] is None:
        lines.append("- no valid enriched candidate found.")
    else:
        lines.extend(
            [
                f"- best step: `{zm_result['best']['step']}`",
                f"- best formula: `{zm_result['best']['expr']}`",
                f"- math: `{zm_result['best']['math']}`",
                f"- composite RMSE: `{zm_result['best']['composite_rmse']:.6e}`",
            ]
        )
    (outdir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    args = parse_args()
    enriched.EXP_CLAMP = float(args.exp_clamp)
    enriched.VALUE_ABS_LIMIT = float(args.value_abs_limit)
    args.outdir.mkdir(parents=True, exist_ok=True)

    bundle = build_city_bundle(args.dataset_url, args.zip_path, args.txt_path, args.sample_points, args.x_low, args.x_high)

    linear_target = bundle["y"] - bundle["linear"]["prediction_sample"]
    linear_result = enriched.run_search(
        bundle["x"],
        linear_target,
        args.beam_width,
        args.max_steps,
        args.keep_all_until_step,
        args.diversity_weight,
        args.constant_variance_threshold,
    )
    if linear_result["best"] is not None:
        linear_result["best"]["composite_rmse"] = enriched.rmse(
            bundle["y"],
            bundle["linear"]["prediction_sample"] + linear_result["best"]["values"],
        )

    zm_target = bundle["y"] - bundle["zm"]["prediction_sample"]
    zm_result = enriched.run_search(
        bundle["x"],
        zm_target,
        args.beam_width,
        args.max_steps,
        args.keep_all_until_step,
        args.diversity_weight,
        args.constant_variance_threshold,
    )
    if zm_result["best"] is not None:
        zm_result["best"]["composite_rmse"] = enriched.rmse(
            bundle["y"],
            bundle["zm"]["prediction_sample"] + zm_result["best"]["values"],
        )

    payload = {
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "corpus": {
            "token_count": bundle["corpus"]["token_count"],
            "unique_words": bundle["corpus"]["unique_words"],
            "entropy": bundle["corpus"]["entropy"],
            "top_cities": bundle["corpus"]["top_cities"],
        },
        "linear_baseline": {
            "intercept": bundle["linear"]["intercept"],
            "slope": bundle["linear"]["slope"],
            "rmse_full": bundle["linear"]["rmse_full"],
            "rmse_sample": bundle["linear"]["rmse_sample"],
        },
        "zm_baseline": {
            "a": bundle["zm"]["a"],
            "b": bundle["zm"]["b"],
            "c": bundle["zm"]["c"],
            "rmse_full": bundle["zm"]["rmse_full"],
            "rmse_sample": bundle["zm"]["rmse_sample"],
        },
        "linear_search": {
            "best": enriched.sanitize_candidate(linear_result["best"]),
            "step_summary": linear_result["steps"],
        },
        "zm_search": {
            "best": enriched.sanitize_candidate(zm_result["best"]),
            "step_summary": zm_result["steps"],
        },
    }
    (args.outdir / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    plot_city_distribution(args.outdir, bundle)
    write_report(args.outdir, bundle, linear_result, zm_result)
    print(f"Saved {args.outdir / 'summary.json'}")
    print(f"Saved {args.outdir / 'report.md'}")


if __name__ == "__main__":
    main()
