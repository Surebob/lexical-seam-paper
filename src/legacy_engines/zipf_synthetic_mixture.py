import importlib.util
import json
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
ZIPF_PATH = ROOT / "eml_zipf_experiment.py"
SEARCH_PATH = ROOT / "eml_zipf_enriched_search.py"
OUTDIR = ROOT / "results" / "zipf_synthetic_mixture"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


zipf = load_module(ZIPF_PATH, "zipf_synthetic_zipf")
search = load_module(SEARCH_PATH, "zipf_synthetic_search")


def power_law(size: int, alpha: float, scale: float):
    ranks = np.arange(1, size + 1, dtype=np.float64)
    return scale / np.power(ranks, alpha)


def build_series(freqs: np.ndarray, x_low: float = 0.05, x_high: float = 1.0):
    freqs = np.asarray(freqs, dtype=np.float64)
    freqs = freqs[np.isfinite(freqs) & (freqs > 0.0)]
    freqs = np.sort(freqs)[::-1]
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)

    linear = zipf.fit_linear_zipf(log_rank, log_freq)
    zm = zipf.fit_zipf_mandelbrot(ranks, log_freq)
    x_full = zipf.normalize_x(log_rank, x_low, x_high)

    return {
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": log_rank,
        "log_freq": log_freq,
        "x": x_full,
        "linear": linear,
        "zm": zm,
        "token_equivalent": float(np.sum(freqs)),
        "vocab_equivalent": int(len(freqs)),
    }


def run_case(name: str, freqs: np.ndarray):
    series = build_series(freqs)

    linear_target = series["log_freq"] - series["linear"]["prediction"]
    linear_result = search.run_search(
        series["x"],
        linear_target,
        beam_width=50,
        max_steps=10,
        keep_all_until_step=2,
        diversity_weight=0.35,
        constant_variance_threshold=1e-10,
    )
    if linear_result["best"] is not None:
        linear_result["best"]["composite_rmse"] = search.rmse(
            series["log_freq"],
            series["linear"]["prediction"] + linear_result["best"]["values"],
        )

    zm_target = series["log_freq"] - series["zm"]["prediction"]
    zm_result = search.run_search(
        series["x"],
        zm_target,
        beam_width=50,
        max_steps=10,
        keep_all_until_step=2,
        diversity_weight=0.35,
        constant_variance_threshold=1e-10,
    )
    if zm_result["best"] is not None:
        zm_result["best"]["composite_rmse"] = search.rmse(
            series["log_freq"],
            series["zm"]["prediction"] + zm_result["best"]["values"],
        )

    step2 = next(item for item in zm_result["steps"] if item["step"] == 2)
    step2_top = step2["top_candidates"][0]

    payload = {
        "name": name,
        "token_equivalent": series["token_equivalent"],
        "vocab_equivalent": series["vocab_equivalent"],
        "zm_baseline": {
            "a": float(series["zm"]["a"]),
            "b": float(series["zm"]["b"]),
            "c": float(series["zm"]["c"]),
            "rmse": float(series["zm"]["rmse"]),
        },
        "step2": {
            "winner": step2_top["expr"],
            "math": step2_top["math"],
            "rmse": float(step2_top["rmse"]),
            "helps": bool(step2_top["rmse"] < series["zm"]["rmse"]),
            "top5": step2["top_candidates"][:5],
        },
        "best_overall": search.sanitize_candidate(zm_result["best"]),
    }
    return payload


def write_report(results):
    lines = [
        "# Synthetic Mixture Test",
        "",
        "All runs use the unchanged enriched search core with beam=50, max_steps=10, keep_all_until_step=2, sample_points=0 equivalent on full synthetic arrays.",
        "",
    ]
    for item in results:
        lines.extend(
            [
                f"## {item['name']}",
                "",
                f"- token equivalent: `{item['token_equivalent']:.6f}`",
                f"- vocab equivalent: `{item['vocab_equivalent']}`",
                f"- ZM params: `a={item['zm_baseline']['a']:.12f}`, `b={item['zm_baseline']['b']:.12f}`, `c={item['zm_baseline']['c']:.12f}`",
                f"- ZM RMSE: `{item['zm_baseline']['rmse']:.12f}`",
                f"- step-2 winner: `{item['step2']['winner']}`",
                f"- step-2 math: `{item['step2']['math']}`",
                f"- step-2 RMSE: `{item['step2']['rmse']:.12f}`",
                f"- step-2 helps ZM: `{item['step2']['helps']}`",
                f"- best overall composite RMSE: `{item['best_overall']['composite_rmse']:.12f}`",
                "",
                "Top 5 step-2 candidates:",
                "",
            ]
        )
        for rank, cand in enumerate(item["step2"]["top5"], start=1):
            lines.append(f"{rank}. `{cand['expr']}` — `{cand['rmse']:.12f}`")
        lines.append("")
    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    a_scale = 10.0
    b_scale = 1.0
    a_size = 50
    b_size = 10000

    alpha_a = 1.5
    alpha_b = 0.8

    component_a = power_law(a_size, alpha_a, a_scale)
    component_b = power_law(b_size, alpha_b, b_scale)
    mix_medium = np.concatenate([component_a, component_b])
    control_bb = np.concatenate([component_b, component_b])

    component_b_small_gap = power_law(b_size, 1.3, b_scale)
    mix_small_gap = np.concatenate([component_a, component_b_small_gap])

    component_b_large_gap = power_law(b_size, 0.5, b_scale)
    mix_large_gap = np.concatenate([component_a, component_b_large_gap])

    cases = [
        ("component_a_only_alpha1.5_n50", component_a),
        ("component_b_only_alpha0.8_n10000", component_b),
        ("mixture_alpha1.5_plus_0.8", mix_medium),
        ("control_two_copies_alpha0.8", control_bb),
        ("mixture_alpha1.5_plus_1.3", mix_small_gap),
        ("mixture_alpha1.5_plus_0.5", mix_large_gap),
    ]

    results = [run_case(name, freqs) for name, freqs in cases]

    payload = {
        "setup": {
            "component_a": {"size": a_size, "alpha": alpha_a, "scale": a_scale},
            "component_b": {"size": b_size, "alpha": alpha_b, "scale": b_scale},
            "control": "two copies of component B merged and reranked",
            "extra_gap_cases": [
                {"alpha_a": 1.5, "alpha_b": 1.3},
                {"alpha_a": 1.5, "alpha_b": 0.5},
            ],
        },
        "results": results,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(results)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
