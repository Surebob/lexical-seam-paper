from __future__ import annotations

import argparse
import csv
import math
import multiprocessing as mp
import os
import statistics
import sys
import time
from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
SHARED_DIR = ROOT / "s2_v3_windows_port" / "shared"
PARAMS_PATH = ROOT / "results" / "windows_s2_decoupled_v3_results_2026-04-18" / "outputs" / "s2_v3_logistic_local_params.csv"
OUTDIR = ROOT / "phase2_addon" / "s2_logistic_synthetic_gate_recovery"

BASE_SEED = 20260415
N_STARTS = 100
MAX_NFEV = 12000
DEFAULT_WORKERS = 6

CFG = None
LOADER = None
MODEL = None


def worker_init() -> None:
    for key in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
    ]:
        os.environ[key] = "1"
    if str(SHARED_DIR) not in sys.path:
        sys.path.insert(0, str(SHARED_DIR))
    global CFG, LOADER, MODEL
    import fit_config as cfg
    import corpus_loader as loader
    import decoupled_smooth_model as model

    CFG = cfg
    LOADER = loader
    MODEL = model


def ensure_imports() -> None:
    if MODEL is None:
        worker_init()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Synthetic logistic-gate recovery check for S2 v3.")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--n-starts", type=int, default=N_STARTS)
    parser.add_argument("--max-nfev", type=int, default=MAX_NFEV)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_generator_params() -> dict[str, dict]:
    return {row["slug"]: row for row in read_csv(PARAMS_PATH)}


def params_array(row: dict):
    import numpy as np

    return np.array(
        [
            float(row["logistic_a1"]),
            float(row["logistic_b1"]),
            float(row["logistic_c1"]),
            float(row["logistic_a2"]),
            float(row["logistic_b2"]),
            float(row["logistic_c2"]),
            float(row["logistic_k"]),
            float(row["logistic_w_gate"]),
            float(row["logistic_w_tail"]),
        ],
        dtype=np.float64,
    )


def fit_synthetic_corpus(spec: dict, generator_row: dict, n_starts: int, max_nfev: int) -> dict:
    ensure_imports()
    import numpy as np

    dataset = LOADER.build_zipf_dataset(CFG.CORPORA_DIR / spec["filename"])
    generating_params = params_array(generator_row)
    synthetic_log_freq = MODEL.decoupled_prediction(
        dataset["ranks"],
        dataset["log_rank"],
        generating_params,
        "logistic",
    )
    synthetic_dataset = dict(dataset)
    synthetic_dataset["log_freq"] = synthetic_log_freq
    synthetic_dataset["freqs"] = np.exp(synthetic_log_freq)

    fit_rows = []
    dispersion_rows = []
    fits = {}
    start_time = time.time()
    for gate in CFG.ORDERED_GATES:
        gate_start = time.time()
        fit = MODEL.fit_with_gate(
            synthetic_dataset,
            gate,
            base_seed=BASE_SEED,
            n_starts=n_starts,
            max_nfev=max_nfev,
        )
        wall_clock = time.time() - gate_start
        fits[gate] = fit
        fit_rows.append(
            {
                "slug": spec["slug"],
                "corpus": spec["name"],
                "gate": gate,
                "bic": repr(float(fit["bic"])),
                "rmse": repr(float(fit["rmse"])),
                "k": repr(float(fit["k"])),
                "w_gate": repr(float(fit["w_gate"])),
                "w_tail": repr(float(fit["w_tail"])),
                "best_start_index": str(int(fit["best_start_index"])),
                "best_nfev": str(int(fit["best_nfev"])),
                "wall_clock_seconds": repr(float(wall_clock)),
                "starts_within_1_bic": str(int(fit["starts_within_1_bic"])),
                "k_hit_lower_bound": str(bool(fit["k_hit_lower_bound"])),
                "k_hit_upper_bound": str(bool(fit["k_hit_upper_bound"])),
                "w_gate_hit_lower_bound": str(bool(fit["w_gate_hit_lower_bound"])),
                "w_gate_hit_upper_bound": str(bool(fit["w_gate_hit_upper_bound"])),
                "w_tail_hit_lower_bound": str(bool(fit["w_tail_hit_lower_bound"])),
                "w_tail_hit_upper_bound": str(bool(fit["w_tail_hit_upper_bound"])),
            }
        )
        dispersion_rows.append(
            {
                "slug": spec["slug"],
                "corpus": spec["name"],
                "gate": gate,
                "min_bic": repr(float(fit["min_bic_across_starts"])),
                "median_bic": repr(float(fit["median_bic_across_starts"])),
                "max_bic": repr(float(fit["max_bic_across_starts"])),
                "starts_within_1_bic": str(int(fit["starts_within_1_bic"])),
            }
        )

    independent_bics = {gate: float(fits[gate]["bic"]) for gate in CFG.INDEPENDENT_GATES}
    all_bics = {gate: float(fits[gate]["bic"]) for gate in CFG.ORDERED_GATES}
    winner_gate = min(independent_bics, key=independent_bics.get)
    worst_gate = max(independent_bics, key=independent_bics.get)
    bic_spread = max(independent_bics.values()) - min(independent_bics.values())
    tanh_bic_diff = float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"])
    tanh_gate_ratio = float(fits["tanh"]["w_gate"]) / (2.0 * float(fits["logistic"]["w_gate"]))
    tanh_pass = abs(tanh_bic_diff) < 1.0 and abs(tanh_gate_ratio - 1.0) < 0.05
    logistic_independent_rank = 1 + sum(v < float(fits["logistic"]["bic"]) for gate, v in independent_bics.items() if gate != "logistic")

    per_corpus_row = {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "vocab_size": str(int(dataset["unique_words"])),
        "generator_gate": "logistic",
        "generator_k": generator_row["logistic_k"],
        "generator_w_gate": generator_row["logistic_w_gate"],
        "generator_w_tail": generator_row["logistic_w_tail"],
        "winner_gate": winner_gate,
        "worst_gate": worst_gate,
        "bic_spread": repr(float(bic_spread)),
        "logistic_independent_rank": str(int(logistic_independent_rank)),
        "bic_erf_minus_logistic": repr(float(fits["erf"]["bic"]) - float(fits["logistic"]["bic"])),
        "bic_tanh_minus_logistic": repr(float(tanh_bic_diff)),
        "w_gate_tanh_over_2w_gate_logistic_ratio": repr(float(tanh_gate_ratio)),
        "tanh_calibration_pass": str(bool(tanh_pass)),
        "total_wall_clock_seconds": repr(float(time.time() - start_time)),
    }
    for gate in CFG.ORDERED_GATES:
        for key in ["bic", "rmse", "k", "w_gate", "w_tail", "starts_within_1_bic"]:
            per_corpus_row[f"{gate}_{key}"] = repr(float(fits[gate][key])) if isinstance(fits[gate][key], float) else str(fits[gate][key])
    return {
        "slug": spec["slug"],
        "fit_rows": fit_rows,
        "dispersion_rows": dispersion_rows,
        "per_corpus_row": per_corpus_row,
    }


def fit_synthetic_corpus_task(task: tuple[dict, dict, int, int]) -> dict:
    spec, generator_row, n_starts, max_nfev = task
    return fit_synthetic_corpus(spec, generator_row, n_starts, max_nfev)


def write_all_outputs(results: list[dict]) -> None:
    ensure_imports()
    order = {slug: idx for idx, slug in enumerate(CFG.SEARCHED_SLUGS)}
    results = sorted(results, key=lambda item: order[item["slug"]])
    fit_rows = [row for result in results for row in result["fit_rows"]]
    dispersion_rows = [row for result in results for row in result["dispersion_rows"]]
    per_corpus_rows = [result["per_corpus_row"] for result in results]

    fit_rows.sort(key=lambda row: (order[row["slug"]], CFG.ORDERED_GATES.index(row["gate"])))
    dispersion_rows.sort(key=lambda row: (order[row["slug"]], CFG.ORDERED_GATES.index(row["gate"])))

    write_csv(
        OUTDIR / "synthetic_gate_recovery_per_fit.csv",
        fit_rows,
        [
            "slug",
            "corpus",
            "gate",
            "bic",
            "rmse",
            "k",
            "w_gate",
            "w_tail",
            "best_start_index",
            "best_nfev",
            "wall_clock_seconds",
            "starts_within_1_bic",
            "k_hit_lower_bound",
            "k_hit_upper_bound",
            "w_gate_hit_lower_bound",
            "w_gate_hit_upper_bound",
            "w_tail_hit_lower_bound",
            "w_tail_hit_upper_bound",
        ],
    )
    write_csv(
        OUTDIR / "synthetic_gate_recovery_per_start_dispersion.csv",
        dispersion_rows,
        ["slug", "corpus", "gate", "min_bic", "median_bic", "max_bic", "starts_within_1_bic"],
    )
    per_corpus_fields = [
        "slug",
        "corpus",
        "vocab_size",
        "generator_gate",
        "generator_k",
        "generator_w_gate",
        "generator_w_tail",
        "winner_gate",
        "worst_gate",
        "bic_spread",
        "logistic_independent_rank",
        "bic_erf_minus_logistic",
        "bic_tanh_minus_logistic",
        "w_gate_tanh_over_2w_gate_logistic_ratio",
        "tanh_calibration_pass",
        "total_wall_clock_seconds",
    ]
    for gate in CFG.ORDERED_GATES:
        per_corpus_fields.extend([f"{gate}_{key}" for key in ["bic", "rmse", "k", "w_gate", "w_tail", "starts_within_1_bic"]])
    write_csv(OUTDIR / "synthetic_gate_recovery_per_corpus.csv", per_corpus_rows, per_corpus_fields)

    aggregate_rows = build_aggregate_rows(per_corpus_rows, fit_rows)
    write_csv(OUTDIR / "synthetic_gate_recovery_aggregate_statistics.csv", aggregate_rows, ["metric_name", "value", "notes"])
    (OUTDIR / "synthetic_gate_recovery_summary_report.md").write_text(build_report(per_corpus_rows, aggregate_rows), encoding="utf-8")


def build_aggregate_rows(per_corpus_rows: list[dict], fit_rows: list[dict]) -> list[dict]:
    ensure_imports()
    rows = []
    winner_counts = {gate: 0 for gate in CFG.INDEPENDENT_GATES}
    for row in per_corpus_rows:
        winner_counts[row["winner_gate"]] += 1
    bic_spreads = [float(row["bic_spread"]) for row in per_corpus_rows]
    erf_minus_logistic = [float(row["bic_erf_minus_logistic"]) for row in per_corpus_rows]
    logistic_rank = [int(row["logistic_independent_rank"]) for row in per_corpus_rows]

    def add(name: str, value, notes: str) -> None:
        rows.append({"metric_name": name, "value": repr(float(value)) if isinstance(value, float) else str(value), "notes": notes})

    add("corpus_count", len(per_corpus_rows), "Number of synthetic corpora generated from known logistic decoupled model.")
    for gate in CFG.INDEPENDENT_GATES:
        add(f"{gate}_bic_wins", winner_counts[gate], f"Count of synthetic corpora where {gate} has lowest BIC among independent gates.")
    add("tanh_calibration_pass_count", sum(row["tanh_calibration_pass"] == "True" for row in per_corpus_rows), "Count passing tanh/logistic calibration.")
    add("logistic_recovered_count", winner_counts["logistic"], "Synthetic generator is logistic; this should be high if gate family is recovered.")
    add("erf_beats_logistic_count", sum(v < 0.0 for v in erf_minus_logistic), "Count where erf BIC is lower than logistic on logistic-generated data.")
    add("median_bic_spread", statistics.median(bic_spreads), "Median independent-gate BIC spread.")
    add("mean_bic_spread", statistics.fmean(bic_spreads), "Mean independent-gate BIC spread.")
    add("median_bic_erf_minus_logistic", statistics.median(erf_minus_logistic), "Negative means erf beats logistic.")
    add("mean_bic_erf_minus_logistic", statistics.fmean(erf_minus_logistic), "Negative means erf beats logistic.")
    add("max_logistic_independent_rank", max(logistic_rank), "Worst rank of logistic among independent gates.")
    for gate in CFG.ORDERED_GATES:
        starts = [int(row["starts_within_1_bic"]) for row in fit_rows if row["gate"] == gate]
        add(f"{gate}_starts_within_1_bic_min", min(starts), f"Minimum starts within 1 BIC for {gate}.")
        add(f"{gate}_starts_within_1_bic_median", statistics.median(starts), f"Median starts within 1 BIC for {gate}.")
    return rows


def build_report(per_corpus_rows: list[dict], aggregate_rows: list[dict]) -> str:
    agg = {row["metric_name"]: row["value"] for row in aggregate_rows}
    lines = [
        "# Synthetic Logistic-Gate Recovery Check",
        "",
        "- synthetic generator: exact decoupled logistic model using local S2 v3 logistic parameters",
        f"- corpora: `{agg.get('corpus_count', 'n/a')}`",
        f"- logistic wins: `{agg.get('logistic_bic_wins', 'n/a')}`",
        f"- erf wins: `{agg.get('erf_bic_wins', 'n/a')}`",
        f"- algebraic wins: `{agg.get('algebraic_bic_wins', 'n/a')}`",
        f"- arctan wins: `{agg.get('arctan_bic_wins', 'n/a')}`",
        f"- tanh calibration pass count: `{agg.get('tanh_calibration_pass_count', 'n/a')}`",
        f"- erf beats logistic count: `{agg.get('erf_beats_logistic_count', 'n/a')}`",
        f"- median BIC(erf - logistic): `{agg.get('median_bic_erf_minus_logistic', 'n/a')}`",
        "",
    ]
    if agg.get("erf_beats_logistic_count") == "0":
        lines.append("Interpretation: erf does not win on logistic-generated synthetic data, so the empirical erf preference is not explained by a generic fitter artifact in this recovery check.")
    else:
        lines.append("Interpretation: erf wins on some logistic-generated synthetic data, so the empirical erf preference may include fitting or parameterization artifact and needs further investigation.")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    progress_path = OUTDIR / "progress.log"
    progress_path.write_text("", encoding="utf-8")

    worker_init()
    params_by_slug = load_generator_params()
    specs = [spec for spec in CFG.SEARCHED_CORPORA if spec["slug"] in params_by_slug]
    tasks = [(spec, params_by_slug[spec["slug"]], args.n_starts, args.max_nfev) for spec in specs]

    results = []
    with mp.Pool(processes=args.workers, initializer=worker_init) as pool:
        for idx, result in enumerate(pool.imap_unordered(fit_synthetic_corpus_task, tasks), start=1):
            results.append(result)
            write_all_outputs(results)
            with progress_path.open("a", encoding="utf-8") as handle:
                handle.write(f"completed {idx}/{len(tasks)}: {result['slug']}\n")


if __name__ == "__main__":
    main()
