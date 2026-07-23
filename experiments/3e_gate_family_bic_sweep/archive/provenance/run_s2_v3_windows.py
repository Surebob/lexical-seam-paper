from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import multiprocessing as mp
import os
import sys
import time
import traceback
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
SHARED_DIR = PACKAGE_ROOT / "shared"
REFERENCE_DIR = PACKAGE_ROOT / "reference"
OUTPUTS_DIR = PACKAGE_ROOT / "outputs"

CFG = None
LOADER = None
MODEL = None


def worker_init(shared_dir: str, blas_threads: int = 2) -> None:
    thread_value = str(blas_threads)
    for key in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
    ]:
        os.environ[key] = thread_value
    if shared_dir not in sys.path:
        sys.path.insert(0, shared_dir)

    global CFG, LOADER, MODEL
    import fit_config as cfg
    import corpus_loader as loader
    import decoupled_smooth_model as model

    CFG = cfg
    LOADER = loader
    MODEL = model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Windows runner for the S2 v3 decoupled five-gate sweep.")
    parser.add_argument("--base-seed", type=int, default=20260415)
    parser.add_argument("--n-starts", type=int, default=100)
    parser.add_argument("--max-nfev", type=int, default=12000)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--blas-threads", type=int, default=2)
    parser.add_argument("--full-sweep-confirmed", action="store_true")
    parser.add_argument("--skip-smoke", action="store_true")
    parser.add_argument("--output-root", type=Path, default=OUTPUTS_DIR)
    return parser.parse_args()


def ensure_shared_imports() -> None:
    if str(SHARED_DIR) not in sys.path:
        sys.path.insert(0, str(SHARED_DIR))
    global CFG
    if CFG is None:
        import fit_config as cfg

        CFG = cfg


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_csv(path: Path, key_field: str | None = None):
    if not path.exists():
        return {} if key_field else []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if key_field is None:
        return rows
    return {row[key_field]: row for row in rows}


def log_message(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")


def preflight_check() -> dict:
    ensure_shared_imports()
    messages = []
    ok = True
    versions = {}
    for package, expected in CFG.EXPECTED_DEPENDENCIES.items():
        try:
            actual = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            actual = None
        versions[package] = actual
        if actual != expected:
            ok = False
            messages.append(f"{package}: expected {expected}, found {actual}")
    missing = [spec["filename"] for spec in CFG.SEARCHED_CORPORA if not (CFG.CORPORA_DIR / spec["filename"]).exists()]
    if missing:
        ok = False
        messages.append(f"missing corpora: {', '.join(missing)}")
    return {"ok": ok, "versions": versions, "messages": messages}


def gate_fit_fieldnames() -> list[str]:
    return [
        "slug",
        "corpus",
        "gate",
        "seed",
        "bic",
        "rmse",
        "k",
        "w_gate",
        "w_tail",
        "best_start_index",
        "best_nfev",
        "runtime_sec",
        "k_hit_lower_bound",
        "k_hit_upper_bound",
        "w_gate_hit_lower_bound",
        "w_gate_hit_upper_bound",
        "w_tail_hit_lower_bound",
        "w_tail_hit_upper_bound",
        "wall_clock_seconds",
        "peak_memory_mb",
        "starts_within_1_bic",
    ]


def dispersion_fieldnames() -> list[str]:
    return [
        "slug",
        "corpus",
        "gate",
        "seed",
        "min_bic",
        "median_bic",
        "max_bic",
        "starts_within_1_bic",
        "notes",
    ]


def per_corpus_fieldnames() -> list[str]:
    ensure_shared_imports()
    fields = ["slug", "corpus", "seed", "vocab_size"]
    for gate in CFG.ORDERED_GATES:
        fields.extend(
            [
                f"{gate}_bic",
                f"{gate}_rmse",
                f"{gate}_k",
                f"{gate}_w_gate",
                f"{gate}_w_tail",
                f"{gate}_best_start_index",
                f"{gate}_best_nfev",
                f"{gate}_runtime_sec",
                f"{gate}_peak_memory_mb",
                f"{gate}_starts_within_1_bic",
                f"{gate}_k_hit_lower_bound",
                f"{gate}_k_hit_upper_bound",
                f"{gate}_w_gate_hit_lower_bound",
                f"{gate}_w_gate_hit_upper_bound",
                f"{gate}_w_tail_hit_lower_bound",
                f"{gate}_w_tail_hit_upper_bound",
            ]
        )
    fields.extend(
        [
            "bic_tanh_minus_logistic",
            "w_gate_tanh_over_2w_gate_logistic_ratio",
            "w_tail_tanh_over_w_tail_logistic_ratio",
            "tanh_calibration_pass",
            "winner_gate",
            "worst_gate",
            "bic_spread",
        ]
    )
    return fields


def tanh_diag_fieldnames() -> list[str]:
    return [
        "slug",
        "corpus",
        "bic_tanh_minus_logistic",
        "w_gate_tanh_over_2w_gate_logistic_ratio",
        "w_tail_tanh_over_w_tail_logistic_ratio",
        "tanh_calibration_pass",
    ]


def corpus_task(spec: dict, base_seed: int, n_starts: int, max_nfev: int) -> dict:
    dataset = LOADER.build_zipf_dataset(CFG.CORPORA_DIR / spec["filename"])
    per_fit_rows = []
    dispersion_rows = []
    fits = {}
    for gate in CFG.ORDERED_GATES:
        start_time = time.time()
        fit = MODEL.fit_with_gate(dataset, gate, base_seed=base_seed, n_starts=n_starts, max_nfev=max_nfev)
        wall_clock = time.time() - start_time
        fits[gate] = fit
        per_fit_rows.append(
            {
                "slug": spec["slug"],
                "corpus": spec["name"],
                "gate": gate,
                "seed": str(base_seed),
                "bic": repr(float(fit["bic"])),
                "rmse": repr(float(fit["rmse"])),
                "k": repr(float(fit["k"])),
                "w_gate": repr(float(fit["w_gate"])),
                "w_tail": repr(float(fit["w_tail"])),
                "best_start_index": str(int(fit["best_start_index"])),
                "best_nfev": str(int(fit["best_nfev"])),
                "runtime_sec": repr(float(wall_clock)),
                "k_hit_lower_bound": str(bool(fit["k_hit_lower_bound"])),
                "k_hit_upper_bound": str(bool(fit["k_hit_upper_bound"])),
                "w_gate_hit_lower_bound": str(bool(fit["w_gate_hit_lower_bound"])),
                "w_gate_hit_upper_bound": str(bool(fit["w_gate_hit_upper_bound"])),
                "w_tail_hit_lower_bound": str(bool(fit["w_tail_hit_lower_bound"])),
                "w_tail_hit_upper_bound": str(bool(fit["w_tail_hit_upper_bound"])),
                "wall_clock_seconds": repr(float(wall_clock)),
                "peak_memory_mb": repr(float(fit["peak_memory_mb"])),
                "starts_within_1_bic": str(int(fit["starts_within_1_bic"])),
            }
        )
        dispersion_rows.append(
            {
                "slug": spec["slug"],
                "corpus": spec["name"],
                "gate": gate,
                "seed": str(base_seed),
                "min_bic": repr(float(fit["min_bic_across_starts"])),
                "median_bic": repr(float(fit["median_bic_across_starts"])),
                "max_bic": repr(float(fit["max_bic_across_starts"])),
                "starts_within_1_bic": str(int(fit["starts_within_1_bic"])),
                "notes": "Across starts within a single seed run.",
            }
        )

    tanh_bic_diff = float(fits["tanh"]["bic"]) - float(fits["logistic"]["bic"])
    tanh_gate_ratio = float(fits["tanh"]["w_gate"]) / (2.0 * float(fits["logistic"]["w_gate"]))
    tanh_tail_ratio = float(fits["tanh"]["w_tail"]) / float(fits["logistic"]["w_tail"])
    tanh_pass = abs(tanh_bic_diff) < CFG.TANH_BIC_THRESHOLD and abs(tanh_gate_ratio - 1.0) < CFG.TANH_GATE_RATIO_THRESHOLD

    independent_bics = {gate: float(fits[gate]["bic"]) for gate in CFG.INDEPENDENT_GATES}
    winner_gate = min(independent_bics, key=independent_bics.get)
    worst_gate = max(independent_bics, key=independent_bics.get)
    bic_spread = max(independent_bics.values()) - min(independent_bics.values())

    per_corpus_row = {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "seed": str(base_seed),
        "vocab_size": str(int(dataset["unique_words"])),
    }
    tanh_diag_row = {
        "slug": spec["slug"],
        "corpus": spec["name"],
        "bic_tanh_minus_logistic": repr(float(tanh_bic_diff)),
        "w_gate_tanh_over_2w_gate_logistic_ratio": repr(float(tanh_gate_ratio)),
        "w_tail_tanh_over_w_tail_logistic_ratio": repr(float(tanh_tail_ratio)),
        "tanh_calibration_pass": str(bool(tanh_pass)),
    }
    for gate in CFG.ORDERED_GATES:
        fit = fits[gate]
        for key in [
            "bic",
            "rmse",
            "k",
            "w_gate",
            "w_tail",
            "best_start_index",
            "best_nfev",
            "peak_memory_mb",
            "starts_within_1_bic",
            "k_hit_lower_bound",
            "k_hit_upper_bound",
            "w_gate_hit_lower_bound",
            "w_gate_hit_upper_bound",
            "w_tail_hit_lower_bound",
            "w_tail_hit_upper_bound",
        ]:
            value = fit[key]
            per_corpus_row[f"{gate}_{key}"] = repr(float(value)) if isinstance(value, float) else str(value)
        runtime = next(row for row in per_fit_rows if row["gate"] == gate)["runtime_sec"]
        per_corpus_row[f"{gate}_runtime_sec"] = runtime
    per_corpus_row.update(tanh_diag_row)
    per_corpus_row["winner_gate"] = winner_gate
    per_corpus_row["worst_gate"] = worst_gate
    per_corpus_row["bic_spread"] = repr(float(bic_spread))
    return {
        "slug": spec["slug"],
        "per_fit_rows": per_fit_rows,
        "dispersion_rows": dispersion_rows,
        "per_corpus_row": per_corpus_row,
        "tanh_diag_row": tanh_diag_row,
    }


def build_aggregate_rows(per_corpus_rows: list[dict]) -> list[dict]:
    ensure_shared_imports()
    if not per_corpus_rows:
        return []

    rows = []
    winner_counts = {gate: 0 for gate in CFG.INDEPENDENT_GATES}
    bic_spreads = []
    tanh_pass_count = 0

    for row in per_corpus_rows:
        if row["tanh_calibration_pass"] == "True":
            tanh_pass_count += 1
        independent_bics = {gate: float(row[f"{gate}_bic"]) for gate in CFG.INDEPENDENT_GATES}
        winner_counts[min(independent_bics, key=independent_bics.get)] += 1
        bic_spreads.append(float(row["bic_spread"]))

    rows.append(
        {
            "metric_name": "tanh_calibration_pass_count",
            "value": str(tanh_pass_count),
            "notes": "Count of corpora satisfying the tanh calibration checks.",
        }
    )
    for gate in CFG.INDEPENDENT_GATES:
        rows.append(
            {
                "metric_name": f"{gate}_bic_wins",
                "value": str(winner_counts[gate]),
                "notes": f"Count of corpora where {gate} has the lowest BIC among the independent gates.",
            }
        )
    rows.extend(
        [
            {
                "metric_name": "gates_indistinguishable_strict_count",
                "value": str(sum(spread < 2.0 for spread in bic_spreads)),
                "notes": "Count of corpora with independent-gate BIC spread < 2.",
            },
            {
                "metric_name": "gates_indistinguishable_positive_count",
                "value": str(sum(spread < 6.0 for spread in bic_spreads)),
                "notes": "Count of corpora with independent-gate BIC spread < 6.",
            },
            {
                "metric_name": "gates_indistinguishable_strong_count",
                "value": str(sum(spread < 10.0 for spread in bic_spreads)),
                "notes": "Count of corpora with independent-gate BIC spread < 10.",
            },
            {
                "metric_name": "median_bic_spread",
                "value": repr(float(sorted(bic_spreads)[len(bic_spreads) // 2] if len(bic_spreads) % 2 == 1 else (sorted(bic_spreads)[len(bic_spreads)//2 - 1] + sorted(bic_spreads)[len(bic_spreads)//2]) / 2.0)),
                "notes": "Median BIC spread across the independent gates.",
            },
        ]
    )
    rows.append(
        {
            "metric_name": "mean_bic_spread",
            "value": repr(float(sum(bic_spreads) / len(bic_spreads))),
            "notes": "Mean BIC spread across the independent gates.",
        }
    )
    for gate in CFG.ORDERED_GATES:
        for param in ["w_gate", "w_tail"]:
            vals = [float(row[f"{gate}_{param}"]) for row in per_corpus_rows]
            vals_sorted = sorted(vals)
            median = vals_sorted[len(vals_sorted) // 2] if len(vals_sorted) % 2 == 1 else (vals_sorted[len(vals_sorted)//2 - 1] + vals_sorted[len(vals_sorted)//2]) / 2.0
            rows.extend(
                [
                    {"metric_name": f"{gate}_{param}_min", "value": repr(float(min(vals))), "notes": f"Minimum fitted {param} for {gate}."},
                    {"metric_name": f"{gate}_{param}_max", "value": repr(float(max(vals))), "notes": f"Maximum fitted {param} for {gate}."},
                    {"metric_name": f"{gate}_{param}_median", "value": repr(float(median)), "notes": f"Median fitted {param} for {gate}."},
                ]
            )
        for bound_name in [
            "k_hit_lower_bound",
            "k_hit_upper_bound",
            "w_gate_hit_lower_bound",
            "w_gate_hit_upper_bound",
            "w_tail_hit_lower_bound",
            "w_tail_hit_upper_bound",
        ]:
            rows.append(
                {
                    "metric_name": f"{gate}_{bound_name}_count",
                    "value": str(sum(row[f"{gate}_{bound_name}"] == "True" for row in per_corpus_rows)),
                    "notes": f"Count of corpora where {gate} hit {bound_name}.",
                }
            )
    return rows


def write_outputs(full_dir: Path, per_fit_rows: list[dict], dispersion_rows: list[dict], per_corpus_rows: list[dict], tanh_rows: list[dict]) -> None:
    write_csv(full_dir / "s2_v3_per_fit_results.csv", per_fit_rows, gate_fit_fieldnames())
    write_csv(full_dir / "s2_v3_per_start_dispersion.csv", dispersion_rows, dispersion_fieldnames())
    write_csv(full_dir / "s2_v3_per_corpus_results.csv", per_corpus_rows, per_corpus_fieldnames())
    write_csv(full_dir / "s2_v3_tanh_calibration.csv", tanh_rows, tanh_diag_fieldnames())
    write_csv(full_dir / "s2_v3_aggregate_statistics.csv", build_aggregate_rows(per_corpus_rows), ["metric_name", "value", "notes"])


def run_smoke_test(base_seed: int, n_starts: int, max_nfev: int, output_root: Path) -> dict:
    ensure_shared_imports()
    smoke_dir = output_root / "smoke"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    worker_init(str(SHARED_DIR), CFG.DEFAULT_BLAS_THREADS)
    shakespeare = next(spec for spec in CFG.SEARCHED_CORPORA if spec["slug"] == "shakespeare")
    result = corpus_task(shakespeare, base_seed=base_seed, n_starts=n_starts, max_nfev=max_nfev)

    expected = json.loads((REFERENCE_DIR / "expected_shakespeare_bic.json").read_text(encoding="utf-8"))
    comparisons = []
    smoke_pass = True
    for row in result["per_fit_rows"]:
        gate = row["gate"]
        actual_bic = float(row["bic"])
        expected_bic = float(expected[gate]["bic"])
        abs_diff = abs(actual_bic - expected_bic)
        if abs_diff > CFG.SMOKE_BIC_TOLERANCE:
            smoke_pass = False
        comparisons.append(
            {
                "gate": gate,
                "expected_bic": repr(float(expected_bic)),
                "actual_bic": repr(float(actual_bic)),
                "abs_diff": repr(float(abs_diff)),
                "passes_1e_6": str(abs_diff <= CFG.SMOKE_BIC_TOLERANCE),
            }
        )

    write_csv(smoke_dir / "s2_v3_smoke_per_fit_results.csv", result["per_fit_rows"], gate_fit_fieldnames())
    write_csv(smoke_dir / "s2_v3_smoke_per_start_dispersion.csv", result["dispersion_rows"], dispersion_fieldnames())
    write_csv(smoke_dir / "s2_v3_smoke_comparison.csv", comparisons, ["gate", "expected_bic", "actual_bic", "abs_diff", "passes_1e_6"])
    (smoke_dir / "s2_v3_smoke_summary.json").write_text(
        json.dumps(
            {
                "passes_smoke": smoke_pass,
                "seed": base_seed,
                "n_starts": n_starts,
                "max_nfev": max_nfev,
                "tolerance": CFG.SMOKE_BIC_TOLERANCE,
                "max_abs_diff": max(float(row["abs_diff"]) for row in comparisons),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return {"passes_smoke": smoke_pass, "comparisons": comparisons, "result": result}


def load_existing_full_rows(full_dir: Path) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    return (
        load_csv(full_dir / "s2_v3_per_fit_results.csv"),
        load_csv(full_dir / "s2_v3_per_start_dispersion.csv"),
        load_csv(full_dir / "s2_v3_per_corpus_results.csv"),
        load_csv(full_dir / "s2_v3_tanh_calibration.csv"),
    )


def run_full_sweep(args: argparse.Namespace) -> None:
    ensure_shared_imports()
    full_dir = args.output_root / "full"
    full_dir.mkdir(parents=True, exist_ok=True)
    runtime_log = full_dir / "s2_v3_runtime_log.txt"
    progress_log = full_dir / "progress.log"

    per_fit_rows, dispersion_rows, per_corpus_rows, tanh_rows = load_existing_full_rows(full_dir)
    completed = {row["slug"] for row in per_corpus_rows}
    pending = [spec for spec in CFG.SEARCHED_CORPORA if spec["slug"] not in completed]

    log_message(runtime_log, f"starting full sweep with {len(pending)} pending corpora")
    tasks = [(spec, args.base_seed, args.n_starts, args.max_nfev) for spec in pending]
    if not tasks:
        write_outputs(full_dir, per_fit_rows, dispersion_rows, per_corpus_rows, tanh_rows)
        log_message(runtime_log, "nothing to do; all corpora already complete")
        return

    with mp.Pool(processes=args.workers, initializer=worker_init, initargs=(str(SHARED_DIR), args.blas_threads)) as pool:
        for payload in pool.starmap(corpus_task, tasks):
            slug = payload["slug"]
            per_fit_rows = [row for row in per_fit_rows if row["slug"] != slug] + payload["per_fit_rows"]
            dispersion_rows = [row for row in dispersion_rows if row["slug"] != slug] + payload["dispersion_rows"]
            per_corpus_rows = [row for row in per_corpus_rows if row["slug"] != slug] + [payload["per_corpus_row"]]
            tanh_rows = [row for row in tanh_rows if row["slug"] != slug] + [payload["tanh_diag_row"]]
            per_fit_rows.sort(key=lambda row: (CFG.SEARCHED_SLUGS.index(row["slug"]), CFG.ORDERED_GATES.index(row["gate"])))
            dispersion_rows.sort(key=lambda row: (CFG.SEARCHED_SLUGS.index(row["slug"]), CFG.ORDERED_GATES.index(row["gate"])))
            per_corpus_rows.sort(key=lambda row: CFG.SEARCHED_SLUGS.index(row["slug"]))
            tanh_rows.sort(key=lambda row: CFG.SEARCHED_SLUGS.index(row["slug"]))
            write_outputs(full_dir, per_fit_rows, dispersion_rows, per_corpus_rows, tanh_rows)
            log_message(progress_log, f"completed {slug}")
            log_message(runtime_log, f"completed {slug}")


def main() -> int:
    args = parse_args()
    preflight = preflight_check()
    if not preflight["ok"]:
        print("Pre-flight check failed:")
        for message in preflight["messages"]:
            print(f"- {message}")
        return 2

    smoke = None
    if not args.skip_smoke:
        smoke = run_smoke_test(args.base_seed, args.n_starts, args.max_nfev, args.output_root)
        max_diff = max(float(row["abs_diff"]) for row in smoke["comparisons"])
        print(f"Smoke test pass: {smoke['passes_smoke']} (max abs diff {max_diff:.12g})")
        if not smoke["passes_smoke"]:
            print("Smoke test failed; not launching full sweep.")
            return 3
    if not args.full_sweep_confirmed:
        print("Smoke test completed. Re-run with --full-sweep-confirmed to launch the full parallel sweep.")
        return 0

    run_full_sweep(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

