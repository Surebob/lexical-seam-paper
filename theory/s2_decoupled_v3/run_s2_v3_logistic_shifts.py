from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
FULL_SCRIPT = ROOT / "phase2_addon" / "s2_decoupled_v3" / "run_s2_v3_full.py"
OUTDIR = ROOT / "results" / "windows_s2_decoupled_v3_results_2026-04-18" / "outputs"
LOGISTIC_PARAMS_PATH = OUTDIR / "s2_v3_logistic_local_params.csv"
SHIFT_PATH = OUTDIR / "s2_v3_parameter_shifts.csv"
PROGRESS_PATH = OUTDIR / "s2_v3_parameter_shifts_progress.log"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


full = load_module(FULL_SCRIPT, "s2_v3_full_module")


def load_existing_rows() -> dict[str, dict]:
    if not LOGISTIC_PARAMS_PATH.exists():
        return {}
    with LOGISTIC_PARAMS_PATH.open(newline="", encoding="utf-8") as handle:
        return {row["slug"]: row for row in csv.DictReader(handle)}


def write_logistic_rows(rows_by_slug: dict[str, dict]) -> None:
    ordered = [rows_by_slug[slug] for slug in full.SEARCHED_SLUGS if slug in rows_by_slug]
    fieldnames = [
        "slug",
        "corpus",
        "vocab_size",
        "logistic_bic",
        "logistic_rmse",
        "logistic_a1",
        "logistic_b1",
        "logistic_c1",
        "logistic_a2",
        "logistic_b2",
        "logistic_c2",
        "logistic_k",
        "logistic_w_gate",
        "logistic_w_tail",
        "logistic_transition_fraction",
        "logistic_best_start_index",
        "logistic_best_nfev",
        "logistic_k_hit_lower_bound",
        "logistic_k_hit_upper_bound",
        "logistic_w_gate_hit_lower_bound",
        "logistic_w_gate_hit_upper_bound",
        "logistic_w_tail_hit_lower_bound",
        "logistic_w_tail_hit_upper_bound",
        "delta_decoupled_minus_existing",
    ]
    full.write_csv(LOGISTIC_PARAMS_PATH, ordered, fieldnames)


def rewrite_shift_file(rows_by_slug: dict[str, dict], threea_by_slug: dict[str, dict]) -> None:
    parameter_rows = full.build_parameter_shift_rows(rows_by_slug, threea_by_slug)
    fieldnames = [
        "row_type",
        "label",
        "slug",
        "corpus",
        *[f"{param}_pct_change" for param in full.SHARED_PARAMS],
        "head_mean_abs_pct_change",
        "tail_mean_abs_pct_change",
        "notes",
    ]
    full.write_csv(SHIFT_PATH, parameter_rows, fieldnames)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if not PROGRESS_PATH.exists():
        PROGRESS_PATH.write_text("", encoding="utf-8")

    threea_by_slug = full.load_threea_rows()
    rows_by_slug = load_existing_rows()
    completed = set(rows_by_slug)
    total = len(full.SEARCHED_CORPORA)

    for idx, spec in enumerate(full.SEARCHED_CORPORA, start=1):
        slug = spec["slug"]
        if slug in completed:
            continue
        print(f"[{idx}/{total}] decoupled logistic shift fit {slug} ...", flush=True)
        dataset = full.build_zipf_dataset(full.DATA_DIR / spec["filename"])
        fit = full.fit_with_gate(dataset, "logistic")

        existing_bic = full.bic_from_rmse(
            float(threea_by_slug[slug]["reranked_rmse"]),
            full.PARAM_COUNT_EXISTING_3A,
            int(threea_by_slug[slug]["vocab_size"]),
        )
        rows_by_slug[slug] = {
            "slug": slug,
            "corpus": spec["name"],
            "vocab_size": str(int(dataset["unique_words"])),
            "logistic_bic": repr(float(fit["bic"])),
            "logistic_rmse": repr(float(fit["rmse"])),
            "logistic_a1": repr(float(fit["a1"])),
            "logistic_b1": repr(float(fit["b1"])),
            "logistic_c1": repr(float(fit["c1"])),
            "logistic_a2": repr(float(fit["a2"])),
            "logistic_b2": repr(float(fit["b2"])),
            "logistic_c2": repr(float(fit["c2"])),
            "logistic_k": repr(float(fit["k"])),
            "logistic_w_gate": repr(float(fit["w_gate"])),
            "logistic_w_tail": repr(float(fit["w_tail"])),
            "logistic_transition_fraction": repr(float(fit["transition_fraction"])),
            "logistic_best_start_index": str(int(fit["best_start_index"])),
            "logistic_best_nfev": str(int(fit["best_nfev"])),
            "logistic_k_hit_lower_bound": str(bool(fit["k_hit_lower_bound"])),
            "logistic_k_hit_upper_bound": str(bool(fit["k_hit_upper_bound"])),
            "logistic_w_gate_hit_lower_bound": str(bool(fit["w_gate_hit_lower_bound"])),
            "logistic_w_gate_hit_upper_bound": str(bool(fit["w_gate_hit_upper_bound"])),
            "logistic_w_tail_hit_lower_bound": str(bool(fit["w_tail_hit_lower_bound"])),
            "logistic_w_tail_hit_upper_bound": str(bool(fit["w_tail_hit_upper_bound"])),
            "delta_decoupled_minus_existing": repr(float(fit["bic"]) - float(existing_bic)),
        }
        write_logistic_rows(rows_by_slug)
        rewrite_shift_file(rows_by_slug, threea_by_slug)
        with PROGRESS_PATH.open("a", encoding="utf-8") as handle:
            handle.write(f"completed {idx}/{total}: {slug}\n")


if __name__ == "__main__":
    main()
