from __future__ import annotations

import argparse
import csv
import importlib.util
import statistics
import time
from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
FULL_SCRIPT = ROOT / "phase2_addon" / "s2_decoupled_v3" / "run_s2_v3_full.py"
DEFAULT_OUTDIR = ROOT / "results" / "s2_v3_seed_check_2026-04-18"
TARGET_SLUGS = ["shakespeare", "aesops_fables", "moby_dick"]
SEEDS = [20260415, 20260515, 20260615]
SPREAD_THRESHOLD = 1.0


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


full = load_module(FULL_SCRIPT, "s2_v3_full_for_seed_check")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step A multi-seed reproducibility check for S2 v3.")
    parser.add_argument("--n-starts", type=int, default=100)
    parser.add_argument("--max-nfev", type=int, default=12000)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_rows(path: Path, key_fields: tuple[str, ...]) -> dict[tuple[str, ...], dict]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {tuple(row[field] for field in key_fields): row for row in csv.DictReader(handle)}


def summarize(seed_rows: list[dict]) -> tuple[list[dict], list[dict], bool]:
    grouped: dict[tuple[str, str], list[dict]] = {}
    for row in seed_rows:
        grouped.setdefault((row["slug"], row["gate"]), []).append(row)

    dispersion_rows = []
    all_pass = True
    for (slug, gate), rows in sorted(grouped.items()):
        bics = [float(row["bic"]) for row in rows]
        seed_count = len(rows)
        spread = max(bics) - min(bics)
        passes = seed_count == len(SEEDS) and spread < SPREAD_THRESHOLD
        if not passes:
            all_pass = False
        dispersion_rows.append(
            {
                "slug": slug,
                "corpus": rows[0]["corpus"],
                "gate": gate,
                "seed_count": str(seed_count),
                "min_bic": repr(float(min(bics))),
                "median_bic": repr(float(statistics.median(bics))),
                "max_bic": repr(float(max(bics))),
                "bic_spread_across_seeds": repr(float(spread)),
                "passes_step_a": str(passes),
                "notes": "Pass criterion is max_bic - min_bic < 1.0 across the three seeds.",
            }
        )

    aggregate_rows = []
    pass_count = sum(row["passes_step_a"] == "True" for row in dispersion_rows)
    aggregate_rows.append(
        {
            "metric_name": "gate_corpus_pass_count",
            "value": str(pass_count),
            "notes": "Count of gate/corpus combinations passing the <1.0 BIC spread criterion.",
        }
    )
    aggregate_rows.append(
        {
            "metric_name": "gate_corpus_total",
            "value": str(len(dispersion_rows)),
            "notes": "Total gate/corpus combinations in the seed check.",
        }
    )
    aggregate_rows.append(
        {
            "metric_name": "step_a_overall_pass",
            "value": str(all_pass),
            "notes": "True only if every gate on every target corpus has BIC spread < 1.0 across seeds.",
        }
    )
    if dispersion_rows:
        aggregate_rows.append(
            {
                "metric_name": "max_bic_spread_across_all_checks",
                "value": repr(float(max(float(row["bic_spread_across_seeds"]) for row in dispersion_rows))),
                "notes": "Worst BIC spread among all gate/corpus combinations.",
            }
        )
    return dispersion_rows, aggregate_rows, all_pass


def write_summary(path: Path, n_starts: int, dispersion_rows: list[dict], aggregate_rows: list[dict]) -> None:
    agg = {row["metric_name"]: row["value"] for row in aggregate_rows}
    failing = [row for row in dispersion_rows if row["passes_step_a"] != "True"]
    lines = [
        "# S2 v3 Step A Multi-Seed Reproducibility Check",
        "",
        f"- corpora: `{', '.join(TARGET_SLUGS)}`",
        f"- seeds: `{', '.join(str(seed) for seed in SEEDS)}`",
        f"- starts per run: `{n_start(s=n_starts)}`",
        f"- overall pass: `{agg.get('step_a_overall_pass', 'False')}`",
        f"- gate/corpus pass count: `{agg.get('gate_corpus_pass_count', '0')}/{agg.get('gate_corpus_total', '0')}`",
        f"- worst BIC spread: `{agg.get('max_bic_spread_across_all_checks', 'n/a')}`",
        "",
    ]
    if failing:
        lines.extend(["## Failing gate/corpus checks", ""])
        for row in failing:
            lines.append(
                f"- {row['corpus']} / {row['gate']}: spread `{float(row['bic_spread_across_seeds']):.6f}`"
            )
    else:
        lines.extend(["## Result", "", "- All gate/corpus combinations passed the <1.0 BIC spread criterion."])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def n_start(s: int) -> str:
    return f"{s}"


def main() -> None:
    args = parse_args()
    outdir = args.outdir / f"{args.n_starts}_starts"
    outdir.mkdir(parents=True, exist_ok=True)

    per_seed_path = outdir / "s2_v3_seed_check_per_seed.csv"
    dispersion_path = outdir / "s2_v3_seed_check_dispersion.csv"
    aggregate_path = outdir / "s2_v3_seed_check_aggregate_statistics.csv"
    summary_path = outdir / "s2_v3_seed_check_summary_report.md"
    progress_path = outdir / "progress.log"

    if not progress_path.exists():
        progress_path.write_text("", encoding="utf-8")

    existing = load_rows(per_seed_path, ("slug", "gate", "seed"))
    rows_by_key = dict(existing)

    corpus_specs = [spec for spec in full.SEARCHED_CORPORA if spec["slug"] in TARGET_SLUGS]
    total = len(corpus_specs) * len(SEEDS) * len(full.GATE_FUNCS)

    done_index = 0
    for spec in corpus_specs:
        dataset = None
        for seed in SEEDS:
            for gate in full.GATE_FUNCS:
                done_index += 1
                key = (spec["slug"], gate, str(seed))
                if key in rows_by_key:
                    continue
                if dataset is None:
                    dataset = full.build_zipf_dataset(full.DATA_DIR / spec["filename"])
                print(
                    f"[{done_index}/{total}] seed-check fit {spec['slug']} / {gate} / seed={seed} ...",
                    flush=True,
                )
                original_seed = full.BASE_SEED
                start_time = time.time()
                try:
                    full.BASE_SEED = seed
                    fit = full.fit_with_gate(dataset, gate, n_starts=args.n_starts, max_nfev=args.max_nfev)
                finally:
                    full.BASE_SEED = original_seed
                runtime_sec = time.time() - start_time
                row = {
                    "slug": spec["slug"],
                    "corpus": spec["name"],
                    "gate": gate,
                    "seed": str(seed),
                    "bic": repr(float(fit["bic"])),
                    "rmse": repr(float(fit["rmse"])),
                    "k": repr(float(fit["k"])),
                    "w_gate": repr(float(fit["w_gate"])),
                    "w_tail": repr(float(fit["w_tail"])),
                    "best_start_index": str(int(fit["best_start_index"])),
                    "best_nfev": str(int(fit["best_nfev"])),
                    "runtime_sec": repr(float(runtime_sec)),
                    "k_hit_lower_bound": str(bool(fit["k_hit_lower_bound"])),
                    "k_hit_upper_bound": str(bool(fit["k_hit_upper_bound"])),
                    "w_gate_hit_lower_bound": str(bool(fit["w_gate_hit_lower_bound"])),
                    "w_gate_hit_upper_bound": str(bool(fit["w_gate_hit_upper_bound"])),
                    "w_tail_hit_lower_bound": str(bool(fit["w_tail_hit_lower_bound"])),
                    "w_tail_hit_upper_bound": str(bool(fit["w_tail_hit_upper_bound"])),
                }
                rows_by_key[key] = row
                seed_rows = [rows_by_key[k] for k in sorted(rows_by_key)]
                dispersion_rows, aggregate_rows, _ = summarize(seed_rows)
                write_csv(
                    per_seed_path,
                    seed_rows,
                    [
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
                    ],
                )
                write_csv(
                    dispersion_path,
                    dispersion_rows,
                    [
                        "slug",
                        "corpus",
                        "gate",
                        "seed_count",
                        "min_bic",
                        "median_bic",
                        "max_bic",
                        "bic_spread_across_seeds",
                        "passes_step_a",
                        "notes",
                    ],
                )
                write_csv(aggregate_path, aggregate_rows, ["metric_name", "value", "notes"])
                write_summary(summary_path, args.n_starts, dispersion_rows, aggregate_rows)
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(f"completed {spec['slug']} / {gate} / seed={seed}\n")


if __name__ == "__main__":
    main()
