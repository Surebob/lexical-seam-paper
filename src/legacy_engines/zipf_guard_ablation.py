import json
import subprocess
from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
SEARCH_SCRIPT = ROOT / "eml_zipf_enriched_search.py"
OUTDIR = ROOT / "results" / "zipf_guard_ablation"

CORPORA = [
    {
        "slug": "shakespeare",
        "name": "Complete Works of Shakespeare",
        "corpus_path": ROOT / "data" / "zipf" / "pg100.txt",
        "baseline_summary": ROOT / "results" / "zipf_enriched_search_full" / "summary.json",
    },
    {
        "slug": "principia_ethica",
        "name": "Principia Ethica",
        "corpus_path": ROOT / "data" / "zipf" / "pg53430.txt",
        "baseline_summary": ROOT / "results" / "zipf_enriched_principia_ethica_full" / "summary.json",
    },
    {
        "slug": "ulysses",
        "name": "Ulysses",
        "corpus_path": ROOT / "data" / "zipf" / "pg4300.txt",
        "baseline_summary": ROOT / "results" / "zipf_enriched_ulysses_full" / "summary.json",
    },
]

PROFILES = [
    {"slug": "tight", "exp_clamp": 20.0, "value_abs_limit": 1.0e4},
    {"slug": "loose", "exp_clamp": 60.0, "value_abs_limit": 1.0e9},
]


def load_summary(path: Path):
    obj = json.loads(path.read_text(encoding="utf-8"))
    step2 = obj["zm_search"]["step_summary"][1]["top_candidates"][0]
    return {
        "summary_path": str(path),
        "zm_c": obj["zm_baseline"]["c"],
        "zm_rmse": obj["zm_baseline"]["rmse_full"],
        "step2_expr": step2["expr"],
        "step2_math": step2["math"],
        "step2_rmse": step2["rmse"],
        "best_expr": obj["zm_search"]["best"]["expr"],
        "best_math": obj["zm_search"]["best"]["math"],
        "best_rmse": obj["zm_search"]["best"]["composite_rmse"],
    }


def run_profile(corpus: dict, profile: dict):
    outdir = OUTDIR / f"{corpus['slug']}__{profile['slug']}"
    cmd = [
        "python3",
        str(SEARCH_SCRIPT),
        "--corpus-path",
        str(corpus["corpus_path"]),
        "--beam-width",
        "50",
        "--max-steps",
        "10",
        "--keep-all-until-step",
        "2",
        "--sample-points",
        "0",
        "--exp-clamp",
        str(profile["exp_clamp"]),
        "--value-abs-limit",
        str(profile["value_abs_limit"]),
        "--outdir",
        str(outdir),
    ]
    subprocess.run(cmd, check=True, cwd=ROOT)
    return load_summary(outdir / "summary.json") | {
        "profile": profile["slug"],
        "exp_clamp": profile["exp_clamp"],
        "value_abs_limit": profile["value_abs_limit"],
        "outdir": str(outdir),
    }


def write_report(rows):
    lines = [
        "# Zipf Guard Ablation",
        "",
        "- Only two search guards were varied: `EXP_CLAMP` and `VALUE_ABS_LIMIT`.",
        "- Baselines are the already-completed original runs with `EXP_CLAMP=30` and `VALUE_ABS_LIMIT=1e6`.",
        "",
    ]
    for corpus in CORPORA:
        baseline = next(
            row for row in rows if row["corpus_slug"] == corpus["slug"] and row["profile"] == "baseline"
        )
        lines.extend(
            [
                f"## {corpus['name']}",
                "",
                "| Profile | EXP_CLAMP | VALUE_ABS_LIMIT | step-2 expr | step-2 RMSE | best overall RMSE |",
                "| --- | ---: | ---: | --- | ---: | ---: |",
            ]
        )
        corpus_rows = [row for row in rows if row["corpus_slug"] == corpus["slug"]]
        for row in corpus_rows:
            lines.append(
                f"| {row['profile']} | {row['exp_clamp']} | {row['value_abs_limit']:.0e} | "
                f"`{row['step2_expr']}` | {row['step2_rmse']:.12f} | {row['best_rmse']:.12f} |"
            )
        lines.extend(
            [
                "",
                f"- baseline step-2: `{baseline['step2_expr']}`",
                "",
            ]
        )
    (OUTDIR / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for corpus in CORPORA:
        baseline = load_summary(corpus["baseline_summary"]) | {
            "corpus_slug": corpus["slug"],
            "corpus_name": corpus["name"],
            "profile": "baseline",
            "exp_clamp": 30.0,
            "value_abs_limit": 1.0e6,
            "outdir": str(corpus["baseline_summary"].parent),
        }
        rows.append(baseline)
        for profile in PROFILES:
            result = run_profile(corpus, profile)
            rows.append(
                result
                | {
                    "corpus_slug": corpus["slug"],
                    "corpus_name": corpus["name"],
                }
            )

    (OUTDIR / "summary.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    write_report(rows)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
