import json
import subprocess
from collections import Counter
from pathlib import Path

import importlib.util


ROOT = Path("/Volumes/External2TB/emlexperiment")
ZIPF_EXPERIMENT_PATH = ROOT / "eml_zipf_experiment.py"
SEARCH_SCRIPT = ROOT / "eml_zipf_enriched_search.py"
SOURCE_CORPUS = ROOT / "data" / "zipf" / "pg100.txt"
DERIVED_DIR = ROOT / "data" / "zipf" / "function_word_test"
OUTDIR = ROOT / "results" / "zipf_function_word_test"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


zipf = load_module(ZIPF_EXPERIMENT_PATH, "zipf_function_word_test_zipf")


def write_tokens(path: Path, tokens):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for idx, token in enumerate(tokens):
            if idx:
                handle.write(" ")
            handle.write(token)


def run_search(corpus_path: Path, outdir: Path):
    cmd = [
        "python3",
        str(SEARCH_SCRIPT),
        "--corpus-path",
        str(corpus_path),
        "--beam-width",
        "50",
        "--max-steps",
        "10",
        "--keep-all-until-step",
        "2",
        "--sample-points",
        "0",
        "--outdir",
        str(outdir),
    ]
    subprocess.run(cmd, check=True)


def summarize_case(name: str, dataset, summary_path: Path):
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    step2 = next(item for item in payload["zm_search"]["step_summary"] if item["step"] == 2)
    step2_top = step2["top_candidates"][0]
    return {
        "case": name,
        "token_count": dataset["token_count"],
        "unique_words": dataset["unique_words"],
        "top_words_preview": dataset["ranked"][:20],
        "zm_a": payload["zm_baseline"]["a"],
        "zm_b": payload["zm_baseline"]["b"],
        "zm_c": payload["zm_baseline"]["c"],
        "zm_rmse": payload["zm_baseline"]["rmse_full"],
        "step2_winner": step2_top["expr"],
        "step2_math": step2_top["math"],
        "step2_rmse": step2_top["rmse"],
        "step2_helps": step2_top["rmse"] < payload["zm_baseline"]["rmse_full"],
        "best_overall": payload["zm_search"]["best"],
        "summary_path": str(summary_path),
    }


def write_report(results, top50, top100):
    lines = [
        "# Shakespeare Function Word Test",
        "",
        "Top-50 proxy function words are defined operationally as the 50 most frequent Shakespeare word types.",
        "Top-100 proxy function words are defined operationally as the 100 most frequent Shakespeare word types.",
        "",
        "## Top 50",
        "",
        "`" + ", ".join(top50) + "`",
        "",
        "## Top 100",
        "",
        "`" + ", ".join(top100) + "`",
        "",
        "## Results",
        "",
    ]
    for item in results:
        lines.extend(
            [
                f"### {item['case']}",
                "",
                f"- token count: `{item['token_count']}`",
                f"- vocab size: `{item['unique_words']}`",
                f"- ZM params: `a={item['zm_a']:.12f}`, `b={item['zm_b']:.12f}`, `c={item['zm_c']:.12f}`",
                f"- ZM RMSE: `{item['zm_rmse']:.12f}`",
                f"- step-2 winner: `{item['step2_winner']}`",
                f"- step-2 math: `{item['step2_math']}`",
                f"- step-2 RMSE: `{item['step2_rmse']:.12f}`",
                f"- step-2 helps ZM: `{item['step2_helps']}`",
                f"- best overall ZM-search result: `{item['best_overall']['expr']}` at RMSE `{item['best_overall']['composite_rmse']:.12f}`",
                "",
            ]
        )
    (OUTDIR / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)

    base = zipf.build_zipf_dataset(SOURCE_CORPUS)
    top100 = [word for word, _ in base["ranked"][:100]]
    top50 = top100[:50]
    top50_set = set(top50)
    top100_set = set(top100)

    cases = [
        (
            "remove_top50",
            [token for token in base["tokens"] if token not in top50_set],
            DERIVED_DIR / "shakespeare_remove_top50.txt",
            OUTDIR / "remove_top50",
        ),
        (
            "remove_top100",
            [token for token in base["tokens"] if token not in top100_set],
            DERIVED_DIR / "shakespeare_remove_top100.txt",
            OUTDIR / "remove_top100",
        ),
        (
            "only_top100",
            [token for token in base["tokens"] if token in top100_set],
            DERIVED_DIR / "shakespeare_only_top100.txt",
            OUTDIR / "only_top100",
        ),
    ]

    results = []
    for name, tokens, corpus_path, case_outdir in cases:
        write_tokens(corpus_path, tokens)
        dataset = zipf.build_zipf_dataset(corpus_path)
        run_search(corpus_path, case_outdir)
        results.append(summarize_case(name, dataset, case_outdir / "summary.json"))

    payload = {
        "source_corpus": str(SOURCE_CORPUS),
        "top50_words": top50,
        "top100_words": top100,
        "results": results,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(results, top50, top100)
    print(f"Saved {OUTDIR / 'summary.json'}")
    print(f"Saved {OUTDIR / 'report.md'}")


if __name__ == "__main__":
    main()
