from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "2a_function_word_ablation"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

FUNCTION_WORD_DIR = ROOT / "results" / "zipf_function_word_test"
FUNCTION_WORD_SUMMARY_JSON = FUNCTION_WORD_DIR / "summary.json"
FUNCTION_WORD_REPORT_MD = FUNCTION_WORD_DIR / "report.md"

ABLATION_LABELS = {
    "remove_top50": "remove_top_50",
    "remove_top100": "remove_top_100",
    "only_top100": "top_100_only",
}
