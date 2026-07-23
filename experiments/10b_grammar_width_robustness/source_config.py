from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "10b_grammar_width_robustness"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLES = {
    "widened_grammar_diagnostic": ROOT / "results" / "zipf_widened_grammar_diagnostic",
    "widened_grammar_extended": ROOT / "results" / "zipf_widened_grammar_extended",
}

SOURCE_FILES = {
    "diagnostic_summary": SOURCE_BUNDLES["widened_grammar_diagnostic"] / "summary.json",
    "diagnostic_report": SOURCE_BUNDLES["widened_grammar_diagnostic"] / "report.md",
    "extended_summary": SOURCE_BUNDLES["widened_grammar_extended"] / "summary.json",
    "extended_table": SOURCE_BUNDLES["widened_grammar_extended"] / "widened_lowc_summary.csv",
    "extended_report": SOURCE_BUNDLES["widened_grammar_extended"] / "report.md",
}

OUTPUTS = {
    "diagnostic": OUTPUT_DIR / "widened_diagnostic_per_mode.csv",
    "diagnostic_steps": OUTPUT_DIR / "widened_diagnostic_step_winners.csv",
    "diagnostic_top10": OUTPUT_DIR / "widened_diagnostic_top10.csv",
    "lowc_summary": OUTPUT_DIR / "widened_lowc_manifold_summary.csv",
    "aggregate": OUTPUT_DIR / "aggregate_statistics.csv",
    "manifest": OUTPUT_DIR / "manifest.json",
}
