from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "10a_search_depth_robustness"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLES = {
    "step10_ablation": ROOT / "results" / "zipf_step10_ablation",
    "head_poly_decomposition": ROOT / "results" / "zipf_head_poly_decomposition",
    "head_poly_transfer": ROOT / "results" / "zipf_head_poly_transfer",
}

SOURCE_FILES = {
    "step10_summary": SOURCE_BUNDLES["step10_ablation"] / "summary.json",
    "step10_report": SOURCE_BUNDLES["step10_ablation"] / "report.md",
    "step10_shakespeare_plot": SOURCE_BUNDLES["step10_ablation"] / "shakespeare_step10_ablation.svg",
    "step10_war_and_peace_plot": SOURCE_BUNDLES["step10_ablation"] / "war_and_peace_step10_ablation.svg",
    "decomposition_summary": SOURCE_BUNDLES["head_poly_decomposition"] / "summary.json",
    "decomposition_report": SOURCE_BUNDLES["head_poly_decomposition"] / "report.md",
    "transfer_summary": SOURCE_BUNDLES["head_poly_transfer"] / "summary.json",
    "transfer_report": SOURCE_BUNDLES["head_poly_transfer"] / "report.md",
}

OUTPUTS = {
    "step10": OUTPUT_DIR / "step10_ablation_per_corpus.csv",
    "decomposition": OUTPUT_DIR / "polynomial_decomposition_per_corpus.csv",
    "transfer": OUTPUT_DIR / "poly_transfer.csv",
    "aggregate": OUTPUT_DIR / "aggregate_statistics.csv",
    "manifest": OUTPUT_DIR / "manifest.json",
}

MISSING_TRANSFER_CLAIMS = [
    "War and Peace -> King James Bible",
    "Moby Dick -> Bible",
]
