from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "9_likelihood_frontier"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLES = {
    "hybrid_headtail_splitfit": ROOT / "results" / "zipf_hybrid_headtail_splitfit",
    "hybrid_mechpenalty": ROOT / "results" / "zipf_hybrid_mechpenalty",
    "hybrid_vs_softk_analysis": ROOT / "results" / "zipf_hybrid_vs_softk_analysis",
    "fully_ours_variants": ROOT / "results" / "zipf_fully_ours_variants",
}

SOURCE_FILES = {
    "hybrid_headtail_table": SOURCE_BUNDLES["hybrid_headtail_splitfit"] / "hybrid_headtail_table.csv",
    "hybrid_headtail_summary": SOURCE_BUNDLES["hybrid_headtail_splitfit"] / "summary.json",
    "hybrid_headtail_report": SOURCE_BUNDLES["hybrid_headtail_splitfit"] / "report.md",
    "mechpenalty_table": SOURCE_BUNDLES["hybrid_mechpenalty"] / "hybrid_mechpenalty_table.csv",
    "mechpenalty_summary": SOURCE_BUNDLES["hybrid_mechpenalty"] / "summary.json",
    "mechpenalty_report": SOURCE_BUNDLES["hybrid_mechpenalty"] / "report.md",
    "analysis_table": SOURCE_BUNDLES["hybrid_vs_softk_analysis"] / "hybrid_vs_softk_table.csv",
    "analysis_summary": SOURCE_BUNDLES["hybrid_vs_softk_analysis"] / "summary.json",
    "analysis_report": SOURCE_BUNDLES["hybrid_vs_softk_analysis"] / "report.md",
    "fully_ours_table": SOURCE_BUNDLES["fully_ours_variants"] / "fully_ours_variants_table.csv",
    "fully_ours_summary": SOURCE_BUNDLES["fully_ours_variants"] / "summary.json",
    "fully_ours_report": SOURCE_BUNDLES["fully_ours_variants"] / "report.md",
    "softk_comparator_source_diagnostic": EXPERIMENT_DIR / "softk_comparator_source_diagnostic.csv",
}

OUTPUTS = {
    "hybrid_headtail": OUTPUT_DIR / "hybrid_headtail_per_corpus.csv",
    "fully_ours": OUTPUT_DIR / "fully_ours_variants.csv",
    "mechpenalty": OUTPUT_DIR / "mechanism_penalty_sweep.csv",
    "hybrid_vs_softk_diagnostic": OUTPUT_DIR / "hybrid_vs_softk_diagnostic.csv",
    "hybrid_structure_summary": OUTPUT_DIR / "hybrid_structure_summary.csv",
    "aggregate": OUTPUT_DIR / "aggregate_statistics.csv",
    "manifest": OUTPUT_DIR / "manifest.json",
}

NESTED_PROTOCOL_CONSTANTS = {
    "nested_seam_free_parameter_count": 8,
    "nested_mini_seam_floor": 5,
    "nested_mini_seam_logV_multiplier": 0.5,
    "nested_mini_seam_width": 0.22,
}
