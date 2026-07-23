from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "7e_hierarchical_k_pooling"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLE = ROOT / "results" / "zipf_angle4_hierk"
SOURCE_TABLE = SOURCE_BUNDLE / "hierk_head_window_table.csv"
SOURCE_SUMMARY = SOURCE_BUNDLE / "summary.json"
SOURCE_REPORT = SOURCE_BUNDLE / "report.md"

SOFTK_CANONICAL_TABLE = ROOT / "experiments" / "7a_canonical_pmf_family" / "outputs" / "splitfit" / "table4_fourway.csv"

CUTOFF_ORDER = ["top50", "top100", "top200", "top500", "top1000", "full"]
