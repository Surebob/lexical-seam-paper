from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "1d_winner_vs_euclidean_gap_analysis"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

GAP_VERIFY_DIR = ROOT / "results" / "zipf_english_gap_verify"
GAP_SUMMARY = GAP_VERIFY_DIR / "summary.json"
GAP_TABLE = GAP_VERIFY_DIR / "english_gap_table.csv"

SOURCE_BUNDLES = [GAP_VERIFY_DIR]

