from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "6_multilingual_extension"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

MULTILANG_ROMANCE_DIR = ROOT / "results" / "zipf_multilang_romance"
MULTILANG_ROMANCE_SUMMARY_JSON = MULTILANG_ROMANCE_DIR / "summary.json"
MULTILANG_ROMANCE_TABLE_CSV = MULTILANG_ROMANCE_DIR / "multilang_table.csv"
MULTILANG_VERIFY_DIR = ROOT / "results" / "zipf_multilang_verify"
MULTILANG_VERIFY_SUMMARY_JSON = MULTILANG_VERIFY_DIR / "summary.json"

NEAR_HALF_LOW = 0.42
NEAR_HALF_HIGH = 0.62
VERY_LOW_C_THRESHOLD = 0.01
LOW_C_THRESHOLD = 5.0
