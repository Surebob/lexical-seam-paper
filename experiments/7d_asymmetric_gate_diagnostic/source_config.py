from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "7d_asymmetric_gate_diagnostic"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLE = ROOT / "results" / "zipf_angle3_asymmetric_gate_splitfit"
SOURCE_TABLE = SOURCE_BUNDLE / "asymmetric_gate_table.csv"
SOURCE_SUMMARY = SOURCE_BUNDLE / "summary.json"
SOURCE_REPORT = SOURCE_BUNDLE / "report.md"

CUTOFF_ORDER = ["top50", "top100", "top200", "top500", "top1000", "full"]
MODEL_ORDER = ["zipf", "zm", "moe", "softk", "asym"]
