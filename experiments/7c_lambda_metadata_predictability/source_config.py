from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "7c_lambda_metadata_predictability"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLE = ROOT / "results" / "zipf_angle2_lambda_metadata"
SOURCE_TABLE = SOURCE_BUNDLE / "lambda_metadata_table.csv"
SOURCE_SUMMARY = SOURCE_BUNDLE / "summary.json"
SOURCE_REPORT = SOURCE_BUNDLE / "report.md"

PREDICTOR_ORDER = [
    "structure_code",
    "log_unit_count",
    "log_author_count",
    "log_era_span_years",
    "heterogeneity_score",
    "token_count",
    "vocab_size",
]

MANUSCRIPT_METADATA_PREDICTORS = [
    "log_unit_count",
    "log_author_count",
    "log_era_span_years",
    "heterogeneity_score",
]
