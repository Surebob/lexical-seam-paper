from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "3d_pos_crossover_scaling"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

POS_ALL_DIR = ROOT / "results" / "zipf_pos_all_corpora"
POS_ALL_POINTS_CSV = POS_ALL_DIR / "pos_all_corpora_points.csv"
POS_ALL_SUMMARY_JSON = POS_ALL_DIR / "summary.json"
POS_ALL_REPORT_MD = POS_ALL_DIR / "report.md"

POS_MANUAL_DIR = ROOT / "results" / "zipf_pos_manual_v2"
POS_MANUAL_POINTS_CSV = POS_MANUAL_DIR / "manual_alpha_points.csv"
POS_MANUAL_SUMMARY_JSON = POS_MANUAL_DIR / "summary.json"
POS_MANUAL_REPORT_MD = POS_MANUAL_DIR / "report.md"

TOP_POS_WINDOW = 500
CROSSOVER_FRACTION = 0.50
MANUAL_VALIDATION_WINDOW = 300

CI_METHOD = "Student-t 95% CI on the forced alpha using historical curve_fit covariance SE and df = n - 1"
PVALUE_METHOD = "Two-sided one-sample t-test on per-corpus alpha values against 0.5"
