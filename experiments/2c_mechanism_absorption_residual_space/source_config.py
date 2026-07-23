from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "2c_mechanism_absorption_residual_space"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

BREAKTHROUGH_DIR = ROOT / "results" / "zipf_breakthrough_probe"
BREAKTHROUGH_SUMMARY_JSON = BREAKTHROUGH_DIR / "summary.json"
BREAKTHROUGH_REPORT_MD = BREAKTHROUGH_DIR / "report.md"

SIMULATION_RECOVERY_DIR = ROOT / "results" / "zipf_simulation_recovery"
SIMULATION_RECOVERY_SUMMARY_JSON = SIMULATION_RECOVERY_DIR / "summary.json"

HELPFUL_STEP2_GAIN_THRESHOLD = 0.001
HEAD_WINDOWS = [50, 100, 200]
