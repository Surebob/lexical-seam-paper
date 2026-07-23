"""Source configuration for experiment 5a migration."""

from pathlib import Path

EXPERIMENT_ID = "5a_simulation_recovery_from_smooth_models"
RESEARCH_QUESTION = "Does synthetic data sampled from fitted smooth models reproduce the empirical symbolic residual signal?"

REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = Path(__file__).resolve().parent
SOURCE_BUNDLE = REPO_ROOT / "results" / "zipf_simulation_recovery"

SOURCE_FILES = {
    "summary": SOURCE_BUNDLE / "summary.json",
    "table": SOURCE_BUNDLE / "simulation_recovery_table.csv",
    "report": SOURCE_BUNDLE / "report.md",
}

OUTPUTS = {
    "per_corpus": EXPERIMENT_DIR / "outputs" / "simulation_recovery_per_corpus.csv",
    "aggregate": EXPERIMENT_DIR / "outputs" / "aggregate_statistics.csv",
    "manifest": EXPERIMENT_DIR / "outputs" / "manifest.json",
}

