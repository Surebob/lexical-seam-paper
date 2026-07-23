"""Source configuration for experiment 5b migration."""

from pathlib import Path

EXPERIMENT_ID = "5b_lowc_manifold_structure"
RESEARCH_QUESTION = "What is the geometry of the low-c manifold, and how does winner identity change under head-weighted scoring?"

REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = Path(__file__).resolve().parent

SOURCE_BUNDLES = {
    "lowc_manifold": REPO_ROOT / "results" / "zipf_lowc_manifold_analysis",
    "phase_coordinate": REPO_ROOT / "results" / "zipf_phase_coordinate",
    "simulation_recovery": REPO_ROOT / "results" / "zipf_simulation_recovery",
}

OUTPUTS = {
    "lowc_per_corpus": EXPERIMENT_DIR / "outputs" / "lowc_manifold_per_corpus.csv",
    "phase_per_corpus": EXPERIMENT_DIR / "outputs" / "phase_coordinate_per_corpus.csv",
    "aggregate": EXPERIMENT_DIR / "outputs" / "aggregate_statistics.csv",
    "manifest": EXPERIMENT_DIR / "outputs" / "manifest.json",
}

