from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "3b_model_family_bic_comparison"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

BIC_COMPARISON_DIR = ROOT / "results" / "zipf_bic_comparison"
MOEZIPF_COMPARISON_DIR = ROOT / "results" / "zipf_moezipf_comparison"
CONTINUOUS_PIECEWISE_DIR = ROOT / "results" / "zipf_continuous_piecewise"
SQRT_V_DIR = ROOT / "results" / "zipf_sqrt_v_all_corpora"
BIC_LANDSCAPE_DIR = ROOT / "results" / "zipf_bic_landscape"

BIC_COMPARISON_SUMMARY = BIC_COMPARISON_DIR / "summary.json"
MOEZIPF_COMPARISON_SUMMARY = MOEZIPF_COMPARISON_DIR / "summary.json"
CONTINUOUS_PIECEWISE_SUMMARY = CONTINUOUS_PIECEWISE_DIR / "summary.json"
SQRT_V_SUMMARY = SQRT_V_DIR / "summary.json"

SOURCE_BUNDLES = [
    BIC_COMPARISON_DIR,
    MOEZIPF_COMPARISON_DIR,
    CONTINUOUS_PIECEWISE_DIR,
    SQRT_V_DIR,
    BIC_LANDSCAPE_DIR,
]

