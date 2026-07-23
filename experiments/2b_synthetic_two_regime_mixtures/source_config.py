from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "2b_synthetic_two_regime_mixtures"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SYNTHETIC_MIXTURE_DIR = ROOT / "results" / "zipf_synthetic_mixture"
SYNTHETIC_MIXTURE_SUMMARY_JSON = SYNTHETIC_MIXTURE_DIR / "summary.json"
SYNTHETIC_MIXTURE_REPORT_MD = SYNTHETIC_MIXTURE_DIR / "report.md"

BREAKTHROUGH_DIR = ROOT / "results" / "zipf_breakthrough_probe"
BREAKTHROUGH_SUMMARY_JSON = BREAKTHROUGH_DIR / "summary.json"

MIXTURE_METADATA = {
    "mixture_alpha1.5_plus_1.3": {
        "mixture_id": "small_gap",
        "alpha_1": 1.5,
        "alpha_2": 1.3,
        "weight_1": 0.5,
        "claimed_by_v4": True,
    },
    "mixture_alpha1.5_plus_0.8": {
        "mixture_id": "medium_gap",
        "alpha_1": 1.5,
        "alpha_2": 0.8,
        "weight_1": 0.5,
        "claimed_by_v4": True,
    },
    "mixture_alpha1.5_plus_0.5": {
        "mixture_id": "large_gap",
        "alpha_1": 1.5,
        "alpha_2": 0.5,
        "weight_1": 0.5,
        "claimed_by_v4": True,
    },
    "control_two_copies_alpha0.8": {
        "mixture_id": "historical_same_exponent_control_alpha0p8",
        "alpha_1": 0.8,
        "alpha_2": 0.8,
        "weight_1": 0.5,
        "claimed_by_v4": False,
    },
}
