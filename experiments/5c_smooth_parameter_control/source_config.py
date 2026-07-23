"""Source configuration for experiment 5c migration."""

from pathlib import Path

EXPERIMENT_ID = "5c_smooth_parameter_control"
CANONICAL_MODEL = "decoupled_erf"
CANONICAL_STATUS = "blocked_missing_decoupled_erf_parameter_sweep"

REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = Path(__file__).resolve().parent

LEGACY_MODEL = "coupled_logistic_parameter_sweep"
LEGACY_ARCHIVE_DIR = EXPERIMENT_DIR / "archive" / "legacy_coupled_logistic"
LEGACY_SOURCE_BUNDLE = REPO_ROOT / "results" / "zipf_smooth_parameter_sweep"

CANONICAL_MISSING_SOURCES = [
    "decoupled-erf low-c smooth-parameter sweep rows",
    "decoupled-erf parameter-control correlations for w_gate and w_tail",
    "decoupled-erf synthetic winner-count aggregates under the parameter sweep",
]

