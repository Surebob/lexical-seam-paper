"""Source configuration for experiment 4 migration.

Canonical Phase 2 model status:
    blocked, because the available seam-theory bundles derive the historical
    coupled-logistic model, while the current canonical smooth model is
    decoupled-erf.

The coupled-logistic bundles are legitimate prior experimental outputs and are
preserved under archive/legacy_coupled_logistic/.
"""

from pathlib import Path


EXPERIMENT_ID = "4_analytic_seam_theory"
CANONICAL_MODEL = "decoupled_erf"
CANONICAL_STATUS = "blocked_missing_decoupled_erf_seam_theory"

REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = Path(__file__).resolve().parent

LEGACY_MODEL = "coupled_logistic"
LEGACY_ARCHIVE_DIR = EXPERIMENT_DIR / "archive" / "legacy_coupled_logistic"

LEGACY_SOURCE_BUNDLES = {
    "sign_theory": REPO_ROOT / "results" / "zipf_seam_sign_theory",
    "projection_theory": REPO_ROOT / "results" / "zipf_seam_projection_theory",
    "second_order_theory": REPO_ROOT / "results" / "zipf_seam_second_order_theory",
}

CANONICAL_MISSING_SOURCES = [
    "decoupled-erf raw seam sign-law outputs",
    "decoupled-erf tangent-space projection outputs",
    "decoupled-erf second-order seam correction outputs",
]

