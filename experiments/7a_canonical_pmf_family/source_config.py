"""Source configuration for experiment 7a migration."""

from pathlib import Path

EXPERIMENT_ID = "7a_canonical_pmf_family"
RESEARCH_QUESTION = "What is the canonical held-out tradeoff among PMF variants under the 80/20 binomial split protocol?"

TRAIN_FRACTION = 0.8
TEST_FRACTION = 0.2
SPLIT_TYPE = "binomial_per_token"
SOFTK_PARAMETER_COUNT = 6

REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = Path(__file__).resolve().parent

SOURCE_BUNDLES = {
    "pmf_free": REPO_ROOT / "results" / "zipf_seam_mandelbrot_pmf",
    "regularized": REPO_ROOT / "results" / "zipf_seam_mandelbrot_regularized",
    "softk_legacy": REPO_ROOT / "results" / "zipf_seam_mandelbrot_softk",
    "softk_splitfit": REPO_ROOT / "results" / "zipf_seam_mandelbrot_softk_splitfit",
    "softkw": REPO_ROOT / "results" / "zipf_seam_mandelbrot_softkw",
    "v4_verification": REPO_ROOT / "results" / "zipf_v4_verification",
}

CANONICAL_TABLE4_SOURCE = SOURCE_BUNDLES["v4_verification"] / "table_a_fourway_pmf.csv"

OUTPUTS = {
    "table4": EXPERIMENT_DIR / "outputs" / "splitfit" / "table4_fourway.csv",
    "table4_provenance": EXPERIMENT_DIR / "outputs" / "splitfit" / "table4_provenance.csv",
    "softk_diagnostic": EXPERIMENT_DIR / "outputs" / "splitfit" / "softk_source_diagnostic.csv",
    "variant_per_corpus": EXPERIMENT_DIR / "outputs" / "splitfit" / "pmf_variant_per_corpus.csv",
    "aggregate": EXPERIMENT_DIR / "outputs" / "splitfit" / "aggregate_statistics.csv",
    "fullrefit_per_corpus": EXPERIMENT_DIR / "outputs" / "fullrefit" / "fourway_per_corpus.csv",
    "fullrefit_aggregate": EXPERIMENT_DIR / "outputs" / "fullrefit" / "aggregate_statistics.csv",
    "manifest": EXPERIMENT_DIR / "manifest.json",
}

