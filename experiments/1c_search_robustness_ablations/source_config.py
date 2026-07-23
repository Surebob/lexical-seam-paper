from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "1c_search_robustness_ablations"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

BOUNDARY_DIR = ROOT / "results" / "zipf_boundary_ablation"
GUARD_DIR = ROOT / "results" / "zipf_guard_ablation"
WLS_DIR = ROOT / "results" / "zipf_wls_test"
ONE_A_DIR = ROOT / "experiments" / "1a_per_corpus_enriched_search"

BOUNDARY_SUMMARY = BOUNDARY_DIR / "summary.json"
GUARD_SUMMARY = GUARD_DIR / "summary.json"
WLS_SUMMARY = WLS_DIR / "summary.json"
ONE_A_TABLE = ONE_A_DIR / "outputs" / "table1_per_corpus.csv"

SOURCE_BUNDLES = [BOUNDARY_DIR, GUARD_DIR, WLS_DIR]

