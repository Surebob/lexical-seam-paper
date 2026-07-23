from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "1b_non_linguistic_control"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

CITY_FULL_DIR = ROOT / "results" / "zipf_enriched_city_populations"
CITY_STEP2_DIR = ROOT / "results" / "zipf_enriched_city_populations_step2_only"

CITY_FULL_SUMMARY = CITY_FULL_DIR / "summary.json"
CITY_STEP2_SUMMARY = CITY_STEP2_DIR / "summary.json"

SOURCE_BUNDLES = [CITY_FULL_DIR, CITY_STEP2_DIR]

