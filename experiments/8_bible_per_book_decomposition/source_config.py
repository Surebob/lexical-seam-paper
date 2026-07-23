from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
EXPERIMENT_DIR = ROOT / "experiments" / "8_bible_per_book_decomposition"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ARCHIVE_DIR = EXPERIMENT_DIR / "archive"

SOURCE_BUNDLE = ROOT / "results" / "zipf_angle6_bible_books"
SOURCE_TABLE = SOURCE_BUNDLE / "bible_books_table.csv"
SOURCE_SUMMARY = SOURCE_BUNDLE / "summary.json"
SOURCE_REPORT = SOURCE_BUNDLE / "report.md"
SOURCE_FIGURE = SOURCE_BUNDLE / "step2_gains_by_book.png"

PMF_VARIANT_TABLE = ROOT / "experiments" / "7a_canonical_pmf_family" / "outputs" / "splitfit" / "pmf_variant_per_corpus.csv"
