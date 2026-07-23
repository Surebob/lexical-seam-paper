"""Source configuration for J1 paper-scope rollup."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

UPSTREAMS = {
    "1a_aggregate": ROOT / "experiments/1a_per_corpus_enriched_search/outputs/aggregate_statistics.csv",
    "3a_aggregate": ROOT / "experiments/3a_smooth_two_regime_fits/outputs/aggregate_statistics.csv",
    "6_aggregate": ROOT / "experiments/6_multilingual_extension/outputs/aggregate_statistics.csv",
    "6_table3": ROOT / "experiments/6_multilingual_extension/outputs/table3_multilang.csv",
}

OUTPUTS = {
    "summary": ROOT / "joins/J1_paper_scope_rollup/outputs/paper_scope_summary.csv",
    "manifest": ROOT / "joins/J1_paper_scope_rollup/outputs/manifest.json",
}
