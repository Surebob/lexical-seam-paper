"""Source configuration for J2 scaling comparison."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

UPSTREAMS = {
    "3e_erf_fits": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/per_fit_results.csv",
    "3e_per_corpus": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/per_corpus_results.csv",
    "3d_pos": ROOT / "experiments/3d_pos_crossover_scaling/outputs/aggregate_statistics.csv",
}

OUTPUTS = {
    "summary": ROOT / "joins/J2_scaling_comparison/outputs/scaling_comparison_summary.csv",
    "manifest": ROOT / "joins/J2_scaling_comparison/outputs/manifest.json",
}
