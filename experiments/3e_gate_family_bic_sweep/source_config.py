"""Source configuration for experiment 3e gate-family BIC sweep migration."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

SOURCES = {
    "real_gate_sweep_dir": ROOT / "results/s2_v3_windows_full_outputs_2026-04-18",
    "synthetic_recovery_dir": ROOT / "results/s2_synthetic_recovery_outputs_workers32_2026-04-18/outputs",
    "logistic_parameter_source": ROOT
    / "results/windows_s2_decoupled_v3_results_2026-04-18/outputs/s2_v3_logistic_local_params.csv",
}

PROVENANCE_SCRIPTS = {
    "mac_full_sweep_script": ROOT / "phase2_addon/s2_decoupled_v3/run_s2_v3_full.py",
    "windows_full_sweep_script": ROOT / "s2_v3_windows_port/run_s2_v3_windows.py",
    "synthetic_recovery_script": ROOT
    / "phase2_addon/s2_logistic_synthetic_gate_recovery/run_synthetic_gate_recovery.py",
}

OUTPUTS = {
    "per_fit": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/per_fit_results.csv",
    "per_corpus": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/per_corpus_results.csv",
    "aggregate": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/aggregate_statistics.csv",
    "synthetic": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/synthetic_recovery.csv",
    "synthetic_aggregate": ROOT
    / "experiments/3e_gate_family_bic_sweep/outputs/synthetic_recovery_aggregate.csv",
    "manifest": ROOT / "experiments/3e_gate_family_bic_sweep/outputs/manifest.json",
}

INDEPENDENT_GATES = ["logistic", "erf", "algebraic", "arctan"]
ORDERED_GATES = ["logistic", "tanh", "erf", "algebraic", "arctan"]
