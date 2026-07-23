from pathlib import Path


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTPUT_DIR = ROOT / "experiments" / "3c_kstat_scaling" / "outputs"
SMOOTH_FIT_CSV = ROOT / "experiments" / "3a_smooth_two_regime_fits" / "outputs" / "smooth_fit_per_corpus.csv"
HISTORICAL_SUMMARY_JSON = ROOT / "results" / "zipf_kstat_scaling" / "summary.json"
HISTORICAL_REPORT_MD = ROOT / "results" / "zipf_kstat_scaling" / "report.md"

FORCED_SCALING_MODEL = "log(k_stat) = alpha * log(V)"
HALF_REFERENCE = 0.5
CI_METHOD = "Student-t 95% CI on the forced-fit alpha using curve_fit covariance SE and df = n - 1"
PVALUE_METHOD = "Two-sided one-sample t-test on per-corpus alpha_stat values against 0.5"
