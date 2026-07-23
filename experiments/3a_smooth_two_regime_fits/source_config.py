from pathlib import Path
import sys


EXPERIMENT_DIR = Path(__file__).resolve().parent
ROOT = EXPERIMENT_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from zipf_analysis_common import SEARCHED_CORPORA


OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
N_RANDOM_STARTS = 100
MAX_NFEV = 12000
SMOOTH_MODEL_PARAMETER_COUNT = 8
PIECEWISE_BREAKPOINT_RANK = 500
SIGMOID_GATE_DEFINITION = "sigma_r = 1 / (1 + exp[(log r - log k) / w])"
RERANK_SHIFT_FRACTION = 0.8

BOUND_A = (5.0, 50.0)
BOUND_B = (0.5, 3.0)
BOUND_C = (0.0, 1000.0)
BOUND_K = (20.0, 2000.0)
BOUND_W = (0.1, 3.0)

RELAXED_BOUND_A = (-100.0, 100.0)
RELAXED_BOUND_B = (0.5, 5.0)
RELAXED_BOUND_C = (0.0, 5000.0)
RELAXED_BOUND_K = (20.0, 2000.0)
RELAXED_BOUND_W = (0.1, 5.0)
RELAXED_SEED_BASE = 20260800

HISTORICAL_CONSTRAINED_SUMMARY = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
HISTORICAL_RELAXED_SUMMARY = ROOT / "results" / "zipf_reranked_model_all_corpora_relaxed" / "summary.json"
