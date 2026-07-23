from pathlib import Path
import sys


EXPERIMENT_DIR = Path(__file__).resolve().parent
ROOT = EXPERIMENT_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from zipf_analysis_common import DATA_DIR, RESULTS_DIR, SEARCHED_CORPORA, STEP2_OURS_EXPR, STEPL2_LOWC_EXPR


OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
BEAM_WIDTH = 50
EMITTED_STEP = 2
REPLAY_MAX_STEPS = 2
HISTORICAL_MAX_STEPS = 10
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
X_LOW = 0.05
X_HIGH = 1.0
CONSTANT_VARIANCE_THRESHOLD = 1e-10
SAMPLE_POINTS = 0
EXP_CLAMP = 30.0
VALUE_ABS_LIMIT = 1.0e6

DUMMY_CORPUS_URL = "unused://local-corpus"

XPow_EXPR = "sub[pow[x,x],sqrt[x]]"
EUCLIDEAN_EXPR = "mul[sub[1,x],sub[1,x]]"

DISPLAY_OVERRIDES = {
    STEP2_OURS_EXPR: "(x - 1) - log(x)",
    STEPL2_LOWC_EXPR: "exp(x - 1) - x",
    XPow_EXPR: "x^x - sqrt(x)",
    EUCLIDEAN_EXPR: "(1 - x)^2",
}

FAMILY_LABELS = {
    STEP2_OURS_EXPR: "is",
    STEPL2_LOWC_EXPR: "exp",
    XPow_EXPR: "xpow",
    EUCLIDEAN_EXPR: "euclidean",
}


def corpus_path(entry: dict) -> Path:
    return DATA_DIR / entry["filename"]


def historical_summary_path(entry: dict) -> Path:
    return RESULTS_DIR / entry["enriched_dir"] / "summary.json"


def display_expression(expr: str, math_expr: str) -> str:
    return DISPLAY_OVERRIDES.get(expr, math_expr)


def winner_family(expr: str) -> str:
    return FAMILY_LABELS.get(expr, "other")
