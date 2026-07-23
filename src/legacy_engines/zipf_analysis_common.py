import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path("/Volumes/External2TB/emlexperiment")
DATA_DIR = ROOT / "data" / "zipf"
RESULTS_DIR = ROOT / "results"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = [
    "*** START OF THE PROJECT GUTENBERG EBOOK",
    "*** START OF THIS PROJECT GUTENBERG EBOOK",
]
END_MARKERS = [
    "*** END OF THE PROJECT GUTENBERG EBOOK",
    "*** END OF THIS PROJECT GUTENBERG EBOOK",
]

STEP2_OURS_EXPR = "sub[sub[x,1],log[x]]"
STEPL2_LOWC_EXPR = "eml[sub[x,1],eml[x,1]]"

SEARCHED_CORPORA = [
    {
        "slug": "shakespeare",
        "name": "Complete Works of Shakespeare",
        "filename": "pg100.txt",
        "enriched_dir": "zipf_enriched_search_full",
    },
    {
        "slug": "war_and_peace",
        "name": "War and Peace",
        "filename": "pg2600.txt",
        "enriched_dir": "zipf_enriched_war_and_peace_full_seq",
    },
    {
        "slug": "moby_dick",
        "name": "Moby Dick",
        "filename": "pg2701.txt",
        "enriched_dir": "zipf_enriched_moby_dick_full_seq",
    },
    {
        "slug": "king_james_bible",
        "name": "King James Bible",
        "filename": "pg10.txt",
        "enriched_dir": "zipf_enriched_bible_full",
    },
    {
        "slug": "federalist_papers",
        "name": "Federalist Papers",
        "filename": "pg1404.txt",
        "enriched_dir": "zipf_enriched_federalist_full",
    },
    {
        "slug": "grimms_fairy_tales",
        "name": "Grimm's Fairy Tales",
        "filename": "pg2591.txt",
        "enriched_dir": "zipf_enriched_grimms_full",
    },
    {
        "slug": "don_quixote",
        "name": "Don Quixote",
        "filename": "pg996.txt",
        "enriched_dir": "zipf_enriched_quixote_full",
    },
    {
        "slug": "pride_and_prejudice",
        "name": "Pride and Prejudice",
        "filename": "pg1342.txt",
        "enriched_dir": "zipf_enriched_pride_full",
    },
    {
        "slug": "canterbury_tales",
        "name": "Canterbury Tales",
        "filename": "pg2383.txt",
        "enriched_dir": "zipf_enriched_canterbury_full",
    },
    {
        "slug": "arabian_nights_vol1",
        "name": "Arabian Nights (Vol 1)",
        "filename": "pg3435.txt",
        "enriched_dir": "zipf_enriched_arabian_full",
    },
    {
        "slug": "aesops_fables",
        "name": "Aesop's Fables",
        "filename": "pg11339.txt",
        "enriched_dir": "zipf_enriched_aesop_full",
    },
    {
        "slug": "complete_sherlock_holmes",
        "name": "Complete Sherlock Holmes",
        "filename": "pg1661.txt",
        "enriched_dir": "zipf_enriched_sherlock_full",
    },
    {
        "slug": "jane_eyre",
        "name": "Jane Eyre",
        "filename": "pg1260.txt",
        "enriched_dir": "zipf_enriched_janeeyre_full",
    },
    {
        "slug": "dubliners",
        "name": "Dubliners",
        "filename": "pg2814.txt",
        "enriched_dir": "zipf_enriched_dubliners_full",
    },
    {
        "slug": "the_iliad",
        "name": "The Iliad",
        "filename": "pg6130.txt",
        "enriched_dir": "zipf_enriched_iliad_full",
    },
    {
        "slug": "democracy_in_america",
        "name": "Democracy in America",
        "filename": "pg815.txt",
        "enriched_dir": "zipf_enriched_democracy_full",
    },
    {
        "slug": "origin_of_species",
        "name": "Origin of Species",
        "filename": "pg1228.txt",
        "enriched_dir": "zipf_enriched_origin_full",
    },
    {
        "slug": "wealth_of_nations",
        "name": "Wealth of Nations",
        "filename": "pg3300.txt",
        "enriched_dir": "zipf_enriched_wealth_full",
    },
    {
        "slug": "les_miserables",
        "name": "Les Miserables",
        "filename": "pg135.txt",
        "enriched_dir": "zipf_enriched_lesmis_full",
    },
    {
        "slug": "decline_and_fall_vol1",
        "name": "Decline and Fall Vol 1",
        "filename": "pg731.txt",
        "enriched_dir": "zipf_enriched_decline_full",
    },
    {
        "slug": "emile",
        "name": "Emile",
        "filename": "pg5427.txt",
        "enriched_dir": "zipf_enriched_emile_full",
    },
    {
        "slug": "ulysses",
        "name": "Ulysses",
        "filename": "pg4300.txt",
        "enriched_dir": "zipf_enriched_ulysses_full",
    },
    {
        "slug": "collected_poe",
        "name": "Collected Poe",
        "filename": "pg2147.txt",
        "enriched_dir": "zipf_enriched_poe_full",
    },
    {
        "slug": "principia_ethica",
        "name": "Principia Ethica",
        "filename": "pg53430.txt",
        "enriched_dir": "zipf_enriched_principia_ethica_full",
    },
    {
        "slug": "critique_of_pure_reason",
        "name": "Critique of Pure Reason",
        "filename": "pg4280.txt",
        "enriched_dir": "zipf_enriched_kant_full",
    },
]

HIGH_C_BETA_CORPORA = [
    "federalist_papers",
    "democracy_in_america",
    "origin_of_species",
    "the_iliad",
    "king_james_bible",
    "war_and_peace",
    "shakespeare",
]


def strip_gutenberg_boilerplate(text: str) -> str:
    start = 0
    end = len(text)
    for marker in START_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            line_end = text.find("\n", idx)
            start = line_end + 1 if line_end != -1 else idx
            break
    for marker in END_MARKERS:
        idx = text.find(marker)
        if idx != -1:
            end = idx
            break
    return text[start:end]


def tokenize_text(text: str):
    return TOKEN_RE.findall(text.lower())


def build_zipf_dataset(corpus_path: Path):
    raw_text = corpus_path.read_text(encoding="utf-8", errors="ignore")
    clean_text = strip_gutenberg_boilerplate(raw_text)
    tokens = tokenize_text(clean_text)
    counts = Counter(tokens)
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    freqs = np.array([freq for _, freq in ranked], dtype=np.float64)
    ranks = np.arange(1, len(freqs) + 1, dtype=np.float64)
    log_rank = np.log(ranks)
    log_freq = np.log(freqs)
    return {
        "tokens": tokens,
        "counts": counts,
        "ranked": ranked,
        "freqs": freqs,
        "ranks": ranks,
        "log_rank": log_rank,
        "log_freq": log_freq,
        "token_count": len(tokens),
        "unique_words": len(freqs),
    }


def normalize_x(values: np.ndarray, low: float = 0.05, high: float = 1.0):
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if vmax <= vmin:
        return np.full_like(values, 0.5 * (low + high))
    scaled = (values - vmin) / (vmax - vmin)
    return low + (high - low) * scaled


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    diff = np.asarray(y_true, dtype=np.float64) - np.asarray(y_pred, dtype=np.float64)
    return float(math.sqrt(float(np.mean(diff * diff))))


def weighted_rmse(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray) -> float:
    diff = np.asarray(y_true, dtype=np.float64) - np.asarray(y_pred, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    return float(math.sqrt(float(np.sum(weights * diff * diff) / np.sum(weights))))


def step2_formula(x: np.ndarray):
    return (x - 1.0) - np.log(x)


def lowc_formula(x: np.ndarray):
    return np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x


def _solve_affine(z: np.ndarray, y: np.ndarray):
    design = np.column_stack([np.ones_like(z), z])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ coeffs
    mse = float(np.mean((pred - y) ** 2))
    return float(coeffs[0]), float(coeffs[1]), pred, mse


def fit_zipf_mandelbrot(ranks: np.ndarray, log_freq: np.ndarray):
    max_rank = float(np.max(ranks))
    c_grid = np.concatenate(
        [
            np.array([0.0], dtype=np.float64),
            np.geomspace(1e-6, max_rank, 4096, dtype=np.float64),
        ]
    )
    best = None
    for c in c_grid:
        z = np.log(ranks + c)
        intercept, slope, pred, mse = _solve_affine(z, log_freq)
        candidate = {
            "a": intercept,
            "b": -slope,
            "c": float(c),
            "mse": mse,
            "rmse": float(math.sqrt(max(mse, 0.0))),
            "prediction": pred,
        }
        if best is None or candidate["mse"] < best["mse"]:
            best = candidate
    return best


def get_corpus_spec(slug: str):
    for spec in SEARCHED_CORPORA:
        if spec["slug"] == slug:
            return spec
    raise KeyError(f"Unknown corpus slug: {slug}")


def corpus_path(spec: dict) -> Path:
    return DATA_DIR / spec["filename"]


def summary_path(spec: dict) -> Path:
    return RESULTS_DIR / spec["enriched_dir"] / "summary.json"


def report_path(spec: dict) -> Path:
    return RESULTS_DIR / spec["enriched_dir"] / "report.md"


def load_enriched_summary(spec: dict):
    return json.loads(summary_path(spec).read_text(encoding="utf-8"))


def get_step2_candidate(summary: dict):
    return summary["zm_search"]["step_summary"][1]["top_candidates"][0]


def get_step2_expr(summary: dict) -> str:
    return get_step2_candidate(summary)["expr"]


def step2_is_ours(summary: dict) -> int:
    return int(get_step2_expr(summary) == STEP2_OURS_EXPR)


def zm_prediction(summary: dict, ranks: np.ndarray):
    zm = summary["zm_baseline"]
    return zm["a"] - zm["b"] * np.log(ranks + zm["c"])

