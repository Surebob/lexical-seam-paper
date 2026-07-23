from __future__ import annotations

from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
CORPORA_DIR = PACKAGE_ROOT / "corpora"
REFERENCE_DIR = PACKAGE_ROOT / "reference"
OUTPUTS_DIR = PACKAGE_ROOT / "outputs"

EXPECTED_DEPENDENCIES = {
    "numpy": "2.4.2",
    "scipy": "1.17.1",
    "psutil": "7.1.3",
}

SEARCHED_CORPORA = [
    {"slug": "shakespeare", "name": "Complete Works of Shakespeare", "filename": "pg100.txt"},
    {"slug": "war_and_peace", "name": "War and Peace", "filename": "pg2600.txt"},
    {"slug": "moby_dick", "name": "Moby Dick", "filename": "pg2701.txt"},
    {"slug": "king_james_bible", "name": "King James Bible", "filename": "pg10.txt"},
    {"slug": "federalist_papers", "name": "Federalist Papers", "filename": "pg1404.txt"},
    {"slug": "grimms_fairy_tales", "name": "Grimm's Fairy Tales", "filename": "pg2591.txt"},
    {"slug": "don_quixote", "name": "Don Quixote", "filename": "pg996.txt"},
    {"slug": "pride_and_prejudice", "name": "Pride and Prejudice", "filename": "pg1342.txt"},
    {"slug": "canterbury_tales", "name": "Canterbury Tales", "filename": "pg2383.txt"},
    {"slug": "arabian_nights_vol1", "name": "Arabian Nights (Vol 1)", "filename": "pg3435.txt"},
    {"slug": "aesops_fables", "name": "Aesop's Fables", "filename": "pg11339.txt"},
    {"slug": "complete_sherlock_holmes", "name": "Complete Sherlock Holmes", "filename": "pg1661.txt"},
    {"slug": "jane_eyre", "name": "Jane Eyre", "filename": "pg1260.txt"},
    {"slug": "dubliners", "name": "Dubliners", "filename": "pg2814.txt"},
    {"slug": "the_iliad", "name": "The Iliad", "filename": "pg6130.txt"},
    {"slug": "democracy_in_america", "name": "Democracy in America", "filename": "pg815.txt"},
    {"slug": "origin_of_species", "name": "Origin of Species", "filename": "pg1228.txt"},
    {"slug": "wealth_of_nations", "name": "Wealth of Nations", "filename": "pg3300.txt"},
    {"slug": "les_miserables", "name": "Les Miserables", "filename": "pg135.txt"},
    {"slug": "decline_and_fall_vol1", "name": "Decline and Fall Vol 1", "filename": "pg731.txt"},
    {"slug": "emile", "name": "Emile", "filename": "pg5427.txt"},
    {"slug": "ulysses", "name": "Ulysses", "filename": "pg4300.txt"},
    {"slug": "collected_poe", "name": "Collected Poe", "filename": "pg2147.txt"},
    {"slug": "principia_ethica", "name": "Principia Ethica", "filename": "pg53430.txt"},
    {"slug": "critique_of_pure_reason", "name": "Critique of Pure Reason", "filename": "pg4280.txt"},
]

CORPUS_NAME_BY_SLUG = {spec["slug"]: spec["name"] for spec in SEARCHED_CORPORA}
SEARCHED_SLUGS = [spec["slug"] for spec in SEARCHED_CORPORA]
STEP_A_TARGET_SLUGS = ["shakespeare", "aesops_fables", "moby_dick"]

ORDERED_GATES = ["logistic", "tanh", "erf", "algebraic", "arctan"]
INDEPENDENT_GATES = ["logistic", "erf", "algebraic", "arctan"]

DEFAULT_WORKERS = 20
DEFAULT_N_STARTS = 100
DEFAULT_MAX_NFEV = 12000
DEFAULT_BASE_SEED = 20260415
DEFAULT_BLAS_THREADS = 2
SMOKE_BIC_TOLERANCE = 1e-6

PARAM_COUNT_DECOUPLED = 9
K_LOWER = 20.0
K_UPPER = 1000.0
W_GATE_LOWER = 0.05
W_GATE_UPPER = 10.0
W_TAIL_LOWER = 0.05
W_TAIL_UPPER = 10.0
BOUND_TOL = 1e-6

TANH_BIC_THRESHOLD = 1.0
TANH_GATE_RATIO_THRESHOLD = 0.05

