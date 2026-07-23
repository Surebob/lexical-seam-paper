from __future__ import annotations

import math

import numpy as np
from scipy import special


def logistic_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = np.clip((log_ranks - math.log(k)) / w_gate, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(z))


def tanh_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = (log_ranks - math.log(k)) / w_gate
    return 0.5 * (1.0 - np.tanh(z))


def erf_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = (log_ranks - math.log(k)) / w_gate
    return 0.5 * (1.0 - special.erf(z))


def algebraic_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = log_ranks - math.log(k)
    return 0.5 * (1.0 - z / np.sqrt((w_gate * w_gate) + (z * z)))


def arctan_gate(log_ranks: np.ndarray, k: float, w_gate: float) -> np.ndarray:
    z = log_ranks - math.log(k)
    return 0.5 - np.arctan(z / w_gate) / math.pi


GATE_FUNCS = {
    "logistic": logistic_gate,
    "tanh": tanh_gate,
    "erf": erf_gate,
    "algebraic": algebraic_gate,
    "arctan": arctan_gate,
}


def gate_width_scale(gate_name: str) -> float:
    return 2.0 if gate_name == "tanh" else 1.0

