from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy import special
from scipy.optimize import brentq, least_squares


OUTDIR = Path("/Volumes/External2TB/emlexperiment/phase2_addon/s2_alternative_gates_v2")
Z_GRID = np.arange(-10.0, 10.0001, 0.01, dtype=float)
W = 1.0
DERIVATIVE_COSINE_THRESHOLD = 0.99
REPARAM_RESIDUAL_THRESHOLD = 1e-4


def sigma_logistic(z: np.ndarray, w: float = W) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(z / w))


def dsigma_logistic(z: np.ndarray, w: float = W) -> np.ndarray:
    ez = np.exp(z / w)
    return -(ez / (w * (1.0 + ez) ** 2))


def sigma_tanh(z: np.ndarray, w: float = W) -> np.ndarray:
    return 0.5 * (1.0 - np.tanh(z / w))


def dsigma_tanh(z: np.ndarray, w: float = W) -> np.ndarray:
    return -0.5 * (1.0 / w) * (1.0 / np.cosh(z / w) ** 2)


def sigma_erf(z: np.ndarray, w: float = W) -> np.ndarray:
    return 0.5 * (1.0 - special.erf(z / w))


def dsigma_erf(z: np.ndarray, w: float = W) -> np.ndarray:
    return -(1.0 / (math.sqrt(math.pi) * w)) * np.exp(-((z / w) ** 2))


def sigma_algebraic(z: np.ndarray, w: float = W) -> np.ndarray:
    return 0.5 * (1.0 - z / np.sqrt((w * w) + (z * z)))


def dsigma_algebraic(z: np.ndarray, w: float = W) -> np.ndarray:
    return -0.5 * (w * w) / np.power((w * w) + (z * z), 1.5)


def sigma_arctan(z: np.ndarray, w: float = W) -> np.ndarray:
    return 0.5 - np.arctan(z / w) / math.pi


def dsigma_arctan(z: np.ndarray, w: float = W) -> np.ndarray:
    return -(1.0 / (math.pi * w)) * (1.0 / (1.0 + (z / w) ** 2))


GATES = {
    "logistic": (sigma_logistic, dsigma_logistic),
    "tanh": (sigma_tanh, dsigma_tanh),
    "erf": (sigma_erf, dsigma_erf),
    "algebraic": (sigma_algebraic, dsigma_algebraic),
    "arctan": (sigma_arctan, dsigma_arctan),
}


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(a, b) / denom)


def reparameterization_sse(func_a, func_b) -> tuple[float, float, float]:
    starts = [
        (math.log(0.5), 0.0),
        (0.0, 0.0),
        (math.log(2.0), 0.0),
        (math.log(0.5), -1.0),
        (math.log(0.5), 1.0),
        (0.0, -1.0),
        (0.0, 1.0),
        (math.log(2.0), -1.0),
        (math.log(2.0), 1.0),
    ]
    best = None
    for log_alpha0, beta0 in starts:
        result = least_squares(
            lambda p: func_a(np.exp(p[0]) * Z_GRID + p[1]) - func_b(Z_GRID),
            x0=np.array([log_alpha0, beta0], dtype=float),
            bounds=(np.array([-8.0, -20.0]), np.array([8.0, 20.0])),
            max_nfev=20000,
        )
        log_alpha, beta = [float(x) for x in result.x]
        alpha = float(math.exp(log_alpha))
        residual = float(np.sum((func_a(alpha * Z_GRID + beta) - func_b(Z_GRID)) ** 2))
        if best is None or residual < best[2]:
            best = (alpha, beta, residual)
    return best


def quarter_tail_position(func) -> float:
    return float(brentq(lambda zz: float(func(np.array([zz], dtype=float))[0] - 0.25), 0.0, 100.0))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    pair_summary: list[tuple[str, str, float, float]] = []

    for gate_name, (sigma_func, deriv_func) in GATES.items():
        central_slope_abs = abs(float(deriv_func(np.array([0.0], dtype=float))[0]))
        rows.append(
            {
                "row_type": "gate_summary",
                "gate": gate_name,
                "gate_a": "",
                "gate_b": "",
                "derivative_cosine": "",
                "passes_derivative_cosine": "",
                "reparameterization_alpha": "",
                "reparameterization_beta": "",
                "reparameterization_residual": "",
                "passes_reparameterization": "",
                "central_slope_abs": repr(central_slope_abs),
                "value_z1": repr(float(sigma_func(np.array([1.0], dtype=float))[0])),
                "value_z3": repr(float(sigma_func(np.array([3.0], dtype=float))[0])),
                "value_z10": repr(float(sigma_func(np.array([10.0], dtype=float))[0])),
                "z_star_sigma_eq_0p25": repr(quarter_tail_position(sigma_func)),
            }
        )

    gate_names = list(GATES.keys())
    for i, gate_a in enumerate(gate_names):
        for gate_b in gate_names[i + 1 :]:
            sigma_a, deriv_a = GATES[gate_a]
            sigma_b, deriv_b = GATES[gate_b]
            derivative_cosine = cosine_similarity(deriv_a(Z_GRID), deriv_b(Z_GRID))
            alpha, beta, residual = reparameterization_sse(sigma_a, sigma_b)
            pair_summary.append((gate_a, gate_b, derivative_cosine, residual))
            rows.append(
                {
                    "row_type": "pairwise_check",
                    "gate": "",
                    "gate_a": gate_a,
                    "gate_b": gate_b,
                    "derivative_cosine": repr(derivative_cosine),
                    "passes_derivative_cosine": str(derivative_cosine < DERIVATIVE_COSINE_THRESHOLD),
                    "reparameterization_alpha": repr(alpha),
                    "reparameterization_beta": repr(beta),
                    "reparameterization_residual": repr(residual),
                    "passes_reparameterization": str(residual > REPARAM_RESIDUAL_THRESHOLD),
                    "central_slope_abs": "",
                    "value_z1": "",
                    "value_z3": "",
                    "value_z10": "",
                    "z_star_sigma_eq_0p25": "",
                }
            )

    csv_path = OUTDIR / "s2_pre_fit_sanity_check_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "row_type",
                "gate",
                "gate_a",
                "gate_b",
                "derivative_cosine",
                "passes_derivative_cosine",
                "reparameterization_alpha",
                "reparameterization_beta",
                "reparameterization_residual",
                "passes_reparameterization",
                "central_slope_abs",
                "value_z1",
                "value_z3",
                "value_z10",
                "z_star_sigma_eq_0p25",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    failing_derivative = [(a, b, c) for a, b, c, _ in pair_summary if c >= DERIVATIVE_COSINE_THRESHOLD]
    failing_reparam = [(a, b, r) for a, b, _, r in pair_summary if r <= REPARAM_RESIDUAL_THRESHOLD]

    lines = [
        "# S2 v2 Pre-fit Sanity Check",
        "",
        "- Stage 1: derivative-cosine check on `z in [-10, 10]` with step `0.01` and `w = 1.0`.",
        f"- Stage 1 threshold: cosine `< {DERIVATIVE_COSINE_THRESHOLD}` for every pair.",
        "- Stage 2: affine reparameterization fit minimizing sum-squared error between `A(alpha*z + beta)` and `B(z)`.",
        f"- Stage 2 threshold: residual `> {REPARAM_RESIDUAL_THRESHOLD}` for every non-equivalent pair.",
        "",
        "## Pairwise results",
        "",
    ]
    for gate_a, gate_b, derivative_cosine, residual in pair_summary:
        lines.append(
            f"- `{gate_a}` vs `{gate_b}`: derivative cosine `{derivative_cosine:.12f}`, reparameterization residual `{residual:.12e}`"
        )
    lines.extend(
        [
            "",
            f"- derivative-threshold failures: `{len(failing_derivative)}`",
            f"- reparameterization-threshold failures: `{len(failing_reparam)}`",
            "",
        ]
    )
    if failing_reparam:
        lines.append(
            "- Interpretation: at least one pair is reparameterization-equivalent (or numerically indistinguishable under affine reparameterization), so the fit sweep should not proceed until the lineup is adjusted."
        )
    else:
        lines.append(
            "- Interpretation: no pair collapsed to machine-precision affine reparameterization, so the lineup passes the equivalence test."
        )
    (OUTDIR / "s2_summary_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
