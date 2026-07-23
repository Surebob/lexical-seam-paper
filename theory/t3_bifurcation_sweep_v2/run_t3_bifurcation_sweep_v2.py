from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:  # pragma: no cover
    plt = None


ROOT = Path("/Volumes/External2TB/emlexperiment")
OUTDIR = ROOT / "phase2_addon" / "t3_bifurcation_sweep_v2"
SWEEP_DIR = OUTDIR / "break_rank_sweep"
BEAMS_DIR = SWEEP_DIR / "per_config_beams"
SEARCH_PATH = ROOT / "eml_zipf_enriched_search.py"

V = 10_000
B_HEAD = 0.8
B_TAIL = 1.5
C_CORE = 10.0
K_VALUES = [50, 100, 200, 500, 1000, 2000, 5000]
BEAM_WIDTH = 50
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
CONSTANT_VARIANCE_THRESHOLD = 1e-10
X_LOW = 0.05
X_HIGH = 1.0


def load_search_module():
    spec = importlib.util.spec_from_file_location("t3v2_search", SEARCH_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


search = load_search_module()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def normalized_log_rank_x(n: int) -> np.ndarray:
    ranks = np.arange(1, n + 1, dtype=np.float64)
    return X_LOW + (X_HIGH - X_LOW) * np.log(ranks) / np.log(float(n))


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0 or not math.isfinite(denom):
        return float("nan")
    return float(np.dot(a, b) / denom)


def _solve_affine(z: np.ndarray, y: np.ndarray):
    design = np.column_stack([np.ones_like(z), z])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ coeffs
    mse = float(np.mean((pred - y) ** 2))
    return float(coeffs[0]), float(coeffs[1]), pred, mse


def fit_zipf_mandelbrot(ranks: np.ndarray, log_freq: np.ndarray):
    max_rank = float(np.max(ranks))
    c_grid = np.concatenate([np.array([0.0]), np.geomspace(1e-6, max_rank, 4096)])
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
            "hit_c_upper_bound": bool(c >= 0.999999 * max_rank),
        }
        if best is None or candidate["mse"] < best["mse"]:
            best = candidate
    return best


def structured_two_regime(K: int):
    ranks = np.arange(1, V + 1, dtype=np.float64)
    A = 1.0
    B = A * (float(K) + C_CORE) ** (-B_HEAD) * float(K) ** B_TAIL
    head = A * (ranks + C_CORE) ** (-B_HEAD)
    tail = B * ranks ** (-B_TAIL)
    freqs = np.where(ranks <= float(K), head, tail)
    continuity_left = A * (float(K) + C_CORE) ** (-B_HEAD)
    continuity_right = B * float(K) ** (-B_TAIL)
    return {
        "ranks": ranks,
        "freqs": freqs,
        "A": A,
        "B": B,
        "continuity_left": continuity_left,
        "continuity_right": continuity_right,
        "continuity_abs_error": abs(continuity_left - continuity_right),
    }


def eval_node(node, x: np.ndarray):
    if isinstance(node, str):
        if node == "x":
            return x
        if node == "1":
            return np.ones_like(x)
        raise ValueError(f"Unsupported terminal {node!r}")
    op = node[0]
    if op in search.UNARY_FUNCS:
        child = eval_node(node[1], x)
        if child is None:
            return None
        return search.UNARY_FUNCS[op](child)
    if op in search.BINARY_FUNCS:
        left = eval_node(node[1], x)
        right = eval_node(node[2], x)
        if left is None or right is None:
            return None
        return search.BINARY_FUNCS[op](left, right)
    raise ValueError(f"Unsupported op {op!r}")


def bregman_boundary_diagnostics(candidate: dict):
    x_dense = np.linspace(X_LOW, X_HIGH, 5000, dtype=np.float64)
    values = eval_node(candidate["node"], x_dense)
    if values is None or not np.all(np.isfinite(values)):
        return {
            "f_at_1": float("nan"),
            "fprime_at_1": float("nan"),
            "fsecond_min": float("nan"),
            "bregman_conditions_pass": False,
        }
    first = np.gradient(values, x_dense, edge_order=2)
    second = np.gradient(first, x_dense, edge_order=2)
    f_at_1 = float(values[-1])
    fprime_at_1 = float(first[-1])
    fsecond_min = float(np.min(second))
    return {
        "f_at_1": f_at_1,
        "fprime_at_1": fprime_at_1,
        "fsecond_min": fsecond_min,
        "bregman_conditions_pass": bool(abs(f_at_1) <= 1e-5 and abs(fprime_at_1) <= 1e-4 and fsecond_min > -1e-5),
    }


def reference_generators(x: np.ndarray):
    return {
        "IS": (x - 1.0) - np.log(x),
        "exp": np.exp(np.clip(x - 1.0, -700.0, 700.0)) - x,
        "Euclidean": (1.0 - x) ** 2,
        "xpow_sqrt": np.power(x, x) - np.sqrt(x),
    }


def run_search_depth2(x: np.ndarray, residual: np.ndarray):
    current = search.initial_vocabulary(x, residual)
    steps = []
    for step in range(1, 3):
        generated = search.generate_candidates(current, residual, step)
        generated = search.dedupe_candidates(generated)
        generated = search.filter_candidates(generated, CONSTANT_VARIANCE_THRESHOLD)
        generated = sorted(generated, key=lambda item: (item["rmse"], item["expr"]))
        if not generated:
            raise RuntimeError(f"No SR candidates generated at step {step}.")
        diverse = search.select_diverse_beam(generated, BEAM_WIDTH, DIVERSITY_WEIGHT)
        steps.append({"step": step, "generated": generated, "diverse": diverse})
        if step < KEEP_ALL_UNTIL_STEP:
            current = current + generated
        else:
            current = diverse
    return steps


def candidate_summary(candidate: dict, x: np.ndarray):
    refs = reference_generators(x)
    values = np.asarray(candidate["values"], dtype=np.float64)
    cosines = {name: cosine(values, ref) for name, ref in refs.items()}
    closest = max(cosines, key=lambda name: abs(cosines[name]))
    diag = bregman_boundary_diagnostics(candidate)
    return cosines, closest, diag


def run_config(K: int, write_beam: bool = True):
    x = normalized_log_rank_x(V)
    data = structured_two_regime(K)
    log_freq = np.log(data["freqs"])
    fit = fit_zipf_mandelbrot(data["ranks"], log_freq)
    if not np.isfinite(fit["c"]) or fit["hit_c_upper_bound"]:
        raise RuntimeError(f"ZM fit failed or hit c upper bound at K={K}: c={fit['c']}")
    residual = log_freq - fit["prediction"]
    steps = run_search_depth2(x, residual)
    step2 = steps[-1]["generated"]
    winner = step2[0]
    cosines, closest, diag = candidate_summary(winner, x)
    top20 = []
    for rank, candidate in enumerate(step2[:20], start=1):
        candidate_cosines, candidate_closest, candidate_diag = candidate_summary(candidate, x)
        top20.append(
            {
                "rank": rank,
                "expression": candidate["expr"],
                "math": candidate["math"],
                "rmse": float(candidate["rmse"]),
                "cosine_vs_IS": candidate_cosines["IS"],
                "cosine_vs_exp": candidate_cosines["exp"],
                "cosine_vs_Euclidean": candidate_cosines["Euclidean"],
                "cosine_vs_xx_sqrtx": candidate_cosines["xpow_sqrt"],
                "closest_generator": candidate_closest,
                "bregman_conditions_pass": candidate_diag["bregman_conditions_pass"],
                "f_at_1": candidate_diag["f_at_1"],
                "fprime_at_1": candidate_diag["fprime_at_1"],
                "fsecond_min": candidate_diag["fsecond_min"],
            }
        )
    if write_beam:
        write_json(BEAMS_DIR / f"K_{K}_top20_beam.json", {"K": K, "top20_step2": top20})
    top5_packed = [
        {"rank": rank, "expression": cand["expr"], "math": cand["math"], "rmse": float(cand["rmse"])}
        for rank, cand in enumerate(step2[:5], start=1)
    ]
    return {
        "K": K,
        "A": data["A"],
        "B": data["B"],
        "continuity_abs_error": data["continuity_abs_error"],
        "c": float(fit["c"]),
        "b": float(fit["b"]),
        "zm_rmse": float(fit["rmse"]),
        "residual_rmse": float(np.sqrt(np.mean(residual * residual))),
        "winner_expr": winner["expr"],
        "winner_math": winner["math"],
        "winner_rmse": float(winner["rmse"]),
        "top5": top5_packed,
        "cosine_vs_IS": cosines["IS"],
        "cosine_vs_exp": cosines["exp"],
        "cosine_vs_Euclidean": cosines["Euclidean"],
        "cosine_vs_xx_sqrtx": cosines["xpow_sqrt"],
        "closest_generator": closest,
        "bregman_conditions_pass": diag["bregman_conditions_pass"],
        "top20_rmse_spread": float(step2[19]["rmse"] - step2[0]["rmse"]) if len(step2) >= 20 else float("nan"),
    }


def preflight():
    BEAMS_DIR.mkdir(parents=True, exist_ok=True)
    endpoint_results = [run_config(50, write_beam=True), run_config(5000, write_beam=True)]
    c_delta = abs(endpoint_results[1]["c"] - endpoint_results[0]["c"])
    payload = {
        "status": "pass" if c_delta > 50.0 else "fail",
        "pass_condition": "abs(c_K5000 - c_K50) > 50",
        "c_delta": c_delta,
        "winner_differs": endpoint_results[0]["winner_expr"] != endpoint_results[1]["winner_expr"],
        "endpoints": endpoint_results,
    }
    write_json(SWEEP_DIR / "preflight_endpoints.json", payload)
    return payload


def analyze(rows: list[dict]):
    winner_changes = []
    closest_changes = []
    for prev, cur in zip(rows[:-1], rows[1:]):
        if prev["winner_expr"] != cur["winner_expr"]:
            winner_changes.append(
                {
                    "from_K": prev["K"],
                    "to_K": cur["K"],
                    "from_c": prev["c"],
                    "to_c": cur["c"],
                    "from_winner": prev["winner_expr"],
                    "to_winner": cur["winner_expr"],
                    "from_closest_generator": prev["closest_generator"],
                    "to_closest_generator": cur["closest_generator"],
                }
            )
        if prev["closest_generator"] != cur["closest_generator"]:
            closest_changes.append(
                {
                    "from_K": prev["K"],
                    "to_K": cur["K"],
                    "from_c": prev["c"],
                    "to_c": cur["c"],
                    "from_closest_generator": prev["closest_generator"],
                    "to_closest_generator": cur["closest_generator"],
                }
            )
    c_band_hits = [
        {"K": row["K"], "c": row["c"], "winner_expr": row["winner_expr"], "closest_generator": row["closest_generator"]}
        for row in rows
        if 66.0 <= row["c"] <= 79.0
    ]
    if not winner_changes:
        transition_type = "no_winner_identity_change"
    elif any({change["from_closest_generator"], change["to_closest_generator"]} == {"IS", "exp"} for change in winner_changes + closest_changes):
        transition_type = "sharp_IS_exp_bifurcation_candidate"
    else:
        transition_type = "winner_change_not_aligned_to_IS_exp"
    return {
        "winner_identity_changes": winner_changes,
        "closest_generator_changes": closest_changes,
        "c_band_66_79_hits": c_band_hits,
        "transition_type": transition_type,
        "winner_unchanged_entire_sweep": len(winner_changes) == 0,
        "c_min": min(row["c"] for row in rows),
        "c_max": max(row["c"] for row in rows),
    }


def plot_rows(rows: list[dict]):
    if plt is None:
        return []
    paths = []
    K = np.asarray([row["K"] for row in rows], dtype=float)
    c = np.asarray([row["c"] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(K, c, marker="o")
    ax.axhspan(66, 79, color="orange", alpha=0.2, label="empirical 66-79 c band")
    ax.set_xscale("log")
    ax.set_xlabel("break rank K")
    ax.set_ylabel("fitted ZM c")
    ax.set_title("T3v2 fitted c vs break rank")
    ax.legend()
    fig.tight_layout()
    path = SWEEP_DIR / "fitted_c_vs_K.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(K, [row["cosine_vs_IS"] for row in rows], marker="o", label="IS")
    ax.plot(K, [row["cosine_vs_exp"] for row in rows], marker="o", label="exp")
    ax.set_xscale("log")
    ax.set_xlabel("break rank K")
    ax.set_ylabel("winner cosine")
    ax.set_title("Winner cosine vs Bregman references")
    ax.legend()
    fig.tight_layout()
    path = SWEEP_DIR / "winner_cosines_vs_K.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))
    return paths


def run_primary():
    rows = [run_config(K, write_beam=True) for K in K_VALUES]
    csv_rows = []
    for row in rows:
        csv_rows.append(
            {
                "K": row["K"],
                "c": f"{row['c']:.17g}",
                "b": f"{row['b']:.17g}",
                "zm_rmse": f"{row['zm_rmse']:.17g}",
                "residual_rmse": f"{row['residual_rmse']:.17g}",
                "winner_expr": row["winner_expr"],
                "winner_math": row["winner_math"],
                "winner_rmse": f"{row['winner_rmse']:.17g}",
                "top5_json": json.dumps(row["top5"], separators=(",", ":")),
                "cosine_vs_IS": f"{row['cosine_vs_IS']:.17g}",
                "cosine_vs_exp": f"{row['cosine_vs_exp']:.17g}",
                "cosine_vs_Euclidean": f"{row['cosine_vs_Euclidean']:.17g}",
                "cosine_vs_xx_sqrtx": f"{row['cosine_vs_xx_sqrtx']:.17g}",
                "closest_generator": row["closest_generator"],
                "bregman_conditions_pass": row["bregman_conditions_pass"],
                "top20_rmse_spread": f"{row['top20_rmse_spread']:.17g}",
            }
        )
    write_csv(
        SWEEP_DIR / "per_config_results.csv",
        csv_rows,
        [
            "K",
            "c",
            "b",
            "zm_rmse",
            "residual_rmse",
            "winner_expr",
            "winner_math",
            "winner_rmse",
            "top5_json",
            "cosine_vs_IS",
            "cosine_vs_exp",
            "cosine_vs_Euclidean",
            "cosine_vs_xx_sqrtx",
            "closest_generator",
            "bregman_conditions_pass",
            "top20_rmse_spread",
        ],
    )
    analysis = analyze(rows)
    analysis["plots"] = plot_rows(rows)
    write_json(SWEEP_DIR / "bifurcation_analysis.json", analysis)
    return rows, analysis


def write_methodology_note():
    lines = [
        "# T3v2 Methodology Note",
        "",
        "T3v2 replaces the fixed-exponent two-Pareto mixture from T3 with a c-generating two-regime rank-frequency construction.",
        "",
        "Synthetic data:",
        "",
        "`f(r) = A * (r + c_core)^(-b_head)` for `r <= K`",
        "",
        "`f(r) = B * r^(-b_tail)` for `r > K`",
        "",
        "`B` is chosen so the two pieces are continuous at `K`.",
        "",
        f"- V: `{V}`",
        f"- b_head: `{B_HEAD}`",
        f"- b_tail: `{B_TAIL}`",
        f"- c_core: `{C_CORE}`",
        f"- K values: `{K_VALUES}`",
        "- single-ZM fit: grid search over `c` plus affine least-squares solve for log-amplitude and slope at each `c`",
        "- residual: `log(f_observed) - log(f_ZM_predicted)`",
        "- SR coordinate: normalized log-rank `x = 0.05 + 0.95 log(r)/log(V)`",
        "- SR grammar: live Section 2.4 grammar from `eml_zipf_enriched_search.py`: terminals `x, 1`; unary `neg, inv, sqr, sqrt, exp, log`; binary `eml, add, sub, mul, div, pow`",
        f"- depth: `2`; beam width: `{BEAM_WIDTH}`; diversity weight: `{DIVERSITY_WEIGHT}`",
    ]
    (OUTDIR / "t3v2_methodology_note.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(preflight_payload: dict, rows: list[dict] | None, analysis: dict | None):
    lines = [
        "# T3v2 Aggregate Report",
        "",
        "T3v2 tests bifurcation structure using a break-rank sweep in an explicitly two-regime rank-frequency construction.",
        "",
        "## Preflight",
        "",
        f"- status: `{preflight_payload['status']}`",
        f"- c delta K=50 to K=5000: `{preflight_payload['c_delta']:.6g}`",
        f"- endpoint winners differ: `{preflight_payload['winner_differs']}`",
    ]
    for item in preflight_payload["endpoints"]:
        lines.append(
            f"- K `{item['K']}`: c `{item['c']:.6g}`, b `{item['b']:.6g}`, winner `{item['winner_expr']}`, closest `{item['closest_generator']}`"
        )
    if preflight_payload["status"] != "pass":
        lines.extend(["", "Preflight failed; full sweep was not run."])
        (OUTDIR / "t3v2_aggregate_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        return
    assert rows is not None and analysis is not None
    lines.extend(
        [
            "",
            "## Primary Sweep",
            "",
            f"- transition type: `{analysis['transition_type']}`",
            f"- winner identity changes: `{len(analysis['winner_identity_changes'])}`",
            f"- closest-generator changes: `{len(analysis['closest_generator_changes'])}`",
            f"- fitted c range: `{analysis['c_min']:.6g}` to `{analysis['c_max']:.6g}`",
            f"- c-band 66-79 hits: `{len(analysis['c_band_66_79_hits'])}`",
            "",
            "## Winner Sequence",
            "",
        ]
    )
    for row in rows:
        lines.append(
            f"- K `{row['K']}`: c `{row['c']:.6g}`, winner `{row['winner_expr']}`, closest `{row['closest_generator']}`, cos(IS) `{row['cosine_vs_IS']:.4f}`, cos(exp) `{row['cosine_vs_exp']:.4f}`"
        )
    lines.extend(["", "## Interpretation", ""])
    if analysis["transition_type"] == "sharp_IS_exp_bifurcation_candidate":
        lines.append("Primary shows a discrete IS/exp-related transition. Dense sweep around the transition is warranted next.")
    elif analysis["transition_type"] == "no_winner_identity_change":
        lines.append("Primary produced c variation but no winner change. Per the prompt, the next conditional test would be a direct ZM c sweep, but this run stops after primary.")
    else:
        lines.append("Primary shows winner changes, but not aligned with the hypothesized empirical c-band IS/exp bifurcation.")
    if analysis["c_band_66_79_hits"]:
        lines.extend(["", "## c-band hits", ""])
        for item in analysis["c_band_66_79_hits"]:
            lines.append(f"- K `{item['K']}`: c `{item['c']:.6g}`, winner `{item['winner_expr']}`, closest `{item['closest_generator']}`")
    (OUTDIR / "t3v2_aggregate_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    BEAMS_DIR.mkdir(parents=True, exist_ok=True)
    write_methodology_note()
    preflight_payload = preflight()
    if preflight_payload["status"] != "pass":
        write_report(preflight_payload, None, None)
        print(json.dumps({"outdir": str(OUTDIR), "status": "preflight_failed", "c_delta": preflight_payload["c_delta"]}, indent=2))
        return
    rows, analysis = run_primary()
    write_report(preflight_payload, rows, analysis)
    print(json.dumps({"outdir": str(OUTDIR), "status": "primary_complete", "transition_type": analysis["transition_type"]}, indent=2))


if __name__ == "__main__":
    main()
