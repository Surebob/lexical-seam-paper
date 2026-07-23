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
OUTDIR = ROOT / "phase2_addon" / "t3_bifurcation_sweep"
SWEEP_DIR = OUTDIR / "mixing_weight_sweep"
BEAMS_DIR = SWEEP_DIR / "per_config_beams"
SEARCH_PATH = ROOT / "eml_zipf_enriched_search.py"

V = 10_000
ALPHA1 = 1.5
ALPHA2 = 1.3
BEAM_WIDTH = 50
KEEP_ALL_UNTIL_STEP = 2
DIVERSITY_WEIGHT = 0.35
CONSTANT_VARIANCE_THRESHOLD = 1e-10
X_LOW = 0.05
X_HIGH = 1.0

IS_EXPR = "sub[sub[x,1],log[x]]"
EXP_EXPR = "eml[sub[x,1],eml[x,1]]"
EUCLIDEAN_EXPR = "mul[sub[1,x],sub[1,x]]"
XPOW_EXPR = "sub[pow[x,x],sqrt[x]]"


def load_search_module():
    spec = importlib.util.spec_from_file_location("t3_search", SEARCH_PATH)
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


def pareto_mixture(w1: float):
    ranks = np.arange(1, V + 1, dtype=np.float64)
    comp1 = ranks ** (-ALPHA1)
    comp2 = ranks ** (-ALPHA2)
    comp1 /= float(np.sum(comp1))
    comp2 /= float(np.sum(comp2))
    mixture = w1 * comp1 + (1.0 - w1) * comp2
    return ranks, comp1, comp2, mixture


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


def endpoint_two_component_check(w1: float, comp1: np.ndarray, comp2: np.ndarray):
    weighted1 = w1 * comp1
    weighted2 = (1.0 - w1) * comp2
    frac1_head = float(weighted1[0] / (weighted1[0] + weighted2[0]))
    frac1_tail = float(weighted1[-1] / (weighted1[-1] + weighted2[-1]))
    tail_cross = bool((frac1_head - 0.5) * (frac1_tail - 0.5) < 0)
    return {
        "w1": w1,
        "component1_mass": w1,
        "component2_mass": 1.0 - w1,
        "component1_fraction_at_rank1": frac1_head,
        "component1_fraction_at_rankV": frac1_tail,
        "head_tail_dominance_crosses": tail_cross,
        "meaningfully_two_component": bool(min(w1, 1.0 - w1) >= 0.05 and tail_cross),
    }


def analyze_transitions(rows: list[dict]):
    changes = []
    for prev, cur in zip(rows[:-1], rows[1:]):
        if prev["winner_expr"] != cur["winner_expr"]:
            changes.append(
                {
                    "from_w1": float(prev["w1"]),
                    "to_w1": float(cur["w1"]),
                    "from_c": float(prev["c"]),
                    "to_c": float(cur["c"]),
                    "from_winner": prev["winner_expr"],
                    "to_winner": cur["winner_expr"],
                    "from_closest_generator": prev["closest_generator"],
                    "to_closest_generator": cur["closest_generator"],
                }
            )
    c_band_hits = [
        {"w1": float(row["w1"]), "c": float(row["c"]), "winner_expr": row["winner_expr"], "closest_generator": row["closest_generator"]}
        for row in rows
        if 66.0 <= float(row["c"]) <= 79.0
    ]
    closest_changes = []
    for prev, cur in zip(rows[:-1], rows[1:]):
        if prev["closest_generator"] != cur["closest_generator"]:
            closest_changes.append(
                {
                    "from_w1": float(prev["w1"]),
                    "to_w1": float(cur["w1"]),
                    "from_c": float(prev["c"]),
                    "to_c": float(cur["c"]),
                    "from_closest_generator": prev["closest_generator"],
                    "to_closest_generator": cur["closest_generator"],
                }
            )
    if not changes:
        transition_type = "no_winner_identity_change"
    elif any(
        {change["from_closest_generator"], change["to_closest_generator"]} == {"IS", "exp"}
        for change in changes + closest_changes
    ):
        transition_type = "sharp_IS_exp_bifurcation_candidate"
    else:
        transition_type = "sharp_noncanonical_expression_changes"
    return {
        "winner_identity_changes": changes,
        "closest_generator_changes": closest_changes,
        "c_band_66_79_hits": c_band_hits,
        "transition_type": transition_type,
        "winner_unchanged_entire_sweep": bool(len(changes) == 0),
    }


def plot_summary(rows: list[dict]):
    if plt is None:
        return []
    paths = []
    w = np.asarray([float(row["w1"]) for row in rows])
    c = np.asarray([float(row["c"]) for row in rows])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(w, c, marker="o")
    ax.axhspan(66, 79, color="orange", alpha=0.2, label="empirical 66-79 band")
    ax.set_xlabel("w1")
    ax.set_ylabel("fitted ZM c")
    ax.set_title("T3 fitted c vs mixing weight")
    ax.legend()
    fig.tight_layout()
    path = SWEEP_DIR / "fitted_c_vs_w1.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(w, [float(row["cosine_vs_IS"]) for row in rows], marker="o", label="IS")
    ax.plot(w, [float(row["cosine_vs_exp"]) for row in rows], marker="o", label="exp")
    ax.set_xlabel("w1")
    ax.set_ylabel("winner cosine")
    ax.set_title("Winner cosine to IS and exp generators")
    ax.legend()
    fig.tight_layout()
    path = SWEEP_DIR / "winner_cosines_vs_w1.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))
    return paths


def run_primary_sweep():
    SWEEP_DIR.mkdir(parents=True, exist_ok=True)
    BEAMS_DIR.mkdir(parents=True, exist_ok=True)
    x = normalized_log_rank_x(V)
    output_rows = []
    endpoint_checks = []
    for i in range(1, 20):
        w1 = round(0.05 * i, 2)
        ranks, comp1, comp2, mixture = pareto_mixture(w1)
        if w1 in {0.05, 0.95}:
            endpoint_checks.append(endpoint_two_component_check(w1, comp1, comp2))
        log_freq = np.log(mixture)
        fit = fit_zipf_mandelbrot(ranks, log_freq)
        if not np.isfinite(fit["c"]) or fit["hit_c_upper_bound"]:
            raise RuntimeError(f"ZM fit failed or hit c upper bound at w1={w1}: c={fit['c']}")
        residual = log_freq - fit["prediction"]
        steps = run_search_depth2(x, residual)
        step2 = steps[-1]["generated"]
        if len(step2) < 5:
            raise RuntimeError(f"Too few step-2 candidates at w1={w1}: {len(step2)}")
        top20 = []
        for rank, candidate in enumerate(step2[:20], start=1):
            cosines, closest, diag = candidate_summary(candidate, x)
            top20.append(
                {
                    "rank": rank,
                    "expression": candidate["expr"],
                    "math": candidate["math"],
                    "rmse": float(candidate["rmse"]),
                    "cosine_vs_IS": cosines["IS"],
                    "cosine_vs_exp": cosines["exp"],
                    "cosine_vs_Euclidean": cosines["Euclidean"],
                    "cosine_vs_xx_sqrtx": cosines["xpow_sqrt"],
                    "closest_generator": closest,
                    "bregman_conditions_pass": diag["bregman_conditions_pass"],
                    "f_at_1": diag["f_at_1"],
                    "fprime_at_1": diag["fprime_at_1"],
                    "fsecond_min": diag["fsecond_min"],
                }
            )
        write_json(BEAMS_DIR / f"w1_{w1:.2f}_top20_beam.json", {"w1": w1, "top20_step2": top20})
        winner = step2[0]
        cosines, closest, diag = candidate_summary(winner, x)
        top5_packed = [
            {"rank": rank, "expression": cand["expr"], "math": cand["math"], "rmse": float(cand["rmse"])}
            for rank, cand in enumerate(step2[:5], start=1)
        ]
        top20_spread = float(step2[19]["rmse"] - step2[0]["rmse"]) if len(step2) >= 20 else float("nan")
        if top20_spread < 0.01 and not any(item["bregman_conditions_pass"] for item in top20):
            raise RuntimeError(f"Step-2 beam appears noise-like at w1={w1}: top20 RMSE spread {top20_spread}")
        output_rows.append(
            {
                "w1": f"{w1:.2f}",
                "c": f"{fit['c']:.17g}",
                "b": f"{fit['b']:.17g}",
                "zm_rmse": f"{fit['rmse']:.17g}",
                "winner_expr": winner["expr"],
                "winner_math": winner["math"],
                "winner_rmse": f"{float(winner['rmse']):.17g}",
                "top5_json": json.dumps(top5_packed, separators=(",", ":")),
                "cosine_vs_IS": f"{cosines['IS']:.17g}",
                "cosine_vs_exp": f"{cosines['exp']:.17g}",
                "cosine_vs_Euclidean": f"{cosines['Euclidean']:.17g}",
                "cosine_vs_xx_sqrtx": f"{cosines['xpow_sqrt']:.17g}",
                "closest_generator": closest,
                "bregman_conditions_pass": diag["bregman_conditions_pass"],
                "top20_rmse_spread": f"{top20_spread:.17g}",
            }
        )
    write_csv(
        SWEEP_DIR / "per_config_results.csv",
        output_rows,
        [
            "w1",
            "c",
            "b",
            "zm_rmse",
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
    analysis = analyze_transitions(output_rows)
    analysis["endpoint_two_component_checks"] = endpoint_checks
    analysis["plots"] = plot_summary(output_rows)
    write_json(SWEEP_DIR / "bifurcation_analysis.json", analysis)
    return output_rows, analysis


def write_methodology_note():
    lines = [
        "# T3 Methodology Note",
        "",
        "Synthetic data are rank-frequency curves generated from equal-support mixtures of two normalized discrete Pareto components over ranks `1..10000`.",
        "",
        f"- exponents: alpha1 `{ALPHA1}`, alpha2 `{ALPHA2}`",
        "- sweep: `w1 = 0.05, 0.10, ..., 0.95`",
        "- single-ZM fit: nonlinear least squares on log-frequencies, implemented as grid search over `c` with affine least-squares solve for intercept and slope at each `c`",
        "- residual: `log(f_observed) - log(f_ZM_predicted)`",
        "- SR coordinate: normalized log-rank `x = 0.05 + 0.95 log(r)/log(V)`",
        "- grammar: live Section 2.4 grammar from `eml_zipf_enriched_search.py`: terminals `x, 1`; unary `neg, inv, sqr, sqrt, exp, log`; binary `eml, add, sub, mul, div, pow`",
        f"- depth: `2`; beam width: `{BEAM_WIDTH}`; diversity weight: `{DIVERSITY_WEIGHT}`; keep-all-until-step: `{KEEP_ALL_UNTIL_STEP}`",
    ]
    (OUTDIR / "t3_methodology_note.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(rows: list[dict], analysis: dict):
    winners = [(row["w1"], row["winner_expr"], row["closest_generator"], row["c"]) for row in rows]
    lines = [
        "# T3 Aggregate Report",
        "",
        "T3 swept mixing weight in a two-Pareto rank-frequency mixture and ran the canonical depth-2 symbolic residual search on each single-ZM residual.",
        "",
        "## Primary Outcome",
        "",
        f"- transition type: `{analysis['transition_type']}`",
        f"- winner identity changes: `{len(analysis['winner_identity_changes'])}`",
        f"- closest-generator changes: `{len(analysis['closest_generator_changes'])}`",
        f"- configs with fitted `c` in empirical 66-79 band: `{len(analysis['c_band_66_79_hits'])}`",
        "",
        "## Winner Sequence",
        "",
    ]
    for w1, expr, closest, c in winners:
        lines.append(f"- w1 `{w1}`: c `{float(c):.6g}`, winner `{expr}`, closest `{closest}`")
    lines.extend(["", "## Transition Details", ""])
    if analysis["winner_identity_changes"]:
        for item in analysis["winner_identity_changes"]:
            lines.append(
                f"- winner changes between w1 `{item['from_w1']}` and `{item['to_w1']}`: `{item['from_winner']}` -> `{item['to_winner']}`; c `{item['from_c']:.6g}` -> `{item['to_c']:.6g}`"
            )
    else:
        lines.append("- no winner identity changes across the primary sweep")
    if analysis["c_band_66_79_hits"]:
        lines.append("")
        lines.append("## Empirical c-band hits")
        lines.append("")
        for item in analysis["c_band_66_79_hits"]:
            lines.append(f"- w1 `{item['w1']}`: c `{item['c']:.6g}`, winner `{item['winner_expr']}`, closest `{item['closest_generator']}`")
    lines.extend(["", "## Interpretation", ""])
    if analysis["transition_type"] == "sharp_IS_exp_bifurcation_candidate":
        lines.append("The primary sweep shows a discrete IS/exp-related transition, consistent with subclaim 2's bifurcation picture.")
    elif analysis["transition_type"] == "no_winner_identity_change":
        lines.append("The primary sweep does not show a bifurcation: winner identity is unchanged across all 19 configurations. Under the stated stop rule, this suggests the exponent gap or synthetic construction is not producing the empirical transition.")
    else:
        lines.append("The primary sweep shows discrete expression changes, but not the hypothesized exp-to-IS bifurcation.")
    (OUTDIR / "t3_aggregate_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_methodology_note()
    rows, analysis = run_primary_sweep()
    write_report(rows, analysis)
    if analysis["winner_unchanged_entire_sweep"]:
        print(json.dumps({"outdir": str(OUTDIR), "status": "stopped_after_primary_no_bifurcation", "transition_type": analysis["transition_type"]}, indent=2))
    else:
        print(json.dumps({"outdir": str(OUTDIR), "status": "primary_complete", "transition_type": analysis["transition_type"]}, indent=2))


if __name__ == "__main__":
    main()
