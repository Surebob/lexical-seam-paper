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
PRIMARY_SCRIPT = ROOT / "phase2_addon" / "t3_bifurcation_sweep_v2" / "run_t3_bifurcation_sweep_v2.py"
OUTDIR = ROOT / "phase2_addon" / "t3_bifurcation_sweep_v2" / "dense_sweeps"


def load_primary_module():
    spec = importlib.util.spec_from_file_location("t3v2_primary_for_dense", PRIMARY_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


t3 = load_primary_module()


IS_EXPR = "sub[sub[x,1],log[x]]"
EXP_EXPR = "eml[sub[x,1],eml[x,1]]"


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def run_config_dense(K: int, beam_dir: Path):
    x = t3.normalized_log_rank_x(t3.V)
    data = t3.structured_two_regime(K)
    log_freq = np.log(data["freqs"])
    fit = t3.fit_zipf_mandelbrot(data["ranks"], log_freq)
    residual = log_freq - fit["prediction"]
    steps = t3.run_search_depth2(x, residual)
    step2 = steps[-1]["generated"]
    winner = step2[0]
    cosines, closest, diag = t3.candidate_summary(winner, x)
    is_entry = None
    exp_entry = None
    top20 = []
    for rank, candidate in enumerate(step2[:20], start=1):
        candidate_cosines, candidate_closest, candidate_diag = t3.candidate_summary(candidate, x)
        payload = {
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
        top20.append(payload)
    for rank, candidate in enumerate(step2, start=1):
        if candidate["expr"] == IS_EXPR and is_entry is None:
            is_entry = {"rank": rank, "rmse": float(candidate["rmse"]), "math": candidate["math"]}
        if candidate["expr"] == EXP_EXPR and exp_entry is None:
            exp_entry = {"rank": rank, "rmse": float(candidate["rmse"]), "math": candidate["math"]}
        if is_entry is not None and exp_entry is not None:
            break
    beam_dir.mkdir(parents=True, exist_ok=True)
    write_json(beam_dir / f"K_{K}_top20_beam.json", {"K": K, "top20_step2": top20})
    top5 = [
        {"rank": rank, "expression": cand["expr"], "math": cand["math"], "rmse": float(cand["rmse"])}
        for rank, cand in enumerate(step2[:5], start=1)
    ]
    is_rmse = is_entry["rmse"] if is_entry is not None else float("nan")
    exp_rmse = exp_entry["rmse"] if exp_entry is not None else float("nan")
    return {
        "K": K,
        "c": float(fit["c"]),
        "b": float(fit["b"]),
        "zm_rmse": float(fit["rmse"]),
        "residual_rmse": float(np.sqrt(np.mean(residual * residual))),
        "winner_expr": winner["expr"],
        "winner_math": winner["math"],
        "winner_rmse": float(winner["rmse"]),
        "top5": top5,
        "cosine_vs_IS": cosines["IS"],
        "cosine_vs_exp": cosines["exp"],
        "cosine_vs_Euclidean": cosines["Euclidean"],
        "cosine_vs_xx_sqrtx": cosines["xpow_sqrt"],
        "closest_generator": closest,
        "bregman_conditions_pass": diag["bregman_conditions_pass"],
        "is_rank": is_entry["rank"] if is_entry else None,
        "is_rmse": is_rmse,
        "exp_rank": exp_entry["rank"] if exp_entry else None,
        "exp_rmse": exp_rmse,
        "is_minus_exp_rmse": is_rmse - exp_rmse if math.isfinite(is_rmse) and math.isfinite(exp_rmse) else float("nan"),
        "both_is_exp_in_top5": bool(is_entry is not None and exp_entry is not None and is_entry["rank"] <= 5 and exp_entry["rank"] <= 5),
    }


def write_results_csv(path: Path, rows: list[dict]) -> None:
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
                "is_rank": row["is_rank"] or "",
                "is_rmse": f"{row['is_rmse']:.17g}",
                "exp_rank": row["exp_rank"] or "",
                "exp_rmse": f"{row['exp_rmse']:.17g}",
                "is_minus_exp_rmse": f"{row['is_minus_exp_rmse']:.17g}",
                "both_is_exp_in_top5": row["both_is_exp_in_top5"],
            }
        )
    write_csv(
        path,
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
            "is_rank",
            "is_rmse",
            "exp_rank",
            "exp_rmse",
            "is_minus_exp_rmse",
            "both_is_exp_in_top5",
        ],
    )


def transition_windows(rows: list[dict]):
    changes = []
    for prev, cur in zip(rows[:-1], rows[1:]):
        if prev["winner_expr"] != cur["winner_expr"]:
            changes.append(
                {
                    "from_K": prev["K"],
                    "to_K": cur["K"],
                    "from_c": prev["c"],
                    "to_c": cur["c"],
                    "from_winner": prev["winner_expr"],
                    "to_winner": cur["winner_expr"],
                    "from_gap": prev["is_minus_exp_rmse"],
                    "to_gap": cur["is_minus_exp_rmse"],
                }
            )
    sign_changes = []
    for prev, cur in zip(rows[:-1], rows[1:]):
        if prev["is_minus_exp_rmse"] == 0:
            sign_changes.append({"K": prev["K"], "gap": 0.0})
        elif prev["is_minus_exp_rmse"] * cur["is_minus_exp_rmse"] < 0:
            sign_changes.append(
                {
                    "from_K": prev["K"],
                    "to_K": cur["K"],
                    "from_gap": prev["is_minus_exp_rmse"],
                    "to_gap": cur["is_minus_exp_rmse"],
                    "linear_interpolated_zero_K": prev["K"]
                    + (cur["K"] - prev["K"]) * (0.0 - prev["is_minus_exp_rmse"]) / (cur["is_minus_exp_rmse"] - prev["is_minus_exp_rmse"]),
                }
            )
    is_wins = [row["K"] for row in rows if row["winner_expr"] == IS_EXPR]
    exp_wins = [row["K"] for row in rows if row["winner_expr"] == EXP_EXPR]
    return {
        "winner_changes": changes,
        "gap_sign_changes": sign_changes,
        "is_win_K_values": is_wins,
        "exp_win_K_values": exp_wins,
        "is_win_count": len(is_wins),
        "exp_win_count": len(exp_wins),
        "K_min": min(row["K"] for row in rows),
        "K_max": max(row["K"] for row in rows),
        "c_min": min(row["c"] for row in rows),
        "c_max": max(row["c"] for row in rows),
    }


def plot_sweep(name: str, rows: list[dict], outdir: Path):
    if plt is None:
        return []
    paths = []
    K = np.asarray([row["K"] for row in rows], dtype=float)
    c = np.asarray([row["c"] for row in rows], dtype=float)
    gap = np.asarray([row["is_minus_exp_rmse"] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(K, c, marker="o")
    ax.set_xlabel("K")
    ax.set_ylabel("fitted ZM c")
    ax.set_title(f"{name}: fitted c vs K")
    fig.tight_layout()
    path = outdir / "fitted_c_vs_K.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axhline(0.0, color="black", linewidth=1)
    ax.plot(K, gap, marker="o")
    ax.set_xlabel("K")
    ax.set_ylabel("IS RMSE - exp RMSE")
    ax.set_title(f"{name}: signed IS/exp RMSE gap")
    fig.tight_layout()
    path = outdir / "is_minus_exp_rmse_vs_K.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(K, [row["cosine_vs_IS"] for row in rows], marker="o", label="winner vs IS")
    ax.plot(K, [row["cosine_vs_exp"] for row in rows], marker="o", label="winner vs exp")
    ax.set_xlabel("K")
    ax.set_ylabel("cosine")
    ax.set_title(f"{name}: winner cosine trajectory")
    ax.legend()
    fig.tight_layout()
    path = outdir / "winner_cosines_vs_K.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(str(path))
    return paths


def run_sweep(name: str, values: list[int], outdir: Path):
    beam_dir = outdir / "per_config_beams"
    outdir.mkdir(parents=True, exist_ok=True)
    beam_dir.mkdir(parents=True, exist_ok=True)
    rows = [run_config_dense(K, beam_dir) for K in values]
    write_results_csv(outdir / "per_config_results.csv", rows)
    analysis = transition_windows(rows)
    analysis["plots"] = plot_sweep(name, rows, outdir)
    write_json(outdir / "transition_analysis.json", analysis)
    return rows, analysis


def classify_transition(analysis: dict):
    if not analysis["winner_changes"]:
        return "no_winner_transition"
    spans = [abs(item["to_K"] - item["from_K"]) for item in analysis["winner_changes"] if "to_K" in item]
    max_span = max(spans) if spans else 0
    if max_span <= 100:
        return "sharp_discrete_transition_width_le_100"
    if max_span <= 500:
        return "moderately_broad_transition_width_le_500"
    return "broad_or_irregular_transition"


def write_report(sweep_a_rows, sweep_a_analysis, sweep_b_rows, sweep_b_analysis):
    class_a = classify_transition(sweep_a_analysis)
    class_b = classify_transition(sweep_b_analysis)
    lines = [
        "# T3v2 Dense Sweep Report",
        "",
        "Dense sweeps refine the two transition windows observed in the primary T3v2 break-rank sweep.",
        "",
        "## Sweep A: K=1000..2000 step 100",
        "",
        f"- classification: `{class_a}`",
        f"- IS wins: `{sweep_a_analysis['is_win_count']}` / `{len(sweep_a_rows)}`",
        f"- exp wins: `{sweep_a_analysis['exp_win_count']}` / `{len(sweep_a_rows)}`",
        f"- c range: `{sweep_a_analysis['c_min']:.6g}` to `{sweep_a_analysis['c_max']:.6g}`",
        f"- winner changes: `{sweep_a_analysis['winner_changes']}`",
        f"- gap sign changes: `{sweep_a_analysis['gap_sign_changes']}`",
        "",
        "## Sweep B: K=2000..5000 step 300",
        "",
        f"- classification: `{class_b}`",
        f"- IS wins: `{sweep_b_analysis['is_win_count']}` / `{len(sweep_b_rows)}`",
        f"- exp wins: `{sweep_b_analysis['exp_win_count']}` / `{len(sweep_b_rows)}`",
        f"- c range: `{sweep_b_analysis['c_min']:.6g}` to `{sweep_b_analysis['c_max']:.6g}`",
        f"- winner changes: `{sweep_b_analysis['winner_changes']}`",
        f"- gap sign changes: `{sweep_b_analysis['gap_sign_changes']}`",
        "",
        "## Interpretation",
        "",
    ]
    if sweep_a_analysis["is_win_count"] > 1 or sweep_b_analysis["is_win_count"] > 1:
        lines.append("IS is not a single isolated point; it occupies a contiguous window in the dense sweeps.")
    else:
        lines.append("IS appears only at one or two points, suggesting a narrow-window phenomenon.")
    if class_a.startswith("sharp") or class_b.startswith("sharp"):
        lines.append("At least one transition is sharp at the current grid resolution.")
    if class_a.startswith("moderately") or class_b.startswith("moderately") or class_a.startswith("broad") or class_b.startswith("broad"):
        lines.append("At least one transition spans multiple grid intervals, so the signed RMSE gap should be read as a continuous crossing rather than a literal discontinuity.")
    (OUTDIR / "dense_sweep_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sweep_a_values = list(range(1000, 2000 + 1, 100))
    sweep_b_values = list(range(2000, 5000 + 1, 300))
    sweep_a_rows, sweep_a_analysis = run_sweep("sweep_A_1000_to_2000", sweep_a_values, OUTDIR / "sweep_A_1000_to_2000")
    sweep_b_rows, sweep_b_analysis = run_sweep("sweep_B_2000_to_5000", sweep_b_values, OUTDIR / "sweep_B_2000_to_5000")
    aggregate = {
        "sweep_A_1000_to_2000": sweep_a_analysis,
        "sweep_B_2000_to_5000": sweep_b_analysis,
        "overall": {
            "sweep_A_classification": classify_transition(sweep_a_analysis),
            "sweep_B_classification": classify_transition(sweep_b_analysis),
            "is_win_total": sweep_a_analysis["is_win_count"] + sweep_b_analysis["is_win_count"],
            "exp_win_total": sweep_a_analysis["exp_win_count"] + sweep_b_analysis["exp_win_count"],
        },
    }
    write_json(OUTDIR / "dense_sweep_analysis.json", aggregate)
    write_report(sweep_a_rows, sweep_a_analysis, sweep_b_rows, sweep_b_analysis)
    print(json.dumps({"outdir": str(OUTDIR), "sweep_A": aggregate["overall"]["sweep_A_classification"], "sweep_B": aggregate["overall"]["sweep_B_classification"]}, indent=2))


if __name__ == "__main__":
    main()
