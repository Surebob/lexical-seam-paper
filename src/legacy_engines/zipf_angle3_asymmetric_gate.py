from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
BASE_PATH = ROOT / "zipf_seam_mandelbrot_pmf.py"
SOFTK_PATH = ROOT / "zipf_seam_mandelbrot_softk.py"
ANGLE1_PATH = ROOT / "zipf_angle1_head_windows.py"
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_angle3_asymmetric_gate"

CUTOFFS = [50, 100, 200, 500, 1000, None]
TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
TRAIN_STARTS = 8
FULL_STARTS = 5
STEP2_VARIANCE_THRESHOLD = 1e-10


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_angle3_common")
base = load_module(BASE_PATH, "zipf_angle3_base")
softk_mod = load_module(SOFTK_PATH, "zipf_angle3_softk")
angle1 = load_module(ANGLE1_PATH, "zipf_angle3_angle1")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_angle3_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def parse_args():
    p = argparse.ArgumentParser(description="Asymmetric gate comparison with explicit split-fit soft-k baseline")
    p.add_argument("--softk-summary", type=Path, default=SOFTK_SUMMARY_PATH)
    p.add_argument("--base-summary", type=Path, default=BASE_SUMMARY_PATH)
    p.add_argument("--outdir", type=Path, default=OUTDIR)
    return p.parse_args()


def cutoff_label(cutoff: int | None) -> str:
    return "full" if cutoff is None else f"top{cutoff}"


def load_map(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def sigma_curve_asym(ranks: np.ndarray, k: float, w_left: float, w_right: float) -> np.ndarray:
    z = np.log(ranks) - math.log(k)
    out = np.empty_like(z, dtype=np.float64)
    mask = z <= 0.0
    out[mask] = 1.0 / (1.0 + np.exp(np.clip(z[mask] / w_left, -60.0, 60.0)))
    out[~mask] = 1.0 / (1.0 + np.exp(np.clip(z[~mask] / w_right, -60.0, 60.0)))
    return out


def asym_log_scores(ranks: np.ndarray, params: np.ndarray) -> np.ndarray | None:
    b1, c1, b2, c2, k, w_left, w_right = [float(v) for v in params]
    if b1 <= 0.0 or b2 <= 0.0 or c1 < 0.0 or c2 < 0.0 or k <= 0.0 or w_left <= 0.0 or w_right <= 0.0:
        return None
    sigma = sigma_curve_asym(ranks, k, w_left, w_right)
    head = -b1 * np.log(ranks + c1)
    tail_rank = base.smooth_tail_local_rank(ranks, k, w_right)
    tail_at_k = float(base.smooth_tail_local_rank(np.array([k], dtype=np.float64), k, w_right)[0])
    head_at_k = -b1 * math.log(k + c1)
    d = head_at_k + b2 * math.log(tail_at_k + c2)
    tail = d - b2 * np.log(tail_rank + c2)
    scores = sigma * head + (1.0 - sigma) * tail
    return scores if np.all(np.isfinite(scores)) else None


def asym_logpmf(params: np.ndarray, ranks: np.ndarray) -> np.ndarray | None:
    return base.normalized_logpmf_from_scores(asym_log_scores(ranks, params))


def asym_starts(vocab_size: int, reranked_row: dict | None, base_row: dict, soft_row: dict, prev_params: list[float] | None, seed: int, n_starts: int) -> list[np.ndarray]:
    k_lo, k_hi = base.k_bounds(vocab_size)
    lower = np.array(
        [base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], base.SEAM_B_BOUNDS[0], base.SEAM_C_BOUNDS[0], k_lo, base.SEAM_W_BOUNDS[0], base.SEAM_W_BOUNDS[0]],
        dtype=np.float64,
    )
    upper = np.array(
        [base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], base.SEAM_B_BOUNDS[1], base.SEAM_C_BOUNDS[1], k_hi, base.SEAM_W_BOUNDS[1], base.SEAM_W_BOUNDS[1]],
        dtype=np.float64,
    )
    starts = []
    if prev_params is not None:
        starts.append(np.asarray(prev_params, dtype=np.float64))

    for source in [
        base_row["models"]["seam"]["params"],
        protocol.split_fit_params(soft_row),
        reranked_row["params"] if reranked_row is not None else None,
    ]:
        if source is None:
            continue
        w = float(source["w"])
        k = float(min(max(source["k"], k_lo), k_hi))
        core = [float(source["b1"]), float(source["c1"]), float(source["b2"]), float(source["c2"]), k]
        starts.append(np.array(core + [w, w], dtype=np.float64))
        starts.append(np.array(core + [max(base.SEAM_W_BOUNDS[0], 0.5 * w), min(base.SEAM_W_BOUNDS[1], 1.8 * w)], dtype=np.float64))
        starts.append(np.array(core + [min(base.SEAM_W_BOUNDS[1], 1.8 * w), max(base.SEAM_W_BOUNDS[0], 0.5 * w)], dtype=np.float64))

    rng = np.random.default_rng(seed)
    while len(starts) < n_starts:
        starts.append(
            np.array(
                [
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    rng.uniform(*base.SEAM_B_BOUNDS),
                    rng.uniform(*base.SEAM_C_BOUNDS),
                    float(np.exp(rng.uniform(math.log(k_lo), math.log(k_hi)))),
                    float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                    float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                ],
                dtype=np.float64,
            )
        )
    return [np.clip(start, lower, upper) for start in starts[:n_starts]]


def fit_asym_softk(
    train_counts: np.ndarray,
    ranks: np.ndarray,
    vocab_size: int,
    reranked_row: dict | None,
    base_row: dict,
    soft_row: dict,
    lam: float,
    seed: int,
    n_starts: int,
    prev_params: list[float] | None = None,
) -> dict:
    bounds = [
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.SEAM_B_BOUNDS,
        base.SEAM_C_BOUNDS,
        base.k_bounds(vocab_size),
        base.SEAM_W_BOUNDS,
        base.SEAM_W_BOUNDS,
    ]
    lower = np.array([b[0] for b in bounds], dtype=np.float64)
    upper = np.array([b[1] for b in bounds], dtype=np.float64)
    target = softk_mod.SMOOTH_ALPHA * math.log(vocab_size)
    n_tokens = max(float(np.sum(train_counts)), 1.0)
    best = None

    def objective(theta: np.ndarray) -> float:
        params = np.clip(theta, lower, upper)
        logpmf = asym_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        if not math.isfinite(ll):
            return math.inf
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens
        return float(avg_nll + lam * dev * dev)

    for start_index, start in enumerate(asym_starts(vocab_size, reranked_row, base_row, soft_row, prev_params, seed, n_starts), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = asym_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens if math.isfinite(ll) else math.inf
        penalized = avg_nll + lam * dev * dev if math.isfinite(avg_nll) else math.inf
        wratio = float(params[6] / params[5])
        record = {
            "params": {
                "b1": float(params[0]),
                "c1": float(params[1]),
                "b2": float(params[2]),
                "c2": float(params[3]),
                "k": float(params[4]),
                "w_left": float(params[5]),
                "w_right": float(params[6]),
                "transition_fraction": float(np.log(params[4]) / np.log(vocab_size)),
                "width_ratio": wratio,
            },
            "params_vec": [float(v) for v in params],
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "train_avg_nll": float(avg_nll),
            "penalized_objective": float(penalized),
            "dev_logk": float(dev),
            "success": bool(result.success),
            "message": result.message,
            "nfev": int(result.nfev),
            "start_index": int(start_index),
            "p": 7,
            "lambda": float(lam),
        }
        if best is None or record["penalized_objective"] < best["penalized_objective"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def full_context_asym(dataset: dict, reranked_row: dict | None, base_row: dict, soft_row: dict, lam: float, seed: int, prev_params: list[float] | None = None) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    fit = fit_asym_softk(counts, ranks, int(dataset["unique_words"]), reranked_row, base_row, soft_row, lam, seed, FULL_STARTS, prev_params=prev_params)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {
        "fit": fit,
        "rmse": base.full_rank_rmse(counts, fit["logpmf"]),
        "prediction_log": pred_log,
    }


def analyze_corpus(spec: dict, base_row: dict, soft_row: dict, reranked_row: dict | None, corpus_index: int) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    ranks = np.arange(1, len(train_counts) + 1, dtype=np.float64)
    vocab_size = split["vocab_size"]

    lambda_rows = []
    prev_params = None
    for lambda_index, lam in enumerate(softk_mod.LAMBDAS, start=1):
        fit = fit_asym_softk(
            train_counts,
            ranks,
            vocab_size,
            reranked_row,
            base_row,
            soft_row,
            lam,
            seed=SPLIT_SEED + 11000 + corpus_index * 50 + lambda_index,
            n_starts=TRAIN_STARTS,
            prev_params=prev_params,
        )
        prev_params = fit["params_vec"]
        evals = base.evaluate_test_fit(fit, test_counts)
        lambda_rows.append(
            {
                "lambda": float(lam),
                "fit": fit,
                "test_loglike": float(evals["test_loglike"]),
                "test_avg_nll": float(evals["test_avg_nll"]),
            }
        )

    best_lambda_row = min(lambda_rows, key=lambda row: row["test_avg_nll"])
    best_full = full_context_asym(
        dataset,
        reranked_row,
        base_row,
        soft_row,
        lam=float(best_lambda_row["lambda"]),
        seed=SPLIT_SEED + 14000 + corpus_index,
        prev_params=best_lambda_row["fit"]["params_vec"],
    )

    model_logpmfs = {
        "zipf": angle1.logpmf_from_params("zipf", base_row["models"]["zipf"]["params"], ranks),
        "zm": angle1.logpmf_from_params("zm", base_row["models"]["zm"]["params"], ranks),
        "moe": angle1.logpmf_from_params("moe", base_row["models"]["moe"]["params"], ranks),
        "softk": angle1.logpmf_from_params("softk", protocol.split_fit_params(soft_row), ranks),
        "asym": best_lambda_row["fit"]["logpmf"],
    }

    heldout = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        vals = {name: angle1.restricted_avg_nll(test_counts, logpmf, cutoff) for name, logpmf in model_logpmfs.items()}
        winner = min(vals, key=vals.get)
        heldout[label] = {
            "avg_nll": vals,
            "winner": winner,
            "asym_minus_moe": float(vals["asym"] - vals["moe"]),
            "asym_minus_zm": float(vals["asym"] - vals["zm"]),
            "asym_minus_softk": float(vals["asym"] - vals["softk"]),
            "asym_minus_zipf": float(vals["asym"] - vals["zipf"]),
        }

    step2 = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        end = len(dataset["ranks"]) if cutoff is None else min(int(cutoff), len(dataset["ranks"]))
        step2[label] = angle1.exact_step2_window(
            dataset["log_rank"][:end],
            dataset["log_freq"][:end],
            best_full["prediction_log"][:end],
        )

    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "best_lambda": {
            "lambda": float(best_lambda_row["lambda"]),
            "selection": {
                "test_avg_nll": float(best_lambda_row["test_avg_nll"]),
                "test_loglike": float(best_lambda_row["test_loglike"]),
                "train_avg_nll": float(best_lambda_row["fit"]["train_avg_nll"]),
                "params": best_lambda_row["fit"]["params"],
            },
            "full_refit": {
                "bic": float(best_full["fit"]["bic"]),
                "rmse": float(best_full["rmse"]),
                "params": best_full["fit"]["params"],
            },
            "protocol_note": "selection.* are split-fit quantities and are canonical for held-out comparisons; full_refit.* are full-corpus diagnostics only.",
        },
        "heldout": heldout,
        "step2": step2,
        "sweep": [
            {
                "lambda": row["lambda"],
                "test_avg_nll": row["test_avg_nll"],
                "train_avg_nll": row["fit"]["train_avg_nll"],
                "width_ratio": row["fit"]["params"]["width_ratio"],
            }
            for row in lambda_rows
        ],
    }


def summarize(rows: list[dict]) -> dict:
    cutoff_summary = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        winner_counts = {"zipf": 0, "zm": 0, "moe": 0, "softk": 0, "asym": 0}
        asym_beats_moe = 0
        asym_beats_zm = 0
        asym_beats_softk = 0
        asym_beats_zipf = 0
        deltas_moe = []
        deltas_zm = []
        deltas_softk = []
        step2_help = 0
        for row in rows:
            info = row["heldout"][label]
            winner_counts[info["winner"]] += 1
            asym_beats_moe += info["avg_nll"]["asym"] < info["avg_nll"]["moe"]
            asym_beats_zm += info["avg_nll"]["asym"] < info["avg_nll"]["zm"]
            asym_beats_softk += info["avg_nll"]["asym"] < info["avg_nll"]["softk"]
            asym_beats_zipf += info["avg_nll"]["asym"] < info["avg_nll"]["zipf"]
            deltas_moe.append(info["asym_minus_moe"])
            deltas_zm.append(info["asym_minus_zm"])
            deltas_softk.append(info["asym_minus_softk"])
            step2_help += int(row["step2"][label]["helps"])
        cutoff_summary[label] = {
            "winner_counts": winner_counts,
            "asym_beats_moe": int(asym_beats_moe),
            "asym_beats_zm": int(asym_beats_zm),
            "asym_beats_softk": int(asym_beats_softk),
            "asym_beats_zipf": int(asym_beats_zipf),
            "median_asym_minus_moe": float(np.median(deltas_moe)),
            "median_asym_minus_zm": float(np.median(deltas_zm)),
            "median_asym_minus_softk": float(np.median(deltas_softk)),
            "step2_help_count": int(step2_help),
        }
    width_ratios = [row["best_lambda"]["selection"]["params"]["width_ratio"] for row in rows]
    material_asym = sum(max(ratio, 1.0 / ratio) > 1.5 for ratio in width_ratios)
    return {
        "rows": rows,
        "cutoffs": cutoff_summary,
        "width_ratio": {
            "median": float(np.median(width_ratios)),
            "mean": float(np.mean(width_ratios)),
            "material_asymmetry_count": int(material_asym),
        },
    }


def write_csv(rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "lambda",
        "width_ratio",
        "w_left",
        "w_right",
        "full_asym",
        "full_softk",
        "full_moe",
        "full_zm",
        "full_step2_help",
        "full_step2_gain",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            params = row["best_lambda"]["selection"]["params"]
            writer.writerow(
                {
                    "slug": row["slug"],
                    "name": row["name"],
                    "lambda": row["best_lambda"]["lambda"],
                    "width_ratio": params["width_ratio"],
                    "w_left": params["w_left"],
                    "w_right": params["w_right"],
                    "full_asym": row["heldout"]["full"]["avg_nll"]["asym"],
                    "full_softk": row["heldout"]["full"]["avg_nll"]["softk"],
                    "full_moe": row["heldout"]["full"]["avg_nll"]["moe"],
                    "full_zm": row["heldout"]["full"]["avg_nll"]["zm"],
                    "full_step2_help": int(row["step2"]["full"]["helps"]),
                    "full_step2_gain": row["step2"]["full"]["gain"],
                }
            )


def build_report(summary: dict) -> str:
    lines = [
        "# Angle 3: Asymmetric Transition Gate",
        "",
        "- The Seam-Mandelbrot PMF is extended with a two-sided logistic gate: separate left and right widths `w_left`, `w_right` around the seam.",
        "- The rest of the PMF is unchanged, and the same soft-k lambda sweep is used to select the per-corpus prior strength.",
        f"- median width ratio w_right / w_left: `{summary['width_ratio']['median']:.6f}`",
        f"- mean width ratio w_right / w_left: `{summary['width_ratio']['mean']:.6f}`",
        f"- corpora with material asymmetry (>1.5x width imbalance): `{summary['width_ratio']['material_asymmetry_count']}` / 25",
        "",
    ]
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        info = summary["cutoffs"][label]
        lines.extend(
            [
                f"## {label}",
                "",
                f"- asymmetric seam beats MOE: `{info['asym_beats_moe']}` / 25",
                f"- asymmetric seam beats ZM: `{info['asym_beats_zm']}` / 25",
                f"- asymmetric seam beats symmetric soft-k: `{info['asym_beats_softk']}` / 25",
                f"- asymmetric seam beats Zipf: `{info['asym_beats_zipf']}` / 25",
                f"- median asym minus MOE held-out avg NLL: `{info['median_asym_minus_moe']:.12f}`",
                f"- median asym minus ZM held-out avg NLL: `{info['median_asym_minus_zm']:.12f}`",
                f"- median asym minus symmetric soft-k held-out avg NLL: `{info['median_asym_minus_softk']:.12f}`",
                f"- winner counts: `{json.dumps(info['winner_counts'])}`",
                f"- asymmetric seam step-2 help count: `{info['step2_help_count']}` / 25",
                "",
            ]
        )
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    base_map = load_map(args.base_summary)
    soft_map = load_map(args.softk_summary)
    reranked_map = load_map(RERANKED_SUMMARY_PATH)
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(
            analyze_corpus(
                spec,
                base_row=base_map[spec["slug"]],
                soft_row=soft_map[spec["slug"]],
                reranked_row=reranked_map.get(spec["slug"]),
                corpus_index=corpus_index,
            )
        )
    summary = summarize(rows)
    summary["metadata"] = {
        "softk_summary": str(args.softk_summary),
        "base_summary": str(args.base_summary),
        "canonical_heldout_protocol": "soft-k comparator uses split-fit params only",
    }
    (args.outdir / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (args.outdir / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, args.outdir / "asymmetric_gate_table.csv")


if __name__ == "__main__":
    main()
