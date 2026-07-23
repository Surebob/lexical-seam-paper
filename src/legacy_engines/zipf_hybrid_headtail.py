from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
BASE_PATH = ROOT / "zipf_seam_mandelbrot_pmf.py"
SOFTK_PATH = ROOT / "zipf_seam_mandelbrot_softk.py"
ANGLE1_PATH = ROOT / "zipf_angle1_head_windows.py"
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_hybrid_headtail"

TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
CUTOFFS = [50, 100, 200, 500, 1000, None]
TRAIN_STARTS = 8
FULL_STARTS = 5


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_hybrid_common")
base = load_module(BASE_PATH, "zipf_hybrid_base")
softk_mod = load_module(SOFTK_PATH, "zipf_hybrid_softk")
angle1 = load_module(ANGLE1_PATH, "zipf_hybrid_angle1")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_hybrid_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def parse_args():
    p = argparse.ArgumentParser(description="Hybrid head-tail comparison with explicit split-fit soft-k baseline")
    p.add_argument("--softk-summary", type=Path, default=SOFTK_SUMMARY_PATH)
    p.add_argument("--base-summary", type=Path, default=BASE_SUMMARY_PATH)
    p.add_argument("--outdir", type=Path, default=OUTDIR)
    return p.parse_args()


def cutoff_label(cutoff: int | None) -> str:
    return "full" if cutoff is None else f"top{cutoff}"


def load_map(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def hybrid_logpmf(params: np.ndarray, ranks: np.ndarray) -> np.ndarray | None:
    b_head, c_head, alpha_tail, beta_tail, k, w = [float(v) for v in params]
    if b_head <= 0.0 or c_head < 0.0 or alpha_tail <= 1.000001 or beta_tail <= 0.0 or k <= 0.0 or w <= 0.0:
        return None
    head_logpmf = base.zm_logpmf(b_head, c_head, ranks)
    tail_logpmf = base.moe_logpmf(alpha_tail, beta_tail, ranks)
    if head_logpmf is None or tail_logpmf is None:
        return None
    sigma = base.sigma_curve(ranks, k, w)
    raw = sigma * np.exp(head_logpmf) + (1.0 - sigma) * np.exp(tail_logpmf)
    if np.any(~np.isfinite(raw)) or np.any(raw <= 0.0):
        return None
    log_raw = np.log(raw)
    return log_raw - logsumexp(log_raw)


def hybrid_starts(vocab_size: int, base_row: dict, soft_row: dict, prev_params: list[float] | None, seed: int, n_starts: int) -> list[np.ndarray]:
    k_lo, k_hi = base.k_bounds(vocab_size)
    lower = np.array(
        [base.ZM_B_BOUNDS[0], base.ZM_C_BOUNDS[0], base.MOE_ALPHA_BOUNDS[0], base.MOE_BETA_BOUNDS[0], k_lo, base.SEAM_W_BOUNDS[0]],
        dtype=np.float64,
    )
    upper = np.array(
        [base.ZM_B_BOUNDS[1], base.ZM_C_BOUNDS[1], base.MOE_ALPHA_BOUNDS[1], base.MOE_BETA_BOUNDS[1], k_hi, base.SEAM_W_BOUNDS[1]],
        dtype=np.float64,
    )
    starts = []
    if prev_params is not None:
        starts.append(np.asarray(prev_params, dtype=np.float64))

    base_seam = base_row["models"]["seam"]["params"]
    base_moe = base_row["models"]["moe"]["params"]
    soft_params = protocol.split_fit_params(soft_row)
    starts.append(
        np.array(
            [
                float(base_seam["b1"]),
                float(base_seam["c1"]),
                float(base_moe["alpha"]),
                float(base_moe["beta"]),
                float(base_seam["k"]),
                float(base_seam["w"]),
            ],
            dtype=np.float64,
        )
    )
    starts.append(
        np.array(
            [
                float(soft_params["b1"]),
                float(soft_params["c1"]),
                float(base_moe["alpha"]),
                float(base_moe["beta"]),
                float(soft_params["k"]),
                float(soft_params["w"]),
            ],
            dtype=np.float64,
        )
    )
    starts.append(
        np.array(
            [
                float(base_row["models"]["zm"]["params"]["b"]),
                float(base_row["models"]["zm"]["params"]["c"]),
                float(base_moe["alpha"]),
                float(base_moe["beta"]),
                float(min(max(vocab_size ** softk_mod.SMOOTH_ALPHA, k_lo), k_hi)),
                1.0,
            ],
            dtype=np.float64,
        )
    )
    rng = np.random.default_rng(seed)
    while len(starts) < n_starts:
        starts.append(
            np.array(
                [
                    rng.uniform(*base.ZM_B_BOUNDS),
                    rng.uniform(*base.ZM_C_BOUNDS),
                    rng.uniform(*base.MOE_ALPHA_BOUNDS),
                    float(np.exp(rng.uniform(math.log(base.MOE_BETA_BOUNDS[0]), math.log(base.MOE_BETA_BOUNDS[1])))),
                    float(np.exp(rng.uniform(math.log(k_lo), math.log(k_hi)))),
                    float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                ],
                dtype=np.float64,
            )
        )
    return [np.clip(start, lower, upper) for start in starts[:n_starts]]


def fit_hybrid(
    train_counts: np.ndarray,
    ranks: np.ndarray,
    vocab_size: int,
    base_row: dict,
    soft_row: dict,
    lam: float,
    seed: int,
    n_starts: int,
    prev_params: list[float] | None = None,
) -> dict:
    bounds = [
        base.ZM_B_BOUNDS,
        base.ZM_C_BOUNDS,
        base.MOE_ALPHA_BOUNDS,
        base.MOE_BETA_BOUNDS,
        base.k_bounds(vocab_size),
        base.SEAM_W_BOUNDS,
    ]
    lower = np.array([b[0] for b in bounds], dtype=np.float64)
    upper = np.array([b[1] for b in bounds], dtype=np.float64)
    target = softk_mod.SMOOTH_ALPHA * math.log(vocab_size)
    n_tokens = max(float(np.sum(train_counts)), 1.0)
    best = None

    def objective(theta: np.ndarray) -> float:
        params = np.clip(theta, lower, upper)
        logpmf = hybrid_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        if not math.isfinite(ll):
            return math.inf
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens
        return float(avg_nll + lam * dev * dev)

    for start_index, start in enumerate(hybrid_starts(vocab_size, base_row, soft_row, prev_params, seed, n_starts), start=1):
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = hybrid_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens if math.isfinite(ll) else math.inf
        penalized = avg_nll + lam * dev * dev if math.isfinite(avg_nll) else math.inf
        record = {
            "params": {
                "b_head": float(params[0]),
                "c_head": float(params[1]),
                "alpha_tail": float(params[2]),
                "beta_tail": float(params[3]),
                "k": float(params[4]),
                "w": float(params[5]),
                "transition_fraction": float(np.log(params[4]) / np.log(vocab_size)),
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
            "p": 6,
            "lambda": float(lam),
        }
        if best is None or record["penalized_objective"] < best["penalized_objective"]:
            best = record

    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def full_context_hybrid(dataset: dict, base_row: dict, soft_row: dict, lam: float, seed: int, prev_params: list[float] | None = None) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    fit = fit_hybrid(counts, ranks, int(dataset["unique_words"]), base_row, soft_row, lam, seed, FULL_STARTS, prev_params=prev_params)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {
        "fit": fit,
        "rmse": base.full_rank_rmse(counts, fit["logpmf"]),
        "prediction_log": pred_log,
    }


def analyze_corpus(spec: dict, base_row: dict, soft_row: dict, corpus_index: int) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    ranks = np.arange(1, len(train_counts) + 1, dtype=np.float64)
    vocab_size = split["vocab_size"]

    lambda_rows = []
    prev_params = None
    for lambda_index, lam in enumerate(softk_mod.LAMBDAS, start=1):
        fit = fit_hybrid(
            train_counts,
            ranks,
            vocab_size,
            base_row=base_row,
            soft_row=soft_row,
            lam=lam,
            seed=SPLIT_SEED + 15000 + corpus_index * 50 + lambda_index,
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
    best_full = full_context_hybrid(
        dataset,
        base_row=base_row,
        soft_row=soft_row,
        lam=float(best_lambda_row["lambda"]),
        seed=SPLIT_SEED + 18000 + corpus_index,
        prev_params=best_lambda_row["fit"]["params_vec"],
    )
    best_step2 = base.exact_step2_on_residual(dataset, best_full["prediction_log"])

    model_logpmfs = {
        "zipf": angle1.logpmf_from_params("zipf", base_row["models"]["zipf"]["params"], ranks),
        "zm": angle1.logpmf_from_params("zm", base_row["models"]["zm"]["params"], ranks),
        "moe": angle1.logpmf_from_params("moe", base_row["models"]["moe"]["params"], ranks),
        "softk": angle1.logpmf_from_params("softk", protocol.split_fit_params(soft_row), ranks),
        "hybrid": best_lambda_row["fit"]["logpmf"],
    }
    heldout = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        vals = {name: angle1.restricted_avg_nll(test_counts, logpmf, cutoff) for name, logpmf in model_logpmfs.items()}
        winner = min(vals, key=vals.get)
        heldout[label] = {
            "avg_nll": vals,
            "winner": winner,
            "hybrid_minus_moe": float(vals["hybrid"] - vals["moe"]),
            "hybrid_minus_zm": float(vals["hybrid"] - vals["zm"]),
            "hybrid_minus_softk": float(vals["hybrid"] - vals["softk"]),
            "hybrid_minus_zipf": float(vals["hybrid"] - vals["zipf"]),
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
                "step2_gain": float(best_step2["gain"]),
                "step2_helps": bool(best_step2["helps"]),
                "params": best_full["fit"]["params"],
            },
            "protocol_note": "selection.* are split-fit quantities and are canonical for held-out comparisons; full_refit.* are full-corpus diagnostics only.",
        },
        "heldout": heldout,
        "step2": step2,
        "delta_best_vs_softk_full": float(heldout["full"]["hybrid_minus_softk"]),
    }


def summarize(rows: list[dict]) -> dict:
    cutoff_summary = {}
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        winner_counts = {"zipf": 0, "zm": 0, "moe": 0, "softk": 0, "hybrid": 0}
        hybrid_beats_moe = 0
        hybrid_beats_zm = 0
        hybrid_beats_softk = 0
        hybrid_beats_zipf = 0
        deltas_moe = []
        deltas_zm = []
        deltas_softk = []
        step2_help = 0
        for row in rows:
            info = row["heldout"][label]
            winner_counts[info["winner"]] += 1
            hybrid_beats_moe += info["avg_nll"]["hybrid"] < info["avg_nll"]["moe"]
            hybrid_beats_zm += info["avg_nll"]["hybrid"] < info["avg_nll"]["zm"]
            hybrid_beats_softk += info["avg_nll"]["hybrid"] < info["avg_nll"]["softk"]
            hybrid_beats_zipf += info["avg_nll"]["hybrid"] < info["avg_nll"]["zipf"]
            deltas_moe.append(info["hybrid_minus_moe"])
            deltas_zm.append(info["hybrid_minus_zm"])
            deltas_softk.append(info["hybrid_minus_softk"])
            step2_help += int(row["step2"][label]["helps"])
        cutoff_summary[label] = {
            "winner_counts": winner_counts,
            "hybrid_beats_moe": int(hybrid_beats_moe),
            "hybrid_beats_zm": int(hybrid_beats_zm),
            "hybrid_beats_softk": int(hybrid_beats_softk),
            "hybrid_beats_zipf": int(hybrid_beats_zipf),
            "median_hybrid_minus_moe": float(np.median(deltas_moe)),
            "median_hybrid_minus_zm": float(np.median(deltas_zm)),
            "median_hybrid_minus_softk": float(np.median(deltas_softk)),
            "step2_help_count": int(step2_help),
        }
    return {"rows": rows, "cutoffs": cutoff_summary}


def write_csv(rows: list[dict], path: Path):
    fieldnames = ["slug", "name", "lambda", "k", "w", "top100_hybrid", "top100_moe", "top100_softk", "full_hybrid", "full_moe", "full_softk", "full_step2_help", "full_step2_gain"]
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
                    "k": params["k"],
                    "w": params["w"],
                    "top100_hybrid": row["heldout"]["top100"]["avg_nll"]["hybrid"],
                    "top100_moe": row["heldout"]["top100"]["avg_nll"]["moe"],
                    "top100_softk": row["heldout"]["top100"]["avg_nll"]["softk"],
                    "full_hybrid": row["heldout"]["full"]["avg_nll"]["hybrid"],
                    "full_moe": row["heldout"]["full"]["avg_nll"]["moe"],
                    "full_softk": row["heldout"]["full"]["avg_nll"]["softk"],
                    "full_step2_help": int(row["best_lambda"]["full_refit"]["step2_helps"]),
                    "full_step2_gain": row["best_lambda"]["full_refit"]["step2_gain"],
                }
            )


def build_report(summary: dict) -> str:
    lines = [
        "# Hybrid Head-Tail PMF",
        "",
        "- Head component: normalized ZM law over ranks.",
        "- Tail component: normalized MOEZipf law over ranks.",
        "- The two are blended by a rank-dependent seam gate and then renormalized into one discrete PMF.",
        "- The same soft-k lambda sweep is used to regularize seam location `k`.",
        "",
    ]
    for cutoff in CUTOFFS:
        label = cutoff_label(cutoff)
        info = summary["cutoffs"][label]
        lines.extend(
            [
                f"## {label}",
                "",
                f"- hybrid beats MOE: `{info['hybrid_beats_moe']}` / 25",
                f"- hybrid beats ZM: `{info['hybrid_beats_zm']}` / 25",
                f"- hybrid beats soft-k seam: `{info['hybrid_beats_softk']}` / 25",
                f"- hybrid beats Zipf: `{info['hybrid_beats_zipf']}` / 25",
                f"- median hybrid minus MOE held-out avg NLL: `{info['median_hybrid_minus_moe']:.12f}`",
                f"- median hybrid minus ZM held-out avg NLL: `{info['median_hybrid_minus_zm']:.12f}`",
                f"- median hybrid minus soft-k held-out avg NLL: `{info['median_hybrid_minus_softk']:.12f}`",
                f"- winner counts: `{json.dumps(info['winner_counts'])}`",
                f"- hybrid step-2 help count: `{info['step2_help_count']}` / 25",
                "",
            ]
        )
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    base_map = load_map(args.base_summary)
    soft_map = load_map(args.softk_summary)
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(analyze_corpus(spec, base_map[spec["slug"]], soft_map[spec["slug"]], corpus_index))
    summary = summarize(rows)
    summary["metadata"] = {
        "softk_summary": str(args.softk_summary),
        "base_summary": str(args.base_summary),
        "canonical_heldout_protocol": "soft-k comparator uses split-fit params only",
    }
    (args.outdir / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (args.outdir / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, args.outdir / "hybrid_headtail_table.csv")


if __name__ == "__main__":
    main()
