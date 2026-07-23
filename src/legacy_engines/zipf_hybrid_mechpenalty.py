from __future__ import annotations

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
HYBRID_PATH = ROOT / "zipf_hybrid_headtail.py"
ANGLE1_PATH = ROOT / "zipf_angle1_head_windows.py"
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"
HYBRID_SUMMARY_PATH = ROOT / "results" / "zipf_hybrid_headtail" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_hybrid_mechpenalty"

TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
LAMBDA_MECHS = [0.01, 0.1, 1.0, 10.0, 100.0]
TRAIN_STARTS = 2
FULL_STARTS = 2
POS_MIN_COUNT = 1


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_hybrid_mech_common")
base = load_module(BASE_PATH, "zipf_hybrid_mech_base")
hybrid_mod = load_module(HYBRID_PATH, "zipf_hybrid_mech_hybrid")
angle1 = load_module(ANGLE1_PATH, "zipf_hybrid_mech_angle1")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_hybrid_mech_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_map(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def precompute_step2_values(log_rank: np.ndarray) -> list[dict]:
    x = common.normalize_x(log_rank)
    target = np.zeros_like(x, dtype=np.float64)
    vocab0 = base.eml.initial_vocabulary(x, target)
    step1 = base.eml.generate_candidates(vocab0, target, 1)
    step1 = base.eml.dedupe_candidates(step1)
    step1 = base.eml.filter_candidates(step1, angle1.STEP2_VARIANCE_THRESHOLD)
    vocab = vocab0 + step1
    step2 = base.eml.generate_candidates(vocab, target, 2)
    step2 = base.eml.dedupe_candidates(step2)
    step2 = base.eml.filter_candidates(step2, angle1.STEP2_VARIANCE_THRESHOLD)
    return [{"expr": row["expr"], "values": np.asarray(row["values"], dtype=np.float64)} for row in step2]


def best_step2_gain(y: np.ndarray, prediction_log: np.ndarray, candidates: list[dict]) -> dict:
    baseline_rmse = common.rmse(y, prediction_log)
    best = None
    best_rmse = float("inf")
    for cand in candidates:
        comp = base.eml.rmse(y, prediction_log + cand["values"])
        if comp < best_rmse:
            best_rmse = comp
            best = cand["expr"]
    gain = float(baseline_rmse - best_rmse)
    return {"expr": best, "gain": gain, "helps": bool(gain > 1e-12)}


def fit_mech_hybrid(
    train_counts: np.ndarray,
    ranks: np.ndarray,
    vocab_size: int,
    hybrid_row: dict,
    base_row: dict,
    softk_row: dict,
    lambda_mech: float,
    train_log_rank_pos: np.ndarray,
    train_log_freq_pos: np.ndarray,
    train_step2_candidates: list[dict],
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
    lam_k = protocol.selection_lambda(hybrid_row)
    target = hybrid_mod.softk_mod.SMOOTH_ALPHA * math.log(vocab_size)
    n_tokens = max(float(np.sum(train_counts)), 1.0)

    starts = []
    if prev_params is not None:
        starts.append(np.asarray(prev_params, dtype=np.float64))
    hparams = protocol.split_fit_params(hybrid_row)
    starts.append(np.array([hparams["b_head"], hparams["c_head"], hparams["alpha_tail"], hparams["beta_tail"], hparams["k"], hparams["w"]], dtype=np.float64))
    sparams = protocol.split_fit_params(softk_row)
    mparams = base_row["models"]["moe"]["params"]
    starts.append(np.array([sparams["b1"], sparams["c1"], mparams["alpha"], mparams["beta"], sparams["k"], sparams["w"]], dtype=np.float64))

    best = None

    def objective(theta: np.ndarray) -> float:
        params = np.clip(theta, lower, upper)
        logpmf = hybrid_mod.hybrid_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        if not math.isfinite(ll):
            return math.inf
        pred_train_log = math.log(float(np.sum(train_counts))) + logpmf[train_counts > 0]
        gain = best_step2_gain(train_log_freq_pos, pred_train_log, train_step2_candidates)["gain"]
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens
        mech_pen = lambda_mech * max(0.0, gain) ** 2
        return float(avg_nll + lam_k * dev * dev + mech_pen)

    for start_index, start in enumerate(starts[:TRAIN_STARTS], start=1):
        start = np.clip(start, lower, upper)
        result = minimize(
            objective,
            x0=start,
            method="Powell",
            bounds=bounds,
            options={"maxiter": 80, "maxfev": 2500, "xtol": 1e-3, "ftol": 1e-4},
        )
        params = np.clip(result.x, lower, upper)
        logpmf = hybrid_mod.hybrid_logpmf(params, ranks)
        ll = base.multinomial_loglike(train_counts, logpmf)
        pred_train_log = math.log(float(np.sum(train_counts))) + logpmf[train_counts > 0]
        step2 = best_step2_gain(train_log_freq_pos, pred_train_log, train_step2_candidates)
        dev = math.log(float(params[4])) - target
        avg_nll = -ll / n_tokens if math.isfinite(ll) else math.inf
        penalized = avg_nll + lam_k * dev * dev + lambda_mech * max(0.0, step2["gain"]) ** 2 if math.isfinite(avg_nll) else math.inf
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
            "train_step2_gain": float(step2["gain"]),
            "train_step2_expr": step2["expr"],
            "dev_logk": float(dev),
            "success": bool(result.success),
            "message": str(result.message),
            "nfev": int(getattr(result, "nfev", 0)),
            "start_index": int(start_index),
            "lambda_k": lam_k,
            "lambda_mech": float(lambda_mech),
        }
        if best is None or record["penalized_objective"] < best["penalized_objective"]:
            best = record
    return best


def analyze_corpus(
    spec: dict,
    corpus_index: int,
    base_row: dict,
    softk_row: dict,
    hybrid_row: dict,
) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    ranks = np.arange(1, len(split["train_counts"]) + 1, dtype=np.float64)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    train_pos = train_counts >= POS_MIN_COUNT
    test_pos = test_counts >= POS_MIN_COUNT

    train_log_rank_pos = np.log(ranks[train_pos])
    train_log_freq_pos = np.log(train_counts[train_pos])
    test_log_rank_pos = np.log(ranks[test_pos])
    test_log_freq_pos = np.log(test_counts[test_pos])
    train_step2_candidates = precompute_step2_values(train_log_rank_pos)
    test_step2_candidates = precompute_step2_values(test_log_rank_pos)

    mech_rows = []
    prev_params = None
    for lam in LAMBDA_MECHS:
        fit = fit_mech_hybrid(
            train_counts=train_counts,
            ranks=ranks,
            vocab_size=split["vocab_size"],
            hybrid_row=hybrid_row,
            base_row=base_row,
            softk_row=softk_row,
            lambda_mech=lam,
            train_log_rank_pos=train_log_rank_pos,
            train_log_freq_pos=train_log_freq_pos,
            train_step2_candidates=train_step2_candidates,
            prev_params=prev_params,
        )
        prev_params = fit["params_vec"]
        evals = base.evaluate_test_fit(fit, test_counts)
        pred_test_log = math.log(float(np.sum(test_counts))) + fit["logpmf"][test_pos]
        heldout_step2 = best_step2_gain(test_log_freq_pos, pred_test_log, test_step2_candidates)
        fit_public = {
            key: val
            for key, val in fit.items()
            if key not in {"logpmf", "params_vec"}
        }
        mech_rows.append(
            {
                "lambda_mech": float(lam),
                "fit": fit_public,
                "test_avg_nll": float(evals["test_avg_nll"]),
                "test_loglike": float(evals["test_loglike"]),
                "heldout_step2_gain": float(heldout_step2["gain"]),
                "heldout_step2_helps": bool(heldout_step2["helps"]),
                "heldout_step2_expr": heldout_step2["expr"],
            }
        )

    return {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": int(dataset["unique_words"]),
        "rows": mech_rows,
    }


def summarize(corpus_rows: list[dict], base_map: dict[str, dict], softk_map: dict[str, dict]) -> dict:
    lambda_summary = {}
    for lam in LAMBDA_MECHS:
        four_way = {"zipf": 0, "zm": 0, "moe": 0, "mech": 0}
        beats_moe = 0
        beats_softk = 0
        beats_zm = 0
        deltas_moe = []
        deltas_zm = []
        deltas_softk = []
        step2_help = 0
        for crow in corpus_rows:
            slug = crow["slug"]
            mech_row = next(row for row in crow["rows"] if abs(row["lambda_mech"] - lam) < 1e-12)
            mech_nll = mech_row["test_avg_nll"]
            zipf_nll = float(base_map[slug]["models"]["zipf"]["test_avg_nll"])
            zm_nll = float(base_map[slug]["models"]["zm"]["test_avg_nll"])
            moe_nll = float(base_map[slug]["models"]["moe"]["test_avg_nll"])
            softk_nll = float(softk_map[slug]["best_lambda"]["test_avg_nll"])
            vals4 = {"zipf": zipf_nll, "zm": zm_nll, "moe": moe_nll, "mech": mech_nll}
            four_way[min(vals4, key=vals4.get)] += 1
            beats_moe += mech_nll < moe_nll
            beats_softk += mech_nll < softk_nll
            beats_zm += mech_nll < zm_nll
            deltas_moe.append(mech_nll - moe_nll)
            deltas_zm.append(mech_nll - zm_nll)
            deltas_softk.append(mech_nll - softk_nll)
            step2_help += int(mech_row["heldout_step2_helps"])
        lambda_summary[f"{lam:g}"] = {
            "lambda_mech": float(lam),
            "beats_moe": int(beats_moe),
            "beats_softk": int(beats_softk),
            "beats_zm": int(beats_zm),
            "winner_counts": four_way,
            "median_minus_moe": float(np.median(deltas_moe)),
            "median_minus_zm": float(np.median(deltas_zm)),
            "median_minus_softk": float(np.median(deltas_softk)),
            "heldout_step2_help_count": int(step2_help),
        }
    candidates = [
        info
        for info in lambda_summary.values()
        if info["heldout_step2_help_count"] <= 4 and info["beats_moe"] >= 18
    ]
    best_candidate = max(candidates, key=lambda info: (info["beats_moe"], info["beats_softk"], -info["median_minus_moe"])) if candidates else None
    return {"corpora": corpus_rows, "lambdas": lambda_summary, "replacement_candidate": best_candidate}


def write_csv(corpus_rows: list[dict], path: Path):
    fieldnames = [
        "slug",
        "name",
        "lambda_mech",
        "test_avg_nll",
        "train_step2_gain",
        "heldout_step2_gain",
        "heldout_step2_help",
        "k",
        "w",
        "alpha_tail",
        "beta_tail",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for crow in corpus_rows:
            for row in crow["rows"]:
                p = row["fit"]["params"]
                writer.writerow(
                    {
                        "slug": crow["slug"],
                        "name": crow["name"],
                        "lambda_mech": row["lambda_mech"],
                        "test_avg_nll": row["test_avg_nll"],
                        "train_step2_gain": row["fit"]["train_step2_gain"],
                        "heldout_step2_gain": row["heldout_step2_gain"],
                        "heldout_step2_help": int(row["heldout_step2_helps"]),
                        "k": p["k"],
                        "w": p["w"],
                        "alpha_tail": p["alpha_tail"],
                        "beta_tail": p["beta_tail"],
                    }
                )


def build_report(summary: dict) -> str:
    lines = [
        "# Experiment A: Mechanism-Preserving Hybrid",
        "",
        "- Hybrid head-ZM + tail-MOE PMF refit with an added in-sample step-2 gain penalty.",
        "- Objective adds `lambda_mech * max(0, train_step2_gain)^2` on top of the existing hybrid k-regularization.",
        "- Held-out step-2 help is evaluated on held-out residuals using the train-fitted PMF and held-out positive-count ranks.",
        "",
    ]
    for key in [f"{lam:g}" for lam in LAMBDA_MECHS]:
        info = summary["lambdas"][key]
        lines.extend(
            [
                f"## lambda_mech = {key}",
                "",
                f"- held-out wins vs MOE: `{info['beats_moe']}` / 25",
                f"- held-out wins vs soft-k: `{info['beats_softk']}` / 25",
                f"- held-out wins vs ZM: `{info['beats_zm']}` / 25",
                f"- four-way winner counts: `{json.dumps(info['winner_counts'])}`",
                f"- median mech minus MOE held-out avg NLL: `{info['median_minus_moe']:.12f}`",
                f"- median mech minus ZM held-out avg NLL: `{info['median_minus_zm']:.12f}`",
                f"- median mech minus soft-k held-out avg NLL: `{info['median_minus_softk']:.12f}`",
                f"- held-out step-2 help count: `{info['heldout_step2_help_count']}` / 25",
                "",
            ]
        )
    if summary["replacement_candidate"] is None:
        lines.append("- No lambda_mech satisfied the decision rule `(held-out step-2 help <= 4/25) AND (beats MOE >= 18/25)`.")
    else:
        info = summary["replacement_candidate"]
        lines.append(
            f"- Replacement candidate found at lambda_mech=`{info['lambda_mech']}` with beats-MOE `{info['beats_moe']}` and held-out step-2 help `{info['heldout_step2_help_count']}`."
        )
    lines.append("")
    return "\n".join(lines)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    base_map = load_map(BASE_SUMMARY_PATH)
    softk_map = load_map(SOFTK_SUMMARY_PATH)
    hybrid_map = load_map(HYBRID_SUMMARY_PATH)
    corpus_rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        corpus_rows.append(
            analyze_corpus(
                spec=spec,
                corpus_index=corpus_index,
                base_row=base_map[spec["slug"]],
                softk_row=softk_map[spec["slug"]],
                hybrid_row=hybrid_map[spec["slug"]],
            )
        )
    summary = summarize(corpus_rows, base_map, softk_map)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(corpus_rows, OUTDIR / "hybrid_mechpenalty_table.csv")


if __name__ == "__main__":
    main()
