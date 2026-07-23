from __future__ import annotations

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
PROTOCOL_UTILS_PATH = ROOT / "zipf_eval_protocol_utils.py"
OUTDIR = ROOT / "results" / "zipf_fully_ours_variants"
BASE_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_pmf" / "summary.json"
SOFTK_SUMMARY_PATH = ROOT / "results" / "zipf_seam_mandelbrot_softk" / "summary.json"

TRAIN_FRACTION = 0.8
SPLIT_SEED = 20260416
TRAIN_STARTS = 6
FULL_STARTS = 4


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_ours_common")
base = load_module(BASE_PATH, "zipf_ours_base")
softk_mod = load_module(SOFTK_PATH, "zipf_ours_softk")
protocol = load_module(PROTOCOL_UTILS_PATH, "zipf_ours_protocol")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_map(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
    return {row["slug"]: row for row in rows}


def sigma_curve(ranks: np.ndarray, k: float, w: float) -> np.ndarray:
    return base.sigma_curve(ranks, k, w)


def zm_scores(ranks: np.ndarray, b: float, c: float) -> np.ndarray | None:
    if b <= 0.0 or c < 0.0:
        return None
    scores = -b * np.log(ranks + c)
    return scores if np.all(np.isfinite(scores)) else None


def nested_headtail_logpmf(params: np.ndarray, ranks: np.ndarray, vocab_size: int) -> np.ndarray | None:
    b1, c1, b2, c2, b3, c3, k, w = [float(v) for v in params]
    if min(b1, b2, b3) <= 0.0 or min(c1, c2, c3) < 0.0 or k <= 0.0 or w <= 0.0:
        return None
    z1 = zm_scores(ranks, b1, c1)
    z2 = zm_scores(ranks, b2, c2)
    z3 = zm_scores(ranks, b3, c3)
    if z1 is None or z2 is None or z3 is None:
        return None
    k0 = max(5.0, 0.5 * math.log(vocab_size))
    w0 = 0.22
    sigma0 = sigma_curve(ranks, k0, w0)
    sigma_main = sigma_curve(ranks, k, w)
    head_comp = sigma0 * z1 + (1.0 - sigma0) * z2
    raw = sigma_main * np.exp(head_comp) + (1.0 - sigma_main) * np.exp(z3)
    if np.any(~np.isfinite(raw)) or np.any(raw <= 0.0):
        return None
    log_raw = np.log(raw)
    return log_raw - logsumexp(log_raw)


def three_regime_logpmf(params: np.ndarray, ranks: np.ndarray) -> np.ndarray | None:
    b1, c1, b2, c2, b3, c3, k1, w1, k2, w2 = [float(v) for v in params]
    if min(b1, b2, b3) <= 0.0 or min(c1, c2, c3) < 0.0 or min(k1, w1, k2, w2) <= 0.0 or not (k1 < k2):
        return None
    z1 = zm_scores(ranks, b1, c1)
    z2 = zm_scores(ranks, b2, c2)
    z3 = zm_scores(ranks, b3, c3)
    if z1 is None or z2 is None or z3 is None:
        return None
    s1 = sigma_curve(ranks, k1, w1)
    s2 = sigma_curve(ranks, k2, w2)
    w_head = s1
    w_mid = (1.0 - s1) * s2
    w_tail = 1.0 - s2
    raw = w_head * np.exp(z1) + w_mid * np.exp(z2) + w_tail * np.exp(z3)
    if np.any(~np.isfinite(raw)) or np.any(raw <= 0.0):
        return None
    log_raw = np.log(raw)
    return log_raw - logsumexp(log_raw)


def fit_variant(
    variant: str,
    train_counts: np.ndarray,
    ranks: np.ndarray,
    vocab_size: int,
    base_row: dict,
    softk_row: dict,
    seed: int,
    n_starts: int,
) -> dict:
    rng = np.random.default_rng(seed)
    soft = protocol.split_fit_params(softk_row)
    zm = base_row["models"]["zm"]["params"]
    k_lo, k_hi = base.k_bounds(vocab_size)

    if variant == "nested":
        bounds = [
            base.SEAM_B_BOUNDS,
            base.SEAM_C_BOUNDS,
            base.SEAM_B_BOUNDS,
            base.SEAM_C_BOUNDS,
            base.SEAM_B_BOUNDS,
            base.SEAM_C_BOUNDS,
            (k_lo, k_hi),
            base.SEAM_W_BOUNDS,
        ]
        starts = [
            np.array([soft["b1"], soft["c1"], zm["b"], zm["c"], soft["b2"], soft["c2"], soft["k"], soft["w"]], dtype=np.float64),
            np.array([zm["b"], zm["c"], soft["b1"], soft["c1"], soft["b2"], soft["c2"], max(k_lo, vocab_size ** softk_mod.SMOOTH_ALPHA), soft["w"]], dtype=np.float64),
        ]
        while len(starts) < n_starts:
            starts.append(
                np.array(
                    [
                        rng.uniform(*base.SEAM_B_BOUNDS),
                        rng.uniform(*base.SEAM_C_BOUNDS),
                        rng.uniform(*base.SEAM_B_BOUNDS),
                        rng.uniform(*base.SEAM_C_BOUNDS),
                        rng.uniform(*base.SEAM_B_BOUNDS),
                        rng.uniform(*base.SEAM_C_BOUNDS),
                        float(np.exp(rng.uniform(math.log(k_lo), math.log(k_hi)))),
                        float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                    ],
                    dtype=np.float64,
                )
            )
        scorer = lambda p: nested_headtail_logpmf(p, ranks, vocab_size)
        pcount = 8
    else:
        early_hi = min(200.0, max(20.0, 0.25 * k_hi))
        bounds = [
            base.SEAM_B_BOUNDS,
            base.SEAM_C_BOUNDS,
            base.SEAM_B_BOUNDS,
            base.SEAM_C_BOUNDS,
            base.SEAM_B_BOUNDS,
            base.SEAM_C_BOUNDS,
            (5.0, early_hi),
            (0.05, 3.0),
            (max(20.0, early_hi), k_hi),
            base.SEAM_W_BOUNDS,
        ]
        starts = [
            np.array([zm["b"], zm["c"], soft["b1"], soft["c1"], soft["b2"], soft["c2"], max(5.0, 0.5 * math.log(vocab_size)), 0.25, soft["k"], soft["w"]], dtype=np.float64),
            np.array([soft["b1"], soft["c1"], zm["b"], zm["c"], soft["b2"], soft["c2"], max(5.0, 0.5 * math.log(vocab_size)), 0.25, soft["k"], soft["w"]], dtype=np.float64),
        ]
        while len(starts) < n_starts:
            k1 = float(np.exp(rng.uniform(math.log(5.0), math.log(early_hi))))
            k2 = float(np.exp(rng.uniform(math.log(max(20.0, early_hi)), math.log(k_hi))))
            if k2 <= k1:
                k2 = min(k_hi, k1 + 10.0)
            starts.append(
                np.array(
                    [
                        rng.uniform(*base.SEAM_B_BOUNDS),
                        rng.uniform(*base.SEAM_C_BOUNDS),
                        rng.uniform(*base.SEAM_B_BOUNDS),
                        rng.uniform(*base.SEAM_C_BOUNDS),
                        rng.uniform(*base.SEAM_B_BOUNDS),
                        rng.uniform(*base.SEAM_C_BOUNDS),
                        k1,
                        float(np.exp(rng.uniform(math.log(0.05), math.log(3.0)))),
                        k2,
                        float(np.exp(rng.uniform(math.log(base.SEAM_W_BOUNDS[0]), math.log(base.SEAM_W_BOUNDS[1])))),
                    ],
                    dtype=np.float64,
                )
            )
        scorer = lambda p: three_regime_logpmf(p, ranks)
        pcount = 10

    lower = np.array([b[0] for b in bounds], dtype=np.float64)
    upper = np.array([b[1] for b in bounds], dtype=np.float64)
    best = None

    def objective(theta: np.ndarray) -> float:
        params = np.clip(theta, lower, upper)
        logpmf = scorer(params)
        ll = base.multinomial_loglike(train_counts, logpmf)
        return -ll if math.isfinite(ll) else math.inf

    for start_index, start in enumerate(starts[:n_starts], start=1):
        start = np.clip(start, lower, upper)
        result = minimize(objective, x0=start, method="L-BFGS-B", bounds=bounds)
        params = np.clip(result.x, lower, upper)
        logpmf = scorer(params)
        ll = base.multinomial_loglike(train_counts, logpmf)
        record = {
            "params_vec": [float(v) for v in params],
            "logpmf": logpmf,
            "train_loglike": float(ll),
            "success": bool(result.success),
            "message": str(result.message),
            "nfev": int(getattr(result, "nfev", 0)),
            "start_index": int(start_index),
            "p": pcount,
        }
        if variant == "nested":
            record["params"] = {
                "b1": float(params[0]),
                "c1": float(params[1]),
                "b2": float(params[2]),
                "c2": float(params[3]),
                "b3": float(params[4]),
                "c3": float(params[5]),
                "k": float(params[6]),
                "w": float(params[7]),
                "transition_fraction": float(np.log(params[6]) / np.log(vocab_size)),
            }
        else:
            record["params"] = {
                "b1": float(params[0]),
                "c1": float(params[1]),
                "b2": float(params[2]),
                "c2": float(params[3]),
                "b3": float(params[4]),
                "c3": float(params[5]),
                "k1": float(params[6]),
                "w1": float(params[7]),
                "k2": float(params[8]),
                "w2": float(params[9]),
                "transition_fraction_2": float(np.log(params[8]) / np.log(vocab_size)),
            }
        if best is None or record["train_loglike"] > best["train_loglike"]:
            best = record
    n = int(np.sum(train_counts))
    best["aic"] = float(2 * best["p"] - 2 * best["train_loglike"])
    best["bic"] = float(math.log(n) * best["p"] - 2 * best["train_loglike"])
    return best


def full_context(variant: str, dataset: dict, base_row: dict, softk_row: dict, seed: int) -> dict:
    counts = dataset["freqs"].astype(np.float64)
    ranks = np.arange(1, len(counts) + 1, dtype=np.float64)
    fit = fit_variant(variant, counts, ranks, int(dataset["unique_words"]), base_row, softk_row, seed, FULL_STARTS)
    pred_log = math.log(float(np.sum(counts))) + fit["logpmf"]
    return {"fit": fit, "rmse": base.full_rank_rmse(counts, fit["logpmf"]), "prediction_log": pred_log}


def analyze_corpus(spec: dict, corpus_index: int, base_row: dict, softk_row: dict) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    split = base.split_counts_by_train_rank(dataset, seed=SPLIT_SEED + corpus_index, train_fraction=TRAIN_FRACTION)
    train_counts = split["train_counts"]
    test_counts = split["test_counts"]
    ranks = np.arange(1, len(train_counts) + 1, dtype=np.float64)

    variants = {}
    for variant, offset in [("nested", 20000), ("three_regime", 30000)]:
        fit = fit_variant(variant, train_counts, ranks, split["vocab_size"], base_row, softk_row, SPLIT_SEED + offset + corpus_index, TRAIN_STARTS)
        evals = base.evaluate_test_fit(fit, test_counts)
        full = full_context(variant, dataset, base_row, softk_row, SPLIT_SEED + offset + 5000 + corpus_index)
        step2 = base.exact_step2_on_residual(dataset, full["prediction_log"])
        variants[variant] = {
            "test_avg_nll": float(evals["test_avg_nll"]),
            "test_loglike": float(evals["test_loglike"]),
            "bic": float(full["fit"]["bic"]),
            "rmse": float(full["rmse"]),
            "step2_gain": float(step2["gain"]),
            "step2_helps": bool(step2["helps"]),
            "params": full["fit"]["params"],
        }
    return {"slug": spec["slug"], "name": spec["name"], "token_count": int(dataset["token_count"]), "vocab_size": int(dataset["unique_words"]), "variants": variants}


def summarize(rows: list[dict], base_map: dict[str, dict], softk_map: dict[str, dict], hybrid_map: dict[str, dict]) -> dict:
    out = {"rows": rows, "variants": {}}
    for variant in ["nested", "three_regime"]:
        beats = {"moe": 0, "zm": 0, "softk": 0}
        winner_counts = {"zipf": 0, "zm": 0, "moe": 0, variant: 0}
        deltas_moe = []
        deltas_zm = []
        deltas_softk = []
        step2_help = 0
        bic_wins_vs_two_regime = 0
        bic_wins_vs_hybrid = 0
        for row in rows:
            slug = row["slug"]
            test_nll = row["variants"][variant]["test_avg_nll"]
            zipf_nll = float(base_map[slug]["models"]["zipf"]["test_avg_nll"])
            zm_nll = float(base_map[slug]["models"]["zm"]["test_avg_nll"])
            moe_nll = float(base_map[slug]["models"]["moe"]["test_avg_nll"])
            softk_nll = float(softk_map[slug]["best_lambda"]["test_avg_nll"])
            vals = {"zipf": zipf_nll, "zm": zm_nll, "moe": moe_nll, variant: test_nll}
            winner_counts[min(vals, key=vals.get)] += 1
            beats["moe"] += test_nll < moe_nll
            beats["zm"] += test_nll < zm_nll
            beats["softk"] += test_nll < softk_nll
            deltas_moe.append(test_nll - moe_nll)
            deltas_zm.append(test_nll - zm_nll)
            deltas_softk.append(test_nll - softk_nll)
            step2_help += int(row["variants"][variant]["step2_helps"])
            bic_wins_vs_two_regime += row["variants"][variant]["bic"] < float(base_map[slug]["models"]["seam"]["bic"])
            bic_wins_vs_hybrid += row["variants"][variant]["bic"] < float(hybrid_map[slug]["best_lambda"]["bic"])
        out["variants"][variant] = {
            "beats_moe": int(beats["moe"]),
            "beats_zm": int(beats["zm"]),
            "beats_softk": int(beats["softk"]),
            "winner_counts": winner_counts,
            "median_minus_moe": float(np.median(deltas_moe)),
            "median_minus_zm": float(np.median(deltas_zm)),
            "median_minus_softk": float(np.median(deltas_softk)),
            "step2_help_count": int(step2_help),
            "bic_wins_vs_two_regime": int(bic_wins_vs_two_regime),
            "bic_wins_vs_hybrid": int(bic_wins_vs_hybrid),
        }
    return out


def write_csv(rows: list[dict], path: Path):
    fieldnames = ["slug", "name", "variant", "test_avg_nll", "bic", "rmse", "step2_gain", "step2_help"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            for variant, info in row["variants"].items():
                writer.writerow(
                    {
                        "slug": row["slug"],
                        "name": row["name"],
                        "variant": variant,
                        "test_avg_nll": info["test_avg_nll"],
                        "bic": info["bic"],
                        "rmse": info["rmse"],
                        "step2_gain": info["step2_gain"],
                        "step2_help": int(info["step2_helps"]),
                    }
                )


def build_report(summary: dict) -> str:
    lines = [
        "# Experiment B: Fully-Ours Variants",
        "",
        "- B1 nested seam: early fixed mini-seam between two ZM components in the head, then a main seam into a tail ZM.",
        "- B2 three-regime: three ZM components blended sequentially by two smooth gates.",
        "",
    ]
    for variant in ["nested", "three_regime"]:
        info = summary["variants"][variant]
        lines.extend(
            [
                f"## {variant}",
                "",
                f"- held-out wins vs MOE: `{info['beats_moe']}` / 25",
                f"- held-out wins vs soft-k: `{info['beats_softk']}` / 25",
                f"- held-out wins vs ZM: `{info['beats_zm']}` / 25",
                f"- four-way winner counts: `{json.dumps(info['winner_counts'])}`",
                f"- median minus MOE held-out avg NLL: `{info['median_minus_moe']:.12f}`",
                f"- median minus ZM held-out avg NLL: `{info['median_minus_zm']:.12f}`",
                f"- median minus soft-k held-out avg NLL: `{info['median_minus_softk']:.12f}`",
                f"- step-2 help count: `{info['step2_help_count']}` / 25",
                f"- BIC wins vs two-regime seam: `{info['bic_wins_vs_two_regime']}` / 25",
                f"- BIC wins vs hybrid head-tail: `{info['bic_wins_vs_hybrid']}` / 25",
                "",
            ]
        )
    return "\n".join(lines)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    base_map = load_map(BASE_SUMMARY_PATH)
    softk_map = load_map(SOFTK_SUMMARY_PATH)
    hybrid_map = load_map(ROOT / "results" / "zipf_hybrid_headtail" / "summary.json")
    rows = []
    for corpus_index, spec in enumerate(common.SEARCHED_CORPORA, start=1):
        rows.append(analyze_corpus(spec, corpus_index, base_map[spec["slug"]], softk_map[spec["slug"]]))
    summary = summarize(rows, base_map, softk_map, hybrid_map)
    (OUTDIR / "summary.json").write_text(json.dumps(sanitize(summary), indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(summary), encoding="utf-8")
    write_csv(rows, OUTDIR / "fully_ours_variants_table.csv")


if __name__ == "__main__":
    main()
