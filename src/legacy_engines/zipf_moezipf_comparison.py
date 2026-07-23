from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import zeta as scipy_zeta


ROOT = Path("/Volumes/External2TB/emlexperiment")
COMMON_PATH = ROOT / "zipf_analysis_common.py"
SQRT_V_SUMMARY_PATH = ROOT / "results" / "zipf_sqrt_v_all_corpora" / "summary.json"
RERANKED_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_moezipf_comparison"

ALPHA_MIN = 1.000001
ALPHA_MAX = 10.0
THETA_ALPHA_BOUNDS = (math.log(1e-4), math.log(ALPHA_MAX - ALPHA_MIN))
THETA_BETA_BOUNDS = (math.log(1e-4), math.log(1e4))
MOBY_GROUP_CUTOFF = 10


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = load_module(COMMON_PATH, "zipf_moezipf_common")


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_reference_maps():
    sqrt_rows = json.loads(SQRT_V_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    reranked_rows = json.loads(RERANKED_SUMMARY_PATH.read_text(encoding="utf-8"))["rows"]
    sqrt_map = {row["slug"]: row for row in sqrt_rows}
    reranked_map = {row["slug"]: row for row in reranked_rows}
    return sqrt_map, reranked_map


def riemann_zeta(alpha: float) -> float:
    return float(scipy_zeta(alpha, 1.0))


def zipf_loglike(freqs: np.ndarray, alpha: float) -> float:
    z = riemann_zeta(alpha)
    if not math.isfinite(z) or z <= 0.0:
        return -math.inf
    return float(-alpha * np.log(freqs).sum() - len(freqs) * math.log(z))


def fit_zipf_mle(freqs: np.ndarray) -> dict:
    def objective(alpha: float):
        if alpha <= ALPHA_MIN or alpha >= ALPHA_MAX:
            return math.inf
        ll = zipf_loglike(freqs, alpha)
        return -ll if math.isfinite(ll) else math.inf

    result = minimize_scalar(objective, bounds=(1.0001, ALPHA_MAX), method="bounded", options={"xatol": 1e-6})
    alpha = float(result.x)
    loglike = float(-result.fun)
    n = int(len(freqs))
    return {
        "alpha": alpha,
        "loglike": loglike,
        "aic": float(2 * 1 - 2 * loglike),
        "bic": float(math.log(n) * 1 - 2 * loglike),
        "success": bool(result.success),
        "message": result.message,
    }


def moe_logpmf(alpha: float, beta: float, x: np.ndarray) -> np.ndarray | None:
    z = riemann_zeta(alpha)
    if not math.isfinite(z) or z <= 0.0 or beta <= 0.0:
        return None
    hz_x = np.asarray(scipy_zeta(alpha, x), dtype=np.float64)
    hz_x1 = np.asarray(scipy_zeta(alpha, x + 1.0), dtype=np.float64)
    beta_bar = 1.0 - beta
    denom_x = z - beta_bar * hz_x
    denom_x1 = z - beta_bar * hz_x1
    if np.any(~np.isfinite(denom_x)) or np.any(~np.isfinite(denom_x1)) or np.any(denom_x <= 0.0) or np.any(denom_x1 <= 0.0):
        return None
    with np.errstate(all="ignore"):
        logp = -alpha * np.log(x) + math.log(beta) + math.log(z) - np.log(denom_x) - np.log(denom_x1)
    if np.any(~np.isfinite(logp)):
        return None
    return logp


def moe_loglike_from_hist(unique_x: np.ndarray, counts: np.ndarray, alpha: float, beta: float) -> float:
    logp = moe_logpmf(alpha, beta, unique_x.astype(np.float64))
    if logp is None:
        return -math.inf
    return float(np.dot(counts.astype(np.float64), logp))


def unpack_theta(theta: np.ndarray) -> tuple[float, float]:
    alpha = ALPHA_MIN + math.exp(float(theta[0]))
    beta = math.exp(float(theta[1]))
    return alpha, beta


def fit_moe_mle(freqs: np.ndarray, zipf_fit: dict) -> dict:
    unique_x, counts = np.unique(freqs.astype(np.int64), return_counts=True)
    starts = []
    for alpha0 in [zipf_fit["alpha"], 1.3, 1.6, 2.0, 2.4]:
        alpha0 = min(max(alpha0, 1.0001), ALPHA_MAX)
        for beta0 in [0.4, 0.7, 1.0, 1.4, 2.0, 4.0]:
            starts.append(np.array([math.log(alpha0 - ALPHA_MIN), math.log(beta0)], dtype=np.float64))

    best = None

    def objective(theta: np.ndarray) -> float:
        alpha, beta = unpack_theta(theta)
        ll = moe_loglike_from_hist(unique_x, counts, alpha, beta)
        return -ll if math.isfinite(ll) else math.inf

    bounds = [THETA_ALPHA_BOUNDS, THETA_BETA_BOUNDS]
    for idx, x0 in enumerate(starts, start=1):
        result = minimize(objective, x0=x0, method="L-BFGS-B", bounds=bounds)
        alpha, beta = unpack_theta(result.x)
        loglike = float(-result.fun) if math.isfinite(result.fun) else -math.inf
        record = {
            "start_index": idx,
            "alpha": alpha,
            "beta": beta,
            "loglike": loglike,
            "success": bool(result.success),
            "status": int(result.status),
            "message": result.message,
            "nfev": int(result.nfev),
        }
        if best is None or record["loglike"] > best["loglike"]:
            best = record

    n = int(len(freqs))
    return {
        **best,
        "aic": float(2 * 2 - 2 * best["loglike"]),
        "bic": float(math.log(n) * 2 - 2 * best["loglike"]),
    }


def zipf_tail_ge(alpha: float, x: np.ndarray) -> np.ndarray:
    z = riemann_zeta(alpha)
    return np.asarray(scipy_zeta(alpha, x), dtype=np.float64) / z


def moe_tail_ge(alpha: float, beta: float, x: np.ndarray) -> np.ndarray:
    base = zipf_tail_ge(alpha, x)
    beta_bar = 1.0 - beta
    return beta * base / (1.0 - beta_bar * base)


def ensure_tail_cutoff(tail_fn, n_obs: int, observed_max: int) -> int:
    xmax = max(int(observed_max), 1)
    while True:
        tail = float(np.asarray(tail_fn(np.array([float(xmax)], dtype=np.float64)))[0])
        if not math.isfinite(tail) or n_obs * tail < 1.0:
            break
        xmax *= 2
        if xmax > 1_000_000:
            break
    return xmax


def rank_curve_from_tail(tail_fn, n_obs: int, observed_max: int) -> np.ndarray:
    xmax = ensure_tail_cutoff(tail_fn, n_obs, observed_max)
    xs = np.arange(1, xmax + 1, dtype=np.float64)
    counts_ge = n_obs * np.asarray(tail_fn(xs), dtype=np.float64)
    if counts_ge[-1] >= 1.0:
        counts_ge = np.append(counts_ge, 0.0)
    ranks = np.arange(1, n_obs + 1, dtype=np.float64)
    predicted = np.searchsorted(-counts_ge, -ranks, side="right").astype(np.float64)
    predicted[predicted < 1.0] = 1.0
    return predicted


def rmse_log_rank(obs_freqs: np.ndarray, pred_freqs: np.ndarray) -> float:
    obs = np.log(np.asarray(obs_freqs, dtype=np.float64))
    pred = np.log(np.asarray(pred_freqs, dtype=np.float64))
    return float(math.sqrt(float(np.mean((obs - pred) ** 2))))


def rmse_bic(rmse_value: float, n_obs: int, p: int) -> float:
    mse = max(float(rmse_value) ** 2, 1e-300)
    return float(p * math.log(n_obs) + n_obs * math.log(mse))


def grouped_chi_square(obs_freqs: np.ndarray, pmf_fn, tail_ge_fn, cutoff: int) -> dict:
    n = int(len(obs_freqs))
    observed_counts = np.array([(obs_freqs == x).sum() for x in range(1, cutoff)], dtype=np.float64)
    observed_tail = float((obs_freqs >= cutoff).sum())
    expected_counts = n * np.asarray(pmf_fn(np.arange(1, cutoff, dtype=np.float64)), dtype=np.float64)
    expected_tail = n * float(tail_ge_fn(np.array([float(cutoff)], dtype=np.float64))[0])
    observed = np.append(observed_counts, observed_tail)
    expected = np.append(expected_counts, expected_tail)
    valid = expected > 0.0
    chi = float(np.sum(((observed[valid] - expected[valid]) ** 2) / expected[valid]))
    return {
        "cutoff": cutoff,
        "observed": observed.tolist(),
        "expected": expected.tolist(),
        "chi_square": chi,
    }


def analyze_corpus(spec: dict, sqrt_ref: dict, reranked_ref: dict) -> dict:
    dataset = common.build_zipf_dataset(common.corpus_path(spec))
    freqs = dataset["freqs"].astype(np.int64)
    n_obs = int(dataset["unique_words"])
    observed_max = int(freqs.max())

    zipf_fit = fit_zipf_mle(freqs)
    moe_fit = fit_moe_mle(freqs, zipf_fit)

    zipf_pred = rank_curve_from_tail(lambda x: zipf_tail_ge(zipf_fit["alpha"], x), n_obs, observed_max)
    moe_pred = rank_curve_from_tail(lambda x: moe_tail_ge(moe_fit["alpha"], moe_fit["beta"], x), n_obs, observed_max)
    zipf_rmse = rmse_log_rank(freqs, zipf_pred)
    moe_rmse = rmse_log_rank(freqs, moe_pred)

    zipf_rank_bic = rmse_bic(zipf_rmse, n_obs, p=1)
    moe_rank_bic = rmse_bic(moe_rmse, n_obs, p=2)
    sqrt_rank_bic = rmse_bic(float(sqrt_ref["sqrt_v_rmse"]), n_obs, p=7)
    reranked_rank_bic = rmse_bic(float(reranked_ref["reranked_rmse"]), n_obs, p=8)

    rank_bics = {
        "zipf": zipf_rank_bic,
        "moezipf": moe_rank_bic,
        "sqrt_v_7param": sqrt_rank_bic,
        "reranked_8param": reranked_rank_bic,
    }
    rank_bic_winner = min(rank_bics, key=rank_bics.get)

    row = {
        "slug": spec["slug"],
        "name": spec["name"],
        "token_count": int(dataset["token_count"]),
        "vocab_size": n_obs,
        "zipf_alpha": float(zipf_fit["alpha"]),
        "zipf_loglike": float(zipf_fit["loglike"]),
        "zipf_aic": float(zipf_fit["aic"]),
        "zipf_bic": float(zipf_fit["bic"]),
        "zipf_rmse": float(zipf_rmse),
        "zipf_rank_bic": float(zipf_rank_bic),
        "moe_alpha": float(moe_fit["alpha"]),
        "moe_beta": float(moe_fit["beta"]),
        "moe_loglike": float(moe_fit["loglike"]),
        "moe_aic": float(moe_fit["aic"]),
        "moe_bic": float(moe_fit["bic"]),
        "moe_rmse": float(moe_rmse),
        "moe_rank_bic": float(moe_rank_bic),
        "sqrt_v_rmse": float(sqrt_ref["sqrt_v_rmse"]),
        "sqrt_v_rank_bic": float(sqrt_rank_bic),
        "reranked_rmse": float(reranked_ref["reranked_rmse"]),
        "reranked_rank_bic": float(reranked_rank_bic),
        "rank_bic_winner": rank_bic_winner,
    }

    if spec["slug"] == "moby_dick":
        zipf_pmf = lambda x: np.exp(-zipf_fit["alpha"] * np.log(x) - math.log(riemann_zeta(zipf_fit["alpha"])))
        moe_pmf = lambda x: np.exp(moe_logpmf(moe_fit["alpha"], moe_fit["beta"], np.asarray(x, dtype=np.float64)))
        zipf_chi = grouped_chi_square(freqs, zipf_pmf, lambda x: zipf_tail_ge(zipf_fit["alpha"], x), MOBY_GROUP_CUTOFF)
        moe_chi = grouped_chi_square(freqs, moe_pmf, lambda x: moe_tail_ge(moe_fit["alpha"], moe_fit["beta"], x), MOBY_GROUP_CUTOFF)
        reduction = 1.0 - (moe_chi["chi_square"] / zipf_chi["chi_square"])
        row["moby_paper_check"] = {
            "paper_reported_reduction": 0.7964,
            "zipf_chi_square_grouped_ge_10": float(zipf_chi["chi_square"]),
            "moe_chi_square_grouped_ge_10": float(moe_chi["chi_square"]),
            "reduction": float(reduction),
        }

    return row


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "slug",
                "name",
                "token_count",
                "vocab_size",
                "zipf_alpha",
                "zipf_loglike",
                "zipf_aic",
                "zipf_bic",
                "zipf_rmse",
                "zipf_rank_bic",
                "moe_alpha",
                "moe_beta",
                "moe_loglike",
                "moe_aic",
                "moe_bic",
                "moe_rmse",
                "moe_rank_bic",
                "sqrt_v_rmse",
                "sqrt_v_rank_bic",
                "reranked_rmse",
                "reranked_rank_bic",
                "rank_bic_winner",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def build_report(rows: list[dict], winner_counts: dict, moby_check: dict | None) -> str:
    lines = [
        "# MOEZipf Prior-Art Comparison",
        "",
        "- MOEZipf source: Pérez-Casany & Casellas (2013), arXiv:1304.4540v2.",
        "- Observations are per-word-type frequencies, matching the paper's `frequency of occurrence of words` setup.",
        "- Likelihood AIC/BIC for Zipf and MOEZipf use the frequency-sample likelihood.",
        "- Cross-model winner comparison uses the same rank-curve RMSE BIC style as the earlier model-comparison runs: `p ln(n) + n ln(MSE)`.",
        "",
        "## Rank-BIC Winner Counts",
        "",
    ]
    for model, count in winner_counts.items():
        lines.append(f"- {model}: `{count}`")

    if moby_check is not None:
        lines.extend(
            [
                "",
                "## Moby Dick Paper Check",
                "",
                f"- paper-reported Pearson chi-square reduction: `{moby_check['paper_reported_reduction']:.4%}`",
                f"- reproduced Zipf chi-square with frequencies `>=10` grouped: `{moby_check['zipf_chi_square_grouped_ge_10']:.6f}`",
                f"- reproduced MOE chi-square with frequencies `>=10` grouped: `{moby_check['moe_chi_square_grouped_ge_10']:.6f}`",
                f"- reproduced reduction: `{moby_check['reduction']:.4%}`",
            ]
        )

    lines.extend(
        [
            "",
            "| Corpus | Zipf RMSE | MOE RMSE | 7-param RMSE | 8-param RMSE | Rank-BIC winner |",
            "| --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['zipf_rmse']:.12f} | {row['moe_rmse']:.12f} | {row['sqrt_v_rmse']:.12f} | {row['reranked_rmse']:.12f} | {row['rank_bic_winner']} |"
        )

    return "\n".join(lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sqrt_map, reranked_map = load_reference_maps()

    rows = []
    moby_check = None
    for spec in common.SEARCHED_CORPORA:
        row = analyze_corpus(spec, sqrt_map[spec["slug"]], reranked_map[spec["slug"]])
        rows.append(row)
        if spec["slug"] == "moby_dick":
            moby_check = row.get("moby_paper_check")

    rows = sorted(rows, key=lambda row: row["moe_bic"])
    winner_counts = {
        "zipf": int(sum(row["rank_bic_winner"] == "zipf" for row in rows)),
        "moezipf": int(sum(row["rank_bic_winner"] == "moezipf" for row in rows)),
        "sqrt_v_7param": int(sum(row["rank_bic_winner"] == "sqrt_v_7param" for row in rows)),
        "reranked_8param": int(sum(row["rank_bic_winner"] == "reranked_8param" for row in rows)),
    }

    summary = {
        "rows": sanitize(rows),
        "winner_counts": winner_counts,
        "moby_paper_check": sanitize(moby_check),
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(rows, winner_counts, moby_check), encoding="utf-8")
    write_csv(rows, OUTDIR / "moezipf_table.csv")


if __name__ == "__main__":
    main()
