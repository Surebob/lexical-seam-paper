from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm, t, ttest_1samp


ROOT = Path("/Volumes/External2TB/emlexperiment")
RERANKED_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_correct_model_reranked_all" / "summary.json"
POS_ALL_SUMMARY_PATH = ROOT / "results" / "zipf_pos_all_corpora" / "summary.json"
OUTDIR = ROOT / "results" / "zipf_kstat_scaling"


def forced_scaling_model(x, alpha):
    return alpha * x


def load_rows():
    stat_summary = json.loads(RERANKED_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))
    pos_summary = json.loads(POS_ALL_SUMMARY_PATH.read_text(encoding="utf-8"))

    pos_map = {row["slug"]: row for row in pos_summary["rows"]}
    rows = []
    for row in stat_summary["rows"]:
        slug = row["slug"]
        V = float(row["vocab_size"])
        k_stat = float(row["k"])
        pos_row = pos_map[slug]
        rows.append(
            {
                "slug": slug,
                "name": row["name"],
                "V": V,
                "k_stat": k_stat,
                "log_V": float(math.log(V)),
                "log_k_stat": float(math.log(k_stat)),
                "alpha_stat_per_corpus": float(math.log(k_stat) / math.log(V)),
                "k_pos": float(pos_row["k_crossover"]),
                "log_k_pos": float(pos_row["log_k"]),
                "alpha_pos_per_corpus": float(pos_row["alpha_per_corpus"]),
            }
        )
    return rows, pos_summary["forced_fit"]


def fit_forced_alpha(rows: list[dict], key: str) -> dict:
    x = np.array([row["log_V"] for row in rows], dtype=float)
    y = np.array([row[key] for row in rows], dtype=float)
    popt, pcov = curve_fit(forced_scaling_model, x, y, p0=(0.5,), maxfev=10000)
    alpha = float(popt[0])
    se = float(math.sqrt(pcov[0, 0]))
    df = len(rows) - 1
    tcrit = float(t.ppf(0.975, df))
    return {
        "alpha": alpha,
        "alpha_se": se,
        "alpha_ci_95": [alpha - tcrit * se, alpha + tcrit * se],
        "df": df,
    }


def summarize_alpha_distribution(rows: list[dict], key: str) -> dict:
    values = np.array([row[key] for row in rows], dtype=float)
    mean_value = float(np.mean(values))
    median_value = float(np.median(values))
    se = float(np.std(values, ddof=1) / math.sqrt(len(values)))
    tcrit = float(t.ppf(0.975, len(values) - 1))
    mean_ci = [mean_value - tcrit * se, mean_value + tcrit * se]
    ttest_result = ttest_1samp(values, popmean=0.5)
    return {
        "mean_alpha": mean_value,
        "median_alpha": median_value,
        "mean_alpha_ci_95": mean_ci,
        "t_statistic_vs_0_5": float(ttest_result.statistic),
        "p_value_vs_0_5": float(ttest_result.pvalue),
        "df": int(len(values) - 1),
    }


def compare_alpha_fits(stat_fit: dict, pos_fit: dict, rows: list[dict]) -> dict:
    stat_alpha = float(stat_fit["alpha"])
    pos_alpha = float(pos_fit["alpha"])
    diff = stat_alpha - pos_alpha
    se_diff_indep = math.sqrt(float(stat_fit["alpha_se"]) ** 2 + float(pos_fit["alpha_se"]) ** 2)
    z = diff / se_diff_indep if se_diff_indep > 0 else float("nan")
    p = float(2.0 * norm.sf(abs(z))) if math.isfinite(z) else float("nan")

    per_corpus_diff = np.array(
        [row["alpha_stat_per_corpus"] - row["alpha_pos_per_corpus"] for row in rows],
        dtype=float,
    )
    mean_diff = float(np.mean(per_corpus_diff))
    median_diff = float(np.median(per_corpus_diff))
    se_paired = float(np.std(per_corpus_diff, ddof=1) / math.sqrt(len(per_corpus_diff)))
    tcrit = float(t.ppf(0.975, len(per_corpus_diff) - 1))
    mean_diff_ci = [mean_diff - tcrit * se_paired, mean_diff + tcrit * se_paired]
    ttest_paired = ttest_1samp(per_corpus_diff, popmean=0.0)

    ci_overlap = not (
        stat_fit["alpha_ci_95"][1] < pos_fit["alpha_ci_95"][0]
        or pos_fit["alpha_ci_95"][1] < stat_fit["alpha_ci_95"][0]
    )

    return {
        "alpha_difference": float(diff),
        "ci_overlap": bool(ci_overlap),
        "independent_fit_z": float(z),
        "independent_fit_p_value": p,
        "paired_mean_difference": mean_diff,
        "paired_median_difference": median_diff,
        "paired_mean_difference_ci_95": mean_diff_ci,
        "paired_t_statistic": float(ttest_paired.statistic),
        "paired_p_value": float(ttest_paired.pvalue),
    }


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "slug",
                "name",
                "V",
                "k_stat",
                "alpha_stat_per_corpus",
                "k_pos",
                "alpha_pos_per_corpus",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_scaling(rows: list[dict], stat_fit: dict, pos_fit: dict) -> None:
    x = np.array([row["log_V"] for row in rows], dtype=float)
    y_stat = np.array([row["log_k_stat"] for row in rows], dtype=float)
    y_pos = np.array([row["log_k_pos"] for row in rows], dtype=float)
    xline = np.linspace(float(x.min()) - 0.05, float(x.max()) + 0.05, 200)

    plt.figure(figsize=(8.8, 6.2))
    plt.scatter(x, y_stat, color="#174a7e", s=46, label="k_stat free-fit")
    plt.scatter(x, y_pos, color="#c23b22", s=34, alpha=0.7, label="k_POS crossover")
    plt.plot(xline, forced_scaling_model(xline, stat_fit["alpha"]), color="#174a7e", linewidth=2.2, label=f"k_stat fit: alpha={stat_fit['alpha']:.3f}")
    plt.plot(xline, forced_scaling_model(xline, pos_fit["alpha"]), color="#c23b22", linewidth=2.0, linestyle="--", label=f"k_POS fit: alpha={pos_fit['alpha']:.3f}")
    plt.plot(xline, forced_scaling_model(xline, 0.5), color="#6a4c93", linewidth=1.8, linestyle=":", label="sqrt(V) reference")
    plt.xlabel("log(V)")
    plt.ylabel("log(k)")
    plt.title("Statistical Seam k vs POS Crossover k")
    plt.legend(loc="best", fontsize=9)
    plt.tight_layout()
    plt.savefig(OUTDIR / "kstat_vs_pos_scaling.png", dpi=220)
    plt.close()


def build_report(stat_fit: dict, stat_dist: dict, pos_fit: dict, compare: dict) -> str:
    lines = [
        "# k_stat Scaling Fit",
        "",
        "## Forced Fit: k_stat = V^alpha",
        "",
        f"- Forced alpha: `{stat_fit['alpha']:.12f}`",
        f"- Forced alpha 95% CI: `[{stat_fit['alpha_ci_95'][0]:.12f}, {stat_fit['alpha_ci_95'][1]:.12f}]`",
        f"- CI {'includes' if stat_fit['alpha_ci_95'][0] <= 0.5 <= stat_fit['alpha_ci_95'][1] else 'excludes'} `0.50`.",
        "",
        "## Per-Corpus alpha_stat Distribution",
        "",
        f"- Mean alpha_stat: `{stat_dist['mean_alpha']:.12f}`",
        f"- Median alpha_stat: `{stat_dist['median_alpha']:.12f}`",
        f"- Mean-alpha 95% CI: `[{stat_dist['mean_alpha_ci_95'][0]:.12f}, {stat_dist['mean_alpha_ci_95'][1]:.12f}]`",
        f"- One-sample t-test vs `0.50`: `t={stat_dist['t_statistic_vs_0_5']:.12f}`, `p={stat_dist['p_value_vs_0_5']:.12f}`",
        "",
        "## Comparison To POS Exponent",
        "",
        f"- POS forced alpha reference: `{pos_fit['alpha']:.12f}`",
        f"- POS forced alpha 95% CI: `[{pos_fit['alpha_ci_95'][0]:.12f}, {pos_fit['alpha_ci_95'][1]:.12f}]`",
        f"- k_stat minus k_POS forced alpha: `{compare['alpha_difference']:.12f}`",
        f"- Forced-fit 95% CIs {'overlap' if compare['ci_overlap'] else 'do not overlap'}.",
        f"- Approx independent-fit z-test: `z={compare['independent_fit_z']:.12f}`, `p={compare['independent_fit_p_value']:.12f}`",
        f"- Paired mean(alpha_stat - alpha_pos): `{compare['paired_mean_difference']:.12f}`",
        f"- Paired median(alpha_stat - alpha_pos): `{compare['paired_median_difference']:.12f}`",
        f"- Paired mean-difference 95% CI: `[{compare['paired_mean_difference_ci_95'][0]:.12f}, {compare['paired_mean_difference_ci_95'][1]:.12f}]`",
        f"- Paired t-test vs zero: `t={compare['paired_t_statistic']:.12f}`, `p={compare['paired_p_value']:.12f}`",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, pos_fit = load_rows()
    rows = sorted(rows, key=lambda row: row["slug"])
    stat_fit = fit_forced_alpha(rows, "log_k_stat")
    stat_dist = summarize_alpha_distribution(rows, "alpha_stat_per_corpus")
    compare = compare_alpha_fits(stat_fit, pos_fit, rows)

    write_csv(rows, OUTDIR / "kstat_scaling_points.csv")
    plot_scaling(rows, stat_fit, pos_fit)

    summary = {
        "rows": rows,
        "kstat_fit": stat_fit,
        "kstat_distribution": stat_dist,
        "pos_fit_reference": pos_fit,
        "comparison": compare,
    }
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTDIR / "report.md").write_text(build_report(stat_fit, stat_dist, pos_fit, compare), encoding="utf-8")


if __name__ == "__main__":
    main()
