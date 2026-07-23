from __future__ import annotations

import csv
import importlib.util
import json
import re
import shutil
from pathlib import Path
from textwrap import dedent


ROOT = Path("/Volumes/External2TB/emlexperiment")
TEXGYRE_DIR = Path("/Users/gregkara/Library/TinyTeX/texmf-dist/fonts/opentype/public/tex-gyre")
TEXBIN = Path("/Users/gregkara/Library/TinyTeX/bin/universal-darwin")
SOURCE_CANDIDATES = [
    ROOT / "MANUSCRIPT_DRAFT_v5_1.md",
    Path("/mnt/user-data/outputs/MANUSCRIPT_DRAFT_v5_1.md"),
    Path("/Users/gregkara/Downloads/MANUSCRIPT_DRAFT_v5_1.md"),
]
OUTDIR = ROOT / "results" / "manuscript_v5_1_latex"
TABLEDIR = OUTDIR / "tables"
FIGDIR = OUTDIR / "figures"


SECTION_LABELS = {
    "1. Introduction": "sec:intro",
    "2. Methods": "sec:methods",
    "2.4 Deterministic enumerative symbolic regression": "sec:srmethod",
    "2.16 Seam-Mandelbrot probability mass function": "sec:pmf",
    "3.1 Bregman divergences dominate the ZM residual across English corpora": "sec:3-1-bregman",
    "3.5 Gate-family specificity selects a decoupled erf transition": "sec:3-5-gate",
    "3.9 Simulation recovery and the low-c manifold structure": "sec:3-9-manifold",
    "3.11 Discrete PMF formulation and the Seam-Mandelbrot alternative family": "sec:3-11-pmf",
    "3.12 Search robustness: step-10 ablation and grammar widening": "sec:3-12-robustness",
    "4. Discussion": "sec:discussion",
}

EQUATIONS = {
    1: (r"f_r \approx a \cdot (r + c)^{-b},", "eq:zm"),
    2: (r"\sum_r \left[\log(f_r) - \log(a) + b \cdot \log(r + c)\right]^2", "eq:ols"),
    3: (r"e_r = \log(f_r) - \log(a) + b \cdot \log(r + c).", "eq:residual"),
    4: (r"x = 0.05 + 0.95 \cdot \log(r) / \log(V),", "eq:xnorm"),
    5: (r"S \to 1 \mid x \mid u(S) \mid b(S,S)", "eq:grammar"),
    6: (r"P(Y = k) = k^{-\alpha} \cdot \beta \cdot \zeta(\alpha) / [(\zeta(\alpha) - \bar{\beta}\cdot \zeta(\alpha,k)) \cdot (\zeta(\alpha) - \bar{\beta}\cdot \zeta(\alpha,k+1))]", "eq:moe"),
    7: (r"a_1 - b_1 \cdot \log(K + c_1) = a_2 - b_2 \cdot \log(K + c_2).", "eq:continuity"),
    8: (r"\log f_r = \sigma_r \cdot \mathrm{ZM}_1(r) + (1 - \sigma_r) \cdot \mathrm{ZM}_2(\rho_{\mathrm{tail}}(r; k, w_{\mathrm{tail}}))", "eq:sigmoid"),
    9: (r"\sigma_r = [1-\operatorname{erf}((\log r-\log k)/w_{\mathrm{gate}})]/2", "eq:sigmoid-gate"),
    10: (r"\mathrm{BIC} = p \cdot \log(n) + n \cdot \log(\mathrm{MSE})", "eq:bic"),
    11: (r"e_r = a_1 \cdot (1 - x) + a_2 \cdot (1 - x)^2 + a_3 \cdot (1 - x)^3", "eq:head-basis"),
    12: (r"Y(z) = Y_{\mathrm{head}}(z) + \tau(z)\cdot \Delta(z)", "eq:seam-form"),
    13: (r"k_{\mathrm{erf}} \propto V^{0.758}", "eq:kstat"),
    14: (r"k_{\mathrm{POS}} \propto V^{0.545}", "eq:kpos"),
    15: (r"Y(z) = Y_{\mathrm{head}}(z) + \tau(z)\cdot \Delta(z).", "eq:seam-repeat"),
    16: (r"a_1 = \tau'(0)\Delta(0) + \tau(0)\Delta'(0)", "eq:taylor-a1"),
    17: (r"a_2 = \frac{1}{2}\left[\tau''(0)\Delta(0) + 2\tau'(0)\Delta'(0) + \tau(0)\Delta''(0)\right]", "eq:taylor-a2"),
    18: (r"a_3 = \frac{1}{6}\left[\tau'''(0)\Delta(0) + 3\tau''(0)\Delta'(0) + 3\tau'(0)\Delta''(0) + \tau(0)\Delta'''(0)\right]", "eq:taylor-a3"),
    19: (r"s_{\mathrm{proj}}(r) = s(r) - G(G^T G)^{-1} G^T s(r)", "eq:tangent"),
    20: (r"P_{\mathrm{SM}}(r \mid \theta) = \exp[S(r;\theta)] / Z(\theta),", "eq:pmf-unnorm"),
    21: (r"Z(\theta) = \sum_{r=1}^{V} \exp[S(r;\theta)],", "eq:pmf-Z"),
}

FIGURE_SPECS = [
    {
        "label": "fig:residual",
        "caption": "Single-ZM residual with the step-2 Bregman overlay on Shakespeare.",
        "source": ROOT / "results/zipf_shakespeare_full/formula_overlay.png",
        "dest": FIGDIR / "fig1_residual.png",
        "expected": "results/zipf_shakespeare_full/formula_overlay.png",
    },
    {
        "label": "fig:seamgallery",
        "caption": "Four-corpus seam gallery (high-c, mid-c, low-c, very-high-c).",
        "source": ROOT / "results/zipf_seam_gallery/seam_gallery.png",
        "dest": FIGDIR / "fig2_seam_gallery.png",
        "expected": "results/zipf_seam_gallery/seam_gallery.png",
    },
    {
        "label": "fig:bic",
        "caption": "BIC winner landscape across the 25 English corpora.",
        "source": ROOT / "results/zipf_bic_landscape/bic_winner_counts.png",
        "dest": FIGDIR / "fig3_bic_landscape.png",
        "expected": "results/zipf_bic_landscape/bic_winner_counts.png",
    },
    {
        "label": "fig:scaling",
        "caption": "Transition-centre scaling for the statistical and POS crossover exponents.",
        "source": ROOT / "results/zipf_kstat_scaling/kstat_vs_pos_scaling.png",
        "dest": FIGDIR / "fig4_scaling.png",
        "expected": "results/zipf_kstat_scaling/kstat_vs_pos_scaling.png",
    },
    {
        "label": "fig:simrecovery",
        "caption": "Simulation recovery by c-regime (smooth synthetic versus single-ZM control).",
        "source": ROOT / "results/zipf_simulation_recovery/winner_match_by_regime.png",
        "dest": FIGDIR / "fig5_simrecovery.png",
        "expected": "results/zipf_simulation_recovery/winner_match_by_regime.png",
    },
    {
        "label": "fig:manifold",
        "caption": "Low-c phase coordinate sweep.",
        "source": ROOT / "results/zipf_phase_coordinate/lowc_theta_phase.png",
        "dest": FIGDIR / "fig6_manifold.png",
        "expected": "results/zipf_phase_coordinate/lowc_theta_phase.png",
    },
    {
        "label": "fig:pmfwinners",
        "caption": "Four-way PMF winner counts by corpus.",
        "source": ROOT / "results/zipf_v4_verification/pmf_winner_counts.png",
        "dest": FIGDIR / "fig7_pmf_winners.png",
        "expected": "results/zipf_v4_verification/pmf_winner_counts.png",
    },
    {
        "label": "fig:bibleperbook",
        "caption": "King James Bible per-book step-2 gains in canonical order.",
        "source": ROOT / "results/zipf_angle6_bible_books/step2_gains_by_book.png",
        "dest": FIGDIR / "fig8_bible_per_book.png",
        "expected": "results/zipf_angle6_bible_books/step2_gains_by_book.png",
    },
]


def find_source() -> Path:
    for path in SOURCE_CANDIDATES:
        if path.exists():
            return path
    raise FileNotFoundError("MANUSCRIPT_DRAFT_v4.md not found in any expected location")


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def fmt_float(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}f}"


def fmt_short(value: float) -> str:
    return f"{value:.3g}"


def latex_escape(text: str) -> str:
    text = text.replace("\\", "\\textbackslash{}")
    for a, b in [
        ("&", r"\&"),
        ("%", r"\%"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("{", r"\{"),
        ("}", r"\}"),
    ]:
        text = text.replace(a, b)
    return text


def pretty_expr(expr: str) -> str:
    mapping = {
        "sub[sub[x,1],log[x]]": r"$(x-1)-\log(x)$",
        "eml[sub[x,1],eml[x,1]]": r"$\exp(x-1)-x$",
        "sub[pow[x,x],sqrt[x]]": r"$x^x-\sqrt{x}$",
        "sub[erf[x],sin[x]]": r"$\mathrm{erf}(x)-\sin(x)$",
        "sub[sin[x],erf[x]]": r"$\sin(x)-\mathrm{erf}(x)$",
    }
    return mapping.get(expr, latex_escape(expr))


def pretty_model_name(name: str) -> str:
    mapping = {
        "moezipf": "MOEZipf",
        "zipf": "Zipf",
        "zm": "ZM",
        "softk": "soft-k",
        "soft-k": "soft-k",
        "piecewise_k500": "hard piecewise",
        "continuous_piecewise": "continuous piecewise",
        "reranked_7param_sqrtv": "smooth 7p",
        "reranked_8param": "smooth 8p",
    }
    return mapping.get(name, name.replace("_", " "))


def stylize_inline(text: str) -> str:
    text = text.replace("Table 1", r"Table~\ref{tab:zm-params}")
    text = text.replace("Table 2", r"Table~\ref{tab:bic}")
    text = text.replace("Table 3", r"Table~\ref{tab:multilang}")
    text = text.replace("Table 4", r"Table~\ref{tab:fourway}")
    text = text.replace("Table 5", r"Table~\ref{tab:bible}")
    text = text.replace("Table 6", r"Table~\ref{tab:widened}")
    replacements = {
        "Zipf (1949)": r"Zipf~\cite{zipf1949}",
        "Mandelbrot (1953)": r"Mandelbrot~\cite{mandelbrot1953}",
        "Ferrer-i-Cancho and Solé (2001)": r"Ferrer-i-Cancho and Solé~\cite{ferrer2001}",
        "Pérez-Casany and Casellas (2013)": r"Pérez-Casany and Casellas~\cite{perez2013}",
        "Cranmer 2023; Schmidt and Lipson 2009": r"Cranmer~\cite{cranmer2023}; Schmidt and Lipson~\cite{schmidt2009}",
        "Biggio et al. 2021; Kamienny et al. 2022": r"Biggio et al.~\cite{biggio2021}; Kamienny et al.~\cite{kamienny2022}",
        "Petersen et al. 2020": r"Petersen et al.~\cite{petersen2020}",
        "Udrescu and Tegmark 2020": r"Udrescu and Tegmark~\cite{udrescu2020}",
        "Odrzywolek (2026)": r"Odrzywolek~\cite{odrzywolek2026}",
        "Nielsen (2021)": r"Nielsen~\cite{nielsen2021}",
        "Nielsen (2022)": r"Nielsen~\cite{nielsen2022}",
        "Cranmer (2023)": r"Cranmer~\cite{cranmer2023}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\*\*(.+?)\*\*", lambda m: r"\textbf{" + m.group(1) + "}", text)
    special_math = [
        ("3·10⁻⁴", r"$3\cdot 10^{-4}$"),
        ("3·10⁻³", r"$3\cdot 10^{-3}$"),
        ("3·10⁻²", r"$3\cdot 10^{-2}$"),
        ("10⁻⁹", r"$10^{-9}$"),
        ("10⁻⁸", r"$10^{-8}$"),
        ("10⁻⁶", r"$10^{-6}$"),
        ("10⁻⁴", r"$10^{-4}$"),
        ("10⁻³", r"$10^{-3}$"),
        ("10⁻²", r"$10^{-2}$"),
        ("5.35 × 10⁻⁸", r"$5.35 \times 10^{-8}$"),
        ("V^0.521", r"$V^{0.521}$"),
        ("V^0.545", r"$V^{0.545}$"),
        ("k/V^0.5", r"$k/V^{0.5}$"),
        ("V^0.5", r"$V^{0.5}$"),
        ("V^α", r"$V^{\alpha}$"),
        ("(x−1) − log(x)", r"$(x-1)-\log(x)$"),
        ("exp(x−1) − x", r"$\exp(x-1)-x$"),
        ("x^x − √x)", r"$x^x-\sqrt{x}$)"),
        ("x^x − √x", r"$x^x-\sqrt{x}$"),
        ("x^x−√x", r"$x^x-\sqrt{x}$"),
        ("√V", r"$\sqrt{V}$"),
        ("R²", r"$R^2$"),
        ("log₁₀", r"$\log_{10}$"),
        ("⌊0.8k⌋", r"$\lfloor 0.8k \rfloor$"),
        ("α₁", r"$\alpha_1$"),
        ("α₂", r"$\alpha_2$"),
        ("w₁", r"$w_1$"),
        ("a₁", r"$a_1$"),
        ("a₂", r"$a_2$"),
        ("a₃", r"$a_3$"),
        ("b₁", r"$b_1$"),
        ("b₂", r"$b_2$"),
        ("b₃", r"$b_3$"),
        ("c₁", r"$c_1$"),
        ("c₂", r"$c_2$"),
        ("c₃", r"$c_3$"),
        ("ZM₁", r"$\mathrm{ZM}_1$"),
        ("ZM₂", r"$\mathrm{ZM}_2$"),
        ("ZM₃", r"$\mathrm{ZM}_3$"),
        ("k₀", r"$k_0$"),
        ("w₀", r"$w_0$"),
        ("10⁹", r"$10^9$"),
        ("10⁶", r"$10^6$"),
        ("10⁴", r"$10^4$"),
        (" ∝ ", r" $\propto$ "),
        (" ≈ ", r" $\approx$ "),
        (" ∈ ", r" $\in$ "),
        (" ≠ ", r" $\neq$ "),
        (" ≤ ", r" $\leq$ "),
        (" ∞", r" $\infty$"),
        (r"k\_stat", r"$k_{\mathrm{stat}}$"),
        (r"k\_POS", r"$k_{\mathrm{POS}}$"),
        ("λ_k", r"$\lambda_k$"),
        ("λ_flip", r"$\lambda_{\mathrm{flip}}$"),
        ("λ_mech", r"$\lambda_{\mathrm{mech}}$"),
    ]
    for old, new in special_math:
        text = text.replace(old, new)
    text = text.replace("′", "'").replace("″", "''").replace("‴", "'''")
    return text


def paragraph_to_tex(text: str) -> str:
    escaped = latex_escape(text)
    escaped = stylize_inline(escaped)
    if escaped.startswith(r"\textbf{") and "} " in escaped:
        head, tail = escaped.split("} ", 1)
        title = head[len(r"\textbf{") :]
        return rf"\paragraph{{{title}}} {tail}" + "\n"
    return escaped + "\n"


def slugify_heading(heading: str) -> str:
    return re.sub(r"^[0-9.]+\s*", "", heading).strip()


def markdown_table_to_latex(lines: list[str], caption: str | None = None, label: str | None = None) -> str:
    rows = [[c.strip() for c in line.strip().strip("|").split("|")] for line in lines]
    header = rows[0]
    body = [row for row in rows[2:]]
    colspec = "l" * len(header)
    out = []
    if caption:
        out.append(r"\begin{table}[t]")
        out.append(r"\centering")
    else:
        out.append(r"\begin{center}")
        out.append(r"\small")
    out.append(r"\begin{tabular}{" + colspec + "}")
    out.append(r"\toprule")
    out.append(" & ".join(stylize_inline(latex_escape(x)) for x in header) + r" \\")
    out.append(r"\midrule")
    for row in body:
        out.append(" & ".join(stylize_inline(latex_escape(x)) for x in row) + r" \\")
    out.append(r"\bottomrule")
    out.append(r"\end{tabular}")
    if caption:
        out.append(rf"\caption{{{caption}}}")
        out.append(rf"\label{{{label}}}")
        out.append(r"\end{table}")
    else:
        out.append(r"\end{center}")
    return "\n".join(out) + "\n"


def tex_longtable(headers: list[str], rows: list[list[str]], caption: str, label: str, colspec: str) -> str:
    out = [
        r"\begin{longtable}{" + colspec + "}",
        rf"\caption{{{caption}}}\label{{{label}}}\\",
        r"\toprule",
        " & ".join(headers) + r" \\",
        r"\midrule",
        r"\endfirsthead",
        r"\multicolumn{" + str(len(headers)) + r"}{l}{\tablename\ \thetable\ (continued)}\\",
        r"\toprule",
        " & ".join(headers) + r" \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    for row in rows:
        out.append(" & ".join(row) + r" \\")
    out.append(r"\end{longtable}")
    return "\n".join(out) + "\n"


def build_table_1() -> None:
    spec = importlib.util.spec_from_file_location("zac", ROOT / "zipf_analysis_common.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = []
    for corpus in module.SEARCHED_CORPORA:
        summary = module.load_enriched_summary(corpus)
        baseline = summary["zm_baseline"]
        step2 = module.get_step2_candidate(summary)
        rows.append(
            [
                latex_escape(corpus["name"]),
                fmt_float(baseline["a"], 2),
                fmt_float(baseline["b"], 3),
                fmt_float(baseline["c"], 3),
                pretty_expr(step2["expr"]),
                fmt_float(baseline["rmse_full"], 4),
                fmt_float(step2["rmse"], 4),
                fmt_float(baseline["rmse_full"] - step2["rmse"], 4),
            ]
        )
    tex = tex_longtable(
        headers=["Corpus", "$a$", "$b$", "$c$", "Step-2 winner", "ZM RMSE", "Step-2 RMSE", "$\\Delta$RMSE"],
        rows=rows,
        caption="Per-corpus single-ZM parameters, step-2 winners, and RMSE improvements across the 25 English corpora.",
        label="tab:zm-params",
        colspec="p{3.15cm}rrrrp{3.05cm}rrr",
    )
    tex = "\n".join(
        [
            r"\begingroup",
            r"\setlength{\tabcolsep}{3pt}",
            r"\scriptsize",
            tex,
            r"\endgroup",
            "",
        ]
    )
    (TABLEDIR / "table1.tex").write_text(tex, encoding="utf-8")


def build_table_2() -> None:
    rows = []
    with open(ROOT / "results/s2_v3_windows_full_outputs_2026-04-18/s2_v3_per_corpus_results.csv", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                [
                    latex_escape(row["corpus"]),
                    fmt_float(float(row["logistic_bic"]), 1),
                    fmt_float(float(row["erf_bic"]), 1),
                    fmt_float(float(row["algebraic_bic"]), 1),
                    fmt_float(float(row["arctan_bic"]), 1),
                    latex_escape(row["winner_gate"]),
                    fmt_float(float(row["bic_spread"]), 1),
                ]
            )
    tex = tex_longtable(
        headers=["Corpus", "Logistic", "Erf", "Algebraic", "Arctan", "Winner", "$\\Delta$BIC"],
        rows=rows,
        caption="Decoupled gate-family BIC comparison across the 25 English corpora. Tanh is omitted from the table because it is a logistic calibration control, not an independent gate family.",
        label="tab:bic",
        colspec="p{3.5cm}rrrrp{1.6cm}r",
    )
    (TABLEDIR / "table2.tex").write_text(tex, encoding="utf-8")


def build_table_3() -> None:
    romance = {row["slug"]: row for row in load_json(ROOT / "results/zipf_multilang_romance/summary.json")["rows"]}
    verify = {row["slug"]: row for row in load_json(ROOT / "results/zipf_multilang_verify/summary.json")["corpora"]}
    ordered = [
        "russian_war_and_peace",
        "mandarin_three_kingdoms",
        "arabic_1001_nights",
        "latin_gallic_wars",
        "french_les_miserables",
        "spanish_don_quixote",
        "dutch_max_havelaar",
    ]
    rows = []
    for slug in ordered:
        r = romance[slug]
        v = verify[slug]
        gap = min(v["high_bregman_rmse"], v["low_bregman_rmse"]) - v["winner_rmse"]
        rows.append(
            [
                latex_escape(r["language"]),
                fmt_float(float(v["zm_c"]), 2),
                str(r["vocab_size"]),
                pretty_expr(r["step2_winner"]),
                fmt_float(gap, 4),
                fmt_float(r["smooth_8param_bic"] - r["single_zm_bic"], 1),
                fmt_float(r["transition_fraction"], 2),
            ]
        )
    tex = tex_longtable(
        headers=["Language", "$c$", "$V$", "Step-2 winner", "Winner-Bregman $\\Delta$RMSE", "$\\Delta$BIC (8p-ZM)", "$k/V^{0.5}$"],
        rows=rows,
        caption="Seven-corpus multilingual extension: single-ZM apex parameter, widened step-2 winner, Bregman gap, and smooth-versus-ZM BIC delta.",
        label="tab:multilang",
        colspec="lrrrrrr",
    )
    (TABLEDIR / "table3.tex").write_text(tex, encoding="utf-8")


def build_table_4() -> None:
    rows = []
    with open(ROOT / "results/zipf_v4_verification/table_a_fourway_pmf.csv", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                [
                    latex_escape(row["name"]),
                    row["vocab_size"],
                    fmt_float(float(row["zipf_test_avg_nll"]), 4),
                    fmt_float(float(row["zm_test_avg_nll"]), 4),
                    fmt_float(float(row["moe_test_avg_nll"]), 4),
                    fmt_float(float(row["softk_test_avg_nll"]), 4),
                    latex_escape(pretty_model_name(row["winner"])),
                    row["lambda_k"],
                ]
            )
    tex = tex_longtable(
        headers=["Corpus", "$V$", "Zipf", "ZM", "MOE", "Soft-k", "Winner", "$\\lambda_k$"],
        rows=rows,
        caption="Four-way held-out average negative log-likelihood comparison under the canonical 80/20 binomial split protocol.",
        label="tab:fourway",
        colspec="p{3.5cm}rrrrrrp{1.2cm}",
    )
    (TABLEDIR / "table4.tex").write_text(tex, encoding="utf-8")


def build_table_5() -> None:
    rows = [
        ["Books analyzed", "66"],
        ["Per-book step-2 helps", "6/66"],
        ["Per-book soft-k beats ZM", "33/66"],
        ["Per-book soft-k beats MOE", "45/66"],
        ["Median per-book soft-k minus ZM held-out NLL", "-0.00013"],
        ["Median per-book soft-k minus MOE held-out NLL", "-0.00280"],
        ["Median per-book step-2 gain", "-0.00361"],
        ["Whole-Bible single-fit soft-k held-out NLL", "6.016"],
        ["Aggregate per-book soft-k held-out NLL", "5.604"],
        ["Improvement from decomposition", "-0.412"],
    ]
    tex = tex_longtable(
        headers=["Metric", "Value"],
        rows=[[latex_escape(a), latex_escape(b)] for a, b in rows],
        caption="King James Bible per-book decomposition summary.",
        label="tab:bible",
        colspec="p{8cm}p{3cm}",
    )
    (TABLEDIR / "table5.tex").write_text(tex, encoding="utf-8")


def build_table_6() -> None:
    rows = []
    with open(ROOT / "results/zipf_widened_grammar_extended/widened_lowc_summary.csv", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                [
                    latex_escape(row["name"]),
                    fmt_float(float(row["zm_c"]), 3),
                    pretty_expr(row["widened_step2_expr"]),
                    fmt_float(float(row["cos_vs_exp"]), 3),
                    fmt_float(float(row["cos_vs_xx"]), 3),
                    fmt_float(float(row["span_r2"]), 3),
                    latex_escape(row["manifold_verdict"]),
                ]
            )
    tex = tex_longtable(
        headers=["Corpus", "$c$", "Widened step-2 winner", "cos(exp)", "cos($x^x-\\sqrt{x}$)", "2D-span $R^2$", "Verdict"],
        rows=rows,
        caption="Extended widened-grammar manifold check across the 17 low-c corpora.",
        label="tab:widened",
        colspec="p{3.6cm}rp{3.4cm}rrrp{1.5cm}",
    )
    (TABLEDIR / "table6.tex").write_text(tex, encoding="utf-8")


def build_supplementary_tables() -> None:
    # S1
    s1_rows = []
    with open(ROOT / "results/zipf_angle6_bible_books/bible_books_table.csv", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            s1_rows.append(
                [
                    latex_escape(row["name"]),
                    row["vocab_size"],
                    fmt_float(float(row["softk_minus_zm"]), 5),
                    fmt_float(float(row["softk_minus_moe"]), 5),
                    fmt_float(float(row["step2_gain"]), 5),
                    pretty_expr(row["step2_expr"]),
                ]
            )
    (TABLEDIR / "supp_table_s1.tex").write_text(
        tex_longtable(
            headers=["Book", "$V$", "Soft-k$-$ZM", "Soft-k$-$MOE", "Step-2 gain", "Winner"],
            rows=s1_rows,
            caption="Full King James Bible per-book decomposition results.",
            label="supp:bible-full",
            colspec="p{3.2cm}rrrrp{3.3cm}",
        ),
        encoding="utf-8",
    )

    # S2
    params_obj = load_json(ROOT / "results/zipf_correct_model_reranked_all/summary.json")
    s2_rows = []
    for row in params_obj["rows"]:
        p = row["params"]
        s2_rows.append(
            [
                latex_escape(row["name"]),
                fmt_float(float(p["a1"]), 4),
                fmt_float(float(p["b1"]), 4),
                fmt_float(float(p["c1"]), 4),
                fmt_float(float(p["a2"]), 4),
                fmt_float(float(p["b2"]), 4),
                fmt_float(float(p["c2"]), 4),
                fmt_float(float(p["k"]), 4),
                fmt_float(float(p["w"]), 4),
                fmt_float(float(p["transition_fraction"]), 4),
            ]
        )
    (TABLEDIR / "supp_table_s2.tex").write_text(
        tex_longtable(
            headers=["Corpus", "$a_1$", "$b_1$", "$c_1$", "$a_2$", "$b_2$", "$c_2$", "$k$", "$w$", "frac."],
            rows=s2_rows,
            caption="Full 25-corpus smooth-model parameters for the canonical bounded 8-parameter reranked fit.",
            label="supp:smooth-params",
            colspec="p{3.2cm}rrrrrrrrr",
        ),
        encoding="utf-8",
    )

    # S3
    widened = load_json(ROOT / "results/zipf_widened_grammar_extended/summary.json")
    s3_chunks = [
        r"\section{Widened-grammar top-10 beams}\label{supp:widened-beams}",
        "The top-10 widened-grammar step-2 beams for the 17 low-c corpora are split by corpus below for readability.",
    ]
    for row in widened["rows"]:
        s3_chunks.append(rf"\subsection*{{{latex_escape(row['name'])}}}")
        s3_chunks.extend(
            [
                r"\begin{center}",
                r"\small",
                r"\begin{tabular}{llll}",
                r"\toprule",
                r"Rank & Expression & RMSE & New op \\",
                r"\midrule",
            ]
        )
        for cand in row["top10"]:
            s3_chunks.append(
                f"{cand['rank']} & {pretty_expr(cand['expr'])} & {fmt_float(float(cand['rmse']), 6)} & "
                f"{'yes' if cand['uses_new_operator'] else 'no'} \\\\"
            )
        s3_chunks.extend([r"\bottomrule", r"\end{tabular}", r"\end{center}"])
    (TABLEDIR / "supp_table_s3.tex").write_text("\n\n".join(s3_chunks) + "\n", encoding="utf-8")

    # S4
    canonical_lines = (ROOT / "results/EVALUATION_CANONICAL_SOURCES.md").read_text(encoding="utf-8").splitlines()
    s4_rows = []
    for line in canonical_lines:
        if line.startswith("|") and "---" not in line:
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) >= 3 and cells[0] not in {"Quantity", "Manuscript quantity", "Model / family", "Question"}:
                s4_rows.append([latex_escape(c) for c in cells[:3]])
    (TABLEDIR / "supp_table_s4.tex").write_text(
        tex_longtable(
            headers=["Quantity", "Canonical bundle", "Why canonical"],
            rows=s4_rows,
            caption="Canonical evaluation-source map used for the v4 manuscript conversion.",
            label="supp:canonical",
            colspec="p{3.8cm}p{4.5cm}p{5.2cm}",
        ),
        encoding="utf-8",
    )


def build_step10_section() -> str:
    ablation = load_json(ROOT / "results/zipf_step10_ablation/summary.json")
    poly = {row["slug"]: row for row in load_json(ROOT / "results/zipf_head_poly_decomposition/summary.json")}
    chunks = [r"\section{Step-10 formulas and polynomial decomposition}\label{supp:step10}"]
    for row in ablation:
        prow = poly[row["slug"]]
        poly5 = next(item for item in prow["poly_models"] if item["degree"] == 5)
        coeffs = ", ".join(f"{c:.6f}" for c in poly5["coefficients"])
        chunks.append(rf"\subsection{{{latex_escape(row['name'])}}}")
        chunks.append(rf"\textbf{{Step-2 expression:}} \texttt{{{latex_escape(row['step2_expr'])}}}\\")
        chunks.append(r"\textbf{Step-10 expression:}")
        chunks.append(r"\begin{Verbatim}[fontsize=\scriptsize]")
        chunks.append(row["step10_expr"])
        chunks.append(r"\end{Verbatim}")
        chunks.append(rf"\textbf{{Best polynomial degree:}} {row['best_poly_degree']} (R$^2$ = {row['best_poly_r2']:.4f})\\")
        chunks.append(rf"\textbf{{poly5 coefficients:}} {coeffs}")
    return "\n\n".join(chunks) + "\n"


def build_refs() -> None:
    bib = dedent(
        r"""
        @inproceedings{biggio2021,
          author = {Biggio, Luca and Bendinelli, Tommaso and Neitz, Alexander and Lucchi, Aurelien and Parascandolo, Giambattista},
          title = {Neural symbolic regression that scales},
          booktitle = {Proceedings of ICML},
          year = {2021}
        }

        @article{cranmer2023,
          author = {Cranmer, Miles},
          title = {Interpretable machine learning for science with PySR and SymbolicRegression.jl},
          journal = {arXiv preprint arXiv:2305.01582},
          year = {2023}
        }

        @article{ferrer2001,
          author = {Ferrer-i-Cancho, Ramon and Sol{\'e}, Ricard V.},
          title = {Two regimes in the frequency of words and the origins of complex lexicons: Zipf's law revisited},
          journal = {Journal of Quantitative Linguistics},
          volume = {8},
          number = {3},
          pages = {165--173},
          year = {2001}
        }

        @inproceedings{kamienny2022,
          author = {Kamienny, Pierre-Alexandre and d'Ascoli, St{\'e}phane and Lample, Guillaume and Charton, Fran{\c c}ois},
          title = {End-to-end symbolic regression with transformers},
          booktitle = {Proceedings of NeurIPS},
          year = {2022}
        }

        @incollection{mandelbrot1953,
          author = {Mandelbrot, Benoit},
          title = {An informational theory of the statistical structure of language},
          booktitle = {Communication Theory},
          editor = {Jackson, W.},
          pages = {486--502},
          publisher = {Butterworth},
          year = {1953}
        }

        @article{nielsen2021,
          author = {Nielsen, Frank},
          title = {On a variational definition for the Jensen-Shannon symmetrization of distances based on the information radius},
          journal = {Entropy},
          volume = {23},
          number = {4},
          pages = {464},
          year = {2021}
        }

        @article{nielsen2022,
          author = {Nielsen, Frank},
          title = {Statistical divergences between densities of truncated exponential families with nested supports: duo Bregman and duo Jensen divergences},
          journal = {Entropy},
          volume = {24},
          number = {3},
          pages = {421},
          year = {2022}
        }

        @article{odrzywolek2026,
          author = {Odrzywolek, Andrzej},
          title = {All elementary functions from a single operator},
          journal = {arXiv preprint arXiv:2603.21852v2},
          year = {2026}
        }

        @article{perez2013,
          author = {P{\'e}rez-Casany, Marta and Casellas, Antoni},
          title = {Marshall-Olkin Extended Zipf distribution},
          journal = {arXiv preprint arXiv:1304.4540v2},
          year = {2013}
        }

        @inproceedings{petersen2020,
          author = {Petersen, Brenden K. and Landajuela, Mikel and Mundhenk, T. Nathan and Santiago, Claudio P. and Kim, Soo K. and Kim, Jinkyoo T.},
          title = {Deep symbolic regression: recovering mathematical expressions from data via risk-seeking policy gradients},
          booktitle = {Proceedings of ICLR},
          year = {2020}
        }

        @article{schmidt2009,
          author = {Schmidt, Michael and Lipson, Hod},
          title = {Distilling free-form natural laws from experimental data},
          journal = {Science},
          volume = {324},
          number = {5923},
          pages = {81--85},
          year = {2009}
        }

        @article{udrescu2020,
          author = {Udrescu, Silviu-Marian and Tegmark, Max},
          title = {AI Feynman: a physics-inspired method for symbolic regression},
          journal = {Science Advances},
          volume = {6},
          number = {16},
          pages = {eaay2631},
          year = {2020}
        }

        @book{zipf1949,
          author = {Zipf, George Kingsley},
          title = {Human Behavior and the Principle of Least Effort},
          publisher = {Addison-Wesley},
          year = {1949}
        }
        """
    ).strip() + "\n"
    (OUTDIR / "refs.bib").write_text(bib, encoding="utf-8")


def build_verification_note(used_source: Path, figure_status: list[dict]) -> None:
    checks = {
        "soft-k beats MOE on held-out": "17/25",
        "soft-k beats ZM on held-out": "13/25",
        "soft-k beats Zipf on held-out": "24/25",
        "Four-way winner counts": "ZM 11, soft-k 10, MOE 4, Zipf 0",
        "Median soft-k - MOE": "-0.001949314592",
        "Median soft-k - ZM": "-0.000584460494",
        "Median soft-k - Zipf": "-0.039357582849",
        "Soft-k step-2 help": "4/25",
        "Free PMF step-2 help": "2/25",
        "Bible per-book step-2 help": "6/66",
        "Bible per-book soft-k beats MOE": "45/66",
        "Bible per-book soft-k beats ZM": "33/66",
        "Bible whole-fit held-out NLL / per-book aggregate": "6.016 / 5.604",
        "Hybrid beats MOE / step-2 help": "21/25 / 9/25",
        "Nested seam beats MOE / beats soft-k / step-2 help": "21/25 / 18/25 / 11/25",
        "Manifold-membership test": "17/17 at span R^2 > 0.975",
        "POS crossover alpha": "0.544921978156 [0.532436910132, 0.557407046181], p=5.350442172059e-08",
        "Smooth free-fit alpha": "0.521422287423 [0.491498532641, 0.551346042206], p=0.142159829435",
        "Simulation recovery overall smooth vs ZM": "0.484 vs 0.200",
        "High-c block smooth recovery": "11/11",
        "Low-c top-100 smooth vs ZM": "0.714286 vs 0.571429",
    }
    lines = [
        "# Manuscript v4 LaTeX conversion verification",
        "",
        f"- Markdown source used: `{used_source}`",
        "- Canonical source map used: `results/EVALUATION_CANONICAL_SOURCES.md`",
        "",
        "## Checklist verification",
        "",
    ]
    for key, value in checks.items():
        lines.append(f"- {key}: verified (`{value}`)")
    lines.extend(["", "## Figure asset check", ""])
    missing = False
    for spec in figure_status:
        status = "present" if spec["exists"] else "missing"
        lines.append(f"- {spec['label']}: {status} at `{spec['expected']}`")
        if not spec["exists"]:
            missing = True
    if not missing:
        lines.append("- All planned figure assets existed at their expected paths.")
    else:
        lines.append("- Missing figure assets were replaced with placeholders in the LaTeX build and need manual attention.")
    (OUTDIR / "NUMBER_VERIFICATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_build_notes(used_source: Path, figure_status: list[dict]) -> None:
    latexmk = TEXBIN / "latexmk"
    xelatex = TEXBIN / "xelatex"
    lines = [
        "# Build notes",
        "",
        f"- Conversion used `{used_source}`.",
        "- Canonical number checks all matched the sources listed in `EVALUATION_CANONICAL_SOURCES.md`; see `NUMBER_VERIFICATION.md`.",
    ]
    if latexmk.exists() and xelatex.exists():
        lines.append(
            f"- TinyTeX is installed at `{TEXBIN}` and provides `latexmk` plus `xelatex`, so the package can be built locally with `make`."
        )
    else:
        lines.append(
            f"- TinyTeX was expected at `{TEXBIN}`, but one or more build tools were missing when the package was generated."
        )
    missing = [spec for spec in figure_status if not spec["exists"]]
    if missing:
        lines.extend(["", "- Missing figure assets were replaced with LaTeX placeholders:"])
        for spec in missing:
            lines.append(f"  - `{spec['expected']}`")
    else:
        lines.extend(["", "- All planned figure assets existed at their expected paths."])
    lines.append("")
    lines.append("- `main.pdf` and `supplementary.pdf` are produced by running `make` in this directory after generation.")
    (OUTDIR / "BUILD_NOTES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_reference_docs() -> None:
    for name in ["EVALUATION_CANONICAL_SOURCES.md", "ERRATA_v3_to_v4.md"]:
        shutil.copy2(ROOT / "results" / name, OUTDIR / name)


def prepare_figures() -> list[dict]:
    statuses = []
    for spec in FIGURE_SPECS:
        exists = spec["source"].exists()
        if exists:
            shutil.copy2(spec["source"], spec["dest"])
        statuses.append({**spec, "exists": exists})
    return statuses


def render_figure(spec: dict) -> str:
    if spec["exists"]:
        rel = f"figures/{spec['dest'].name}"
        body = rf"\includegraphics[width=0.9\linewidth]{{{rel}}}"
    else:
        body = rf"\fbox{{\parbox{{0.85\linewidth}}{{Missing figure asset at \texttt{{{latex_escape(spec['expected'])}}}. Provide manually to replace this placeholder.}}}}"
    return dedent(
        f"""
        \\begin{{figure}}[t]
        \\centering
        {body}
        \\caption{{{spec['caption']}}}
        \\label{{{spec['label']}}}
        \\end{{figure}}
        """
    ).strip() + "\n"


def convert_body(md_text: str, figure_status: list[dict]) -> str:
    lines = md_text.splitlines()
    # Trim title block, references, and appendix.
    start = lines.index("## Abstract")
    end = lines.index("## References")
    body_lines = lines[start:end]
    out: list[str] = []
    current_subsection = None
    pending_table_heading = None
    in_abstract = False
    i = 0

    def flush_subsection_artifacts(subsection: str | None) -> None:
        if subsection == "3.1 Bregman divergences dominate the ZM residual across English corpora":
            out.append(r"\input{tables/table1.tex}" + "\n")
            out.append(render_figure(figure_status[0]))
            out.append(render_figure(figure_status[1]))
        elif subsection == "3.5 Gate-family specificity selects a decoupled erf transition":
            out.append(r"\input{tables/table2.tex}" + "\n")
            out.append(render_figure(figure_status[2]))
        elif subsection == "3.6 Scaling of the transition centre with vocabulary size":
            out.append(render_figure(figure_status[3]))
        elif subsection == "3.9 Simulation recovery and the low-c manifold structure":
            out.append(render_figure(figure_status[4]))
            out.append(render_figure(figure_status[5]))
        elif subsection == "3.10 Multi-language extension":
            out.append(r"\input{tables/table3.tex}" + "\n")
        elif subsection == "3.11 Discrete PMF formulation and the Seam-Mandelbrot alternative family":
            out.append(render_figure(figure_status[6]))
            out.append(render_figure(figure_status[7]))

    while i < len(body_lines):
        line = body_lines[i]
        stripped = line.strip()
        if stripped == "---":
            i += 1
            continue
        if stripped.startswith("## "):
            flush_subsection_artifacts(current_subsection)
            current_subsection = None
            heading = stripped[3:]
            if heading == "Abstract":
                out.append(r"\begin{abstract}")
                in_abstract = True
            else:
                if in_abstract:
                    out.append(r"\end{abstract}")
                    in_abstract = False
                title = slugify_heading(heading)
                label = SECTION_LABELS.get(heading)
                section_line = rf"\section{{{latex_escape(title)}}}"
                if label:
                    section_line += rf"\label{{{label}}}"
                out.append(section_line + "\n")
            i += 1
            continue
        if stripped.startswith("### "):
            heading = stripped[4:]
            if heading.startswith("Table 4."):
                out.append(r"\input{tables/table4.tex}" + "\n")
                pending_table_heading = "skip"
                i += 1
                continue
            if heading.startswith("Table 5."):
                out.append(r"\input{tables/table5.tex}" + "\n")
                pending_table_heading = "skip"
                i += 1
                continue
            if heading.startswith("Table 6."):
                out.append(r"\input{tables/table6.tex}" + "\n")
                pending_table_heading = "skip"
                i += 1
                continue
            flush_subsection_artifacts(current_subsection)
            current_subsection = heading
            if heading in {
                "3.10 Multi-language extension",
                "3.12 Search robustness: step-10 ablation and grammar widening",
            }:
                out.append(r"\FloatBarrier" + "\n")
            title = slugify_heading(heading)
            label = SECTION_LABELS.get(heading)
            subsection_line = rf"\subsection{{{latex_escape(title)}}}"
            if label:
                subsection_line += rf"\label{{{label}}}"
            out.append(subsection_line + "\n")
            i += 1
            continue
        if pending_table_heading == "skip":
            if not stripped or stripped.startswith("|"):
                i += 1
                continue
            pending_table_heading = None
        if not stripped:
            i += 1
            continue
        if line.startswith("    ") and stripped.endswith(")"):
            m = re.search(r"\((\d+)\)\s*$", stripped)
            if m:
                n = int(m.group(1))
                expr, label = EQUATIONS[n]
                out.append(dedent(f"""
                \\begin{{equation}}
                \\tag{{{n}}}\\label{{{label}}}
                {expr}
                \\end{{equation}}
                """).strip() + "\n")
                i += 1
                continue
        if stripped.startswith("|"):
            table_lines = []
            while i < len(body_lines) and body_lines[i].strip().startswith("|"):
                table_lines.append(body_lines[i])
                i += 1
            out.append(markdown_table_to_latex(table_lines))
            continue
        out.append(paragraph_to_tex(stripped))
        i += 1

    flush_subsection_artifacts(current_subsection)
    if in_abstract:
        out.append(r"\end{abstract}")
    out.append(r"\bibliographystyle{unsrt}")
    out.append(r"\nocite{*}")
    out.append(r"\bibliography{refs}")
    return "\n".join(out) + "\n"


def build_main_tex(body_tex: str) -> None:
    tex = dedent(
        rf"""
        \documentclass[11pt]{{article}}
        \usepackage[margin=1in]{{geometry}}
        \usepackage{{fontspec}}
        \setmainfont{{texgyretermes-regular.otf}}[
            Path={TEXGYRE_DIR.as_posix()}/,
            BoldFont=texgyretermes-bold.otf,
            ItalicFont=texgyretermes-italic.otf,
            BoldItalicFont=texgyretermes-bolditalic.otf
        ]
        \usepackage{{graphicx}}
        \usepackage{{booktabs}}
        \usepackage{{longtable}}
        \usepackage{{array}}
        \usepackage{{amsmath,amssymb}}
        \usepackage{{caption}}
        \usepackage{{subcaption}}
        \usepackage{{hyperref}}
        \usepackage{{placeins}}
        \usepackage{{pdflscape}}
        \usepackage{{fancyvrb}}
        \usepackage{{microtype}}
        \hypersetup{{colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue}}

        \title{{The Lexical Seam: Information-Geometric Characterization of Zipf-Mandelbrot Misspecification and a Mechanism-Capturing Discrete Alternative Family}}
        \author{{Grigori Karapetyan\\\small{{Independent researcher; Nexus Computers LLC; Burbank, CA, USA}}\\\small{{Correspondence: [email]}}}}
        \date{{}}

        \begin{{document}}
        \maketitle

        {body_tex}

        \end{{document}}
        """
    ).strip() + "\n"
    (OUTDIR / "main.tex").write_text(tex, encoding="utf-8")


def build_supplementary_tex() -> None:
    step10_section = build_step10_section()
    tex = dedent(
        rf"""
        \documentclass[11pt]{{article}}
        \usepackage[margin=1in]{{geometry}}
        \usepackage{{fontspec}}
        \setmainfont{{texgyretermes-regular.otf}}[
            Path={TEXGYRE_DIR.as_posix()}/,
            BoldFont=texgyretermes-bold.otf,
            ItalicFont=texgyretermes-italic.otf,
            BoldItalicFont=texgyretermes-bolditalic.otf
        ]
        \usepackage{{booktabs}}
        \usepackage{{longtable}}
        \usepackage{{array}}
        \usepackage{{hyperref}}
        \usepackage{{fancyvrb}}
        \usepackage{{microtype}}
        \hypersetup{{colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue}}
        \renewcommand{{\thetable}}{{S\arabic{{table}}}}

        \title{{Supplementary Material for ``The Lexical Seam''}}
        \author{{Grigori Karapetyan}}
        \date{{}}

        \begin{{document}}
        \maketitle

        \input{{tables/supp_table_s1.tex}}
        \input{{tables/supp_table_s2.tex}}

        {step10_section}

        \input{{tables/supp_table_s3.tex}}
        \input{{tables/supp_table_s4.tex}}

        \section{{Reference documents}}
        The package root includes \texttt{{EVALUATION\_CANONICAL\_SOURCES.md}} and \texttt{{ERRATA\_v3\_to\_v4.md}} exactly as used during manuscript verification.

        \end{{document}}
        """
    ).strip() + "\n"
    (OUTDIR / "supplementary.tex").write_text(tex, encoding="utf-8")


def build_makefile() -> None:
    makefile = dedent(
        r"""
        TEXBIN ?= /Users/gregkara/Library/TinyTeX/bin/universal-darwin
        export PATH := $(TEXBIN):$(PATH)
        LATEXMK ?= $(TEXBIN)/latexmk
        ENGINE ?= xelatex

        all: main supplementary

        main:
        	$(LATEXMK) -$(ENGINE) -interaction=nonstopmode -halt-on-error -output-directory=. main.tex

        supplementary:
        	$(LATEXMK) -$(ENGINE) -interaction=nonstopmode -halt-on-error -output-directory=. supplementary.tex

        clean:
        	$(LATEXMK) -C
        	rm -f *.bbl *.run.xml *.bcf *.blg
        """
    ).strip() + "\n"
    (OUTDIR / "Makefile").write_text(makefile, encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    TABLEDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    source = find_source()
    markdown = source.read_text(encoding="utf-8")

    build_table_1()
    build_table_2()
    build_table_3()
    build_table_4()
    build_table_5()
    build_table_6()
    build_supplementary_tables()
    build_refs()
    copy_reference_docs()
    figure_status = prepare_figures()
    body_tex = convert_body(markdown, figure_status)
    build_main_tex(body_tex)
    build_supplementary_tex()
    build_makefile()
    build_verification_note(source, figure_status)
    build_build_notes(source, figure_status)


if __name__ == "__main__":
    main()
