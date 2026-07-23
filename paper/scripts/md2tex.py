"""Convert MANUSCRIPT_v6.md to LaTeX (article class) with figures. Tuned to the
exact markdown subset used in v6; not a general converter."""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
MD = PAPER / "MANUSCRIPT_v6.md"
TEX = PAPER / "MANUSCRIPT_v6.tex"

PHRASES = [
    ("f(r) ∝ (r+c)^(−b)", r"$f(r) \propto (r+c)^{-b}$"),
    ("e^(x−1) = e^(−0.95)·r^(0.95/ln V)", r"$e^{x-1} = e^{-0.95}\,r^{0.95/\ln V}$"),
    ("e^(−0.95)·r^(0.95/ln V)", r"$e^{-0.95}\,r^{0.95/\ln V}$"),
    ("e^(x−1) − x", r"$e^{x-1}-x$"),
    ("e^(x−1)−x", r"$e^{x-1}-x$"),
    ("e^(x−1)", r"$e^{x-1}$"),
    ("(x−1) − log x", r"$(x-1)-\log x$"),
    ("(x−1)−log x", r"$(x-1)-\log x$"),
    ("x^x − √x", r"$x^{x}-\sqrt{x}$"),
    ("x^x−√x", r"$x^{x}-\sqrt{x}$"),
    ("k ≈ √V", r"$k \approx \sqrt{V}$"),
    ("√V", r"$\sqrt{V}$"),
    ("√e", r"$\sqrt{e}$"),
    ("V^0.545", r"$V^{0.545}$"),
    ("V^0.977", r"$V^{0.977}$"),
    ("V^1.003", r"$V^{1.003}$"),
    ("10⁻⁴", r"$10^{-4}$"),
    ("R²", r"$R^2$"),
    ("β ≈ 1/b", r"$\beta \approx 1/b$"),
    ("corr(β, 1/b)", r"$\mathrm{corr}(\beta, 1/b)$"),
    ("f(r) ≈ (r+c)^(−b)", r"$f(r)\approx(r+c)^{-b}$"),
    ("(a, b, c)", r"$(a, b, c)$"),
    ("[a-z]+(?:'[a-z]+)?", r"\texttt{[a-z]+(?:'[a-z]+)?}"),
    ("x = 0.05 + 0.95·log r / log V", r"$x = 0.05 + 0.95\,\log r/\log V$"),
    ("a₁ − b₁log(r+c₁)", r"$a_1 - b_1\log(r+c_1)$"),
    ("a₂ − b₂log(ρ+c₂)", r"$a_2 - b_2\log(\rho+c_2)$"),
    ("ρ(r; k, w_tail)", r"$\rho(r;k,w_{tail})$"),
    ("s = max(1, k·w_tail)", r"$s=\max(1, k\,w_{tail})$"),
    ("σ((log r − log k)/w_gate)", r"$\sigma((\log r-\log k)/w_{gate})$"),
    ("(π_H, μ_H, σ_H, μ_T, σ_T)", r"$(\pi_H,\mu_H,\sigma_H,\mu_T,\sigma_T)$"),
    ("V/(1−P₀)", r"$V/(1-P_0)$"),
    ("n ~ Poisson(λ)", r"$n \sim \mathrm{Poisson}(\lambda)$"),
    ("s = k·w_tail", r"$s = k\,w_{tail}$"),
    ("s ≈ 0.012·V", r"$s \approx 0.012\,V$"),
    ("s ∝ V^1.003", r"$s \propto V^{1.003}$"),
    ("w_tail", r"$w_{tail}$"),
    ("w_gate", r"$w_{gate}$"),
    ("k_POS", r"$k_{POS}$"),
    ("π_H", r"$\pi_H$"),
    ("σ_H", r"$\sigma_H$"),
    ("σ_T", r"$\sigma_T$"),
    ("μ_H", r"$\mu_H$"),
    ("μ_T", r"$\mu_T$"),
    ("P₀", r"$P_0$"),
    ("β = 0.758", r"$\beta = 0.758$"),
    ("β = 0.298", r"$\beta = 0.298$"),
    ("β = 0.476", r"$\beta = 0.476$"),
    ("λ-ZM", r"$\lambda$-ZM"),
]

UNI = {
    "≈": r"$\approx$", "≤": r"$\le$", "≥": r"$\ge$", "±": r"$\pm$", "×": r"$\times$",
    "→": r"$\to$", "∝": r"$\propto$", "·": r"$\cdot$", "−": "--UMINUS--",
    "Δ": r"$\Delta$", "σ": r"$\sigma$", "λ": r"$\lambda$", "β": r"$\beta$",
    "μ": r"$\mu$", "π": r"$\pi$", "ρ": r"$\rho$", "φ": r"$\varphi$",
    "α": r"$\alpha$", "ζ": r"$\zeta$", "τ": r"$\tau$", "Φ": r"$\Phi$",
    "¹": r"$^{1}$", "²": r"$^{2}$", "³": r"$^{3}$", "⁴": r"$^{4}$", "⁶": r"$^{6}$",
    "₀": r"$_{0}$", "₁": r"$_{1}$", "₂": r"$_{2}$",
    "§": r"\S", "—": "---", "–": "--", "’": "'", "‘": "`", "“": "``", "”": "''",
    "…": r"\dots", "≠": r"$\ne$", "⅙": r"$1/6$", "∈": r"$\in$", "∼": r"$\sim$",
    "∞": r"$\infty$",
}


def esc(t: str) -> str:
    t = t.replace("\\", r"\textbackslash{}")
    for a, b in [("{", r"\{"), ("}", r"\}"), ("%", r"\%"), ("&", r"\&"),
                 ("#", r"\#"), ("_", r"\_"), ("$", r"\$"), ("^", r"\^{}"),
                 ("~", r"\~{}")]:
        t = t.replace(a, b)
    return t


def prose(t: str) -> str:
    toks = []

    def stash(s):
        toks.append(s)
        return f"@@T{len(toks)-1}@@"

    t = re.sub(r"`([^`]+)`", lambda m: stash(r"\texttt{" + esc(m.group(1)) + "}"), t)
    for a, b in PHRASES:
        if a in t:
            t = t.replace(a, stash(b))
    for a, b in UNI.items():
        if a in t:
            # stash every replacement so esc() cannot re-escape the LaTeX it contains
            # (un-stashed "\S" became "\textbackslash{}S" in the v6.0 build)
            t = t.replace(a, stash(b))
    # straight-quote pairs -> LaTeX quotes (bare " typesets as a closing quote)
    t = re.sub(r'"([^"\n]+)"', r"``\1''", t)
    t = esc(t)
    t = t.replace("--UMINUS--", "-")
    t = re.sub(r"\*\*([^*]+)\*\*", r"\\textbf{\1}", t)
    t = re.sub(r"\*([^*]+)\*", r"\\emph{\1}", t)
    for i, s in enumerate(toks):
        t = t.replace(f"@@T{i}@@", s.replace("--UMINUS--", "-"))
    return t


def mathify(t: str) -> str:
    t = re.sub(r"\s*\(\d\)\s*$", "", t.strip()).rstrip(",")
    reps = [("log", r"\log "), ("≈", r"\approx "), ("∝", r"\propto "), ("·", r"\cdot "),
            ("λ", r"\lambda "), ("σ", r"\sigma "), ("−", "-"), ("²", "^{2}"),
            ("^(x-1)", "^{x-1}"), ("^(−b)", "^{-b}"), ("^(-b)", "^{-b}"),
            ("V^1.003", "V^{1.003}"), ("%", r"\%"), ("R^{2}", "R^2"),
            ("95\\% CI", r"\ \text{95\% CI}\ "), ("e^(x-1)", "e^{x-1}")]
    for a, b in reps:
        t = t.replace(a, b)
    return t


FIGS = {
    "The residual is reproducible": ("fig1_residual", 0.72,
        "The lexical seam: single-ZM residual on Shakespeare with the step-2 correction overlaid."),
    "The smooth two-regime model absorbs": ("fig3_s_law", 0.72,
        "The width law: seam width $s$ against vocabulary $V$ for 25 classic corpora (circles) and three modern registers (squares), with $s=0.0118\\,V$."),
    "Mandelbrot's c is sampling depth": ("fig2_c_collapse", 0.72,
        "Subsampling collapses the fitted shift $c$ along smooth trajectories."),
    "A generative account": ("fig4_b_prediction", 0.5,
        "The exponent predicted from the count histogram alone versus the exponent fitted to the rank curve (25 corpora)."),
}
EXTRA_FIGS = [
    ("fig5_trajectories", 0.72,
     "Downward $c(T)$ trajectories: binomially thinned data (circles) versus the mixture fitted at full depth (squares)."),
    ("fig6_heaps", 0.5,
     "Heaps' exponent against $1/b$: the measured Zipf--Heaps duality (r = 0.991)."),
]

PRE = r"""\documentclass[11pt]{article}
\usepackage[margin=1.1in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\graphicspath{{figures/}}
\setlength{\parskip}{2pt}
"""


def main():
    lines = MD.read_text(encoding="utf-8").splitlines()
    out = [PRE]
    i = 0
    title = lines[0].lstrip("# ").strip()
    out.append("\\title{" + prose(title) + "}")
    out.append("\\author{Grigori Karapetyan\\\\ \\small Independent researcher; Nexus Computers LLC; Burbank, CA, USA}")
    out.append("\\date{Draft v6.1 --- July 23, 2026}")
    out.append("\\begin{document}\\maketitle")
    while i < len(lines) and not lines[i].startswith("## Abstract"):
        i += 1
    i += 1
    ab = []
    while i < len(lines) and not lines[i].startswith("---"):
        ab.append(lines[i])
        i += 1
    out.append("\\begin{abstract}\n" + prose(" ".join(x.strip() for x in ab if x.strip())) + "\n\\end{abstract}")

    pending_extra = list(EXTRA_FIGS)
    body = lines[i:]
    j = 0
    in_table = []
    in_list = None  # 'i' or 'e'
    par = []

    def flush_par():
        nonlocal par
        if par:
            out.append(prose(" ".join(par)) + "\n")
            par = []

    def close_list():
        nonlocal in_list
        if in_list:
            out.append("\\end{itemize}" if in_list == "i" else "\\end{enumerate}")
            in_list = None

    def flush_table():
        nonlocal in_table
        if not in_table:
            return
        hdr = [c.strip() for c in in_table[0].strip("|").split("|")]
        out.append("\\begin{center}\\small\\begin{tabular}{" + "l" * len(hdr) + "}\\toprule")
        out.append(" & ".join(prose(h) for h in hdr) + r" \\ \midrule")
        for row in in_table[2:]:
            cells = [c.strip() for c in row.strip("|").split("|")]
            cells += [""] * (len(hdr) - len(cells))
            out.append(" & ".join(prose(c) for c in cells[:len(hdr)]) + r" \\")
        out.append("\\bottomrule\\end{tabular}\\end{center}")
        in_table = []

    while j < len(body):
        ln = body[j]
        s = ln.strip()
        if s.startswith("|"):
            flush_par(); close_list()
            in_table.append(s); j += 1; continue
        if in_table:
            flush_table()
        if not s:
            flush_par(); close_list(); j += 1; continue
        if s.startswith("---"):
            flush_par(); close_list(); j += 1; continue
        if ln.startswith("## "):
            flush_par(); close_list()
            h = s[3:].strip()
            if h.lower().startswith("references"):
                out.append("\\section*{References}")
            elif h.lower().startswith("appendix"):
                out.append("\\section*{" + prose(h) + "}")
            else:
                out.append("\\section{" + prose(re.sub(r"^\d+\.\s*", "", h)) + "}")
            j += 1; continue
        if ln.startswith("### "):
            flush_par(); close_list()
            h = re.sub(r"^\d+\.\d+\s*[—-]?\s*", "", s[4:].strip())
            out.append("\\subsection{" + prose(h) + "}")
            for key, (fname, w, cap) in list(FIGS.items()):
                if key.lower() in s.lower():
                    out.append("\\begin{figure}[t]\\centering\\includegraphics[width=%s\\linewidth]{%s}"
                               % (w, fname))
                    out.append("\\caption{" + cap + "}\\end{figure}")
                    if "generative" in key.lower():
                        for fname2, w2, cap2 in pending_extra:
                            out.append("\\begin{figure}[t]\\centering\\includegraphics[width=%s\\linewidth]{%s}"
                                       % (w2, fname2))
                            out.append("\\caption{" + cap2 + "}\\end{figure}")
                        pending_extra = []
                    del FIGS[key]
            j += 1; continue
        if re.match(r"^\s{4,}\S", ln) and not s.startswith("-"):
            flush_par(); close_list()
            out.append("\\begin{equation}" + mathify(prose_math_guard(s)) + "\\end{equation}")
            j += 1; continue
        m = re.match(r"^(\d+)\.\s+(.*)", s)
        if m and (in_list == "e" or looks_like_list(body, j)):
            flush_par()
            if in_list != "e":
                close_list(); out.append("\\begin{enumerate}"); in_list = "e"
            out.append("\\item " + prose(m.group(2)))
            j += 1; continue
        if s.startswith("- "):
            flush_par()
            if in_list != "i":
                close_list(); out.append("\\begin{itemize}"); in_list = "i"
            out.append("\\item " + prose(s[2:]))
            j += 1; continue
        if in_list and ln.startswith(("   ", "  ")):
            out.append(prose(s))
            j += 1; continue
        par.append(s)
        j += 1
    flush_par(); flush_table(); close_list()
    out.append("\\end{document}")
    TEX.write_text("\n".join(out), encoding="utf-8")
    print(f"wrote {TEX} ({TEX.stat().st_size/1024:.0f} KB)")


def prose_math_guard(s: str) -> str:
    return s


def looks_like_list(body, j) -> bool:
    s = body[j].strip()
    if not re.match(r"^\d+\.\s", s):
        return False
    k = j + 1
    seen = 1
    while k < len(body) and seen < 2:
        t = body[k].strip()
        if re.match(r"^\d+\.\s", t):
            seen += 1
        elif t and not body[k].startswith(("   ", "  ")):
            break
        k += 1
    return seen >= 2


if __name__ == "__main__":
    sys.exit(main())
