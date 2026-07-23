"""Build v6 figures (paper/figures/) and bootstrap CIs (docs/v6_bootstrap_cis.md)."""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FIG = REPO / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
E = REPO / "experiments"
rng = np.random.default_rng(20260722)

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
C1, C2, C3 = "#C03422", "#2F5DA8", "#0F7A4D"


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print(f"fig: {name}", flush=True)


# ---------- F1: Shakespeare residual + correction ----------
TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
text = (REPO / "data" / "zipf" / "pg100.txt").read_text(encoding="utf-8", errors="ignore")
i = text.find("*** START OF THE PROJECT GUTENBERG EBOOK")
j = text.find("*** END OF THE PROJECT GUTENBERG EBOOK")
text = text[text.find("\n", i) + 1:j]
freqs = np.array(sorted(Counter(TOKEN_RE.findall(text.lower())).values(), reverse=True), dtype=float)
V = len(freqs)
ranks = np.arange(1, V + 1, dtype=float)
logf = np.log(freqs)
best = None
for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 1024)]):
    d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
    coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
    mse = float(np.mean((d @ coef - logf) ** 2))
    if best is None or mse < best[0]:
        best = (mse, d @ coef)
resid = logf - best[1]
logr = np.log(ranks)
x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
is_curve = (x - 1.0) - np.log(x)
fig, ax = plt.subplots(figsize=(4.6, 3.0))
sel = np.unique(np.geomspace(1, V, 400).astype(int)) - 1
ax.plot(ranks[sel], resid[sel], ".", ms=2.5, color=C2, alpha=.7, label="single-ZM residual")
ax.plot(ranks[sel], is_curve[sel], "-", lw=1.6, color=C1, label=r"$(x-1)-\log x$")
ax.set_xscale("log")
ax.set_xlabel("rank"); ax.set_ylabel("residual (log frequency)")
ax.set_title("The seam: single-ZM residual, Shakespeare", fontsize=9)
ax.legend(frameon=False, fontsize=8)
save(fig, "fig1_residual")

# ---------- F2: c collapse ----------
rows = list(csv.DictReader(open(E / "f6_c_sampling_depth" / "outputs" / "f6_slices.csv", encoding="utf-8")))
fig, ax = plt.subplots(figsize=(4.6, 3.0))
for name, col in [("Complete Works of Shakespeare", C1), ("Les Miserables", C2), ("Don Quixote", C3)]:
    sub = sorted([r for r in rows if r["corpus"] == name], key=lambda r: int(r["tokens"]))
    ax.plot([int(r["tokens"]) for r in sub], [float(r["zm_c"]) + 1 for r in sub],
            "o-", color=col, lw=1.5, ms=4, label=name.replace("Complete Works of ", ""))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("tokens counted"); ax.set_ylabel("Zipf–Mandelbrot $c$ (+1)")
ax.set_title("Mandelbrot's $c$ collapses with sampling depth", fontsize=9)
ax.legend(frameon=False, fontsize=8)
save(fig, "fig2_c_collapse")

# ---------- F3: the s-law ----------
f2rows = list(csv.DictReader(open(E / "f2_k_profile_likelihood" / "outputs" / "f2_per_corpus.csv", encoding="utf-8")))
Ve = np.array([float(r["V"]) for r in f2rows]); Se = np.array([float(r["s_at_min"]) for r in f2rows])
panel = {r["corpus"]: r for r in csv.DictReader(open(E / "f5_expanded_panel" / "outputs" / "f5b_gate_fits.csv", encoding="utf-8"))}
regs = [("brown", "Brown"), ("wikitext_1M", "WikiText"), ("cornell_dialogs", "Dialogue")]
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.plot(Ve, Se, "o", ms=4, color=C2, label="25 classic English corpora")
for key, lab in regs:
    r = panel[key]
    ax.plot(float(r["V"]), float(r["s"]), "s", ms=6, color=C1)
    ax.annotate(lab, (float(r["V"]), float(r["s"])), textcoords="offset points",
                xytext=(5, -3), fontsize=7, color=C1)
vv = np.geomspace(3000, 60000, 50)
ax.plot(vv, 0.0118 * vv, "--", color="gray", lw=1.2, label=r"$s = 0.0118\,V$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("vocabulary size $V$"); ax.set_ylabel("seam width $s = k\\,w_{tail}$")
ax.set_title("The width law: $s \\approx 1.2\\%$ of vocabulary", fontsize=9)
ax.legend(frameon=False, fontsize=8)
save(fig, "fig3_s_law")

# ---------- F4: predicted vs actual b (f7 run log, commit a588ab4) ----------
PAIRS = [("Shakespeare",1.654,1.724),("War and Peace",1.688,1.681),("Moby Dick",1.239,1.182),
("KJ Bible",1.683,1.896),("Federalist",1.572,1.555),("Grimm",1.628,1.580),("Quixote",1.417,1.484),
("Pride & Prejudice",1.470,1.456),("Canterbury",1.191,1.256),("Arabian Nights",1.149,1.188),
("Aesop",1.129,1.147),("Sherlock",1.174,1.253),("Jane Eyre",1.208,1.258),("Dubliners",1.015,1.137),
("Iliad",1.378,1.488),("Democracy",1.402,1.521),("Origin",1.640,1.689),("Wealth",1.949,1.863),
("Les Misérables",1.442,1.526),("Decline & Fall",1.314,1.481),("Émile",1.374,1.464),
("Ulysses",1.038,1.045),("Poe",1.010,1.134),("Principia",1.900,1.815),("Critique",1.955,1.690)]
bp = np.array([p[1] for p in PAIRS]); ba = np.array([p[2] for p in PAIRS])
fig, ax = plt.subplots(figsize=(3.6, 3.4))
ax.plot([0.95, 2.05], [0.95, 2.05], "--", color="gray", lw=1)
ax.plot(bp, ba, "o", ms=5, color=C2, alpha=.85)
ax.set_xlabel("$b$ predicted from count histogram")
ax.set_ylabel("$b$ fitted to rank curve")
ax.set_title(f"Histogram predicts the exponent (r = {np.corrcoef(bp,ba)[0,1]:.3f})", fontsize=9)
save(fig, "fig4_b_prediction")

# ---------- F5: downward trajectories, model vs real ----------
tr = list(csv.DictReader(open(E / "f6b_heavy_tail_extrapolation" / "outputs" / "f6b_trajectories.csv", encoding="utf-8")))
fig, ax = plt.subplots(figsize=(4.6, 3.0))
for name, col in [("Complete Works of Shakespeare", C1), ("War and Peace", C2)]:
    for model, ls, lab in [("real_thinned", "o-", "real"), ("pln_fullfit", "s--", "model")]:
        sub = sorted([r for r in tr if r["corpus"] == name and r["model"] == model],
                     key=lambda r: float(r["scale"]))
        ax.plot([float(r["scale"]) for r in sub], [float(r["c"]) + 1 for r in sub], ls,
                color=col, lw=1.3, ms=4, alpha=.9,
                label=f"{name.replace('Complete Works of ','')} ({lab})")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("fraction of corpus"); ax.set_ylabel("$c$ (+1)")
ax.set_title("The mixture predicts the collapse", fontsize=9)
ax.legend(frameon=False, fontsize=7)
save(fig, "fig5_trajectories")

# ---------- F6: Heaps duality ----------
hp = list(csv.DictReader(open(E / "f11_loose_ends" / "outputs" / "f11_heaps.csv", encoding="utf-8")))
hb = np.array([float(r["heaps_beta"]) for r in hp]); bb = np.array([float(r["zm_b"]) for r in hp])
fig, ax = plt.subplots(figsize=(3.6, 3.4))
ax.plot(1 / bb, hb, "o", ms=5, color=C3, alpha=.85)
lim = [0.35, 1.05]
ax.plot(lim, lim, "--", color="gray", lw=1)
ax.set_xlabel(r"$1/b$ (Zipf)"); ax.set_ylabel(r"Heaps exponent $\beta$")
ax.set_title(f"Zipf–Heaps duality (r = {np.corrcoef(hb, 1/bb)[0,1]:.3f})", fontsize=9)
save(fig, "fig6_heaps")

# ---------- F7: family fingerprints (s vs V across systems) ----------
mixrows = [r for r in csv.DictReader(open(E / "f12_forced_mixing" / "outputs" / "f12_mixtures.csv", encoding="utf-8")) if int(r["m"]) > 1]
concat = list(csv.DictReader(open(E / "f15_deep_multilingual" / "outputs" / "f15c_concat.csv", encoding="utf-8")))
ga = list(csv.DictReader(open(E / "f14_ga_break" / "outputs" / "f14_ga_fits.csv", encoding="utf-8")))
nulls = list(csv.DictReader(open(E / "f14_ga_break" / "outputs" / "f14b_nulls.csv", encoding="utf-8")))
simon_decay = [r for r in nulls if r["variant"] == "pure_simon_decay"]
simon_classic = [r for r in nulls if r["variant"] == "classic_simon"]
sur = panel["census_surnames"]
BELT_V, BELT_S = 1_552_974, 25_749  # x1_asteroid_control/outputs/x1_summary.md

fig, ax = plt.subplots(figsize=(5.2, 3.6))
vv7 = np.geomspace(3000, 2_200_000, 60)
ax.plot(vv7, 0.0118 * vv7, "--", color="gray", lw=1.2, zorder=1)
ax.annotate(r"$s = 0.0118\,V$", (5.5e5, 0.0118 * 5.5e5), textcoords="offset points",
            xytext=(2, -16), fontsize=8, color="gray")
ax.plot(Ve, Se, "o", ms=3.5, color=C2, alpha=.8, label="25 classic English corpora", zorder=3)
for key in ["brown", "wikitext_1M", "cornell_dialogs"]:
    r = panel[key]
    ax.plot(float(r["V"]), float(r["s"]), "s", ms=5, color=C1, zorder=3,
            label="modern registers" if key == "brown" else None)
ax.plot([float(r["V"]) for r in mixrows], [float(r["s"]) for r in mixrows], "D", ms=4,
        mfc="none", mec=C2, label="forced mixtures $m=2$–$14$", zorder=3)
deep7 = [r for r in concat if float(r["tok_per_V"]) >= 15]
shal7 = [r for r in concat if float(r["tok_per_V"]) < 15]
ax.plot([float(r["V"]) for r in deep7], [float(r["s"]) for r in deep7], "^", ms=6, color=C3,
        label="languages at matched depth", zorder=3)
ax.plot([float(r["V"]) for r in shal7], [float(r["s"]) for r in shal7], "^", ms=6, mfc="none",
        mec=C3, label="languages below matched depth", zorder=3)
sims7 = ga + simon_decay
ax.plot([float(r["V"]) for r in sims7], [float(r["s"]) for r in sims7], "x", ms=4.5, color="#666666",
        label="simulated decaying-innovation", zorder=2)
ax.plot([float(r["V"]) for r in simon_classic], [float(r["s"]) for r in simon_classic], "+", ms=7,
        color="#B8860B", label="simulated constant-innovation (0.0101)", zorder=2)
ax.plot(float(sur["V"]), float(sur["s"]), "*", ms=10, color="#7A1F1F", zorder=4)
ax.annotate("surnames\n$s/V = 0.0266$", (float(sur["V"]), float(sur["s"])), textcoords="offset points",
            xytext=(9, 2), fontsize=7, color="#7A1F1F", ha="left")
ax.plot(BELT_V, BELT_S, "p", ms=7, color="black", zorder=4)
ax.annotate("asteroid belt\n$s/V = 0.0166$", (BELT_V, BELT_S), textcoords="offset points",
            xytext=(-10, -4), fontsize=7, ha="right")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("vocabulary / catalogue size $V$"); ax.set_ylabel("seam width $s$")
ax.set_title("Family fingerprints: one law for language, parallel lines for others", fontsize=9)
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
save(fig, "fig7_fingerprints")

# ---------- Bootstrap CIs ----------
def boot_ci(stat, n_items, reps=10000):
    vals = []
    for _ in range(reps):
        idx = rng.integers(0, n_items, n_items)
        vals.append(stat(idx))
    v = np.sort(np.array(vals))
    return float(v[int(0.025 * reps)]), float(v[int(0.975 * reps)])

lines = ["# v6 bootstrap 95% CIs (resampling corpora, 10,000 reps)\n"]
lV, lS = np.log(Ve), np.log(Se)
def slope_stat(idx):
    X = np.column_stack([np.ones(len(idx)), lV[idx]])
    return float(np.linalg.lstsq(X, lS[idx], rcond=None)[0][1])
lo, hi = boot_ci(slope_stat, len(Ve))
lines.append(f"- s-law exponent (log s ~ log V, 25 EN): point 1.003, bootstrap CI [{lo:.3f}, {hi:.3f}]")
sv = Se / Ve
lo, hi = boot_ci(lambda i: float(np.median(sv[i])), len(sv))
lines.append(f"- median s/V: point {np.median(sv):.4f}, CI [{lo:.4f}, {hi:.4f}]")
f1rows = list(csv.DictReader(open(E / "f1_fresh_reproduction" / "outputs" / "f1_per_corpus.csv", encoding="utf-8")))
impr = np.array([float(r["lzm_impr_pct"]) for r in f1rows])
lo, hi = boot_ci(lambda i: float(np.median(impr[i])), len(impr))
lines.append(f"- λ-ZM median RMSE improvement (25 EN): point {np.median(impr):.2f}%, CI [{lo:.2f}%, {hi:.2f}%]")
lo, hi = boot_ci(lambda i: float(np.corrcoef(bp[i], ba[i])[0, 1]) if len(set(i)) > 2 else 0.9, len(bp))
lines.append(f"- corr(b_pred, b_actual): point {np.corrcoef(bp,ba)[0,1]:.3f}, CI [{lo:.3f}, {hi:.3f}]")
lo, hi = boot_ci(lambda i: float(np.corrcoef(hb[i], 1/bb[i])[0, 1]) if len(set(i)) > 2 else 0.99, len(hb))
lines.append(f"- corr(Heaps β, 1/b): point {np.corrcoef(hb,1/bb)[0,1]:.3f}, CI [{lo:.3f}, {hi:.3f}]")
(REPO / "docs" / "v6_bootstrap_cis.md").write_text("\n".join(lines), encoding="utf-8")
print("\n".join(lines).encode("ascii", "replace").decode(), flush=True)
