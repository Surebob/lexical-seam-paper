"""F10 — the gate-vs-labels capstone: does the erf gate fitted to the RANK CURVE
(no identity information) match the labeled grammar-word mixing curve built from
ANNOTATIONS (no curve information)?

For each of the 10 label-complete corpora (EXP03 manual_main):
  1. pi(r) = cumulative share of closed-class types among my top-r ranked words,
     using annotation tags (closed = ADP, AUX, CCONJ, DET, PART, PRON, SCONJ).
  2. Fit an erf curve to pi(r): pi ~ [1 - erf((ln r - ln k_lab)/w_lab)]/2.
  3. Compare (k_lab, w_lab) against the canonical rank-curve erf gate (k, w_gate)
     from the 2026-04-18 sweep — fitted blind to labels.
  4. k_POS from labels (pi crosses 0.5); regress log k_POS ~ log V (claim: 0.545).
  5. Labeled ablation: remove labeled closed types -> fit ZM on content remainder;
     fit ZM on closed types alone; does the unit step-2 correction stop helping?
  6. Inter-annotator agreement on the Principia 3-tagger pilot.

Outputs: ../outputs/f10_per_corpus.csv, f10_summary.md
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
ANN = REPO / "data" / "annotations" / "manual_main" / "english"
PILOT = REPO / "data" / "annotations" / "tagger_pilot"
ARCH = Path(r"C:\Users\Greg Kara\Desktop\temporary\emlexperiment\results\s2_v3_windows_full_outputs_2026-04-18\s2_v3_per_fit_results.csv")
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]
CLOSED = {"ADP", "AUX", "CCONJ", "DET", "PART", "PRON", "SCONJ"}

CORPORA = {
    "aesops_fables": ("Aesop's Fables", "pg11339.txt"),
    "collected_poe": ("Collected Poe", "pg2147.txt"),
    "complete_sherlock_holmes": ("Complete Sherlock Holmes", "pg1661.txt"),
    "critique_of_pure_reason": ("Critique of Pure Reason", "pg4280.txt"),
    "democracy_in_america": ("Democracy in America", "pg815.txt"),
    "dubliners": ("Dubliners", "pg2814.txt"),
    "emile": ("Emile", "pg5427.txt"),
    "federalist_papers": ("Federalist Papers", "pg1404.txt"),
    "grimms_fairy_tales": ("Grimm's Fairy Tales", "pg2591.txt"),
    "principia_ethica": ("Principia Ethica", "pg53430.txt"),
}

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


def load_ranked_types(fname):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    st, e = 0, len(text)
    for m in SM:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            st = nl + 1 if nl != -1 else i
            break
    for m in EMK:
        i = text.find(m)
        if i != -1:
            e = i
            break
    counts = Counter(TOKEN_RE.findall(text[st:e].lower()))
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return [t for t, _ in ranked], np.array([f for _, f in ranked], dtype=float)


def load_labels(slug):
    tags = {}
    for line in (ANN / f"{slug}.jsonl").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            o = json.loads(line)
            tags[str(o["type"]).lower()] = str(o["tag"]).upper()
        except Exception:
            continue
    return tags


def fit_zm(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 768)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, d @ coef, float(c))
    return best[1], math.sqrt(best[0]), best[2]


def step2_helps(freqs):
    pred, rmse, _ = fit_zm(freqs)
    logf = np.log(freqs)
    resid = logf - pred
    logr = np.log(np.arange(1, len(freqs) + 1, dtype=float))
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
    scores = {k: float(np.sqrt(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
    w = min(scores, key=scores.get)
    return w, scores[w] < rmse, rmse


def main():
    erf_fit = {r["corpus"]: r for r in csv.DictReader(open(ARCH, newline="", encoding="utf-8")) if r["gate"] == "erf"}
    rows = []
    for slug, (name, fname) in CORPORA.items():
        types, freqs = load_ranked_types(fname)
        V = len(types)
        tags = load_labels(slug)
        n_lab = len(tags)
        # coverage window: largest R such that >=85% of top-R types are labeled
        labeled_flag = np.array([1 if t in tags else 0 for t in types])
        cum_cov = np.cumsum(labeled_flag) / np.arange(1, V + 1)
        R = 50
        for r in range(50, min(V, 4000)):
            if cum_cov[r - 1] >= 0.85:
                R = r
        closed_flag = np.array([1.0 if tags.get(t, "X") in CLOSED else 0.0 for t in types])
        # INSTANTANEOUS closed-share at rank r (geometric window r/1.4 .. r*1.4),
        # which is what the gate sigma(log r) actually corresponds to. The cumulative
        # share is kept separately only for the paper-protocol k_POS.
        rr = np.arange(1, R + 1, dtype=float)
        pi = np.empty(R)
        for i, r in enumerate(rr):
            lo_i = max(0, int(math.floor(r / 1.4)) - 1)
            hi_i = min(R, int(math.ceil(r * 1.4)))
            pi[i] = float(np.mean(closed_flag[lo_i:hi_i]))
        pi_cum = np.cumsum(closed_flag)[:R] / np.arange(1, R + 1)

        def resid_fn(p):
            lk, lw = p
            z = (np.log(rr) - lk) / math.exp(lw)
            return (0.5 * (1.0 - sp_erf(z))) - pi

        best = None
        for lk0 in (2.0, 3.0, 4.0, 5.0):
            try:
                sol = least_squares(resid_fn, x0=[lk0, math.log(1.0)], method="lm", max_nfev=4000)
            except Exception:
                continue
            r2 = float(np.sqrt(np.mean(sol.fun**2)))
            if best is None or r2 < best[0]:
                best = (r2, sol.x)
        fit_rmse, (lk, lw) = best
        k_lab, w_lab = math.exp(lk), math.exp(lw)
        # k_POS: cumulative share crossing 0.5 (paper protocol)
        below = np.where(pi_cum < 0.5)[0]
        k_pos = float(rr[below[0]]) if len(below) else float(R)
        # gate curve match: evaluate canonical gate vs pi
        kg = float(erf_fit[name]["k"])
        wg = float(erf_fit[name]["w_gate"])
        z = (np.log(rr) - math.log(kg)) / wg
        gate_curve = 0.5 * (1.0 - sp_erf(z))
        gate_vs_pi_rmse = float(np.sqrt(np.mean((gate_curve - pi) ** 2)))
        # labeled ablation
        closed_types = {t for t in types if tags.get(t, "X") in CLOSED}
        content_freqs = np.array([f for t, f in zip(types, freqs) if t not in closed_types], dtype=float)
        closed_freqs = np.array(sorted((f for t, f in zip(types, freqs) if t in closed_types), reverse=True), dtype=float)
        w_full, helps_full, rmse_full = step2_helps(freqs)
        w_cont, helps_cont, rmse_cont = step2_helps(content_freqs)
        w_clos, helps_clos, rmse_clos = step2_helps(closed_freqs) if len(closed_freqs) > 30 else ("na", False, float("nan"))

        rows.append({
            "corpus": name, "V": V, "n_labels": n_lab, "coverage_R": R,
            "k_lab": k_lab, "w_lab": w_lab, "labfit_rmse": fit_rmse,
            "k_gate": kg, "w_gate": wg, "gate_vs_pi_rmse": gate_vs_pi_rmse,
            "log_k_ratio": math.log(kg / k_lab), "w_ratio": wg / w_lab,
            "k_pos": k_pos,
            "n_closed_types": len(closed_types),
            "winner_full": w_full, "helps_full": helps_full, "rmse_full": rmse_full,
            "winner_content": w_cont, "helps_content": helps_cont, "rmse_content": rmse_cont,
            "winner_closed": w_clos, "helps_closed": helps_clos, "rmse_closed": rmse_clos,
        })
        print(f"{name[:26]:27} R={R:5} k_lab={k_lab:7.1f} w_lab={w_lab:.3f} | k_gate={kg:7.1f} "
              f"w_gate={wg:.3f} | gate-vs-pi rmse={gate_vs_pi_rmse:.4f} | kPOS={k_pos:.0f} "
              f"| full:{w_full}/{helps_full} content:{w_cont}/{helps_cont}", flush=True)

    with open(OUT / "f10_per_corpus.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    lk_gate = np.log([r["k_gate"] for r in rows])
    lk_lab = np.log([r["k_lab"] for r in rows])
    wgates = np.array([r["w_gate"] for r in rows])
    wlabs = np.array([r["w_lab"] for r in rows])
    lV = np.log([r["V"] for r in rows])
    lkpos = np.log([r["k_pos"] for r in rows])
    X = np.column_stack([np.ones_like(lV), lV])
    bpos, *_ = np.linalg.lstsq(X, lkpos, rcond=None)
    sse = float(np.sum((lkpos - X @ bpos) ** 2))
    se = math.sqrt(sse / (len(lV) - 2) / float(np.sum((lV - lV.mean()) ** 2)))

    # pilot agreement
    manual = load_labels("principia_ethica")
    agree_lines = []
    for tf in sorted(PILOT.glob("tagger*_principia_ethica.jsonl")):
        tt = {}
        for line in tf.read_text(encoding="utf-8").splitlines():
            try:
                o = json.loads(line)
                tt[str(o["type"]).lower()] = str(o["tag"]).upper()
            except Exception:
                continue
        common = [t for t in tt if t in manual]
        if not common:
            continue
        exact = sum(1 for t in common if tt[t] == manual[t]) / len(common)
        binary = sum(1 for t in common if (tt[t] in CLOSED) == (manual[t] in CLOSED)) / len(common)
        agree_lines.append(f"- {tf.stem} vs manual: exact tag {100*exact:.1f}%, closed/open {100*binary:.1f}% (n={len(common)})")

    lines = ["# F10 — gate vs labels (the identity capstone)\n"]
    lines.append(f"- corpora: {len(rows)} (all EXP03 label-complete)")
    lines.append(f"- corr(log k_gate, log k_lab) = {float(np.corrcoef(lk_gate, lk_lab)[0,1]):.4f}; "
                 f"median k_gate/k_lab = {float(np.exp(np.median(lk_gate - lk_lab))):.3f}")
    lines.append(f"- corr(w_gate, w_lab) = {float(np.corrcoef(wgates, wlabs)[0,1]):.4f}; "
                 f"median w_gate/w_lab = {float(np.median(wgates/wlabs)):.3f}")
    lines.append(f"- median gate-vs-pi curve RMSE = {float(np.median([r['gate_vs_pi_rmse'] for r in rows])):.4f} "
                 f"(pi is a share in [0,1])")
    lines.append(f"- labeled k_POS ~ V^beta: beta = {bpos[1]:.3f} 95% CI [{bpos[1]-1.96*se:.3f}, {bpos[1]+1.96*se:.3f}] "
                 f"(paper claim 0.545; n=10)")
    helps_full_n = sum(1 for r in rows if r["helps_full"])
    helps_cont_n = sum(1 for r in rows if r["helps_content"])
    lines.append(f"- labeled ablation: unit step-2 helps on full corpus {helps_full_n}/10 -> on content-only {helps_cont_n}/10")
    lines.append("\n## Inter-annotator pilot (Principia, 3 taggers vs manual)")
    lines.extend(agree_lines if agree_lines else ["- (no overlap parsed)"])
    (OUT / "f10_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
