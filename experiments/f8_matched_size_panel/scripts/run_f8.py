"""F8 — matched-size cross-language panel (the size-confound killer).

f6 proved ZM's c is dominated by sampling depth, which contaminates every
cross-corpus comparison made at unequal sizes — including the old Table 3 and the
f4e sigma_T gradient. Here every corpus is cut to the SAME depth (first 65,000
tokens) and re-measured on equal footing:

  - ZM (b, c) and the unit-amplitude step-2 winner
  - the PLN mixture signature (pi_H, sd_H, sd_T) at matched depth

Corpora: 8 English classics + 7 canonical multilingual + 6 langext (Polish noted
separately: only 33k tokens available, excluded from the matched tier).

Outputs: ../outputs/f8_matched_panel.csv, f8_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
ML1 = REPO / "data" / "zipf_multilang"
ML2 = REPO / "data" / "zipf_multilang_romance"
LX = REPO / "data" / "langext"
OUT = HERE.parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

T_STAR = 65_000
EN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
UNI_RE = re.compile(r"[^\W\d_]+", re.UNICODE)
SM = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
EMK = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

ENGLISH = [
    ("EN Shakespeare", DATA / "pg100.txt"), ("EN War and Peace (tr)", DATA / "pg2600.txt"),
    ("EN KJ Bible", DATA / "pg10.txt"), ("EN Pride and Prejudice", DATA / "pg1342.txt"),
    ("EN Moby Dick", DATA / "pg2701.txt"), ("EN Don Quixote (tr)", DATA / "pg996.txt"),
    ("EN Ulysses", DATA / "pg4300.txt"), ("EN Federalist", DATA / "pg1404.txt"),
]
MULTI = [
    ("RU War and Peace", ML1 / "russian_war_and_peace" / "combined_clean.txt", "uni"),
    ("ZH Three Kingdoms", ML1 / "mandarin_three_kingdoms" / "combined_clean.txt", "jieba"),
    ("AR 1001 Nights", ML1 / "arabic_1001_nights" / "combined_clean.txt", "uni"),
    ("LA Gallic Wars", ML1 / "latin_gallic_wars" / "combined_clean.txt", "uni"),
    ("FR Les Miserables", ML2 / "french_les_miserables" / "combined_clean.txt", "uni"),
    ("ES Don Quixote", ML2 / "spanish_don_quixote" / "combined_clean.txt", "uni"),
    ("NL Max Havelaar", ML2 / "dutch_max_havelaar" / "combined_clean.txt", "uni"),
    ("IT Promessi Sposi", LX / "italian.txt", "uni"),
    ("FI Seitseman veljesta", LX / "finnish.txt", "uni"),
    ("PT Dom Casmurro", LX / "portuguese.txt", "uni"),
    ("DE Wahlverwandtschaften", LX / "german.txt", "uni"),
    ("SV Gosta Berling", LX / "swedish.txt", "uni"),
    ("PL Pan Tadeusz (33k only)", LX / "polish.txt", "uni"),
]

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}
GH_X, GH_W = hermgauss(48)
LOG_GH_W = np.log(GH_W)


def strip_gutenberg(text):
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
    return text[st:e]


def tokens_for(path, mode):
    text = path.read_text(encoding="utf-8", errors="ignore")
    text = strip_gutenberg(text)
    if mode == "en":
        return EN_RE.findall(text.lower())
    if mode == "jieba":
        import jieba
        return [t for t in jieba.cut(text) if UNI_RE.fullmatch(t)]
    return UNI_RE.findall(text.lower())


def fit_zm_full(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=float)
    logf = np.log(freqs)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 1024)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(-coef[1]), float(c), d @ coef)
    _, b, c, pred = best
    resid = logf - pred
    logr = np.log(ranks)
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
    unit = {k: float(np.sqrt(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
    return b, c, min(unit, key=unit.get)


def pln_log_pn(ns, mu, sd):
    loglam = mu + math.sqrt(2.0) * sd * GH_X
    lam = np.exp(np.clip(loglam, -30, 30))
    n = ns[:, None]
    ll = n * loglam[None, :] - lam[None, :] - gammaln(n + 1.0)
    return logsumexp(ll + LOG_GH_W[None, :], axis=1) - 0.5 * math.log(math.pi)


def mixture_nll(theta, uniq, mult):
    pi_h = 0.001 + 0.299 / (1.0 + math.exp(-theta[0]))
    mu_t = theta[1]
    sd_t = math.exp(theta[2])
    mu_h = mu_t + 0.5 + math.exp(theta[3])
    sd_h = math.exp(theta[4])
    if not (0.03 < sd_t < 6.0 and 0.03 < sd_h < 6.0):
        return 1e12
    lp_h = pln_log_pn(uniq, mu_h, sd_h)
    lp_t = pln_log_pn(uniq, mu_t, sd_t)
    lp = np.logaddexp(math.log(pi_h) + lp_h, math.log(1 - pi_h) + lp_t)
    z = np.array([0.0])
    lp0 = np.logaddexp(math.log(pi_h) + pln_log_pn(z, mu_h, sd_h),
                       math.log(1 - pi_h) + pln_log_pn(z, mu_t, sd_t))[0]
    return -float(np.dot(mult, lp - math.log1p(-min(math.exp(lp0), 1 - 1e-12))))


def fit_pln(freqs):
    uniq, mult = np.unique(freqs, return_counts=True)
    mult = mult.astype(float)
    best = None
    for x0 in [
        [-3.0, 0.3, math.log(1.2), math.log(4.0), math.log(1.0)],
        [-2.0, 0.0, math.log(1.0), math.log(6.0), math.log(1.4)],
        [-4.0, 0.5, math.log(1.5), math.log(3.0), math.log(0.7)],
    ]:
        try:
            r = minimize(mixture_nll, x0=x0, args=(uniq, mult), method="Nelder-Mead",
                         options={"maxiter": 3000, "xatol": 1e-5, "fatol": 1e-7})
        except Exception:
            continue
        if best is None or r.fun < best.fun:
            best = r
    t = best.x
    return {
        "pi_h": 0.001 + 0.299 / (1.0 + math.exp(-t[0])),
        "sd_t": float(math.exp(t[2])), "sd_h": float(math.exp(t[4])),
    }


def main():
    rows = []
    jobs = [(n, p, "en") for n, p in ENGLISH] + list(MULTI)
    for name, path, mode in jobs:
        try:
            toks = tokens_for(path, mode)
        except Exception as ex:
            print(f"FAIL {name}: {ex!r}", flush=True)
            continue
        matched = len(toks) >= T_STAR
        use = toks[:T_STAR] if matched else toks
        freqs = np.array(sorted(Counter(use).values(), reverse=True), dtype=float)
        b, c, winner = fit_zm_full(freqs)
        mix = fit_pln(freqs)
        rows.append({"corpus": name, "matched": matched, "tokens": len(use), "V": len(freqs),
                     "zm_b": b, "zm_c": c, "winner": winner, **mix,
                     "hapax": float(np.sum(freqs == 1.0)) / len(freqs)})
        r = rows[-1]
        print(f"{name:26} T={r['tokens']:6} V={r['V']:6} c={r['zm_c']:7.2f} b={r['zm_b']:.3f} "
              f"win={r['winner']:6} sd_T={r['sd_t']:.3f} pi_H={100*r['pi_h']:.2f}%", flush=True)

    with open(OUT / "f8_matched_panel.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    en = [r for r in rows if r["corpus"].startswith("EN ") and r["matched"]]
    xx = [r for r in rows if not r["corpus"].startswith("EN ") and r["matched"]]
    lines = ["# F8 — matched-size panel (all corpora at 65k tokens)\n"]
    lines.append(f"- English at 65k: c median {float(np.median([r['zm_c'] for r in en])):.2f} "
                 f"(range {min(r['zm_c'] for r in en):.2f}-{max(r['zm_c'] for r in en):.2f}); "
                 f"sd_T median {float(np.median([r['sd_t'] for r in en])):.3f}")
    lines.append(f"- Non-English at 65k: c median {float(np.median([r['zm_c'] for r in xx])):.2f} "
                 f"(range {min(r['zm_c'] for r in xx):.2f}-{max(r['zm_c'] for r in xx):.2f}); "
                 f"sd_T median {float(np.median([r['sd_t'] for r in xx])):.3f}")
    wc_en = Counter(r["winner"] for r in en)
    wc_xx = Counter(r["winner"] for r in xx)
    lines.append(f"- winner map English@65k: {dict(wc_en)}; non-English@65k: {dict(wc_xx)}")
    lines.append("\n| corpus | V | c | b | winner | sd_T | pi_H% | hapax |")
    lines.append("|---|---:|---:|---:|---|---:|---:|---:|")
    for r in sorted(rows, key=lambda r: r["sd_t"]):
        lines.append(f"| {r['corpus']} | {r['V']} | {r['zm_c']:.2f} | {r['zm_b']:.3f} | {r['winner']} | "
                     f"{r['sd_t']:.3f} | {100*r['pi_h']:.2f} | {r['hapax']:.3f} |")
    (OUT / "f8_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
