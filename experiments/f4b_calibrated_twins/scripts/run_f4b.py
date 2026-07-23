"""F4b — Calibrated synthetic twins: does each real corpus's own type-frequency
structure reproduce the erf-gate preference and the seam fingerprint?

For each of 12 representative English corpora:
  1. Fit a 2-component Gaussian mixture (EM) to the corpus's type log-frequencies.
  2. Simulate a "lognormal twin": V type-rates drawn from the fitted mixture,
     Poisson counts, zeros dropped.
  3. Simulate a "pareto twin": two power-law type-rate populations with rank-curve
     exponents matched to the f1 two-population fit (alpha_function on top-100 mass,
     alpha_content on the rest), levels moment-matched to the corpus.
  4. Fit all 5 gates (16 starts) to each twin; record winner, erf margin, BIC spread.
  5. Also record each twin's ZM c and unit-amplitude step-2 winner
     ({IS, exp, euclid, xpow}) to compare the full seam fingerprint with the corpus.

Outputs: ../outputs/f4b_fits.csv, f4b_twins.csv, f4b_summary.md
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf as sp_erf

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DATA = REPO / "data" / "zipf"
OUT = HERE.parent / "outputs"
F1CSV = REPO / "experiments" / "f1_fresh_reproduction" / "outputs" / "f1_per_corpus.csv"

TOKEN_RE = re.compile(r"[a-z]+(?:'[a-z]+)?")
START_MARKERS = ["*** START OF THE PROJECT GUTENBERG EBOOK", "*** START OF THIS PROJECT GUTENBERG EBOOK"]
END_MARKERS = ["*** END OF THE PROJECT GUTENBERG EBOOK", "*** END OF THIS PROJECT GUTENBERG EBOOK"]

CORPORA = [
    ("Complete Works of Shakespeare", "pg100.txt"), ("War and Peace", "pg2600.txt"),
    ("King James Bible", "pg10.txt"), ("Federalist Papers", "pg1404.txt"),
    ("Origin of Species", "pg1228.txt"), ("Wealth of Nations", "pg3300.txt"),
    ("Moby Dick", "pg2701.txt"), ("Pride and Prejudice", "pg1342.txt"),
    ("Dubliners", "pg2814.txt"), ("Ulysses", "pg4300.txt"),
    ("Grimm's Fairy Tales", "pg2591.txt"), ("Don Quixote", "pg996.txt"),
]

SEED = 20260721
N_STARTS = 16
MAX_NFEV = 8000
LOWER = np.array([-100.0, 0.5, 0.0, -100.0, 0.5, 0.0, 20.0, 0.05, 0.05])
UPPER = np.array([100.0, 3.0, 1000.0, 100.0, 3.0, 1000.0, 1000.0, 10.0, 10.0])

GENS = {
    "is": lambda x: (x - 1.0) - np.log(x),
    "exp": lambda x: np.exp(x - 1.0) - x,
    "euclid": lambda x: (1.0 - x) ** 2,
    "xpow": lambda x: np.power(x, x) - np.sqrt(x),
}


def load_freqs(fname: str):
    text = (DATA / fname).read_text(encoding="utf-8", errors="ignore")
    start, end = 0, len(text)
    for m in START_MARKERS:
        i = text.find(m)
        if i != -1:
            nl = text.find("\n", i)
            start = nl + 1 if nl != -1 else i
            break
    for m in END_MARKERS:
        i = text.find(m)
        if i != -1:
            end = i
            break
    counts = Counter(TOKEN_RE.findall(text[start:end].lower()))
    return np.array(sorted(counts.values(), reverse=True), dtype=np.float64)


def gmm2_em(y: np.ndarray, iters=250, seed=0):
    rng = np.random.default_rng(seed)
    mu = np.array([np.quantile(y, 0.97), np.quantile(y, 0.35)])
    sd = np.array([1.0, 1.0])
    pi = np.array([0.05, 0.95])
    for _ in range(iters):
        logp = -0.5 * ((y[:, None] - mu[None, :]) / sd[None, :]) ** 2 - np.log(sd[None, :]) + np.log(pi[None, :])
        logp -= logp.max(axis=1, keepdims=True)
        w = np.exp(logp)
        w /= w.sum(axis=1, keepdims=True)
        nk = w.sum(axis=0) + 1e-12
        pi = nk / len(y)
        mu = (w * y[:, None]).sum(axis=0) / nk
        sd = np.sqrt((w * (y[:, None] - mu[None, :]) ** 2).sum(axis=0) / nk) + 1e-6
    order = np.argsort(-mu)
    return pi[order], mu[order], sd[order]


def simulate_lognormal_twin(freqs, rng):
    y = np.log(freqs)
    pi, mu, sd = gmm2_em(y, seed=int(rng.integers(2**31)))
    V = len(freqs)
    comp = rng.random(V) < pi[0]
    lam = np.where(comp, np.exp(rng.normal(mu[0], sd[0], V)), np.exp(rng.normal(mu[1], sd[1], V)))
    counts = rng.poisson(lam)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(np.float64), {"pi_head": float(pi[0]), "mu_head": float(mu[0]),
                                                      "sd_head": float(sd[0]), "mu_tail": float(mu[1]),
                                                      "sd_tail": float(sd[1])}


def simulate_pareto_twin(freqs, a_func, a_cont, rng):
    V = len(freqs)
    n_h = 130
    a_h = max(1.0 / max(a_func, 0.2), 0.3)
    a_t = max(1.0 / max(a_cont, 0.2), 0.3)
    raw_h = (1.0 + rng.pareto(a_h, n_h))
    raw_t = (1.0 + rng.pareto(a_t, V - n_h))
    raw_h *= np.median(freqs[:100]) / np.median(raw_h)
    tail_ref = np.median(freqs[100:])
    raw_t *= max(tail_ref, 1.2) / np.median(raw_t)
    lam = np.concatenate([raw_h, raw_t])
    counts = rng.poisson(lam)
    counts = counts[counts > 0]
    return np.sort(counts)[::-1].astype(np.float64)


def zm_and_step2(freqs):
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    logf = np.log(freqs)
    logr = np.log(ranks)
    best = None
    for c in np.concatenate([[0.0], np.geomspace(1e-6, V, 1024)]):
        d = np.column_stack([np.ones_like(logf), np.log(ranks + c)])
        coef, *_ = np.linalg.lstsq(d, logf, rcond=None)
        mse = float(np.mean((d @ coef - logf) ** 2))
        if best is None or mse < best[0]:
            best = (mse, float(coef[0]), float(-coef[1]), float(c), d @ coef)
    _, a, b, c, pred = best
    resid = logf - pred
    x = 0.05 + 0.95 * (logr - logr[0]) / (logr[-1] - logr[0])
    unit = {k: float(np.sqrt(np.mean((resid - g(x)) ** 2))) for k, g in GENS.items()}
    return c, min(unit, key=unit.get)


def gate_sigma(gate, z):
    if gate == "logistic":
        return 1.0 / (1.0 + np.exp(np.clip(z, -60, 60)))
    if gate == "tanh":
        return 0.5 * (1.0 - np.tanh(z))
    if gate == "erf":
        return 0.5 * (1.0 - sp_erf(z))
    if gate == "algebraic":
        return 0.5 * (1.0 - z / np.sqrt(1.0 + z * z))
    if gate == "arctan":
        return 0.5 - np.arctan(z) / math.pi
    raise ValueError(gate)


def predict(ranks, log_ranks, params, gate):
    a1, b1, c1, a2, b2, c2, k, w_gate, w_tail = [float(v) for v in params]
    z = (log_ranks - math.log(k)) / w_gate
    sigma = gate_sigma(gate, z)
    head = a1 - b1 * np.log(ranks + max(c1, 0.0))
    scale = max(1.0, k * w_tail)
    zz = np.clip((ranks - k) / scale, -60.0, 60.0)
    tail_rank = 1.0 + scale * np.log1p(np.exp(zz))
    tail = a2 - b2 * np.log(tail_rank + max(c2, 0.0))
    return sigma * head + (1.0 - sigma) * tail


def piecewise_anchor(ranks, logf):
    k = min(500, len(ranks) // 3)

    def zmfit(rr, yy):
        best = None
        for c in np.concatenate([[0.0], np.geomspace(1e-6, rr[-1], 384)]):
            d = np.column_stack([np.ones_like(yy), np.log(rr + c)])
            coef, *_ = np.linalg.lstsq(d, yy, rcond=None)
            mse = float(np.mean((d @ coef - yy) ** 2))
            if best is None or mse < best[0]:
                best = (mse, float(coef[0]), float(-coef[1]), float(c))
        return best[1:]

    a1, b1, c1 = zmfit(ranks[:k], logf[:k])
    rr = np.arange(1, len(ranks) - k + 1, dtype=np.float64)
    a2, b2, c2 = zmfit(rr, logf[k:])
    return [a1, min(max(b1, 0.5), 3.0), min(c1, 1000.0), a2, min(max(b2, 0.5), 3.0), min(c2, 1000.0)]


def fit_task(task):
    corpus, family, gate, freqs_list = task
    freqs = np.asarray(freqs_list, dtype=np.float64)
    V = len(freqs)
    ranks = np.arange(1, V + 1, dtype=np.float64)
    log_ranks = np.log(ranks)
    y = np.log(freqs)
    rng = np.random.default_rng(SEED + hash((corpus, family, gate)) % (2**31))
    anchor6 = piecewise_anchor(ranks, y)

    def resid(p):
        return predict(ranks, log_ranks, p, gate) - y

    best_rmse = None
    for i in range(N_STARTS):
        if i < N_STARTS // 2:
            x0 = np.array(anchor6 + [0.0, 0.0, 0.0])
            x0[:6] += rng.normal(0.0, [6.0, 0.3, 90.0, 6.0, 0.3, 90.0])
            x0[6] = 400.0 + rng.normal(0.0, 150.0)
            x0[7] = 0.5 * math.exp(rng.normal(0.0, 0.45))
            x0[8] = 0.5 * math.exp(rng.normal(0.0, 0.45))
        else:
            x0 = np.array([
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                rng.uniform(y.min() - 3, y.max() + 3), rng.uniform(0.5, 3.0), rng.uniform(0, 1000),
                math.exp(rng.uniform(math.log(20), math.log(1000))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
                math.exp(rng.uniform(math.log(0.05), math.log(10))),
            ])
        x0 = np.clip(x0, LOWER + 1e-8, UPPER - 1e-8)
        try:
            sol = least_squares(resid, x0=x0, bounds=(LOWER, UPPER), method="trf", max_nfev=MAX_NFEV)
        except Exception:
            continue
        r = float(np.sqrt(np.mean(sol.fun**2)))
        if best_rmse is None or r < best_rmse:
            best_rmse = r
    bic = 9 * math.log(V) + V * math.log(best_rmse**2)
    return {"corpus": corpus, "family": family, "gate": gate, "V": V, "rmse": best_rmse, "bic": bic}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f1 = {r["corpus"]: r for r in csv.DictReader(open(F1CSV, newline="", encoding="utf-8"))}
    rng = np.random.default_rng(SEED)

    twins, tasks = [], []
    for name, fname in CORPORA:
        freqs = load_freqs(fname)
        ln, gmm = simulate_lognormal_twin(freqs, rng)
        pa = simulate_pareto_twin(freqs, float(f1[name]["alpha_function"]), float(f1[name]["alpha_content"]), rng)
        for fam, fr in [("lognormal_twin", ln), ("pareto_twin", pa)]:
            c, s2 = zm_and_step2(fr)
            twins.append({"corpus": name, "family": fam, "V": len(fr), "tokens": int(fr.sum()),
                          "zm_c": c, "step2_winner": s2, "real_c": float(f1[name]["zm_c"]),
                          "real_winner": f1[name]["canon_family"], **({f"gmm_{k}": v for k, v in gmm.items()} if fam == "lognormal_twin" else {})})
            print(f"twin {name[:24]:25} {fam:15} V={len(fr):6} c={c:8.2f} (real {float(f1[name]['zm_c']):7.2f}) "
                  f"step2={s2} (real {f1[name]['canon_family']})", flush=True)
            for gate in ["logistic", "tanh", "erf", "algebraic", "arctan"]:
                tasks.append((name, fam, gate, fr.tolist()))

    with open(OUT / "f4b_twins.csv", "w", newline="", encoding="utf-8") as f:
        cols = sorted({k for t in twins for k in t}, key=lambda k: (k != "corpus", k))
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(twins)

    rows = []
    with Pool(processes=14) as pool:
        for i, res in enumerate(pool.imap_unordered(fit_task, tasks, chunksize=1), 1):
            rows.append(res)
            print(f"[{i}/{len(tasks)}] {res['corpus'][:22]:23} {res['family']:15} {res['gate']:9} rmse={res['rmse']:.5f}", flush=True)

    with open(OUT / "f4b_fits.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lines = ["# F4b calibrated twins — summary\n"]
    lines.append("| corpus | family | winner | erf margin (BIC) | indep spread |")
    lines.append("|---|---|---|---:|---:|")
    from collections import Counter as C2
    wc = C2()
    for name, _ in CORPORA:
        for fam in ["lognormal_twin", "pareto_twin"]:
            sub = [r for r in rows if r["corpus"] == name and r["family"] == fam and r["gate"] != "tanh"]
            best = min(sub, key=lambda r: r["bic"])
            erf_row = next(r for r in sub if r["gate"] == "erf")
            others = min(r["bic"] for r in sub if r["gate"] != "erf")
            spread = max(r["bic"] for r in sub) - min(r["bic"] for r in sub)
            wc[(fam, best["gate"])] += 1
            lines.append(f"| {name} | {fam} | {best['gate']} | {erf_row['bic'] - others:.1f} | {spread:.1f} |")
    lines.append("\n## Winner counts")
    for (fam, gate), nwin in sorted(wc.items()):
        lines.append(f"- {fam} -> {gate}: {nwin}")
    tw_by = {(t["corpus"], t["family"]): t for t in twins}
    ln_c_match = [(tw_by[(n, 'lognormal_twin')]["zm_c"], tw_by[(n, 'lognormal_twin')]["real_c"]) for n, _ in CORPORA]
    lines.append("\n## Fingerprint checks (lognormal twins)")
    lines.append(f"- corr(twin c, real c): {float(np.corrcoef([a for a, b in ln_c_match], [b for a, b in ln_c_match])[0, 1]):.3f}")
    s2m = sum(1 for n, _ in CORPORA if tw_by[(n, 'lognormal_twin')]["step2_winner"] == tw_by[(n, 'lognormal_twin')]["real_winner"])
    lines.append(f"- step-2 winner matches real corpus: {s2m}/{len(CORPORA)}")
    (OUT / "f4b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nsummary written", flush=True)


if __name__ == "__main__":
    sys.exit(main())
