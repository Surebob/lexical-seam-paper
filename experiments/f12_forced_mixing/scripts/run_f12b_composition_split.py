"""
Adversarial test of the Williams (2015) text-mixing hypothesis against the s-law.

Williams argues the two-regime structure is a MIXING artifact: heavily aggregated
texts should behave differently from single continuous works. If the seam WIDTH
s = k * w_tail were a mixing artifact, s/V should differ between:
  - AGGREGATES (collections of many independent sub-texts: plays, tales, books, essays)
  - SINGLE continuous works (one author, one continuous narrative/argument)

We test this on the existing 25-EN panel (no new fitting; uses f2_per_corpus.csv).
"""
import csv
import numpy as np
from pathlib import Path

repo = Path(__file__).resolve().parents[3]
with open(repo / "experiments/f2_k_profile_likelihood/outputs/f2_per_corpus.csv", newline="", encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh))

# Classify by degree of aggregation (Williams' notion of "text mixing").
AGGREGATE = {
    "Complete Works of Shakespeare",  # ~37 plays + poems
    "Complete Sherlock Holmes",       # 4 novels + 56 stories
    "Grimm's Fairy Tales",            # ~200 tales
    "Aesop's Fables",                 # hundreds of fables
    "Canterbury Tales",               # 24 tales
    "Arabian Nights (Vol 1)",         # many tales
    "Federalist Papers",              # 85 essays, 3 authors
    "King James Bible",               # 66 books, many authors, translation
    "Collected Poe",                  # many poems/stories
    "Dubliners",                      # 15 short stories
}
# Everything else = single continuous work (one author, one narrative/argument).
recs = []
for r in rows:
    corpus = r["corpus"]
    V = float(r["V"])
    s = float(r["s_at_min"])
    kind = "aggregate" if corpus in AGGREGATE else "single"
    recs.append({"corpus": corpus, "V": V, "s": s, "s_over_V": s / V, "kind": kind})

out = []
def w(x): out.append(x)

w("Per-corpus s/V (sorted):")
for r in sorted(recs, key=lambda r: r["s_over_V"]):
    w(f"  {r['kind']:9s}  s/V={r['s_over_V']:.5f}  V={int(r['V']):6d}  s={r['s']:7.1f}  {r['corpus']}")

w("")
for kind in ("single", "aggregate"):
    v = np.array([r["s_over_V"] for r in recs if r["kind"] == kind])
    w(f"{kind:9s} n={len(v):2d}  mean s/V={v.mean():.5f}  sd={v.std(ddof=1):.5f}  "
      f"median={np.median(v):.5f}  range=[{v.min():.5f},{v.max():.5f}]")

agg = np.array([r["s_over_V"] for r in recs if r["kind"] == "aggregate"])
sng = np.array([r["s_over_V"] for r in recs if r["kind"] == "single"])

# Welch t-test (manual, no scipy dependency needed but use it if present)
try:
    from scipy import stats
    t, p = stats.ttest_ind(agg, sng, equal_var=False)
    u, pu = stats.mannwhitneyu(agg, sng, alternative="two-sided")
    w("")
    w(f"Welch t-test  single vs aggregate:  t={t:.3f}  p={p:.3f}")
    w(f"Mann-Whitney U (rank, distribution-free): U={u:.1f}  p={pu:.3f}")
except Exception as e:
    w(f"(scipy stats unavailable: {e})")

# Bootstrap difference in means (distribution-free), 10000 reps, seeded.
rng = np.random.default_rng(20260722)
diffs = []
for _ in range(10000):
    a = rng.choice(agg, size=len(agg), replace=True)
    s = rng.choice(sng, size=len(sng), replace=True)
    diffs.append(a.mean() - s.mean())
diffs = np.array(diffs)
lo, hi = np.percentile(diffs, [2.5, 97.5])
w("")
w(f"Bootstrap mean(aggregate) - mean(single) = {agg.mean()-sng.mean():+.5f}")
w(f"  95% CI [{lo:+.5f}, {hi:+.5f}]  (contains 0 => indistinguishable)")

# Also: does s/V correlate with V within each group? (Williams would predict a
# composition-driven trend; the s-law predicts flat s/V independent of V.)
for kind in ("single", "aggregate"):
    Vv = np.array([r["V"] for r in recs if r["kind"] == kind])
    sv = np.array([r["s_over_V"] for r in recs if r["kind"] == kind])
    if len(Vv) > 3:
        r = np.corrcoef(np.log(Vv), np.log(sv))[0, 1]
        w(f"  corr(log V, log s/V) [{kind}] = {r:+.3f}  (near 0 => s/V flat in V)")

(Path(__file__).resolve().parent.parent / "outputs" / "f12b_composition_split.txt").write_text("\n".join(out), encoding="utf-8")
print("wrote f12b_composition_split.txt")
