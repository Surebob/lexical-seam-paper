# Experiment 9 Soft-k Comparator Diagnostic

## Question

The `zipf_hybrid_vs_softk_analysis` bundle reports a `16/25` hybrid-vs-soft-k result, while `zipf_hybrid_headtail_splitfit` reports `22/25`. Before selecting a canonical source, we traced the analysis bundle's `softk_full_nll` values against every soft-k-bearing source in the bundle forest.

## Result

The analysis bundle's soft-k values differ from `zipf_hybrid_headtail_splitfit` on all `25/25` rows, not just a small subset.

They match exactly, on all `25/25` rows, the legacy pre-splitfit derivative branch:

- `results/zipf_hybrid_headtail/hybrid_headtail_table.csv`, column `full_softk`
- `results/zipf_angle1_head_windows/head_window_table.csv`, `cutoff=full`, column `softk`
- `results/zipf_angle3_asymmetric_gate/asymmetric_gate_table.csv`, column `full_softk`
- the corresponding `summary.json` full held-out soft-k values in those same legacy bundles

They match `0/25` rows for the current or patched soft-k sources:

- `results/zipf_v4_verification/table_a_fourway_pmf.csv`, column `softk_test_avg_nll`
- `results/zipf_seam_mandelbrot_softk_splitfit/softk_table.csv`, column `softk_test_avg_nll`
- `results/zipf_seam_mandelbrot_softkw/softkw_table.csv`, repeated `softk_test_avg_nll`
- `results/zipf_hybrid_headtail_splitfit/hybrid_headtail_table.csv`, column `full_softk`
- `results/zipf_angle1_head_windows_splitfit/head_window_table.csv`, `cutoff=full`, column `softk`
- `results/zipf_angle3_asymmetric_gate_splitfit/asymmetric_gate_table.csv`, column `full_softk`

## Mechanism

This is not false provenance. It is a known legacy-comparator branch.

The producer script [`/Volumes/External2TB/emlexperiment/zipf_hybrid_vs_softk_analysis.py`](/Volumes/External2TB/emlexperiment/zipf_hybrid_vs_softk_analysis.py) loads:

- `HYBRID_SUMMARY_PATH = results/zipf_hybrid_headtail/summary.json`
- `SOFTK_SUMMARY_PATH = results/zipf_seam_mandelbrot_softk/summary.json`

But for the actual `softk_full_nll` column it uses:

```python
softk_full = float(h["heldout"]["full"]["avg_nll"]["softk"])
```

where `h` is the legacy `zipf_hybrid_headtail` row. The separate soft-k summary is used only for metadata such as lambda, transition fraction, and step-2 help. Therefore the analysis bundle is explicitly analyzing the legacy hybrid-headtail comparator, not the splitfit or v4 canonical soft-k comparator.

This matches the existing audit docs: `results/audit_evaluation_protocols.md` classifies `zipf_hybrid_headtail` as `incorrect-conflation` and `zipf_hybrid_vs_softk_analysis` as `different-protocol-but-labeled-clearly`.

## Interpretation for consolidation

The `16/25` result is legitimate as a legacy pre-splitfit comparator analysis. It should be preserved as research record and may be useful for Paper 2 history.

It should **not** be used as the current canonical Paper 1/Paper 2 headline hybrid-vs-soft-k count without explicitly choosing to revive the legacy comparator convention. The current canonical post-cleanup source remains `zipf_hybrid_headtail_splitfit`, yielding:

- hybrid beats soft-k: `22/25`
- median hybrid minus soft-k: `-0.001677944396615949`

## Detailed row-level diagnostic

See [`softk_comparator_source_diagnostic.csv`](/Volumes/External2TB/emlexperiment/experiments/9_likelihood_frontier/softk_comparator_source_diagnostic.csv).
