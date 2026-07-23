# Experiment 4 README

## Experiment id

`4_analytic_seam_theory`

## Research question

What sign pattern does the smooth seam predict analytically, and how far can tangent-space projection and second-order correction carry that prediction toward the empirical residual?

## Current canonical status

**BLOCKED for canonical decoupled-erf output.**

The current smooth-model canon is the decoupled-erf model. The historical seam-theory bundles available in `results/` are coupled-logistic derivations from the project's prior canonical model. They are legitimate original experimental outputs, not drift or errors, but they do not satisfy the current decoupled-erf canonical output contract.

Per the approved Option C migration decision, the canonical experiment records this blocked state and preserves the historical coupled-logistic results as first-class archived research output.

## Why this is one experiment

The three historical seam-theory bundles are progressive layers of one analytic question:

- `zipf_seam_sign_theory`: raw/local seam sign law
- `zipf_seam_projection_theory`: first-order tangent-space projection
- `zipf_seam_second_order_theory`: second-order correction and quadratic surrogate

They are not duplicates; together they form the analytic seam-theory sequence.

## Upstream dependencies

Canonical decoupled-erf dependencies:

- Missing: decoupled-erf seam sign-law outputs
- Missing: decoupled-erf tangent-space projection outputs
- Missing: decoupled-erf second-order correction outputs

Historical legacy dependencies:

- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_sign_theory/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_seam_sign_theory/summary.json)
- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_projection_theory/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_seam_projection_theory/summary.json)
- [`/Volumes/External2TB/emlexperiment/results/zipf_seam_second_order_theory/summary.json`](/Volumes/External2TB/emlexperiment/results/zipf_seam_second_order_theory/summary.json)

## Outputs produced

### Canonical

The canonical output area contains only:

- `BLOCKED.md`
- `outputs/manifest.json`

No canonical `outputs/seam_sign_checks.csv` or canonical `outputs/aggregate_statistics.csv` is emitted because the decoupled-erf seam theory has not been computed.

### Historical archive

The coupled-logistic prior-canon output is preserved under:

- `archive/legacy_coupled_logistic/README.md`
- `archive/legacy_coupled_logistic/scripts/run_legacy_coupled_logistic.py`
- `archive/legacy_coupled_logistic/outputs/seam_sign_checks.csv`
- `archive/legacy_coupled_logistic/outputs/aggregate_statistics.csv`
- `archive/legacy_coupled_logistic/outputs/source_consistency_checks.csv`
- `archive/legacy_coupled_logistic/outputs/manifest.json`

## Manuscript lines fed

Under the old coupled-logistic manuscript contract, this experiment fed Section 3.8 and discussion restatements of the seam-theory sign counts.

Under the current decoupled-erf canon, those empirical seam-theory validation claims are not satisfied by canonical outputs. The analytic equations themselves may still be manuscript-internal mathematical derivation, but any empirical sign/projection counts require either:

- a fresh decoupled-erf seam-theory computation, or
- explicit manuscript framing as historical coupled-logistic prior-canon evidence.

## Methods

Canonical method status:

- No fresh computation is run in Phase 2 migration.
- The decoupled-erf seam-theory experiment is recorded as blocked because source outputs are absent.

Legacy method status:

- The archive producer reads the three historical coupled-logistic `summary.json` files.
- It joins per-corpus rows by `slug`.
- It emits one consolidated per-corpus table and one aggregate table matching the old claim-map metric names.
- Repeated counts across bundles are written to `source_consistency_checks.csv` to make agreement/disagreement auditable.

## Rerun command

```bash
python3 experiments/4_analytic_seam_theory/run_experiment.py
```

This command records the canonical blocked state and regenerates the coupled-logistic archive tables. It does not run a decoupled-erf seam derivation.

## Verification mapping

Canonical verification:

- `BLOCKED.md` must exist.
- `outputs/manifest.json` must report `canonical_status = blocked_missing_decoupled_erf_seam_theory`.
- No canonical seam-sign CSV should be present unless a decoupled-erf source is later supplied.

Legacy verification:

- `archive/legacy_coupled_logistic/outputs/seam_sign_checks.csv` must contain 25 rows.
- `archive/legacy_coupled_logistic/outputs/aggregate_statistics.csv` must contain the historical coupled-logistic seam-theory counts.
- `archive/legacy_coupled_logistic/outputs/source_consistency_checks.csv` must show whether duplicated metrics agree across source bundles.

## AUDIT

Archive means "superseded by later canonical choice," not "low quality." The coupled-logistic seam-theory outputs are preserved because they were real research outputs under the previous model canon and remain useful as model-canon-comparison reference material.

The current manuscript/claim-map entries that point to experiment 4 need revision if the paper continues to use decoupled-erf as its canonical smooth model.

