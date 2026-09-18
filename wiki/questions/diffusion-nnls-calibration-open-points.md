---
title: Diffusion NNLS Calibration Open Points
type: question
status: draft
updated: 2026-09-10
sources:
  - reports/calibration_restoration_20260909/false_split_definition.md
  - reports/calibration_restoration_20260909/validation_protocol.md
  - reports/calibration_restoration_20260909/README.md
  - reports/calibration_restoration_20260909/external_contract_check.json
  - reports/diffusion_nnls_versions_20260909/README.md
  - reports/diffusion_nnls_versions_20260909/support_review_results.csv
tags:
  - calibration
  - diffusion
  - nnls
  - validation
---

# Diffusion NNLS Calibration Open Points

## Question

What must be defined or established to validate the restored diffusion NNLS
calibration while preserving its signal recovery?

## Current State

The empirical calibration rule is restored and benchmark recovery is recorded
in [[calibration-rule-restoration-20260909]]. The held-out full-pipeline
calibration study has not been run. These open points distinguish study-design
decisions from findings that require evidence. An item closes only when its
answer, applicable conditions and supporting artifact are recorded.

- [x] **CAL-01 — Null models and truth: definition recorded.** A local false
  split is a sibling rejection used by traversal while the generator-defined
  population parameters of its children are equal. The specification separates
  this from a final false partition under the global null and fragmentation on
  signal data, defines mixed-child truth, and lists required data conditions.
  See [false-split definition](../../reports/calibration_restoration_20260909/false_split_definition.md).
  This closes the definition; generator implementation and calibration validation
  remain pending under the subsequent points.

- [ ] **CAL-02 — Error-rate target.** Which guarantee is sought: local sibling
  Type I error, false root splits, any false split under a one-cluster null, or
  FDR among selected splits? Define the tested family, numerator, denominator,
  conditioning and target alpha for each reported quantity. Separate
  conditional-on-eligible outcomes from whole-dataset outcomes.

- [ ] **CAL-03 — Data-dependent selection.** Which selection events affect the
  focal null distribution, and what must be reproduced or conditioned on?
  Specify diffusion fitting, tree construction, NNLS branch lengths, covariance,
  PCA/rank selection, edge selection, calibration-record selection and traversal.
  Rebuild the complete method in each null replicate; use controlled stage
  ablations to identify where the reference distribution changes.

- [ ] **CAL-04 — Calibration-record contamination and dependence.** Do
  apparently null-like or edge-blocked records contain signal, and how does their
  weighting affect the fitted scale? Measure contamination using generator
  truth, observation overlap, stopping-event dependency groups and concentration
  of calibration weights. Establish the independent sampling unit rather than
  treating overlapping nodes as independent calibration observations.

- [ ] **CAL-05 — Adequacy of empirical inflation.** Does the restored
  `chi2.sf(T / (reference_scale * c_hat), df)` produce valid null tails after
  selection, including uncertainty in the estimated `c_hat`? Evaluate whether
  `Pr(p <= t | declared null/selection context) <= t` across relevant tails and
  contexts. Distinguish conservative nonuniformity from excessive rejection;
  matching the mean or failing a uniformity test does not settle this question.

- [ ] **CAL-06 — Necessary support and the unsupported rule.** Which support
  conditions are actually necessary for reliable empirical calibration? Measure
  tail error and signal recovery against group count, effective support,
  weight concentration and leave-one-group sensitivity. Examine zero-support
  cases separately: useful cluster assignments and validated p-values are
  different outcomes. Establish an evidence-based output policy before changing
  the currently retained no-support rule.

- [ ] **CAL-07 — BH and traversal.** Does traversal-aligned sibling BH control
  the intended error rate with overlapping nodes, edge selection and pass-through
  traversal? Evaluate root and descendant decisions and final partitions under
  both global null and mixed null/signal configurations. Identify the assumptions
  needed for a formal guarantee; marginal p-value checks alone do not establish
  method-level FDR.

- [ ] **CAL-08 — Independent held-out evaluation.** How will development and
  evaluation datasets/seeds be separated, and which choices are frozen before
  evaluation? Record disjoint replicate identifiers, generator settings, source
  hashes, geometry, solver settings, weights, bandwidths and support thresholds.
  Any correction chosen using evaluation outcomes needs a fresh held-out test.

- [ ] **CAL-09 — Simulation precision.** How many independent datasets and
  eligible selected contexts are needed at each alpha? Prespecify confidence
  bounds, tolerances, tail resolution and stopping rules. Count matching
  independent simulations, not just selected records; report rare contexts and
  unavailable outcomes. The existing 499-matching-simulation resolution floor
  does not by itself determine the sample size for the desired error-rate precision.

- [ ] **CAL-10 — Signal recovery and the source of improvement.** Which part of
  diffusion, NNLS and calibration produces the observed gain, and what power is
  lost by a proposed correction? Use matched cases/seeds and controlled ablations
  to compare the restored rule with uninflated testing and any candidate correction.
  Measure ARI, recovered K, over-splits and under-splits alongside false splits;
  closing every gate cannot count as successful calibration and recovery.

- [ ] **CAL-11 — Sensitivity and transfer.** Over which alpha values, sample
  sizes, feature dimensions, feature families, noise levels and cluster balances
  do the findings hold? Check bandwidth/support choices and diffusion NNLS
  versions separately. Keep Hamming-compatible representations distinct from
  native continuous data, and state the conditions supported by evidence rather
  than transferring one preset's result to every method version.

- [ ] **CAL-12 — Acceptance and the next correction.** What prespecified
  evidence is sufficient to retain the empirical rule, restrict its scope or
  change calibration? Record error and power criteria before evaluating candidates.
  If the restored tail exceeds its target, determine whether a selected-tail
  correction or a narrower selection-law derivation resolves the failure and
  verify the result on fresh held-out data. State separately what is empirically
  validated and what has a mathematical guarantee.

CAL-01, CAL-02, CAL-03, CAL-08, CAL-09 and the acceptance criteria in CAL-12 must
be specified before the confirmatory simulations. Their outcomes then address
CAL-04 through CAL-07, CAL-10, CAL-11 and the correction decision in CAL-12.

## Evidence

- `reports/calibration_restoration_20260909/validation_protocol.md` defines the
  complete-pipeline replay, held-out evaluation, tail requirements and sampling
  limitations from which this checklist is derived.
- `reports/calibration_restoration_20260909/README.md` records restored benchmark
  behavior and distinguishes it from a selected-tail calibration guarantee.
- `reports/calibration_restoration_20260909/external_contract_check.json` records
  the existing external diagnostic's insufficient matching-simulation counts.
- `reports/diffusion_nnls_versions_20260909/support_review_results.csv` records
  matched support/gate/alpha comparisons showing why support, clustering quality
  and false splits require separate examination.

## Links

- [[calibration-rule-restoration-20260909]]
- [[empirical-null-calibration-reference-law-contract]]
- [[diffusion-nnls-versions-and-support-review-20260909]]
- [[open-mathematical-questions]]
