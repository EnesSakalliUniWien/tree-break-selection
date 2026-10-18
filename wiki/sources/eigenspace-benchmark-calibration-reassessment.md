---
title: Eigenspace Benchmark and Calibration Reassessment
type: source
status: reviewed
updated: 2026-09-28
sources:
  - reports/eigenspace_reassessment_20260928/README.md
  - reports/eigenspace_reassessment_20260928/spectral_blocks/comparison_summary.json
  - reports/eigenspace_reassessment_20260928/wald_rank/counts.csv
  - reports/eigenspace_reassessment_20260928/wald_rank/verification.json
  - reports/eigenspace_reassessment_20260928/full_pipeline/protocol.md
  - reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/manifest.json
  - reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/summary.csv
  - reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/root_tails.csv
  - reports/eigenspace_reassessment_20260928/full_pipeline/evaluation/verification.json
  - reports/eigenspace_reassessment_20260928/signal_recovery_summary.csv
tags:
  - eigenvalues
  - benchmark
  - calibration
---

# Eigenspace Benchmark and Calibration Reassessment

## Summary

The complete-eigenspace rank repair removes the reproduced basis dependence
but does not calibrate contrast-selected rank across distinct eigenvalue groups
or the complete empirical pipeline. Matched benchmark recovery improves and
worsens in different runs. Independent null datasets show excess false
partitions under several supported feature/method combinations.

## Key Points

- The same 11-case, 47-block, eight-tree-method, two-geometry benchmark has
  752 attempts in both versions: 525 successful, 170 failed trees and 57
  unsupported. Eighteen of 582 available trees and 33 of 525 partitions change.
  Twenty-two partition changes occur on unchanged trees. ARI improves in 15
  runs and worsens in 15; single-cluster outputs increase from 436 to 459.
- All raw geometries and 44/47 diffusion geometries are exact. The three
  changed diffusion geometries have disconnected leading eigenspaces of
  dimensions 4, 2 and 6. The below-BBP one-population example changes from
  10/36 multicluster results to 0/36, but these are repeated analyses of one
  dataset, not independent calibration trials.
- The isolated Wald experiment uses six million independent Gaussian
  contrasts, four million held out. At alpha .01, repaired adaptive energy
  rank .90 gives held-out rates .009949/.009993 for `(1,1)`, .014464/.014595
  for `(2,1)` and .015531/.015342 for `(4,2,2,1)`. All 16 held-out current
  distinct/partially tied adaptive cells exceed their targets after Holm
  adjustment across 216 cells. Rotation invariance now holds in the checks.
- A frozen current-implementation pipeline study uses 500 null and 200 signal
  datasets in each of six families, five registered methods and prescribed iid
  Bernoulli alpha sensitivity. Its 4,200 independent datasets yield 23,800
  method/alpha attempts: 19,595 successful, 4,201 skips, four unsupported, zero
  execution errors and zero audit errors. Skips include native continuous data
  with registered Hamming diffusion; they are never counted as null successes.
- At sibling alpha .01, adaptive diffusion NNLS gives false partitions on
  316/500 iid Bernoulli null datasets, Wilson interval [.589,.673], and
  263/500 categorical null datasets, interval [.482,.569]. Finite selected-root
  null p-values also have excess lower-tail mass: 280/482 and 240/481 are at
  most .01. The generic full-pipeline methods disable adaptive energy rank,
  so the isolated optional rank defect does not explain away these failures.
- Adaptive NNLS mean signal ARI is .877 for iid Bernoulli and .649 for
  categorical data; these recovery results coexist with null failures.
  Generator-nonnull calibration records appear in 118/200 and 194/200 signal
  runs. Base support weights are recorded; focal kernel-weight contamination
  and its causal effect are not established.
- High-dimensional Gaussian data have zero observed false partitions but
  also return one cluster on all 200 planted-signal datasets for every
  supported method, with ARI zero. Low error alone is not successful recovery.

## Evidence

The report records matched seeds, inputs, source hashes and settings. Frozen
source hashes remained stable throughout both executions. The independent
pipeline verifier regenerates every input, checks 1,184,337 node records,
reconstructs completed-run aggregates and 4,140 summary/tail values, and passes
39 deterministic production reruns. Per-edge rank aggregates are verified on
the reruns; complete per-edge raw ranks were not retained. The Wald verifier
decodes all saved rejection counts and checks 2,688 actual-kernel comparisons.

Global-null partitions, traversal-used generator-null split FDP, selected-root
tails and signal recovery are distinct endpoints. Null truth uses exact
generator-mixture counts, and uncertainty uses independent datasets rather
than overlapping nodes. Holm families are fixed at 34 global-null partition
cells and 272 root-tail cells, retaining unavailable cells. These finite
experiments do not establish a general selected-tail or whole-tree guarantee.

The full-pipeline study evaluates the current implementation only. The spectral
benchmark and isolated rank experiment supply the before/after comparisons.
No production calibration change, commit or push was made.

## Links

- [[eigenvector-calibration-construction-review]]
- [[eigenvalue-block-tree-comparison]]
- [[diffusion-nnls-calibration-open-points]]
- [[empirical-null-calibration-reference-law-contract]]
- [[current-adaptive-diffusion-subspace-tree-pipeline]]
