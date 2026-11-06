---
title: Conditional Diffusion Monte Carlo Evidence
type: source
status: reviewed
updated: 2026-11-04
sources:
  - raw/inbox/conditional-diffusion-statistics-20261104.md
  - raw/inbox/conditional-diffusion-wolfram-20261104.json
  - raw/inbox/conditional-diffusion-experiment-20261104.py.txt
  - raw/inbox/conditional-diffusion-mc-initial-20261104.py.txt
  - raw/assets/conditional-diffusion-20261104/summary.json
  - raw/assets/conditional-diffusion-20261104/model-control.json
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py
  - tests/statistics/67_test_diffusion_root_conditional_mc.py
  - raw/inbox/conditional-small-benchmarks-review-20261104.md
  - raw/inbox/conditional-small-benchmarks-20261104.py.txt
  - raw/inbox/conditional-small-certificates-20261104.py.txt
  - raw/assets/conditional-small-benchmarks-20261104/model-audit.json
  - raw/assets/conditional-small-benchmarks-20261104/descendant-and-spectrum-audit.json
  - raw/assets/conditional-small-benchmarks-20261104/summary.json
  - reports/selective_rebuild_comparison_20261104/README.md
tags:
  - diffusion
  - selective-inference
  - gaussian
---

# Conditional Diffusion Monte Carlo Evidence

## Summary

The original known-covariance Gaussian root diagnostic conditions on the complete native
diffusion average-linkage merge history. Wolfram checked its chi density and
tail identities. A fixed-proposal conditional rank test completed all 100
nine-row Gaussian cases after handling strict history exclusion before later
near-ties. It had zero rejections under both null and signal scenarios.

The later [[selective-rebuild-comparison]] changes the sampling computation:
random-position reversible-chain ranks detect three of six strong examples
while preserving full history, and five while preserving the root split at
0.05. The ordinary sampler accepts none in these same six examples. These
restricted numerical results do not establish cancer-data or descendant validity.

## Key Points

- The statistic uses full input dimension, known identity covariance and
  frozen nuisance offset and contrast direction. Every proposal reconstructs
  the native diffusion geometry and hierarchy.
- Conditional on the number K accepted from a fixed proposal budget, the
  null observation and accepted radii are exchangeable. The rank p-value
  includes the observation and has minimum 1/(K+1).
- At alpha 0.01 the nominal chi tails rejected 20/60 homogeneous nulls,
  17/20 moderate signals and 20/20 strong signals. The conditional rank test
  rejected none. Each strong signal had zero accepted radii among 999 proposals.
- The 0/60 null interval extends to 0.05963 at 95% confidence. The finite
  panel cannot establish the 0.01 target. Its principal practical finding is
  inadequate sampling resolution for strongly selected histories.
- Whole-history conditioning is stronger than the root-partition event.
  The numerical replay is separate from exact-real membership certification.
  Cancer-data covariance, prior preprocessing and traversal remain outside
  the tested model.

## Evidence

### Existing small benchmark check

The registered 60-by-40 clear and 45-by-30 moderate continuous Gaussian cases
preserve the raw continuous inputs and have known spherical noise scales.
Their root numerical p-values are both 0.001 with 999 accepted proposals.
Global-tree descendant replays tested the two immediate internal children.
The clear mixed 45-row child rejects, both pure 15-row children remain unsplit,
and the moderate mixed 30-row child has zero accepted draws and p=1. The family
is a fixed root-plus-children family conditional on the full history, using
Bonferroni thresholds. Production adaptive traversal remains outside this check.

The registered 200-by-80 continuous null root has nominal p=9.41e-13 and
conditional rank p=1 from zero accepted proposals. Its replay took 403 seconds.
All seven root/child contrasts completed and their ranks were independently
recomputed. These fixed cases do not provide an empirical error-rate bound.

The exact-real certifier gives [0,1] for both signal roots. Their graphs have
four and three connected components, yielding repeated leading eigenvalues
outside its isolated-leading-mode assumptions. A scalar control certifies a
small observed interval but retains almost [0,1] bounds overall. The certificate
code reconstructs the root contrast internally, so it cannot certify the tested
descendant paths. Seven other generated benchmark cases and six count-model
specifications were audited for incompatible observation-model assumptions.

### Original Gaussian panel

The raw derivation includes the exchangeability argument and probability
control. The JSON captures preserve Wolfram input/output and 20,000 exact-region
replicates. The native panel retains all datasets, accepted radii, random seeds,
source hashes and runtimes. Its initial unresolved cases are archived alongside
the paired strict-exclusion rerun. Tests cover the native failure, rank
accounting, unresolved membership, empty acceptance and reproducibility.

## Links

- [[root-conditional-sampler]]
- [[selective-rebuild-comparison]]
- [[diffusion-selected-root-diagnostic]]
- [[analytic-selected-root-laws]]
- [[choir-cell-lines-recovery-plan]]
