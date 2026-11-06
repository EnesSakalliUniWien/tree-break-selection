---
title: Diffusion Topology Branch Rank Isolation
type: analysis
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/diffusion-factorial-20261102.md
  - raw/assets/diffusion-factorial-20261102/runs.csv
  - raw/assets/diffusion-factorial-20261102/paired_contrasts.csv
  - raw/assets/diffusion-factorial-20261102/nnls_reversal_trace.json
  - benchmarks/diagnostics/calibration/diffusion_factorial.py
  - benchmarks/diagnostics/calibration/diffusion_factorial_trace.py
  - tests/core/test_diffusion_factorial.py
tags:
  - diffusion
  - branch-lengths
  - calibration
---

# Diffusion Topology Branch Rank Isolation

## Summary

The paired factorial experiment identifies a Gaussian-null failure from the
diffusion topology intervention even without NNLS or contrast-energy rank,
and a Bernoulli-null interaction between all three components. It does not
identify a universally calibrated replacement or justify threshold tuning.

## Details

With branch time off and nominal rank, original-distance versus diffusion
topology yields 0 versus 29 false splits in 50 Gaussian-null datasets. This
is a comparison of complete topology constructors at identical downstream
policies. The original-distance roots have a child of size 1–2; diffusion
roots are usually more balanced. Selection and covariance-regime changes
remain pathways within the measured topology effect.

For Bernoulli diffusion, false split counts for (NNLS, energy rank) are
(off,off)=0, (off,on)=4, (on,off)=4, (on,on)=18 out of 50. Rank interventions
hold tree and branch lengths fixed; NNLS interventions hold topology fixed.
The downstream empirical-null model is refitted under the same policy.

In one replay, NNLS reduces raw Q from 109.87886 to 99.59536 but reduces
c_hat from 8.60026 to 5.98049, making adjusted p fall from .01242 to .00226
and changing one cluster into two. Thus branch variance inflation does not
imply a conservative final decision when empirical inflation is refitted.

Diffusion improves mean ARI from .80668 to .97340 on the five-seed
three-category signal case in the topology-only comparison. Five other
relatively easy signal cases are recovered perfectly in every condition.
These recovery results do not remove the null failures.

## Evidence

[[diffusion-factorial-evidence]] links the full protocol and all eight-condition
counts. Dataset-level pairing, complete source hashes and invariant checks are
retained. Four parity tests check NNLS diagnostic results against the existing
runner. Fifty null repetitions cannot establish nominal .001 or .01 tail
control; 0/50 has exact two-sided .95 upper bound .0711.

## Links

- [[diffusion-factorial-evidence]]
- [[adaptive-diffusion-nnls-method-library-audit]]
- [[analytic-selected-root-laws]]
- [[mathematical-glossary]]

## Open Questions

- How much of the Gaussian topology effect comes from selected membership
  versus the changed child-size/covariance regime?
- Can empirical-null fitting account for branch-dependent changes in support
  and selected-record distributions without the observed significance reversals?
- What selected-region model covers the adaptive diffusion constructor and
  contrast-energy rank jointly? The restricted squared-Euclidean root laws do not.
