---
title: Diffusion Factorial Evidence
type: source
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/diffusion-factorial-20261102.md
  - raw/assets/diffusion-factorial-20261102/manifest.json
  - raw/assets/diffusion-factorial-20261102/summary.csv
  - raw/assets/diffusion-factorial-20261102/nnls_reversal_trace.json
tags:
  - diffusion
  - calibration
  - benchmarks
---

# Diffusion Factorial Evidence

## Summary

A paired eight-condition experiment separates topology, NNLS branch-time use
and contrast-energy rank in 180 datasets, producing 1,440 pipeline runs.

## Key Points

- Explicit iid nulls use 50 datasets each for Gaussian, Bernoulli and categorical
  observations; six existing signal cases use five seeds each.
- Gaussian diffusion false splits occur without NNLS or adaptive rank.
- Bernoulli diffusion, NNLS and energy rank interact; a traced NNLS reversal
  occurs because fitted inflation falls more than the raw Wald statistic.
- Original-distance Gaussian null trees have root children of size 1–2, so
  zero final splits there does not demonstrate calibration.
- Compact raw outcomes, root annotations, provenance and executable scripts
  preserve results; no production method or threshold was changed.

## Evidence

`raw/inbox/diffusion-factorial-20261102.md` specifies the design, null models,
reproduction commands, paired findings, uncertainty and unresolved mechanisms.
The cited manifest identifies the actual working-tree sources and environment.

## Links

- [[diffusion-topology-branch-rank-isolation]]
- [[analytic-selected-root-laws]]
