---
title: Descendant Chain Family Test
type: source
status: reviewed
updated: 2026-11-06
sources:
  - raw/inbox/descendant-chain-validation-20261106.md
  - raw/inbox/descendant-chain-validation-protocol-20261106.md
  - raw/inbox/descendant-chain-validation-summary-20261106.csv
  - raw/inbox/descendant-chain-validation-results-20261106.jsonl
  - raw/inbox/descendant-chain-validation-panels-script.py.txt
  - raw/inbox/descendant-chain-validation-summary-script.py.txt
  - raw/inbox/descendant-chain-benchmark-results-20261106.jsonl
  - raw/inbox/descendant-chain-benchmark-script.py.txt
  - raw/inbox/conditional-small-benchmarks-20261104.py.txt
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py
  - tests/statistics/69_test_diffusion_descendant_chain.py
tags:
  - selective-inference
  - diffusion
  - traversal
---

# Descendant Chain Family Test

## Summary

`selected_diffusion_family_chain_tests` extends the root conditional chain test
to a pre-set family: the root and its internal children. Each node is tested
conditional on its ancestor path, at a fixed Bonferroni level `alpha / 3`. On
known-noise Gaussian panels its null p-values are calibrated; children's power
is low because the diffusion tree splits them as small peels.

## Key Points

- Each node's mean contrast between its children uses the random-position
  reversible chi chain; only the node's rows move. A proposal counts only if
  every split from the root down to the node recurs. At the root this
  reproduces `selected_diffusion_root_chain_test` exactly.
- In 390 pre-set families (30 to 200 rows), 801 null nodes give p at most
  0.05, 0.0167 and 0.01 in 5.37%, 1.75% and 1.12%. Null children under a real
  root split are rejected 1 of 141 times.
- Family-wise false rejections on global nulls are 13/200: consistent with
  0.05, but the one-sided 95% upper bound is 0.1014, so the target is not yet
  established.
- All 140 hierarchical children were split by the observed tree across the true
  subgroups (median smaller side 4 of 20 rows); their rejection rates were 5/100
  and 6/40. With 32 rows and strong signal the tree splits along subgroups and
  every position rejects.
- On the six benchmark cases inside the model, all six single-population nodes
  are kept and 11 of 12 multi-population nodes rejected; the miss is the
  20,000-feature root (p 0.03, 14 chain moves).
- Sequential joblib keeps merge signatures identical and reduces each replay to
  8 to 20 ms for 30 to 200 rows; a family takes 10 to 33 s.

## Evidence

`raw/inbox/descendant-chain-validation-20261106.md` records the design, the
test-first implementation, the protocol fixed before results, all panel
endpoints, the benchmark case selection and results, runtimes and limits. The
per-dataset and per-case JSON lines and the scripts are stored beside it. The
earlier prototype is `descendant_path` in
`raw/inbox/conditional-small-benchmarks-20261104.py.txt`.

The test covers known-noise Gaussian rows, the root and its immediate children
only; counts, fitted PCA, unknown covariance and production traversal with
pass-through remain outside it.

## Links

- [[root-conditional-sampler]]
- [[conditional-diffusion-monte-carlo-evidence]]
- [[choir-cell-lines-recovery-plan]]
