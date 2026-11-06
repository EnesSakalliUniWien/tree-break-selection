---
title: Adaptive Alpha Tree Procedures
type: source
status: reviewed
updated: 2026-11-06
sources:
  - raw/inbox/adaptive-alpha-tree-20261106.md
  - raw/inbox/adaptive-alpha-tree-protocol-20261106.md
  - raw/inbox/adaptive-alpha-tree-libraries-20261106.md
  - raw/inbox/adaptive-alpha-tree-results-20261106.jsonl
  - raw/inbox/adaptive-alpha-tree-procedures-summary-20261106.csv
  - raw/inbox/adaptive-alpha-tree-assumptions-20261106.json
  - raw/inbox/adaptive-alpha-tree-panels-script.py.txt
  - raw/inbox/adaptive-alpha-tree-procedures-script.py.txt
  - raw/inbox/adaptive-alpha-tree-summary-script.py.txt
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py
tags:
  - alpha
  - adaptive
  - selective-inference
  - multiple-testing
---

# Adaptive Alpha Tree Procedures

## Summary

Three ways to let the chain-test alphas adapt to the opening tree were compared
with fixed Bonferroni on 650 known-noise Gaussian families of depth at most 2.
Alpha-investing rejected 30% more multi-population nodes at mFDR; alphas fixed
by each node's event keep a provable FWER but gain little; weighting by the
diffusion manifold has no validity argument and gains nothing. AdaPT and STAR
cannot reject anything in families this small.

## Key Points

- Alpha-investing (`W0 = B0 = 0.05`, half the wealth per test, breadth-first,
  parent first) rejected 313 of 780 mixed nodes against 240 for the current
  `alpha / 3` family test, with 3/150 global-null families falsely rejecting.
  Its guarantee is mFDR and needs null p-values super-uniform given earlier
  decisions.
- That tree assumption held at the decision level: null children were rejected
  3/92 and 6/190 times after a rejected null or mixed parent. Parent and child
  p-values are negatively correlated (Spearman -0.14 and -0.12), the
  conservative direction, so independence fails while super-uniformity was not
  contradicted.
- Alphas fixed by each node's event, `alpha * c_d * n_v / n` and a learned
  depth/size allocation, keep FWER at most alpha on any tree by the union bound.
  Null p-values are unrelated to depth, size share and child balance (all
  p >= 0.80). They gained 3 mixed nodes over the current test; learning added
  nothing.
- Diffusion merge height, which no event fixes, is related to null p-values
  (rho -0.061, p 0.001), so manifold-weighted alphas lack a validity argument;
  in these panels they matched the current test.
- AdaPT's `(1 + A) / R` and STAR's `(h(1) + sum h) / (1 + |R|)` FDP estimates
  need at least 20 and 19 rejections; a family has at most 7 nodes. AdaPT is
  archived from CRAN; STAR is unlicensed 2019 paper code; neither has a Python
  port. `online-fdr` (Python, BSD-3) is the reference for alpha-investing.
- No alpha rescues peel-shaped child splits: depth-1 power in the 2 x 2 panel
  was 3 to 11 of 100 for every procedure.

## Evidence

`raw/inbox/adaptive-alpha-tree-20261106.md` records the design, endpoints,
assumption tests and limits; the protocol was fixed before any panel ran.
`raw/inbox/adaptive-alpha-tree-libraries-20261106.md` records the package
checks and the update rules and assumptions read from the alpha-investing,
AdaPT and STAR papers. Per-family node rows, the endpoint table, the
assumption-test JSON and the three scripts sit beside them. Evidence covers
known-noise Gaussian rows with 30 to 80 rows only.

## Links

- [[descendant-chain-family-test]]
- [[root-conditional-sampler]]
- [[adaptive-alpha-gate-design-20260628]]
- [[alpha-structure-sweep-20260627]]
