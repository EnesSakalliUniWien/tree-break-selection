---
title: Known-Branch Skew Diagnostic
type: source
status: reviewed
updated: 2026-09-28
sources:
  - raw/inbox/recursive-method-program-20260604.txt
  - benchmarks/diagnostics/calibration/sibling/known_branch_skew_diagnostic.py
  - tests/validation/calibration/195_test_known_branch_skew_diagnostic.py
  - reports/known_branch_skew_20260928/protocol.md
  - reports/known_branch_skew_20260928/README.md
  - reports/known_branch_skew_20260928/evaluation/manifest.json
  - reports/known_branch_skew_20260928/verification.json
  - reports/known_branch_skew_20260928/paired_stage_summary.csv
  - reports/known_branch_skew_20260928/node_recovery_summary.csv
  - reports/passthrough_benchmark_20260928/pre_run_manifest.json
  - reports/passthrough_benchmark_20260928/evaluation/manifest.json
  - reports/passthrough_benchmark_20260928/independent_verification.json
tags:
  - calibration
  - diagnostics
  - traversal
---

# Known-Branch Skew Diagnostic

## Summary

A bounded implementation of the June 4 plan separates covariance,
same-data PCA/MP, selected-root inference, and FDR/traversal effects using
known null and signal sibling contrasts. The 10,000-dataset Gaussian study
finds excess null rejection under same-data projection and especially
same-data root selection. It also demonstrates how a null ancestor can block
a strong descendant despite a successful local test.

## Key Points

- The primary plan is Q10, Q17, Q19 and Q24–Q25 in
  `raw/inbox/recursive-method-program-20260604.txt`. Q19's ablation structure
  is implemented with Gaussian data here; its categorical study remains open.
- Five patterns on `R -> (A,B)`, `A -> (A0,A1)`, `B -> (B0,B1)` use declared
  population means and exact null truth. Each of two covariance profiles has
  1,000 independent datasets per pattern, each with four groups of 32
  eight-dimensional observations. The 50-dataset pilot is excluded.
- At alpha .01, global-null fixed-root known-covariance/full-coordinate
  rejection counts are 12/1000 and 15/1000 (identity/AR). Same-data PCA/MP
  yields 33/1000 and 37/1000. Same-data average-linkage selected roots yield
  765/1000 and 762/1000 even with known covariance and full coordinates.
- Sample pooled covariance includes between-child mean variation. Its
  PCA/MP variant keeps all eight dimensions because whitening produces a
  tied spectrum, and matches full-coordinate rejection decisions. This is
  not evidence for validity of dimension-reducing same-data projection.
- Raw-statistic skew is compared only among null rows against the appropriate
  chi-square df-mixture reference. Changing degrees of freedom changes the
  reference skew; null p-value tails remain a separate calibration endpoint.
- With deliberately open edge prerequisites, any reached false split occurs
  in 12/1000 and 15/1000 global-null datasets without pass-through, and 20/1000
  and 23/1000 with pass-through. Their pointwise Wilson intervals include
  the controlled-design references .01 and .0199 respectively.
- Signal only below A is detected locally in every dataset. Without
  pass-through it is reached in 11/1000 and 12/1000; with pass-through it is
  reached and split in all datasets. With a strong root and weak A, local
  rejection counts 81/1000 and 281/1000 become 51/1000 and 213/1000 actual
  A splits after sibling BH. Local power and traversal recovery are distinct.

## Evidence

All 120,000 local tests, 10,000 selected-root tests and 20,000 traversals
completed with status `ok`. The independent verifier regenerates every input,
recomputes all fixed/selected oracle statistics, checks exact selected group
composition, checks chi-square arithmetic, reconstructs every three-node
BH/traversal outcome, and matches 17 source hashes. Hashes were saved after
execution; a separate pre-run source snapshot was not retained. Matched
stage changes, node recovery, ranks, p-value quantiles, skew, actual reached
splits, pass-through decisions and unavailable denominators are preserved.

`make check` passes all 697 tests plus lint, dependency/dead-code audits and
wiki validation. No production method or default was modified for this
diagnostic. No commit or push was made.

A subsequent requested Luna-agent passthrough benchmark reran the same 10,000
datasets with enabled/disabled traversal. Its independent audit reconstructed
all 20,000 outcomes, verified identical local tests within each pair, matched
the earlier traversal rows exactly, and matched 17 source hashes captured
before, after and at verification. This repeat confirms reproducibility; it
does not add independent simulation evidence. The targeted 22-test baseline
passed. Passthrough requires open binary/edge prerequisites, a closed sibling
gate, an available descendant split and any configured support guard. It
records `pass_through` separately from a significant `split`.

This study does not validate categorical/Bernoulli covariance, diffusion,
NNLS, other tree inference methods, empirical inflation, estimated edge
prerequisites, or general selected-tree error control. The subsequent
[[independent-projection-vandermonde-control]] implements the independent-data
comparison: fixed-rank controls isolate a same-data eigenvector effect,
while canonical global-null MP ranks offer almost no rank-selection contrast.

## Links

- [[passthrough-pvalue-normalisation-audit]]
- [[eigenspace-benchmark-calibration-reassessment]]
- [[independent-projection-vandermonde-control]]
- [[diffusion-nnls-calibration-open-points]]
- [[open-mathematical-questions]]
- [[eigenvector-calibration-construction-review]]
