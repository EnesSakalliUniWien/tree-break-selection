---
title: Sibling Inflation and Traversal Debugging
type: tool
status: draft
updated: 2026-11-06
sources:
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/inflation_adjusted_sibling_tests.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/projection_dimension_annotation.py
  - tree_break_selection/hierarchy_analysis/tree_decomposition.py
  - benchmarks/shared/runners/tbs_runner.py
  - tests/statistics/30_test_gate_annotation_pipeline_parity.py
  - tests/hierarchy_analysis/test_tree_decomposition.py
  - raw/inbox/choir-cell-lines/tbs-inflation-debug-20261106.md
  - raw/inbox/choir-cell-lines/tbs-inflation-debug-script.py.txt
tags:
  - debugging
  - calibration
  - traversal
---

# Sibling Inflation and Traversal Debugging

## Summary

Every TBS run's annotations record what lies behind each sibling p-value, and the
full traversal trace records what each visited node decided. Together they
answer which factor produced a p-value, which decisions a rule changes and
whether the final partition depends on them, without instrumenting the code.

## Usage

1. Read the five inflation columns from the annotations of any tested sibling
   parent: `Sibling_Raw_Test_Statistic`, `Sibling_Reference_Scale`,
   `Sibling_Inflation_Factor` (applied), `Sibling_Local_Inflation_Factor` (kernel
   estimate before the floor at one and the cap) and
   `Sibling_Family_Baseline_Inflation_Factor`. `Sibling_Test_Statistic` is the
   raw statistic already divided by reference scale and applied factor.
2. Recompute a p-value under any factor `c` as
   `chi2.sf(raw / (reference_scale * c), df)`. The cap applied where
   `Sibling_Test_Method` ends in `:family_baseline` and the local factor exceeds
   the family baseline.
3. For node-level decisions, pass `trace_level="full"` to
   `run_tbs_on_distance` (or `TreeDecomposition`) and save
   `result.extra["traversal_trace"]`. Each visited node has its decision (split,
   pass-through or boundary), child-edge raw and BH p-values, sibling p-values,
   the five inflation fields and the pass-through audit.
4. For a counterfactual rule, replace
   `inflation_adjusted_sibling_tests.decide_empirical_null_calibration` with a
   wrapper, rerun on the same distances, and compare sibling decisions,
   traversal counters and labels. `tbs-inflation-debug-script.py.txt` does this
   for the cap and shares one disk-backed distance array between passes.
5. Compare partitions, not only p-values: with pass-through enabled, a closed
   sibling gate above a splitting descendant becomes a pass-through and can
   leave the clusters unchanged.

## Evidence

On the merged cancer rerun every stored sibling p-value reproduces from the
recorded columns. The cap applies to the 39 largest tested nodes (local
factors 382 to 801 against a baseline of 16.1542). Removing it turns 23 splits
into pass-throughs and leaves all 70 clusters identical
(`raw/inbox/choir-cell-lines/tbs-inflation-debug-20261106.md`). Tests check that
annotated p-values reproduce from the columns on a benchmark case and that the
full trace carries the factors.

## Links

- [[choir-singleton-peel-rerun]]
- [[choir-cell-lines-recovery-plan]]
