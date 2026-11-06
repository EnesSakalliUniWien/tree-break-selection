# Sibling inflation factors in the merged cancer rerun

Run on 2026-11-06 on merge `50455316` plus the uncommitted recording change
that writes inflation factors into the sibling annotations.

## Recording change

Sibling annotations now carry `Sibling_Raw_Test_Statistic`,
`Sibling_Reference_Scale`, `Sibling_Inflation_Factor` (applied),
`Sibling_Local_Inflation_Factor` (kernel estimate before the floor at one and
the cap) and `Sibling_Family_Baseline_Inflation_Factor`. Any tested sibling
p-value is `chi2.sf(raw / (reference_scale * factor), df)`. The full traversal
trace (`trace_level="full"`) carries the same five fields per visited node.

## Run

`tbs-inflation-debug-script.py.txt` reran the cancer TBS twice in one process
on one disk-backed distance array, from the same cached PCA and diffusion
coordinates as `tbs-merged-rerun-20261106.md`:

- `applied`: unchanged decision rule.
- `uncapped`: wherever the cap applied, the applied factor is replaced by the
  recorded local factor (floored at one) and the p-value recomputed. Nothing
  else changes.

The script ends with `os._exit(0)` after writing outputs, avoiding the teardown
stall of the two earlier full-data runs. The passes took 110 s and 126 s.

## Results

The applied pass reproduces the merged rerun exactly (ARI 1.0 against its
assignments; 70 clusters, ARI 0.027593 against cell lines). Every stored p-value
reproduces from the recorded raw statistic, reference scale and factor.

| Quantity | Applied | Uncapped |
| --- | ---: | ---: |
| Sibling tests | 78 | 78 |
| `family_baseline` (capped) / `local_kernel` | 39 / 39 | 0 / 39 (39 recomputed) |
| Sibling splits | 45 | 22 |
| Traversal pass-throughs | 27 | 50 |
| Clusters | 70 | 70 |
| Assignments | identical | identical |

- The cap applies to the 39 largest tested nodes (32,679 to 48,877 leaves); all
  local-kernel tests are at nodes of at most 4,069 leaves.
- At the capped nodes the family baseline factor is 16.1542, while the local
  factors are 382 to 801 (23.6 to 49.6 times the baseline). Local-kernel
  factors elsewhere range from 3.911 to 1,305.
- 34 of the 45 splits are capped tests. Without the cap, 23 of them become
  same; 11 remain different. No local-kernel decision changes.
- With pass-through enabled, a node whose sibling gate closes but which has a
  splitting descendant is passed through, so traversal reaches the same 70
  boundaries. In this run the cap changes 23 nodes from split to pass-through
  and no cluster.

The cap therefore does not determine the cancer partition. The partition is
set by the deeper splits and the stops analysed in
`tbs-merged-cell-lines-pvalues-20261106.md`.

## Files

- `tbs-inflation-factors-20261106.csv`: one row per tested node with raw
  statistic, reference scale, df, applied, local and family-baseline factors,
  and p, BH p and decision under both rules.
- `tbs-inflation-debug-applied-summary-20261106.json` and
  `tbs-inflation-debug-uncapped-summary-20261106.json`.
- Local-only outputs: `results/analyses/choir_cell_lines_tbs_inflation_debug_20261106/`.

SHA-256:

- `applied_labels.npy` and `uncapped_labels.npy` 40c18c69bf63520feb426ce3f98d7526de41e1aa39276b56817f4fb99b8bea00
- `applied_sibling_tests.csv` 2b2bafe011da9e4dc04b445a135bd09a01a3a5b6f80dfd997864057b7b6ec394
- `uncapped_sibling_tests.csv` d1926776b120547e467cbedc0798eca2a5421ae9316a9241ed2ec21febed1346

## Limits

The uncapped rule is a counterfactual for this dataset only; the count-null
and benchmark effects of the cap were not separated here. Neither factor makes
the selected-tree sibling p-values calibrated.
