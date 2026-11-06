# TBS statistical logic at the cancer cell-line root

Date 2026-11-04. Each replay uses one saved hierarchy and its original 30-PC
input for all 48,879 cells. The symmetric/average replay uses the saved native
NNLS lengths. The Markov/average replay uses the saved independent diagnostic
lengths because the original runner stopped at the NNLS guard. Instrumentation
reproduced every saved root annotation and traversal counter for both replays.

## From selected children to a traversal decision

1. The diffusion coordinates and linkage select the two root children. For
   each edge the test uses a child-minus-parent mean contrast on the test rows.
   Its covariance uses the factor \((1/n_c-1/n_p)(1+\ell_e/\bar\ell)\),
   with the regularization in the Wald implementation. The parent PCA basis
   uses separate projection rows. The adaptive rank retains 0.90 of the
   contrast energy, with ties retained. The projected, whitened squared norm
   is compared with a nominal chi-square tail at that rank. Both edges compare
   the same child-versus-child mean contrast up to count scaling before their
   edge-specific variance adjustments.
2. Tree-BH opens the next family only along a rejected edge. The code computes
   the raw edge statistics before Tree-BH determines which families are reached.
3. The sibling pass computes raw projected Wald records on the full rows.
   Its two-child covariance multiplier is
   \((1/n_L+1/n_R)(1+(\ell_L+\ell_R)/(2\bar\ell))\). A sibling record
   is focal if either immediate edge rejected. Other records can support the
   internal empirical null. Their role comes from the edge gate, including
   ancestor-blocked records, and supplies no known-null label. When focal
   records exist, the projected-Wald inflation method uses
   \(p=\Pr(\chi^2_k\ge Q/(\mathrm{scale}\,\hat c))\), then applies sibling BH.
4. Traversal consumes the gate annotations. A root with both edge gates
   closed becomes a boundary. Pass-through also requires an open edge gate.

| Root quantity | Symmetric / average | Markov / average |
| --- | --- | --- |
| Child sizes | 1 ASPC1 cell / 48,878 | 151 OAW28 cells / 48,728 |
| Child rows used for edge tests | 1 / 24,439 | 76 / 24,364 |
| Parent projection rows | 24,439 | 24,439 |
| Covariance from all parent rows | Diagonal parent-total, singleton guard | Diagonal pooled within-child |
| Edge Wald statistics | 18.3594 / 30.0196 (28 df) | 1921.5478 / 2377.2637 (25 df) |
| Root edge Tree-BH p-values, alpha 0.001 | 0.916734 / 0.724561 | Both underflow to zero |
| Tree-BH reached out of 97,756 raw edge statistics | 2 | 994 |
| Raw root sibling Wald statistic | 23.3033 (30 df), nominal p = 0.802587 | 3995.6124 (25 df), nominal p underflows |
| Sibling records, focal out of 48,878 | 0 | 273 |
| Sibling gate at root, alpha 0.01 | Skipped, no empirical model fit | Adjusted p = 1.50351e-28 |
| Root decision and clusters in replay | Boundary, 1 | Split, 72 |

Turning off edge branch variance while retaining the same edge data, covariance
and basis gave symmetric/average nominal p-values 0.120617 and 0.120680.
The resulting two-edge BH adjustment is 0.120680 for both. The Markov/average
root still has both edge p-values underflow with branch variance off. The
baseline failure is at the immediate singleton contrast, before the sibling
decision can evaluate any split farther down that hierarchy.

## Holding the Markov / average Wald statistic fixed

The empirical fit used 48,605 calibration records in 257 dependency groups.
There were 224 reached null-like records and 48,381 ancestor-blocked records.
The largest calibration parent has 32,681 cells. The root has 48,879. The local
kernel has effective sample size 1.239 and puts 0.891863 of its weight on one
dependency group. Local support falls below the implementation's thresholds.
The replay reports those diagnostics with threshold enforcement disabled.

| Rule evaluated on identical model and root record | Inflation \(\hat c\) | Nominal calibrated p |
| --- | ---: | ---: |
| Committed `8e0093d3`, local kernel | 746.498315 | 0.999988962 |
| Working tree, family baseline fallback | 20.343353 | 1.503510e-28 |

The working-tree rule uses the family baseline when the root exceeds the
calibration size range and local inflation exceeds the family baseline. The
raw \(Q=3995.612392\), reference scale 1 and 25 degrees of freedom are identical
in this comparison. The working-tree annotation stores the adjusted statistic
\(Q/\hat c=196.408744\). The p-value crosses the nominal sibling alpha of 0.01
when the inflation rule changes. The committed comparison evaluates the root
calibration function only. The 72-cluster count belongs to the full working-tree
replay with diagnostic lengths.

## Independent-data control

The archived fixed-tree Gaussian control rejected 50 of 50 selected-data roots
and 0 of 50 independent-data roots with known identity covariance, all 24
dimensions, zero branch variance and nominal alpha 0.01. Guarded covariance
rejected 8 of 50 independent-data roots. The current replay of replicate 0
reproduced the selected and independent known-covariance statistics 76.835983
and 18.917610, with p-values 1.92637e-7 and 0.756369. Its complete merge-history
fingerprint differs from the archived pre-correction fingerprint. The saved
matched self-distance comparison records a merge-history change with the same
clades and root partition for this replicate. This control isolates use of the
tree-selection data in the nominal reference and is a sibling-statistic check.
On the cell-line hierarchy the gates and empirical null still use the selected
data, so these replays give descriptive gate decisions.

## Evidence and reproduction

- `tbs-statistical-logic-trace-evidence.json` preserves raw edge and sibling
  records, annotations, calibration estimates and the Gaussian control.
- `tbs-statistical-logic-trace-script.py.txt` takes `symmetric_average` or
  `markov_average` and an existing output directory. Run it with
  `uv run --no-sync --offline python` from the repository root with the saved
  run directories referenced in the script. It checks exact annotation and
  traversal parity and independently evaluates the calibrated chi-square tail.
- The edge stages live in `child_parent_divergence_annotation/tree_testing.py`
  and `child_parent_divergence_annotation/child_parent_divergence_tree_bh.py`.
  The sibling stages live in `pair_testing/collection/record_collection.py`
  and `inflation_correction/empirical_null_inflation_estimation.py`.
  The traversal decision is in `decomposition/gates/gate_evaluator.py`, under
  `tree_break_selection/hierarchy_analysis/`.
- Gaussian evidence is in `raw/assets/diffusion-attribution-20261102/base_runs.csv`
  and `raw/assets/diffusion-self-distance-20261102/postfix_operator_v2/matched_topology_changes.csv`.
- Eight focused existing tests passed across gate annotation reuse, gate
  pipeline parity and the beyond-range inflation case.
