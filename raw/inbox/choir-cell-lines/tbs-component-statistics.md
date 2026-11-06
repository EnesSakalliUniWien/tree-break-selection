# TBS statistics across disconnected graph components

The three alternative full-data runs stopped when the runner rejected their
NNLS lengths. The exception in `benchmarks/shared/runners/tbs_runner.py` occurs
before node distributions are populated or the edge and sibling gates are
called. Graph disconnection does not enter this guard. The symmetric/average
baseline reached the edge gate. Both root edges failed at alpha 0.001, with
adjusted p-values 0.9167 and 0.7246. The sibling gate was skipped at the root,
and traversal stopped above the disconnected components. Its prior edge replay
without branch-length variance still gave adjusted p-values of 0.12068.

## Frozen-tree replay

On 2026-11-04 we ran the existing gate annotation and traversal on all four
saved topologies and their common 30-PC input. We reapplied the saved native
lengths for the baseline and the independent diagnostic lengths for the other
three. We used the native linkage ratio to divide the combined root length,
normalized branch-length variance, the 0.90 projection energy fraction, and
the original edge and sibling alpha levels of 0.001 and 0.01. These runs
measure the behavior of the existing statistics on these fixed inputs.

| Tree | Root edge rejections | Adjusted root sibling p-value | Root decision | Diagnostic clusters |
| --- | ---: | ---: | --- | ---: |
| Symmetric / average | 0 of 2 | Skipped | Boundary | 1 |
| Symmetric / Ward | 2 of 2 | Underflow | Split | 711 |
| Markov / average | 2 of 2 | 1.50e-28 | Split | 72 |
| Markov / Ward | 2 of 2 | 2.20e-309 | Split | 804 |

Each alternative traversal recovered all 250 A204_SOFT_TISSUE cells and all
185 IM95_STOMACH cells as two whole clusters. Every other cluster contained
only cells from the main graph component. The main component has 48,444 cells
from 188 cell lines. Its subdivision into 709, 70 or 802 clusters depends
strongly on the tree. A204 and IM95 are intact subtrees in all four saved
hierarchies. In the baseline, the one final cluster mixes all three components.

## Direct component contrasts

We compared each pair of graph components in the original PCA feature space,
using the TBS `fixed_global_chi_square` sibling statistic without a branch
time. We also calculated the same mean contrast with separate component
covariances, using `S_a / n_a + S_b / n_b` as the covariance of the difference.

| Components | Cells | Pooled-covariance Wald statistic | Separate-covariance Wald statistic |
| --- | --- | ---: | ---: |
| Main and A204 | 48,444 and 250 | 403,636 | 67,475 |
| Main and IM95 | 48,444 and 185 | 23,981 | 66,510 |
| A204 and IM95 | 250 and 185 | 50,158 | 50,753 |

All exceed the nominal chi-square threshold of 50.89 at alpha 0.01 with
30 degrees of freedom. The comparisons of sample means show large observed
separation under both covariance choices. The same expression data built the
graph and supplied these contrasts. Their nominal p-values, including the
root replay's adjusted p-values, do not establish post-selection error control.
Ordinary sample splitting with data-dependent assignment can still invalidate
the contrast. [Gao, Bien and Witten](https://arxiv.org/abs/2012.02936)

## Handling disconnected components

1. **Use a forest for within-component analysis.** The diffusion transition
   matrix is block diagonal after ordering cells by connected component.
   Diffusion stays in each block and each block has its own stationary mode.
   Run local diffusion coordinates, tree construction and TBS inside each
   component. Retain the component ID in the outputs. A shared root used for
   display can carry a structural link without an estimated branch time or a
   statistical gate that blocks entry into every component.
2. **Test between-component differences as a separate family.** Define the
   question in measurement space, such as equality of the component mean
   vectors, and use a two-sample contrast with an appropriate covariance.
   Freeze graph membership on discovery measurements and test independently
   measured features for those cells under an explicit data model. Count
   splitting is one option under suitable assumptions when raw counts are
   available. This run used deposited CPM values. [Neufeld et al.](https://arxiv.org/abs/2207.00554)
3. **Calibrate any same-data split claim on the full selection pipeline.** For
   a declared one-population null, one possible global diagnostic replays
   preprocessing, graph fitting, component selection and the maximum pairwise
   statistic on each null dataset. A global-null check leaves pairwise error
   control under mixtures and error control inside the forest to be established
   separately. Keep graph connectivity and confirmed expression differences
   as distinct outputs.

## Evidence

- `tbs-component-statistics-evidence.json` preserves the nominal contrasts,
  root gate flags and final component membership counts. Input hashes, cell
  order, graph components and edge order were checked. The native baseline
  edge statistics and adjusted p-values reproduced. The alternative repeat
  runs reproduced their root annotations and traversal counters exactly.
- `tbs-component-statistics-script.py.txt` records the replay. Its direct
  contrast was checked against an independent covariance solve including the
  native ridge. Seven existing gate tests passed.
- Local outputs are under `results/analyses/choir_component_statistics_20261104_102936/`.
  The component-aware forest and its error control have not been implemented
  or validated by these replays.
