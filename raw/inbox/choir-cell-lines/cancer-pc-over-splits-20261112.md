# The accepted extra splits of today's test on the cancer_pc cases, default tree against Markov/Ward tree

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/cancer_pc_over_splits_20261112` (git-ignored);
protocol, deviations, both scripts and the complete summary output (every split of every run) are in
`cancer-pc-over-splits-scripts-20261112.py.txt`. Working tree on top of `bc044e07` with the uncommitted opt-in
reference; the runs use the registered parameter cells only (no opt-in). No production code was changed.
Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Question

On the 4 `cancer_pc` cases (496 cells, 30 features, 6 true groups of 256, 64, 32, 16, 64, 64 cells) today's
test finds 10, 13, 16 and 15 clusters on the Markov/Ward tree (`tbs_family_diffusion_markov_ward`, mean ARI
0.3033) and 6, 5, 6 and 6 on the default tree (`tbs_family_diffusion`, mean ARI 0.8218). Which splits are
accepted on the Markov/Ward tree, and what differs from the default tree?

## Setup

4 cases x 2 methods through `benchmark_cluster_algorithm`. The tree, the annotation table and the cluster
labels were read by wrapping `tbs_runner.TreeDecomposition` and `tbs_runner.labels_and_report_from_decomposition`
at run time. `measure.py` wrote one record per internal node (3,960 records): cells per side, true-group counts
per side, the sibling and edge columns of the annotation table, whether the node is split in the final clusters,
the reference (`best_flat_cut_reference`) for the same split, and `rho` in all 30 features for nodes with at
least 40 cells. Self-checks passed: found clusters and ARI equal the 11 November run in 8 of 8, and the number
of split nodes equals found clusters minus 1 plus the pass-through nodes in 8 of 8. The protocol was written
before the run; the changes made after the first output are listed under Deviations.

## Part 1: the splits in the final clusters

"between" = accepted by the sibling test and the two sides have different majority true groups and no divided
group; "inside" = accepted and a true group is divided (both sides at least 10% and 5 of its cells) or both
sides have the same majority true group; "pass-through" = not accepted by the sibling test.

| Case | Tree | Clusters | ARI | between | inside | pass-through |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| hub_disconnected | default | 6 | 0.9982 | 5 | 0 | 1 |
| hub_disconnected | Markov/Ward | 10 | 0.4459 | 5 | 4 | 0 |
| hub_connected | default | 5 | 0.5934 | 3 | 1 | 2 |
| hub_connected | Markov/Ward | 13 | 0.2651 | 3 | 9 | 0 |
| nohub_disconnected | default | 6 | 0.9982 | 5 | 0 | 0 |
| nohub_disconnected | Markov/Ward | 16 | 0.2677 | 5 | 10 | 0 |
| nohub_connected | default | 6 | 0.6974 | 3 | 2 | 0 |
| nohub_connected | Markov/Ward | 15 | 0.2344 | 4 | 10 | 0 |

Totals: Markov/Ward 50 split nodes, 17 between and 33 inside; default 22 split nodes, 16 between, 3 inside and
3 pass-through. The reference refuses (p > 0.01) all 33 inside splits of the Markov/Ward tree and all 3 of the
default tree; it also refuses 3 of the 17 between splits of the Markov/Ward tree and 6 of the 16 of the default
tree.

The inside splits of the Markov/Ward tree:

- Disconnected cases: group 0 (256 cells) is cut into 5 and 8 pieces (22|234, 62|172, 94|78, 32|46 with the
  hub; 43|213, 56|157, 30|127, 30|97, 58|39, 22|34, 13|17 without it). Without the hub the groups of 64 cells
  are also halved: group 1 27|38, group 4 29|35, group 5 28|36.
- Connected cases: groups 4 and 5 (offset 3 from group 0) are not separated cleanly (with the hub 47|337 with
  43 cells of group 5 against 21 left behind and 47|290 with 45 cells of group 4 against 19; without it 61|323
  with 52 against 12 and 39|284 with 38 against 26), and the remaining mixture of group 0 with the rest of
  groups 4 and 5 is then cut 7 times in each case. Without the hub group 1 is also halved (27|38).
- 32 of the 33 accepted inside splits have 30 degrees of freedom, an inflation factor of 1.76 to 2.70, a
  statistic of 51 to 183 and a raw p-value of 0.0095 to 8.0e-24; one (13|17) has 4 degrees of freedom, an
  inflation factor of 6.83 and p = 0.0041.

The inside splits of the default tree are all in the connected cases: 53|329, 336|48 and 273|63, each taking
42 to 49 cells of group 4 or 5 together with 6 to 14 cells of group 0.

## Part 2: nodes with at least 40 cells, all from one true group

| | default | Markov/Ward |
| --- | ---: | ---: |
| Nodes | 48 | 24 |
| Share of the smaller side, median (quartiles) | 0.124 (0.069 to 0.180) | 0.382 (0.262 to 0.441) |
| Smaller side under 10 cells | 21 | 1 |
| `rho` in 30 features, median (min, max) | 1.214 (0.088, 7.901) | 3.872 (0.730, 17.368) |
| An edge of the node was tested | 6 | 21 |
| Edge significant for at least one child | 0 | 17 (both children in all 17) |
| Sibling p-value present | 0 | 17 |
| `Sibling_Divergence_Skipped` | 48 | 7 |
| Raw sibling p <= 0.01 | 0 | 12 |
| `Sibling_BH_Different` | 0 | 12 |
| Split in the final clusters | 0 | 12 |
| Reference p <= 0.01 | 0 | 0 |

The reference's reasons: default 48 `no_invariant_directions`; Markov/Ward 13 `no_invariant_directions` and 11
`tested`. At the 17 Markov/Ward nodes with a sibling p-value the statistic per degree of freedom has median
2.240 (1.370 to 4.918) with 30 degrees of freedom, and the inflation factor has median 1.996 (1.757 to 2.474).

The largest node of a true group with all cells from that group, same cells in both trees:

| Case | Group, cells | default: split, `rho` in 30 features, edge p (BH) | Markov/Ward: split, `rho` in 30 features, edge p (BH), sibling p |
| --- | --- | --- | --- |
| hub_disconnected | 0, 255 or 256 | 244\|11, 0.100, 0.789 | 22\|234, 0.730, 8.5e-05 and 4.9e-07, 0.0079 (accepted) |
| hub_disconnected | 4, 64 | 10\|54, 0.334, 0.992 | 29\|35, 4.51, 0.00021, 0.0257 (not accepted) |
| hub_disconnected | 5, 64 | 6\|58, 0.355, 0.904 | 28\|36, 4.05, 6.7e-05 and 1.5e-05, 0.0437 (not accepted) |
| nohub_disconnected | 0, 256 | 6\|250, 0.088, 0.996 | 43\|213, 1.05, 9.8e-16 and 4.1e-20, 4.0e-16 (accepted) |
| nohub_disconnected | 4, 64 | 10\|54, 0.334, 0.991 | 29\|35, 4.51, 0.00021, 0.0012 (accepted) |
| nohub_disconnected | 5, 64 | 6\|58, 0.355, 0.904 | 28\|36, 4.05, 7.2e-05 and 1.5e-05, 0.0028 (accepted) |

Per run, a sibling p-value is present at 4 or 5 nodes of the default tree and at 12 to 17 nodes of the
Markov/Ward tree; these are exactly the nodes with a significant edge for at least one child.

## Reading

- The over-splitting on the Markov/Ward tree is over-splitting inside one true group: 33 of its 50 split nodes.
- The hypothesis of the protocol is supported: inside one true group the Markov/Ward tree's splits are more
  balanced (median smaller side 0.382 against 0.124) and separate the cells more strongly (`rho` in 30 features
  3.87 against 1.21; 4.05 to 4.51 against 0.33 to 0.36 on the same 64 cells), and today's tests reject them.
- The difference between the trees is made at the edge gate, before the sibling test. On the default tree no
  edge inside a true group is significant (edge p 0.79 to 0.996 at the largest nodes), so the sibling test is
  never run there (0 of 48 nodes). On the Markov/Ward tree both edges are significant at 17 of 24 nodes, and
  the sibling test then rejects 12 of the 17.
- So the default tree's 6 clusters do not show that today's sibling test holds its level inside one Gaussian
  group; it is not asked there. Where it is asked on such nodes (Markov/Ward), it rejects 12 of 17 at 0.01.
- The reference refuses every one of the 24 nodes on the Markov/Ward tree and all 33 accepted inside splits.

## Deviations

The classification rule was changed after the first output (a 22|234 split of group 0 alone had been counted
as between groups; pass-through nodes are now listed apart), and Part 2 gained the edge columns, the skipped
count and the table of the largest nodes, because the first output showed that the default tree has no sibling
p-value inside a true group. The columns were already recorded; `measure.py` was rerun only after an import
fix and wrote an identical `nodes.jsonl`.

## Not checked

The two methods differ in the diffusion operator and in the linkage; which of the two makes the balanced
splits was not separated. Why the edge test is significant for the Markov/Ward tree's splits of one Gaussian
group was not examined (branch lengths, selection of the split by Ward). One seed per case, 4 synthetic cases,
nodes with at least 40 cells only; `rho` in 30 features is not comparable between nodes of different size.
