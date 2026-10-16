# tree/

Core tree representation, construction adapters, rooting, and branch-length
handling. Distance and embedding geometry belongs in
`tree_break_selection/space_separation/`; benchmark registration belongs in
`benchmarks/shared/runners/`.

## Construction method map

Tree construction has five distinct stages. Keeping them separate prevents a
branch-time or gate configuration from being mistaken for a new topology
method.

```text
features -> geometry/distance -> topology builder -> rooting -> PosetTree
                                                        |
                                                        v
                                      branch-length fit -> gates/traversal
```

| Topology route | Input geometry | Builder | Root | Active use |
| --- | --- | --- | --- | --- |
| Linkage | Any valid condensed distance | SciPy `linkage`; average, complete, weighted, single, or Ward (centroid and median are rejected: their merge heights are nonmonotone) | Final linkage merge | Canonical TBS, diffusion TBS, applications, and the adaptive-K tree grid |

The registered geometry routes feeding the linkage builder are direct
feature/precomputed distance, fixed Hamming-neighbor diffusion, adaptive
pydiffmap diffusion, and optional graphtools kernel diffusion. Adaptive-cosine
subspace applications add block-coordinate Euclidean or adaptive-diffusion
distance before using average linkage.

`linkage_ultrametric` and `fixed_topology_nnls` are branch-length strategies;
they do not construct different topologies. Sibling gates, edge gates,
traversal profiles, and the benchmark consensus selector are also downstream
of tree construction.

The durable, cross-application audit is in
`wiki/analyses/tree-construction-method-map.md`.

## Construction and representation modules

| Module | Responsibility |
| --- | --- |
| `construction/build.py` | Validates linkage inputs, builds the tree rooted at the final merge, and returns construction evidence |
| `construction/defaults.py` | Immutable canonical binary-distance and linkage defaults; callers pass alternatives explicitly |
| `construction/hierarchical.py` | Promote monotone linkage merge output into `PosetTree`; non-monotone heights fail closed |
| `poset_tree.py` | Stable rooted-tree representation and hierarchy operations |
| `branch_lengths.py` | Convert monotone linkage heights to normalized ultrametric edge lengths |
| `optimized_branch_lengths.py` | Refit non-negative edge lengths on an already fixed topology |

`build_tree` builds linkage trees only and takes no builder or rooting
argument. Neighbor joining, MAD rooting and IQ-TREE were removed, as were the unused sklearn-agglomerative and arbitrary
edge-list adapters and the pass-through `PosetTree.from_*` constructors;
live callers use the construction interface directly.

## poset_tree.py — `PosetTree`

NetworkX `DiGraph` subclass. Central data structure for the entire pipeline.

| Method                                      | What it does                                                                                   |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| `root()`                                    | Return the root node (in-degree 0), cached after first call.                                   |
| `get_leaves(node, return_labels)`           | Collect leaf labels globally or under a subtree.                                               |
| `compute_descendant_sets()`                 | Map every node → frozenset of its descendant leaf labels.                                      |
| `find_lca_for_set(nodes)`                   | LCA for a non-empty collection of nodes.                                                       |
| `populate_node_divergences(leaf_data, feature_space=None)` | Populate flat raw-coordinate distributions and leaf counts. Categorical and continuous data require an explicit feature-space contract. Stores result in `annotations_df`. |

## topology.py

Strict rooted-tree topology helpers used by `PosetTree`, the gates, and Tree-BH.

| Function                             | What it does                                                                 |
| ------------------------------------ | ---------------------------------------------------------------------------- |
| `is_leaf(tree, node)`                | Read the explicit `is_leaf` node attribute.                                   |
| `get_leaf_values(tree, ...)`         | Collect leaf labels or ids without label/id fallback coercion.                |
| `compute_descendant_leaf_sets(tree)` | Bottom-up descendant leaf-set aggregation.                                    |
| `lowest_common_ancestor(tree, a, b)` | Depth-based LCA for nodes in the same rooted tree.                            |
| `lowest_common_ancestor_for_set(...)` | LCA for a non-empty node collection.                                          |
| `tree_root(tree)`                    | Cached `graph["root"]` or the single parentless node; raises otherwise.      |
| `bottom_up_nodes(tree)`              | Deterministic leaves-first, root-last traversal order.                        |
| `compute_node_depths(tree)`          | BFS depth from the root(s), used by Tree-BH and LCA.                          |

## distributions.py

Bottom-up distribution population (called by `populate_node_divergences`).

| Function                                             | What it does                                                                                                                                          |
| ---------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `populate_distributions(tree, leaf_data, feature_space=None)` | Postorder traversal: leaves get validated flat raw-coordinate vectors under the active `FeatureSpace`, and internal nodes get empirical subtree barycenters (leaf-count-weighted means). Sets `distribution` and `leaf_count` on each node; each continuous block gets an empirical descendant covariance under `continuous_covariance_by_block`: full dense covariance when feasible, or a diagonal variance vector in high-dimensional regimes. |
| `_calculate_leaf_distribution(tree, node, data)`     | Set distribution for a single leaf from `leaf_data`.                                                                                                  |
| `_calculate_hierarchy_node_distribution(tree, node)` | Compute an internal node's distribution as the empirical subtree barycenter of its children's distributions; also updates the parent's total descendant leaf count. |

## continuous_distance.py

Continuous tree-distance helpers for Gaussian/Brownian data.

| Function | What it does |
| -------- | ------------ |
| `continuous_time_distance_condensed(data, feature_space, ...)` | Returns condensed pairwise Brownian-time dissimilarities \((x_a-x_b)^T\widehat\Sigma^{-1}(x_a-x_b)/d\) for pure continuous feature spaces. |
| `estimate_continuous_covariance_by_block(data, feature_space, ...)` | Estimates shrinkage empirical covariance blocks with an identity target and scale-relative ridge for stable Mahalanobis-time distances. |
