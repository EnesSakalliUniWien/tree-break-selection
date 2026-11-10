# Fixed labeled subtree geometry

2026-11-09, Europe/Vienna. This diagnostic enumerates cuts and describes geometry on a supplied fixed subtree. The module is `benchmarks/diagnostics/calibration/selected/hierarchy/subtree_topology_geometry.py`. The executable examples and hashed outputs are in `results/analyses/subtree_topology_reference_20261109/geometry/`.

The cut universe contains every nonempty unordered bipartition, canonicalized to put row zero on the left. The module bounds cut enumeration at 10 leaves and abstract rooted topology enumeration at 6. It uses the positive centered feature Gram spectrum, groups numerical eigenvalue ties with the existing repository rule, and returns every complete positive eigenvalue block and each distinct leading cumulative block union. For each space it exposes its rank, descending half-open eigenblock ranges, row-space projector, and centered unwhitened projected coordinates. All cuts are evaluated by the existing SciPy strict-affine feasibility LP with unit signed margin. LP column rescaling preserves feasibility. Constant rows yield zero positive-rank candidates.

The returned geometric distances use the projected coordinates. Ward linkage uses distances rounded relative to their maximum at 1e-12 with fixed row order. This numerical tie convention can affect near-ties. Condensed distances retain their unrounded values. Rooted Robinson–Foulds compares nontrivial clades. Edge splits include pendant edges, with the complementary root edges counted once. Cophenetic distances retain Ward linkage heights when those are supplied. An abstract rooted topology uses height equal to the longest edge-count path from each node to a descendant leaf.

| Leaves | Unordered nonempty cuts | Rooted labeled binary topologies |
|---:|---:|---:|
| 2 | 1 | 1 |
| 3 | 3 | 3 |
| 4 | 7 | 15 |
| 5 | 15 | 105 |
| 6 | 31 | 945 |

All enumerated rooted clade families are laminar. The union of their edge splits equals the cut universe at each tested size. This combinatorial family exists independently of the data. A Ward tree is one output per tested spectral space.

- In a four-point rectangle, the two rank-one blocks each admit one strict cut, but those cuts cross. The Ward trees have rooted RF distance 4, condensed distance RMSE 2.828, and Ward cophenetic RMSE 4. The cumulative rank-two space admits six strict cuts and shares the first block's Ward topology, while its cophenetic distances still differ (RMSE 1.155).
- In the square, the two-dimensional tied block is retained whole and admits six strict cuts. The projector, projected distances, strict cuts, Ward topology, and cophenetic distances agree under the tested coordinate rotation. The diagonal cut is infeasible.
- In the coincident-pair example only the cut between distinct positions is feasible. Four centered simplex vertices span rank three and admit all seven cuts.

These examples distinguish spectral rank and multiplicity, affine cut membership, and hierarchy membership. The enumeration of abstract trees is combinatorial, and does not describe which trees a data-selected algorithm can return or establish a probability reference. The finite-precision LP and rounded Ward policy bound the numerical scope of the examples.

Verification used Python 3.12.11, NumPy 2.3.4 and SciPy 1.16.3 on macOS arm64. The 17 focused tests and Ruff checks on both new Python files passed. The demonstration also checks all 1,069 abstract trees up to six leaves for laminarity and checks the tied-block rotation at tolerance 1e-12.

## Ward alternatives supported by distances

The expanded diagnostic enumerates every merge sequence admitted by the Euclidean Ward criterion for two through six leaves. For disjoint current clusters A and B, the increase in within-cluster SSE is |A||B|/(|A|+|B|) times the squared distance between their centroids. The implementation obtains this value from the original squared pair distances. At each step it retains every merge whose cost is at most the minimum plus rtol times the maximum original squared pair distance. The default is rtol=1e-10. The zero setting uses floating-point equality. Canonical ordering by smallest leaf and deduplication of forests remove duplicate histories. The leaf bound is enforced before exploration.

| Geometry | Ward topologies, rtol=0 | Ward topologies, rtol=1e-10 |
|---|---:|---:|
| Rectangle first axis | 1 | 1 |
| Rectangle second axis | 1 | 1 |
| Square | 2 | 2 |
| Rotated square | 2 | 2 |
| Square with both right points shifted right by 1e-11 | 1 | 2 |
| Generic one-dimensional four-point example | 1 | 1 |
| Four-point regular simplex | 15 | 15 |

The rectangle's two axis families have union size 2 and empty intersection. From leaf 0, the first axis places leaf 1 at distance 0 and leaf 2 at distance 4. The second axis places leaf 1 at distance 2 and leaf 2 at distance 0. The first axis distinguishes leaves 1 and 3 at distances 0 and 4 from leaf 0, while the second axis makes both distance 2. Over the 15 comparisons of its six pair distances, the axes reverse four strict distance orders and give equality on one axis only in eight comparisons. The archived output preserves every comparison and every returned topology. Equal Ward costs in the full square admit two distinct rooted topologies. Equal costs in the simplex admit all 15 abstract rooted topologies.

These counts are conditional on the supplied subtree membership, subspace, Ward objective, and numerical tolerance. An empty intersection demonstrates conflicting Ward outcomes across those spaces. It does not identify a true population topology. The number of topology alternatives does not calibrate a sibling test.

## Leaf distances fitted on each topology

The inverse comparison fixes a rooted topology and uses all unordered leaf pairs to fit nonnegative additive edge lengths by least squares. The design marks an edge when exactly one leaf of a pair lies on its canonical split side. Complementary root edges share one column, so the fitted length is their identifiable sum. Root placements with identical unrooted splits yield the same fit. The targets are the supplied Euclidean distances in their original units. The production solver's target uses squared distances, so this small diagnostic uses SciPy NNLS directly.

Each of 15 rooted four-leaf topologies was fitted to each geometry, for 75 fits. The table shows the two axis Ward topologies. Entries give RMSE, with maximum absolute residual in parentheses. Rows retain their own geometry's units.

| Supplied distances | AB \| CD topology | AC \| BD topology | Rooted fits tied for best RMSE |
|---|---:|---:|---:|
| Rectangle first axis | 0 (0) | 1.8856 (2.6667) | 5 |
| Rectangle second axis | 0.9428 (1.3333) | <1e-12 (<1e-12) | 5 |
| Full rectangle | 0.1927 (0.2361) | 1.0716 (1.4907) | 5 |
| Full square | 0.3382 (0.4142) | 0.3382 (0.4142) | 10 |
| Regular simplex | <1e-12 (<1e-12) | <1e-12 (<1e-12) | 15 |

The rectangle axes favor opposing unrooted splits, while the full rectangle favors AB \| CD with positive residual. The full square ties two unrooted splits at positive residual. Triangle inequalities hold for all five input distance matrices and all 75 fitted distance matrices. The simplex has an exact additive star representation and all 15 binary representations contract their internal edge to zero within floating-point tolerance. Near-zero fitted branches leave an unresolved split in the fitted metric. These comparisons measure additive distance fit and do not give a statistical reference.

A fixture with known nonnegative lengths reproduces its six distances to numerical precision, gives the same fit after rerooting, and has positive residual on the incompatible topology. The 28 focused geometry, spectral, rule-out, and production branch-length tests passed. The expanded `geometry.json` contains all fitted lengths, distances, residuals, and source hashes. The pre-fit-extension source is archived as `subtree_topology_geometry_before_distance_fit.py` in the geometry output directory, with the existing integration source history unchanged.
