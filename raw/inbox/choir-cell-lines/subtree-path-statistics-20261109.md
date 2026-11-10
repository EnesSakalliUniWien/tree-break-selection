# Paths between surviving nodes and their missing neighborhoods

Date: 2026-11-09. The confirmed extension covers both ancestor–descendant and
cross-branch paths, including disappearing neighboring parents, children and
siblings, and relates them to the existing native statistical records.

**When both endpoint clades survive, a connecting path always exists. What can
disappear is the original sequence of intermediate clades and its neighborhood.**
The fixed full-data endpoint mean contrast survives those changes. The native
edge and sibling tests attached to the path need not survive or remain equal.

## Measured scope

Reused the preceding 160 saved fits on 16 Gaussian matrices of 64 rows and six
features, with two reference
trees compared against eight alternatives per matrix. All 256 tree pairs were
analysed without refitting clustering. Exact leaf sets identify nodes; native
linkage numbers are used only to locate records in their original files.

There are 569,954 paths between shared endpoint clades, including singletons.
The full output retains all of them. The primary interpretation below uses the
1,890 paths with two internal endpoints; the 516,096 leaf/leaf paths would
otherwise dominate. Paths, nodes and tree alternatives reuse observations and
are not independent replicates or hypothesis tests.

| Internal endpoint paths | Ancestor–descendant | Cross-branch |
| --- | ---: | ---: |
| Total | 816 | 1,074 |
| At least one original intermediate clade absent | 635 | 993 |
| All original intermediates absent, when originally nonempty | 401 | 205 |
| At least one original neighboring clade absent | 777 | 1,067 |
| Same entire labeled connecting path | 123 | 78 |
| Same connecting path and same neighborhood union | 8 | 3 |
| Endpoint traversal record changes, including reachability changes | 46 | 104 |

Among cross-branch paths, 563 change LCA: 416 lose the original LCA clade;
147 retain it elsewhere, above the new connecting path. In a nested path,
a surviving intermediate clade remains between the same endpoints by set
inclusion, so no corresponding off-path movement occurs.

The neighborhood is the parents, children and siblings of path vertices, excluding
the path itself. On internal-endpoint paths, at least one neighboring parent
disappears in 341 comparisons, a neighboring child in 1,821, and a neighboring
sibling in 1,241. Roles overlap. Globally missing clades are distinguished from
clades surviving elsewhere but no longer belonging to this neighborhood.

## Relationship to the statistics

For a path C0,...,Cm the full-data means obey

```text
sum_j (mean_Cj - mean_C(j-1)) = mean_Cm - mean_C0.
```

If a_C is the normalized leaf indicator of C, the covariance of two increments
with weights w_i,w_j is `(w_i^T w_j) Sigma` for independent rows with common
covariance and fixed memberships. Full-data increments along a fixed nested
chain are orthogonal in these weights. Across branches, the two increments
around a binary LCA express the same split and are positively collinear.

**Neither relationship makes the Wald scores additive.** In an exact four-leaf
example, edge scores sum to 42 but the endpoint score is 4.5. Wolfram verifies
the contrast and covariance algebra in [the saved query](wolfram_evidence.json).
The [mathematical explanation](mathematics.md) states the assumptions and shows
why fixed-membership identities do not give selected-path error control.

The native edge tests use distinct child test halves and pooled parent test
halves. At an intermediate node, the mean used on its incoming edge can differ
from the mean used on its outgoing edge. Parent bases, covariance, fitted branch
factors and ranks can also differ. The native half-row increments fail to
telescope to their own endpoint context means on 657/816 reference nested paths
and 971/1,074 reference cross-branch paths with internal endpoints. These are
diagnostic identities with tolerance 1e-10, not rejection rates. Individual
context changes are retained even when their signed contributions cancel.

Among 608 unique matched directed edges across tree pairs, raw edge T changes
in 91, while df and the Tree-BH-tested flag do not change in this panel. All
paired raw edge T values are finite. This edge-level count is separate from
the many connecting paths that reuse those edges.

Endpoint traversal records change in 150/1,890 internal-path comparisons;
16 involve a change in which endpoint is visited. The connecting path remains
exactly the same in 201 comparisons, but six of those still change an endpoint's
boundary/split decision and sibling gate state. All six have a changed
neighborhood. This association does not isolate neighborhood changes from
global refitting and empirical-reference changes. No decision changes occur in
the 11 comparisons retaining both path and neighborhood union; that finite
observation is not a general invariance result.

## Concrete example: same path, different comparison below its endpoint

In crossed replicate 0, the symmetric/average reference and the full-coordinate
Ward alternative preserve the same direct path from a 32-leaf clade to its
16-leaf child. The parent above the 32-leaf clade changes from 48 to 64 leaves.
Below the 16-leaf endpoint, its child partition changes from sizes 14|2 to 8|8.

| At the same 16-leaf endpoint | Reference | Alternative |
| --- | ---: | ---: |
| Collected raw sibling T | 4.2517 | 78.7587 |
| df | 6 | 6 |
| Sibling branch factor | 1.0000 | 1.0347 |
| Active sibling gate | Unavailable | Reject |
| Traversal decision | Boundary | Split |

Both endpoints remain visited. The original unavailable gate is not a
non-rejection, and the two sibling scores compare different child partitions.
All 16 observations belong to the same known joint group in this simulation.
The alternative fragments that group; whole-partition AMI changes from 1.0000
to 0.9366. Thus a larger raw statistic need not indicate better clustering.
Labels are used only for this evaluation, not to construct the candidate.
The exact clades, original records, input hashes and contrast are retained in
[the example artifact](panel/same_path_example.json).

A separate recorded contraction removes a 48-leaf intermediate between the
64-leaf root and a preserved 32-leaf clade. The fitted path length changes from
.34035 to .27893, and the 32-leaf incoming-edge T from 388.45 to 663.68. The
new comparison is directly against the root; summing or copying the old edge
statistics would not produce that test.

## Production consequence

The diagnostic now records the ordered original and alternative paths, their
LCA, missing and displaced clades, fitted length sums, actual half-row contexts,
edge gates and endpoint traversal records. This extends the earlier node-local
comparison to the relationships between surviving nodes.

Three questions remain separate: whether the same endpoint contrast survives;
which native hypotheses and fitted quantities are now attached to the path and
its neighborhood; and whether root-to-node traversal reaches those endpoints.
Ancestors above the connecting path and empirical-null records elsewhere can
affect the last two. A direct endpoint test would need its own declared reference
covering selection, covariance and rank. No sum of scores, p-values, branches,
missing-node counts or eigenvalues is introduced as a replacement test.

No production preset or statistical correction changed. Existing calibration,
runtime and frozen-validation requirements remain in force.

## Evidence and reproduction

- `protocol.md` fixes the scope and reuse of the 16 datasets.
- `run.py` reads the prior native panel and writes all paths to `panel/paths.csv.gz`,
  plus matched edges, counts, examples and input/source hashes. Its existing
  manifest is protected against overwrite.
- `summarize_endpoints.py` reads those paths and reports endpoint decisions and
  missing-neighborhood roles for all internal/internal pairs.
- `benchmarks/diagnostics/calibration/selected/hierarchy/subtree_path_context.py`
  implements exact clade and path comparison; test 81 covers missing chains,
  displaced LCAs, invalid trees, invariant child order, contrast algebra and
  reconstruction from actual native masks and branch lengths.
- The focused context/path suite passes 31 tests. Final style, wiki and
  provenance checks are recorded in `verification.json`.

Run both report modules with the existing project Python environment and a new
output location for the main `run(output=...)` function. Original native
artifacts and their manifests remain unchanged. The Wolfram semantic-context
endpoint returned an internal error; its separate evaluator successfully returned
the exact expressions preserved here.
