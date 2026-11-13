# Paths, missing neighborhoods and statistical comparisons

## Node identity and path identity

Match each node by its exact labeled descendant leaf set, on the same leaf
universe. If two matched clades are nested, their ancestor–descendant relation
cannot reverse in another rooted hierarchy. If they are disjoint, they remain
on different branches. Intermediate clades, the lowest common ancestor (LCA),
edge lengths and the comparisons attached to those edges can nevertheless change.

An original intermediate clade has three different outcomes: it remains on the
new path; it survives elsewhere in the tree but leaves that path; or it is absent
from the alternative tree. Across branches, an original LCA can survive above
a newly introduced, smaller LCA and therefore leave the connecting path. On an
ancestor–descendant path, a surviving intermediate clade must remain on the path
because its set lies strictly between the two endpoint sets.

The measured neighborhood is the union of parents, children and siblings of path
vertices, excluding the path itself. Roles can overlap. A neighbor can disappear
as a clade, or survive elsewhere without remaining in that neighborhood. Leaves
do not disappear when both trees have the same leaf universe.

## Mean contrasts telescope; Wald scores do not

For fixed clade C define a_C as its normalized indicator vector: entry i is
1/|C| for i in C and zero otherwise. Its sample mean is a_C^T X. Along an
oriented path C0,...,Cm, let w_j=a_Cj-a_C(j-1). Then exactly

```text
sum_j w_j = a_Cm - a_C0
sum_j (mean_Cj - mean_C(j-1)) = mean_Cm - mean_C0.
```

Changing intermediate nodes cannot change this full-data endpoint contrast
when the endpoint clades and data are held fixed. It can change the resolution
of that contrast into different comparisons.

For independent rows with common feature covariance Sigma and fixed clades,

```text
Cov(w_i^T X, w_j^T X) = (w_i^T w_j) Sigma
Var(mean_V - mean_U)
  = [1/n_U + 1/n_V - 2 |U intersect V|/(n_U n_V)] Sigma.
```

For nested U containing V the coefficient is 1/n_V-1/n_U; for disjoint U,V it is
1/n_U+1/n_V. The variance of the sum uses the whole increment covariance matrix,
including its cross terms, not just the sum of diagonal variances.

There is an instructive exception to blanket claims of path dependence. On a
fixed nested chain, full-data mean increments are mutually uncorrelated under
the common-covariance model, and independent under Gaussian rows. Their sampling
variance coefficients telescope. Their squared standardized scores still do not
equal the direct endpoint score: squaring discards the directions and signs
needed to combine contrasts. A sum would test a different collection of mean
contrasts even in that idealized fixed-path model.

On a cross-branch path, the two increments touching a binary LCA are positively
collinear when both are oriented along the path. They express the same child
partition with size-dependent scaling. For the four-leaf path
`0 -> 01 -> 0123 -> 23 -> 2`, the increment covariance coefficient matrix is

```text
[[1/2,   0,   0,   0],
 [  0, 1/4, 1/4,   0],
 [  0, 1/4, 1/4,   0],
 [  0,   0,   0, 1/2]].
```

With scalar values (0,2,3,7) and known variance one, the individual edge scores
sum to 42. The endpoint score is 9/2. Summing diagonal variance coefficients
gives 3/2; the correct endpoint coefficient is 2. The exact
`wolfram_evidence.json` and deterministic fixture verify these expressions and
the nested-chain diagonal coefficients (1/12,1/4,1/2). This is fixed-membership
algebra, not a test calibrated after selecting the tree or path.

## Why the native edge chain differs

An incoming edge uses the child's test half and the pooled test halves of the
parent's immediate children. The same intermediate node can therefore use one
mean as a child on its incoming edge and a different mean as a parent on its
outgoing edges. The saved masks identify both means. The path diagnostic sums
the oriented native half-row increments and measures the discrepancy from its
first and last edge-context means. It also counts intermediate context changes,
because different changes can cancel in the sum.

Each native edge additionally has its own parent basis, covariance, ridge,
branch multiplier and selected rank. In the normalized branch policy, the edge
factor is 1+b_child/mean_positive_branch. The sibling factor instead uses the
sum of two child branches divided by twice that mean. Eigenvalue weighting in
diffusion, covariance whitening, and PCA ordering/tie/rank selection remain the
separate mechanisms established in the preceding context study.

Consequently, equality of the full endpoint contrast establishes none of these:
equality of native half-row contrasts; equality of projected scores or df;
additivity of Wald statistics; a valid combined p-value; or preserved traversal.
The current analysis constructs no new p-value.

## What a disappearing neighborhood changes

| Missing or changed object | Statistical consequence |
| --- | --- |
| Intermediate parent clade | The original incoming edge hypothesis disappears; replacement edges use different parent means, complements, sample counts and spectral contexts. |
| Original child clade | Its incoming edge and the old sibling partition may disappear; the surviving parent can compare a new child partition. |
| Original sibling clade | The original parent/child context changes, including its pooled covariance and empirical-null contribution. The leaves remain but are regrouped. |
| Original LCA | The central comparison between the two branches is replaced. A surviving LCA moved off the path has a different structural status from an absent LCA. |
| An ancestor outside the connecting path | Its gate can change whether either endpoint is reached from the root, even if the endpoint-to-endpoint path survives. |

Empirical inflation also uses records from elsewhere in the tree, and global
positive branch normalization can change after NNLS refitting. These dependencies
extend beyond the defined local neighborhood. Path preservation alone is not a
guarantee of equal numeric evidence.

For graph distance alone, suppressing a degree-two vertex while summing its
adjacent lengths preserves path distance. Disappearing binary clades generally
have off-path branches, and alternative topologies reattach their leaf groups
and refit NNLS. Thus the new fitted length is not assumed to equal the old path's
length sum. The saved-data comparison measures both explicitly.

## Implication for testing and production

Keep three separate comparisons: a fixed full-data endpoint contrast; native
edge/sibling records attached to each actual path; and reachability under that
tree's root-to-node traversal. A direct endpoint diagnostic can describe what
survives structural changes, but a calibrated test must state its null, covariance
and rank treatment and selection event. Neither the number of deleted nodes,
summed edge scores, path length nor eigenvalue magnitude supplies that reference.
The earlier fixed-union Gaussian Beta and declared-event Ward diagnostics keep
their original scope; they do not automatically cover newly selected paths.
