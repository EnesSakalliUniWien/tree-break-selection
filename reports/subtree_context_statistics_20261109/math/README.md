# Native statistics across matched clade contexts

Date: 2026-11-09. Scope is the live projected-Wald edge and native
empirical-inflation sibling route for continuous inputs. The optional
split-count, fixed-subspace, sparse and dense evidence paths are separate.
The exact algebra checks and their limitations are recorded alongside this
audit.

## Match the hypothesis before comparing statistics

Match a focal node C by its exact set of labeled leaves across the two
trees. Match its parent by the parent's leaf set and its two children by
their unordered pair of leaf sets. Numeric linkage node IDs can change
when merge order changes and are unsuitable for cross-tree matching.
A root has no incoming edge. A leaf has no binary outgoing sibling test.

| Matched C | Incoming C vs parent | Outgoing edges and sibling test at C |
| --- | --- | --- |
| Same parent, same children | Same nested group contrast and held-out memberships. | Same child contrasts and local samples. |
| Different parent, same children | Different parent and sibling complement, so different nested comparison and potentially different held-out parent basis and covariance. | Same local child contrasts and samples, with potentially changed branch normalization, tree-BH reachability, calibration and traversal. |
| Same parent, different children | Same C and parent complement in a binary tree, so same incoming held-out memberships and same incoming comparison. | New outgoing child groups and held-out rows. The sibling contrast, guarded covariance, projection dimension inputs and gate evidence can change. |
| Different parent, different children | Incoming comparison changes. | Outgoing comparisons change. |

These are statements about which groups and samples are compared.
They do not assert numerical equality in the same-context cell. A
different whole tree can change fitted branches, their positive global
mean, edge reachability, empirical-null support and pass-through. For
matched comparisons, equal statistics require equality of the actual
inputs to those statistics. For a clade absent from one tree, record
absence and no corresponding native comparison. Do not substitute a
p-value of one, a failed test, a non-rejection, or a final boundary. If a
clade is present but unreachable, record that separately.

## Incoming and outgoing edge statistics

The edge route sorts each child's leaves by their string labels and
alternates positions, starting with a test row. It uses the same rule
for each immediate child, with no linkage-node ID seed. The test
parent mean pools the test halves of all immediate children. The basis
parent rows pool the complementary halves. For the usual distinct
string labels, the split depends on child membership and label order,
not input row order or numerical node IDs. Equal string keys for
distinct exotic label types are outside that simple invariance claim.

For a binary parent, let m and n be the counts of the test halves of
the two children and let d be the difference of their test means.
The first child's nested contrast is n d/(m+n). Its ordinary sampling
variance multiplier is 1/m−1/(m+n). The second child has nested
contrast −m d/(m+n) and multiplier 1/n−1/(m+n). Under iid rows with
one fixed common covariance, both squared normalized contrasts equal
the square for the two test-half child means with factor 1/m+1/n.
This algebra does not make their native p-values interchangeable.
The optional normalized edge branch factor is 1+L_child/mean_positive_L,
with an edge-specific length. Ridge regularization is applied after
variance scaling. Native rank can adapt to the projected contrast.

The edge statistic uses a parent basis recomputed from the basis halves.
Continuous basis whitening can fall back to the full parent-total
covariance if the half-row guard fails. The tested edge contrast is
whitened using the separately resolved full-node covariance, under
the configured covariance policy. Thus a split of basis rows from
contrast rows does not by itself make every fitted quantity independent.
Tree selection on the same data remains another dependence. Edge raw
statistics are calculated before Tree-BH. A finite raw edge p-value can
exist for an ancestor-blocked edge. An existing edge whose parent has
zero-dimensional basis gets the explicit numerical convention T=0,
df=0, p=1. This is a different situation from an absent edge.

## Sibling statistic, eigenvalue magnitudes and gate decisions

For an existing binary split of C, the sibling route uses full child
distributions and sample sizes. In continuous raw coordinates write
d=mean_L−mean_R. Under the continuous full covariance path the
unregularized contrast covariance is

```text
V0 = (1/n_L + 1/n_R) * Sigma_hat
g  = 1 + (L_L + L_R) / (2 * mean_positive_L)
V  = g * V0 + ridge * I
```

The guarded covariance pools within-child residual covariance when
supported and otherwise uses the parent-total covariance. Diagonal
blocks and the full-block jitter policy are handled by the native
whitener. Let z be the whitened contrast and Q the selected orthonormal
PCA rows in its native coordinates. The raw statistic is
T=||Q z||², with nominal df equal to the retained rank and
reference_scale=1. The nominal chi-square law requires a standard
Gaussian contrast with a fixed or appropriately independent subspace
and rank. These requirements are not established by native tree
selection, estimated covariance, and contrast-adaptive rank.

Three spectra play different roles.

| Spectrum | How magnitude enters |
| --- | --- |
| Diffusion-operator eigenvalues | The retained embedding coordinates are weighted by max(lambda,0)^t. This can change Euclidean geometry, the tree, branch fit and selected comparisons. Diffusion lambda is not inserted as a Wald variance. |
| Raw covariance eigenvalues | If that covariance is used in the whitener, its magnitude changes standardized contrast. In the full-rank unprojected no-ridge case, an aligned coordinate contributes d_j²/((1/n_L+1/n_R) g sigma_j²). Native guarded covariance, ridge, and projected Cholesky-whitened coordinates complicate this simple directional picture. |
| Native null-whitened tangent covariance eigenvalues | They order PCA directions, determine positive direction availability and numerical ties, and enter the MP count above (1+sqrt(p_active/n_threshold))² before rank floors/caps. Given the same chosen rows Q, these eigenvalues do not divide the contrast, weight T, or change reference_scale. |

The sibling test uses the all-leaf spectral context returned by the edge
pipeline. Its nonadaptive target rank is derived from the children's own
all-leaf MP dimensions, using the rounded geometric mean when both are
positive and capping at the parent's available test dimension. These
are not the df of the actually tested incoming edges. In the adaptive
route the kernel finds the smallest prefix reaching .90 of the
observed squared contrast energy in the available parent basis, then
retains the complete boundary tie group. The energy is cumulative
sum of (q_i^T z)^2, with no eigenvalue weighting. The kernel's absolute
near-zero energy rule can also affect rank under rescaling.

For example, z=(3,1) in a fixed ordered orthonormal basis reaches
.90 at rank one and gives T=9. Increasing an already strictly leading
PCA eigenvalue from 4 to 100 has no effect if basis and rank inputs
stay fixed. Tying the two eigenvalues can extend the retained block to
rank two and give T=10. Moving the contrast to z=(1,3), with the same
untied spectrum and basis, gives only .10 in the first component.
An energy change is not a covariance eigenvalue change.

Global positive rescaling of all fitted lengths preserves g when the
positive edge set is held fixed. Setting g=1 is different. In the
unregularized fixed context, T scales inversely with g. The production
fixed ridge and adaptive rank can break that identity. Rescaling the
raw data or refitting a topology can also change the branch fit itself.

A raw sibling record is not the final decision. Child-edge evidence
assigns focal and null-like roles, support weights and stopping-event
dependency groups. A contextual empirical factor c_hat, where supported,
produces a displayed T/(reference_scale*c_hat) and plug-in chi-square
tail. Failed calibration can preserve a finite raw p-value in the
annotation while marking the gate invalid/skipped and withholding an
active calibrated decision. The sibling BH frontier further conditions
which adjusted tails enter comparisons. Traversal then records split,
boundary or supported pass-through. A boundary at an ancestor can
make a present node unreachable.

## Extractable trace

- `split_edge_test_rows(tree, leaf_data)` supplies contrast and basis
  membership indices. Record leaf labels as well as counts. To recover
  actual edge-half eigenvalues and basis, call
  `compute_child_parent_spectral_context` with those parent basis rows
  and the same spectral options. The returned
  `edge_gate_result.spectral_context` is the all-leaf sibling context.
- Edge annotations retain `Child_Parent_Divergence_Test_Statistic`,
  `_df`, `_P_Value`, `_P_Value_BH`, `_Tested`,
  `_Ancestor_Blocked` and `_Significant` under the same prefix.
- `collect_sibling_pair_records` supplies native raw T, reference scale,
  df, p, branch sum, focal/null-like role, edge-blocked flag, calibration
  dependency group/weight and spectral support. Annotated values include
  `Sibling_Inflation_Factor`, `Sibling_Test_Statistic`,
  `Sibling_Divergence_P_Value_Corrected`, `Sibling_BH_Different`,
  skip/invalid flags and `Sibling_Gate_P_Value_Role`.
- Use the full traversal trace and actual visited-clade list to tell a
  structural absence, unreachable existing node, closed gate and final
  boundary apart. Use branch sum and the fitted tree-wide positive mean
  to reconstruct the dimensionless factor.

## Connection to earlier evidence

The earlier paired Gaussian null panel had 47/50 root rejections at
nominal .01 in both native NNLS configurations. On its first fixture,
changing branch source opened previously closed outgoing edge gates.
Another saved reversal changed empirical inflation while root raw T,
df and outgoing edge p-values were unchanged. That is direct evidence
that a changed final decision need not come from eigenvalue magnitude,
or even from a changed raw sibling statistic. The earlier covariance
review similarly distinguishes raw covariance from node-local
whitening, and the crossing-partition experiment shows how one
topology can structurally omit a compatible-with-data alternative.
An omitted clade provides no native p-value. Geometry, spectral
magnitude, topology count and a nominal Wald tail each require a
selection-aware reference before supporting calibrated inference.
