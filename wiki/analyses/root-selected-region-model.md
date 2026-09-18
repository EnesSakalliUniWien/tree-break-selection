---
title: Root Selected Region Model
type: analysis
status: reviewed
updated: 2026-09-17
sources:
  - reports/latest_changes_benchmark_review_20260909.md
  - wiki/analyses/method-proof-web.md
  - wiki/analyses/selected-hierarchy-selection-geometry.md
  - wiki/analyses/selected-hierarchy-geometric-law-map.md
  - wiki/sources/edge-selection-null-audit-20260601.md
  - wiki/sources/selected-hierarchy-null-audit-20260601.md
  - wiki/sources/selected-ratio-tail-law-diagnostic-20260602.md
  - wiki/sources/root-selected-region-margins-20260603.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/selected_hierarchy_geometry_covariates.py
  - benchmarks/diagnostics/calibration/root/selected/root_selected_region_margins.py
  - benchmarks/shared/tbs_tree_context.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_radial.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_hierarchy.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_diffusion.py
  - tests/statistics/51_test_selected_gaussian_radial.py
  - tests/statistics/52_test_selected_gaussian_hierarchy.py
  - tests/statistics/53_test_selected_gaussian_branch_lengths.py
  - tests/statistics/54_test_selected_gaussian_diffusion.py
tags:
  - analysis
  - selection
  - geometry
  - proof
---

# Root Selected Region Model

## Summary

This page defines the first tractable selected-region object for Tree-Break Selection: a root
sibling context with a fixed hierarchy-construction procedure, fixed feature
chart, fixed projected-Wald kernel, no sibling FDR, and no traversal. It does
not solve the full production selected-hierarchy law. It defines the event
that must be conditioned on before a fixed-object Wald reference can be used
honestly after same-data hierarchy selection.

The target object is
\[
\mathcal L\!\left(R_\rho\mid X\in\mathcal S_\rho\right),
\qquad
R_\rho=\frac{W_\rho}{a_\rho\nu_\rho},
\]
where \(\rho\) is the root parent, \(W_\rho\) is the raw root sibling
projected-Wald statistic, and \(\mathcal S_\rho\) is the event that the same
data select the observed root split and open the child-parent edge path.

## Details

### Simplified Setting

Fix the following objects before analysis:

1. a feature-space chart and null covariance contract;
2. a deterministic hierarchy-construction algorithm \(H\);
3. a linkage score \(D_X(A,B)\) for every candidate pair of current clusters;
4. a deterministic tie-breaking rule;
5. an edge alpha \(\alpha_{\mathrm{edge}}\);
6. a fixed projected-Wald edge and sibling kernel.

Ignore sibling FDR, empirical inflation, pass-through traversal, and downstream
non-root focal selection. The only selection target is the observed root split
\[
\rho \to (L_\rho,R_\rho).
\]

### Hierarchy Selection Event

Agglomerative hierarchy construction is a sequence of merges. At merge step
\(t\), let \(\mathcal P_t(X)\) be the set of available cluster pairs and let
\((A_t,B_t)\) be the pair merged by the observed run. With deterministic
tie-breaking, the event that the same merge sequence occurs is
\[
\mathcal H_\rho
=
\bigcap_t
\bigcap_{(A,B)\in\mathcal P_t(X)\setminus\{(A_t,B_t)\}}
\left\{
D_X(A_t,B_t)\le D_X(A,B)
\right\}.
\]

This is already a selected region in data space. The sets
\(\mathcal P_t(X)\) are determined recursively by earlier merge inequalities,
but for a fixed observed merge sequence they become fixed candidate sets
inside that sequence cell.

For Euclidean average linkage, each linkage score is an average of squared
pairwise distances or Euclidean distances. In the squared-distance form,
\(D_X(A,B)\) is quadratic in \(X\), so each merge comparison is a quadratic
inequality. Therefore the fixed-sequence root-selection region is
semi-algebraic:
\[
\mathcal H_\rho
=
\{X:\ g_i(X)\le 0,\ i=1,\ldots,m_H\}
\]
for finitely many polynomial inequalities \(g_i\). For Hamming/binary
distances, the region is piecewise linear after fixing the sign/coordinate
cells of the absolute-value terms.

### Edge-Opening Event

Let \(Q_{\rho\to L}(X)\) and \(Q_{\rho\to R}(X)\) be the child-parent
projected-Wald edge statistics for the two root children, with corresponding
reference thresholds \(q_{L,\alpha}\) and \(q_{R,\alpha}\). The root edge-path
opening event is
\[
\mathcal E_\rho
=
\left\{
Q_{\rho\to L}(X)\ge q_{L,\alpha}
\right\}
\cap
\left\{
Q_{\rho\to R}(X)\ge q_{R,\alpha}
\right\},
\]
or the Tree-BH-adjusted analogue when the root path is embedded in the full
tree correction. In this simplified root model, the unadjusted form is the
minimal conditioning object.

Each \(Q\) is a quadratic form in the null-whitened child-parent contrast when
the projection is fixed. When the projection is selected from local spectral
rows, the event also includes the spectral-selection cell that fixes the
chosen projection dimension and basis.

For a fixed projected edge subspace, write
\[
Q_{\rho\to c}(X)=\|P_{\rho\to c}z_{\rho\to c}(X)\|_2^2,
\qquad
q_{c,\alpha}=\chi^2_{k_c,1-\alpha_{\mathrm{edge}}}.
\]
The local radial signed distance from the observed projected edge vector to
the opening boundary is
\[
d_{\mathrm{edge},c}(X)
=
\sqrt{Q_{\rho\to c}(X)}-\sqrt{q_{c,\alpha}}.
\]
The root edge path requires both child edges to open, so the path distance is
\[
d_{\mathrm{edge},\rho}(X)
=
\min_{c\in\{L_\rho,R_\rho\}} d_{\mathrm{edge},c}(X).
\]
This is an actual fixed-subspace boundary coordinate. It is still not the
full selected law because the production hierarchy also selects the tree,
local spectral basis, and downstream sibling context from the same data.

### Edge-Sibling Algebra At A Binary Parent

For a binary parent \(u\) with children \(L\) and \(R\), the implementation
stores the parent distribution as the leaf-count barycenter
\[
p_u
=
\frac{n_Lp_L+n_Rp_R}{n_u},
\qquad
n_u=n_L+n_R.
\]
Therefore the raw child-parent contrasts are not new directions:
\[
p_L-p_u
=
\frac{n_R}{n_u}(p_L-p_R),
\qquad
p_R-p_u
=
-\frac{n_L}{n_u}(p_L-p_R).
\]
In the no-branch-time limit, the variance scales are
\[
s_{\mathrm{sib}}
=
\frac1{n_L}+\frac1{n_R}
=
\frac{n_u}{n_Ln_R},
\]
and
\[
s_{L,u}
=
\frac1{n_L}-\frac1{n_u}
=
\frac{n_R}{n_Ln_u},
\qquad
s_{R,u}
=
\frac1{n_R}-\frac1{n_u}
=
\frac{n_L}{n_Rn_u}.
\]
The current runtime multiplies these sampling scales by normalized branch-time
factors, \(1+t_{u,c}/\bar t\) for child-parent edges and
\(1+(t_{u,L}+t_{u,R})/(2\bar t)\) for siblings. Without those multipliers, both
edge and sibling contrasts use the same parent/null covariance at a binary
parent, so these scales cancel the barycentric coefficients:
\[
z_{L,u}=z_{L,R},
\qquad
z_{R,u}=-z_{L,R}.
\]
Thus, locally, the root child-parent edge tests and the root sibling test
measure the same null-whitened barycentric direction. They differ through the
projection dimension, the chi-square opening threshold, Tree-BH edge
correction, sibling FDR/inflation, and same-data selection conditioning.

### Fixed-Projection Edge-Conditioned Sibling Law

The barycentric identity makes the first edge/sibling conditioning calculation
explicit in the fixed-subspace root model. Let
\[
X\sim\chi^2_{k_s}
\]
be the sibling projected energy and let
\[
Y\sim\chi^2_{k_e-k_s}
\]
be the additional parent-edge projection energy when the edge projection
contains the sibling projection as a prefix. The raw edge-opening event is
\[
X+Y\ge q_{1-\alpha_{\mathrm{edge}},k_e}.
\]
For an observed sibling statistic \(w\), the fixed-projection tail conditional
only on this edge opening is
\[
p_{\mathrm{edge-cond}}(w)
=
\Pr\left(
X\ge w
\mid
X+Y\ge q_{1-\alpha_{\mathrm{edge}},k_e}
\right).
\]
When \(k_s=k_e=k\), this reduces to the truncated chi-square tail
\[
p_{\mathrm{edge-cond}}(w)
=
\frac{
\Pr\!\left(\chi^2_k\ge \max(w,q_{1-\alpha_{\mathrm{edge}},k})\right)
}{
\Pr\!\left(\chi^2_k\ge q_{1-\alpha_{\mathrm{edge}},k}\right)
}.
\]
This law is mathematically useful because it isolates the action of the root
edge-opening boundary. It is not the production selected-hierarchy law: it
does not condition on hierarchy construction, Tree-BH selection cells,
selected PCA, sibling FDR, empirical inflation, traversal, or non-root focal
selection.

The schema `v6` diagnostic shows that this restricted edge conditioning is
not enough to explain the current root blocker cases. For the diffuse
dimensional root examples, raw and edge-conditioned sibling p-values remain
far below `SIBLING_ALPHA = 0.01`; the block appears only after the internal
empirical-inflation layer. Therefore the next object is not another
edge-opening truncation formula. It is the relationship between internal
empirical inflation and the fuller selected-hierarchy law.

### Root Sibling Selected Region

The root selected region is
\[
\mathcal S_\rho=\mathcal H_\rho\cap\mathcal E_\rho.
\]

The correct root selected law is therefore
\[
\mathcal L\!\left(
W_\rho
\mid
X\in \mathcal H_\rho\cap\mathcal E_\rho
\right),
\]
not the unconditional fixed-object law.

### Nuisance-Derived Continuous Gaussian Radial Prototype

The 2026-09-09 review found numerical counterexamples to the implementation's
exactness claim: a narrow hierarchy interval gives selected p-value 0.8561354
instead of approximately 0.8, and absolute chi-tail underflow rejects a positive
selected event. The 2026-09-17 fix preserves the supplied floating-point path
factors as exact rational values during polynomial construction, linkage
updates, and discriminant evaluation, then solves roots at 80 decimal digits
before returning floating-point boundaries. Tail integration now normalizes
log interval masses, with scaled density integration for narrow intervals.
Regressions cover the counterexamples and independent reconstructed-distance
checks. This repairs the reported numerical failures without certifying all
floating-point boundary cases or establishing production selective validity.
The original review is preserved in [[latest-changes-benchmark-review-20260909]].

The internal prototype in
`tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_radial.py`
implements the smallest exact conditional object for a continuous Gaussian
contrast with known fixed feature covariance. Let \(X\in\mathbb R^{n\times q}\),
let \(\eta\in\mathbb R^n\) be the tested row contrast, and set

\[
h=\lVert\eta\rVert_2,
\qquad
X_\perp
=
\left(I-\frac{\eta\eta^\top}{h^2}\right)X.
\]

For a fixed positive-definite covariance \(\Sigma=LL^\top\), define the
whitened contrast

\[
z=L^{-1}\frac{X^\top\eta}{h}.
\]

The prototype learns an orthonormal projection \(P\) only from the whitened
nuisance matrix \(X_\perp L^{-\top}\). It then conditions on \(X_\perp\), \(P\),
the projected direction

\[
u=\frac{Pz}{\lVert Pz\rVert_2},
\]

and the projection-orthogonal coordinate

\[
z_\perp=(I-P^\top P)z.
\]

Only the projected radius \(r\ge 0\) varies along

\[
X(r)
=
X_\perp
+
\frac{\eta}{h}
\left[L\left(P^\top(ru)+z_\perp\right)\right]^\top.
\]

Under the fixed-covariance Gaussian contrast null, \(R=\lVert Pz\rVert_2\)
has a chi law with \(k=\operatorname{rank}(P)\) degrees of freedom after the
conditioned coordinates are fixed. The Gaussian score/likelihood-ratio
statistic is \(R^2\). If an exact hierarchy replay supplies the selected radial
region as an interval union \(S\subseteq[0,\infty)\), the selected p-value is

\[
p_{\mathrm{sel}}
=
\frac{
\int_{S\cap[r_{\mathrm{obs}},\infty)} f_{\chi_k}(r)\,dr
}{
\int_S f_{\chi_k}(r)\,dr
}.
\]

The implementation evaluates these interval integrals deterministically with
log-domain chi CDF and survival functions, using scaled quadrature when finite
interval tails nearly cancel. Absolute probability fields may underflow to
zero while the normalized p-value remains representable. Focused tests verify exact reconstruction of
the observed matrix, invariance of the contrast-orthogonal nuisance and
projection-orthogonal contrast along the path, nuisance-only projection
selection, the unconditional chi limit, and a disconnected selected region.

The companion implementation in
`tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_hierarchy.py`
now constructs the hierarchy component of \(S\) for one explicit model:
deterministic average linkage over pairwise squared-Euclidean distances along
the same radial path. Since every row has the affine form
\[
x_i(r)=a_i+r b_i,
\]
every leaf-pair squared distance is quadratic:
\[
\lVert x_i(r)-x_j(r)\rVert_2^2
=
\lVert a_i-a_j\rVert_2^2
+2r(a_i-a_j)^\top(b_i-b_j)
+r^2\lVert b_i-b_j\rVert_2^2.
\]
Average linkage preserves that form. At each observed merge step \(t\), the
selected pair \((A_t,B_t)\) is compared with every active competitor
\((C,D)\), giving the complete constraint family
\[
D^{(2)}_{A_t,B_t}(r)-D^{(2)}_{C,D}(r)\le 0.
\]
The constructor enumerates every nonnegative real root of these quadratic
differences, evaluates the constant-sign intervals between roots, and replays
the full merge signature in every retained interval interior. It preserves
disconnected cells and does not merge distinct nearby roots. Deterministic
tie-breaking can differ at finitely many equality points; these boundary
points have zero probability under the continuous chi-radial law and do not
change the selected p-value.

### Branch-Length Compatibility Boundary

The fixed-covariance radial construction supports a fixed scalar branch-time
variance multiplier. If

\[
\Sigma'=m\Sigma,
\qquad m>0,
\]

then the adjusted path satisfies

\[
X_{\Sigma'}(r/\sqrt m)=X_\Sigma(r),
\]

and its squared-Euclidean hierarchy interval union is the original interval
union divided by \(\sqrt m\). A focused check with the runtime multiplier
\(m=1+0.4/0.2=3\) verifies reconstruction, interval-boundary scaling, and the
independently derived truncated-\(\chi_1\) tail. This covers a branch-time
multiplier held fixed along the path; it does not cover branch lengths refitted
from each candidate matrix.

A 100-case random positive-definite covariance challenge initially exposed a
floating-point violation of this scaling law. Computing the radial slope as
`reconstruct(1) - reconstruct(0)` perturbed a mathematically cancelled
quadratic merge coefficient to about machine epsilon; its sign then created a
false remote interval near radius \(10^{15}\). The hierarchy implementation now
constructs the shared raw radial direction directly from the conditioned path
factors, so equal slope-energy coefficients cancel before root solving. The
captured counterexample and a fresh 1,000-case challenge pass after this
correction; the random challenge is numerical evidence, not a universal proof.

Fixed-topology NNLS branch-length refitting preserves the observed merge
signature by construction. A four-leaf integration check verifies that the
analytic merge signature agrees with the SciPy squared-Euclidean average-linkage
tree before and after NNLS changes its edge lengths.

The stronger compatibility assumption fails when branch lengths and
branch-length-adjusted internal nodes are recomputed along the radial path. At
two radii inside one retained merge-signature interval, both
`linkage_ultrametric` and `fixed_topology_nnls` produce different branch
lengths, and `branch_length_state` produces different one-dimensional root
spectral projectors. Therefore the hierarchy merge event \(\mathcal H_\rho\)
does not fix the branch lengths or the branch-length-adjusted internal-node
projection. The current chi-radial law cannot be claimed for either recomputed
internal-node variant until that additional data-dependent selection is
conditioned on, or the branch lengths and projection are derived only from
coordinates already fixed by the conditioning construction.

### Diffusion-Hierarchy Compatibility Boundary

The two production diffusion families require different radial treatment.
Adaptive pydiffmap can accept each continuous candidate matrix, but it refits
the neighbor graph, local bandwidths, kernel eigenspace, diffusion distances,
and average-linkage hierarchy at every radius. A deterministic 12-leaf
continuous fixture was replayed at \(0.2r_{\mathrm{obs}}\),
\(r_{\mathrm{obs}}\), and \(2r_{\mathrm{obs}}\). The first merges were,
respectively, leaves \(9/10\), \(6/11\), and \(0/11\); all three complete merge
signatures and diffusion-distance vectors differed. The observed pydiffmap
epsilon also changed across those radii. Thus the squared-Euclidean quadratic
constraint constructor cannot be reused for adaptive diffusion. An exact
selected event would additionally have to condition on or solve every
data-dependent diffusion fit and its resulting hierarchy cell.

Fixed Hamming-neighbor diffusion and adaptive pydiffmap with Hamming metric
have a stricter support boundary. They both replay the observed binary matrix,
but a nearby Gaussian radial candidate such as \(0.5r_{\mathrm{obs}}\) is no
longer binary or one-hot. The implementation rejects such candidates instead
of silently applying Hamming distance to continuous values. Consequently,
these variants do not have a continuous chi-radial selection interval under
this Gaussian construction; they need a discrete selected-region law or a
different conditioning path that remains inside the binary feature space.

This closes hierarchy interval construction only for the declared
squared-Euclidean average-linkage prototype. It does not cover ordinary
Euclidean average linkage, exact adaptive-diffusion selection cells, Hamming
selected-region laws, standardized-Euclidean or estimated-Mahalanobis
transformations, other tree builders, edge-path opening, focal sibling
selection, covariance estimation, recomputed branch-length-adjusted internal-node
projections, discrete feature families, or production gates. The next
conditional boundary is to express and intersect the fixed-projection
edge-opening event with this hierarchy interval union.

### Proposition: The Root Selected Law Is Generally Not Chi-Square

Assume the fixed-object root sibling statistic satisfies
\[
W_\rho\sim\chi^2_k
\]
when the root split, edge path, and projection are fixed independently of the
tested data. If \(\mathcal S_\rho\) depends on the same data and has nonzero
association with \(W_\rho\), then
\[
\mathcal L(W_\rho\mid X\in\mathcal S_\rho)
\ne
\chi^2_k
\]
in general.

Proof. This is an application of the conditioning counterexample in
[[method-proof-web]]. Since \(\mathcal E_\rho\) contains lower-bound
constraints on child-parent Wald energies, and \(\mathcal H_\rho\) selects
clusters by small between-cluster linkage or large induced barycentric
separation, the event is a nontrivial function of the same empirical
coordinates that enter \(W_\rho\). Unless the event is independent of
\(W_\rho\), conditioning changes the law. Independence is not implied by the
method construction and is contradicted by the edge-selection null audit.

### Differential-Geometric Objects

Inside a fixed smooth sequence cell, write
\[
\mathcal S_\rho=\{X:\ g_i(X)\le 0,\ i=1,\ldots,m\}.
\]
At an observed point \(x\), the active set is
\[
\mathcal A(x)=\{i:\ g_i(x)=0\}.
\]
The tangent cone is
\[
T_{\mathcal S_\rho}(x)
=
\{h:\ \nabla g_i(x)^\top h\le 0\ \text{for all } i\in\mathcal A(x)\}.
\]
The signed distance in the null-whitened metric is
\[
d_{\mathcal S_\rho}(x)
=
\inf_{y\in\partial\mathcal S_\rho}
\left\|
\Sigma_0^{-1/2}(x-y)
\right\|.
\]
Curvature is carried by the Hessians \(\nabla^2 g_i(x)\) of active smooth
constraints after projecting onto the local tangent boundary. These are the
differential-geometric objects suggested by the selected-inference literature:
selected region, boundary, signed distance/action, tangent cone, and curvature.

The edge p-value action
\[
A_{\mathrm{edge},\rho}
=
\min_{c\in\{L_\rho,R_\rho\}} -\log_{10}p_{\rho\to c,\mathrm{BH}}
\]
is therefore best interpreted as a tail-probability action for the root edge
path. The radial distance \(d_{\mathrm{edge},\rho}\) above is the direct
fixed-subspace geometric boundary coordinate.

### First-Order Merge-Boundary Law

For a fixed merge step, write the selected average-linkage pair as
\((A_t,B_t)\) and the nearest competitor pair as \((C_t,D_t)\). The observed
merge inequality is
\[
g_t(X)
=
D_X(A_t,B_t)-D_X(C_t,D_t)
\le 0.
\]
The observed score-space margin is
\[
m_t(X)
=
D_X(C_t,D_t)-D_X(A_t,B_t)
=
-g_t(X).
\]
When \(D_X\) is average Euclidean linkage and all involved pairwise distances
are nonzero, this cell is smooth. The first-order signed distance from the
observed point \(x\) to the nearest merge boundary in the ambient Euclidean
data metric is
\[
d_t^{(1)}(x)
=
\frac{m_t(x)}{\|\nabla g_t(x)\|_2}.
\]
This is not yet the null-whitened selected-region distance and not a
calibration law. It is the first local geometric object that can be computed
from the current tree builder without changing production inference.

For continuous Euclidean cases with independent leaves and root
empirical-Gaussian feature covariance \(\widehat\Sigma_\rho\), schema `v5`
also records the first-order null-whitened scale
\[
\sigma_{g,t}^2
=
\nabla g_t(x)^\top
\left(I_n\otimes\widehat\Sigma_\rho\right)
\nabla g_t(x),
\]
and the corresponding diagnostic distance
\[
d_{\Sigma,t}^{(1)}(x)
=
\frac{m_t(x)}{\sigma_{g,t}}.
\]
This is still a diagnostic geometry coordinate. It is not a production
selected-tail law because it conditions only on the local merge inequality and
does not include edge-opening, spectral selection, sibling FDR, or traversal.

For Hamming/discrete metrics, or when the selected merge is tied with another
minimum pair, the smooth formula is not the right object. The local selected
region is a discrete or nonsmooth tie cell. In that regime the next law is not
a curvature correction to \(d_t^{(1)}\); it is a discrete selected-region law
or a polyhedral/tie-breaking tangent-cone model.

### Connection To Existing Diagnostics

The selected-hierarchy null audit simulates from the conditional region
approximately by regenerating null data and reapplying \(H\) and the edge
tests. The geometry covariate diagnostic then records proxies for the
mathematical objects above:

- edge action as a signed-distance or large-deviation proxy;
- eigenvalue excess and effective rank as local spectral-mode descriptors;
- \(\cos^2\theta\) as selected-subspace angular alignment;
- parent size and child balance as finite-sample fluctuation descriptors.

The next proof-level development is not to add another scalar correction. It
is to decide whether these proxies can be replaced by actual functions of
\(\mathcal S_\rho\): active constraints, distances to boundaries, tangent
cones, and curvature terms.

### Observed Root Margin Diagnostic

`benchmarks/diagnostics/calibration/root/selected/root_selected_region_margins.py` now
extracts one concrete part of \(\mathcal H_\rho\): for each average-linkage
merge step, it checks that the observed merge is a minimum active pair and
records the nearest-competitor margin. The final root merge has no competitor,
so the diagnostic summarizes the inequalities that construct the two root
child clusters.

The 2026-06-03 representative run separates two geometries. In Hamming,
discretized, and categorical cases, many root-child construction merges are
exact or numerical ties, and schema `v2` classifies those rows as requiring
discrete tie-cell geometry. In the continuous Euclidean diffuse representative,
all `178` root-child construction constraints have smooth first-order
geometry. The minimum merge margin is about `4.05e-4`, and the minimum
first-order signed distance after gradient normalization is about `2.24e-4`.
This means the selected-region boundary is not one uniform object: tie-heavy
discrete hierarchy cells and positive-margin continuous cells should not be
collapsed into a single scalar selected-tail explanation.

The same run records large root selected sibling ratios in all representatives
(`32.7` to `755`). Therefore merge margins alone do not explain the selected
tail. They are one component of the selected region and must be read together
with edge-opening strength, spectral excess, feature family, and covariance
geometry. Curvature remains explicit missing work: the continuous diagnostic
marks the high-dimensional Hessian operator as not materialized.

The supported eight-case continuous rerun makes that limitation sharper. Raw
merge margin, ambient first-order distance, and null-whitened first-order
merge distance have weak descriptive relationships with the log root sibling
selected ratio. The schema `v5` edge-opening fields
`root_edge_path_radial_distance`, `root_edge_path_statistic_margin`, and
`root_edge_path_bh_action` each have Spearman correlation `1.0` in this small
panel. Therefore the next selected-region object is the edge-opening
boundary/action, not another adjustment to merge-margin geometry.

Schema `v5` then checks the edge/sibling relationship directly. In all five
representative rows and all eight supported continuous rows, the default
no-branch-scaling barycentric z-identity is verified up to numerical
residuals below `4.17e-8`. In the supported continuous root panel, the sibling
projection dimension equals the parent edge projection dimension in every
row, so the edge and sibling statistics are the same projected energy before
different thresholds/FDR/inflation are applied. In categorical representative
rows, the parent edge projection can retain additional spectral energy beyond
the sibling projection prefix; this is a projection-policy difference, not a
different raw barycentric contrast.

Schema `v6` adds the fixed-projection edge-conditioned sibling tail and the
current internal empirical-inflation sibling tail. The result is a clean
separation: edge conditioning alone does not turn strong diffuse-dimensional
root sibling signals into blockers, while empirical inflation does. The
high-cardinality categorical representative remains explicitly unsupported by
the internal empirical-null contract. This keeps the next mathematical target
honest: diagnose the empirical-inflation/selected-hierarchy law, not the
already isolated root edge-opening law.

### Minimal Diagnostic Contract

A focused root selected-region diagnostic should record, for each selected
root context:

```text
root merge sequence id
root child leaf sets
root sibling W
root selected ratio R
root sibling raw p-value
root sibling edge-conditioned p-value
root sibling empirical-inflation p-value or unsupported status
root edge-path radial distance
root edge-path statistic margin
root edge-path action
edge/sibling z-identity residual
edge extra parent-projection energy
active or near-active merge inequalities
smooth first-order signed distances or discrete tie-cell status
null-whitened first-order merge distances for continuous Euclidean cells
edge-statistic margins to threshold
lambda_k / lambda_plus
selected eigenvalue mass
effective rank
cos^2 theta
parent size
child balance
diagnostic status
```

The diagnostic remains descriptive unless it can estimate the conditional tail
\[
\Pr(R_\rho\ge r\mid X\in\mathcal S_\rho)
\]
with a declared support and precision contract.

## Evidence

- [[method-proof-web]] records the fixed-object projected-Wald proof and the
  selected-region proof gap.
- [[selected-hierarchy-selection-geometry]] records why same-data hierarchy
  selection creates selected barycentric contrasts.
- [[selected-hierarchy-geometric-law-map]] records the current edge-action,
  spectral, angular, parent-size, branch-length, and barycentric variables.
- `wiki/sources/edge-selection-null-audit-20260601.md` records that same-data
  hierarchy selection opens nearly all tested null edges, while fixed-tree
  permutations do not.
- `wiki/sources/selected-hierarchy-null-audit-20260601.md` records selected
  root and non-root null simulations with large selected-hierarchy ratios.
- `benchmarks/diagnostics/calibration/selected/hierarchy/selected_hierarchy_geometry_covariates.py`
  records the current diagnostic proxy variables.
- `tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_radial.py`
  implements the nuisance-derived Gaussian radial path and deterministic
  selected chi-radial interval integration without production wiring.
- `tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_hierarchy.py`
  replays squared-Euclidean average linkage and analytically constructs the
  complete observed merge-sequence interval union from quadratic constraints.
- `tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_diffusion.py`
  refits adaptive-pydiffmap or Hamming-neighbor geometry and replays average
  linkage at one radius without claiming an exact interval union.
- `tests/statistics/51_test_selected_gaussian_radial.py` protects the path
  invariants and analytic selected-tail examples.
- `tests/statistics/52_test_selected_gaussian_hierarchy.py` protects full
  merge-signature replay, every active competitor constraint, analytic
  boundaries, the two-leaf full region, and narrow disconnected cells.
- `tests/statistics/53_test_selected_gaussian_branch_lengths.py` protects fixed
  branch-time covariance scaling, fixed-topology NNLS merge-signature
  preservation, the cancelled-quadratic remote-interval regression, and the
  observed failure of fixed projection under recomputed linkage/NNLS branch
  lengths with `branch_length_state`.
- `tests/statistics/54_test_selected_gaussian_diffusion.py` protects adaptive
  diffusion-distance and topology changes across the continuous path, observed
  binary replay for both Hamming variants, and fail-closed rejection away from
  binary support.
- [[root-selected-region-margins-20260603]] records the first concrete replay
  of root merge-selection inequalities and their margins for representative
  benchmark contexts.

## Links

- [[method-proof-web]]
- [[selected-hierarchy-selection-geometry]]
- [[selected-hierarchy-geometric-law-map]]
- [[selected-hierarchy-null-support-contract]]
- [[selected-ratio-tail-law-diagnostic-20260602]]
- [[root-selected-region-margins-20260603]]
- [[open-mathematical-questions]]

## Open Questions

- For the current Hamming and Euclidean benchmark metrics, which constraints
  dominate the null-whitened distance from \(x\) to
  \(\partial\mathcal S_\rho\)?
- How should exact/tie-heavy Hamming cells be represented geometrically:
  tangent cones of polyhedral cells, random tie-breaking cells, or a discrete
  selected-region law?
- How does edge-path p-value action relate to the radial boundary distance
  for selected edge opening once the projection and Tree-BH selection cells
  are included?
- Does the tangent cone explain the selected-ratio tail better than the
  current edge-plus-spectral proxy equation?
- Can this root model be extended to non-root focal sibling contexts without
  making exact depth matching too sparse?
