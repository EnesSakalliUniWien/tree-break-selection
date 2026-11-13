# Spectral alignment and topology compatibility

Date: 2026-11-09. This review distinguishes operator identities, descriptive
agreement in fixed coordinates, and evidence for an inferential claim.
The live sources are `tree_break_selection/space_separation/diffusion.py`,
`tree_break_selection/eigenvalue_groups.py`, and the node covariance and
spectral estimator under `tree_break_selection/tree/distributions.py` and
`tree_break_selection/hierarchy_analysis/statistics/projection/spectral/`.
Primary source locations and exact Wolfram checks are stored alongside this
note. No implementation or production default changed.

## Symmetric diffusion and Markov left/right vectors

Fix the same symmetric nonnegative weight matrix W on the same ordered
leaves, with strictly positive degrees D. Define
\[
P=D^{-1}W,\qquad S=D^{-1/2}WD^{-1/2}=D^{1/2} P D^{-1/2}.
\]
For an orthonormal eigenbasis U of S, the corresponding right eigenvectors
of P are R=D^{-1/2}U. Its left eigenvectors, written as columns, are
L=D^{1/2}U=DR. In these normalizations LᵀR=I,
RᵀDR=I, and LᵀD⁻¹L=I. A Markov spectral projector is R_B L_Bᵀ,
the similarity transport of the symmetric block projector U_B U_Bᵀ.
Right and left agree after transport to the symmetric coordinates.
Raw Euclidean vector, projector or distance equality is a separate claim.
The stationary weights are proportional to D. Using stationary weights
instead of D rescales right and left vectors reciprocally by a common
constant. [Coifman and Lafon, Appendix A](https://doi.org/10.1016/j.acha.2006.04.006)

The live code first clears the diagonal and symmetrizes by maximum, then
uses a zero-degree floor. With isolated rows the floor is a numerical
extension and an isolated row of D⁻¹W is still zero. The usual stochastic
and positive stationary measure interpretation requires positive degrees
on the analyzed support. The similarity itself uses the same degree
diagonal in both expressions. The live Markov coordinates are right
vectors, with the same mode selection and weights as the symmetric
coordinates. It retains complete numerical eigenvalue groups at the
rank boundary, removes a leading mode only if isolated, and clips
negative retained eigenvalues before applying diffusion time. For
positive time, a mode clipped to zero contributes no embedding
direction. These policies matter when comparing actual embeddings to
unweighted eigenspaces or the exact full-spectrum diffusion distance.

For a positive nontrivial mode, take the four-node unweighted path.
P has eigenvalue 1/2 with right vector (2,1,−1,−2)/sqrt(12).
At time 1 the contribution of this mode to squared distance between the
first two leaves is (3−2 sqrt(2))/24 in symmetric coordinates and 1/48
in Markov coordinates. Thus identical eigenvalues and exactly transported
modes need not give the same Euclidean geometry. A global orthogonal
change of embedding coordinates preserves all distances, whereas degree
transport is a row scaling. On a regular graph it reduces to one global
positive scale.

## Common coordinates and a child graph

Restricting a parent's diffusion coordinates to the child's leaf IDs
is an inherited embedding. It solves neither the induced nor the
refitted child eigenproblem in general. For a parent eigenvector u,
\[
S_{CC} u_C = \lambda u_C - S_{C\bar C}u_{\bar C}.
\]
Even with no leakage term, the induced operator generally renormalizes
with degrees computed inside C. A refitted child graph can additionally
change neighbor support, variable bandwidth, kernel scale, and weights.
Even if C is a disconnected parent component and all its degrees
match, parent and child rank selection need not retain the same modes.
On the four-node path above, restricting the right vector to the first
three leaves gives an induced-child residual
P_C r_C−(1/2)r_C=(0,0,sqrt(3)/4). This difference persists with the
same leaf identities and no kernel refit.

A defensible comparison states the target explicitly. For two operators
on the same graph, use their known degree transport. For inherited vs
induced/refitted geometry, match and order the shared leaf IDs, state
whether a common metric or a degree-coordinate transport is applied,
then recompute the numerical rank and orthonormalize the restricted
basis. Restriction need not preserve dimension or orthonormality.
A child basis is not automatically the intersection of the parent
eigenspace and the child row space.

## Covariance, projectors and reflections

A real covariance matrix C is symmetric positive semidefinite, so its
left and right operator eigenvectors can be the same orthonormal
vectors. The left and right singular vectors of a rectangular
centered data matrix are different objects, indexed by observations
and features respectively. The feature vectors from its right singular
space diagonalize feature covariance. Equal covariance eigenvalues
alone do not identify the eigenspaces in a common feature system.
For example diag(4,1) and diag(1,4) have identical spectra and
orthogonal leading eigenspaces.

In shared Euclidean coordinates, write Q_A and Q_B for orthonormal
bases after accounting for rank. Compare the projectors
Q_A Q_Aᵀ and Q_B Q_Bᵀ. The cosines of the principal angles are the
singular values of Q_AᵀQ_B. The squared Frobenius projector gap is
r_A+r_B−2 sum cos²(theta). All reported angles can be zero for
nested spaces of unequal dimension, so equality also requires
equal rank. For equal ranks, a simple isolated eigenline is
unchanged by a sign flip. A repeated block permits any orthogonal
basis change within the block. Compare the entire block projector,
and document the numerical tie rule and any excluded or unmatched
directions. A near-tie grouping is an analysis convention, not an
assertion of exact eigenvalue multiplicity. If embedding weights
are unequal within a numerical group, rotating the basis while
leaving the diagonal weights fixed does not generally preserve
the weighted geometry. For exact multiplicity the weights are
equal and the rotation commutes with them.
[Björck and Golub, Theorem 1](https://www.ams.org/mcom/1973-27-123/S0025-5718-1973-0348991-3/S0025-5718-1973-0348991-3.pdf)

A prescribed ambient reflection H=I−2nnᵀ for a fixed unit normal n
changes the projector to H Q_A Q_Aᵀ Hᵀ. Report direct and mapped
comparisons separately. Central inversion H=−I leaves every linear
subspace projector unchanged. Choosing a reflection after seeing two
unit vectors can always map one vector to the other when they differ:
take the normal parallel to their difference. Allowing an arbitrary
ambient orthogonal map sends any r-dimensional subspace to any other
of the same dimension. Exact orthogonal conjugacy of real symmetric
matrices with matching spectra is likewise automatic when the ambient
map is unrestricted. In the diagonal example, flipping the first
feature leaves the original covariance unchanged, whereas reflecting
across the feature-swap hyperplane maps it to the second matrix.
The first map and the second map answer different questions.

An orthogonal alignment *within columns* of two embeddings is a different
operation. It preserves the matched-row Gram and every pair distance in
each embedding, and cannot alter its column space. A zero orthogonal
Procrustes residual for two equally shaped matched-row embeddings
expresses equality up to a global coordinate rotation/reflection,
equivalently equality of their row Grams. If translations are allowed
first, apply the statement to their centered row Grams. An ambient
map chosen to force subspace equality cannot be substituted for this
test of matched-leaf geometry.

## Parent and sibling covariance

Use exactly the same ordered feature variables and units for raw
parent, left and right covariances. Let T_v be centered scatter
\(X_v^\top H_v X_v\), with H_v the centering matrix, and let
\(\delta=\mu_L-\mu_R\), \(n=n_L+n_R\). For a partition into the two
children,
\[
T_P=T_L+T_R+\frac{n_L n_R}{n}\delta\delta^\top.
\]
Expand each parent-centered row as its child-centered row plus the
child mean minus the parent mean. The cross terms sum to zero in each
child. Divide by n−1 to compare unbiased sample covariances, which
weights child covariance matrices by n_L−1 and n_R−1. Scatter also
covers singleton children with zero within scatter. The corresponding
population law adds within covariance weighted by the mixture shares
and the between term \(\pi(1-\pi)\delta\delta^\top\).
The between term is positive semidefinite and rank at most one.
Adding it can rotate the parent directions or change their rank.
Equality, inclusion, and preservation of leading order require
additional conditions. Even equality of the two child eigenspaces
alone does not remove the between term.

This identity refers to raw total covariance. The guarded sibling
covariance pools within-child covariance when its support conditions
hold, with parent-total fallback. The native node spectral estimator
decomposes covariance of its node-specific null-whitened tangent
rows. A parent/child comparison of those returned directions needs a
common coordinate system or an explicit transport of the loadings.
For a linear whitened coordinate z=A(x−mu), a projected linear
functional qᵀz has raw covector Aᵀq. Transport of a displacement
vector instead uses A⁻¹ where defined. Direct comparison of raw
feature PCA, node-local whitened PCA, and sample-space diffusion
vectors therefore compares different objects.

## Incompatible alternatives and statistical scope

Rooted clades of one tree form a laminar family. Two nodes are either
ancestor/descendant, so their clades are nested, or lie in disjoint
branches. Crossing non-nested clades cannot both occur. For unrooted
binary splits on the same leaf set, a conflict requires all four
intersections of the two sides to be nonempty. To see necessity,
remove the first edge. The second lies in one component, so one side
of its split is contained in that component and disjoint from the
other. For two splits with one empty intersection, the surviving
parts can be assembled on a common two-edge backbone, so the
criterion is also sufficient for the pair.

A four-corner rectangle, with horizontal and vertical coordinates,
has distinct axis eigendirections if the side lengths differ. The
two axis cuts give the two balanced side pairings. All four
intersections contain one leaf, so one tree cannot display both
cuts. On a square the covariance block is repeated. Rotating the
chosen basis inside that block can change the coordinate-axis cuts,
while the whole plane and its complete affine cut family are
unchanged. Whole-block and cumulative-span cut families avoid
turning arbitrary axes into scientific alternatives. Geometry can
enumerate or prioritize competing splits, including some omitted by
one fitted topology. Matching eigenvalues, projector gaps,
alternative counts, and rank- or reflection-selected matches are
descriptive. None is a selected-null p-value. A reference would have
to account for the complete data-dependent comparison and selection.

The Wolfram check verifies the reversible similarity and three
eigen-equations, the biorthogonal normalization, both distances in
the four-path example, its child residual, the two prescribed
reflections, a basis rotation of a repeated covariance, and the
exact within-plus-between scatter identity on one unequal-size
fixture. The exact query and output are in `wolfram_evidence.json`.

## Diagnostic implementation cross-check

I checked the new `subtree_spectral_alignment.py` contracts against
these derivations. Its dense, unclipped reversible operator diagnostic
requires positive degrees, while the production embedding retains its
own mode policy. Its reflection input is a prescribed orthogonal
involution, which includes both plane reflections and central
inversion. Its directional residuals distinguish inclusion from
equal span and yield the Frobenius projector gap without constructing
dense projectors. Its parent restriction recomputes rank relative to an
orthonormal parent basis and compares against the supplied child span.
The 20 focused diagnostic tests passed. This cross-check is mathematical
and numerical on the covered fixtures, with no selected-inference claim.
