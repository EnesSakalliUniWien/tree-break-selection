# Spectral alignment in a subtree

## Equality of spans and prescribed reflections

The diagnostic represents bases by columns in a declared common ambient space.
A rank-revealing SVD removes redundant columns. If (Q_A,Q_B) are the resulting
orthonormal bases, the singular values of (Q_A^TQ_B) are cosines of the principal
angles. The number of such angles is the smaller rank, so zero angles do not
establish equality when the ranks differ. In particular a line in a plane
has one angle zero, but the projector distance is one. Zero-rank inputs
have an empty angle vector, with containment reported explicitly.

The two directional residual norms are
(|(I-Q_BQ_B^T)Q_A|_F) and (|(I-Q_AQ_A^T)Q_B|_F).
Their squares sum to the squared Frobenius distance of the projectors.
Equal numerical rank and small residuals in both directions define the
reported numerical span match. Sign flips, permutations, and orthogonal
changes of basis within an exact eigenspace preserve its projector.

Covariance comparisons operate in the same **raw feature coordinates**.
Matching eigenvalues, even with multiplicities, establish spectral agreement
only. For example, covariance diagonals proportional to ((4,4,1)) and
((1,4,4)) have the same eigenvalues. The leading planes share one direction,
have angles zero and pi/2, and have projector distance sqrt(2). Numerical
eigenvalue groups follow the project's adjacent-tie convention and do not
assert exact population multiplicity. Native spectral context can use
node-specific whitening and tangent coordinates. Comparing those contexts
requires explicit coordinate alignment, even when their ranks agree.

For a prescribed reflection (F=F^T=F^{-1}), the separate checks compare
(FQ_A) with (Q_B) and (F C_A F^T) with (C_B). The line through
((1,1)) is changed by reflecting the first coordinate, while the line
through ((1,0)) is preserved. Under that same reflection a covariance
with nonzero cross-covariance changes the sign of its off-diagonal entries.
A paired observed span or covariance can match its prescribed reflected
counterpart without either object being invariant by itself. No rotation or
reflection is estimated to manufacture a match.

## Parent covariance and restricted parent row spans

With raw coordinates (x_i), group sizes (n_g), group means (ar x_g)
and parent mean (ar x), write
[
 (n-d)C =
 sum_gsum_{iin g} (x_i-ar x_g)(x_i-ar x_g)^T +
 sum_g n_g(ar x_g-ar x)(ar x_g-ar x)^T.
]
The parent divisor is (n-d) in both terms, with (d=0) for population-form
empirical moments or (d=1) for sample covariance. A child of size one
contributes zero within scatter. For non-singleton child sample covariances
the within numerator weights each child by (n_g-1), and the between
component need not share its leading eigenspace with any child covariance.
For a population law, the same decomposition is the conditional covariance
identity with probability weights.

For a row-based diffusion eigenspace, first restrict the orthonormalized
parent basis to the child row IDs in the declared order, then take another
rank-revealing SVD. The reported rank loss is relative to the original
parent span. This operation maps rows, whereas covariance feature
directions live in a shared feature space. A child diffusion graph rebuilt
on those rows can differ because edges or degrees change. Even restricting
an existing weight matrix and renormalizing its degrees gives a new
operator. Agreement or disagreement between the resulting spans is a
diagnostic about those specified constructions.

## Reversible diffusion

For finite nonnegative weights the diagnostic removes self-loops and uses
the production elementwise-maximum symmetrization. Its supported model
requires positive degree at every row and uses a dense decomposition up
to the production dense size limit. Unlike production it rejects isolates
rather than inserting the degree floor. If (D) is the degree matrix,
[
 S=D^{-1/2}WD^{-1/2},qquad P=D^{-1}W,qquad
 P=D^{-1/2}S D^{1/2}.
]
From (SU=ULambda) with (U^TU=I), the Markov right vectors are
(R=D^{-1/2}U), and the left vectors stored as columns are
(L=D^{1/2}U). Thus (PR=RLambda), (P^TL=LLambda), and (L^TR=I).
The degree-transformed bases (D^{1/2}R) and (D^{-1/2}L) match the
symmetric span block by block. Without those transforms, the Euclidean
spans can differ for unequal degrees. The diagonal transform is not
generally an isometry for row distances, so equivalence of the operator
spectrum does not imply identical clustering distances. A whole-space
comparison of full-rank bases is uninformative, hence the diagnostic
compares complete numerical eigenvalue blocks individually.

This is the full unclipped operator spectrum, including zero and negative
eigenvalues. Production can drop an isolated leading mode, keep numerical
ties together at its chosen rank, clip negative eigenvalues, and weight
the retained vectors by diffusion time. A test fixture checks the span of
both production embeddings after those steps. A weighted star fixture
keeps its two-dimensional numerical zero block and its negative mode.
Neither signed operator eigenvalues nor clipped coordinates establish
that a subtree contains more than one population.

## Topology interpretation and checks

Two overlapping views may preserve different crossing cuts. On the square,
the two coordinate axes yield the crossing root splits AC|BD and AB|CD.
The full tied plane is preserved under basis rotation. Examining only one
basis axis or the realized binary hierarchy can hide the other cut.
Compare cuts and candidate trees using an explicit distance and linkage rule
on the same labels. A projector match, a reflected match, or similar
eigenvalues do not establish validity of a topology or a p-value.

The irregular-degree fixture has degrees (4, 4, 6, 2). Its leading
symmetric/right projector distance is 0.259938, and the right/left distance
is 0.471405. The corresponding symmetric/degree-transformed distances are
about 2.48e-16. All operator and duality residuals in this fixture are
below 4e-15, and both production embeddings passed the span check. A separate test obtains a
different leading symmetric span after rebuilding a child graph, while
its leading right Markov span remains constant. The focused new suite
covers sign changes, within-block basis rotation, unequal and zero rank,
matching eigenvalues with different eigenspaces, covariance decomposition,
reflection, and row restriction with rank loss. No production defaults
or gate policies are changed.

The new 22 checks and five existing feature-family diffusion checks passed
together under Python 3.12.11 and single-thread BLAS. Ruff passed.
