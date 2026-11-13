# Subspaces, parent–child alignment and clustering alternatives

The constructed examples and finite panels support the proposed mechanism:
different subspaces can expose different cluster partitions that one rooted
topology cannot represent simultaneously. They also distinguish omission of a
factor grouping from loss of the joint groups. No production configuration or
statistical gate changed in this study.

## What was implemented and checked

The new `subtree_spectral_alignment.py` diagnostic compares numerical ranks,
principal angles, directional containment and whole-block projectors. It checks
prescribed reflections, raw covariance decompositions, Markov left/right
eigenvectors and degree transport, and parent bases restricted to child rows.
Equal angles without equal rank are not reported as equal subspaces.

Three Sol agents covered alignment implementation, mathematical/Wolfram review,
and alternative-clustering experiments. The coordinator connected those results
to actual native pipeline replays. Existing geometry, metrics, branch fitting,
diffusion construction and tests were reused. The initial 30-test baseline
passed; the final focused run passed 48 tests. Source and input verification is
recorded separately from statistical claims.

## Equality, reflections and eigenvalues

For orthonormal bases Q_A and Q_B in the same ambient coordinates,

\[
\cos\theta_j=\sigma_j(Q_A^TQ_B),\qquad
\|Q_AQ_A^T-Q_BQ_B^T\|_F^2=r_A+r_B-2\sum_j\cos^2\theta_j.
\]

Sign changes and orthogonal basis changes inside an exact repeated eigenspace
preserve its projector. Compare the entire block rather than individual
eigenvectors. Numerical near-ties use the existing grouping convention and do
not establish population multiplicity. Unequal weights within a near-tied block
can still make a basis rotation change weighted embedding distances.

The covariance matrices diag(4,1) and diag(1,4) have identical eigenvalues but
orthogonal leading eigenspaces. The prescribed feature-swap reflection maps
one covariance into the other; a sign flip of the first feature does not.
Allowing an arbitrary fitted orthogonal map would make equal-rank subspace
alignment automatic. The diagnostic therefore keeps direct comparison and
comparison under a declared reflection separate. It probes the six feature
sign reflections and the first-two-feature swap in the finite panel; it does
not claim to enumerate every possible ambient reflection.

## Left and right eigenvectors

For the production-symmetrized nonnegative kernel W with positive degrees D,

\[
S=D^{-1/2}WD^{-1/2},\quad P=D^{-1}W,\quad
R=D^{-1/2}U,\quad L=D^{1/2}U,
\]

where S U=U Lambda. Consequently P R=R Lambda, P^T L=L Lambda,
L^T R=I, and both bases agree with U after the corresponding degree transport.
This is an exact relationship. Their unweighted Euclidean distances need not
agree because degree transport rescales rows.

Across the 160 native replays, the maximum degree-transported block projector
gap was about 1.65e-15. Median per-dataset maximum raw left/right block gaps were
.146–.166 in the signal regimes and .392 in the null.
The full unclipped spectrum is retained separately from production mode
selection, nonnegative clipping and diffusion time. Reconstructing the actual
production distances had maximum absolute error 1.39e-14.

For a real covariance matrix, left and right operator eigenvectors may be the
same orthonormal vectors. The left/right singular vectors of a rectangular data
matrix are different objects, indexed by observations/features. The report does
not conflate these comparisons. Exact Wolfram requests and primary sources are
in `math/wolfram_evidence.json` and `math/primary_sources.md`.

## Parents and children

In fixed feature coordinates with sample covariances,

\[
(n-1)C_P=(n_L-1)C_L+(n_R-1)C_R+
\frac{n_Ln_R}{n}(\mu_L-\mu_R)(\mu_L-\mu_R)^T.
\]

The last term can dominate a parent direction even when it is absent within
each child. Therefore matching parent and child eigenvalues is not required for
a valid hierarchy. The numerical decomposition error was at most 7.15e-15.
The raw-covariance comparisons use common features and units; actual native
Wald/MP contexts are captured separately because they use node-specific
null-whitened coordinates. No untransported native direction comparison is
claimed.

The 80 shared matrices were each replayed under symmetric/average and
Markov/Ward, giving 160 native fits and 320 child comparisons. Every native and
both induced/refitted child graph attempts returned `ok` in this panel.
The primary parent–child comparisons use unweighted eigenvectors from the
selected, unclipped spectral windows. Restricting that parent Markov-right span
lost numerical rank in 239/320 children. The induced parent-weight graph and newly fitted child graph
differed in symmetric-operator Frobenius norm by more than 1e-10 in 316/320
children. The observed range was 2.09e-16 to .8165. These counts measure distinct
geometric operations; they are not independent-data rejection rates.

A separate check uses the actual clipped, diffusion-time-three coordinates on
the same saved memberships, with no new native fits. Nonpositive modes receive
zero weight in 154/320 induced and 160/320 refitted child windows; none of the
parent windows has a zero weight. The inherited and child production-coordinate
spans still differ in all 320 comparisons under the stated numerical tolerance.
The median projector gap is about 1.414 in both operator arms. Reconstructed
child-coordinate distances agree with the native geometry. The distinction
between an unweighted spectral window and a surviving production span is
retained in `parent_child/production_span_summary.json`.

Representative traces preserve the actual W, S, P and raw covariance matrices,
the native spectra/projections, branches and gate/traversal records. Native
recovery remained high on these simple signals: joint AMI was 1 except the
Markov tied-population mean .9968. Both configurations rejected the null root
in 20/20 datasets and returned mean 3.60/4.55 groups for symmetric/Markov.
Recovering geometry has not repaired the existing null calibration failure.

## Alternative partitions and lost representation

Each of the 80 six-feature matrices generated 63 nonempty unions of its six
positive numerical covariance blocks. Together with the full-coordinate
control this gives 5,120 candidate rows. One canonical Ward tree is constructed
per union; exhaustive tie exploration is limited to the exact four-leaf examples.

| Regime, 20 datasets each | Full tree has A | Full tree has B | Spectral family has both | One family tree has both |
| --- | ---: | ---: | ---: | ---: |
| Crossed means, amplitudes 3 and 2 | 20 | 0 | 20 | 0 |
| Crossed means, equal population amplitudes 2.5 | 13 | 7 | 11 | 0 |
| Compatible hierarchy, nested B | 20 | 20 | 20 | 20 |

In the crossed example, the same four groups support AB|CD and AC|BD. All four
intersections are nonempty, so both splits cannot be displayed by one tree.
The full tree still recovered the four joint clades in all 20 crossed datasets.
That is an omitted factor grouping, not destruction of the underlying data.
There is a separate loss from projection: in the compatible panel, the leading
rank-two span retained all four joint clades in only 6/20 datasets, compared
with 20/20 for full coordinates.

Selecting by minimum common-target NNLS residual over the declared shortlist
improved residuals in 6/20 crossed, 7/20 tied and 1/20 compatible datasets but
preserved the full tree's factor-coverage counts. Better fitting did not recover
the omitted second crossed factor. Conversely, all 20 common-mean null datasets
also generated 63 distinct Ward topologies. A topology count is not a population
count or a correction for selection.

The tied population does not imply an exactly repeated sample eigenvalue. The
exact square has one intrinsic tied rank-two block, whose Ward tie family has
two crossing outputs with equal additive RMSE .3382. Prescribing separate axes
inside that block adds external structure; the spectrum alone does not choose
those axes. The retrospective best AMI across unions is labeled an oracle
bound and is never used to select a reported deployable method.

## Consequences for the production research plan

1. Record subtree alternative partitions and their compatibility, alongside
   eigenblocks, distances and fitted branches. Evaluate the grouping factors
   and their joint refinement separately; a single reference partition cannot
   describe every meaningful grouping relation.
2. Keep inherited parent coordinates, induced child operators and refitted
   child operators distinguishable. Compare ranks and projectors only after
   the appropriate row restriction, degree transport or feature-coordinate map.
3. Measure whether rank truncation removes a group distinction before choosing
   a faster candidate policy. Repeated eigenblocks cannot supply uniquely
   identifiable individual directions without additional information.
4. Calibrate the complete selection of spectral blocks/unions, topology, branch
   lengths and cuts. Neither a low residual, a symmetry nor many alternatives
   supplies the missing selective reference. Previous promotion gates remain.

These are measured diagnostics and requirements for a next candidate. They are
not an implemented multi-partition production output or a valid new p-value.
The finite panel is small, Gaussian, balanced and low dimensional; it does not
establish behavior on CHOIR or general population identifiability.

## Reproduction and verification

`alternatives/protocol.md` and `protocol.md` specify the cases and comparisons.
`summary.json` verifies all 80 generated input hashes against both panels and
records the scalar results. `parent_child/rows.jsonl.gz` retains all 160 records.
The source manifests were stable during the native panel. The alternatives
script's exact executed version is archived; the current version differs only
in imports, with an unchanged non-import AST. No original outputs were replaced.

Run `python -m reports.subtree_spectral_alignment_20261109.run --output NEW_PATH`
with one BLAS thread for native replay; existing manifests cannot be overwritten.
Run the adjacent `summarize.py` as a module to reconstruct the summaries.
`production_spans.py` independently checks clipping and production-coordinate
spans using the same inputs and saved child memberships.
`figures/plot.py` produces PNG/PDF/SVG. Focused test, style, artifact and wiki
verification is retained in `verification.json`.
