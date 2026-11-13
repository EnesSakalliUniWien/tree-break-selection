# Subtree spectral alignment and alternative clustering protocol

Confirmed scope: equal subspaces despite eigenvector signs or basis choices,
prescribed reflections, matching eigenvalues, diffusion and covariance matrices,
left/right eigenvectors, parent/child comparisons, and clusters represented by
different combinations of subspaces but obscured by one topology.

## Scope and provenance

This is a diagnostic study, not a production preset or a calibrated significance
test. Preserve existing working changes and previous reports. Use the existing
Python 3.12.11 environment, NumPy 2.3.4, SciPy 1.16.3 and pydiffmap 0.2.0.1;
single BLAS thread for experiments. No installation, commits, pushes or services.
The narrow existing geometry/candidate/feature-diffusion baseline passed 30 tests.
The repository has no .github directory; no remote CI result is asserted.

The shared geometry and diffusion source, eigenvalue grouping, covariance
spectral-context source, callers and tests were inspected before edits. Numerical
eigenvalue ties use the existing rule; this rule is not a test of population
eigenvalue equality. Exact fixtures and noisy finite panels are separate.

## Mathematical comparisons

Compare orthogonal projectors, rank, principal angles and directional containment
after specifying a common ambient coordinate system. Equal nonzero spectra do
not imply equal subspaces. Sign and within-block basis changes preserve a span.
Reflection comparisons require a declared orthogonal involution; an arbitrary
fitted transformation does not establish an intrinsic data symmetry.

For covariance, retain identical feature coordinates and units across parent and
child. Check the exact within-plus-between covariance decomposition. The native
Wald/MP covariance additionally uses node-specific whitening; a raw covariance
comparison must not be labeled as that native comparison.

For diffusion, capture the actual kernel passed to the existing production
coordinate function without altering that function. Preserve zero-diagonal and
maximum-symmetrization conventions. Examine unclipped eigenvalues of S and P,
both left and right eigenvectors, and the degree transformations linking them.
Keep the production retained/clipped/time-weighted coordinates separate.
For a child, distinguish the row-restricted parent eigenvectors, the induced
parent weight subgraph, and a newly fitted child kernel. Compare only after
restriction into the child's ordered row coordinates, recording rank loss.

## Finite experiments

Exact examples: rectangle, tied square, equal spectra with rotated eigenspaces,
prescribed reflection preserving versus changing a covariance/subspace, and an
irregular-degree graph with known left/right transformation and child restriction.

The alternatives panel defines 20 seeds in each of four cases: crossed means
(amplitudes 3,2), crossed means with equal population amplitudes (2.5,2.5),
hierarchy-compatible groups, and a common-mean Gaussian null. Each has 64 rows,
six features and independent Gaussian noise with standard deviation .5. The
alternatives protocol owns exact means, seeds and truth definitions; truth is
used only for evaluation. Examine every nonempty union of whole positive sample
covariance blocks, at most 63, using deterministic Ward reconstruction. Refit a
declared truth-free shortlist against a common all-pair Euclidean target.

Reuse those same 80 matrices for parent/child comparisons under the native
symmetric/average and Markov/Ward configurations. Inspect the root and its two
children when supported; retain small-child/kernel failures explicitly. Record
native clustering outcomes and actual spectral contexts separately from raw
covariance geometry. Never count nodes or repeats as independent datasets.

Record AMI, matched group recovery, clade/split representation, spectral ranks,
eigenvalues, projector differences, containment and kernel changes. Any maximum
over subspace recovery scores is a labeled descriptive oracle bound, not an
implementable truth-free method. Crossing partitions need not be representable
as clades of one rooted tree; their joint refinement is a different partition.

## Outputs and completion

Keep every attempt/status, source and data hashes, finite thresholds and original
outputs. Verify invariance and analytical identities, run focused tests, and
record findings in the existing wiki with raw evidence captures. A numerical
symmetry, eigenvalue match, recovered alternative or topology count does not
establish a p-value, true population clustering, or final-partition error control.
