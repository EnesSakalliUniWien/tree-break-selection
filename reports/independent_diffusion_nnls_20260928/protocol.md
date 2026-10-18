# Independent discovery diffusion/NNLS benchmark protocol

Evaluation seed namespace 2026092813; pilot namespace 2026092812. Freeze all
settings before evaluation. Generate 250 null and 50 signal pairs for each of
identity covariance in 2 dimensions, AR1(.6) covariance in 8 dimensions, and
identity covariance in 40 dimensions. There are 32 corresponding leaves with
independent Gaussian noise across rows and between discovery/inference matrices.
Signal leaf means are balanced groups differing by two in the first two
coordinates. Covariance is externally known from the generator.

Matched discovery geometries: known-covariance Euclidean average linkage;
adaptive pydiffmap Euclidean diffusion plus average linkage; the same diffusion
and linkage followed by fixed-topology NNLS on original standardized discovery
columns. Diffusion k=10, time=3, components=30, bandwidth=-1/(d+2), epsilon=median.
NNLS uses all 496 pairs (cap 50,000), seed zero, tolerance 1e-5, at most 1000
iterations; nonconverged fits fail explicitly. Save optimality and residuals,
not just the optimizer success flag. These are the existing continuous-Euclidean
capabilities, not a validation of the binary-Hamming preset.

Freeze discovery topology, branch lengths, diffusion coordinates and PCA before
seeing inference measurements. Candidate PCA is fitted to discovery measurements
whitened with known covariance. Energy fraction=.9; nominal versus conditional
energy tails share identical test contrasts. Edge alpha=.001; sibling alphas
.005,.01,.05, passthrough=True. Root-tail cutoffs .001,.005,.01,.05. No empirical
inflation or fitted branch-time covariance is applied under known independent
row covariance. Therefore NNLS alone is predicted to have no test/partition
effect when discovery topology is fixed; check this directly for every matched
case. This is a model implication, not a reason to omit the NNLS fit.

Also run the current same-matrix adaptive diffusion and adaptive diffusion NNLS
through their existing dispatcher, overriding only metric=euclidean for these
Gaussian inputs and using their default branch-time/inflation policies. These
whole-method comparisons change data quantity, covariance treatment and rank
selection as well as the reference law. Record actual configurations and do not
attribute their differences solely to the new equation.

Save dataset hashes, source hashes, solver/geometry metadata, before/after model
hashes, per-node statistics/ranks/p-values/truth, root tails, partitions, ARI and
reached false splits. Unsupported/failed/skipped cases remain explicit and are
excluded only from successful paired comparisons with denominators retained.
Pointwise Wilson intervals use independent datasets, not overlapping nodes or
repeated alpha variants. Check conditional root tails with one-sided binomial
tests and Holm correction over the planned 3 profiles x 3 geometries x 4 cutoffs
(36 comparisons). No threshold/equation tuning on evaluation outcomes. A finite
non-rejection does not prove calibration; rare-tail precision is limited.
