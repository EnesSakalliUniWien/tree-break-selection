# Branch distributions and variance: frozen comparison protocol

This diagnostic changes no production defaults. It reuses the 900 Gaussian
dataset pairs from `independent_diffusion_nnls_20260928`: 250 null and 50 signal
replicates for identity2, ar1_8 and identity40, 32 leaves, evaluation seed
namespace 2026092813. A separate pilot uses 2026092820. These are matched
extensions of an existing evaluation, not a new independent confirmatory study.

## Current single-matrix pipeline

Cross adaptive Euclidean diffusion linkage / diffusion plus fixed-topology NNLS,
three spectral contexts (leaf_only, empirical_barycenter, branch_length_state),
and edge variance policies (none, normalized_branch_length). Run sibling alpha
.005/.01/.05, edge alpha .001, passthrough true. All other settings match the
previous comparison. Raw linkage variance is explicitly a diagnostic override.
Use the existing distance runner with the actual diffusion distances because
the adaptive diffusion dispatcher does not forward the internal-state options.
Verify resolved configuration and compare baseline rows to the saved evaluation.

The edge flag does not disable sibling branch variance. Sibling tests always
use 1 + (left_length + right_length)/(2 * positive_mean_length). Edge tests use
1 + edge_length/positive_mean_length only when the edge flag is enabled.
Do not describe the none setting as globally branch-free.

Save statuses, actual configurations, data/topology/branch hashes, labels,
node p-values, rank, spectrum, empirical-inflation metadata, independent row
counts and final partitions. Compare matched variants, including unavailable
outcomes. Null false partition means more than one final cluster. Signal ARI
and false splits whose child population means are equal are separate outcomes.

## Independent-equation compatibility diagnostic

Fit discovery-only diffusion/NNLS using the existing independent API. For
internal-row variants, use the existing branch-state helper or empirical
barycenters on discovery whitened by the externally known covariance. Refit
only the candidate PCA on the augmented discovery rows, preserving the
leaf-only candidate-rank cap. Internal rows are deterministic combinations,
not additional independent observations. Freeze the basis before inference.
This report-local adapter is experimental, not a new production flag. It
does not reproduce the current pipeline's estimated whitening or MP selection.
Run the existing conditional energy law and full correction/traversal at all
three sibling alphas. Record orthonormality and frozen-context checks.

Separately compute the uncorrected plug-in tail Fbar(T/m) for the sibling branch
variance multiplier m. This is an explicitly incompatible negative control
under the known independent-row covariance generator, not a production gate.
The correct scaled reference is Fbar(m * (T/m)) = Fbar(T). Record both tails and
check their equality after reference correction. Do not invent additional
independent observations or a branch-correlated generator to justify this scale.

## Assessment

Report root null tails at .001/.005/.01/.05 with Wilson intervals. Independent
conditional root tests form one prespecified 72-cell family (3 profiles x 2
branch sources x 3 contexts x 4 cutoffs); use Holm for excess-tail tests.
Current plug-in and incompatible scaled tails remain descriptive. Nodes within
a dataset are dependent; tail denominators count datasets, not pooled nodes.
250 null replicates have limited rare-tail resolution; neither zero false
partitions nor non-rejection proves calibration. No tuning after evaluation.
