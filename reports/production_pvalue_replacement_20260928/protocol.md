# Independent Gaussian replacement benchmark protocol

Frozen before evaluation: base seed 2026092811 (pilot used 2026092810); 250 null and 50 signal paired
datasets per covariance profile; 32 leaves; dimensions 2 (identity), 8 (AR1 .6),
and 40 (identity, high-dimensional). Discovery and inference use distinct
SeedSequence children. Leaf means coincide across the paired measurements;
noise is independent across all rows and between matrices. Covariance is known
from the generator, not estimated. Signal means differ by two in the first two
coordinates (balanced groups). All methods share the same paired datasets.

Infer seven linkage trees and neighbor joining with MAD rooting from discovery
measurements only. Fit positive PCA candidate spaces on discovery data after
known-covariance whitening. Freeze topology, directions, eigenvalue groups and
energy fraction .9. Compare nominal chi-square after energy selection against
the exact conditional energy law on the SAME independent test contrasts.
Reuse production Tree-BH, traversal-aligned sibling BH and TreeDecomposition.
No empirical inflation is estimated from testing data. Edge alpha .001; sibling
alpha .01, plus .005 and .05 for average linkage and neighbor joining. Root-tail
cutoffs .001, .005, .01, .05. Record raw node statistics/ranks/p-values, inferred
node generator nulls, actual reached false splits, partitions, ARI, input hashes,
matched method statuses, source hashes and runtimes. Pointwise Wilson intervals
use independent dataset root/partition counts; do not pool correlated methods.

Also run the existing single-matrix production dispatcher on discovery data for
average/complete linkage and neighbor joining at sibling alpha .01, explicitly
using standardized Euclidean distances. These comparisons change data quantity,
covariance and projection fitting as well as tails; they are not an isolated
p-value replacement comparison. Existing dispatcher unsupported/error outcomes
remain separate. The separate 420-dataset baseline covers all prior six feature
families and five methods, but does not validate the new Gaussian-only API on
Bernoulli/categorical inputs. Diffusion on binary inputs is outside this contract.

No benchmark outcome will be used to fit parameters or select equations. The
conditional formula is fixed before execution. The benchmark does not establish
universal calibration with estimated covariance, dependent paired observations,
or other distribution families. Known covariance and genuine independent
measurements are explicit required production inputs, not manufactured by IDs.

The implementation pilot (two null/one signal per profile) found that existing
centroid/median construction rejects nonmonotone merge heights. Those cases
are explicitly classified unsupported_tree, not successful null results. No
linkage-height repair or method removal is applied. Pilot outcomes do not tune
the equation or thresholds; pilot data are not part of the evaluation stream.
