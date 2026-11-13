# Subtree candidate geometry and fit checks

The bounded candidate selector keeps the current rooted topology and proposes
one canonical Ward topology per positive whole eigenvalue block and leading
union of blocks. It deduplicates trees by their labeled unranked topology. The
spectral candidate computations use centered SVD without materializing a row
projector for every span. Whole-block distances are invariant to basis rotation
within an exact tied block. The deterministic Ward computation uses relative
distance rounding to 12 decimal places and stable leaf-label order. This
resolves numerical ties by a declared convention and does not exhaust all
tied-Ward alternatives. The existing exhaustive cut and topology limits remain
10 and 6 leaves, respectively. The new bounded selector supports 2–256 leaves
and fails explicitly outside this range. Its optional exhaustive arm includes
all rooted binary topologies only at up to six leaves.

Every candidate is fit to the same full-subtree target and every unordered
leaf pair. The two arms are ordinary Euclidean distance and production mean
squared standardized Euclidean distance, including constant features in the
standardized denominator. Compare residuals only within one arm. The selector
retains the current topology when its RMSE is within 1e-10 times the maximum
target of the minimum, and otherwise uses canonical order for that tie. The
returned selected linkage preserves the original current linkage when supplied,
or the candidate's Ward linkage and merge heights in the chosen spectral span.
An exhaustive-only candidate has synthetic monotone heights that represent
its topology, with no interpretation as observed merge heights or fitted
branch times. Gate or branch-time evaluation needs to refit separately.

## Compressed all-pair NNLS objective

For each edge split j, let M_ij mark whether leaf i is on its canonical side.
An incidence row in A marks the splits that separate that leaf pair. For
splits j and k, let n_ab count the leaves with side indicators a and b.
Then (A'A)_jk = n_11 n_00 + n_10 n_01, counting the unordered pairs crossing
both splits. The jth component of A'b is the sum of the supplied target
distances between the two sides of split j. Removing the duplicate pair
incidence columns of the complementary root edges leaves a nonsingular
split Gram for a binary topology. The coefficients of a represented
additive distance are unique on those splits. Pendant-edge coefficients
follow from a three-leaf difference, and internal-edge coefficients follow
from a quartet that takes one leaf in each of the four incident subtrees.

If R'R=A'A is its Cholesky factorization, minimizing ||R l−R'^(-1)A'b||² for
l≥0 differs from the original all-pair objective only by a constant. The
solver uses SciPy NNLS in that compressed system and reports RMSE and
maximum residual back in the input target units. It stores no
pair-by-edge matrix. This solves the same objective as the dense fit up to
floating arithmetic, with no sampled pair set. Both target arms agreed
with the existing dense NNLS on all 945 rooted six-leaf trees in the
deterministic fixture.

A root split has only a fitted sum of its two root-edge lengths. Rerooted
topologies with the same unrooted splits share the entire fit. Different
unrooted topologies can also be equivalent as fitted metrics if internal
edges collapse to zero. The result records zero internal edges and flags
fitted-metric equivalence to the selected topology with the stated numerical
tolerance. An RMSE tie alone does not imply equality of the fitted metrics.

## Bounded coverage and runtime

All eight frozen runs completed on single-thread BLAS. The figures below are
one measured runtime per shape and arm and exclude downstream gates.

| Leaves | Features | Unique candidates | Pair count | Euclidean seconds | Standardized squared seconds |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 128 | 12 | 23 | 8,128 | 0.400 | 0.373 |
| 240 | 12 | 23 | 28,680 | 1.850 | 1.734 |
| 240 | 72 | 143 | 28,680 | 11.259 | 11.127 |
| 256 | 4 | 7 | 32,640 | 0.664 | 0.554 |

The current topology was retained in two Euclidean runs and was changed in
all four standardized-squared runs. Changes concern the declared residual
alone. At 240 leaves and 72 features, Euclidean RMSE changed from 0.677202
to 0.676911 and standardized-squared RMSE changed from 0.225963 to 0.225850.
A separate 240-by-12 profile attributed 1.48 of 1.85 seconds to the compiled
NNLS routine. High-rank candidate search is therefore a material runtime
cost. The deterministic fixture and coverage panel do not establish a
production runtime ratio or clustering accuracy improvement.

## Six-leaf perturbation check

The geometric family contains eight distinct trees. It is the union of
tied-Ward output families over the observed whole blocks and leading unions
under the declared relative tolerance 1e-10. The bootstrap
and noise scenarios use one canonical full-coordinate Ward result per
draw. The Gaussian orbit changes row orientation while preserving the
observed feature Gram and mean. These arms describe support and sensitivity.
The results are not p-values or estimates of which hierarchy is true.

| Perturbation | Draws | Distinct trees | Draws outside observed geometric family |
| --- | ---: | ---: | ---: |
| Feature bootstrap | 2,000 | 5 | 2 |
| Small Gaussian noise | 2,000 | 4 | 0 |
| Large Gaussian noise | 2,000 | 28 | 111 |
| Canonical row-order control | 2,000 | 1 | 0 |
| Gaussian orthogonal orbit | 2,000 | 665 | 1,956 |

The fixed protocol, source, and machine-readable outputs are beside this
report. There were no failed coverage fits. Focused geometry and candidate
checks passed 16 tests, including square, rectangle, simplex, rotation,
row-label relabeling, exact small fits and the declared bounds. These
experiments are confined to the selected leaves and do not establish
performance of replacing the topology within the full TBS procedure.
