# Frozen geometry and perturbation checks

This isolated geometry arm evaluates same-subtree topology candidates. It does
not fit TBS gates or test a clustering hypothesis. No outside rows enter these
diagnostics.

## Bounded candidate coverage and timing

Use fixed seed 2026110911 with a separate seed sequence per shape
(128,12), (240,12), (240,72), and (256,4). Draw iid standard-normal coordinates,
generate a current full-space Ward linkage and include it in the candidate
family. Evaluate one canonical Ward result per whole spectral block and leading
union. Deduplicate labeled trees. Run both full-subtree targets separately,
ordinary Euclidean distance and production mean squared standardized Euclidean
distance. Use every unordered leaf pair. Record candidate count, spectral span
count, elapsed time, selected versus current change, residual and fit failures.
No timing threshold is learned from these data. Use single-thread BLAS and
one worker.

The compressed path-incidence objective is checked against the existing dense
NNLS on all 945 rooted six-leaf candidates in each distance arm. Deterministic
checks cover rectangle, square, simplex, tied-block rotation, row relabeling,
root-edge sum, zero internal edges, and failure above the stated bounds.

## Six-leaf perturbation comparison

The fixed feature matrix has the six rows of a 2-by-3 grid. Its columns are
the two row-factor copies with values -3 and 3, two column-factor copies with
values -2, 0 and 2, and one alternating value -0.15 or 0.15. The baseline
geometric family is the union of all tied-Ward output families over the
unchanged whole spectral blocks and leading unions, using the existing
relative tie tolerance 1e-10. This is a geometric diagnostic family.

For each of five regimes use 2,000 draws from a separate stream with seed
sequence [2026110912, regime_index]. Recompute one canonical full-coordinate
Ward tree on every draw, rounding relative distances to 12 decimal places
as in the bounded geometry diagnostic. Count trees outside the baseline
family, their support frequencies, and the number of distinct observed
topologies. These are support/stability comparisons.

- Feature bootstrap samples five columns with replacement.
- Small and large Gaussian perturbations add independent entries with
  standard deviation 0.05 and 0.5 times the RMS centered entry, respectively.
  These are noise sensitivity scenarios around the observed geometry.
- Row control permutes rows and then restores stable leaf-label order.
- Gaussian orthogonal reference draws a Haar orientation in the
  centered row subspace, preserving feature Gram and mean. It changes row
  distances. It is a separate model-based perturbation; outside-family
  frequency here is a diagnostic, not a p-value.

Retain counts and source for rerun. No tune, retry, or filter by outcome.
