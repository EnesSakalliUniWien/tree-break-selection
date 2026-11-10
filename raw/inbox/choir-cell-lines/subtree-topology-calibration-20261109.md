# Subtree Gaussian reference calibration

Date: 2026-11-09, Europe/Vienna. The opt-in core and its tests are in
`benchmarks/diagnostics/calibration/selected/hierarchy/subtree_gaussian_reference.py`
and `tests/statistics/71_test_subtree_gaussian_reference.py`.
Runnable protocols, complete attempts, summaries and source hashes are in
`results/analyses/subtree_topology_reference_20261109/calibration/`.

## Conditional target

For a specified subtree with m independent common-mean Gaussian rows and
one unknown common covariance, fix its column mean, centered feature Gram,
and all outside rows. In Helmert coordinates the centered matrix has m−1
rows, and its conditional law is invariant to left orthogonal rotations.
Conditional on the Gram, matrices with that Gram form one orbit. A
caller-supplied selection event restricts the orbit. Selected membership
requires replay of the actual selection and any score-defining choices
that were selected. Membership alone fixes the Helmert basis through
canonical row-ID order.

The chain proposes a symmetric random orthogonal transform on a coordinate
pair, including reflections, and stays in place when a proposal leaves the
event. It uses an independently uniform observed position and counts all
upper-tail ties. This is a stationary reversible random-position rank
argument. Movement is a power and resolution diagnostic. Failures give an
unresolved result, and event membership is evaluated numerically. Outside
rows are unchanged. The common-mean null is stronger than an equal aggregate
contrast under heterogeneous row means.

For a fixed contrast, the projection onto a Gram-determined rank-k left
subspace has conditional Beta(k/2,(m−1−k)/2) law. Rank zero, a full centered
span and a cutoff through a numerically tied Gram eigenspace are unresolved.
A Bonferroni factor of 2^(m−1)−1 covers all unordered nonempty binary labels
on a prespecified subtree. It does not cover selecting the subtree.

## Fixed, selected and full-span checks

The primary panels used four disjoint eight-row candidate subtrees,
a known simulation covariance with condition number 16, and 119 steps.
The reference itself used no covariance parameter. At nominal .05:

| Panel | Naive chain | Event-conditioned chain | Analytic reference |
| --- | ---: | ---: | --- |
| Fixed membership null, p=3 | 56/1200 | | One-feature Beta 55/1200, all-feature rank-3 Beta 59/1200 |
| Selected winner null, p=3 | 219/1200 | 60/1200 | Naive Beta 228/1200, selected truncated Beta 59/1200 |
| Full-span fixed null, p=12 | 14/300 | | Full-span Beta unresolved 300/300, explicit rank-2 Beta 12/300 |
| Selected signal, p=3 | 149/200 | 68/200 | Exact selected truncated Beta 89/200 |

The selected-null chain exact 95% interval was [.0384,.0639]. Its median
was 79 accepted moves and 2.29 ms per dataset. There were no unresolved
chain attempts. The true signal subtree won in 183/200 datasets and was
both selected and rejected in 67/200 by the conditional chain and 88/200
by the exact selected reference. The short budget and choice of statistic
matter for power. The p=12 run used the same first-feature score in the
full-row-span orbit. The all-feature fixed-rank Beta degeneracy is a
separate reference issue.

## All-label maximum and a selected descendant

In an independent p=2 Gaussian panel, each dataset had 16 rows. The score
was the maximum squared projection over all 127 unordered labels on an
eight-row subtree. The chain used 199 steps.

- Fixed membership. At .05 the unadjusted winning-label Beta rejected
  320/320, all-label Bonferroni rejected 14/320, and the max-score chain
  rejected 13/320, exact 95% interval [.0218,.0685].
- Selected lower child of a predeclared root median split. The conditional
  event required the same eight rows to remain the lower child, with
  outside rows fixed. At .05 the naive max-score chain rejected 11/160
  and the event-conditioned chain rejected 8/160, exact 95% interval
  [.0218,.0961]. At .01 the counts were 0/160 and 1/160, respectively.
  The conditioned chain had median 143 accepted moves and 5.48 ms.
  No chain was unresolved.

An initial descendant run supplied the selected rows in observed feature
order. That order would make the proposal basis depend on additional data
beyond the event. The core now orders row IDs canonically. The initial
extension output and code are preserved, and the corrected extension ran
with the same seeds and budget. The four main panels already used canonical
row order. Their first replicate each reproduced exactly against the final
core, excluding elapsed time. The main runner's initial JSON summary write
failed on a NumPy scalar after all individual attempts were saved. The
summary was regenerated from those records without resimulation.

Twenty-six focused core tests pass. Forty-five existing root and
descendant sampler tests passed as a baseline. Ruff and diff checks passed.
These panels concern the stated Gaussian null and simple fixed or replayed
membership events. They do not validate native topology events, production
A inflation and eligibility, complete traversal, or final partitions.
