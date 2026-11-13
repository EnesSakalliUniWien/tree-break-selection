# Score attribution and shared-event Ward family reference

Development report, 2026-11-09, Europe/Vienna. See the locked
`protocol.md` and `variant_protocol.md`, the 300 dataset records,
four matched native-orbit illustrations, exact intervals, and source
hashes. Run from the repository root with the existing `.venv`, setting
BLAS thread counts to one. The runner requires a fresh output directory
at the same repository depth.

## Native statistic controls

The native raw record is the continuous sibling Wald statistic from the
Ward/NNLS fit, guarded continuous covariance, native edge projection
context, and .90 contrast-energy rank. For a fixed fitted context the
variance multiplier is

\[
g = 1 + (L_L+L_R)/(2\bar L_+),
\]

where the denominator averages strictly positive lengths across the
whole fitted tree. Missing or zero sibling time uses one. Multiplying
all branch lengths by the same positive constant preserves the positive
edge set and both numerator and denominator scale, leaving g unchanged.
This is an invariance of the branch factor with the fitted topology and
other score inputs held fixed. It is not a statement about refitting to
changed coordinates, selection, or an underlying branch noise law. Native
NNLS estimates branch lengths from pair distances across the full
tree. Its default target uses feature standard deviations fit to the
full matrix, so changing rows in a subtree can change the global
normalizer and edges outside it. The default child-parent edge branch
variance policy is none, and the sibling factor is still applied.

Holding the unregularized base covariance and projection fixed,
multiplying contrast variance by g multiplies the whitened contrast
by g^(-1/2) and the projected quadratic by 1/g. The production covariance
adds a fixed ridge after scaling, and may add jitter to a full covariance
block. The exact unit-branch control reruns the native score with g=1.
It can differ from g times the native score and can in principle change
the contrast-energy rank. The initial development assertion incorrectly
required the ridge-free scaling at relative tolerance 1e-12 and stopped
on dataset four. The first three records and original protocol and
runner are preserved. The corrected comparison kept the same draws,
cells, estimands, and budgets.

All 100 observed native fits completed. Median g was 1.656 in the null
cell and 1.544 in the signal cell. The largest absolute difference from
the ridge-free identity was 4.20e-9, with no rank changes under the
unit-branch control in these cells. The oracle replaces the covariance
with the externally known simulation covariance while retaining the
native projection context and rank rule. The Gram control replaces the
projection and rank by the leading raw subtree Gram direction at fixed
rank one. Neither control supplies an analytic p-value for the guarded
native statistic.

The displayed statistic uses raw T divided by reference scale and applied
empirical inflation, where supported. The observed display factor was
available for 4/50 null and 48/50 signal fits. Its availability depends
on the native pipeline's eligibility, so it cannot serve as a universal
replay score without defining a common selection event or a score on
every required state. The matched illustration uses one accepted
trajectory per input/budget for all variants and gives the display rank
only if every required state has its native factor.

| Input | Steps | Native raw p | Unit branch p | Oracle covariance p | Gram rank 1 p | Display p | Seconds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Null | 99 | .94 | .86 | .94 | .98 | | 1.01 |
| Null | 999 | .994 | .101 | .982 | .984 | | 7.48 |
| Shift 6 | 99 | .03 | .07 | .18 | .03 | .10 | 1.92 |
| Shift 6 | 999 | .006 | .003 | .23 | .003 | | 13.90 |

These two inputs illustrate score sensitivity. They neither select a
variant nor estimate comparative power.

## Conditional family argument

Let E=e be the common event that a prespecified path selects the observed
labeled subtree with its entire Ward topology and all ancestral child
pairs. This fixes the family of every internal node within the subtree.
For each true node j, let C_j contain its column mean, centered feature
Gram, and outside data. For fixed row labels under its null, the tested
rows are independent Gaussian with one common mean and covariance,
independent of the outside rows. The conditional distribution given
C_j and E=e is the restricted left orthogonal orbit. The reversible
random-position rank is superuniform under this distribution in exact
arithmetic. Its score choices and kernel must be prespecified or fixed
by the same conditioning. Integrate the node-specific nuisance variables
over their distribution given E=e. Each true-node p remains
superuniform conditional on this same event e. Holm then gives strong
familywise homogeneity control conditional on e without requiring joint
exchangeability of the nodes' ranks or independence. Replacing
unresolved ranks by one for adjustment keeps the family size fixed and
preserves the bound. Unresolved node decisions remain visible.

The null requires common row means within the node. The heterogeneous
means (-a,a,-a,a) have equal first-half and second-half means. For a
nonzero a their Helmert mean has squared norm 4a². Left orthogonal
rotation changes its mean direction, invalidating the orbit argument
for the weaker equal-aggregate-contrast null. The family error bound
concerns homogeneity tests, with no implication for partition recovery,
BH over traversal, or error rates under covariance heterogeneity.

The first Holm threshold is alpha/m. A chain needs at least
ceil(m/alpha)-1 steps to reach it. At seven hypotheses and alpha=.05,
the minimum is 139 steps, and at 127 hypotheses and alpha=.01 it is
12,699 steps. This is necessary resolution, with additional power
constraints from the restricted event, score, and movement.

## Development results

The exact fixed-Gram Beta helper and orbit score agreed on the observed
statistics for all 50 fixed-contrast datasets. The exact Beta helper
rejected 4/50, the 99-step chain 2/50, and the 999-step chain 3/50.
Median absolute difference between orbit and exact p fell from .02625
to .00977. The finite randomized ranks need not equal the exact
conditional tail.

| Family | Steps | False family rejections | True family rejections | Unresolved | Median seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Four-row homogeneous root | 99 | 1/50 | 0/50 | 0/50 | .011 |
| Four-row homogeneous root | 999 | 0/50 | 0/50 | 0/50 | .105 |
| Selected homogeneous child | 99 | 1/50 | 0/50 | 0/50 | .012 |
| Selected homogeneous child | 999 | 3/50 | 0/50 | 0/50 | .112 |
| Eight-row mixed root | 99 | 0/50 | 0/50 | 50/50 | .026 |
| Eight-row mixed root | 999 | 1/50 | 0/50 | 0/50 | .251 |

The 99-step mixed-root panel cannot resolve the first Holm threshold,
so the zero count there gives no calibration evidence. The 999-step
mixed panel demonstrates neither useful power nor practical scale-up.
The exact intervals are too wide for promotion, and the native replay
timings leave the runtime gate unestablished. No candidate from these
development panels is ready for promotion.

Wolfram independently returned zero for the positive global branch
scaling difference, zero for the ridge-free quadratic scaling
difference, and zero for the cancelling contrast, with Helmert mean
(-sqrt(2)a, sqrt(2/3)a, -2a/sqrt(3)) and squared norm 4a².
It also returned both minimum-step counts above. The scalar ridge
difference simplifies to
lambda(1-g)/((v+lambda)(g v+lambda)) for
1/(v+lambda)-g/(g v+lambda).

Fifty-six focused tests passed. A final failure-path change makes
observed Ward numerical failures unresolved. The exact family source
used in the development run is preserved alongside the manifest, and
the first dataset in each of the three family cells reproduced after
the change, excluding runtime. The stored rank comparison is numerical
and does not certify exact-real Ward membership or ties.
