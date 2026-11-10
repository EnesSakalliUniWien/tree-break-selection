# Independent audit of the three sibling-test constructions

Date: 2026-11-09, Europe/Vienna. Scope: the supplied A/B/C description, current
implementation, independent reruns, and mathematical/statistical checks through
Wolfram. No production statistical code was changed.

## Finding

The numerical gains reported for B are reproducible, but the claimed general
error control is not established. The statement that every cut of a unimodal
population has rho <= 3 is false. B also has an incorrect boundary variance,
and its flat-cut count does not turn that reference into a valid p-value. C has
correct restricted-model constants and a useful change of scale, but its
uniform-marginal model, smooth-field argument, and Euler-characteristic envelope
do not supply a general selected-tree reference law. A's supplied description
omits current context weighting and contrast-dependent dimension selection.

## Reproducibility and scope

The runnable record is `results/analyses/three_tests_audit_20261109/`: `protocol.md`,
`manifest.json`, `rerun.py`, per-replicate outputs, and separate fit outputs.
Original proposal scripts and previous results were not overwritten. The source
hash manifest records the live production files used. The supplied description
is captured in `three-tests-supplied-construction-20261109.txt`; exact Wolfram
inputs and returned outputs are in `three-tests-wolfram-20261109.json`.
The executable source, protocol, execution notes, and original-derivation replay
are archived in `three-tests-audit-scripts-20261109.py.txt`.

Baseline: 84 focused production tests passed, both original proposal self-checks
passed, and wiki-lint passed on 259 pages. Environment: Python 3.12.11, NumPy
2.3.4, SciPy 1.16.3, pandas 2.3.3, scikit-learn 1.7.2, pytest 8.4.2. No dependency
installation, commit, push, or persistent service was needed.

Fresh seed 20261109: 15,000 low-dimensional and 2,400 high-dimensional C checks,
6,000 unimodal counterexample draws, 10,000 fixed-uniform-cut draws, and 100,000
standard-normal contrasts for A's rank-selection control. The C scripts supply
k directly: they do not validate estimating invariant directions from the same
data, choosing descendants, or a complete tree's final error rate. Searched
maxima at k > 1 are lower bounds, not certified global maxima.

The original Gaussian-bound simulation, exhaustive flat-cut check (140/140),
and Gaussian median-cut check (0.0095 rejection at nominal 0.01) also reproduce.
That last result concerns the Gaussian boundary and does not validate ceiling 3.

## A: current implementation and the chi-square reference

The branch-time multiplier in the supplied standardization agrees with
`sibling_divergence_test.py`. Three qualifications matter:

1. `resolve_node_continuous_covariance_by_block` defaults to guarded pooled
   within-child covariance, with parent-total covariance as its fallback.
2. The 90% rule in `run_projected_wald_kernel` retains the leading PCA directions
   until their squared **observed standardized contrast** reaches 90% of its
   total energy; it is not a 90% PCA-eigenvalue variance rule. Numerical
   eigenvalue ties are retained together.
3. `decide_empirical_null_calibration` normally estimates a local factor using
   kernel weights in log projection dimension and log parent size, within the
   feature family. It also has a family-baseline fallback. Thus A does not
   generally use one constant factor for the entire continuous-feature tree.

For a known correctly standardized Gaussian contrast and a fixed orthonormal
basis P of fixed rank d, ||P' z||^2 is chi-square_d. Fixing group membership alone
does not also fix the covariance estimate, P, d, or branch-time multiplier.
For fixed Gaussian groups with an estimated pooled covariance, the exact
unprojected two-sample reference is Hotelling's scaled F:

    (n-k-1)/(k(n-2)) * n rho ~ F[k, n-k-1]

under the common-covariance null. The noncentral F version requires the same
Gaussian-group model under an alternative. Threshold-selected groups are
truncated samples and do not satisfy that model merely because their parent
population was Gaussian.

A direct control with independent z ~ N(0,I_30), a fixed identity basis, and
distinct supplied eigenvalues gave 1,604/100,000 rejections at nominal 0.01 after
the 90% contrast-energy rule. The first 20 calculations matched the production
kernel exactly. This isolates rank selection; it is not a measurement of the
complete production pipeline's rejection rate. Estimating inflation from the
same selected tree does not itself prove that this or the other selection
effects have been calibrated.

Current BH is applied over eligible top-down frontiers in
`apply_traversal_aligned_sibling_bh_results`, rather than one global family of
all sibling tests. B/C's diagnostic replacement occurs after empirical-null
eligibility is computed, so that inherited eligibility remains part of their
implemented procedure.

## B: which algebra survives, and where the reference fails

### The statistic and the best-flat-cut argument

For positive-definite S_w and both sides nonempty,

    rho = pi(1-pi) d' S_w^(-1) d

is the largest between/within ratio over linear directions, using the stated
sample divisors. For any fixed side size, take a direction attaining that
ratio. Assigning the largest projected observations to one side maximizes its
mean contrast. Total variance is unchanged, so the between/within ratio cannot
decrease. Consequently a maximum over all partitions is attained by a flat cut
(with the usual treatment of ties). This deterministic comparison does not
identify the sampling distribution of the maximum.

### Ceiling 3 requires a restriction on the cut

For a scalar unimodal variable cut at its mode, Khintchine's representation
implies conditional variances at least one third of squared conditional means
relative to that mode. If their magnitudes are a,b and the left probability is
pi, then

    rho <= 3 pi(1-pi)(a+b)^2 / [pi a^2+(1-pi)b^2] <= 3,

because the denominator minus pi(1-pi)(a+b)^2 equals
[pi a-(1-pi)b]^2. Wolfram confirmed the last identity. The representation is
stated and proved in [Ion, Klaassen and van den Heuvel, Appendix A.1](https://arxiv.org/html/2208.08813v2).
This is a statement about a mode cut; it does not cover arbitrary cuts chosen
by a tree or arbitrary multivariate projections without further assumptions.

A concrete continuous unimodal counterexample is

    X ~ 0.9 Uniform(0,0.001) + 0.1 Uniform(0,1), cut at 1/3.

Its density is nonincreasing on its support. Wolfram gives the population ratio
1367002729/184734185 = 7.3998363053. Shrinking the short component and the long
tail's mass yields a ratio approaching 8 at that cut. The original 8 November
derivation already documented this limitation; the supplied summary omitted it.

At n=400, B rejects 1,996/2,000 such samples at nominal 0.01; C rejects
1,997/2,000. At n=1,600, each rejects 2,000/2,000. These are tests of the proposed
tails with k=1 and the specified fixed threshold, not complete tree fits.

### The noncentral chi-square has the wrong boundary variance

For a centered uniform population, the influence function of the R statistic
can be written, with Y uniform on [-1,1], as

    h(Y) = sqrt(3) (|Y| - 1/4 - 3Y^2/4).

Wolfram gives E[h]=0 and E[h^2]=1/40. Since

    d[r^2/(1-r^2)]/dr at r=sqrt(3)/2 = 16 sqrt(3),

the actual leading variance of sqrt(n)(rho_hat-3) is 96/5 = 19.2. B's
chi-square_k(3n) instead has variance 2k+12n, giving a leading value of 12.
This mismatch persists as n grows.

The fresh fixed-threshold uniform check at n=1,600 measured n Var(rho_hat)
=18.9018. The unmultiplied noncentral tail rejects 337/10,000 at nominal 0.01
(95% exact binomial interval 0.03025-0.03743), agreeing with the asymptotic
prediction 0.032947. C's fixed-direction tail gives 82/10,000 in this control.

The complete B count also fails on a boundary case from the existing validation
script: uniform on the radius-sqrt(3) sphere in three dimensions, n=100, where
every one-dimensional projection is uniform and the population best-cut ratio
is 3. B rejects 40/1,000 at nominal 0.01 (95% interval 0.02873-0.05407).
The maximum is searched, so this observed rejection rate is a lower bound for
using the exact maximum. This directly contradicts a general conservative-tail
claim, even for this specific ceiling-3 boundary model.

### Cover's count is a combinatorial fact

[Cover's original Theorem 1](https://isl.stanford.edu/~cover/papers/paper2.pdf)
gives 2 sum[j=0..d-1] binomial(n-1,j) homogeneous labeled dichotomies in general
position. Lifting affine k-dimensional cuts to k+1 dimensions and identifying
complementary labels gives the proposed sum[j=0..k], including one empty-side
partition. Removing that partition subtracts one. Retaining it is harmless as
an upper count.

The count does not give the tail probability of a cut selected from the sample.
A union bound needs valid tails for a suitable indexed family, or a conditional
bound that remains valid after its geometry is selected. Counting the distinct
partitions realized on each sample does not, by itself, supply that argument.
It also cannot repair the false population ceiling or incorrect per-cut law.

## C: correct identities, restricted constants, missing guarantees

### Exact change of scale

Let S_T be total covariance with divisor n. The identity

    S_T = (n-2)/n S_w + pi(1-pi) dd'

and Sherman-Morrison give

    R_hat^2 = pi(1-pi) d' S_T^(-1) d = rho'/(1+rho'),
    rho' = n rho/(n-2).

This was checked numerically independently of the proposal. The transformation
is correct. Replacing the whole best-cut problem by a centered absolute-moment
statistic is only a first-order argument for an appropriate symmetric null
whose population optimum is centered and balanced. In general the side-share
factor and threshold cannot be omitted.

### The constants are model-specific

The preceding influence function verifies variance 1/40. For X uniform on the
three-dimensional sphere of radius sqrt(3), conditioning on its first
coordinate gives tangent-coordinate variance (3-X_1^2)/2. Integrating the
squared derivative of the influence function yields 13/40, hence lambda=13.
Wolfram checked these integrals. Gaussian influence variance is 1-3/pi, also
confirmed by Wolfram.

There cannot be a k-dimensional population with a Uniform(-sqrt(3),sqrt(3))
marginal in **every** unit direction when k>3. Those marginals imply

    E[X_i^4]=9/5, E[X_i^2 X_j^2]=3/5,
    Var(||X||^2)=2k(3-k)/5.

The last quantity would be negative. Thus carrying the uniform-marginal field
model into k=6 or 9 is not just an unestimated smoothness parameter. It requires
a different population model and an argument for its set of maximizing cuts.
Products of uniforms have maximizing coordinate directions, not the complete
sphere of equivalent optima.

### The cited smooth-field theorem does not apply directly

For the actual sphere-3 influence field, direct moment calculations give

    Cov(h(v'X),h(w'X))
      = (2/pi)(sqrt(1-r^2)+r asin(r)) - (51+27r^2)/80,
    r=v'w.

Wolfram expands its normalized angular correlation as

    1 - (13/2) theta^2 + [80/(3pi)] |theta|^3 + O(theta^4).

The absolute cubic term obstructs twice mean-square differentiability. An
independent two-million-point sphere simulation checked the covariance formula
at four angles within five Monte Carlo standard errors. The Gaussian variant
also has a nonzero absolute cubic term.

[Taylor, Takemura and Adler, Assumption 2.1](https://arxiv.org/pdf/math/0507442)
requires a square-integrable C2 process, along with further regularity. Their
result establishes a high-threshold approximation with a controlled error
under those assumptions. It is not a blanket finite-sample upper bound for
these fields. Working on real projective space handles the antipodal equality;
halving the sphere curvatures is consistent with that quotient, but it does not
resolve the smoothness problem.

The L_j/EC-density implementation passes the original unhalved chi-distribution
identity checks. Taking a nonincreasing envelope removes numerical sign swings;
it does not prove domination of the true survival probability.

### Fresh calibration results

Across the 15 original low-dimensional cases, C's 0.01 rejection counts range
from 0 to 16 out of 1,000 per case. The Gaussian cases use C's separate Gaussian
constants, not the ceiling-3 constants. The ceiling-3 cases cover k=1 and k=3;
k=2 was tested only with Gaussian constants. At nominal 0.05, the one-dimensional
uniform case with n=100 rejects 75/1,000 (95% interval 0.05945-0.09311), so a
blanket finite-sample calibration claim is false. Uniform sphere-3, n=100 gives
1/1,000 at 0.01 and 14/1,000 at 0.05: a conservative observation on this panel.

At k=6 and 9, n=400 and 1,600, the fresh searched-maxima check gives 0-1/300
rejections per case at 0.01. Even 0/300 has a two-sided 95% exact upper bound of
0.01222, and the search may miss larger cuts. This is evidence about that panel,
not proof of calibration or conservatism.

## Rerun results on saved splits and the continuous benchmark

The saved 1,838 Markov/Ward split calculations reproduce exactly. At raw p<=0.01,
B/C retain 14/1 splits at 32-159 cells, 65/62 at 160-999 cells, and 42/43 at
1,000+ cells; neither retains any of the 363 splits with <=31 cells. None of the
27 hidden-boundary top splits passes C. This replay alone performs no refit,
BH, or traversal.

All 54 original continuous-benchmark rows return status `ok` and match the old
cluster counts and ARIs exactly. Mean ARI for B is 0.695011 for both trees;
A gives 0.880198 on default and 0.852483 on Markov/Ward. The rule-out comparator
also matches. The n<=p+1 early return in `invariant_directions` forces k=0,
regardless of signal strength; this explains some of the wide/small-node losses.
The MP edge is an asymptotic noise-spectrum rule, not a finite-sample test that
all retained directions are signal. Selecting k and the subspace is an
additional part of the procedure that supplied-k simulations do not validate.

All six full-cell fits completed with status `ok`, using 48,879 cells, 30 saved
PCs, the original saved diffusion coordinates, sibling alpha 0.01, edge alpha
0.001, and the registered tree/branch settings:

| Tree | Test | Clusters | ARI against cell lines | Lines recovered at 90/90 |
| --- | --- | ---: | ---: | ---: |
| Default | A | 75 | 0.027502 | 27 |
| Default | B | 64 | 0.027521 | 29 |
| Default | C | 57 | 0.027748 | 28 |
| Markov/Ward | A | 954 | 0.548774 | 41 |
| Markov/Ward | B | 217 | 0.832948 | 119 |
| Markov/Ward | C | 173 | 0.815201 | 111 |

All four A/B fits reproduce the earlier cell assignments exactly, not just the
aggregate scores. C's end-to-end runs are new diagnostic applications: its
smaller Markov/Ward partition has lower ARI and fewer recovered lines than B.
The default tree remains near ARI 0.028 under all three constructions.

B and C use identical split IDs, dimensions and rho values (90 default-tree
records and 1,838 Markov/Ward records). All 1,928 C tail values were independently
recomputed and matched. All recorded production/input hashes remained unchanged,
and temporary distance files were removed. Full results and checks are archived
in `three-tests-rerun-evidence-20261109.json`; local assignments remain beside
the rerun script. These agreement metrics do not measure Type-I error on CHOIR.

## Interpretation and unresolved work

Keep the exact statistic identities and the deterministic best-flat-cut
comparison. Treat B and C as diagnostic scoring rules under stated population
models. A valid test needs an explicit null class, a justified tail for that
class with estimated covariance/subspace included, and treatment of selected
descendants. Improving clustering agreement on CHOIR does not establish those
properties.

BH cannot manufacture valid p-values from invalid input tails. Its usual
guarantees require suitable null p-values and dependence conditions; separately
testing selected frontiers also needs a corresponding hierarchical argument.
See [Benjamini and Yekutieli (2001)](https://www.math.tau.ac.il/~ybenja/depApr27.pdf).
No final-partition FDR/FWER theorem was established for A, B, or C here.

The earlier 935-dataset end-to-end panel was not rerun in this audit. Diffusion
coordinates and PCs were reused, as in the original full-cell experiment;
their preprocessing was not regenerated. No production replacement was made.
