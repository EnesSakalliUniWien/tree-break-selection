---
title: Analytic Selected Root Laws
type: analysis
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/analytic-selected-root-derivation-20261102.md
  - raw/inbox/sibling-null-mechanism-isolation-20261102.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/gaussian_root_region.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/gaussian_root_unknown_variance.py
  - tests/statistics/59_test_gaussian_root_region.py
  - tests/statistics/60_test_gaussian_root_unknown_variance.py
  - tree_break_selection/hierarchy_analysis/statistics/child_parent_divergence/child_parent_divergence_annotation/child_parent_divergence_annotation.py
  - tree_break_selection/tree/distributions.py
tags:
  - mathematics
  - calibration
  - selection
  - gaussian
---

# Analytic Selected Root Laws

## Summary

Two executable analytic diagnostics account for root selection under explicitly
restricted Gaussian models. They condition on the complete observed merge
sequence of average linkage on **squared Euclidean** distances. The first uses
a truncated chi law with known isotropic covariance; the second uses a
truncated F law for scalar observations with unknown common variance and a
global common-mean null. They implement the analytic direction chosen by the
user; they do not replace the general production edge/sibling gates.

## Details

### Known covariance

For observed root children L and R, let v_i=1/n_L on L and -1/n_R on R,
e=v/||v||, z=X' e/sigma, u=z/||z|| and A=(I-ee')X. For independent Gaussian
rows with covariance sigma² I and contrast null e'E[X]=0, condition on A and u.
The remaining radius r=||z|| has chi_d law before selection. The path

\[
X(r)=A+\sigma e r u^\top
\]

preserves those conditioned objects. Fix every cluster pair in the observed
merge sequence. Each requirement that its chosen pair has no larger average
squared distance than every available competitor is a quadratic inequality
in r. Intersect all such regions with r>0 to obtain S, possibly disconnected.

\[
p=\frac{\int_{S\cap[r_{\rm obs},\infty)}f_{\chi_d}(r)\,dr}
         {\int_S f_{\chi_d}(r)\,dr}.
\]

This conditions more strongly than merely retaining the final root partition.
`selected_gaussian_root_test(values, sigma=...)` recomputes the specified
hierarchy; it does not accept an arbitrary tree and infer its construction.

### Unknown scalar variance and small groups

For scalar iid N(mu,sigma²) observations under the **global common-mean null**,
let a=e'x, w=x-mean(x)1-ea, nu=n-2, eta=w/||w||, s=sign(a) and
R=||x-mean(x)1||. Condition on mean, R, eta and s. Then

\[
h=|a|/\|w\|,\qquad \nu h^2\sim F_{1,\nu},\qquad
x(h)=\bar x1+\frac{R(\eta+seh)}{\sqrt{1+h^2}}.
\]

The common factor R²/(1+h²) cancels from all squared-distance linkage
comparisons. Their quadratic interval region therefore truncates the F law
without estimating a plug-in variance scale. The residual sum of squares is
within the observed root children; no between-child signal is inserted in the
variance denominator. `selected_scalar_gaussian_root_test(values)` requires
n>=3 and positive residual variation. Two observations have no residual df.

This is not a solution for arbitrary mean heterogeneity satisfying only e'mu=0:
the residual distribution can then be noncentral. Nor does the derivation cover
multivariate unknown covariance, adaptive PCA or arbitrary supplied distances.

### Numerical and inferential boundaries

- Quadratic coefficients use exact rational arithmetic on the supplied floating
  path factors; boundaries use Decimal arithmetic. This does not make the
  original floating data or all downstream probability integration exact.
- Log-tail differences and scaled quadrature avoid ordinary tail subtraction
  where possible. Unresolved masses and collapsed intervals raise errors.
- Observed minimum-linkage ties, zero contrast and unsupported residual states
  are explicit errors. Under the continuous nondegenerate Gaussian model these
  degeneracies have probability zero; rounded/duplicate data can encounter them.
- Full merge-sequence conditioning does not condition on the current adaptive
  edge-opening, empirical-null support, sibling BH or traversal events.
- The method's default precomputed distance can be ordinary Euclidean, Hamming,
  diffusion or another geometry. Squaring distances can change average-linkage
  trees. These functions therefore remain explicitly named diagnostics.

### Sibling mechanism isolation

The source capture calls the production projected-Wald kernel on 3,000 iid
40-by-8 Gaussian datasets, fixed 20/20 membership and known covariance. At
alpha .05, fixed rank-two coordinate projection rejected 132 times; same-data
PCA rejected 467; fixed basis with contrast-energy rank rejected 225. Applying
branch variance multiplier two to the PCA case reduced the count to 59.
Topology, MP selection and covariance estimation were held out of this control.
Thus a branch factor can conceal projection miscalibration rather than repair it.

A separate 200,000-draw independent scale control at focal df2, true scale3
and calibration df4 rejected 9,878 times with known scale and 26,477 with the
plug-in chi-square tail (including the production floor). An exact independent
F control rejected 9,896 times. The discrepancy decreased with larger calibration
df. Ordinary F cannot simply replace the production tail: its weighted,
dependent, selected calibration records do not satisfy this control's assumptions.

### Integration decision

The production edge docstrings were corrected: row splitting within an already
selected child does not establish independent membership or a calibrated law.
Production covariance, alpha, sibling projection, inflation and traversal were
preserved. Removing the small-node guard alone or using an ordinary F tail on
half-row contrasts with all-row covariance would be unjustified.

Before integration, specify the actual inferential null, distance/linkage,
projection/rank conditioning and whether the selected event includes edge-path
opening. The two analytic diagnostics provide a verified starting point for
that work; they are not advertised as general TBS calibration fixes.

## Evidence

- Full ordered `make test` verification passed 905 tests on 2026-11-02,
  including the new diagnostic tests. One scikit-learn MDS future-default warning
  appeared. Ruff and wiki structural checks passed.
- New focused tests cover closed-form 2/3-row cases, interval intersections,
  disconnected regions, path reconstruction and independent SciPy linkage
  replay, narrow/extreme tails, invalid inputs and scalar scale/location invariance.
- Independent 300-dataset n4 global-null checks rejected 20 times for known
  sigma and 13 for unknown sigma at alpha .05. This is finite verification,
  not proof of universal numerical calibration.
- Two distinct observations per child near -100 and +100 gave selected
  scalar p=3.1601693694910593e-6, residual df2. This resolves the variance
  construction for this declared model, not the general production path.
- Primary selective-inference references, Wolfram density check, exact commands
  and outputs are preserved in the two raw captures and summarized in
  [[analytic-root-and-sibling-evidence]].

## Links

- [[analytic-root-and-sibling-evidence]]
- [[mathematical-glossary]]
- [[mathematical-construction-wolfram-review]]
- [[root-selected-region-model]]
- [[open-mathematical-questions]]

## Open Questions

- Extend the event to actual edge opening and focal descendant selection.
- Treat adaptive projection and rank jointly with the conditional law.
- Derive a law for selected, weighted empirical-null scale estimation.
- Establish discrete-feature and unknown multivariate covariance contracts.
- Assess power loss from conditioning on the entire merge sequence and scalability
  beyond small diagnostic trees before production integration.
