---
title: Local Marchenko-Pastur Rule
type: analysis
status: reviewed
updated: 2026-11-04
sources:
  - wiki/sources/selected-geometry-mp-integral-literature-20260602.md
  - wiki/sources/local-mp-identity-law-diagnostic-20260602.md
  - wiki/sources/legacy-internal-spectral-comparison-panel-20260616.md
  - tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/decomposition.py
  - tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/operators.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projection_dimension_estimation/projection_dimension_estimators.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/spectral/spectral_decomposition_result.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/spectral/marchenko_pastur.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/spectral/tree_estimator.py
  - manuscript/sections/method/edge_test.tex
  - manuscript/sections/method/assumptions_validation.tex
  - reports/statistical_construction_audit_20261104/mp_audit.md
  - reports/statistical_construction_audit_20261104/mp_results.json
  - benchmarks/validation/contracts/method_constants_manifest.py
  - raw/assets/mp-dimension-rule-analysis/20260525-controlled-spectrum-summary.md
  - raw/assets/mp-dimension-rule-analysis/20260525-dimension-contract-subset.csv
  - raw/assets/mp-dimension-rule-analysis/20260525-dimension-contract-finite-null-smoke.csv
  - raw/assets/mp-dimension-rule-analysis/20260526-threshold-policy-targeted-smoke.csv
  - raw/assets/mp-dimension-rule-analysis/20260601-leaf-only-mp-contract-targeted.csv
  - raw/assets/mp-dimension-rule-analysis/20260601-regression-gate-leaf-only-floor-contract.csv
  - raw/assets/mp-dimension-rule-analysis/20260601-regression-gate-leaf-only-spectral.csv
  - raw/assets/mp-dimension-rule-analysis/20260601-full-kl-only-leaf-only-spectral.csv
  - raw/assets/mp-dimension-rule-analysis/20260601-full-kl-only-leaf-only-spectral-failure-report.md
tags:
  - method
  - spectral
  - validation
---

# Local Marchenko-Pastur Rule

## Summary

The primal and dual backend spectra have the same scale. Both centered
scatter matrices divide by the actual spectral row count \(N_u\), and they
share the same positive eigenvalues. The implemented identity-population
Marchenko--Pastur edge uses a separate leaf-row count \(m_u^{\mathrm{MP}}\):
\[
b_u^{\mathrm{MP}} = \left(1+\sqrt{d_u/m_u^{\mathrm{MP}}}\right)^2.
\]
The dual branch does not transpose this aspect ratio. Optional internal rows
increase \(N_u\) without increasing the current worker's MP leaf count.

That does not make the rule a validated calibration-and-power guarantee. The
formula is an asymptotic bulk edge for its specified independent, unit-scale
reference ensemble, not a finite-sample quantile. The edge spectral context
uses the assigned basis halves; the sibling context uses all descendant leaves
and supports opt-in internal rows. `effective_independent_rows` records a leaf
count, not an estimated number of independent observations. Raw MP count,
nominal rank, available positive rank and final retained rank are different.
When contrast-energy selection is enabled and a positive basis is available,
the kernel replaces the nominal request using the tested contrast. The
reported chi-square df is the final retained row count, whose null-law
justification remains separate. See [[statistical-construction-audit]].

A deformed MP law can describe certain fixed nonidentity population spectra
\(H_u\) under its own assumptions. It does not automatically describe
whitening fitted from the same rows, internal averaging operators or selected
tree membership. Deriving the correct ensemble precedes choosing a spectral
replacement; changing the nominal edge alone does not calibrate the Wald test.

The 2026-06-02 local identity-law screen shows that this question is active,
not cosmetic. Bernoulli and discretized Gaussian selected node spectra often
sit near or above the identity MP edge. Continuous empirical-covariance spectra
sit far below the identity MP positive support under the current active-feature
aspect-ratio calculation because they follow a different self-whitening law:
after each node's rows are whitened by the same node's empirical covariance,
positive covariance eigenvalues concentrate at the backend scale
\((m_u-1)/m_u\). Therefore \(H=\delta_1\) is not a uniform empirical
description across production feature families and selected contexts.

The 2026-06-01 regression-gate check shows how the stricter leaf-only contract
reacts in practice. With \(k_{\min}=2\), the regression gate completed with
mean ARI \(0.6017\), median ARI \(0.6803\), exact \(K\) in 4 of 17 rows in the
shared benchmark output, and six explicit unsupported-calibration skips. A
floor diagnostic on the same 17 cases found that \(k_{\min}=1\) reduced errors
from six to two but lowered mean ARI from \(0.6017\) to \(0.5182\) and median
ARI from \(0.6803\) to \(0.5004\). Therefore changing the floor is not a
validated fix; the exposed problem is missing empirical-null calibration
support when the stricter leaf-only spectral basis opens all local sibling
contexts as selected non-null.

The 2026-06-01 full TBS-only benchmark with plots and relationship analysis
disabled completed 110 cases under the same contract. It produced 85 `ok` rows
and 25 explicit `skip` rows, with valid-row mean ARI \(0.8612\), valid-row
median ARI \(1.0\), and exact \(K\) in 63 of 110 rows. All but one skip were
strict sibling-calibration-support failures; the remaining skip was the
intentional dense continuous covariance limit for the 20,000-feature continuous
stress case.

The 2026-06-16 copied-commit comparison records a then-current reconstruction of
the old internal-barycenter row path. On a one-replicate overlap panel, the
legacy diagnostic increases raw MP signal counts by about `191` nodes on
average and adds thousands of MP threshold rows, but the completed partitions
match the current leaf-only path exactly. This confirms that internal rows
alter the local spectral algebra but are not, by themselves, a validated
anti-fragmentation rule.

## Details

### What The Code Calculates

For each internal node, the spectral worker maps descendant rows into the
null-whitened tangent coordinates used by the projected-Wald tests. It then
removes constant active features, centers the rows, and eigendecomposes either
the feature covariance or the smaller dual Gram matrix. Both paths divide by
the physical spectral row count \(N_u\), so the positive dual eigenvalues remain on the same
scale as the feature covariance.

The dimension estimator counts eigenvalues above
\((1+\sqrt{d_u/m_u^{\mathrm{MP}}})^2\), then floors the returned dimension to
the configured minimum. Here \(m_u^{\mathrm{MP}}\) is the row count used by the
threshold; it is recorded separately from the field
\(m_u^{\mathrm{eff}}\), which currently also counts supplied leaf rows. In the edge-gate context, the default minimum dimension is
two. The raw MP signal count and the returned test dimension are therefore
different mathematical objects:
\[
k_u^{\mathrm{raw}}
=
\#\{\lambda_{u,i}>b_u^{\mathrm{MP}}\},
\qquad
k_u^{\mathrm{test}}
=\operatorname{extendGroups}\!\left(\min\{r_+,\max\{2,k_u^{\mathrm{raw}}\}\}\right).
\]

Here \(r_+\) is the available positive rank and zero remains zero. Numerical
eigenvalue-group completion can expand the nominal rank beyond the floor.
With energy selection, the final test rank can be below or above this nominal
value; the available candidate basis still limits it.

For a leaf-only centered Gaussian matrix with fixed known whitening, the backend divides by the original row count
rather than by the centered rank. A scale-matched centered Gaussian edge would
be \(\frac{m_u-1}{m_u}(1+\sqrt{d_u/(m_u-1)})^2\), which is slightly lower than
the current asymptotic edge and increases finite-sample false spike rates in
controlled null probes. The more conservative unscaled \(m_u-1\) edge is not
scale-matched to the backend covariance. This reinforces that local
finite-sample calibration needs a null quantile or asymptotic fluctuation
correction, not only an \(m\) versus \(m-1\) substitution.

The code exposes these objects through a typed spectral decomposition result:
`raw_mp_signal_counts_by_node`, `test_projection_dimensions_by_node`,
`effective_independent_rows_by_node`, and `mp_threshold_rows_by_node`.
These dimensions are nominal inputs. The projected-Wald kernel may replace a
positive request with contrast-energy rank and then extend the cutoff through
ties. It uses the actual retained row count as the nominal reference df.
The 2026-11-04 replay shows zero MP spikes with nominal rank 30 after full
same-row covariance whitening creates one numerical tie group; separate
contrast-energy controls return rank 1 or 30 from nominal requests 1, 2 or 30.

### What The Controlled Probe Shows

Under pure independent Gaussian noise, the raw asymptotic upper edge produced
at least one false spike in roughly 8--13% of the small local matrices checked.
Using \(m_u-1\) after centering lowered some false-spike rates, but it did not
turn the asymptotic edge into an exact finite-sample test. This is expected:
the largest noise eigenvalue has finite-sample fluctuation, and the MP edge is
an asymptotic support boundary rather than a calibrated finite-\(m_u\)
quantile.

The rank-one spiked probe behaved qualitatively as expected: weak spikes below
the high-dimensional edge were often missed, while stronger spikes were counted
reliably. This supports the rule as a sensible signal-dimension heuristic, not
as a proof of optimal power for the tree-selected projected-Wald pipeline.

The most important failure mode is internal spectral rows. When deterministic
internal distribution summaries were stacked with leaf rows and counted as
additional observations, high-dimensional null spectra could cross the MP edge
far too often. In a controlled \(n_{\mathrm{leaf}}=40,d=100\) case with about
62 augmented rows, the augmented-row threshold produced at least one null spike
in every run and an average of about 2.4 counted spikes. Counting only the leaf
row count for the threshold was much more conservative in that probe.

The other important failure mode is null-whitening scale error. For
\(n=40,d=100\), scaling pure noise by \(1.2\) made the fixed unit-scale rule
count null spikes in every run, with a mean of about five counted components. A
bulk-median scale estimator corrected that isolated scale-misspecification
probe, but it did not solve the dependent/internal-row setting. Therefore a
median-scale rule is a candidate diagnostic or validation arm, not yet a
production replacement.

### Dimension-Contract Diagnostics

Historical diagnostics compared internal-row PCA directions, leaf-only spectra,
projection floors \(0,1,2\), and finite-null upper-edge smoke variants. Those
runs are still useful for failure attribution: the old augmented-row behavior
could improve selected benchmark ARI, especially in high-cardinality
categorical settings, but it did so through a spectral basis that is not
calibrated by the fixed-subspace projected-Wald reference. The historical
diagnostic script compared leaf-only MP floor and finite-null variants.
Current edge basis construction is leaf-only; the sibling context still
supports an explicit internal-row configuration. The archived comparison does
not establish validity for that option.

The historical floor-zero experiment encountered the kernel's rejection of a
zero-dimensional request with nonzero contrast. The current edge wrapper
instead emits \(Q=0,df=0,p=1\) before reaching the kernel when its basis rank
is zero. This sentinel means no edge test; it does not establish equality or
validate a full-method floor-zero benchmark. Floor \(1\) was runnable in that
historical experiment, but it over-split at least one
clear Gaussian case that floor \(2\) recovered exactly in the representative
subset.

The finite-null smoke with 50 null replications on four cases matched the
floor-two clustering outcomes of the current leaf-count threshold, but it is
not a validation study. A real finite-null calibration would need enough
replications and case coverage to estimate local upper-edge quantiles
reliably.

A targeted 2026-05-26 threshold-policy smoke showed why this is a real method
decision rather than a mechanical cleanup. On `binary_many_features`, the old
augmented behavior, leaf-count behavior, and finite-null floor-two variants all
recovered four clusters with ARI \(1.0\). On `cat_highcard_20cat_4c`, the old
augmented behavior found two clusters with ARI about \(0.36\), while the
leaf-count and finite-null variants collapsed to one cluster. After the
selected-PCA validation, that benchmark advantage is treated as evidence of an
unresolved calibration-support problem, not as justification for retaining the
internal-row spectral basis.

A 2026-06-01 leaf-only diagnostic sharpened the failure mode. For
`binary_many_features`, `cat_highcard_20cat_4c`, and
`overlap_heavy_4c_med_feat`, \(k_{\min}=1\) produced admissible supported
calibration records, while \(k_{\min}=2\) produced zero supported records and
only selected non-null positive-weight records. The same diagnostic showed
that floor 1 is not a general solution: on the regression gate it reduced
unsupported-calibration errors but weakened the aggregate clustering score.
The correct next mathematical object is therefore not an internal-row fallback
or an unvalidated floor change, but a calibrated rule for empirical-null
support or an external conditional-null calibration model.

The 2026-06-16 reconstruction appended internal barycenters and counted them
for the MP edge to match the historical comparison, while recording the leaf
count separately. The current worker counts only supplied leaf rows for the
MP threshold. In the historical compact overlap benchmark, augmentation produced much
larger raw MP signal and threshold-row sums, but no partition change among the
four completed paired rows. Thus internal rows are best treated as a spectral
sensitivity diagnostic and possible local covariate, not as a rescue condition
for split/pass-through decisions.

Selection-aware thresholding remains a different mathematical object from the
finite-null smoke. A finite-null edge samples fixed local null matrices. A
selection-aware edge would condition on hierarchy construction, edge-gate
selection, and the selected spectral context. The current evidence only
supports keeping selection-aware MP calibration as an open research question,
not adding it as a production fallback.

The literature-backed analytic alternative is to replace the white-noise edge
only after deriving or estimating a local spectral law \(H_u\). In the
Silverstein--Choi convention, a sample covariance spectrum with aspect ratio
\(c\) and population spectral distribution \(H_u\) has limiting Stieltjes
transform
\[
m_u(z)=-
\left(
z-c\int \frac{t}{1+t\,m_u(z)}\,dH_u(t)
\right)^{-1}.
\]
The local density is recovered from
\[
f_u(x)=\pi^{-1}\operatorname{Im}m_u(x+i0),
\]
and support edges are determined through
\[
z_u(m)=-\frac{1}{m}+c\int\frac{t}{1+t\,m}\,dH_u(t).
\]
For \(H_u=\delta_1\), this collapses to the production edge
\((1+\sqrt{d_u/m_u})^2\). Therefore the unresolved mathematical task is to
test whether \(H_u=\delta_1\) is a defensible local null-whitened model, and
if not, to derive a deformed-MP edge with a locked validation study.

The first local identity-law diagnostic makes the split concrete. In
`gauss_null_large`, `binary_many_features`, and `dim_diffuse_6c_136f`, raw MP
spikes occur in about `59%`, `75%`, and `75%` of evaluated nodes, with median
top eigenvalue over MP edge near one. In the continuous empirical-covariance
versions, raw spike frequency is `0%` and median top-edge ratios are about
`0.081` and `0.048`, while median top-eigenvalue and median positive-eigenvalue
ratios against \((m_u-1)/m_u\) are exactly `1.0`. That is not evidence for
changing production defaults by itself, but it rules out treating the identity
MP edge as a family-independent validated law.

The selected-tree spectral-law rerun separates the next diagnostic targets.
For Bernoulli and discretized Gaussian screens, exceeding the finite
identity-null top edge remains a selected-hierarchy inflation phenomenon:
`gauss_null_large`, `binary_many_features`, and `dim_diffuse_6c_136f` exceed
the 95% finite identity-null edge in about `54.5%`, `72.1%`, and `70.5%` of
evaluated nodes. Categorical screens should be handled at the selected
extreme-node level instead. `cat_highcard_20cat_4c` has only about `7.0%`
finite-null exceedances, and its log top-edge exceedance has strong
node-size/aspect-ratio relationships. `cat_highd_3cat_500feat` has about
`12.0%` exceedances but includes root/half-tree extremes with ratios far above
the finite identity-null edge plus many tiny selected extreme nodes. Therefore
the categorical open problem is not the same object as broad
Bernoulli/discretized selected spectral inflation; it is selected extreme-node
geometry coupled to categorical covariance, projection dimension, and sibling
testing.

The finite-sample identity-null comparison separates ordinary MP edge
fluctuation from selected-hierarchy spectral inflation. With 80 simulated
identity-null replicates per node dimension, the finite-null 95% top-edge
quantile is close to the asymptotic edge in these high-aspect-ratio settings.
`gauss_null_large`, `binary_many_features`, and `dim_diffuse_6c_136f` still
exceed that finite-null envelope in about `54.5%`, `72.1%`, and `70.5%` of
evaluated nodes. By contrast, `cat_highcard_20cat_4c` and
`cat_highd_3cat_500feat` exceed it in about `7.0%` and `12.0%` of evaluated
nodes. The Bernoulli/discretized Gaussian issue is therefore selected spectral
inflation beyond ordinary finite fluctuation; the categorical issue is more
localized and should be analyzed through covariance/projection/FDR interactions
and selected large-node extremes.

## Evidence

- `projection_dimension_estimators.py` implements the fixed unit-scale MP
  upper edge and returns an explicit dimension estimate with raw count, test
  dimension, effective independent row count, and MP threshold row count.
- `spectral_decomposition_result.py` defines the typed tree-level spectral
  output contract consumed by edge gate and sibling gate.
- `decomposition.py` and `operators.py` show that both primal and dual
  eigendecomposition paths divide by the same row count, preserving the
  feature-covariance eigenvalue scale.
- `marchenko_pastur.py` maps descendant leaf rows through the null-whitened
  tangent chart and records the row count used for the current MP threshold.
- `tree_estimator.py` documents the current leaf-only spectral orchestration.
- `edge_test.tex` states the MP edge and the minimum-dimension floor in the
  manuscript method.
- `assumptions_validation.tex` marks the MP rule and minimum spectral dimension
  as validation gaps, and records leaf-only spectral rows as the active
  inferential contract.
- `method_constants_manifest.py` lists the MP upper-edge threshold as an
  explicit method constant requiring dimension-selection, null-calibration,
  and planted-signal validation outputs.
- `20260525-controlled-spectrum-summary.md` records the controlled null,
  spiked, internal-row, and noise-scale probes from this audit.
- `compare_mp_dimension_contracts.py` runs leaf-only MP contract variants
  without adding production config.
- `20260525-dimension-contract-subset.csv` records the 14-case representative
  contract benchmark.
- `20260525-dimension-contract-finite-null-smoke.csv` records the finite-null
  smoke on four representative cases.
- `20260526-threshold-policy-targeted-smoke.csv` records the targeted
  threshold-policy smoke showing that finite-null thresholding is not a
  drop-in production replacement for the old augmented-row behavior.
- `20260601-leaf-only-mp-contract-targeted.csv`,
  `20260601-regression-gate-leaf-only-floor-contract.csv`, and
  `20260601-regression-gate-leaf-only-spectral.csv` record the benchmark
  reaction after removing internal spectral rows from the production path.
- `20260601-full-kl-only-leaf-only-spectral.csv` and
  `20260601-full-kl-only-leaf-only-spectral-failure-report.md` record the
  full-suite TBS-only benchmark reaction under the same production contract.
- `wiki/sources/selected-geometry-mp-integral-literature-20260602.md` records
  the Stieltjes-transform integral route for general MP spectra and the
  no-bootstrap selected-region geometry boundary.
- `wiki/sources/local-mp-identity-law-diagnostic-20260602.md` records the
  representative production-spectrum screen against the identity MP law.
- `wiki/sources/legacy-internal-spectral-comparison-panel-20260616.md`
  records the copied-commit internal-barycenter spectral diagnostic and its
  small overlap benchmark comparison against the current leaf-only path.
- `benchmarks/diagnostics/spectral/mp/compare_mp_dimension_contracts.py` was retired on 2026-10-16 (last present at commit 7cc4e1c0).

## Links

- [[open-mathematical-questions]]
- [[projected-wald-statistic]]
- [[oracle-gate-path-diagnostic]]
- [[legacy-internal-spectral-comparison-panel-20260616]]

## Open Questions

1. What finite-sample upper-edge calibration, if any, should replace the
   asymptotic edge if local Type-I behavior is the target? The targeted smoke
   does not support the tested finite-null edge as a drop-in replacement.
2. Does the local null-whitened tangent spectrum have \(H=\delta_1\), or is a
   deformed MP law with local population spectrum \(H_u\) required?
3. Is a local bulk-scale estimator valid for Bernoulli, categorical, and
   continuous null-whitened tangent spectra, or only as a diagnostic?
4. What explicit semantics would make \(k_{\min}=0\) a valid full-method
   experiment rather than a projected-Wald contract violation?
5. How much of the remaining benchmark weakness is caused by MP dimension
   selection versus hierarchy recoverability, sibling FDR, or empirical-null
   inflation support?
6. Can calibration support be defined from selected leaf-only contexts without
   reusing selected non-null records as empirical-null evidence?
