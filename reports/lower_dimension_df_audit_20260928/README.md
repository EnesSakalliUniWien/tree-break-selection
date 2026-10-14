# Lower dimensions and degrees of freedom

Date: 2026-09-28. Two Luna agents review the current implementation and prior
research. No production statistic, degrees of freedom or default is changed.

## Finding

The concern is real, but there are two distinct questions:

1. Does the recorded degrees of freedom count the coordinates actually used
   in the statistic? In the inspected kernel, yes: it uses the length of the
   projected vector, including complete numerical eigenvalue groups.
2. Does that statistic actually have the assumed chi-square distribution?
   Not generally when the projection or its dimension is chosen using the
   tested contrast, or when its null covariance is not the identity.

Thus a correctly recorded lower dimension can still produce an incorrectly
calibrated test. Replacing it with the original feature count is not a general
repair. Passthrough acts later and does not resolve this issue.

## What the degrees of freedom means

For a null vector `z ~ N(0,I_d)` and fixed orthonormal projection rows R of
rank k, `T = ||R z||²` has law `chi_square(k)`. The rank counts retained
independent directions, not distinct eigenvalue values or the original feature
count. An eigenvalue repeated m times can contribute m retained directions.
For the symmetric covariance operators used here, algebraic and geometric
multiplicities agree in exact arithmetic; numerical ties are tolerance-defined
groups rather than a proof of exact multiplicity.

If R is learned from independent data and the evaluation vector is still
isotropic Gaussian conditional on those data, the same statement holds
conditional on R and k. Learning R from the tested data does not preserve that
condition merely because its rows are orthonormal.

An exact counterexample makes the distinction explicit: choose the one-row
projection `R(z) = z.T / ||z||`. Its rank is one, yet
`||R(z)z||² = ||z||² ~ chi_square(d)`. The nominal rank-one reference is wrong.
This illustrates selection; it is not the repository's PCA algorithm, and it
does not prove that every adaptive test should use df=d.

For a fixed R but `Cov(z)=C`, the eigenvalues w of `R C R.T` give
`T ~ sum_i w_i chi_square_1` under a Gaussian null. These are **null quadratic
form weights**, not automatically the parent PCA eigenvalues. Matching two
moments gives scale `a=sum(w²)/sum(w)` and
`nu_eff=sum(w)²/sum(w²)`, an approximation that generally needs both scale
and degrees of freedom. It does not account for selecting R using z.

These distinctions follow the fixed quadratic-form result in
[Lockhart's notes](https://www.sfu.ca/~lockhart/richard/350/08_2/lectures/Theory/slides.pdf)
and the selected-test example in
[Tibshirani's presentation](https://www.stat.cmu.edu/~ryantibs/talks/ssc-2016.pdf).
The hierarchy adds a separate selection problem, as studied by
[Gao, Bien and Witten](https://arxiv.org/abs/2012.02936v3). These sources do not
establish a null law for this repository's full pipeline.

## New controlled check

`check_reference_law.py` uses 100,000 independent eight-dimensional standard
Gaussian contrasts, seed 2026092803, and alpha .01. All method comparisons use
the same contrasts. Covariance is known; there is no tree or fitted inflation.
The vectorized formulas were checked against 1,024 actual production kernel
calls. A separate tie case verifies nominal target one expands to actual df two.

| Projection/reference | Rejections / 100,000 |
| --- | ---: |
| Fixed coordinate, df 1 | 957 |
| Independent random direction, df 1 | 1,012 |
| Full fixed coordinates, df 8 | 985 |
| Production contrast-energy rank rule at .8, df equal to selected k | 2,148 |
| Illustrative contrast-aligned direction, nominal df 1 | 57,371 |
| Same illustrative statistic, its exact df 8 | 985 |

The energy rule selects enough ordered coordinates to retain .8 of the
observed contrast energy. Its reported df is exactly the number retained;
nevertheless, its rejection rate exceeds the .01 target in this controlled
example. This isolates contrast-dependent rank selection with a fixed basis.
It does not estimate the production pipeline's error rate or show how much of
its error comes from this optional rule.
The orchestrator and shared TBS runner default this energy fraction to `None`;
the rule is enabled only when configured. Same-data PCA/MP basis selection is
a separate concern even with the energy rule disabled.
The pancreas benchmark explicitly enables fraction .90 in
`applications/scrna/pancreas_benchmark.py::_base_tbs_params`; the audit did not
run that application or quantify its resulting calibration error.

A fixed two-dimensional correlated Gaussian example has null weights 1.8 and
.2. Using unscaled df 2 gives 2,424 rejections; the Satterthwaite approximation
(scale 1.64, df 1.219512) gives 1,048. Whitening with the known full covariance
gives 946 using df 2. This is a separate covariance control, not evidence for
installing an effective-df correction in production. Pointwise Wilson intervals
for every count are in `reference_law_summary.csv`.

## Previous research to retain

The June 1 locked Gaussian fixed-membership study and June 5 rerun are still
available as primary CSVs. With mean projection dimension two in every setting,
the locked leaf-only cases reject 45–50/1000 at alpha .05, whereas appending
deterministic child-mean rows gives 648–1000/1000. This changes the spectral
input while leaving the df count unchanged. The near-nominal leaf-only result
is limited to that empirical-covariance design, rank floor and sample sizes;
it is not a theorem that same-data PCA is always valid, nor does later failure
in a different design imply that it must always fail.

The June 13 statistic-distribution-shape and null-law-decomposition panels
already distinguish df-mixture skew, covariance-derived scale/effective df,
and same-data projection selection. Their source audits are in
`prior_research_audit.md`; historical smoke results must not be treated as
newly rerun evidence.

The earlier isolated rank reassessment is retained at
`../eigenspace_reassessment_20260928/wald_rank/README.md`. At energy fraction
.90 and alpha .01, two held-out seeds gave rejection rates .014464/.014595
for distinct eigenvalues (2,1), and .015531/.015342 for the partially tied
spectrum (4,2,2,1), using one million contrasts per seed and dimension.
The fully tied (1,1) case became near nominal after complete-group retention.
That study reports excess in all 16 held-out adaptive distinct/partially tied
cells after Holm correction across its 216-cell family. The new .80 control
is a bounded confirmation of an already identified issue, not its discovery.
The earlier full-pipeline comparisons disabled the energy rule, so their
miscalibration cannot be attributed solely to this optional rank rule.

The September 28 independent-projection control strengthens the distinction:
holding rank at one, same-data PCA gave 33/1000 and 37/1000 global-null root
rejections, versus 11/1000 and 10/1000 for independent-data PCA. Because rank
was identical, changing the integer df alone cannot explain that comparison.
The canonical MP runs had almost no independent-rank variation, so they did
not establish that rank selection is harmless in other regimes. The new
contrast-energy example supplies a separate rank-selection counterexample.

## What remains unresolved

- A selected null law or independently validated calibration for the actual
  combination of tree inference, projection, rank, covariance and inflation.
- Whether covariance scaling/dependence or selection dominates in each actual
  feature family and node regime; no universal replacement df is established.
- Calibration of the optional contrast-energy rank rule after selection.
- Finite-sample covariance estimation and categorical boundary effects, which
  are absent from the known-Gaussian controls.

`implementation_audit.md` gives the caller and test audit. The controlled check
records before/after source hashes in `reference_law_verification.json`.
The narrow baseline passed 30 tests in projected-Wald statistics and eigenspace
invariance. These tests validate implemented behavior, not the adaptive null law.
The implementation agent separately passed 41 tests across projected-Wald,
contrast covariance and reference-law reporting; coverage overlaps with the
root baseline. Report scripts pass Ruff. Final wiki and report-link checks
are recorded in the wiki log.
