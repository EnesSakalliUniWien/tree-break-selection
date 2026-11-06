# Conditional inference on the existing small benchmarks

## Scope and frozen inputs

This check uses the registered `gauss_clear_medium_continuous` (seed 101,
60 by 40), `gauss_moderate_3c_continuous` (seed 200, 45 by 30) and
`gauss_null_large_continuous` (seed 401, 200 by 80) cases at their original sizes.
The existing nine-row scalar reentry fixture is an interval-certification
control. Data come through `prepare_case_inputs` without altering their case
definitions. The diagnostic uses one native symmetric-diffusion average-linkage
hierarchy for each dataset. All radial proposals reconstruct the global graph
and full merge history. The proposal budget is 999 per tested contrast.

## Covariance and preprocessing

The three continuous blob generators have known within-population covariance
\(\sigma^2 I\), with sigma 0.6, 1.2 and 1.0 respectively. These declared scales
are used in the conditional statistic. The generated continuous data and the
matrix entering diffusion are exactly equal. Distance preparation creates a
separate array and leaves the observations intact. These cases have no fitted
upstream PCA or variable-gene selection.

The model audit also generated seven existing cases outside the prototype's
current observation model. `gauss_clear_small` and `gauss_null_small` are
median-binary observations. The dimensional, spike and low-rank cases have
nonspherical covariance. The outlier case has unequal noise laws. The Brownian
case has marginal dependence between rows from shared taxon effects. The audit
records these differences explicitly and assigns no Gaussian selective p-value
to those inputs.

The earlier count benchmarks, including `signal_balanced`, `signal_imbalanced`,
`signal_singleton` and `signal_density_variation`, use negative-binomial counts
followed by normalization, log transformation, variable-gene selection, scaling
and fitted PCA. Their preprocessing is outside the tested Gaussian path. Their
model audit is preserved separately from the continuous benchmark replays.

The production covariance resolver returned diagonal arrays for both signal
roots. Their average entries were 54.7261 and 9.7034 times the known Gaussian
noise variances. These roots contain several true populations, so their pooled
within-child variation also contains remaining mean structure. These ratios
describe these two signal cases. The conditional prototype uses the declared
noise covariance, and their measured ratios do not validate substitution of a
data-estimated covariance into the known-covariance law.

## Root and descendant tests

For each signal tree, the tested family consists of its root and its two
immediate internal children. For each descendant define a unit contrast e with
entries proportional to 1/n_left and -1/n_right inside its two children and
zero outside that node. The path

\[
A=(I-ee^T)X,\quad u=X^Te/\|X^Te\|,\quad
X(r)=A+\sigma eru^T
\]

leaves all outside rows fixed and tests the original selected descendant in
the full original hierarchy. Assertions verify reconstruction of X, nuisance
orthogonality, fixed outside rows and observed-history recovery. Fitting a new
tree on only the descendant's rows would use a different selection event.

Under the known spherical Gaussian model, the same conditional chi law applies
to this full-data descendant contrast under e^T mu=0. Conditional on the entire
selected merge history M, the family of three node identities is fixed. If each
conditional p-value is valid, the union bound gives familywise error at most
alpha by using alpha/3 at each node. The experiment applies this simultaneous
family rule. It does not implement the production adaptive Tree-BH and
pass-through traversal.

| Case and node | Populations present | Children | Accepted radii | Numerical p |
| --- | --- | --- | ---: | ---: |
| Clear root, 60 rows | Four populations | 15 / 45 | 999 | 0.001 |
| Clear child, 15 rows | One population | 1 / 14 | 0 | 1 |
| Clear child, 45 rows | Three populations | 15 / 30 | 999 | 0.001 |
| Moderate root, 45 rows | Three populations | 15 / 30 | 999 | 0.001 |
| Moderate child, 15 rows | One population | 1 / 14 | 12 | 0.923077 |
| Moderate child, 30 rows | Two populations | 15 / 15 | 0 | 1 |
| Registered null root, 200 rows | One population | Root split | 0 | 1 |

The root and mixed 45-row child of the clear case reject at 0.01/3. The
moderate root rejects. The pure children remain unsplit in these diagnostic
decisions. The mixed 30-row moderate child remains unresolved in substance
because no proposals enter its selected event. Its returned rank p=1 records
the available Monte Carlo resolution. No final cluster-count or ARI claim is
made from this one-level test family.

The registered 200-row null root has nominal chi p=9.40775e-13. Its conditional
rank p=1 results from zero accepted proposals and consequently has no power at
the current resolution. Its production covariance diagonal averages 0.985554
times the declared variance. This single case isolates a severe selection
effect despite a close average variance estimate. The root replay took 403.22
seconds. Individual 45-row and 60-row contrasts took approximately 12 to 20
seconds. These timings include the instrumented full-history comparisons.

## Numerical event certification

The existing interval certifier was applied to the two signal roots and the
scalar control with observed half-width 1e-10. Both signal cases returned
uncertified observed membership and conditional p bounds [0,1]. Their reported
reason is that the leading or retained spectral group is not uniformly
separated. The clear graph has four connected components and four eigenvalues
equal to one to floating precision. The moderate graph has three components
and three such eigenvalues. The certifier requires an isolated leading mode.
The native constructor retains a tied leading group, so these cases exceed
that certificate's present spectral assumptions.

The scalar control certifies the small observed interval. Its complete
conditional p bound is still approximately [3.56e-12, 0.999999999718] because
the remaining intervals and infinite tail retain unresolved probability mass.
This is a local membership certificate with an undecided statistical test.

The current certificate functions reconstruct e from the final root merge
inside `_exact_path` and `_exact_observed_radius_bounds`. Passing a descendant
path to them would therefore reconstruct a root contrast. The audit verifies
this explicitly. For example the pure clear child has radius 7.85504, while
the root-only reconstruction returns 244.94110. Descendant paths were used
only in numerical replay, with no root certificate assigned to them.

## Reproduction and limits

`conditional-small-benchmarks-20261104.py.txt` accepts `audit` or a registered
continuous case name. `conditional-small-certificates-20261104.py.txt` accepts
either signal case or `existing_scalar_reentry`. Run with
`uv run --no-sync --offline python`. Outputs use exclusive creation and live
under `raw/assets/conditional-small-benchmarks-20261104/`.

The three native continuous cases are single fixed benchmark datasets. Their
results assess compatibility and behavior, and supply no empirical false-split
rate bound. The numerical p-values remain subject to the event-certification
limitations above. The 200-row null has a root replay only. The one-level
descendant family and certification control are explicitly identified.

All seven numerical tests completed. Every reported rank p-value was
independently recomputed from its saved accepted radii. Prototype source hashes
matched across the runs. Twenty existing focused tests passed, both diagnostic
scripts compiled, and wiki lint passed over 232 pages. The summary is saved in
`raw/assets/conditional-small-benchmarks-20261104/summary.json`.

## Sources

- `benchmarks/shared/cases/gaussian.py`,
  `benchmarks/shared/generators/gaussian_cases.py` and
  `benchmarks/shared/util/case_inputs.py` define the actual cases and inputs.
- `raw/inbox/choir-cell-lines/tbs-recovery-protocol.json` defines the separate
  count models and fitted preprocessing.
- `benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py`
  provides the numerical conditional rank statistic reused by the scripts.
- `benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_kernel.py`
  and `diffusion_root_cells.py` define the existing exact-real certificate scope.
