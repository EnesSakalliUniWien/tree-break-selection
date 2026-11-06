# Contrast, whitening, projection and raw reference laws

2026-11-04. Audit of the live working tree. Scope covers raw child-parent and
sibling test statistics. Counts from selected cancer surrogate nodes describe
one dependent tree. The Gaussian controls use fixed groups and repeated
independent data sets.

## Trace through production

1. The tree stores leaf-count weighted node means and continuous sample
   covariances. Child-parent testing alternates each child's leaves by
   stringified label, starting with a test row, then concatenates the test
   rows from all siblings into the parent mean. These are nested samples.
   Sibling testing uses each child's full stored mean and count. When groups
   are fixed and rows independent with a common covariance Σ, the variance
   factors are 1/n_child−1/n_parent and 1/n_left+1/n_right, respectively.
   For a binary parent with mean weighted by counts,
   Δ_edge=(n_other/n_parent)Δ_sib. The two standardized edge vectors are
   opposites if using the same whitening and branch multiplier. They are
   dependent evidence. In heterogeneous fixed groups, use
   Σ_left/n_left + Σ_right/n_right for the sibling covariance and the squared
   corresponding weight for the edge. With arbitrary independent rows the
   actual covariance is the sum of squared contrast coefficients times
   their row covariances. Selected membership changes the sampling law of
   the contrast itself.
   Sources: `child_parent_divergence_annotation.py:84–119`,
   `spectral_context.py:119–148`, `contrast_covariance.py:167–191,757–771`,
   `pair_observations.py:38–63`, `tree/distributions.py:315–359`.

2. The continuous test policy defaults to guarded within-child with minimum
   8 leaves per child. The guard requires a binary parent, both count minima,
   and residual df at least the width of each dense block. It pools child
   sample covariances with weights n_i−1. It falls back to full subtree
   parent-total covariance if the guard fails. Memory limits and small
   child subtrees can force diagonal storage. Pooling degrades to variances
   if either child's covariance is diagonal. Parent-total includes the
   between-child term. Guard success does not establish a well-conditioned
   covariance estimate or independence from the contrast.
   Sources: `tree/distributions.py:25–32,94–112,267–300,381–476`.

3. Edge PCA takes only basis rows for the scatter and centering mean.
   Its continuous whitening uses pooled basis-row residuals, with the
   helper's default guard and full parent-total fallback. This policy is
   separate from the configurable Wald covariance policy. The edge Wald
   contrast always resolves covariance from the full tree. Thus a row split
   alone neither isolates the topology nor makes all nuisance quantities
   independent of the edge contrast. Sibling PCA uses all parent rows and
   full parent-total covariance, regardless of the test's pooled policy.
   Basis and Wald coordinates can therefore use different Cholesky maps.
   An orthoprojector in truly standardized test coordinates still gives the
   fixed-basis law for any independent orientation. Different maps do
   change which original-coordinate contrast the PCA directions select.
   A two-leaf binary parent has no held-out basis rows. Its edges receive
   Q=0, df=0, p=1, while its full-data sibling test can have rank one.
   Sources: `projection/spectral/tree_estimator.py:224–233,361–387,401–424`,
   `tree/distributions.py:479–518`, `tree_testing.py:63–80`,
   `record_collection.py:78–135`.

4. Contrast whitening is a block Cholesky solve or division by diagonal
   standard deviations. Bernoulli uses p(1−p). Grouped categorical uses
   diag(p)−ppᵀ after dropping one coordinate, preserving within-block
   multinomial dependence. Cross-block covariances are absent from the
   model. Exact Gaussian chi-square requires correct covariance on the
   retained subspace. For Bernoulli or categorical plug-in whitening the
   comparable claim is asymptotic under a fixed-group CLT, adequately
   populated categories, and consistent nuisance estimates. Boundary
   probabilities, correlated blocks, tiny samples, and selected groups
   need separate laws. The fixed 1e−12 ridge and dense jitter stabilize
   the solve but change its covariance. Even with an otherwise known Σ,
   whitening with L Lᵀ=V+ridge·I gives projected null covariance
   R L⁻¹ V L⁻ᵀ Rᵀ. Its Gaussian squared norm is a weighted chi-square
   with those covariance eigenvalues, not generally chi-square with a
   changed df. Near-singular directions can also be dominated by the ridge.
   Sources: `contrast_covariance.py:49–75,194–233,840–895,900–973`.

5. Normalized edge variance multiplies by 1+ℓ_edge/mean(ℓ). Default edge
   variance omits this. Projected siblings always use
   1+(ℓ_L+ℓ_R)/(2 mean(ℓ)), and fixed-subspace siblings omit it.
   Under an iid fixed-group null with truly known sampling covariance,
   fixed external multiplier b≥1 and negligible ridge, Q becomes
   χ²_k/b. A genuine variance model with covariance b times the baseline
   would restore χ²_k. Data-estimated linkage or NNLS length creates
   dependence and the factor alone specifies no such model. Pointwise Q
   decreases under scalar inflation for fixed whitening geometry in the
   zero-ridge limit. That does not establish conditional calibration.
   Sources: `contrast_covariance.py:774–803`, `tree_testing.py:59–61,81–97`,
   `sibling_divergence_test.py:34–55`,
   `fixed_subspace_annotation.py:201–234`.

6. The kernel reports Q=‖Rz‖², scale 1, df=retained row count. Eigenvalues
   validate and delimit the basis. The kernel does not divide projected
   coordinates by their PCA eigenvalues. Ties extend a prefix to the
   end of a numerical eigenvalue group, preserving the subspace under
   rotations within that group. The sibling nominal rank is derived from
   full-data per-node spectral dimensions: rounded geometric mean of the
   two child dimensions, one positive child directly, or the parent's
   dimension when both are zero. It is capped by the parent test dimension
   before any adaptive energy replacement. Adaptive energy can override
   the nominal target in either direction using the whole returned
   positive candidate basis. Rank is bookkeeping, and the same number is
   then used as the chi-square parameter. Direct low-level callers must
   satisfy orthonormality themselves. The basis helper checks availability
   and eigenvalues but accepts scaled or duplicated rows. Internal calls
   use eigensolver rows. No malformed production PCA output was observed.
   Sources: `projection_dimensions.py:44–110`,
   `projected_wald_kernel.py:60–85,102–170`,
   `projected_wald_projection_basis.py:10–52`,
   `projected_wald_reference_distribution.py:58–96`,
   `tree_break_selection/eigenvalue_groups.py:8–30`.

7. The gate orchestrator requests the feature-space candidate width and
   forwards the same energy fraction to both projected tests when enabled.
   Its default energy fraction is unset. The cancer surrogate diagnostic
   uses 0.90 and projected-Wald inflation, as does the base pancreas
   application. The registered continuous guarded benchmark uses the
   fixed-coordinate BH sibling path. Branch-time NNLS presets explicitly
   enable normalized edge time. Fixed-coordinate and block alternatives
   bypass PCA and adaptive sibling rank, while keeping their own
   estimated covariance, selection and multiplicity requirements.
   Sources: `orchestrator.py:315,450–520`,
   `benchmarks/shared/runners/method_registry.py:68–78,113–125`,
   `applications/scrna/pancreas_benchmark.py:240–264`,
   `reports/cancer_pc_statistics_debug_20261104/run_diagnostic.py:117–130`.

## Exact restricted null controls

The formulas in this section use zero ridge and no extra branch variance.
The Gaussian groups and row assignments are fixed before observing values.

**Known covariance, fixed or independent projection and rank.** For z
standard Gaussian conditional on the selected R and k, with R Rᵀ=I_k,
Q is χ²_k. A random rank chosen independently of the test contrast is
permitted conditional on that rank and basis. For any fixed or independent
projection and Gaussian z with conditional mean zero and arbitrary known
covariance, the weights in the generalized chi-square law are the
eigenvalues of the projected covariance. A genuine rank deficiency requires
an explicit support and covariance rule.

**Same-row parent-total covariance.** For fixed iid Gaussian groups, take the
centered residual coordinate matrix Y with ν=N−1 rows and p≤ν columns.
The normalized mean contrast is cᵀY for a fixed unit vector c in centered
row space. With YᵀY/ν=LLᵀ, the columns of Y L⁻ᵀ/√ν are orthonormal
and their frame is Haar conditional on YᵀY by left rotational invariance.
Projection of that frame onto any k orthonormal directions fixed conditional
on YᵀY yields Q/ν~Beta(k/2,(ν−k)/2), for 0<k<ν, and a point mass at 1
if k=ν. This applies to a fixed nested test subset coefficient embedded in
the full parent data, but does not cover projection depending on the
residual frame, adaptive contrast rank, ridge, or diagonal approximations.
For one column at N=2 the limit is a point mass at 1. In p independent
columns with coordinatewise variances and N=2, unregularized whitening
produces signs. Fitting the leading PC from those same rows captures all
p units of energy before branch scaling.

**Within-child covariance and full-group Gaussian contrasts.** Within-child
residuals are orthogonal to the full-group mean contrast. For a dense
single block, W=ν S is Wishart_p(Σ,ν) and independent of the mean contrast
standardized for its sample sizes. After standardization by true Σ, rotate
the numerator to the first axis independently of the isotropic Wishart.
The inverse-Wishart quadratic form is its χ²_p length divided by an
independent χ²_(ν−p+1) Schur complement. Therefore the full projected
statistic is νp/(ν−p+1)·F(p,ν−p+1). A fixed k-dimensional original-coordinate
projection followed by estimating and inverting its k-dimensional
covariance has the analogous k formula. An arbitrary k-row projection
after estimated full Cholesky whitening is a different operation and does
not inherit that k-dimensional Hotelling law in general.

**Held-out edge contrast with full-data within-child denominator.** Even
with fixed groups, the test subset contrast coefficient has a component
inside the full residual subspace. For two groups of size 6, 3 test rows
each, squared residual overlap is ρ²=1/2. In a scalar Gaussian control
Q=ν(√(1−ρ²)G₀+ρG₁)²/(G₁²+W), with independent standard normals G₀,G₁ and
W~χ²_(ν−1). Cov(Z², residual SS)=2ρ². A within-child covariance estimated
from residuals around the test subset means could satisfy mean/residual
independence under the Gaussian fixed-group model. The current denominator
includes full-child residuals.

**Same-data direction and rank.** For a frozen two-row comparison with
known identity covariance, its rank-one PCA direction is z/‖z‖ and
Q=‖z‖²~χ²_p. Even a fixed orthonormal basis with adaptive energy rank
breaks the nominal law. If m candidates have iid standard normal projected
coordinates, their squared energy fractions have Dirichlet(1/2,…,1/2)
distribution independent of the total H~χ²_m. For selected prefix K
after tie extension, Q=H·∑_{j≤K} fraction_j, so the exact restricted law
is the mixture over the Dirichlet fractions that lead to that K. Conditioning
on K changes the distribution of the fractions. The exact special-case
law gives no calibration for a fitted basis or selected hierarchy.

Primary research notes are in
[raw/inbox/tbs-reference-law-primary-sources-20261104.md](../../raw/inbox/tbs-reference-law-primary-sources-20261104.md).

For Bernoulli siblings under a fixed-group common success probability,
conditional on the pooled count, the left count is hypergeometric. For a
fixed categorical block the left category counts are multivariate
hypergeometric given the pooled category counts. The same applies to a
fixed binary edge's test subset and its complement. The plug-in probability
is the pooled count fraction, and conditioning fixes it. The tangent
dimension of a positive K-category block is K−1 because the probabilities
sum to one, with smaller covariance rank when some categories have zero
probability. For a scalar Bernoulli block it is one at interior p and zero
at p∈{0,1}. The ridge restores invertibility, without restoring random
variation in directions that are structurally absent. These conditional
count laws are discrete even with correctly specified probabilities.
They are distinct from the multiblock projected law after fitting,
correlation, adaptivity, or tree selection.

## Reproducible checks

Run from the repository root with the existing environment.

```sh
uv run --no-sync --offline python reports/statistical_construction_audit_20261104/reference_probe.py
uv run --no-sync --offline pytest -q tests/statistics/{22*,34*,36*,39*,42*,45*,46*,49*,55*,58*}
```

The probe uses 20,000 datasets per control and seed 20261104. It replays
the production rank/statistic kernel for all adaptive draws and for 25
rank-one fitted directions. Covariance controls evaluate the stated
Gaussian formulas without ridge. The forced within-child controls use
6+6 leaves to expose small residual df. The production default guard would
select parent-total covariance at this size. A lowered guard would still
pool only variances for p=6 because each child has residual df 5. The test
suite separately exercises
the production contrast covariance and Cholesky paths.

| Fixed-group control | Nominal chi-square rejects at .05 | Restricted exact-law rejects at .05 |
|---|---:|---:|
| Known covariance and fixed k=3 | 989 | 989 |
| Parent-total and fixed k=3, ν=11 | 305 | 972 |
| Within-child full covariance, full p=6, ν=10 | 9,872 | 1,007 |
| Known covariance, rank-one direction fitted to the two tested rows | 20,000 | 1,001 |
| Known covariance, fixed basis, 0.90 energy rank with m=8 | 1,562 | — |

The corresponding adaptive count at .01 is 349. Conditional rank results
include 105/437 rejections at .05 for k=4 and 408/8,255 for k=8. In the
scalar held-out edge, observed Cov(Z², residual SS)=1.0554, with exact
expectation 1. The nominal χ² and inapplicable F references rejected
1,331 and 760 times at .05. For the fitted-pair direction, all 20,000 draws
rejected at .01 using χ²₁ and 194 rejected using χ²₃₀. These contrasts
separate covariance, basis fitting and rank selection mechanisms.

A second 20,000-draw control (seed 20261105) uses 8+8 fixed leaves and
p=6, which passes the default guard and preserves dense child covariance.
For 25 datasets, native distribution population and covariance resolution
returned the analytic pooled matrix. Both native test wrappers with
a fixed full basis and no branch multiplier matched the unregularized
quadratic forms to maximum absolute Q error 7.42e−10. The edge replay
uses native alternating row selection. For full sibling means, the
nominal chi-square rejected 6,540 times at .05 and the exact
14·6/9-scaled F(6,9) rejected 1,027 times. With the alternating edge
means and the same pooled denominator, the corresponding counts are
4,060 and 241. The full-means Hotelling reference is inapplicable to the
held-out edge. Its mean Q was 9.021, compared with 12.042 for full means.
Analytically the squared residual overlap is again 1/2. Conditional on
the residual scatter, the contrast's residual part has mean whitened
energy p by row rotational invariance, and the independent full-means
part has mean energy νp/(ν−p−1)=12. The cross term has mean zero. Thus
the expected held-out Q is (6+12)/2=9. This is a fixed-basis replay of
the default covariance path, with PCA, branch and hierarchy selection
excluded by construction.

The saved selected cancer tree records contain 183 mean-null two-leaf
pairs with Q·branch factor in [29.9999647804, 29.99999999997] and 182
raw p≤.01. In the accompanying saved intervention, known noise with
native basis yields 84/183 and known noise with a fixed coordinate
or full ambient norm yields 0/183. Branch scaling cannot explain
the unscaled energy near 30. None of these 183 became an eligible
native sibling rejection. The selected-tree counts are descriptive.

The focused baseline passed 92 tests. It covers nested overlap scaling,
branch-policy routing, simplex covariance, regularization, projection
width/df, ties and adaptive option forwarding. The edge calibration test
in `tests/statistics/58_test_edge_null_calibration.py` includes an independent
topology Bernoulli control and an upper bound of twice nominal on its Gaussian
case. This is different evidence from exact raw reference validity.
Environment: Darwin arm64, uv 0.9.5, Python 3.12.11, NumPy 2.3.4,
SciPy 1.16.3 and pytest 8.4.2. The probe passed Ruff.

## Claim verdicts and correction requirements

| Claim | Verdict |
|---|---|
| Nested 1/n_child−1/n_parent and disjoint 1/n_left+1/n_right | Algebra correct under fixed membership, independent rows and common covariance. Selected or heterogeneous samples need their own law. |
| Full-data fallback is interchangeable with pooled covariance | Assumption mismatch. It contains between-child variation and has a different dependence structure and reference even under Gaussian fixed groups. |
| Held-out edge PCA makes inference independent | Assumption mismatch. The topology is already selected, the fallback may use test rows, the Wald covariance always uses full rows, and energy rank may reuse the tested contrast. |
| Numerical eigenvalue groups make retained eigenspaces rotation-invariant | Correct for complete groups, up to floating-point ties. This is a geometric property. |
| Original MP count equals Wald df | Incorrect. MP count is the number above the estimated edge. Floor and positive availability determine a nominal target, the sibling route can derive a different target, contrast energy can replace it, and the final tied prefix determines the realized rank. The kernel uses that final rank as nominal df. |
| Rank k establishes χ²_k | Assumption mismatch. Exact only for standard Gaussian test coordinates conditional on an orthonormal projection and rank. Scaled-Beta, Hotelling F and adaptive Dirichlet controls establish distinct restricted laws. |
| Low-level kernel validates the orthonormal basis contract | Input-validation gap. It computes with scaled or duplicate rows supplied by a direct caller. Production eigensolver output has a separate orthonormality contract. |
| Calibrated p-value under branch scaling, ridge, or estimated whitening | Unverified as a general law. Explicitly identify and justify the full projected covariance and test whether rank and covariance are independent of the test contrast. |

A production calibration claim requires the null target and the precise
conditioning regime, including whether groups are fixed, externally
discovered with independent test data and assignments, or selected from the
same observations. It also requires a reference for the actual covariance
path, branch factor, PCA construction and adaptive rank. Correcting df alone
would not address those dependencies. The restricted controls are
diagnostics for individual mechanisms. They supply no replacement
selected-tree p-value or traversal guarantee. No production changes or
large reruns were made.
