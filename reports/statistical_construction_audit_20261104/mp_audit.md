# Marchenko–Pastur rank and PCA basis audit

2026-11-04. Current working tree at HEAD `8e0093d3`, including existing uncommitted edge-split documentation and eigenspace tests. The audit scope is spectral construction and its interface with the Wald rank.

## Findings

1. **At energy fraction 0.9 the MP count does not determine the tested rank when the basis is available.** The gate asks for the entire contrast dimension (30 for the cancer PC configuration). The worker still computes the raw MP count and nominal test dimension, then exposes every numerically positive PCA direction if the nominal dimension is positive. The kernel replaces that positive nominal rank with the shortest eigenvalue-ordered prefix covering 90% of the squared projected *tested contrast*, then extends through the cutoff tie group. This can go below the floor or above the MP count. Fixed-subspace sibling gates take their own path. Rank zero is handled before the adaptive selection.
2. **Continuous fitted whitening changes the spectrum being compared to the unit MP edge.** Full same-row empirical covariance nearly flattens the whitened spectrum, and complete tie retention can expand the default floor to the full available rank even with zero raw MP count. Diagonal same-row whitening gives a scaled sample correlation matrix in the single-block case, with off-diagonal structure and estimated scales. The two-leaf case makes the dependence exact up to ridge. Correcting the Gaussian centering convention would leave these mechanisms intact.
3. **The current edge basis and sibling basis differ.** Edge tests alternate leaf-label-sorted rows within each child between contrast and PCA halves. A missing edge basis produces the sentinel `Q=0, df=0, p=1` before the projected test. Siblings consume the all-leaf spectral context. The edge split separates PCA fitting from contrast rows for a fixed assignment, while hierarchy selection and the covariance used for the edge contrast still involve the tested data. The split alone establishes neither independent rank selection under energy mode nor an exact reference distribution.

## Source trace

Paths below are relative to the repository root. Within `tree_break_selection/hierarchy_analysis`, `statistics/projection` contains the spectral worker, dimension estimator and Wald kernel.

| Step | Source and observed behavior |
| --- | --- |
| Gate configuration | `decomposition/gates/orchestrator.py:288–303,433–471,499–512` validates the minimum and adaptive fraction. The edge floor defaults to 2 (`statistics/child_parent_divergence/child_parent_divergence_annotation/spectral_context.py:21`), while the generic MP estimator defaults to 1. The gate's adaptive fraction defaults to unset. `benchmarks/shared/runners/dispatch.py:346,406–411` forwards opt-in parameters. The internal-filter benchmark variants enable internal rows (`benchmarks/shared/runners/method_registry.py:285–308`). The factorial diagnostic explicitly sets 0.9 when energy is on (`benchmarks/diagnostics/calibration/diffusion_factorial.py:111–125`). The cancer PC reconstruction uses 0.9 (`reports/cancer_pc_statistics_debug_20261104/run_diagnostic.py:120–128,146–161`). |
| Row ownership | `child_parent_divergence_annotation.py:84–101,109–135` builds the returned all-leaf sibling context, then an independent *call* for the split-row edge context. `spectral_context.py:119–148` defines the label split. Internal rows can enter the former via the opt-in flag. `tree_estimator.py:175–236,360–388,416–425` supplies actual basis leaf indices and optional descendant internal state vectors. For the edge basis it replaces the working row list and computes its null mean and guarded within-child continuous covariance from the basis groups. |
| Whitening | `statistics/contrast_covariance.py:314–349,398–417,518–557,900–955` applies block maps at the task's empirical null state using ridge `1e-12` from `marchenko_pastur.py:98–124`. Bernoulli coordinates divide by the square root of `p(1-p)+ridge`. Categorical blocks drop one coordinate, solve through the simplex covariance Cholesky, and use contrast width. Continuous blocks use a diagonal variance divisor or full Cholesky solve. There is no additional cross-block whitening. The subsequent sample centering removes a common location shift while leaving the fitted scale in place. |
| Continuous covariance source | `tree_break_selection/tree/distributions.py:267–296` stores parent-total empirical scatter divided by `max(leaf_count−1,1)` and falls back to diagonal for undersized full blocks. The normal all-leaf spectral context retrieves that stored covariance (`tree_estimator.py:225–234`). The edge split basis can use pooled within-child residual scatter divided by `sum(n_i−1)` under the guard, with fallback to the original all-leaf parent total (`distributions.py:381–400,479–518`). The edge contrast separately resolves configured covariance from tree distributions (`tree_testing.py:67–99`). This can use a different estimated whitening map from the basis. |
| Spectrum and aspect ratio | `decomposition/backends/eigen/preparation.py:21–40` counts the actual post-whitening coordinates with finite strictly positive sample variance as `d`. Thus `d` can depend on the rows and can exceed the rank. `operators.py:8–24` centers columns and divides primal or dual centered scatter by physical row count `N`. `marchenko_pastur.py:125–159` passes basis leaf count `m` as both `effective_independent_rows` and `mp_threshold_rows`, even with augmentation. The label records a count, not a statistical effective sample size estimate. There is no independent estimate of noise scale. `projection_dimension_estimators.py:84–98` uses unit scale and strict `lambda > (1+sqrt(d/m))²`. |
| Available basis and ties | `projection_dimension_estimators.py:146–167` floors the raw count, caps by positive available directions, then extends the prefix through its numerical tie group. `decomposition/backends/eigen/projection.py:56–109` uses numerical positivity tolerance `eps × eigenvalue_count × max(lambda_max,1)`, caps the dual at `N−1`, and recovers normalized feature directions from dual sample vectors. It restores zeros on removed coordinates (`projection.py:112–121`). `eigenvalue_groups.py:8–30` joins adjacent ties transitively using relative tolerance `1e-10` and absolute tolerance `1e-12` times the spectral scale. The positive cap can put the actual dimension below the nominal floor. |
| Adaptive override | `marchenko_pastur.py:162–197` exposes the requested basis independently of the nominal rank when the latter is positive. `projected_wald_kernel.py:64–85,91–139` handles zero nominal rank first, then recomputes positive adaptive rank from ordered squared scores. It sets target 1 when total energy is at most `1e-24`; the final cutoff extends through ties (`projected_wald_projection_basis.py:10–36`). Basis eigenvalues do not scale the squared scores. The resulting reference df is the retained row count. |
| Zero-rank metadata | `marchenko_pastur.py:134–148` returns `mp_threshold_rows=N` for an augmented matrix with no active features, while successful augmented evaluation uses `m` (150–159). No threshold is evaluated in the zero-rank branch. A reproduced 2-leaf plus 1-internal constant matrix reports threshold rows 3 and effective rows 2. This is a metadata inconsistency. |

## Mathematical checks

**Fixed known Gaussian whitening.** Suppose a predetermined subtree contains `m` iid Gaussian leaf rows, the whitening map is fixed and makes the population covariance identity, and the active feature set of width `d` is fixed. Orthogonal rotation of the sample axis transforms the centered scatter into a sum of `m−1` independent standard Gaussian outer products. Thus the backend spectrum has the exact matrix law `W_d(m−1,I)/m`. The corresponding scale-matched asymptotic edge convention is `(sqrt(m−1)+sqrt(d))²/m`. The implemented edge uses `sqrt(m)` in the numerator and is higher. Replacing only the aspect-ratio denominator by `m−1` produces a different, higher edge because it ignores the backend divisor. Neither asymptotic edge is a fixed-level quantile. The MP bulk limit, Gaussian largest-eigenvalue fluctuations, and spiked-model detectability are separate results. [Primary source notes](../../raw/inbox/tbs-mp-primary-sources-20261104.md)

**Fitted full whitening.** Let `S=XcᵀXc/(m−1)` be the covariance of the very rows used for PCA, and let `LLᵀ=S+rho I`. Then the PCA covariance for a full-rank single continuous block is `((m−1)/m)L⁻¹SL⁻ᵀ`, whose eigenvalues are `((m−1)/m)lambda_i(S)/(lambda_i(S)+rho)`. At zero ridge this is exactly a scalar identity. An asymptotic deformed-MP law with a fixed nonidentity population covariance does not account for this fitted transformation. With diagonal same-row whitening instead, the exact zero-ridge expression is `((m−1)/m)R_sample`. For two Gaussian leaves all nonconstant standardized columns have identical centered row patterns up to sign. The rank is one and its lone PCA eigenvalue is `d/2`, below the implemented edge for `m=2`. Nevertheless the native rank is one, and the tested sibling contrast lies in that direction. The previous cancer PC pair reference mismatch therefore has a direct spectral mechanism independent of a raw MP spike.

**Augmentation and selected membership.** With fixed memberships and fixed internal averages expressible as linear combinations of the leaves, stack a row operator `A` whose first block is identity and whose remaining rows are the averaging weights. With fixed known whitening the centered scatter is `Xᵀ Aᵀ H_N A X / N`, where `H_N` removes the augmented mean. Its row weights matter even if `rank(Aᵀ H_N A)=m−1`. A leaf count alone cannot recover its spectrum. Branch-length-dependent weights and selected memberships add further conditioning. Even in leaf-only mode, selected labels, data-dependent active columns, grouped features, covariance fitting and ridge can change the null ensemble. Independent Gaussian rows, known isotropic covariance and growing dimensions with stable aspect ratio are a sufficient control, not conditions established by leaf-only input.

## Reproduced controls

Run from the repository root with the existing environment.

```bash
uv run --no-sync --offline python reports/statistical_construction_audit_20261104/mp_probe.py
uv run --no-sync --offline pytest -q tests/statistics/29_test_method_k_estimators_parity.py tests/statistics/55_test_eigenspace_invariance.py
```

The focused tests passed (37). Environment was macOS arm64, uv 0.9.5, Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, pytest 8.4.2. [Probe output](generated/mp/probe.json)

| Fixed Gaussian control | Verified result |
| --- | --- |
| Same 60×30 matrix, known identity whitening | Raw MP count 0, nominal rank 2, adaptive basis rank 30, eigenvalue range 0.100416–2.658316. |
| Same 60×30 matrix, fitted full whitening | Raw MP count 0, nominal rank 30 through one numerical tie group, eigenvalues within `1e-9` of `59/60`. |
| Same 60×30 matrix, fitted diagonal whitening | Raw MP count 0, nominal rank 2, adaptive basis rank 30, eigenvalue range 0.100213–2.488276. |
| Two rows and fitted diagonal whitening | Raw MP count 0, nominal and available rank 1, eigenvalue 14.9999999854. |
| Same 30-row PCA basis, positive nominal rank changed across 1, 2, 30 | Fraction 0.9 selects ranks 1, 1, 1 for contrast along the first PCA direction and 30, 30, 30 along the last. |

Top eigenvalue exceedances among 3,000 fixed-known-whitening Gaussian draws per shape were:

| `m × d` | Implemented edge | Centered scale-matched edge | Centered unscaled edge |
| --- | ---: | ---: | ---: |
| 6 × 30 | 158 | 294 | 22 |
| 30 × 30 | 274 | 373 | 187 |
| 60 × 30 | 286 | 368 | 243 |

These are repeatable bounded controls of the rank heuristic. They do not estimate false-split rates, calibration of the selected tree, or operating characteristics of the cancer PC pipeline.

## Accumulated claim verdicts

| Claim and location | Verdict |
| --- | --- |
| Primal and dual positive eigenvalues have identical scale (`wiki/analyses/local-marchenko-pastur-rule.md:38–48`). | Verified. Both divide by the physical spectral row count. |
| Centering convention and finite-sample MP exceedances are distinct (`local-marchenko-pastur-rule.md:133–158`). | Verified. The two centered conventions produce different exceedance counts. |
| The inferential basis uses leaf rows and MP determines the dimension (`local-marchenko-pastur-rule.md:50–59,148–151`; `manuscript/sections/method/edge_test.tex:145–176`). | Accurate for the default input type and nominal rank only. Current edge tests use the basis halves, siblings use all leaves, and opt-in internal rows still exist for sibling spectral input. Adaptive energy replaces the nominal rank. The manuscript's `d` cap also needs positive-rank and tie qualifications. |
| Continuous self-whitening can flatten the spectrum near `(m−1)/m` (`local-marchenko-pastur-rule.md:70–79`). | Verified for a dense same-row full block. Diagonal and guarded pooled covariance have distinct behavior. |
| The MP rule is a heuristic and finite local use is unvalidated (`manuscript/sections/method/assumptions_validation.tex:30–38,72`). | Supported. Known-white asymptotic results do not provide a selected, fitted, finite-subtree reference. |
| Controlled leaf-only selected PCA was calibrated in the settings recorded in `wiki/analyses/selected-pca-projected-wald-validation.md`. | Preserve as a historical bounded result. A same-pair basis can align fully with the tested contrast even for a fixed pair with identity noise, as shown by the existing cancer PC control. |

## Correction targets

- Document the nominal MP rank, available positive rank, realized contrast-energy rank and zero-rank sentinel separately. Describe the split edge basis and all-leaf sibling basis as distinct.
- Make augmented zero-rank threshold metadata consistent with the declared leaf-row convention, or represent a threshold that was never evaluated explicitly if the contract is revised.
- Treat calibrated spectral rank and calibrated Wald inference as separate questions. Any proposed law should match the actual centering divisor, whitening estimator, row operator, membership selection and rank choice. Establish control cases for fixed known, diagonal fitted, full fitted, and split-basis pooled or fallback covariance before making selected-subtree claims.

Production and test code were left untouched. The full suite and real-data regimes were outside this probe.
