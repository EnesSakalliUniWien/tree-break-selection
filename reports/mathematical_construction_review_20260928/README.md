# Mathematical construction review — 2026-09-28

This is the review snapshot **before** the subsequent
[rank-invariance repair](/Users/berksakalli/Projects/kl-te-cluster/reports/mathematical_construction_review_20260928/rank_invariance_repair.md).
Its original numerical results and hashes remain preserved. The follow-up records
which findings were repaired and which mathematical questions remain open.

Three agents reviewed eigenvector geometry, calibration, and documentation/module reachability. The coordinating agent independently reproduced the numerical findings. **The construction has real basis-invariance defects and unresolved selected-calibration questions. Missing left eigenvectors of the current symmetric operators is not the explanation.**

This reviews the working tree at commit `f34f59cacd09b19c34fafa3a4a6e783dd3973f72`, including the existing uncommitted adaptive-cosine tie fix. Production code was not changed. Existing benchmark work was preserved. The results below distinguish demonstrated defects, model choices needing a contract, and questions requiring new statistical evidence.

## Reproducible findings

[reproduce_findings.py](/Users/berksakalli/Projects/kl-te-cluster/reports/mathematical_construction_review_20260928/reproduce_findings.py) produces [findings.json](/Users/berksakalli/Projects/kl-te-cluster/reports/mathematical_construction_review_20260928/findings.json), including source hashes, package versions, exact matrices, seeds, ranks, angles and numerical results:

```bash
.venv/bin/python reports/mathematical_construction_review_20260928/reproduce_findings.py
```

### P1 — Adaptive Wald rank changes the null law and can change a rejection under eigenspace rotation

[projected_wald_kernel.py:114](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py:114) chooses the shortest prefix capturing a requested fraction of the **tested contrast's** squared component energy. It then uses the ordinary rank-based chi-square tail at line 150. The reference helper explicitly requires a projection not selected from that contrast. See [projected_wald_reference_distribution.py:10](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py:10).

For `z=(3,0)`, eigenvalues `(1,1)`, candidate rank two, and energy fraction `0.90`:

| Basis for the same repeated eigenspace | Returned rank | Statistic | p-value |
|---|---:|---:|---:|
| Identity | 1 | 9 | 0.002699796 |
| Rotated 45 degrees | 2 | 9 | 0.011108997 |

The rejection changes at alpha `0.01`. Preserving the operator and its full eigenspace is insufficient when a basis-dependent prefix is selected.

One million independent Gaussian `N(0,I₂)` contrasts, seed `1907`, isolate the selection effect with a fixed identity basis:

| Nominal alpha | Adaptive-rank rejections | Adaptive rate | Monte Carlo standard error | Fixed-rank-two rejections |
|---|---:|---:|---:|---:|
| 0.001 | 1,525 | 0.001525 | 0.00003902 | 1,027 |
| 0.01 | 14,720 | 0.014720 | 0.00012043 | 10,070 |
| 0.05 | 67,820 | 0.067820 | 0.00025144 | 49,963 |

The vectorized calculation was checked against 128 calls to the actual kernel. These are independent contrast-level trials, **not** an estimate of the full pipeline's false-split rate. Core and generic benchmark defaults leave this option disabled; [pancreas_benchmark.py:239](/Users/berksakalli/Projects/kl-te-cluster/applications/scrna/pancreas_benchmark.py:239) enables `0.90` for edge and sibling testing.

Required correction: select a whole subspace independently of the tested contrast, or derive/validate the law for the actual adaptive selection rule. Keeping tied eigenvalues together alone does not establish a selected-null law for other data-dependent rank choices.

### P1 — Diffusion component truncation cuts repeated eigenspaces

[diffusion.py:184](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/space_separation/diffusion.py:184) drops the first eigenvector and retains exactly `n_components`, without preserving a repeated eigenvalue group. On the six-cycle, the normalized spectrum is `[-1,-0.5,-0.5,0.5,0.5,1]`. Permuting sample order by `[2,5,0,4,1,3]` and mapping the result back gives:

| Components, diffusion time 2 | Maximum change in matched pairwise distances |
|---|---:|
| One, cutting the 0.5 eigenspace | 0.21650635094611065 |
| Two, retaining that eigenspace | 9.71445146547012e-16 |

This can alter inputs to every downstream tree inference method using the shared diffusion function. It does not establish how often the earlier benchmark panel encountered the defect. A correction must preserve complete numerical groups and specify stationary-space removal for disconnected graphs; simply dropping one vector can also split a repeated eigenvalue-one space.

### P1 — Local PCA rank truncation has the same defect, including the rank floor

[projection.py:40](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py:40) slices eigenvectors by count. The [MP dimension rule:155](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projection_dimension_estimation/projection_dimension_estimators.py:155) floors a zero signal count to the configured minimum without checking ties. The [Wald basis helper:22](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_projection_basis.py:22) can truncate again.

For `X=[[1,0],[-1,0],[0,1],[0,-1]]`, covariance eigenvalues are `(0.5,0.5)`, the raw MP count is zero, and the default dimension-estimator floor is one. Exchanging feature columns and mapping back changes the physical projector from `diag(0,1)` to `diag(1,0)`: Frobenius difference `sqrt(2)`. This precedes optional adaptive Wald rank selection. A configured floor of two avoids this particular two-dimensional fixture, but does not generally protect a cutoff inside a larger tied group.

The existing adaptive-cosine tie fix is working in its own module. Its protection has not propagated to these separate PCA and diffusion cutoffs.

### P2 — The principal-angle diagnostic truncates unequal-rank spaces incorrectly

[recursive_pvalue_geometry.py:155](/Users/berksakalli/Projects/kl-te-cluster/benchmarks/diagnostics/path_b/recursive_pvalue_geometry.py:155) truncates both row bases to the smaller rank before forming their cross-Gram matrix. With parent basis `I₂` and child basis `[0,1]`, it reports cosine zero and normalized distance one. Swapping the parent's rows changes those answers to one and zero, although its subspace and projector are identical.

The full rectangular cross-Gram has singular value one: the child lies inside the parent. SciPy's full-subspace calculation returns angle zero. This affects diagnostic interpretation, not production gates. The live consumer is [mixed_null_signal_geometry_validation.py:1264](/Users/berksakalli/Projects/kl-te-cluster/benchmarks/diagnostics/path_b/mixed_null_signal_geometry_validation.py:1264).

Use full orthonormal bases for containment/principal angles. Also retain both ranks: zero principal angles do not imply equality when dimensions differ. The projector distance satisfies

`||Π_parent − Π_child||_F² = rank_parent + rank_child − 2 Σ cos²(theta_i)`.

Thus this example has containment angle zero but projector distance one. A diagnostic must state which question it measures.

### P2 — Hamming `k_neighbors=1` produces empty geometry

[diffusion.py:224](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/space_separation/diffusion.py:224) requests exactly `k` neighbors while querying the training observations, so the neighbor list includes self. The shared diffusion function erases the diagonal. For distinct rows `[[1,0,0],[1,1,0],[0,1,1],[0,0,1]]`, `k_neighbors=1`, time two and three components, all six distances are zero. More generally, the requested count includes self. The neighbor-count contract needs to count non-self neighbors or explicitly reject unsupported settings.

## A structural calculation question: whitening and MP

For one dense, full-rank continuous block with leaf-only spectral rows, the node's empirical covariance is used to whiten those same observations. The covariance uses divisor `m−1`, and the eigen backend uses divisor `m`. With zero ridge, if `X_c` is centered and `C=X_cᵀX_c/(m−1)`, whitening gives a subsequent covariance exactly `(m−1)/m I`. With positive ridge `epsilon`, its eigenvalues become `(m−1)/m * c_i/(c_i+epsilon)`, where `c_i` are empirical covariance eigenvalues. The production ridge `1e-12` has a small effect in these fixtures, but its effect can be large when a covariance eigenvalue is comparable to or below the ridge.

The actual `_process_node` path gives the following for 20 observations and two features, seed `104`:

| Input | Raw covariance eigenvalues | Whitened eigenvalues | MP signal count | Test dimension with floor 1 |
|---|---|---|---:|---:|
| Gaussian observations | 0.43633, 1.17425 | 0.95, 0.95 | 0 | 1 |
| First ten x-values shifted −10; remaining ten +10 | 0.49453, 102.79165 | 0.95, 0.95 | 0 | 1 |

Sources: [distributions.py:284](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/tree/distributions.py:284), [marchenko_pastur.py:106](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/projection/spectral/marchenko_pastur.py:106), and [contrast_covariance.py:937](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py:937).

This spectrum cannot discriminate the demonstrated covariance changes. It also makes nearly repeated eigenvalues structural rather than exceptional. The question is what information MP should identify after same-observation dense whitening, and what reference model supports that interpretation. This result does not directly extend to diagonal covariance, multiple blocks or augmented internal rows, and does not prove failure of the entire clustering procedure.

## Left/right eigenvectors, subspaces and multiplicity

The current operators must be distinguished before comparing their eigenvectors:

| Construction | Operator | Relationship |
|---|---|---|
| Adaptive cosine | `C = Y Yᵀ`, row-normalized input `Y` | Symmetric positive semidefinite; left/right eigenvectors can coincide |
| Local covariance | `ZᵀZ/m`, or dual `ZZᵀ/m` | Symmetric positive semidefinite; primal/dual recovery supplies feature/sample directions |
| Implemented diffusion | `S = D^(-1/2) W D^(-1/2)` after symmetrization | Symmetric; not necessarily positive semidefinite |
| Associated Markov operator | `P = D^(-1) W`, positive degrees | Generally nonsymmetric but similar to `S`; separate paired left/right vectors |

For these symmetric operators, algebraic multiplicity and geometric multiplicity agree. The associated `P` is also diagonalizable because it is similar to `S`; its nonsymmetry alone does not imply a multiplicity gap. A tolerance-based numerical eigenvalue group is not a proof of exact repeated algebraic multiplicity. A future general nonsymmetric operator could be defective and would require additional analysis.

SciPy's [symmetric eigensolver documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html) and [general left/right eigenvector definitions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html) support the operator distinction; the following relations are directly derived and checked in the saved probe.

Let `S U = U Λ` and `UᵀU=I`. Then

`R=D^(-1/2)U`, `L=D^(1/2)U`, `P R=R Λ`, `LᵀP=Λ Lᵀ`, `LᵀR=I`.

All paired left-left relationships are `LᵀL=UᵀDU`; all right-right relationships are `RᵀR=UᵀD^(-1)U`. Neither need be identity in the Euclidean metric. In their respective degree-weighted metrics, `RᵀDR=I` and `LᵀD^(-1)L=I`. Comparing only eigenvalues misses these metric differences.

![All right-right, left-left and left-right acute line angles for an illustrative irregular-degree Markov operator](/Users/berksakalli/Projects/kl-te-cluster/reports/mathematical_construction_review_20260928/paired_eigenvector_angles.png)

The figure uses a four-node weighted path with edge weights 1, 2 and 3, not benchmark data. Columns use ascending eigenvalue order. Absolute normalized dot products remove eigenvector sign ambiguity. Corresponding left/right line angles are approximately 25.24, 32.31, 32.31 and 25.24 degrees; off-diagonal left/right pairs are orthogonal because of biorthogonality. The JSON includes every Gram entry and both paired-subspace principal angles. Residuals for the eigenvector and duality identities are below `2e-15` in this run.

For a complete block `B`, distinguish the symmetric orthogonal projector `Π_B=U_B U_Bᵀ` from the Markov spectral projector `E_B=R_B L_Bᵀ`. The latter is generally oblique in the Euclidean metric. For isotropic Gaussian `z`, `||E_B z||²` generally has a weighted chi-square law determined by `E_BᵀE_B`, not simply chi-square with block rank. Changing projectors therefore requires revisiting calibration.

The current diffusion coordinates use `U`, not Markov right vectors `R`, and clip negative eigenvalues to zero. This is an explicit symmetric spectral geometry, but it is not ordinary Markov-power geometry. A three-node path at time two has implemented distances below `2e-30`, while two Markov-power row distances are `1.224744871`. The intended operator, metric, signed-spectrum treatment and stationary-space treatment need a clear contract; this comparison alone does not justify replacing the implementation.

GO feature loadings already compute `V=YᵀU Λ^(-1/2)` for positive eigenvalues at [adaptive_diffusion_go_ic_terms.py:30](/Users/berksakalli/Projects/kl-te-cluster/applications/endotypes/analysis/go_ic/adaptive_diffusion_go_ic_terms.py:30). These are paired **sample/feature singular vectors**, distinct from left/right eigenvectors of one square operator. That construction is correct. Individual axis signs and rankings inside a repeated space remain nonunique; complete-block feature energy or projectors would avoid that interpretation ambiguity.

Principal angles between disjoint eigenspaces of one symmetric operator are already 90 degrees. More informative comparisons concern matched subspaces across fitted operators, perturbations or tree contexts, with sample/feature alignment and an explicit inner product. The standard [principal-angle API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.subspace_angles.html) compares full column spaces.

**“Shadows” has no mathematical definition in the inspected project.** Coordinates `U_Bᵀz`, a projected vector `Π_Bz`, squared projection energy, spectral contribution `U_BΛ_BU_Bᵀ`, and an oblique projected vector `E_Bz` are distinct objects. The review does not silently choose one as the meaning of that term.

## Why documentation is ambiguous

| Conflict | Live evidence and effect |
|---|---|
| Active empirical p-values versus unsupported selected law | [README opening:15](/Users/berksakalli/Projects/kl-te-cluster/README.md:15) and line 47 describe the restored plug-in rule correctly. [README:108](/Users/berksakalli/Projects/kl-te-cluster/README.md:108), line 120 and [TreeDecomposition:56](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/tree_decomposition.py:56) still imply a validated law is required and unresolved cases always close. Executable code computes supported plug-in tails and applies BH while labeling the reference unresolved. |
| Annotation ownership | [hierarchy_analysis/README.md:10](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/README.md:10) assigns annotation preparation to `TreeDecomposition`; its [constructor:76](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/tree_decomposition.py:76) requires completed annotations. |
| “Invariant/equivariant” without a declared transformation | `invariant_equivariant.py:3` (renamed `leading_axes.py`, 2026-10-13) labels the first standardized PCA axis invariant and following axes equivariant. The computation is a PCA partition; it does not establish invariance for a declared group/action. This is an interpretation gap, not evidence that SVD is numerically wrong. |
| Missing coverage cited as present | [recursive-pvalue-geometry source:79](/Users/berksakalli/Projects/kl-te-cluster/wiki/sources/recursive-pvalue-geometry-20260606.md:79) cites nonexistent `tests/validation/80_test_recursive_pvalue_geometry.py`. No current test calls the angle helper. Wiki lint checks frontmatter paths, not every path in prose. |
| Retired code described as current | [wiki index:716](/Users/berksakalli/Projects/kl-te-cluster/wiki/index.md:716) describes an importable/registered legacy package; its [source page:20](/Users/berksakalli/Projects/kl-te-cluster/wiki/sources/legacy-c2ef9a69-method-package-20260616.md:20) records retirement, and executable searches find no package or registration. |

The common cause is incomplete propagation of behavior changes among README sections, class documentation, manuscript/wiki interpretations and historical index entries. Passing structural wiki lint does not verify mathematical claims or executable coverage. A repair should distinguish current behavior, assumptions, proof status and historical evidence explicitly. This review records conflicts without rewriting their history.

## Dead-module assessment

| Candidate | Reachability evidence | Assessment |
|---|---|---|
| [Diagnostic distributional-action enrichment](/Users/berksakalli/Projects/kl-te-cluster/benchmarks/diagnostics/analysis/distributional_action.py:23) | No incoming static import, symbol call, registry entry, CLI or test found; wiki advertises an opt-in helper | Unused internally; decide public retention before deletion |
| `selected_gaussian_radial.py`, `selected_gaussian_hierarchy.py`, `selected_gaussian_diffusion.py` | Focused tests, mutual imports, explicit non-production status | Intentional tested prototypes |
| `bootstrap_consensus.py` | Documented public API and focused tests; no application/benchmark consumer found | Optional API, not proven dead |
| `invariant_equivariant.py`, core `distributional_action.py` | Application and test consumers | Live |
| `legacy_methods/commit_c2ef9a69` | Absent executable package/registry; historical wiki/raw material retained | Already retired |

The existing Vulture command exited zero with no findings. A separate AST import scan across package, applications, benchmarks, scripts and tests found no core module lacking incoming imports. These bounded static checks cannot prove the absence of external consumers. No deletion is justified solely by the evidence above.

## What is correct, and what remains unproved

- Under independent observations and the stated covariance model, sibling mean-contrast scaling `1/n_left+1/n_right` and nested child-parent scaling `1/n_child−1/n_parent` are consistent. Their implementation is at [contrast_covariance.py:757](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py:757).
- For an isotropic **Gaussian** standardized contrast and a fixed orthonormal rank-`k` row basis, squared projected norm has a chi-square-`k` law. PCA eigenvalues do not need to be used as extra null weights under those assumptions. Isotropy alone is insufficient for an exact Gaussian chi-square law.
- The independent common-scale F API is separate and guarded by fixed unit weights, common scales and disjoint observation ownership. The caller must still justify independence and distributional assumptions. These checks do not certify selected same-tree calibration.
- The historical selected-Gaussian radius `40.01`/`[40,infinity)` underflow and narrow hierarchy-interval errors are repaired in live source and pass the relevant regressions. They should not remain on the unresolved numerical-defect list.
- The active supported empirical sibling rule remains `chi2.sf(T/(reference_scale*c_hat), df)`. [empirical_null_inflation_estimation.py:642](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py:642) distinguishes support availability from selective validity. No-support positive-dimensional cases are withheld, while support-threshold enforcement defaults to false. A scalar inflation estimate does not establish the selected tail when topology, covariance, projection, rank, edge selection and calibration records use the same observations.
- [Traversal-aligned BH:156](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/fdr_annotation.py:156) operates separately by depth and explores edge-reachable candidates. Per-depth BH is not a general whole-tree guarantee: six independent all-null families tested at `0.01` have any-rejection probability `1−0.99^6=0.0585198506`. This is a mathematical counterexample to an unrestricted guarantee, not an estimated repository error rate.

## Open questions and calibration improvements

1. **Operator and metric:** should each diffusion analysis measure the implemented symmetric geometry or Markov-power geometry? How should negative modes, disconnected stationary spaces and all left/right pairs be treated? Which precisely defined quantity is meant by “shadows”?
2. **Subspace selection:** what is the intended MP signal after same-observation dense whitening? How should rank floors, nearly repeated eigenvalues and complete tied spaces interact across every stage? What group/action supports the invariant/equivariant terminology?
3. **Selected reference law:** what law accounts for covariance estimation, regularization, PCA/rank selection, tree construction, branch-time factors and edge-based calibration selection? Eigenvector angles can describe these contexts but do not supply a valid tail by themselves.
4. **Calibration records:** are stopped/null-like records generator-null? How are contamination, overlapping observations, weights, estimated inflation and tail uncertainty accounted for? A large number of dependent nodes is not a large number of independent replicates.
5. **Error target:** distinguish local false sibling splits (`theta_left=theta_right`), a global-null final partition with more than one cluster, and false discoveries among traversal-used splits. `estimated_K > true_K` alone does not identify all erroneous splits; mixed children need generator-defined truth.
6. **Validation:** add sample/feature permutations and rotations within repeated spaces, unequal-rank principal-angle cases, and fixed/independent/same-data projection ablations. Then use independent held-out datasets and seeds to rebuild every fitted stage, reporting null error and signal recovery with uncertainty. Keep unsupported and failed outcomes separate.

The existing [calibration validation protocol](/Users/berksakalli/Projects/kl-te-cluster/reports/calibration_restoration_20260909/validation_protocol.md) is a useful starting point, not completed confirmatory evidence. Geometric correctness should be repaired before promoting extra angle-based calibration strata, because choosing strata on the same tested data adds another selection event.

For the earlier [11-case eigenvalue-block comparison](/Users/berksakalli/Projects/kl-te-cluster/reports/eigenvalue_block_comparison_20260926/README.md), retain the separation between eigenvalue geometry, rooted trees and final partitions. Its 752 attempts and frequent one-cluster outcomes are descriptive matched-case evidence, not calibration evidence. This review did not rerun that panel or attribute each collapsed result to the newly reproduced defects. A later matched rerun should quantify which blocks and inference methods change after concrete repairs.

## Verification and scope

The agents passed **184 existing targeted tests** across eigendecomposition, diffusion, GO loadings, Wald/reference calculations, empirical inflation, covariance, traversal, branch lengths and selected-Gaussian replay. The coordinating probe independently reproduced all five numerical defects/edge cases and the dense-whitening observation. Passing existing tests therefore exposes missing adversarial coverage rather than disproving the counterexamples.

Commands, all with `.venv/bin/python -m pytest -q -p no:cacheprovider`:

```text
# 33 passed
tests/statistics/48_test_space_separation.py
tests/core/test_adaptive_diffusion_runner.py
tests/applications/endotypes/pipelines/test_adaptive_diffusion_go_ic.py

# 4 passed
tests/statistics/38_test_eigen_backend.py

# 25 passed
tests/statistics/51_test_selected_gaussian_radial.py
tests/statistics/52_test_selected_gaussian_hierarchy.py
tests/statistics/49_test_projected_wald_reference_distribution.py

# 122 passed
tests/statistics/35_test_empirical_null_inflation_estimation.py
tests/statistics/39_test_projected_wald_statistics.py
tests/statistics/42_test_contrast_covariance.py
tests/statistics/43_test_traversal_aligned_sibling_fdr.py
tests/statistics/45_test_edge_gate_math_contract.py
tests/statistics/46_test_continuous_covariance_numerical_psd.py
tests/statistics/53_test_selected_gaussian_branch_lengths.py
tests/statistics/54_test_selected_gaussian_diffusion.py
```

Environment: Darwin arm64, Python 3.11.12, uv 0.9.5, NumPy 2.3.4, SciPy 1.16.3, pytest 8.4.2. Manifests, lockfile, Makefile, test setup and installed executables were inspected before writing review artifacts. The probe was linted; wiki lint was run before and after adding synthesis. The figure was rendered and visually inspected. No full `make check`, full benchmark rerun, production correction, deletion, dependency installation, commit or push was performed in this review.

External reference pages were captured as concise source notes under [raw/inbox/eigensystem-reference-notes-20260928.md](/Users/berksakalli/Projects/kl-te-cluster/raw/inbox/eigensystem-reference-notes-20260928.md). The browsed current SciPy manual is newer than the installed runtime; all numerical claims above were reproduced with the recorded installed versions.
