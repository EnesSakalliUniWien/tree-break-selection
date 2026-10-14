# Repository naming and misleading-claims review

The findings below describe the **pre-fix review snapshot**. The subsequent
[corrections and verification](fixes.md) address all seven through clarified
claims, visible labels and metadata, preserving formulas and API identifiers.

Three user-requested Luna agents reviewed statistical terminology, geometry and
solver terminology, and documentation claims. The coordinating agent checked
each retained finding against then-current source and relevant callers. Seven
actionable issues were recorded. This is a targeted repository review, not a proof
that every name or statement in the repository is correct.

Snapshot: HEAD `f34f59cacd09b19c34fafa3a4a6e783dd3973f72`, including the existing
dirty and untracked work. No runtime source, public symbol, default or numerical
result was changed. Only this review and its wiki record were added. Historical
reports remain evidence for their original snapshots.

The code-review skill supplied review guidance, but its fixed-point diff and
issue-spec workflow did not apply to the explicitly requested repository-wide
naming audit. The agents reviewed present meaning against implementation.
Priorities indicate the risk of scientific or operational misunderstanding;
they do not imply a reproduced runtime failure.

## Statistical names and claims

### N1 — P1: fixed projection is described as data-independent testing

Claim: [fixed_subspace_annotation.py:1](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/fixed_subspace_annotation.py#L1)
calls the module “Data-independent fixed-subspace sibling gate annotation.”
The [hierarchy README:34–37](../../tree_break_selection/hierarchy_analysis/README.md#L34)
calls `fixed_global_chi_square` the full fixed-subspace chi-square reference.

Actual behavior: the coordinates/rank avoid learned parent PCA, but
[fixed_subspace_annotation.py:251–264](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/fixed_subspace_annotation.py#L251)
resolves covariance from the parent. The default `parent_total` policy returns
empirical subtree covariance, computed from the observed scatter with divisor
`n-1` in [distributions.py:265–289](../../tree_break_selection/tree/distributions.py#L265).
The [global test:219–220](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/fixed_subspace_annotation.py#L219)
then applies `chi2.sf(dot(z,z), df=dimension)`.

Consequence: readers may infer that covariance estimation and hierarchy
selection have also been made independent of the tested data. They have not.
This source trace establishes the terminology mismatch; it is not a fresh
numerical calibration experiment.

Suggested clarification: qualify independence as applying to projection/rank
choices, and describe the chi-square reference as nominal where its complete
null assumptions have not been established. Known covariance alone would not
resolve selection of the tested groups. Confidence: high.

## Documentation claims

### N2 — P1: README and class docstring describe a validation gate that is absent

The [README:107–124](../../README.md#L107) says sibling FDR is applied only when
the reference law is validated, and positive-dimensional same-data tests are
unsupported. [TreeDecomposition:56–58](../../tree_break_selection/hierarchy_analysis/tree_decomposition.py#L56)
similarly says unresolved calibration fails closed.

Actual behavior: [inflation_adjusted_sibling_tests.py:54–98](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/inflation_adjusted_sibling_tests.py#L54)
accepts `internal_admissible` and returns adjusted p-values. The
[decision implementation:841–851](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py#L841)
sets that status while explicitly recording
`same_selected_hierarchy_reference_law_unvalidated`.

Suggested clarification: internal support permits the empirical plug-in rule
and subsequent BH; absent support closes the focal gate. Reference-law validity
remains an open question rather than an implemented approval condition. Retain
the correctly qualified README discussion at lines 47–50. This contradiction
was already recorded in the earlier mathematical construction review and is
revalidated here as still unresolved. Confidence: high.

### N3 — P2: a current report calls available plug-in p-values calibrated

The [passthrough report:57–60 and 70–71](../passthrough_pvalue_normalisation_20260928/README.md#L57)
calls 345 values “calibrated sibling p-values” and describes a “calibrated
pre-BH” value. Its [p-value audit:7](../passthrough_pvalue_normalisation_20260928/pvalue_audit.md#L7)
uses the same description.

Actual behavior: [empirical_null_inflation_estimation.py:642–645](../../tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py#L642)
states that internal admissibility means availability and the plug-in value
does not establish selective error control. Numerical results and the report's
no-population-truth caveat are not challenged by this finding.

Suggested clarification: “empirical-inflation plug-in p-values” or
“inflation-adjusted pre-BH p-values.” Reserve “calibrated” as a validity claim
for a stated target and supporting evidence, or explicitly identify it as the
name of an adjustment operation. Confidence: high.

## Geometry, units and provenance

### N4 — P2: invariant/equivariant labels lack a declared transformation

`InvariantEquivariantSpace` (renamed `LeadingAxesSpace` in `leading_axes.py`, 2026-10-13),
its properties, and `decompose_invariant_equivariant_space` label the first
axis invariant and later axes equivariant. The
implementation lines 73–78
centers the matrix, performs SVD and selects the leading axes. It does not
declare a group action or transformation under which these particular claims
are to hold. The module disclaims a formal Cartan decomposition, which does
not define the remaining invariant/equivariant terminology.

Suggested clarification: state exactly which transformation and property are
intended, or explicitly identify these as operational labels for leading and
following PCA axes. No API rename is authorized or performed by this review.
This finding does not assert that PCA has no equivariance properties.
Confidence: high.

### N5 — P2: the NNLS target name hides feature-count normalization

`squared_standardized_euclidean` is the registered target name, but
[_squared_standardized_targets:182–195](../../tree_break_selection/tree/optimized_branch_lengths.py#L182)
computes `sum(standardized_difference**2) / n_features`.

The division is material when comparing branch-length units across feature
dimensions or reproducing the target from its metadata. The result is a mean
squared standardized difference, rather than the unnormalized sum suggested
by the name alone.

Suggested clarification: document the exact normalization in target metadata
and captions, using “mean squared standardized Euclidean distance” where prose
is appropriate. Preserve numerical behavior and existing identifiers until a
specific compatibility-aware rename is authorized. Confidence: high.

### N6 — P2: a plotting caption gives fitted branch lengths time semantics

[selected_nnls_fit_summary.py:341–342](../../applications/scrna/plots/selected_nnls_fit_summary.py#L341)
calls the fitted value “summed branch time along the real tree path.”
The same script calls `_squared_standardized_targets` at line 112 and sums
saved fitted edge lengths in `_tree_pair_distances` at lines 68–97. The
[NNLS module](../../tree_break_selection/tree/optimized_branch_lengths.py#L1)
explicitly states that these fitted lengths are not calibrated stochastic time.

Suggested clarification: “summed fitted branch length,” with the normalized
squared-distance target and its units described. This is a caption issue,
not evidence that the fitting objective itself is incorrect. Confidence: high.

### N7 — P2: backend metadata omits the coordinate implementation

[diffusion.py:308–350](../../tree_break_selection/space_separation/diffusion.py#L308)
fits and evaluates a pydiffmap kernel, then calls repository
`compute_diffusion_coordinates`. That helper constructs the symmetric
normalized operator and uses SciPy `eigh` at lines 188–195. Metadata records
only `backend: pydiffmap`; the benchmark runner forwards it in
`adaptive_diffusion` metadata at line 616.

This is incomplete provenance rather than a wholly false backend label.
The [space-separation README:26–31](../../tree_break_selection/space_separation/README.md#L26)
already distinguishes these coordinates from Markov right-eigenvector
coordinates, but standalone output metadata does not.

Suggested clarification: record the kernel backend, coordinate operator and
eigensolver separately. Confidence: medium-high.

## Qualified terms to retain

- `internal_admissible` explicitly means support-policy availability, not
  selective-tail validity. Improve surrounding contradictions rather than
  treating the status itself as an unqualified guarantee.
- The exact-F API states restricted independence, common-scale and weight
  requirements. The review found no reason to strip “exact” from that
  conditional mathematical statement.
- Independent Gaussian annotations explicitly distinguish retained rank from
  a fitted reference df (`independent_gaussian.py:338–344`). The legacy
  `Sibling_Degrees_of_Freedom` transport field alone is not an unqualified
  chi-square-law claim on that path.
- Numerical eigenvalue groups explicitly disclaim exact multiplicity
  certification. Current tie-expansion and symmetric-coordinate explanations
  should remain.
- “Fixed-topology NNLS” and “optimization” describe the implemented fit;
  neither alone asserts an exact optimality certificate or calibrated time.
- Recent independent Gaussian reports state matched denominators, unavailable
  outcomes and limits. This audit found no denominator overclaim in the
  inspected current reports and does not reclassify old evidence as current.

## Verification and boundary

All findings were checked by reading current code and the relevant data flow;
no new benchmark, mathematical proof or runtime regression is claimed. The
pre-edit documentation baseline passed `make wiki-lint` for 208 pages. The
review adds a wiki source page and append-only log entry; post-edit wiki lint
and local report-link checks verify those artifacts. No corrective rename,
implementation change, commit, push or dependency installation was performed.
