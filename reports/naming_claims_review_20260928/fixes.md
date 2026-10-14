# Corrections to the naming and claims review

All seven review findings are addressed through documentation, visible labels
and explicit metadata. The existing numerical formulas, method identifiers,
function names and data-column names are preserved. Retained PCA identifiers
are explicitly documented as operational labels, not transformation claims.
These corrections do not establish a new statistical calibration guarantee.

The user requested Sol with appropriate reasoning effort and the cavecrew
skill. A Sol agent at high effort performed a two-file cavecrew-builder task
for diffusion provenance and NNLS target metadata. The coordinating agent
handled the connected documentation, captions and existing tests. The same Sol
agent reviewed the combined changes using the cavecrew-reviewer format; this
was not an independent review of its own two-file implementation.

## Resolution by finding

| Finding | Correction |
|---|---|
| N1: data-independent fixed-subspace claims | Module and hierarchy README now limit the fixed choices to projection/rank and explicitly retain same-data covariance and selected-tree caveats. Chi-square tails are described as nominal until their full assumptions are established. |
| N2: validated-reference gating claims | Root README and TreeDecomposition docstring now describe supported empirical plug-in p-values and actual gate behavior. Internal support does not certify the reference law. |
| N3: calibrated passthrough values | The passthrough report and audit now say empirical-inflation plug-in or inflation-adjusted pre-BH p-values. Reported numbers remain unchanged. |
| N4: invariant/equivariant PCA labels | Core docstrings explicitly define the retained identifiers as operational labels with no declared transformation/group action. scRNA help, plot labels and manifest descriptions identify leading and following PCA axes. |
| N5: hidden target normalization | NNLS documentation and output metadata identify mean squared standardized Euclidean distance, sample standard deviation with ddof=1, and division by total input feature count, including constant columns. Prepared squared Euclidean targets have divisor one. |
| N6: branch time caption | The scRNA fit caption describes summed fitted branch length and its normalized target, and states that it is not calibrated stochastic time. |
| N7: incomplete diffusion provenance | Adaptive metadata separately records pydiffmap as kernel backend, the repository coordinate function, the symmetric normalized operator, SciPy eigh and negative-eigenvalue clipping. |

## Metadata compatibility

`BranchLengthOptimizationResult` gains optional trailing fields
`target_definition`, `target_normalization` and
`target_normalization_divisor`. Actual fits populate all three in both
`to_dict()` and `tree.graph["branch_length_optimization"]`. Existing manual
constructors leave them `None`, rather than inventing provenance. Existing
runner serialization forwards the added result fields.

For `squared_standardized_euclidean`, the normalization is
`divide_by_total_feature_count` and the divisor is the number of input columns.
For `squared_euclidean`, normalization is `none` and the divisor is one.
Neither method identifier nor target calculation changes.

Adaptive diffusion retains `backend="pydiffmap"` and adds `kernel_backend`,
`coordinate_implementation`, `coordinate_operator`,
`coordinate_operator_definition`, `coordinate_eigensolver` and
`coordinate_eigenvalue_policy`. Thus existing consumers keep their identifiers
while new exports describe which implementation supplied the coordinates.

## Verification

- The pre-edit baseline passed 34 relevant tests; a separate 13-test check
  covered branch-time constructor compatibility and independent Gaussian NNLS.
- The final focused run passed 71 tests across space separation, NNLS,
  adaptive diffusion, branch-time integration and independent Gaussian tests.
  Existing tests now check normalization metadata against numerical target
  properties and verify that coordinates come from the reported helper.
- The Sol agent compared two NNLS fits against the pre-edit snapshot: four
  leaves, three features including one constant, both target metrics. All
  preexisting result fields except execution time and all fitted edge lengths
  were exactly equal. One seed-42 adaptive diffusion case, 16 by 4 normal
  observations, also had exactly equal coordinates, condensed distances and
  preexisting metadata. These are bounded regression checks, not a benchmark.
- The scRNA diagnostics CLI help renders the corrected descriptions. Changed
  Python files pass Ruff; wiki lint and diff whitespace checks validate the
  documentation changes.

Saved benchmark measurements and their original source hashes were not
regenerated. They describe the source snapshot used for those runs; the newly
added metadata is emitted by subsequent executions. No commit, push,
installation or new calibration benchmark was performed.
