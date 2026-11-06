# Space Separation

This package owns reusable methodological space separation. It has no
dataset-specific file paths, labels, plotting, or report writing.

- `decompose_leading_axes_space(...)` separates the leading PCA axis from the
  following orthogonal PCA axes. No group action or transformation property is
  declared.
- `separate_adaptive_cosine_space(...)` weights a feature matrix, builds its
  sample cosine eigensystem, and segments the ordered spectrum into adaptive
  coordinate blocks.
- `hamming_knn_diffusion_geometry(...)`, `adaptive_diffusion_geometry(...)`,
  `block_diffusion_geometry(...)`, and `block_adaptive_diffusion_geometry(...)`
  return one explicit result containing coordinates, condensed distances, and
  construction evidence without importing benchmark runners.

Dataset adapters belong under `applications/`; benchmark evaluation and null
diagnostics belong under `benchmarks/`.

Cosine and diffusion rank arguments are targets. Rank selection extends through
the last complete numerical eigenvalue group, using the shared rule in
`tree_break_selection/eigenvalue_groups.py`: adjacent values within relative
tolerance `1e-10` and absolute tolerance `1e-12` times the spectral scale join
transitively. Returned dimensions can therefore exceed the requested target;
diffusion metadata records the actual dimension.

Diffusion removes its leading mode only when that eigenvalue is isolated. A
repeated leading group is retained whole, including on disconnected graphs.
This preserves the group rather than choosing an arbitrary direction to remove.
The coordinates use the symmetric normalized operator by default and clip
negative eigenvalues. `coordinate_operator="markov"` (on
`compute_diffusion_coordinates`, `adaptive_diffusion_geometry` and
`feature_family_diffusion_geometry`) instead returns `D^(-1/2)` times them:
the right eigenvectors of the random walk `D^(-1) W`, without the per-row
root-degree factor. Metadata then records
`coordinate_operator="markov_right_eigenvectors"`. The opt-in benchmark method
`tbs_family_diffusion_markov_ward` uses these coordinates with Ward linkage.

Adaptive diffusion metadata records `backend="pydiffmap"`, `coordinate_implementation`,
`coordinate_operator`, `coordinate_operator_definition`,
`coordinate_eigensolver`, and `coordinate_eigenvalue_policy`. Pydiffmap supplies
the kernel; repository code computes the symmetric normalized coordinates
using SciPy eigensolvers. The backend label alone does not identify the coordinates.

Large sparse operators use SciPy `eigsh`. Disconnected graphs are solved by
connected component, using `eigh` only for small blocks or a block requiring its
full spectrum. The component spectra are merged before the global rank and
complete-tie rule is applied. This avoids a full dense graph allocation while
preserving repeated leading modes. Solver metadata identifies the mixed
`eigsh`/`eigh` component path. Recovery of repeated nontrivial eigenvalues inside
a large connected block retains the existing ARPACK limitation documented in code.

Adaptive diffusion sets the scaled-distance identity diagonal to zero before
epsilon selection and kernel computation, so positive floating-point self
distances cannot enter epsilon's positive-distance sample. Off-diagonal entries,
including stored zero distances between duplicate observations, are preserved.
Metadata records
`self_distance_policy="zero_identity_diagonal_before_epsilon_and_kernel"`;
feature-family and block adaptive diffusion retain this policy.
