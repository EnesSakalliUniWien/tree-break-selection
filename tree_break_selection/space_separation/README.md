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
The coordinates still use the symmetric normalized operator and clip negative
eigenvalues; this change does not redefine them as Markov right-eigenvector
coordinates.

Adaptive diffusion metadata retains `backend="pydiffmap"` for compatibility
and separately records `kernel_backend`, `coordinate_implementation`,
`coordinate_operator`, `coordinate_operator_definition`,
`coordinate_eigensolver`, and `coordinate_eigenvalue_policy`. Pydiffmap supplies
the kernel; repository code computes the symmetric normalized coordinates
using SciPy `eigh`. The backend label alone does not identify the coordinates.
