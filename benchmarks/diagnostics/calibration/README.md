# Calibration Diagnostics

This package contains investigation and validation tools for calibration
behavior. It is not a production calibration layer. Production tree
construction, testing, and plotting interfaces remain under
`tree_break_selection/`.

## Categories

- `edge/`: edge-null behavior and selected-edge post-run analysis.
- `selected/hierarchy/`: hierarchy selection null-audit and stratification studies.
- `sibling/gates/`: data-independent sibling-gate profiles.
- `sibling/known_branch_skew_diagnostic.py`: known-tree Gaussian ablation of
  covariance, same-data PCA/MP projection, hierarchy selection, and traversal.
- `sibling/nulls/`: sibling null-weight rule validation.
- `statistics/`: covariance, test-statistic, and internal-calibration studies.
- `traversal/`: traversal CLI and production-admissibility checks.

The matching tests live under `tests/validation/calibration/` with the same
category path. Import modules from their owning category; the package roots do
not re-export moved names.

Shared diagnostic contracts remain at this package root:

- `reporting.py` writes named table bundles and their manifest envelope.

## Maintained runners

- `sibling/independent_projection_control.py`: matched independent basis/rank
  controls with pre-run source freezing, known covariance and fixed branches.
- `sibling/known_branch_skew_diagnostic.py`
