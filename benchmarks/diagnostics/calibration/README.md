# Calibration Diagnostics

This package contains investigation and validation tools for calibration
behavior. It is not a production calibration layer. Production tree
construction, testing, and plotting interfaces remain under
`tree_break_selection/`.

## Categories

- `selected/hierarchy/`: hierarchy selection null-audit and stratification studies.
- `sibling/gates/`: data-independent sibling-gate profiles.
- `sibling/known_branch_skew_diagnostic.py`: known-tree Gaussian ablation of
  covariance, same-data PCA/MP projection, hierarchy selection, and traversal.
- `statistics/`: covariance and test-statistic studies.
- `traversal/`: production-admissibility checks.

The matching tests live under `tests/validation/calibration/` with the same
category path. Import modules from their owning category; the package roots do
not re-export moved names.

Shared diagnostic contracts remain at this package root:

- `reporting.py` writes named table bundles and their manifest envelope.

## Maintained runners

- `sibling/known_branch_skew_diagnostic.py`
