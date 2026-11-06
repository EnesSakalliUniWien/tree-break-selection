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

## Conditional root decisions

`selected/hierarchy/diffusion_root_conditional_mc.py` provides the opt-in
`selected_diffusion_root_chain_test` API for one native diffusion/average root:

```python
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_conditional_mc import (
    selected_diffusion_root_chain_test,
)

result = selected_diffusion_root_chain_test(
    values, sigma=known_noise_standard_deviation, steps=999, seed=17, alpha=0.05
)
```

Supply independent Gaussian rows with known common covariance `sigma**2 * I`
in the actual input features. The function does not estimate that scale or
validate fitted-PC, transformed-count or correlated-row inputs. It preserves
the observed root children and uses a random-position reversible chi chain;
other merges may differ in the simulated reference data.

Read `decision` (`reject`, `do_not_reject`, `unresolved`) together with
`status`, `sampling_status`, `accepted_moves` and `resolution_sufficient`.
`do_not_reject` does not establish equal means. A stalled chain retains its
conservative rank p-value of one but returns an unresolved decision. A step
budget unable to reach `alpha` also returns an unresolved decision. Failed
geometry gives no p-value. `chain_radii`, `anchor_index`, `root_children` and
`source_hashes` support replay and auditing.

Set the event, budget and seed before evaluating the result; choosing a seed
or rerunning until rejection invalidates the single-test interpretation.
Run this diagnostic serially because native geometry replay temporarily
instruments module functions. Numerical root membership is not an exact-real
certificate. This API does not replace production TBS gates or provide a
descendant/multiple-testing rule. Validation and the fixed protocol are in
`reports/root_conditional_sampler_20261104/`.
