# NNLS branch-length solver: switch to column-scaled L-BFGS-B

Work on 2026-11-06, uncommitted, on merge `50455316`.

## Why

The fixed-topology NNLS fit used SciPy `lsq_linear(method="trf",
lsmr_tol="auto")` and accepted any `solution.success`, including stops on
relative cost change. `tbs-nnls-convergence-diagnosis.md` had shown these fits
end far above the optimum. A cold-start timing on the saved Markov/Ward system
(56,414 by 97,756, 2,252,160 non-zeros; 4 threads;
`tbs-nnls-solver-benchmark-script.py.txt`, `tbs-nnls-solver-benchmark-20261106.json`):

| Solver | Seconds | Cost | Scaled optimality | Relative gap |
| --- | ---: | ---: | ---: | ---: |
| `lsq_linear` TRF, production settings | 24.24 (1,000-iteration limit) | 1494.83 | 2.17 | 151.8 |
| L-BFGS-B, column-scaled | 22.04 (3.66 to within 1e-4 of final cost) | 482.41 | 2.1e-6 | 8.4e-6 |
| scikit-learn positive Lasso, coordinate descent | 66.76 (limit) | 482.41 | 2.6e-5 | 1.8e-5 |

scikit-learn `Ridge(positive=True)` (unscaled L-BFGS-B) took 126.6 s at
alpha 1e-6 for cost 482.44. `LinearRegression(positive=True)` rejects sparse
input. JAX, PyTorch, OSQP, Clarabel and cvxpy are not installed and were not
timed.

## Change

`tree_break_selection/tree/optimized_branch_lengths.py` now minimizes
`0.5 ||A tau - b||^2` over `tau >= 0` with SciPy L-BFGS-B on unit-norm columns
from a zero start. A callback checks every 25 iterations and stops once a
certified relative duality gap is at most `solver_tolerance`. The fit's status
is `ok` only when that gap holds for the returned lengths; the gap is recorded
as `relative_duality_gap`. `optimality` is the first-order bound measure
`max |g * where(g > 0, tau, 1)|`.

The certificate uses the dual point `y = r + w` with `w_i` the largest
`max(-g_j, 0) / colsum_j` over the columns in row `i`, which makes
`A.T @ y >= 0` for a nonnegative design without empty columns. A single global
shift (as in the diagnosis) is valid but about ten times looser here: for the
same solutions it gave gaps near 1e-4. On a 400-by-25 test it floors near 4e-8
at a floating-point optimum.

The registered budget `branch_length_optimization_max_iterations` rises from
1,000 to 5,000 in the method registry and both pancreas scRNA configurations,
because L-BFGS-B iterations differ from TRF iterations. The registered
tolerance 1e-5 now means the relative gap.

## Measured effects

On the four saved cancer systems at tolerance 1e-5: 975 to 1,850 iterations,
8 to 18 s, certified gaps 7.9e-6 to 9.1e-6, and costs equal to the diagnosis's
independent optima (535.35, 545.70, 477.20, 482.41).

Full cancer data (`tbs-lbfgsb-*-summary-20261106.json`):

| Method | NNLS cost (before → after) | Clusters | ARI vs cell lines | Lines at 90%/90% |
| --- | --- | --- | --- | --- |
| `tbs_family_diffusion` | 1856.37 → 535.35 | 70 → 75 | 0.0276 → 0.0275 | 27 → 27 |
| `tbs_family_diffusion_markov_ward` (5,000) | 1451.18 → 481.03 | 1,015 → 954 | 0.534 → 0.549 | 45 → 41 |

Frozen count nulls, same datasets, 2,400 fits all certified (300 to 475
iterations, 0.5 to 1.2 s each; `tbs-lbfgsb-null-summary-20261106.csv`); false
partitions out of 300, before → after:

| Variant | Low depth | High depth |
| --- | --- | --- |
| Symmetric / average | 2 → 2 | 264 → 259 |
| Symmetric / Ward | 295 → 295 | 228 → 220 |
| Markov / average | 41 → 40 | 298 → 298 |
| Markov / Ward | 297 → 298 | 300 → 300 |

Accurate branch lengths change neither calibration nor cell-line recovery
materially. The sampled design still leaves at least 41,342 null directions,
so individual edge lengths remain non-identified.

## Checks

New tests cover the certificate against dense `nnls`, a badly scaled solve and
the non-converged policy; the full ordered suite passed 1,045 tests, and Ruff
and the dependency/dead-code audit passed.
