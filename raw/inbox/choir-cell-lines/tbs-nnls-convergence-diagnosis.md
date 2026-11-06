# Full-data TBS NNLS convergence diagnosis

Diagnostic completed on 2026-11-04 using the four saved cancer cell-line
linkages, their common PCA and the installed SciPy 1.16.3. No production
solver setting was changed, and diagnostic edge lengths were not used in TBS.

## Findings

The three alternative fits reproduce their nonconvergence at 1,000 outer
iterations. A controlled change to the accuracy of the inner LSMR solves
removes that iteration-limit termination on all three systems. The native
`success` flag also accepts small relative cost changes when the first-order
criterion is far from its requested tolerance. Separately, the sampled
distance equations do not identify individual edge lengths in these fits.

All rows below use the same frozen problem for each topology: minimize
`0.5 * ||A tau - b||^2` over nonnegative lengths. The sampled pair indices and
targets are identical across topologies; only the path incidence differs.

| Variant | Native stop / iterations | Native cost | Fixed inner tolerance: iterations / cost | Independent cost | Independent first-order residual |
| --- | --- | ---: | --- | ---: | ---: |
| Symmetric / average | Relative cost / 796 | 1856.3662 | — | 535.3519 | 3.960e-6 |
| Symmetric / Ward | Iteration limit / 1,000 | 1570.6670 | 364 / 784.8352 | 545.7035 | 4.530e-5 |
| Markov / average | Iteration limit / 1,000 | 1595.5600 | 388 / 743.5476 | 477.1999 | 4.737e-6 |
| Markov / Ward | Iteration limit / 1,000 | 1494.8301 | 311 / 642.2954 | 482.4050 | 3.674e-6 |

## What causes the iteration-limit failures

Native `lsq_linear` uses trust-region reflective outer iterations and
`lsmr_tol="auto"` for the sparse inner solves. The installed source adapts the
inner tolerance from the current scaled gradient and caps it at 0.1. In the
three failures the inner solves ran at most two LSMR iterations on 994, 992
and 992 of 1,000 outer steps, respectively. The native sparse designs have
between one and approximately 28,000 observed paths per edge; this uneven
support is measured, but its independent effect on solver conditioning was
not isolated.

Holding the initial unconstrained solve, topology, pairs, targets, outer
tolerance and 1,000-iteration budget fixed, replacing only the inner LSMR
tolerance with 0.001 gives the third numeric column above. All three runs
stop on relative cost change and still have scaled first-order residuals
0.3714, 0.3796 and 0.3688. This isolates inner-solve accuracy as a contributor
to the original slow progress; the tolerance change does not establish
first-order convergence. Two exploratory runs with inner tolerance `1e-7`
were interrupted after several minutes of CPU work and had no final fit.

Increasing only the Markov/Ward outer budget to 5,000 preserves its first
1,000 logged iterations and yields status 2 after 1,028 iterations, cost
1473.7089 and scaled optimality 1.4118. The native Symmetric/average baseline
also has status 2 and scaled optimality 2.8945. Both are well above `1e-5`,
the separate first-order criterion. Installed SciPy accepts status 2 for
relative cost change; the TBS wrapper accepts `solution.success` without an
additional optimality check. All 2,540 earlier simulations also stopped on
relative cost change with reported optimality above `1e-5`. Their clustering
results describe the procedure that ran, including its stopping rule.

## What a numerical convergence check should establish

With `r = A tau - b` and `g = A.T @ r`, the NNLS KKT conditions are primal
feasibility `tau >= 0`, dual feasibility `g >= 0`, and complementarity
`tau * g = 0`. Recompute them at the returned lengths in the original units;
the reported gradient can describe the preceding iterate at the iteration
limit. The scaled residual used here takes the maximum absolute component
of `g * where(g > 0, tau, 1)`, matching the nonnegative-bound TRF test.

For an independent check, column-normalize the variables and solve the
identical objective with L-BFGS-B, starting from each native returned vector.
Map the solution back to original length units before evaluating the residual
and cost. To bound objective error, shift `r` upward by
`max(0, max(-g / (A.T @ 1)))` plus a small floating-point margin. The resulting
`y` has numerically verified `A.T @ y >= 0`; its dual objective
`-0.5 * ||y||^2 - b.T @ y` is a lower bound. Each independent fit's relative
primal-dual gap is below `2e-5`, and three of four have the recomputed scaled
residual below `1e-5`. These are numerical checks, not interval arithmetic.
The independent solver's own stop flag is not the evidence of accuracy.

A reliable acceptance rule should use feasibility, scale-appropriate
optimality/complementarity and a checked objective gap, with tolerances
explicitly justified for the intended use. The diagnostic thresholds above
are measured comparisons; no production acceptance threshold was validated.

## Why an accurate objective does not identify edge lengths

Each final design has **56,414 rows and 97,756 columns**, hence at least 41,342
linear null directions. Coverage augmentation removed all zero columns; the
previously logged zero-column count refers to the design before augmentation.
The two root columns are identical, and the native linkage-ratio convention
addresses only that pair.

For a two-leaf subtree whose within-subtree pair was not sampled, adding
length to its incoming edge and subtracting the same amount from both
pendant edges changes no sampled path length. Such missing pairs occur in
16,515–17,625 two-leaf subtrees per topology. For every independent fit, an
explicit perturbation away from the root pair changes three lengths, remains
nonnegative and preserves predictions to within `1.78e-15`. For Markov/Ward
the perturbation magnitude is 1.3411 length units. Thus even near-minimal
sampled cost does not establish stable individual lengths. An identifying
pair design or an explicit rule selecting among equivalent lengths is still
needed, and the downstream effect on TBS gates has not been tested.

## Evidence and verification

- Local run: `results/analyses/choir_nnls_convergence_diagnosis_20261103_185503/`.
  It retains the four sparse designs, original and diagnostic length vectors,
  outer/inner traces, exact perturbations and `verification.json`.
- `tbs-nnls-convergence-evidence.json` and `tbs-nnls-convergence-comparison.csv`
  preserve numerical outputs and hashes. Script snapshots are adjacent.
- Exact replay of original baseline and Markov/Ward cost, iterations and
  reported optimality was verified. The other two original runs saved only
  nonconvergence status; their detailed diagnostics are from these replays.
- Seven focused NNLS tests, tiny dense-reference checks, KKT/dual-bound
  self-checks and explicit feasible null-direction checks passed. Native and
  installed solver hashes stayed fixed. All workers ended.
- [SciPy `lsq_linear` reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.lsq_linear.html)
  explains the status and tolerance distinction. The linked manual is newer
  than 1.16.3; installed source was checked for the executed semantics.
