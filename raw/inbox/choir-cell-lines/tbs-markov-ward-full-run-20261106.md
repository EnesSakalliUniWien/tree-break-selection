# Markov/Ward TBS on all cancer cells, and the NNLS solver

Runs on 2026-11-06 of the registered opt-in method
`tbs_family_diffusion_markov_ward` (merge `50455316` plus uncommitted working-tree
changes adding Markov coordinates and the method).

## Setup

`tbs-markov-ward-run-script.py.txt` reuses the 2026-11-03 kernel and 30-PC
input. Markov coordinates come from the production
`compute_diffusion_coordinates(..., coordinate_operator="markov")`; they match
the 2026-11-03 comparison's Markov coordinates to maximum absolute difference
6.94e-18. All TBS settings come from the method registry: Ward linkage,
fixed-topology NNLS branch lengths, edge alpha 0.001, sibling alpha 0.01,
pass-through on, energy fraction 0.90. The script records a failure instead of
retrying with other settings and exits with `os._exit(0)` after writing outputs.

## Registered 1,000-iteration budget

The NNLS fit stopped at the iteration limit, so the runner refused to apply
unrefit branch lengths and ran no gates (`ValueError`,
`status='solver_not_converged'`), as on 2026-11-03. Elapsed 111 s.
Summary: `tbs-markov-ward-1000-summary-20261106.json`.

## One-off 5,000-iteration budget

The only change is `branch_length_optimization_max_iterations = 5000` (two
added lines in the script; recorded in the summary's `params`). The fit
stopped after 1,017 iterations on relative cost change (cost 1451.18, reported
optimality 1.19) and was applied. Elapsed 171 s.

| Quantity | Symmetric / average (merged) | Markov / Ward, 5,000 | Published `CHOIR_1` |
| --- | ---: | ---: | ---: |
| Clusters | 70 | 1,015 | 367 |
| ARI against cell lines | 0.0276 | 0.5337 | 0.6640 |
| NMI against cell lines | 0.5482 | 0.8602 | — |
| Lines recovered at 90%/90% | 27 | 45 | 79 |
| Cells in clusters at least 90% one line | — | 83.7% | 99.3% |
| Median clusters covering 90% of a line | — | 4 | 2 |

- Traversal visited 2,051 nodes: 596 splits, 429 pass-throughs, 1,015 boundaries.
- Cluster sizes: median 27 cells, 27 singletons, 95 clusters of at most 5 cells,
  largest 1,451.
- 811 of 1,015 clusters are at least 90% one line; 130 of 190 lines have at
  least 90% of their cells in their own pure clusters. 52 lines are in one
  cluster; 87 need five or more to cover 90% of their cells.
- ARI against `CHOIR_1` is 0.4115.

The tree separates cell lines far better than the symmetric/average tree but
fragments them, consistent with its 297/300 and 300/300 false partitions on the
frozen count nulls (`tbs-recovery-comparison.md`).

## NNLS solver

Production fits branch lengths with SciPy `lsq_linear(method="trf",
lsmr_tol="auto")` (`tree_break_selection/tree/optimized_branch_lengths.py`) and
treats `solution.success` as convergence, which includes SciPy status 2
(relative cost change). `tbs-nnls-convergence-diagnosis.md` showed for the same
sampled systems:

- The automatic inner LSMR tolerance (capped at 0.1) ran at most two LSMR
  iterations on about 99% of outer steps, slowing progress.
- An independent L-BFGS-B fit (also SciPy) of the identical objective reached
  primal-dual gaps below 2e-5 with much lower cost: 535.35 against 1856.37 for
  the accepted symmetric/average fit, and 482.41 against 1494.83 at 1,000
  iterations for Markov/Ward. The 5,000-budget Markov/Ward fit (1451.18) is
  still about three times that reference cost.
- Each design has 56,414 sampled paths and 97,756 edges, so at least 41,342
  null directions: even an exact optimum does not identify individual edge
  lengths.

The library is not the root problem: the solver configuration and acceptance
rule stop far from the optimum, and the sampled design does not identify the
lengths. Gate decisions on these runs use those approximate, non-identified
branch lengths.

## Files

- Local-only outputs:
  `results/analyses/choir_cell_lines_tbs_markov_ward_20261106/` and
  `results/analyses/choir_cell_lines_tbs_markov_ward_5000_20261106/`.
- SHA-256: Markov coordinates
  fb5aa8c30e1008e7945a1b3f9a0019a6ea2a10f1c10fcdf432d2f5f14cfcdfd2;
  5,000-run `cluster_assignments.csv`
  ef92f067f5d651efb6a29b3774aaac1ef7af125fd677c9cef4986554f7724bc8,
  `linkage.npy` e4c7c7e259190a7a67afd9f091914ffffa9647c1657e92181c51ea22051d6a7b,
  `annotations.csv` 208cded462cb8c6b31ba7b2c49ac2357b90ddab2ce69449e40ae7758aa8208fc.
