# The prototype rule for k with the opt-in reference on the CHOIR cells

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/k_margin_choir_20261112` (git-ignored);
protocol, deviations, the three scripts and their outputs are in `k-margin-choir-scripts-20261112.py.txt`.
Follow-up of `k-margin-prototype-20261112.md`. Working tree on top of `bc044e07` with the uncommitted opt-in
reference. No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Setup

The harness of `results/analyses/spr_refit_20261111` on the unchanged Markov/Ward tree (48,879 cells, 30 PCs,
`tbs_family_diffusion_markov_ward`, sibling alpha 0.01, edge alpha 0.001, opt-in
`sibling_best_flat_cut_reference`). One variable: `invariant_directions` inside `best_flat_cut_reference` is
replaced at run time by the prototype rule `directions_with_margin` (alpha 0.01). The baseline with the current
rule on the same tree was recorded on 11 November (`k15_base`). The protocol with three predictions was written
before the run. The rule was called 1,838 times; the run took 320 s (current rule 126 s).

## Result

| | Current rule | Prototype rule |
| --- | ---: | ---: |
| Found clusters | 230 | 64 |
| ARI against the cell lines | 0.8219 | 0.1285 |
| NMI | 0.9358 | 0.7336 |
| Lines recovered | 116 | 36 |
| Lines split only / merged only / both | 21 / 38 / 15 | 2 / 143 / 9 |
| Largest cluster | 1,990 | 12,884 |
| Calls with p <= 0.01 | 149 | 26 |
| Accepted by BH | 134 | 9 |

ARI between the two assignments: 0.1544.

Call reasons (1,838 calls, the same nodes under both rules):

| | Current rule | Prototype rule |
| --- | ---: | ---: |
| `tested` | 1,475 | 1,268 |
| `no_invariant_directions` | 363 | 2 |
| `too_many_invariant_directions` (k above 15) | 0 | 417 |
| `too_few_cells` (k above 10 with fewer than 50 cells, or fewer than 15 cells) | 0 | 151 |

By size of the node:

| Cells | Calls | Current rule: tested, accepted | Prototype rule: tested, accepted |
| --- | ---: | --- | --- |
| fewer than 15 | 78 | 0, 0 | 0, 0 |
| 15 to 31 | 285 | 0, 0 | 267, 0 |
| 32 to 99 | 606 | 606, 8 | 517, 1 |
| 100 to 999 | 766 | 766, 76 | 484, 8 |
| 1,000 and more | 103 | 103, 50 | 0, 0 |

k is higher under the prototype rule at 1,815 of the 1,838 nodes and equal at 23; at nodes with at least 32
cells the difference is 0 to 19 directions. k of the calls: 3 to 15 (current rule, tested calls) and 1 to 28
(prototype rule). 125 splits accepted under the current rule are not accepted under the prototype rule; 9 are
accepted under both; none is accepted only under the prototype rule.

Predictions: P1 held (267 of the 285 calls with 15 to 31 cells are tested; 2 have p <= 0.01 and none is accepted
by BH). P2 not held (64 clusters against 230). P3 not held (ARI 0.1285 against 0.8219).

## Cause (`cause.py`, added after the result)

- Root, 48,879 cells: the eigenvalues of the 30 PCs fall from 53.09 to 5.12. Eigenvalue j over the median of
  the remaining ones is 6.38 at step 1 and still 1.04 at step 27; the 99% point for one spherical Gaussian
  population of this size is 1.055 to 1.023. 27 steps pass. The current rule gives k = 14 at the root.
- Single cell lines with at least 100 cells (175 lines, median 238 cells), the line's own cells:
  - current rule: k from 5 to 12 (k = 8, 9 or 10 in 139 lines); no line has k = 0;
  - prototype rule: k from 6 to 27; above 10 in 166 lines; above 15 in 71 lines;
  - the largest eigenvalue is 7.0 to 94.1 times the median eigenvalue within a line (median 25.6); the 99% point
    for one spherical Gaussian population at these sizes is 1.30 to 2.93.

## Reading

- The prototype rule fails on the real data. With k above 15 the reference refuses the split, and that happens
  at all 103 nodes with at least 1,000 cells, where the current rule accepted 50 splits.
- The cause is the null, not the margin or the Gram route: the cells of one cell line are far from spherical,
  so a count of eigenvalues above a spherical edge is large inside one population under both rules (5 to 12
  with the current rule). The margin makes the rule more sensitive, and with many cells the 99% point is close
  to 1, so nearly every PC counts.
- The synthetic benchmark hid this because its noise is spherical by construction.
- For step A2 of the plan this means: k cannot be defined as "eigenvalues above the edge of one spherical
  population". On real cells that count does not separate one population from several. The current rule's
  values (3 to 15) come from the same kind of count and are not derived either; they only happen to stay
  within the reference's limit of 15.

## Not checked

Only the Markov/Ward tree and one dataset. Whether the eigenvalue profile inside one cell line comes from the
30 PCs of the whole dataset, from sequencing depth or from biology was not examined. No alternative definition
of k was tried. The reference with a limit above 15 directions was not tried.
