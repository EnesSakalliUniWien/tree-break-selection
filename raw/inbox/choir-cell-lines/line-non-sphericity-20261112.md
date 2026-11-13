# What makes the cells of one cell line non-spherical

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/line_non_sphericity_20261112` (git-ignored);
the protocol, the deviations of both parts, the scripts, their outputs and two independent spot checks are in
`line-non-sphericity-scripts-20261112.py.txt`. Follow-up of `k-margin-choir-20261112.md`. No production code was
changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Question

On the CHOIR cells (48,879 cells, 30 PCs, 190 cell lines) the largest eigenvalue inside one cell line is 7.0 to
94.1 times the median eigenvalue (175 lines with at least 100 cells; median 25.6), against a 99% point of 1.30 to
2.93 for one spherical Gaussian population. Where does this come from, and can it be removed without the
cell-line labels?

The protocol with definitions, self-checks and six predictions was written before any run. The work was split
in two parts, built separately and checked before use: Part A (the lines, with labels) and Part B (without
labels, from the Markov/Ward tree). T1 is the largest eigenvalue over the median eigenvalue of a set of cells;
W is the pooled within-line covariance over all 190 lines.

## Part A: the lines, with labels

Self-checks passed: the k of the 175 lines reproduces the earlier run under both rules; one Gaussian population
with covariance W, whitened by W, gives a direction in 33 of 2,000 draws (allowed range 8 to 36).

| Step | T1 over the 175 lines, median (min, max) | k current rule | k prototype rule |
| --- | --- | --- | --- |
| PC coordinates | 25.6 (7.0, 94.1) | 5 to 12 | 6 to 27 |
| Whitened by W without the line | 13.8 (3.7, 117.1) | 4 to 10 | 6 to 28 |
| Depth fit removed | 25.4 (7.0, 85.4) | 5 to 12 | 7 to 27 |
| Depth fit removed, then whitened | 13.7 (3.7, 108.6) | 4 to 11 | 6 to 28 |

- Profile in PC coordinates: the median within-line variance per PC is 0.72 to 3.43 (largest at PC18, 3.43,
  which is 45% of that PC's total variance; PC1 1.62, 3% of its total variance). The share of off-diagonal
  energy in a line's covariance is 0.658 at the median (0.435 to 0.846).
- Sequencing depth (recovered from the CPM matrix as on 7 November; log depth median 9.694, sd 0.549): a linear
  fit of every PC on log depth explains 3.0% of the within-line variance at the median (0.7% to 8.6%); the first
  principal score of a line correlates with log depth at 0.154 in absolute value at the median (0.002 to 0.580).
- Pools: no line has cells from more than one `pool_id`; all 190 lines have cells from more than one
  `source_library`, whose means explain 0.8% of the within-line variance at the median (0.1% to 11.4%).
- After whitening by W without the line, T1 exceeds the spherical 99% point in 175 of 175 lines (T1 over the
  99% point: median 6.7, 1.7 to 63.9), and the prototype rule gives at least 2 directions in 175 of 175 lines.

## Part B: without labels, from the tree

Self-checks passed: the root holds all 48,879 cells; both split files list the same 1,838 called nodes; the
recorded k of 20 random called nodes is reproduced under both rules; W has trace 47.240 and largest eigenvalue
6.274 in both parts.

W_tree(m) is the pooled covariance of the maximal subtrees with at most m cells and at least 5:

| m | Subtrees | Cells covered | Subtrees with at least 90% from one line | Eigenvalues of W^-1/2 W_tree W^-1/2: smallest, median, largest |
| ---: | ---: | ---: | --- | --- |
| 30 | 2,610 | 48,270 | 2,274 (87.1%) | 0.243, 0.781, 0.985 |
| 50 | 1,663 | 48,499 | 1,402 (84.3%) | 0.260, 0.808, 1.027 |
| 100 | 928 | 48,675 | 747 (80.5%) | 0.284, 0.851, 1.082 |

The 175 lines whitened by W_tree(m): T1 median 14.6, 14.5 and 14.3 for m = 30, 50, 100 (25.6 in PC coordinates);
k current rule 4 to 10, k prototype rule 5 to 28.

The 1,838 called nodes whitened by W_tree(50), k above 15 under the prototype rule: 290 nodes (417 in PC
coordinates); at the 103 nodes with at least 1,000 cells k is 18 to 28 (18 to 28 in PC coordinates as well).
Median k by size class, PC coordinates to whitened: current rule 6 to 4 (32 to 99 cells), 9 to 7 (100 to 999),
14 to 13 (1,000 and more); prototype rule 10 to 7, 14 to 12, 27 to 28.

Called nodes with at least 40 cells, by composition:

| | Nodes | With k prototype >= 1 | Median k prototype, PC coordinates / whitened | Median k current, PC coordinates / whitened |
| --- | ---: | --- | --- | --- |
| At least 90% of the cells from one line | 1,017 | 1,017 in both | 12 / 9 | 7 / 6 |
| All other nodes | 342 | 342 in both | 18 / 17 | 10 / 8 |

## Predictions

| Prediction | Result |
| --- | --- |
| PA1 after whitening by W without the line the median T1 is below 5 | not held (13.8) |
| PA2 after that, k prototype >= 1 in more than half of the lines | held (175 of 175) |
| PA3 after removing the depth fit the median T1 stays above 10 | held (25.4) |
| PB1 W_tree(50) agrees with W within a factor 2 in every direction | not held (smallest eigenvalue 0.260) |
| PB2 whitened by W_tree(50), k prototype is at most 15 at all 103 nodes with at least 1,000 cells | not held (0 of 103) |
| PB3 whitened, nodes inside one line have k prototype >= 1 less often than the other nodes | not held (1,017 of 1,017 against 342 of 342) |

## Independent spot checks (own code, `verify_a.py`, `verify_b.py`)

Whitened by W without the line: T1 median 13.774 (3.716 to 117.054). Log depth median 9.694, sd 0.549; depth
share 0.030 (0.007 to 0.086). W_tree(50): 1,663 subtrees, 48,499 cells, eigenvalues 0.260, 0.808, 1.027; lines
whitened by it: T1 median 14.5 (3.89 to 112). All equal the values of Parts A and B. The direction in which the
small subtrees hold least of the within-line variance loads most on PC2, PC10, PC9, PC3 and PC1.

## Reading

- The non-sphericity is specific to each cell line. A covariance common to all lines removes about half of the
  ratio (25.6 to 13.8) and leaves every line far above the spherical 99% point.
- Sequencing depth (3.0% of the within-line variance) and pools (0.8%) are not the cause.
- A common covariance can be estimated without labels from small subtrees and is close to the label-based one in
  most directions (median ratio 0.81), but whitening by it changes the picture as little as the label-based
  one. In one direction the small subtrees hold only a quarter of the within-line variance (0.26).
- A count of eigenvalues above a spherical edge does not tell nodes inside one line from nodes with several
  lines, in PC coordinates or after whitening: every such node has a direction. The count is higher at nodes
  with several lines (median 17 against 9 whitened), but no threshold was examined.
- For step A2 of the plan: a common whitening does not rescue a definition of k by a spherical edge. Each line
  has its own covariance.

## Deviations

Part A: depth recomputed from the CPM matrix because no per-cell depth file existed; the depth fit is ordinary
least squares of all 30 PCs on log depth within each line; no line spans more than one `pool_id`, so the pool
means are the `source_library` means. Part B: a check that both split files list the same nodes was added; "in
PC coordinates" at the called nodes uses the recorded k; W_tree(30) and W_tree(100) were not applied to the
called nodes. No assert or prediction was changed in either part.

## Not checked

Which directions carry the structure inside a line (cell cycle or other gene programs). Whether a
line-specific covariance could be estimated without labels. A threshold on the count of directions. Only the
Markov/Ward tree, one dataset, lines with at least 100 cells.
