# Distances in the different subspaces and cells equidistant in two spaces

Analysis on 2026-11-08. Local folder: `results/analyses/subspace_distance_ties_20261108` (git-ignored); protocol,
code and summary output in `subspace-distance-ties-scripts-20261108.py.txt`. Descriptive; no deviations.

## Spaces and measure

For a split of subtree v into children A and B: tree (the coordinates the tree was built in: symmetric diffusion
for the default tree, Markov diffusion for Markov/Ward), pcs (30 test PCs), test (PCs whitened by the pooled
within-children covariance, the sibling test's metric), invariant and fuzzy (PCs projected on the eigenvectors
of v's covariance above and below its noise band). Margin of a cell m = (distance to the other child's centroid
- distance to its own) / distance between the centroids: 1 at its own centroid, 0 equidistant, negative on the
wrong side. Tie: |m| < 0.1 (stated convention). Separation: squared centroid distance / mean squared distance to
the own centroid. Split types from true labels: between groups, within one group, mixed. Real trees: their saved
coordinates reproduce the first leaf-leaf merge heights (49/49, 16/16).

Sanity check: two groups 12 SD apart give 0 ties and 0 wrong-side cells; one population cut at its mean gives
37% ties and 8% wrong-side cells in 30-dimensional Euclidean distance.

## Simulated (185 datasets; median over splits >= 60 cells)

| Split type | Tied cells: tree / test / invariant / fuzzy | Wrong side: tree / test | Separation: tree / test | Tied in tree and test | Side differs tree vs test | Side differs invariant vs fuzzy |
| --- | --- | --- | --- | ---: | ---: | ---: |
| Between groups and mixed (group data) | 0 / 0 / 0 / 34-40% | 0 / 0 | 2-29 / 25-79 | 0 | 0 | 40-51% |
| Within one group (no-group data) | 7-11% / 11-16% / 2-11% / 20-29% | 4-6% / 5-11% | 0.52-1.37 / 0.19-0.34 | 1.6-4.0% | 7-12% | 36-51% |
| Within one group (group data) | 7-10% / 12-21% / 3-10% / 25-30% | 4-6% / 4-11% | 0.66-1.26 / 0.19-0.37 | 2.0-3.6% | 7-12% | 39-45% |

Exception: the mixed splits of the density-variation data on the default tree (29 splits) have test separation
0.76, 1% tied cells in the tree and test spaces and 5.4% wrong-side cells in the test space. The default tree's
within-group splits have more tied cells in the test space (14-21%) than Markov/Ward's (11-13%).

## Real CHOIR data (splits >= 500 cells with a smaller child >= 50; medians)

| Tree | Splits | Tied: tree / test / invariant / fuzzy | Wrong side: tree / test | Separation: tree / test | Side differs tree vs test |
| --- | ---: | --- | --- | --- | ---: |
| Markov/Ward | 171 | 0.6% / 1.3% / 5.7% / 26.3% | 0.4% / 0.7% | 7.8 / 1.2 | 1.2% |
| Default | 172 | 3.5% / 3.5% / 7.4% / 25.9% | 1.4% / 2.1% | 4.3 / 0.7 | 3.7% |

Splits that cut at least one line (each child holds >= 10% of the line's cells in the subtree): 47 of 171 cut 70
lines (Markov/Ward), 85 of 172 cut 423 (default). In those splits, cells of the cut lines are tied in the test
space at a median 22.8% (Markov/Ward) and 39.6% (default), other cells at 7.6% and 6.7%.

## Reading

Real separations leave no cell on the fence in the tree, test or invariant space; the fuzzy directions carry no
side information at all (the side there agrees with the invariant side about half the time). Splits inside one
population put 7-21% of cells on the fence and 7-12% of cells on different sides in the tree's space and the
test space: the tree decides them in its own space, where its split looks stronger (separation 0.5-1.4) than in
the test space (0.2-0.4). In the real data, the cell lines a split cuts are cut where their cells are on the
fence in the test space (23-40% tied against 7% for other cells), so those cuts are forced choices rather than
supported ones.
