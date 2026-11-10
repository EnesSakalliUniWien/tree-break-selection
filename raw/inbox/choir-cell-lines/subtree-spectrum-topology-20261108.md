# Splits, leaves, eigenvectors, eigenvalues, geometric multiplicity and separations outside the topology

Analysis on 2026-11-08. Local folder: `results/analyses/subtree_spectrum_topology_20261108` (git-ignored);
protocol, scripts, calibration table and summary output in `subtree-spectrum-topology-scripts-20261108.py.txt`.
Descriptive; no decision rule. No deviations from the protocol.

## Definitions

Per subtree, in the 30 PCs the tests use: n leaves; 2^(n-1) - 1 two-group splits (n - 1 bits); (2n - 3)!! rooted
binary topologies (log10: 985 for 400 cells, 2,864 for 1,000, 222,681 for 48,879). Eigenvalues of the subtree's
cell covariance; noise edge from the eigenvalue bulk (median matched to the Marchenko-Pastur median); k =
eigenvectors above the edge (invariant directions), p - k below (fuzzy). Geometric multiplicity: consecutive
eigenvalues above the edge stay in one group while a block likelihood-ratio test of equality,
(n - 1)[m log(mean) - sum log l_i] against chi2((m - 1)(m + 2)/2), does not reject at 0.01 (derived here;
attribution to Anderson 1963 or Lawley 1956 not verified). A separation fits the rooted topology below a
subtree when one side is a clade there (95% coverage and purity). A rooted binary tree over g groups holds
exactly 2g - 3 of the 2^(g-1) - 1 group bipartitions (3 of 3 for g = 3, 5 of 7 for g = 4, 7 of 15 for g = 5).

## E2 Calibration of the multiplicity rule (500 runs each, p = 30, unit bulk)

| Truth | n = 100 | n = 400 | n = 1,000 |
| --- | ---: | ---: | ---: |
| (20, 20, 20): one eigenspace, multiplicity 3 | 98.6% | 99.2% | 99.0% |
| (20): one spike | 94.4% | 93.2% | 95.8% |
| (30, 20, 10): three distinct | 5.8% | 85.2% | 99.8% |

Multiplicity is a statement at a given number of cells. Failure modes for (30, 20, 10), same seed: at 100
cells the top two merge into one eigenspace of multiplicity 2 in 78% of runs ((2, 1) 391/500, (1, 2) 62,
(3) 13, resolved 33); at 400 cells (2, 1) remains in 63/500.

## E1 Invariant directions and multiplicity (simulated; 50 replicates per no-group model, 35 group datasets)

- No-group roots: k = 0 in all 50 at log-depth sd 0.25; k = 1 in 60% at sd 0.42 and 100% at sd 1.0 (the
  sequencing-depth direction stands above the noise band).
- Group roots: balanced 4 groups k = 3 with multiplicity 3 (10/10, a degenerate eigenspace); imbalanced and
  singleton k = 3 with 1+1+1 (10/10); density variation k = 4 with 1+1+1+1 (5/5).
- One-group subtrees chosen by the tree: Markov/Ward subtrees of no-group data show a median 3 (sd 0.25) to 7
  (sd 0.42) directions above their noise band; default-tree subtrees mostly 0 (86% at sd 0.25). A subtree cut
  out by the tree is not a random sample of one population, so its spectrum is already shaped by the tree's
  choice.

## E3 Separations the tree cannot hold (simulated nodes with 4 groups)

| Model | Multiplicity at node | Bipartitions | Held by the tree (median) | Dropped, strength relative to best held (median) |
| --- | --- | ---: | ---: | --- |
| Balanced (250 x 4) | 3 | 7 | 5 | the two other 2+2 splits, 0.95-0.97 |
| Imbalanced (700/200/80/20) | 1+1+1 | 7 | 6 | {700, 80} vs {200, 20}, 0.92 |
| Singleton | 1+1+1 | 7 | 6 | {700, 80} vs {200, 20}, 0.92 |
| Density variation | 1+1+1+1 | 7 | 5-6 | group of 600 cells alone, 1.005 (5/5); {2100, 240} vs rest 0.97 |

(6 held for g = 4 exceeds 2g - 3 = 5 because the 95% tolerance lets a clade without the 20-cell group stand
for a side that contains it.) In the balanced data every 1-vs-3 separation is held and one of the three 2+2
separations; the other two are almost as strong and cannot be in the same tree. With distinct eigenvalues the
dropped separations are near copies of a held one with a small group moved, 8% weaker. In the density-variation
data the strongest separation of all (the 600-cell group alone) is not a clade in either tree.

## E4 Real CHOIR data (48,879 cells, 190 lines; nodes >= 500 cells whose smaller child has >= 50)

| Tree | Nodes | Median k | Median lines present | Nodes cutting >= 1 line | Lines cut (sum) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Default (symmetric/average) | 172 | 13 | 13 | 49% | 423 |
| Markov/Ward | 171 | 12 | 6 | 27% | 70 |

Markov/Ward root: 14 eigenvectors above the noise band (top eigenvalues 53.1, 35.7, 24.8, 20.4, 18.2, 16.2;
grouping 1+1+1+1+1+1+1+1+3+2+1), 190 lines, splits off 267 cells, cuts no line. The default tree's root is a
peel and is outside the analysed set. Most analysed nodes have some eigenvalue group of multiplicity >= 2
(162/172 default, 149/171 Markov/Ward). Nodes with one cut 51% (default) and 27% (Markov/Ward) against 30%
and 27% for nodes without; with 10 and 22 such nodes the comparison is weak.

## Reading

The number of splits is astronomical (2^(n-1) - 1) at every size, but the data support only the separations
between groups, and a tree can hold just 2g - 3 of their 2^(g-1) - 1 bipartitions. When the top eigenvalues
form one eigenspace (balanced groups, multiplicity 3), several separations are about equally strong and the
tree keeps one of them arbitrarily; the others are alternative clusterings it cannot show. When the eigenvalues
are distinct the hierarchy is close to determined and what is dropped is weaker. Whether eigenvalues count as
equal depends on the number of cells. Inside subtrees the tree has already chosen, invariant directions appear
even in pure noise, so this reading holds at roots or with known groups, not as a test inside selected
subtrees.
