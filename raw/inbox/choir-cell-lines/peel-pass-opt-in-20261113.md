# The opt-in's current rule with the peel pass on one Gaussian population

Run on 2026-11-13, Europe/Vienna. Local folder: `results/analyses/peel_pass_opt_in_20261113` (git-ignored);
protocol, scripts and outputs are in `peel-pass-opt-in-scripts-20261113.py.txt`. Follow-up of
`peel-pass-default-test-20261113.md`. No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4,
SciPy 1.16.3, uv 0.9.5.

## Setup

The runs of `peel-pass-default-test-20261113.md` with one change: the method parameters carry
`sibling_best_flat_cut_reference=True`, so the sibling test is the opt-in's current rule (k from
`invariant_directions`, `rho` against the table with k up to 15). One Gaussian population, method
`tbs_family_diffusion` with operator symmetric or markov and linkage average or ward. Two arms on the same data
and tree: unchanged code, and the peel pass installed at run time (a node whose smaller child has fewer than 15
cells passes the edge gate like a peel of one cell). Grid `e3` ((200, 80), (200, 30), (256, 30), (496, 30),
seeds 401 to 405) and four grids with 30 features and seeds 1 to 50 (256 or 496 cells, identity covariance or
the unequal variances of the CHOIR root). 880 fits per arm. The protocol with three predictions was written
before the runs.

Self-check passed: the unchanged arm on grid `e3` equals `cancer_pc_systematic_debugging_20261112/e4.json` in 80
of 80 fits.

## Result: fits with more than 1 cluster

| Grid | Fits | Unchanged | Peel pass | Paired under peel pass: more / equal / fewer |
| --- | ---: | ---: | ---: | --- |
| `e3`, seeds 401 to 405 | 80 | 0 | 0 | 0 / 80 / 0 |
| (256, 30) identity, seeds 1 to 50 | 200 | 0 | 0 | 0 / 200 / 0 |
| (496, 30) identity, seeds 1 to 50 | 200 | 0 | 0 | 0 / 200 / 0 |
| (256, 30) unequal, seeds 1 to 50 | 200 | 0 | 0 | 0 / 200 / 0 |
| (496, 30) unequal, seeds 1 to 50 | 200 | 1 | 1 | 0 / 200 / 0 |
| All | 880 | 1 | 1 | 0 / 880 / 0 |

The one fit with more than 1 cluster is the same in both arms: 496 cells, unequal variances, seed 24,
markov/ward, 4 clusters. The found clusters are equal in the two arms in 880 of 880 fits.

The three sibling rules with the peel pass on these 880 fits:

| Sibling rule | More than 1 cluster, unchanged | More than 1 cluster, peel pass |
| --- | ---: | ---: |
| Today's default test (860 fits in 30 features) | 703 | 825 |
| Opt-in, current rule (table, k up to 15) | 1 of 880 | 1 of 880 |
| Opt-in, third version of the draft k | 0 of 880 | 0 of 880 |

## Predictions

- U1 held: identity covariance, 0 of 480 fits with more than 1 cluster in both arms.
- U2 held: unequal variances, 1 of 400 fits under the peel pass.
- U3 held: 0 fits gain clusters and 0 lose clusters.

## Reading

- With the opt-in's current rule the peel pass does not change the result on one Gaussian population in any of
  the 880 fits. The opt-in's current rule holds one Gaussian population at 1 cluster with and without it.
- The over-splitting under the peel pass seen on 13 November belongs to today's default sibling test, not to
  the peel pass.
- The opt-in's current rule is not at 0: 1 of 400 fits with unequal variances is split into 4 clusters, in both
  arms.

## Not checked

- The calls of the reference were not recorded, so the number of tested nodes and their p-values under the
  peel pass are not known for the current rule.
- Limits other than 15 cells. Populations that are not Gaussian. The benchmark cases with groups under the
  current rule with the peel pass (on the CHOIR default tree it gave 432 clusters and ARI 0.6210 on 12
  November). Why seed 24 is split. The peel pass exists at run time only.
