# Today's default sibling test with the peel pass on one Gaussian population

Run on 2026-11-13, Europe/Vienna. Local folder: `results/analyses/peel_pass_default_test_20261113`
(git-ignored); protocol, scripts and outputs are in `peel-pass-default-test-scripts-20261113.py.txt`. Follow-up
of `peel-pass-20261112.md`. No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy
1.16.3, uv 0.9.5.

## Setup

One Gaussian population (registered case `gauss_null_large_continuous` with other sizes and seeds), method
`tbs_family_diffusion` with operator symmetric or markov and linkage average or ward, registered parameters, no
opt-in: today's default sibling test. Two arms on the same data and tree: unchanged code, and the peel pass of
12 November installed at run time (a node whose smaller child has fewer than 15 cells passes the edge gate like
a peel of one cell). Grid `e3`: (200, 80), (200, 30), (256, 30), (496, 30), seeds 401 to 405. Four grids with 30
features and seeds 1 to 50: 256 or 496 cells, identity covariance or the unequal variances of the CHOIR root
with one fixed rotation. 880 fits per arm. The protocol with four predictions was written before the runs.

Self-check passed: the unchanged arm on grid `e3` equals `cancer_pc_systematic_debugging_20261112/e3.json` in 80
of 80 fits.

## Result: fits with more than 1 cluster, and mean found clusters

| Grid | Fits | More than 1 cluster, unchanged | More than 1 cluster, peel pass | Mean found, unchanged | Mean found, peel pass | Paired under peel pass: more / equal / fewer |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| (200, 80), seeds 401 to 405 | 20 | 3 | 7 | 1.20 | 1.55 | 4 / 16 / 0 |
| (200, 30), seeds 401 to 405 | 20 | 16 | 19 | 3.95 | 4.40 | 5 / 14 / 1 |
| (256, 30), seeds 401 to 405 | 20 | 17 | 20 | 5.25 | 7.20 | 9 / 10 / 1 |
| (496, 30), seeds 401 to 405 | 20 | 13 | 20 | 7.50 | 12.60 | 13 / 6 / 1 |
| (256, 30) identity, seeds 1 to 50 | 200 | 151 | 195 | 4.75 | 6.49 | 83 / 108 / 9 |
| (496, 30) identity, seeds 1 to 50 | 200 | 115 | 172 | 7.23 | 11.03 | 89 / 99 / 12 |
| (256, 30) unequal, seeds 1 to 50 | 200 | 197 | 199 | 6.56 | 6.92 | 51 / 128 / 21 |
| (496, 30) unequal, seeds 1 to 50 | 200 | 194 | 200 | 12.53 | 14.77 | 108 / 69 / 23 |
| All 30-feature fits | 860 | 703 | 825 | 7.61 | 9.68 | 358 / 434 / 68 |

By operator and linkage, identity covariance, seeds 1 to 50, fits with more than 1 cluster of 50 (unchanged,
then peel pass):

| Cells | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | --- | --- | --- | --- |
| 256 | 32, 46 | 50, 50 | 19, 49 | 50, 50 |
| 496 | 0, 28 | 50, 50 | 16, 45 | 49, 49 |

Mean found clusters in the same cells (unchanged, then peel pass): 256 cells 4.34 and 6.86, 5.70 and 5.70, 3.08
and 7.48, 5.86 and 5.92; 496 cells 1.00 and 6.84, 11.48 and 12.42, 5.02 and 13.56, 11.42 and 11.28.

For comparison, the third version of the draft definition of k with the peel pass (12 November): 0 of 880 fits
with more than 1 cluster.

## Predictions

- W1 held: 825 of 860 30-feature fits with more than 1 cluster under the peel pass against 703.
- W2 held: 358 fits gain clusters under the peel pass, 68 lose clusters.
- W3 not held: at (200, 80) the unchanged arm has more than 1 cluster in 3 of 20 fits and the peel pass in 7 of
  20. The prediction "1 cluster in 20 of 20 in both arms" contradicted the recorded result of 12 November (3 of
  20); that was an error in the protocol.
- W4 held in total (9.68 against 7.61, a factor 1.27) and not held for 496 cells with identity covariance
  (11.03 against 7.23, a factor 1.53; seeds 401 to 405: 12.60 against 7.50).

## Reading

- With today's default sibling test the peel pass makes the over-splitting of one Gaussian population worse:
  825 of 860 fits against 703.
- The increase is in the cells with average linkage. With Ward linkage nearly every fit is split in both arms.
- The stop at small peels was the reason the default cell (symmetric operator, average linkage) kept 496 cells
  of one Gaussian population at 1 cluster in 50 of 50 fits; with the peel pass it is 22 of 50.
- The peel pass is safe only together with a sibling rule that holds its level on the method's own trees. The
  third version did (0 of 880); today's default test does not.

## Not checked

- Limits other than 15 cells. Today's default test with the peel pass on the benchmark cases with groups and
  on the CHOIR cells. The current rule of the opt-in (table, k up to 15) with the peel pass on one Gaussian
  population. Populations that are not Gaussian. The peel pass exists at run time only.
