# Alternative spectral subspaces can expose crossing splits

Each of the 80 fixed subtrees yielded 63 distinct Ward topologies from 63 unions of its six positive singleton sample eigenvalue blocks. The full-coordinate Ward topology matched the all-block union for every dataset. The experiment examined one deterministic Ward result per space at n=64. The four-leaf examples additionally enumerate all Ward tie histories under the declared tolerance.

## Signal recovery

Counts are out of 20 datasets per regime. A split is represented when either side is a clade. The two-factor union counts cases where both binary contrasts occur somewhere in the spectral family. The same-tree column requires both in one tree. In the compatible regime, contrast B is a nested 16-cell group versus its complement, while A separates the two 32-cell halves.

| Regime | Full split A | Full split B | Full joint four clades | Union has both | One tree has both |
|---|---:|---:|---:|---:|---:|
| Crossed | 20 | 0 | 20 | 20 | 0 |
| Tied-population crossed | 13 | 7 | 20 | 11 | 0 |
| Compatible | 20 | 20 | 20 | 20 | 20 |

The crossed factors each partition all 64 leaves into two groups, with four nonempty intersections. Both splits cannot occur in one rooted tree. Full coordinates recovered the four joint groups in all 20 crossed datasets, while the missing second factor appeared under another subspace. For the tied-population simulation, finite-sample spectra had six distinct numerical blocks, and the unions exposed both factor splits in 11 datasets. The compatible example preserved both contrasts in the full tree, while the leading rank-two span retained B in 12 datasets and all four clades in six.

At the fixed two-cluster cut in the crossed panel, full-coordinate mean AMI was 1.000 for factor A and −.012 for B. The retrospective best-over-unions mean AMI was 1.000 for each. In the tied-population panel the corresponding full means were .646 and .342, and the oracle means were .946 and .797. The oracle evaluates a different best union for each target and dataset, using truth labels. It gives a descriptive upper bound for this searched family. Matched F1 and the fixed four-cluster results are in the row table. All simulated signal groups are balanced, so the rare-label summary reflects a deterministic label tie break.

## Distance fit and exact examples

The common-target NNLS shortlist contained eight distinct unrooted trees per dataset. All fits used the same original Euclidean target for a dataset. Relative to the full tree, minimum residual in this shortlist improved in 6/20 crossed, 7/20 tied-population, 1/20 compatible, and 0/20 null datasets, using an absolute improvement threshold of 1e-9. The selected trees preserved the same exact factor-split and joint-clade counts as the full trees. In particular, a smaller residual on the full-space target did not recover the second crossed contrast.

For the exact rectangle, the two rank-one blocks expose opposite factor splits. The full-coordinate target has additive RMSE .1927 for its Ward topology and 1.0716 for the crossing alternative. For the exact square, the only intrinsic nonempty spectral space is the whole tied rank-two block. Its Ward tie family contains two crossing alternatives with equal additive RMSE .3382 to numerical precision. Prescribing the coordinate axes externally gives one of those trees per axis. This axis choice does not follow from the tied spectrum.

The common-mean null still yields 63 distinct topologies per dataset. It has no factor labels or recovery scores. The number of alternative trees describes the algorithm's sensitivity to its distance geometry. It does not count populations or calibrate any test.

## Evidence

- `protocol.md` records the design frozen before scoring.
- `source_during_run.py.txt` preserves the exact executed source. `run.py` exposes
  `REGIMES` and `dataset(regime, replicate)` with import-only cleanup documented
  in `source_change_review.json`; the non-import AST is unchanged.
- `rows.csv.gz` records all 5,120 candidate outcomes, including seeds, dataset hashes, topology hashes, recovery and common-target residuals where fitted.
- `summary.csv` gives the 80 dataset-level summaries.
- `evidence.json` records source hashes and the exact rectangle and square examples. `verification.json` hashes outputs and records the verified row counts.

Run with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python reports/subtree_spectral_alignment_20261109/alternatives/run.py`. All 15 focused baseline geometry and profile tests passed. No production source or default changed.
