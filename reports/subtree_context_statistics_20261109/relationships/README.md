# Eigenvalue magnitude and what a subtree preserves

This analysis reuses the 80 fixed datasets, 5,120 saved spectral Ward outcomes, and 160 native replays from the spectral alignment study. It regenerates the inputs to verify data hashes and recover covariance eigenvalues. It refits no clustering, topology, branch lengths, or diffusion kernel. The four simulation regimes are distinct from the four node-context classifications in the accompanying study.

## Raw covariance magnitude and factor coverage

For a union of complete positive sample covariance blocks, its raw variance share is the sum of its eigenvalues divided by the covariance trace. Eigenvalues are in squared feature units. We join each saved union's exact block ranges to the regenerated covariance eigenspectrum. We summarize paired contrasts within each dataset. The lowest-share covering union and maximum AMI are retrospective label-based descriptions of the searched family.

| Regime and factor omitted in the full Ward tree | Datasets with factor recovered in a union | Median lowest covering variance share | Range |
|---|---:|---:|---:|
| Crossed, B | 20/20 | .2885 | .2580–.3144 |
| Tied-population crossed, A | 5/7 | .4454 | .4366–.4750 |
| Tied-population crossed, B | 6/13 | .4452 | .4231–.4948 |

In the crossed regime, PC2 alone recovers B in all 20 datasets with two-cluster AMI 1. Full coordinates and the leading two-block span both omit B and have B AMI −.0117, even though the leading span retains median variance share .9300. In replicate 0, the leading span retains .9274 of raw covariance trace. PC2 alone retains .2652 and gives the crossing factor split. The full tree still preserves the four joint groups. The loss here concerns the factor split displayed by a single binary tree.

The B-covering alternatives also fit the full-space Euclidean target worse than the full tree in this shortlist. Across 20 crossed datasets, the best B-covering shortlisted tree has a median additive RMSE increase of 1.3962 distance units, with range 1.3035–1.4473. In replicate 0 the full tree RMSE is .5568 and the best B-covering shortlist RMSE is 1.9638. Selecting by lower full-target residual therefore answers a different question from preserving a crossing factor.

In the compatible regime, the leading span retains median raw covariance share .8022, yet preserves all four joint clades in only 6/20 datasets. Full coordinates preserve them in 20/20. The nested B split is absent from the leading span in eight datasets. Their lowest covering union has median share .3190. This is compatible with losing a grouping through projection even when the two target splits can coexist in a single rooted tree.

## Diffusion weights and native test scales

The diffusion window reports the retained operator eigenvalues. At time three, each retained mode's coefficient in squared distance is `max(lambda, 0)^6`. The leading whole block's coefficient share sums those coefficients over the block. This is an operator-weight summary. Actual pair distance energy also includes the spread of the mode's coordinate over rows. For a coordinate column (z), the sum of squared distances across unordered pairs is (nsum_i(z_i-ar z)^2). Exact tied blocks are summarized whole, so the main comparisons do not depend on a choice of individual vector within such a block.

In the crossed replicate 0, the leading diffusion block has three modes and .7628 of coefficient mass in both operator arms. Its share of actual squared pair distances is .6820 for symmetric/average and .6807 for Markov/Ward. The eight representative traces reproduce the full sum of saved squared distances using this decomposition. The two operators have identical leading block coefficient shares on every paired dataset. Nevertheless, their root splits differ in 54/80 pairs. In crossed replicate 6, both have leading block coefficient share .7569, while recorded root p-values are (5.79\times10^{-20}) and (3.51\times10^{-285}). The comparison changes coordinate operator and linkage together. It does not isolate either as the cause.

Native node-context PCA eigenvalues measure a separate, whitened construction. In the 160 roots, all six are near .984375 and the largest within-root spread is (8.38\times10^{-12}). The raw MP count is zero and the context test dimension is six throughout. All 160 saved root decisions reject, including both configurations on each of the 20 common-mean null datasets. Similar native context spectra therefore coexist with both genuine signal and selected null rejection in this panel. Nineteen stored root p-values are numerical zeros and are retained as such. The eight representative traces record the raw sibling statistic, reference scale, degrees of freedom, fitted inflation, and branch variance multiplier separately. In crossed replicate 0, the two raw statistics are 95.19 and 92.66, with branch multipliers 1.209 and 1.242. In null replicate 0, they are 145.89 and 97.44, with multipliers 1.191 and 1.582 and applied inflation factors 2.740 and 2.282.

These comparisons do not supply an eigenvalue threshold for testing or a probability guarantee. The 63 unions for a dataset reuse its same observations and count as one dataset for between-dataset interpretation. The native configurations are paired on the same 80 datasets. This descriptive analysis evaluates no hypothesis test over the candidate rows.

## Reproduce

From the project root, run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m reports.subtree_context_statistics_20261109.relationships.run`. The manifest records input/source and output SHA-256 hashes. `alternatives_with_magnitudes.csv.gz` preserves the 5,120 joined rows, `factor_pairs.csv` gives 120 dataset-by-factor summaries, and `native_pairs.csv` aligns the two native configurations on each of 80 datasets. The eight trace-level variance records are in `root_trace_scales.csv`.
