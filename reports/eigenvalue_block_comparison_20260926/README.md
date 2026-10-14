# Eigenvalue block and tree inference comparison

All 11 `method_proof` examples were evaluated on the same generated data and seeds across
seven linkage methods (average, complete, weighted, single, centroid, median, Ward) and
neighbor joining with MAD rooting. Each eigenvalue block was evaluated using raw Euclidean
coordinates and adaptive diffusion. There are 47 blocks and 752 attempted runs.

## Findings

- Non-common-mode adjacent boundaries have a median eigenvalue ratio of
  1.035, ranging from 1.003 to
  9.921 across 27 boundaries. These are mostly
  changes in decay regime, not isolated groups separated by large empty eigenvalue intervals.
- Between blocks, the median normalized rooted tree distance is 1.000
  over 1080 available pairs. Identical final labels need separate
  interpretation: 661 of
  931 available partition pairs are both single-cluster outputs.
- Changing the tree inference method within the same block and geometry has median rooted
  tree distance 0.483 over 1527 pairs.
- 525 runs returned clusters; 436 returned one cluster.
  Status counts: `{'ok': 525, 'failed_tree': 170, 'unsupported': 57}`. Failure and unsupported rows are
  excluded from partition comparisons, never converted to one-cluster assignments.
- Single-population examples produced more than one cluster in
  10/36 completed
  runs. These are dependent method/block evaluations, not an estimated false-positive rate.
- Raising the target rank from 40 to 80 changes grouping among the already-retained modes
  in 5/11 cases. Rank-80 tree inference was not run. See `rank_sensitivity.csv`.

## Figures and tables

- [Eigenvalue blocks](eigenvalue_blocks.png): every spectrum with existing block boundaries.
- [Benchmark cluster results](benchmark_clusters.pdf): one page per case, with spectra,
  TBS cluster counts and adjusted Rand index against generator labels for every method.
- [First binary example](binary_example.png).
- `summary_by_case.csv`: per-case results and denominators.
- `block_pairs.csv`: all pairwise eigenvalue separations; `adjacent_boundaries.csv`: adjacent gaps only.
- `geometry_pairs.csv`: matched sample-pair distance comparisons.
- `run_pairs_with_collapse.csv`: tree and partition comparisons with single-cluster matches flagged.
- `assignments.csv`, `trees.csv`, `runs.csv`: sample assignments, rooted clades and explicit outcomes.

## Per-case comparison

The tree distance is zero for identical rooted clades and one for no shared nontrivial clades.
Partition agreement is adjusted Rand index: one means identical partitions. All medians below
are descriptive; multiple pairs from one dataset are not independent observations.

| case | geometry | ok_runs | total_runs | median_between_block_rooted_rf | pairs_with_at_least_one_multicluster_output | median_ari_excluding_both_single_cluster |
| --- | --- | --- | --- | --- | --- | --- |
| bar_binary_balanced_4c | raw | 20 | 24 | 1.000 | 10 | 0.000 |
| bar_binary_balanced_4c | adaptive_diffusion | 18 | 24 | 1.000 | 2 | 0.000 |
| bar_binary_unbalanced_95_5 | raw | 20 | 24 | 1.000 | 8 | 0.000 |
| bar_binary_unbalanced_95_5 | adaptive_diffusion | 18 | 24 | 1.000 | 10 | 0.000 |
| mp_spike_below_bbp_continuous | raw | 18 | 24 | 1.000 | 10 | 0.000 |
| mp_spike_below_bbp_continuous | adaptive_diffusion | 18 | 24 | 1.000 | 10 | 0.000 |
| mp_spike_above_bbp_continuous | raw | 20 | 24 | 1.000 | 4 | 0.000 |
| mp_spike_above_bbp_continuous | adaptive_diffusion | 18 | 24 | 1.000 | 6 | 0.000 |
| cat_simplex_face_rare_20cat | raw | 50 | 64 | 1.000 | 34 | 0.000 |
| cat_simplex_face_rare_20cat | adaptive_diffusion | 48 | 64 | 1.000 | 53 | 0.000 |
| cat_dirichlet_overdispersed_20cat | raw | 32 | 40 | 1.000 | 19 | 0.000 |
| cat_dirichlet_overdispersed_20cat | adaptive_diffusion | 29 | 40 | 1.000 | 11 | 0.000 |
| cont_lowrank_pggn_shrinkage | raw | 17 | 32 | 1.000 | 0 | — |
| cont_lowrank_pggn_shrinkage | adaptive_diffusion | 15 | 32 | 1.000 | 0 | — |
| support_selected_only_fail_closed | raw | 19 | 24 | 1.000 | 10 | 0.000 |
| support_selected_only_fail_closed | adaptive_diffusion | 17 | 24 | 1.000 | 8 | 0.000 |
| phylo_brownian_null_16taxa | raw | 38 | 48 | 1.000 | 29 | 0.000 |
| phylo_brownian_null_16taxa | adaptive_diffusion | 36 | 48 | 1.000 | 15 | 0.000 |
| traversal_deep_signal_under_same_parent | raw | 26 | 32 | 1.000 | 12 | 0.000 |
| traversal_deep_signal_under_same_parent | adaptive_diffusion | 24 | 32 | 1.000 | 6 | 0.000 |
| traversal_deep_branch_recovery_stress | raw | 14 | 40 | 1.000 | 9 | 0.000 |
| traversal_deep_branch_recovery_stress | adaptive_diffusion | 10 | 40 | 1.000 | 4 | 0.000 |

### Single-population examples

| case | geometry | attempts | completed | multi_cluster_outputs | maximum_clusters |
| --- | --- | --- | --- | --- | --- |
| mp_spike_below_bbp_continuous | adaptive_diffusion | 24 | 18 | 5 | 5.000 |
| mp_spike_below_bbp_continuous | raw | 24 | 18 | 5 | 5.000 |

## Measurement contract

The existing cosine operator is formed from row-normalized raw feature values. The name
`binary` in the weighting API means raw values, including signed continuous features.
Blocks are the existing contiguous, tie-preserving BIC segments of the log eigenspectrum.
The target rank is 40, minimum block length 4, maximum segments 8. Full-operator energy
fractions use trace = number of nonzero sample rows; retained-energy fractions are separate.

For an earlier block A and later block B, boundary ratio = min(A)/max(B).
Log-Wasserstein distance gives equal mass to each eigenvalue and compares log eigenvalues;
for disjoint ordered blocks it equals the difference of their mean log eigenvalues. It is
unchanged by a common rescaling. Within-block widths are reported separately.
Euclidean coordinates are Q_B sqrt(Lambda_B); distance Spearman correlation and RMS-normalized
distance error compare the same sample pairs. Unequal block dimensions remain visible.
Disjoint eigenspaces are already orthogonal, so principal angles would not measure useful
differences in their sample clustering.

Tree distance = size of the symmetric difference of nontrivial rooted clade sets divided by
the sum of their sizes. Neighbor joining includes its MAD-rooting choice. Cluster comparisons
use all original samples in their original order. Ground-truth labels are used only for
evaluation; no block, method, alpha or geometry was selected using them.

TBS gates use the original feature matrix and its FeatureSpace, edge alpha 0.001 and sibling
alpha 0.01. This follows the existing runner defaults with support-threshold enforcement off
and branch-length variance policy `none`. It does not establish selected-tail calibration,
false-positive control, or a winning method. `calibration_counts` remain in `runs.csv`.
Adaptive diffusion uses K=15, time=3, 30 components, Euclidean distance, bandwidth
`-1/(d+2)` and median epsilon. Failed geometry is recorded explicitly.

IQ-TREE takes state alignments rather than these block distances and is excluded. NNLS changes
branch lengths on a fixed topology; branch-time variants and other gate profiles were not run.
This is the complete method_proof panel, not the repository's full benchmark suite, and there
is one seed per example. The rank sensitivity check does not establish resampling stability.

## Reproduction and verification

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg uv run --no-sync python -m benchmarks.diagnostics.spectral.adaptive_cosine.compare_spectral_blocks --suite method_proof --output-dir <new-output-directory>
```

`manifest.json` records full case definitions, package versions, Git state and source hashes.
Verification regenerated all 11 datasets, matched their hashes and retained eigenvalues,
checked every successful assignment's sample coverage, cluster count and adjusted Rand index,
and confirmed source hashes did not change during inference. Analytic checks covered
log-distance values, common-scale invariance, proportional and constant distance vectors,
and rooted-clade distance. The 18 space-separation and 16 construction/rooting tests passed.
No production method, inference default, commit or push was changed by this analysis.
