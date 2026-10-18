---
title: Eigenvalue Block Differences Across Tree Inference Methods
type: analysis
status: draft
updated: 2026-10-13
sources:
  - benchmarks/diagnostics/spectral/adaptive_cosine/compare_spectral_blocks.py
  - tree_break_selection/space_separation/adaptive_cosine.py
  - tree_break_selection/tree/construction/build.py
  - reports/eigenvalue_block_comparison_20260926/README.md
  - reports/eigenvalue_block_comparison_20260926/manifest.json
  - reports/eigenvalue_block_comparison_20260926/adjacent_boundaries.csv
  - reports/eigenvalue_block_comparison_20260926/rank_sensitivity.csv
  - reports/eigenvalue_block_comparison_20260926/summary_by_case.csv
  - reports/eigenspace_reassessment_20260928/spectral_blocks/comparison_summary.json
tags:
  - eigenvalues
  - benchmark
  - tree-inference
---

# Eigenvalue Block Differences Across Tree Inference Methods

## Summary

The complete 11-case `method_proof` panel shows mostly narrow eigenvalue
boundaries between non-common `SpectralBlock` objects, despite substantial
differences in the trees inferred from their coordinates. Tree differences
must be separated from TBS partition agreement because most completed outputs
contain one cluster. The grouping also depends on the retained-rank target.

The numerical counts below describe the original 2026-09-26 snapshot. The
2026-09-28 rank-repair rerun preserves all 752 statuses but changes 18/582
available trees and 33/525 successful partitions. Single-cluster outputs rise
from 436 to 459; ARI improves in 15 runs and worsens in 15. See
[[eigenspace-benchmark-calibration-reassessment]] for the matched comparison and
independent null studies.

## Details

The comparison uses identical generated data and seeds across average,
complete, weighted, single, centroid, median and Ward linkage, and neighbor
joining with MAD rooting. Raw block Euclidean distances and adaptive diffusion
are evaluated separately. Gates use the original feature matrix and its
`FeatureSpace`, edge alpha 0.001 and sibling alpha 0.01. The standard runner's
support-threshold enforcement remains off; no calibration claim follows.

The existing algorithm constructs contiguous BIC segments of the log cosine
eigenspectrum, preserving numerical ties. It does not group eigenvalues by
proximity. Tree inference takes place after this eigendecomposition, so
changing the tree inference method cannot change the input eigenvalues.

Measurements separate four questions:

- Eigenvalue separation: adjacent boundary ratio, within-block width,
  log-Wasserstein distance and full versus retained operator energy.
- Sample geometry: matched pairwise-distance rank correlation and
  RMS-normalized distance error.
- Tree difference: normalized symmetric difference of nontrivial rooted
  clades, with zero for identical sets and one for disjoint sets.
- Partition difference: adjusted Rand index on the same samples, with
  both-single-cluster pairs flagged separately.

Disjoint spectral blocks already have orthogonal eigenspaces. Principal
angles between those spaces would therefore not answer whether their sample
distances, inferred trees or TBS partitions differ usefully.

## Evidence

- Executed 47 blocks and 752 block/geometry/tree-method combinations:
  525 returned clusters, 170 failed tree construction and 57 were unsupported.
  All construction failures were rejected nonmonotone centroid/median heights.
  Unavailable outcomes were never converted to single-cluster assignments.
- Across 27 non-common adjacent boundaries, the median eigenvalue ratio is
  1.035, with range 1.003–9.921. Two ratios exceed 2: the Brownian example's
  modes 11–16 versus 17–26, and the deep-branch example's modes 2–6 versus 7–10.
- Between blocks, median normalized rooted distance is 1.000 across 1,080
  available tree pairs. Within the same block and geometry, changing the tree
  inference method gives median distance 0.483 across 1,527 pairs.
- Of 525 completed runs, 436 return one cluster; 661 of 931 available
  between-block partition comparisons are both single-cluster outputs.
  These dependent pair counts are descriptive, not independent replicates.
- The single-generator-cluster continuous example produces multiple clusters
  in 10 of 36 completed method/block runs. This is not an estimated error rate.
- Increasing the target rank from 40 to 80 changes grouping among the first
  40 modes in 5 of 11 cases. Rank-80 trees were not rerun.
- Regenerated all 11 inputs and matched their hashes and retained spectra;
  verified successful assignment coverage, counts, adjusted Rand indices and
  unchanged inference source hashes. The 18 space-separation tests and 16
  construction/rooting tests passed. The report includes an 11-page figure deck.

## Links

- [[current-adaptive-diffusion-subspace-tree-pipeline]]
- [[adaptive-cosine-kak-benchmark-probe-20260605]]
- [[tree-construction-method-map]]
- [[benchmark-pipeline-contract]]
- [[eigenspace-benchmark-calibration-reassessment]]

## Open Questions

- Which block boundaries and partitions persist under independent data
  replicates and rank changes?
- How much do the five changed rank-80 groupings alter tree and TBS results?
- How do branch-time and other gate profiles affect these comparisons?
  This panel fixes the branch-length variance policy to `none`; it does not
  compare NNLS branch-time variants.
