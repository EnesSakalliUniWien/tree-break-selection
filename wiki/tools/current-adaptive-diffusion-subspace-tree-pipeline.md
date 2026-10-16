---
title: Adaptive Diffusion Cosine Subspace Clustering and GO-IC Analysis
type: tool
status: reviewed
updated: 2026-10-16
sources:
  - tree_break_selection/eigenvalue_groups.py
  - tree_break_selection/space_separation/diffusion.py
  - reports/mathematical_construction_review_20260928/rank_invariance_repair.md
  - reports/eigenspace_reassessment_20260928/README.md
  - tree_break_selection/space_separation/adaptive_cosine.py
  - applications/endotypes/pipelines/run_adaptive_diffusion_cosine_subspace_clustering_go_ic_analysis.py
  - applications/endotypes/pipelines/adaptive_diffusion_go_ic_workflow.py
  - applications/endotypes/analysis/go_ic/adaptive_diffusion_go_ic_terms.py
  - applications/endotypes/analysis/go_ic/adaptive_diffusion_go_ic_ranking.py
  - applications/endotypes/plots/terms/axis_term_bar_plot.py
  - applications/endotypes/plots/terms/axis_term_heatmap.py
  - applications/endotypes/plots/trees/dendrogram_plot.py
  - applications/endotypes/plots/trees/tree_cluster_plot.py
  - applications/endotypes/plots/embeddings/subspace_embedding_plot.py
  - applications/endotypes/plots/embeddings/diffusion_embedding_plot.py
  - applications/endotypes/plots/embeddings/embedding_comparison_plot.py
  - applications/endotypes/plots/shared/_cluster_legend.py
  - applications/endotypes/plots/ranking/subspace_ranking_plot.py
  - applications/endotypes/io/adaptive_diffusion_dataframe_io.py
  - applications/endotypes/reports/go_ic/adaptive_diffusion_plot_exports.py
  - applications/endotypes/reports/go_ic/adaptive_diffusion_subspace_reports.py
  - applications/endotypes/reports/go_ic/adaptive_diffusion_subspace_artifacts.py
tags:
  - pipeline
  - allgo
  - adaptive-diffusion
  - subspace
---

# Adaptive Diffusion Cosine Subspace Clustering and GO-IC Analysis

## Summary

Use `applications/endotypes/pipelines/run_adaptive_diffusion_cosine_subspace_clustering_go_ic_analysis.py`
to construct adaptive diffusion cosine subspace trees, apply Tree-Break
Selection clustering, and run extended GO information-criterion analysis,
including coherence, TF-IDF quality, and specificity-aware ranking. It creates
the same directory structure for every input matrix: a timestamped experiment
root under `results/analyses/`, method-separated PDFs, connected manifests,
rankings, and one artifact-complete directory per cosine subspace.

## Usage

Run the pipeline with an explicit feature matrix:

```bash
MPLBACKEND=Agg python applications/endotypes/pipelines/run_adaptive_diffusion_cosine_subspace_clustering_go_ic_analysis.py \
  --input /path/to/feature_matrix.tsv
```

The script name describes its responsibilities. Output prefixes and directory
names retain the existing artifact contract for report readers and saved runs.

The entry point owns CLI parsing and delegates execution to
`pipelines/adaptive_diffusion_go_ic_workflow.py`. The workflow uses separate
modules for GO-term metadata and loadings (`analysis/go_ic/adaptive_diffusion_go_ic_terms.py`),
GO-IC and specificity ranking (`analysis/go_ic/adaptive_diffusion_go_ic_ranking.py`),
independent plotting types under `plots/`, PDF assembly
and experiment figures (`reports/go_ic/adaptive_diffusion_subspace_reports.py`),
and artifact paths, indexes, manifests, and README export
(`reports/go_ic/adaptive_diffusion_subspace_artifacts.py`).

The default output root is:

```text
results/analyses/<matrix-stem-without-feature_matrix_>_current_adaptive_diffusion_subspace_tree_<timestamp>/
```

Every full run should write these top-level artifacts:

- `rankings/current_adaptive_diffusion_subspace_tree_ranking.csv`
- `rankings/current_adaptive_diffusion_subspace_tree_specificity_aware_ranking.csv`
- `rankings/current_adaptive_diffusion_subspace_tree_subspace_blocks.csv`
- `rankings/current_adaptive_diffusion_subspace_tree_spectrum.csv`
- `ARTIFACT_INDEX.md`
- `artifact_index.csv`
- `subspace_plot_index.csv`
- `connected_results_manifest.json`
- `<artifact-prefix>_quality_aware_go_ic_by_method/current__adaptive_diffusion_cosine_subspace/`
- `<artifact-prefix>_quality_aware_go_ic_plots/`
- `subspaces/<weighting>/<block_name>/`

Each `subspaces/<weighting>/<block_name>/` directory should contain the
subspace coordinates, diffusion metadata, linkage tree, cluster assignments or
explicit failure status, GO-IC quality summary, coherence tables, TF-IDF
quality tables, axis term-loading CSVs, per-axis term plots, combined term
plot, and plain plus term-annotated embeddings.

The preferred reading order is `specificity_aware_rank`: quality tier,
specificity score, specific-cluster fraction, weighted specificity delta, then
GO-BIC active per gene. `raw_go_ic_rank` preserves the raw GO-IC-only order
for audit. The `old_display_rank` snapshot column was removed on 2026-10-16;
retained artefacts written before then still carry it.

## Plot and DataFrame interfaces

Each plot type has its own module, grouped by responsibility:
`plots/terms/` contains bars and heatmaps; `plots/trees/` contains dendrograms
and cluster strips; `plots/embeddings/` contains subspace and diffusion
embeddings; `plots/ranking/` contains ranking bars. `plots/shared/` contains
formatting shared across plot types. Their `create_*` functions return Matplotlib figures
without reading or writing files. Callers own closing the returned figures.
Shared draw functions preserve the PNG/PDF chart semantics.

Embedding comparison pages use full-width unannotated plot panels and separate
12-point GO-term summaries, retaining native text in PDFs. Standalone annotated
embeddings use 11-point summaries. Tree and embedding figures share the existing
discrete cluster palette and explicit integer legends, with room allocated for
every cluster ID. Ranking figures show quality tiers, specificity scores, and
GO-BIC on aligned rows, preserving the computed specificity-aware order.

`build_axis_loadings` returns `(all_loadings, top_loadings)` DataFrames.
`select_axis_term_heatmap` returns the selected pivot table. Neither persists
files. `io/adaptive_diffusion_dataframe_io.py` owns table input/output;
`reports/go_ic/adaptive_diffusion_plot_exports.py` owns PNG saving and figure closure.
The workflow orchestrates these boundaries, and report assembly owns PDFs.

See [[go-plot-library-research-20260919]] for the Context7 research and library
adoption decisions. Existing numerical routines and artifact names are retained.

## Evidence

### Shared rank-group contract (2026-09-28)

The numerical-tie rule now also covers diffusion, local PCA, MP test dimensions
and Wald projection; see [[eigenvector-calibration-construction-review]]. Both
`--max-rank` and `--diffusion-components` are targets that can expand through a
complete group. Diffusion removes an isolated leading mode, but retains a tied
leading group whole. Actual dimensions are recorded in diffusion metadata.
Full `make check` passed 693 tests for the rank repair. The matched 11-case
spectral-block benchmark and a separate calibration study of the registered
methods have since completed; see
[[eigenspace-benchmark-calibration-reassessment]]. The calibration study found
excess false splitting in the tested settings. The earlier saved GO application
outputs below remain historical and were not regenerated; the registered-method
study does not establish calibration of this complete GO application workflow.

### Eigenvalue numbering and grouping review (2026-09-26)

Cosine eigenvalues and eigenvector columns are sorted together in descending
order. Mode numbers are one-based; block endpoints are inclusive and extraction
uses `[block_start - 1:block_end]`. Block IDs are zero-based identifiers.
The algorithm groups contiguous log-spectrum decay regimes using linear-fit
SSE and a BIC penalty; it does not cluster eigenvalues by proximity. Each block
supplies diffusion distances for an average-linkage tree; TBS tests use the
original feature matrix to obtain gene cluster assignments.

Replaying the four saved September 19 spectra reproduced all 27 block ranges
with complete coverage of retained modes 1–80, without overlaps. None of these
saved internal boundaries has equal eigenvalues. Forty small random spectra
also matched exhaustive enumeration of the implemented BIC objective, and the
six existing space-separation tests passed.

The review reproduced a defect where segmentation cuts and `max_rank` truncation
did not preserve repeated-eigenvalue eigenspaces. For the descending spectrum
`[10,9,8,7,6,6,6,6,6,6,2,1]`, minimum length 4 and maximum 8 segments produce
blocks 1–8 and 9–12, splitting the eigenvalue-6 eigenspace. Rotating eigenvectors
8 and 9 by pi/4 preserved the reconstructed operator to about 1e-16 but changed
first-block pairwise distances by up to 2.44949. Thus downstream geometry can
depend on the eigensolver's basis inside a repeated eigenspace. The saved
rank-80 spectra alone cannot establish the gap to discarded mode 81.

The fix uses adjacent numerical ties (relative tolerance 1e-10 and absolute
tolerance 1e-12 times the largest retained eigenvalue), joined transitively.
Segmentation and common-mode isolation cannot cut these groups. `max_rank` is
now a target that extends through the whole boundary group; an entirely tied
positive spectrum is retained in full. Minimum segment length still applies,
so fewer blocks may be selected. Regression tests cover exact and numerical
ties, rotation-invariant block distances, common-mode isolation, and rank
extension. Replaying all 27 saved real-analysis blocks leaves their boundaries
unchanged. This fixes cosine separation; it does not change diffusion-component
truncation or certify the resulting gene clusters.

### Matched benchmark after the tie fix (2026-09-26)

`results/analyses/eigenvalue_tie_benchmark_20260926/report.md` records a fresh
comparison of committed `f34f59ca` cosine separation against the fixed working
module, sharing all downstream code. The run includes all 42 binary-suite cases,
both real GO matrices with binary/TF-IDF weighting, and a tied-rank stress case:
45 datasets, 47 dataset/weighting configurations and 342 block executions.

All 169 ordinary matched blocks have identical coordinates. All 167 available
distance/tree pairs and all 155 successful label pairs are exactly equal,
including all 27 real-analysis blocks. Each version has 12 unsupported binary
blocks and two binary diffusion failures. Synthetic truth ARI is unchanged on
all 128 jointly successful binary blocks. The rank stress extends 10 modes to
40, preserving the 39-dimensional repeated-eigenvalue tail; both tails yield
one cluster and both common modes fail diffusion because of duplicate
coordinates. Its different tail ranges are not counted as matched blocks.

The separate tied-spectrum rotation check reduces block-distance sensitivity
from 2.44949 to 4.44e-16. Total block times are 1568.51 seconds before and
1562.92 after; single-run timing does not establish a performance improvement.
Input hashes, module snapshots, scripts, per-block outputs and verification are
retained beside the report. GO-IC metrics and plots were not regenerated.

- `applications/endotypes/pipelines/run_adaptive_diffusion_cosine_subspace_clustering_go_ic_analysis.py`
  derives the experiment root and artifact prefix from the input matrix name
  and writes the full connected artifact structure in one run.

## Links

- [[julia-allgo-new-current-adaptive-diffusion-subspace-tree-20260618]]
- [[allgo-new-interactome-current-adaptive-diffusion-subspace-tree-20260618]]
