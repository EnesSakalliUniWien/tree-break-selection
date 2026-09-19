---
title: Adaptive Diffusion Cosine Subspace Clustering and GO-IC Analysis
type: tool
status: reviewed
updated: 2026-09-19
sources:
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
GO-BIC active per gene. `old_display_rank` preserves the previous
quality-tiered GO-IC order, and `raw_go_ic_rank` preserves the raw GO-IC-only
order for audit.

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

- `applications/endotypes/pipelines/run_adaptive_diffusion_cosine_subspace_clustering_go_ic_analysis.py`
  derives the experiment root and artifact prefix from the input matrix name
  and writes the full connected artifact structure in one run.

## Links

- [[julia-allgo-new-current-adaptive-diffusion-subspace-tree-20260618]]
- [[allgo-new-interactome-current-adaptive-diffusion-subspace-tree-20260618]]
