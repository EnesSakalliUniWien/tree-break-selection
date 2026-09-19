---
title: GO Plot Library Research 2026-09-19
type: source
status: reviewed
updated: 2026-09-19
sources:
  - raw/inbox/2026-09-19-go-plot-library-research.md
  - applications/endotypes/plots/trees/dendrogram_plot.py
  - applications/endotypes/plots/terms/axis_term_bar_plot.py
  - applications/endotypes/plots/terms/axis_term_heatmap.py
  - applications/endotypes/io/adaptive_diffusion_dataframe_io.py
  - applications/endotypes/reports/go_ic/adaptive_diffusion_plot_exports.py
tags:
  - plotting
  - context7
  - modularity
---

# GO Plot Library Research 2026-09-19

## Summary

Context7 research supports reusable Matplotlib Figure/Axes functions for the
adaptive diffusion GO-IC plots. Independent plot modules accept arrays or
DataFrames and return figures; file export and DataFrame serialization have
separate owners.

## Key Points

- Reuse Matplotlib Axes renderers for GO-term bars, heatmaps and embedding
  scatter plots across PNG and PDF output.
- Reuse SciPy dendrogram geometry and its returned leaf order for plain trees
  and cluster strips; the pipeline does not need a custom tree renderer.
- Seaborn heatmaps could replace labeling/colorbar boilerplate but introduce an
  undeclared dependency and different rendering. They were researched, not adopted.
- SciPy special functions could express binary entropy, but that numerical
  change is outside this plotting refactor and was not adopted.
- `build_axis_loadings` and `select_axis_term_heatmap` are pure DataFrame
  transformations. CSV/TSV operations live in
  `io/adaptive_diffusion_dataframe_io.py`.

## Evidence

- `raw/inbox/2026-09-19-go-plot-library-research.md` records Context7 IDs,
  official source URLs, checked local versions, and the adoption decisions.
- `plots/trees/dendrogram_plot.py` owns the shared SciPy dendrogram call.
- `plots/terms/axis_term_bar_plot.py` and `plots/terms/axis_term_heatmap.py` own the shared
  term renderers; `reports/go_ic/adaptive_diffusion_plot_exports.py` handles PNG export.

## Links

- [[current-adaptive-diffusion-subspace-tree-pipeline]]
