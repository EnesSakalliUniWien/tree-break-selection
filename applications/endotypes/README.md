# Endotype and GO-Annotation Application

This application groups all user-facing feature-matrix, GO-annotation,
subspace, endotype/reference, and report commands.

Start with the orchestrator:

```bash
python applications/endotypes/pipelines/run_go_annotation_feature_matrix_pipeline.py --help
```

The adaptive diffusion cosine subspace clustering and GO-IC analysis pipeline is
`pipelines/run_adaptive_diffusion_cosine_subspace_clustering_go_ic_analysis.py`. It consumes the
main `tree_break_selection.space_separation` interface. It constructs diffusion
trees, applies Tree-Break Selection clustering, then evaluates GO information
criteria, coherence, TF-IDF quality, and specificity-aware ranking.

The entry point parses CLI arguments and calls
`pipelines/adaptive_diffusion_go_ic_workflow.py`. GO-term loadings and ranking
live in `analysis/go_ic/adaptive_diffusion_go_ic_terms.py` and
`analysis/go_ic/adaptive_diffusion_go_ic_ranking.py`; subspace plots live in
separate modules grouped under `plots/terms/`, `plots/trees/`,
`plots/embeddings/`, and `plots/ranking/`. PDF assembly
lives in `reports/go_ic/adaptive_diffusion_subspace_reports.py`, while artifact paths,
indexes, manifests, and the experiment README live in
`reports/go_ic/adaptive_diffusion_subspace_artifacts.py`.

Plot constructors return Matplotlib figures from DataFrames/arrays without file
I/O. `reports/go_ic/adaptive_diffusion_plot_exports.py` saves and closes PNG figures;
`io/adaptive_diffusion_dataframe_io.py` reads and writes tables. Analysis
functions such as `build_axis_loadings` return DataFrames for reuse in scripts.

For example, once `all_loadings` is available:

```python
from pathlib import Path
from applications.endotypes.plots.terms.axis_term_heatmap import create_axis_term_heatmap
from applications.endotypes.io.adaptive_diffusion_dataframe_io import write_dataframe
from applications.endotypes.reports.go_ic.adaptive_diffusion_plot_exports import write_figure

fig = create_axis_term_heatmap(all_loadings)
write_figure(fig, Path("terms.png"), dpi=220)
write_dataframe(all_loadings, Path("loadings.csv"))
```

Responsibilities are explicit:

- `pipelines/` owns executable workflows and orchestration.
- `analysis/` owns data-quality, correctness, and biological interpretation.
- `analysis/go_ic/` owns GO-term loading calculations and GO-IC ranking.
- `io/` owns DataFrame file input/output for the modular GO-IC workflow.
- `reports/` owns tabular and PDF artifact export.
- `reports/go_ic/` owns GO-IC figure exports, PDF assembly, and artifact indexes.
- `plots/` owns endotype-specific visualization commands.
- `plots/terms/` owns GO-term bar charts and heatmaps.
- `plots/trees/` owns dendrograms and cluster-annotated trees.
- `plots/embeddings/` owns subspace and diffusion embeddings and their shared scatter renderer.
- `plots/ranking/` owns subspace ranking figures.
- `plots/shared/` owns formatting used by multiple plot types.
- `_shared.py` is the single internal seam for matrix naming, loading, reference
  endotype parsing, and symbol-to-Entrez resolution shared across those command
  categories.

Reference endotype tables belong in `data/reference/`. Canonical input feature
matrices belong in `data/feature_matrices/`.

The `plots/kak_signal_adaptive_umap_tree_page.py` command assembles KAK/cosine
geometry and tree pages from matrix-probe outputs. Application commands may
compose domain reports; reusable separation and plotting primitives belong in
`tree_break_selection/`.
