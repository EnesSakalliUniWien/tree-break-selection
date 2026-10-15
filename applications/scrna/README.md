# Single-Cell RNA Applications

The adult pancreas and Goncalves fetal-pancreas adapters live here. Dataset
loading, labels, manifests, and application reports stay local to this folder;
reusable space separation and plotting primitives come from
`tree_break_selection/`.

- `analysis/` contains dataset-specific biological comparisons and audit
  commands.
- `plots/` contains dataset-specific Python and R figure/report composition.
- The four top-level Python files below are the primary application entry
  points and orchestration commands.

```bash
python applications/scrna/pancreas_benchmark.py --help
python applications/scrna/goncalves_benchmark.py --help
python applications/scrna/analyze_space_decomposition.py --help
python applications/scrna/plot_pipeline.py --help
```

The two benchmarks, `analyze_space_decomposition.py`, the selected-NNLS reports and the
audits in `analysis/` default to a new `results/analyses/` folder. `plot_pipeline.py`, the
UMAP and R tree renderers, `analysis/analyze_goncalves_progenitors.py`,
`analysis/analyze_pancreas_inner_nodes.py`, `analysis/compare_pancreas_progenitor_signatures.py`
and `plots/goncalves_progenitor_coherence.py` work in place on one benchmark run and
require `--output-dir`; pass a `raw/assets/benchmark-results/` folder only to refresh
retained evidence.

Specialized scRNA figure generators used by `plot_pipeline.py` live in
`applications/scrna/plots/`; dataset-specific follow-up analyses live in
`applications/scrna/analysis/`. They are application code, not reusable
method or plotting-engine code. The R tree renderers share their application
tree-construction, branch-display, cluster-metadata, palette, and MRCA-audit
behavior through `plots/tree_plot_helpers.R`; dataset-specific composition
stays in each renderer.
