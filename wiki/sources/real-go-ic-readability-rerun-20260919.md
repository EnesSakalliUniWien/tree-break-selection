---
title: Real GO-IC Analyses with Readable Reports 2026-09-19
type: source
status: reviewed
updated: 2026-09-19
sources:
  - data/feature_matrices/feature_matrix_allGO_new_interactome.tsv
  - data/feature_matrices/feature_matrix_julia_allGO_new.tsv
  - results/analyses/allgo_new_interactome_go_ic_readability_20260919_170423/analysis_summary.md
  - results/analyses/allgo_new_interactome_go_ic_readability_20260919_170423/verification.json
  - results/analyses/julia_allGO_new_go_ic_readability_20260919_170423/analysis_summary.md
  - results/analyses/julia_allGO_new_go_ic_readability_20260919_170423/verification.json
  - applications/endotypes/plots/embeddings/embedding_comparison_plot.py
tags:
  - source
  - allgo
  - go-ic
  - readability
---

# Real GO-IC Analyses with Readable Reports 2026-09-19

## Summary

Reran both canonical real matrices with the saved rank-80, binary/TF-IDF
configuration and the current checkout. Both processes exited successfully.
These are new numerical analyses, not exact copies of the June outputs.

## Key Points

- Interactome: 339 genes by 5873 terms; all 12 subspaces completed, versus
  seven completed and five gate failures previously. TF-IDF modes 06–11
  remains top-ranked with 13 clusters, identical assignments, specificity
  0.753376, and GO-BIC/gene 1796.282946. The five formerly failed subspaces
  now return one cluster and are ranked in the degenerate tier.
- Julia allGO: 602 genes by 6368 terms; all 15 subspaces completed, versus
  fourteen completed and one gate failure previously. Binary modes 12–21
  remains top-ranked, now with 32 clusters versus 38, specificity 0.733687,
  and GO-BIC/gene 1872.654051. The formerly failed binary modes 02–06 block
  returns one cluster and is ranked in the degenerate tier.
- Other assignments also changed. Per-subspace comparisons are in each
  run's summary and verification JSON. The comparison does not establish
  the cause of changes or biological validity.
- Workflow PDFs contain 36 interactome pages and 45 Julia pages. Plot review
  confirmed readable top-ranked comparisons and explicit ranking inputs.
  Dense pages with hundreds of cluster IDs still require zoom.
- Full real-data layout checks exposed insufficient annotation-panel height
  in 35–39-mode blocks. Increased that height and regenerated workflow reports.
  All 27 embedding pages then passed bounds checks, with 12-point native PDF
  annotations. All 118 interactome and 145 Julia CSV hashes stayed unchanged
  during the report refresh.

## Evidence

Each run retains its exact command, input hash, saved settings, runtime versions,
and source hashes in `run_manifest.json`; execution output in `run.log` and
`exit_status.json`; numerical comparisons in `verification.json`; and annotation
bounds checks in `plot_layout_checks.json`. The layout regression test now
covers 39-mode summaries. Eighteen plot tests and four pipeline tests passed.

## Links

- [[current-adaptive-diffusion-subspace-tree-pipeline]]
- [[allgo-new-interactome-current-adaptive-diffusion-subspace-tree-20260618]]
- [[julia-allgo-new-current-adaptive-diffusion-subspace-tree-20260618]]
