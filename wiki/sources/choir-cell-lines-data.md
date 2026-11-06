---
title: CHOIR Pooled Cancer Cell-Line Data
type: source
status: reviewed
updated: 2026-11-06
sources:
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-20261106.md
  - raw/inbox/choir-cell-lines/README.md
  - raw/inbox/choir-cell-lines/GSE157220.series.soft.txt
  - raw/inbox/choir-cell-lines/GSM4759276.sample.soft.txt
  - raw/inbox/choir-cell-lines/GSE157220_Pool_composition.xlsx
  - raw/inbox/choir-cell-lines/zenodo-14641222.json
  - raw/inbox/choir-cell-lines/Kinker_2020_cancer_cell_lines.R
  - raw/inbox/choir-cell-lines/provenance.json
  - raw/inbox/choir-cell-lines/pool-count-validation.json
  - raw/inbox/choir-cell-lines/validation.json
  - raw/inbox/choir-cell-lines/tbs-preparation.json
  - raw/inbox/choir-cell-lines/tbs-method-selection.json
  - raw/inbox/choir-cell-lines/tbs-preparation-script.py.txt
  - raw/inbox/choir-cell-lines/tbs-run.md
  - raw/inbox/choir-cell-lines/tbs-summary.json
  - raw/inbox/choir-cell-lines/tbs-diffusion-verification.json
  - raw/inbox/choir-cell-lines/tbs-run-verification.json
  - tree_break_selection/space_separation/diffusion.py
  - tests/statistics/55_test_eigenspace_invariance.py
  - raw/inbox/choir-cell-lines/tbs-diagnosis.md
  - raw/inbox/choir-cell-lines/tbs-root-variance-diagnostic.json
  - raw/inbox/choir-cell-lines/tbs-structure-diagnostic.json
  - raw/inbox/choir-cell-lines/tbs-recovery-comparison.md
  - raw/inbox/choir-cell-lines/tbs-recovery-protocol.json
  - raw/inbox/choir-cell-lines/tbs-recovery-comparison.csv
  - raw/inbox/choir-cell-lines/tbs-recovery-validation-null.csv
  - raw/inbox/choir-cell-lines/tbs-recovery-verification.json
  - raw/inbox/choir-cell-lines/published-choir-results.md
tags:
  - scrna
  - data
  - choir
  - cancer
---

# CHOIR Pooled Cancer Cell-Line Data

## Summary

The pooled cancer cell-line dataset from CHOIR Figures 3–4 is loaded locally
as sparse AnnData: 48,879 cells × 22,722 genes, with 190 cell lines and 22 cancer
types. Expression is the deposited GSE157220 CPM matrix; cell metadata and
published clustering labels come from Zenodo 14641222.
The subsequent full-data family-diffusion TBS run returned one cluster because
neither root edge passed its gate; it did not recover the known cell lines.
The controlled coordinate/linkage alternatives subsequently failed full-data
NNLS convergence, and independent count simulations found calibration failures.
After the singleton-peel merge, the same tree gives 70 clusters (ARI 0.027593);
see [[choir-singleton-peel-rerun]].

## Key Points

- Ready-to-load file:
  `local_data/choir/kinker_2020/kinker_2020_choir_cells_cpm.h5ad`.
  The source note contains an absolute-path, backed-mode loading example.
- All 109 published metadata columns are preserved. `CHOIR_1` is the imported
  default result with 367 clusters; no clustering was performed in this load.
- All 48,879 authors' cell IDs map uniquely to GEO. Independent GEO
  pool-composition counts agree for every one of the 190 cell lines.
- The original 53,513-cell GEO download remains intact. Its extra 4,634 cells
  belong to the eight custom-pool cell lines absent from the CHOIR label table.
- Full expression and metadata round-trip checks passed; source/output hashes
  and environment versions are recorded. The Zenodo workbook's published MD5
  was verified. Bulk files are local-only and absent from fresh clones.
- The authors' reprocessed count matrix is different from this deposited CPM
  input. The original load applied no transformations; the subsequent TBS run
  used a separate transformed PCA representation and retained the CPM matrix.

## Evidence

### Published CHOIR results

The published article reports that CHOIR distinguished all 190 lines in 48,879
cells, uniquely among the six compared methods under the evaluated settings
(Figure 3). It also performed better when each line was separately downsampled
to 50 cells. In 74.3% of lines divided into multiple clusters, proliferation
module scores differed significantly; A375 and T47D state differences were
supported in independent datasets (Figure 4).

The cancer benchmark uses entropy of cluster accuracy and a primary-cluster
majority/uniqueness criterion. Recovery of every line does not mean perfect
per-cell purity or exactly 190 clusters. The imported default `CHOIR_1` has 367
clusters; ARI 0.664014 against line IDs is our local rescoring, distinct from
the paper's cancer benchmark. Within-line states can reduce that ARI. CHOIR was
not rerun, and the deposited CPM workflow differs from the authors' processed
counts. `published-choir-results.md` records the source, sections and access
limits, keeping published findings separate from local calculations.

### TBS execution preparation, 2026-11-03

The user selected the highest performing available diffusion TBS method.
The latest full benchmark selects `tbs_family_diffusion` by mean ARI:
0.825843 across 122 supported cases. It ties `tbs_diffusion_adaptive_nnls`
at 0.823085 on their 101 common supported cases. This selection does not
establish performance on the cancer cell-line dataset.

All 48,879 cells were prepared with `log1p(CPM / 100)`, 2,000 Seurat-flavor
variable genes, scaling clipped at 10, and 30 PCs with seed 0. No library-total
renormalization or cell subsampling was applied. The selected preset's
15-neighbor variable-bandwidth kernel has three connected components of
48,444, 250, and 185 cells. Current `compute_diffusion_coordinates` therefore
falls back to dense eigendecomposition. One full float64 square array is
17.8006 GiB, versus about 12 GiB available at preparation; multiple such arrays
are required. This initially blocked execution; the SciPy component solution
below subsequently enabled the full run.

Reusable local checkpoints are under
`results/analyses/choir_cell_lines_tbs_20261103_163359/`: `pca.npy`,
`pca_loadings.npy`, `gene_selection.csv`, `cell_metadata.csv`,
`diffusion_kernel.npz`, and `graph_components.npy`. The three compact
`tbs-*` source captures above preserve settings, hashes, graph evidence,
and the preparation script. Nineteen focused diffusion/dispatcher tests passed.

### Full TBS result and library solution

Installed SciPy now solves disconnected diffusion graphs by component and
merges their spectra before applying the global complete-tie rule. Dense-array
reference tests include repeated cross-component eigenvalues, an isolated
vertex and permuted rows. All 985 project tests passed. The full-data eigensolve
took 3.415 s with maximum residual 1.81e-14; Dask was not installed.

The unchanged `tbs_family_diffusion` preset completed numerical execution in
129.035 s on all cells, with NNLS branch times applied. It returned one cluster
(ARI/NMI 0 against 190 cell-line labels). Root `N97756` separates one cell from
48,878; its child edge Tree-BH p-values are 0.9167335 and 0.7245607, exceeding
0.001. Descendant edge tests remain ancestor-blocked. Published `CHOIR_1`
assignments, rescored without rerunning CHOIR, have 367 clusters and ARI 0.664014.

`cancer_cell_lines_tbs.h5ad` and `cluster_assignments.csv` in the local run folder
contain the result. Full matrix and metadata export checks passed. The primary
run report and verification JSON record source/output hashes, solver limits,
and termination of a stalled temporary diagnostic cleanup after successful
numerical output verification. The one-cluster result is not a biological
homogeneity finding, and no parameters were tuned to force additional clusters.

### Why the one-cluster result occurs

Follow-up diagnosis confirms that the root tests a singleton ASPC1 cell against
the other 48,878 observations. All 88 of its kernel neighbors share its cell-line
label, and 987,220 of 1,030,508 stored kernel edges connect matching cell lines.
The singleton fails the eight-leaf within-child covariance guard, so the root
uses parent-total covariance. Its test sample size remains one despite the
large total dataset. Tree-BH schedules no descendant families after both root
edges fail, and pass-through does not override those edge prerequisites.

A frozen-tree replay reproduces the original statistics. Removing branch-length
variance alone changes both adjusted root p-values to 0.120680, still above
0.001; zero edges are significant. Thus that scaling alone cannot fix the
collapse. A known-K global height cut at 190 clusters has ARI 0.076053, but is
neither TBS nor an upper bound on adaptive tree cuts. The ranking benchmark's
122 cases have only 30–1,000 observations. Local graph structure, tree structure,
hierarchical gating, and benchmark transfer must be distinguished.

The primary diagnosis report and its two JSON captures document these checks,
the discarded initial replay with a mismatched basis dimension, and unresolved
geometry and calibration questions. No production behavior or original result
was changed during diagnosis.

### Recovery comparison and independent validation

The frozen symmetric/average, symmetric/Ward, Markov/average and Markov/Ward
comparison retained all cells and the same PCA, kernel and remaining settings.
All three new NNLS fits failed to converge within the fixed budget and were
refused by the native runner. They provide no new TBS assignments. Their root
child sizes are 250/48,629, 151/48,728 and 267/48,612 respectively, removing the
singleton without establishing successful recovery.

Independent count simulations rebuilt the complete pipeline on 600 homogeneous
nulls and 35 signal mixtures, with all four variants (2,540 successful method
evaluations). On the high-library-depth-variation null, false final partitions
occurred in 262/300, 228/300, 298/300 and 300/300 datasets respectively. No variant
met the prespecified P(K > 1) target 0.01 across both null models. Balanced
mixture recovery improved with symmetric/Ward but its false partitions preclude
adoption from this evidence. Null sample size was 400; signal sizes were
1,000–3,000. No 48,879-cell calibration or descendant-FDR guarantee follows.

The primary recovery report, frozen protocol and tables document model details,
separate one-sided 95% exact binomial bounds, complete fingerprints, failed
real-data fits and verification. The local comparison outputs are under
`results/analyses/choir_tbs_recovery_comparison_20261103_174622/`. Production
defaults and dependencies were unchanged; temporary distance files and workers
were cleaned up. [[choir-cell-lines-recovery-plan]] retains the open questions.

- `raw/inbox/choir-cell-lines/README.md` documents sources, selection, loading,
  identifier mapping, and the distinction between CPM and raw counts.
- `raw/inbox/choir-cell-lines/provenance.json` captures source hashes and extraction.
- `raw/inbox/choir-cell-lines/pool-count-validation.json` records count agreement.
- `raw/inbox/choir-cell-lines/validation.json` records the 181,398,276 nonzero
  float64 values, full matrix round-trip equality, output size, and SHA-256.

## Links

- [[choir-cell-lines-recovery-plan]]
- [[choir-singleton-peel-rerun]]
- [[index]]
- [[project-overview]]
