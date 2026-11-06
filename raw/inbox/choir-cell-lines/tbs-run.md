# Full cancer cell-line family-diffusion TBS run

Run on 2026-11-03 after the user selected the highest performing available
diffusion TBS method and suggested a library solution for the full dataset.

## Method and execution

The latest full benchmark selects `tbs_family_diffusion` by mean ARI
(0.825843 on 122 supported cases). Its ARI equals adaptive-diffusion NNLS on
their 101 common supported cases. This is a selection from prior benchmarks,
not evidence that the method performs best on this dataset.

All 48,879 cells were retained. Preparation used `log1p(CPM / 100)` without
library-total renormalization, 2,000 Seurat-flavor variable genes, scaling clipped
at 10, 30 PCs, and seed 0. The original CPM input remains unchanged.

The selected continuous preset uses 15-neighbor Euclidean variable-bandwidth
diffusion, median epsilon, diffusion time 3, 30 effective diffusion coordinates,
average linkage, fixed-topology NNLS branch times targeting mean squared
standardized Euclidean distance in the original PCA geometry, and normalized
branch-length variance. The registered 50,000-pair sample, seed 0, tolerance
1e-5 and maximum 1,000 NNLS iterations were retained. Coverage augmentation
produced 56,414 fitted pairs. Edge alpha is 0.001, sibling alpha is 0.01,
projection energy fraction is 0.90, and pass-through is enabled.

The three-component kernel was handled with installed SciPy: solve each
connected component, merge eigenvalues globally, and retain complete numerical
tie groups. Large blocks use `eigsh`; small blocks use `eigh`. This avoids a
full dense graph allocation. The actual coordinate computation took 3.415 s,
with maximum eigen residual 1.81e-14 and orthogonality error 1.55e-14. The
existing ARPACK limitation for repeated nontrivial eigenvalues inside a large
connected block remains documented. Dask installation was authorized but was
not needed or performed. No dependency or lockfile changes were made.

## Results

| Assignment | Cells | Clusters | ARI against cell-line labels |
| --- | ---: | ---: | ---: |
| TBS family diffusion | 48,879 | 1 | 0.000000 |
| Published CHOIR_1 labels, rescored | 48,879 | 367 | 0.664014 |

CHOIR was not rerun. TBS NMI against the 190 cell-line labels is also zero.
The one-cluster result does not recover the known cell lines and is not evidence
that they are biologically homogeneous.

The final root split is one cell (`L43204`) versus 48,878 (`N97755`), under
root `N97756`. Root-child edge p-values after Tree-BH are 0.9167335 and
0.7245607, both above 0.001. No edge is significant; 97,754 descendant nodes
are ancestor-blocked. Traversal visits only the root and returns one boundary.
This describes the gate decision; it does not establish why the statistical
reference or branch-time model fits this dataset poorly.

The native runner reports status `ok`; NNLS branch times were applied after
796 iterations. SciPy stopped on relative cost change, with reported optimality
2.8945; no stronger optimality claim is made. Numerical TBS execution took
129.035 s, including distance computation. No parameters were tuned to the
published cell-line or CHOIR labels.

## Outputs and verification

Local-only artifacts are under
`results/analyses/choir_cell_lines_tbs_20261103_163359/`:

- `cancer_cell_lines_tbs.h5ad`: original CPM matrix/metadata plus `TBS_cluster`,
  `X_pca`, and `X_diffusion`.
- `cluster_assignments.csv`, `cluster_composition.csv`, and
  `cell_line_fragmentation.csv`.
- `linkage.npy`, `annotations.csv`, `summary.json`, and preparation checkpoints.
- `prepare.py`, `coordinates.py`, `run_tbs.py`, and `export.py` preserve the
  execution steps. The 8.9 GiB temporary distance array was removed.

Every original CPM entry, cell/gene ID and metadata column was compared with
the annotated AnnData export. All assignments, output metrics, method preset
and root-gate decisions were checked independently. Hashes and compact results
are in `tbs-run-verification.json`, `tbs-summary.json`, and
`tbs-diffusion-verification.json`. Source/script snapshots are captured beside
this report.

All 985 tests passed through `make test`, including six new cases comparing
disconnected sparse and dense results, global ties, an isolated vertex and
permuted rows. One existing sklearn MDS default-value FutureWarning was emitted.
Focused Ruff checks and formatting also passed.

The temporary faulthandler diagnostic logging added for execution monitoring
stalled during cleanup after the numerical summary had been written. Once the
annotated export and all numerical outputs were independently verified, only
that completed run process was terminated (exit 143). No TBS numerical error
was reported; this was not an uninterrupted process exit. `runtime.json` records
the cleanup, and no execution process remains running.
