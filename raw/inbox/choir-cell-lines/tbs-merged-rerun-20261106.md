# Cancer cell-line TBS rerun after the singleton-peel merge

Run on 2026-11-06 after merging branch commit `cbbdf1fd` ("fix: keep singleton
peels from ending the edge descent") into `dev` as merge commit `50455316`.

## What changed and what was held fixed

The rerun reused the 2026-11-03 prepared inputs unchanged: `pca.npy`,
`diffusion_coordinates.npy`, `cell_metadata.csv`, `preparation.json` and
`diffusion-verification.json` were copied from
`results/analyses/choir_cell_lines_tbs_20261103_163359/`. The run script is
byte-identical to `tbs-run-script.py.txt`. The 2026-11-03 preparation already
zeroed the scaled-distance diagonal before epsilon and solved the kernel per
connected component, so the cached coordinates match the committed diffusion
code.

The resulting linkage is byte-identical to the 2026-11-03 linkage, with the
same 1 versus 48,878 root (`L43204` | `N97755`). The only difference is the
merged gate and traversal code: Tree-BH families merged through singleton
peels, edge-gate and sibling-frontier pass-through, and the sibling inflation
cap at the family baseline beyond the largest calibration parent.

The working tree held uncommitted benchmark-registration, documentation and
wiki edits; none is on the `run_tbs.py` import path. Versions: SciPy 1.16.3,
NumPy 2.3.4, pandas 2.3.3, scikit-learn 1.7.2.

## Results

| Quantity | 2026-11-03 | 2026-11-06 merged |
| --- | ---: | ---: |
| Clusters | 1 | 70 |
| Live nodes visited | 1 | 145 |
| Splits / pass-throughs / boundaries | 0 / 0 / 1 | 45 / 27 / 70 |
| Significant edges after Tree-BH | 0 | 148 |
| Root edge BH p, `L43204` / `N97755` | 0.916734 / 0.724561 | 0.916734 / 0.434736 |
| ARI against 190 cell lines | 0.000000 | 0.027593 |
| NMI against 190 cell lines | 0.000000 | 0.548219 |
| ARI against published `CHOIR_1` | 0.000000 | 0.013880 |

Both root edges remain non-significant; traversal now passes through the
singleton root instead of stopping. Published `CHOIR_1` labels score ARI
0.664014 against the cell lines (not rerun).

The largest cluster holds 28,611 cells (58.5%) from 144 cell lines; 120 of
those lines have at least 90% of their cells in it. There are four singleton
clusters and the median cluster size is 175. Of 190 cell lines, 27 are
recovered under this definition: at least 90% of the line's cells lie in its
dominant cluster, and that cluster is at least 90% that line. Of 70 clusters,
46 are at least 90% one cell line. These counts were computed from
`cluster_assignments.csv` with `pandas.crosstab`.

Numerical TBS execution took 145.5 s. As on 2026-11-03, the process did not
exit after writing its outputs and removing the temporary distance file; it
was terminated after 34 minutes, once the outputs were checked (48,879
assignment rows). `runtime.json` still reads `running`.

## Scope

The opened tree is not evidence of calibrated selected-tree inference. The
same production pipeline produced false final partitions on many independent
count nulls (`tbs-recovery-comparison.md`); this rerun did not repeat those
nulls. Passed-through peel leaves whose edges are not significant join the
nearest cluster below the peel, so the result is not an exact tree cut.

## Files

- `tbs-merged-rerun-summary-20261106.json`: the run's `summary.json`.
- `tbs-merged-rerun-composition-20261106.csv`: per-cluster composition.
- `tbs-merged-rerun-fragmentation-20261106.csv`: per-cell-line fragmentation.
- Local-only outputs: `results/analyses/choir_cell_lines_tbs_merged_20261106/`.

SHA-256:

- `pca.npy` a2ef1e696acada626d996e7fda7e40522054cd3ebb1a855a5cb83be58a984283
- `diffusion_coordinates.npy` c6844134eb7b8257650b67c2108d7ed7d06aeab2e7f8b19f98a5d5a88e18c79e
- `cell_metadata.csv` 744d353a2aa4308b04e67b9c4e8e487d65084900791c1b6b0faf9d7ad705c9ad
- `linkage.npy` b1c6e8811cc9f1251135c1f20243c8177ccc4909b18c0372379ee7f1cb86da36
- `summary.json` 1ec2846763ca77df4c8ef1cf825cf4b14456a8f1765ee0a33ac5de4e4efcf0df
- `cluster_assignments.csv` 7033e8d16f76a83fb1501257425d50f8717c2786c6f42d101e8f08fb63eb631e
- `annotations.csv` 5af9188050071d8ce7a7fb26145c739914457fbd0b5f3d174f937fc3e3e2654a
