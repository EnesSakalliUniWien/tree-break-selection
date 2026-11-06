# CHOIR pooled cancer cell-line data

Loaded on 2026-11-03 after the user selected the pooled cancer cell lines from
Figures 3–4 of [Sant et al. 2025](https://www.nature.com/articles/s41588-025-02148-8).

## Local data

The ready-to-load file is
`local_data/choir/kinker_2020/kinker_2020_choir_cells_cpm.h5ad`:

- 48,879 cells × 22,722 genes; 190 cell lines and 22 cancer types.
- Sparse CSR `float64` expression, with the published GEO CPM values preserved.
- All 109 columns from the authors' `Cell lines - Cell labels` sheet in `obs`,
  including `CellID`, `Cell_line`, `Cancer_type`, and the published clustering
  labels. `CHOIR_1` contains 367 published clusters; these were imported, not rerun.
- `geo_cell_id`, `source_library`, and `pool_id` record the identifier mapping.
  `obs_names` use the authors' `CellID`; `var_names` retain GEO gene names.
- 507,757,282 bytes on disk. The bulk data and the local conversion script
  `local_data/choir/kinker_2020/load_geo.py` are Git-ignored, local-only files.

From the existing project environment:

```python
import anndata as ad

adata = ad.read_h5ad(
    "/Users/berksakalli/Projects/kl-te-cluster/local_data/choir/kinker_2020/kinker_2020_choir_cells_cpm.h5ad",
    backed="r",
)
# adata.obs["Cell_line"], adata.obs["Cancer_type"], adata.obs["CHOIR_1"]
# adata.file.close() when finished using backed mode.
```

## Sources and selection

- [GEO GSE157220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157220):
  original `GSE157220_CPM_data.txt.gz` retained unchanged under `local_data/choir/kinker_2020/`.
  The complete download has 53,513 cells. Its series record, one representative
  sample record, and the pool-composition workbook are captured in this directory.
- [Zenodo 14641222](https://zenodo.org/records/14641222):
  `SourceData1_RealData.xlsx`, retained unchanged beside the GEO download.
  Extracted its `Cell lines - Cell labels` sheet, header row 18, to
  `local_data/choir/kinker_2020/choir_cell_labels.csv` using standard-library
  ZIP/XML parsing. Workbook shared strings and spreadsheet column positions
  were resolved explicitly; no dependency was installed.
- [Authors' analysis script](https://github.com/corceslab/2025_Sant_CHOIR/blob/3d057a88e87ef318dc5206b815023820dddf5954/additional_analysis_scripts/real_datasets/Kinker_2020_cancer_cell_lines.R):
  unchanged public source snapshot `Kinker_2020_cancer_cell_lines.R`.

The observed identifier relationship is
`Pool{library}_{pool}#{barcode}-1` → `{barcode}-{library}-{pool}`.
All 48,879 mapped IDs are unique and present in GEO. Every cell-line count also
matches the independent GEO pool-composition workbook. The remaining 4,634
GEO cells belong to eight cell lines in the custom head-and-neck pool, absent
from the CHOIR benchmark label sheet. They remain in the original download.

## Verification and limits

`provenance.json` records source URLs, byte sizes, SHA-256 hashes, the verified
Zenodo MD5, runtime versions, and extraction details. `pool-count-validation.json`
records the independent count comparison; `validation.json` records the output
hash, dimensions, and full matrix read-back equality.

All source values were checked for finite, nonnegative numeric values. The gzip
was read through its end, validating its CRC. All expression entries, cell IDs,
gene IDs, and metadata survived an AnnData write/read check without changes.

This loads the cells and published labels used for the CHOIR benchmark with
GEO's processed CPM expression. The authors' script uses separately
preprocessed raw counts; that count matrix is not reproduced here. No raw SRA
reads, alignment, normalization, log transformation, clustering, or figure
reproduction was run. CPM column totals are below one million in this deposited
gene matrix, as recorded in `validation.json`; no rescaling was applied.
