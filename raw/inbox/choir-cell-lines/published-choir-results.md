# Published CHOIR cancer cell-line results

Retrieved 2026-11-03 from the indexed full text of the published article:
Sant, Mucke and Corces, Nature Genetics 57, 1309–1319 (2025).
DOI: https://doi.org/10.1038/s41588-025-02148-8.
Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC12108546/.
Direct publisher access redirected to authentication and direct PMC access
showed a browser challenge; the indexed published PMC text supplied the
Results, Methods and figure captions. The earlier preprint was not used for
the following findings.

## Published findings, paraphrased

- Figure 3 uses 48,879 cells from 190 cancer cell lines. Among the six methods
  compared, only CHOIR distinguished all 190 lines under the evaluated settings.
  The alternatives were Cytocipher, GiniClust3, SCCAF, sc-SHC and Seurat.
- Figure 3j reports stronger recovery by CHOIR when each line was separately
  downsampled to 50 cells while the other lines remained present.
- Figure 4 reports significantly different proliferation module scores in
  74.3% of the lines that CHOIR divided into multiple clusters. A375 and T47D
  state differences were supported in independent datasets from those lines.
- The cancer benchmark uses entropy of cluster accuracy and a line-recovery
  criterion: the cluster containing the most cells from a line must have that
  line as its majority and cannot also be another line's primary cluster.
  Distinguishing every line therefore does not assert perfect per-cell purity
  or exactly 190 clusters. Within-line states are permitted.

## Existing local rescoring, distinct from the paper's reported metric

The imported default labels `CHOIR_1` contain 367 clusters. Existing local
rescoring gives ARI 0.6640144797441466 against the 190 line IDs, versus 0 for
the one-cluster TBS baseline. Sources: `README.md`, `tbs-run.md`, and
`tbs-summary.json` in this directory. This ARI was calculated locally, not
reported as the paper's cancer benchmark. CHOIR was not rerun. Our deposited
CPM-derived workflow differs from the authors' separately processed counts.
