# Which gene programs carry the directions inside a cell line

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/line_gene_programs_20261112` (git-ignored);
protocol, deviations, the gene list, both scripts and their outputs are in
`line-gene-programs-scripts-20261112.py.txt`. Follow-up of `line-non-sphericity-20261112.md`. No production code
was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, anndata, uv 0.9.5.

## Looked-up sources (retrieved 2026-11-12)

1. G. S. Kinker et al., "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular
   heterogeneity", Nature Genetics 52(11), 1208-1218, 2020. https://doi.org/10.1038/s41588-020-00726-6
   (PubMed record 33128048, abstract). Confirmed: the paper of this dataset profiles 198 cancer cell lines and
   reports 12 expression programs that are recurrently heterogeneous within multiple cell lines, associated
   with cell cycle, senescence, stress and interferon responses, epithelial-mesenchymal transition and protein
   metabolism.
2. I. Tirosh et al., "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
   Science 352(6282), 189-196, 2016. https://doi.org/10.1126/science.aad0501 (PubMed record 27124452).
   Confirmed from the abstract: heterogeneity of malignant cells associated with the cell cycle. The scanpy
   documentation (https://scanpy.readthedocs.io/en/1.11.x/how-to/cell-cycle.html) names this group's work as the
   source of a cell-cycle list of 97 genes, 43 for S phase and 54 for G2M.
3. Gene symbols: the original list file could not be retrieved (the documented GitHub paths return "not
   found"). The symbols used are a secondary copy from
   https://natian-and-ryabhatta.web.app/cell-cyle-scores.html, which attributes its lists to the Seurat package:
   42 S phase and 51 G2M symbols as returned by the lookup, 4 fewer than the documented 97. 91 of the 93 are
   gene names of the data (not found: UHRF1, CASP8AP2). Unverified: that this copy equals the original list.

## Setup

175 cell lines with at least 100 cells. Expression log1p(CPM / 100). Per line: the three leading eigenvectors of
the line's covariance in the 30 PCs (directions 1 to 3) and the cells' scores; per gene detected in at least
10% of the line's cells its correlation with each score. Cell-cycle scores per cell: the mean of the
within-line standardized expression over the S list and over the G2M list. The protocol with three predictions
was written before the run; the self-checks passed (median T1 25.64; cell order equal; 91 symbols found).

## Results

Genes of the directions without any gene list. A gene drawn at random is among the 50 genes of a direction in
1.04 of 175 lines on average. Most recurrent genes (lines; * = in the cell-cycle list):

- direction 1: RPL29 45, RPL41 45, CDC20* 42, CENPF* 41, CCNB1 41, UBE2S 40, RPLP1 39, RPS2 39, RPLP0 36,
  HSP90AA1 36, AURKA* 35, S100A6 35, UBC 34, HLA-A 33, RPS14 33, BIRC5* 33; 9 of the top 40 are in the list.
- direction 2: AURKA* 84, CCNB1 82, CENPF* 78, KPNA2 72, TPX2* 72, CDC20* 71, UBE2C* 70, ARL6IP1 69, ASPM 66,
  PLK1 66, CKS2* 65, TOP2A* 65, TUBB4B* 61, MKI67* 60; 21 of the top 40 are in the list.
- direction 3: CENPF* 98, TPX2* 97, AURKA* 96, TOP2A* 91, UBE2C* 81, MKI67* 76, CCNB1 74, ARL6IP1 74; 22 of
  the top 40 are in the list.
- any of the three: AURKA* 159, CENPF* 154, TPX2* 153, CCNB1 150, TOP2A* 145, CDC20* 144, PLK1 140; 303 genes
  are among the 50 in at least 20 lines.

Cell cycle:

| | Median over 175 lines (min, max) |
| --- | --- |
| Share of the within-line variance explained by the S and G2M scores | 13.6% (4.8%, 40.2%) |
| Control: 20 random pairs of gene sets of the same sizes, per-line median share | 3.5% (1.0%, 20.1%) |
| Largest correlation of directions 1 to 3 with the S or G2M score | 0.679 (0.219, 0.934) |
| T1 in PC coordinates | 25.64 (7.03, 94.13) |
| T1 after removing the cell-cycle fit | 24.29 (4.95, 87.54) |

One of the directions 1 to 3 correlates with a cell-cycle score at 0.5 or more in 155 of 175 lines: direction
1 in 34 lines, direction 2 in 72, direction 3 in 74. Of the 50 genes of direction 1 the median number in the
cell-cycle list is 0 (at least 10 in 26 lines). After removing the cell-cycle fit, k of the current rule is 4
to 12 (5 to 12 before) and k of the prototype rule 7 to 27 (6 to 27 before).

Direction 1 after removing the cell-cycle fit, most recurrent genes: RPL41 44, RPS2 42, RPLP0 42, RPL29 41, B2M
41, ANXA2 40, HLA-C 40, S100A11 39, HLA-A 38, RPLP1 36, HLA-B 36, UBC 36, SLC3A2 36, PRDX1 35, TMSB10 35, FXYD3
35, S100A6 35; none of the top 40 is in the cell-cycle list. These genes are reported by symbol; no looked-up
gene list assigns them to a program here.

Predictions: P1 held (13.6% against the 10% stated); P2 held (155 of 175 lines); P3 held (median T1 24.29
after removing the cell cycle).

## Added after the result (`after.py`)

The 86 genes whose symbols are RPL or RPS with a number (and RPLP0 to RPLP2, RPSA), scored like the cell-cycle
lists, and all measured factors together:

| Factor removed | Share of the within-line variance explained, median (min, max) | T1 afterwards, median |
| --- | --- | ---: |
| Cell cycle (S and G2M scores) | 13.6% (4.8%, 40.2%) | 24.29 |
| RPL/RPS score | 6.5% (1.5%, 43.0%) | 23.49 |
| Log depth | 3.0% (0.7%, 8.6%) | 25.41 |
| All three together | 23.6% (11.2%, 57.9%) | 19.44 |

Direction 1 correlates with the RPL/RPS score at 0.5 or more in 37 of 175 lines (median 0.286). The RPL/RPS
score correlates with log depth at 0.219 at the median (-0.147 to 0.599).

## Reading

- The cell cycle is the most recurrent content of the leading directions inside a line: its genes are among
  the 50 of a direction in up to 159 of 175 lines against 1.04 expected, and it explains 13.6% of the
  within-line variance, four times a random gene set and more than depth (3.0%).
- It sits mostly in directions 2 and 3. Direction 1, the largest, follows the cell cycle in 34 lines; its most
  recurrent genes are RPL and RPS genes and HLA-A, HLA-B, HLA-C, B2M.
- Removing the cell cycle, the RPL/RPS score and depth together removes 23.6% of the within-line variance and
  leaves T1 at 19.4. The lines stay far from spherical.
- This agrees with the dataset's own paper (source 1): cell lines are heterogeneous within themselves along
  several recurring programs, the cell cycle among them. The within-line structure is biology of the cells,
  not noise.
- For step A2 of the plan: inside one population there are real directions of variation. A count of
  directions with excess variance cannot be zero for one population on such data, whatever the edge.

## Not checked

The cell-cycle symbols are a secondary copy (93 of the documented 97; 91 in the data). No gene lists for the
other programs of source 1 were looked up, so the RPL/RPS and HLA genes carry no program label here. Linear
fits only. The RPL/RPS set was chosen after seeing the result. One dataset; lines with at least 100 cells.
