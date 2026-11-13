---
title: What Makes the Cells of One Cell Line Non-Spherical
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/line-non-sphericity-20261112.md
  - raw/inbox/choir-cell-lines/line-non-sphericity-scripts-20261112.py.txt
  - raw/inbox/choir-cell-lines/line-gene-programs-20261112.md
  - raw/inbox/choir-cell-lines/line-gene-programs-scripts-20261112.py.txt
tags:
  - choir
  - real-data
  - invariant-directions
  - covariance
  - sequencing-depth
---

# What Makes the Cells of One Cell Line Non-Spherical

## Summary

Inside one CHOIR cell line the largest eigenvalue is 25.6 times the median
eigenvalue (median of 175 lines). The cause is specific to each line: a
covariance common to all lines brings the ratio to 13.8, sequencing depth
explains 3.0% of the within-line variance and pools 0.8%. A common covariance
estimated without labels from small subtrees acts like the label-based one.
After either whitening every line and every called node of at least 40 cells
still has directions above a spherical edge. The leading directions inside a
line are carried by genes: the cell cycle explains 13.6% of the within-line
variance, and cell cycle, RPL/RPS genes and depth together 23.6%.

## Key Points

- Whitened by the pooled within-line covariance of the other lines: ratio
  13.8 (3.7 to 117.1); above the spherical 99% point in 175 of 175 lines;
  at least 2 directions under the prototype rule in 175 of 175.
- Depth fit removed: ratio 25.4. Depth explains 0.7% to 8.6% of the
  within-line variance; `source_library` means 0.1% to 11.4%; no line spans
  more than one `pool_id`.
- In PC coordinates 65.8% of the energy of a line's covariance is off the
  diagonal at the median; the largest median within-line variance is at PC18
  (3.43, 45% of that PC's total variance).
- Label-free covariance from maximal subtrees of at most 50 cells (1,663
  subtrees, 84.3% of them at least 90% one line): against the label-based
  one its ratio per direction is 0.260 to 1.027, median 0.808. Lines
  whitened by it: ratio 14.5.
- Called nodes whitened by it: k above 15 under the prototype rule at 290 of
  1,838 nodes (417 in PC coordinates); 18 to 28 at all 103 nodes with at
  least 1,000 cells.
- Nodes with at least 90% of their cells from one line (1,017) and all other
  nodes (342) have a direction in every case; median k prototype 9 against
  17 whitened, 12 against 18 in PC coordinates.
- Predictions: 2 of 6 held. "The common covariance carries most of it" and
  all three label-free predictions did not hold.
- Independent recomputation with separate code reproduces the ratios 13.774
  and 14.5, the depth share 0.030 and the eigenvalues 0.260, 0.808, 1.027.
- Genes of the directions (follow-up): cell-cycle genes are among the 50
  genes of one of the three leading directions in up to 159 of 175 lines
  (AURKA 159, CENPF 154, TPX2 153, TOP2A 145) against 1.04 expected.
- The S and G2M scores explain 13.6% of the within-line variance (4.8% to
  40.2%; random gene sets 3.5%); one of directions 1 to 3 correlates with
  them at 0.5 or more in 155 of 175 lines, direction 1 in 34, direction 2
  in 72, direction 3 in 74.
- Direction 1 is carried mostly by other genes: RPL and RPS genes, HLA-A,
  HLA-B, HLA-C, B2M. An RPL/RPS score explains 6.5%; cell cycle, RPL/RPS and
  depth together 23.6%, and the eigenvalue ratio is still 19.4 afterwards.
- According to PubMed, the paper of this dataset (Kinker et al. 2020,
  https://doi.org/10.1038/s41588-020-00726-6) reports 12 programs that are
  recurrently heterogeneous within cell lines, the cell cycle among them.
  The within-line structure is biology of the cells, not noise.
- Limits: the cell-cycle symbols are a secondary copy (93 of the documented
  97); no gene lists for other programs; the RPL/RPS set was chosen after
  seeing the result; no threshold on the count was tried; one dataset and
  one tree; lines with at least 100 cells.

## Evidence

`raw/inbox/choir-cell-lines/line-non-sphericity-20261112.md` holds the
question, both parts with their tables, the six predictions, the independent
spot checks, the reading, the deviations and the limits; the protocol, the
scripts and their outputs are in
`raw/inbox/choir-cell-lines/line-non-sphericity-scripts-20261112.py.txt`.
The genes of the directions, the looked-up sources and the cell-cycle
measurements are in
`raw/inbox/choir-cell-lines/line-gene-programs-20261112.md` with protocol,
gene list, scripts and outputs in
`raw/inbox/choir-cell-lines/line-gene-programs-scripts-20261112.py.txt`.

## Links

- [[k-margin-prototype]]
- [[non-heuristical-statistic-plan]]
- [[marchenko-pastur-multiplicity-literature]]
- [[choir-cell-lines-data]]
- [[sequencing-depth-normalization-literature]]
