---
title: Sequencing Depth Normalization Literature
type: source
status: draft
updated: 2026-11-07
sources:
  - raw/inbox/choir-cell-lines/sequencing-depth-literature-20261107.md
tags:
  - scrna
  - normalization
  - sequencing-depth
  - literature
---

# Sequencing Depth Normalization Literature

## Summary

Nine works, checked in PubMed or bioRxiv on 2026-11-07, on how per-cell
sequencing depth enters the expression table. Lun (2018, preprint) reports that
taking the log of depth-normalized counts creates artificial population
structure in simulations when cells differ in depth, which is the table step
TBS currently uses. The field's other answers are count models (NB or Poisson
residuals, GLM-PCA) or changes to the log step.

## Key Points

- Lun (2018, bioRxiv, not peer reviewed): the log transform of normalized
  values gives spurious differences and artificial structure when coverage
  varies; derives a larger pseudo-count that bounds the error.
  [10.1101/404962](https://doi.org/10.1101/404962)
- Townes et al. (2019): log CPM plus highly-variable-gene selection produces
  false variability in dimension reduction; UMI counts follow multinomial
  sampling. [10.1186/s13059-019-1861-6](https://doi.org/10.1186/s13059-019-1861-6)
- Hafemeister and Satija (2019) and Lause et al. (2021): depth as a covariate
  in an NB or Poisson model, using Pearson residuals; Lause et al. find
  negative-control data close to Poisson.
  [10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1),
  [10.1186/s13059-021-02451-7](https://doi.org/10.1186/s13059-021-02451-7)
- Ahlmann-Eltze and Huber (2023): in benchmarks, log with a pseudo-count plus
  PCA performs as well as or better than the alternatives.
  [10.1038/s41592-023-01814-1](https://doi.org/10.1038/s41592-023-01814-1)
- Choudhary and Satija (2022), Bacher et al. (2017), Vallejos et al. (2017)
  and Booeshaghi et al. (2022, preprint) are listed with DOIs in the raw note.
- Project measurements: default-tree false partitions rise with within-group
  depth spread (2, 50, 238, 259 of 300 at log-depth sd 0.25, 0.42, 0.61, 1.0).
  Equal-depth downsampling removes most of them; a linear depth regression does
  not at sd 1.0 (231). This is consistent with the log-step artefact, not a
  proof of it.
- Not verified: Lun's pseudo-count formula and recommendations, and the PFlog
  step order; both need the full texts.

## Evidence

`raw/inbox/choir-cell-lines/sequencing-depth-literature-20261107.md` holds the
records, the retrieval route and the comparison with the local depth-confound
runs (`results/analyses/depth_confound_20261107`, git-ignored).

## Links

- [[choir-cell-lines-recovery-plan]]
- [[tree-break-selection]]
