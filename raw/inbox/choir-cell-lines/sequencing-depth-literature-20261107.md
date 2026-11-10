# Sequencing depth in the table-building step: literature lookup

Retrieved 2026-11-07. Peer-reviewed records were retrieved from PubMed
(`lookup_article_by_citation`, `get_article_metadata`); the two preprints from
the bioRxiv API (`get_preprint`). Only what the retrieved abstracts state is
recorded here; nothing is filled in from memory.

## Records

| Work | Venue | DOI | What the abstract states |
| --- | --- | --- | --- |
| Lun A (2018), Overcoming systematic errors caused by log-transformation of normalized single-cell RNA sequencing data | bioRxiv preprint, not peer reviewed | [10.1101/404962](https://doi.org/10.1101/404962) | The log-transformation of normalized values can introduce systematic errors when cells vary in sequencing coverage, giving spurious non-zero differences and artificial population structure in simulations; similar effects in real data; derives a larger pseudo-count that controls the transformation-induced error to a set threshold. |
| Hafemeister C, Satija R (2019), Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression | Genome Biol 20:296 | [10.1186/s13059-019-1874-1](https://doi.org/10.1186/s13059-019-1874-1) | The number of molecules detected per cell can confound biological heterogeneity; Pearson residuals of a regularized NB regression with sequencing depth as a covariate remove the technical influence. |
| Townes FW, Hicks SC, Aryee MJ, Irizarry RA (2019), Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model | Genome Biol 20:295 | [10.1186/s13059-019-1861-6](https://doi.org/10.1186/s13059-019-1861-6) | UMI counts follow multinomial sampling with no zero inflation in negative controls; log of counts per million and highly-variable-gene selection produce false variability in dimension reduction; proposes GLM-PCA and deviance-based selection. |
| Lause J, Berens P, Kobak D (2021), Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data | Genome Biol 22:258 | [10.1186/s13059-021-02451-7](https://doi.org/10.1186/s13059-021-02451-7) | A parsimonious NB model has an analytic solution equal to rank-one Poisson GLM-PCA; negative-control UMI data are close to Poisson with moderate overdispersion; analytic Pearson residuals outperform other methods for finding variable genes. |
| Choudhary S, Satija R (2022), Comparison and evaluation of statistical error models for scRNA-seq | Genome Biol 23:27 | [10.1186/s13059-021-02584-9](https://doi.org/10.1186/s13059-021-02584-9) | Across 59 datasets, Poisson fits sparse data; genes with enough depth show overdispersion; its degree varies across datasets, systems and gene abundance. |
| Ahlmann-Eltze C, Huber W (2023), Comparison of transformations for single-cell RNA-seq data | Nat Methods 20:665-672 | [10.1038/s41592-023-01814-1](https://doi.org/10.1038/s41592-023-01814-1) | In simulated and real benchmarks, the logarithm with a pseudo-count followed by PCA performs as well as or better than the delta-method, residual, latent-state and factor-analysis alternatives. |
| Bacher R et al. (2017), SCnorm: robust normalization of single-cell RNA-seq data | Nat Methods 14:584-586 | [10.1038/nmeth.4263](https://doi.org/10.1038/nmeth.4263) | Assumptions of most normalization methods do not hold for single cells, and applying them introduces artifacts that bias downstream analyses. |
| Vallejos CA, Risso D, Scialdone A, Dudoit S, Marioni JC (2017), Normalizing single-cell RNA sequencing data: challenges and opportunities | Nat Methods 14:565-571 | [10.1038/nmeth.4292](https://doi.org/10.1038/nmeth.4292) | Commonly used normalization approaches can produce misleading results; gives recommendations. |
| Booeshaghi AS, Hallgrimsdottir IB, Galvez-Merchan A, Pachter L (2022), Depth normalization for single-cell genomics count data | bioRxiv preprint, not peer reviewed | [10.1101/2022.05.06.490859](https://doi.org/10.1101/2022.05.06.490859) | A variance-stabilizing transform for NB counts followed by a proportional fitting step (PFlog) is the only feature-relabeling-equivariant method meeting their three requirements; equivalent to a shifted centred-log-ratio transform. |

Not verified: the exact pseudo-count formula of Lun (2018), whether it
recommends downsampling, and the exact order of steps in PFlog. These need the
full texts.

## Relation to the project's measurements (2026-11-07)

The project's table step is per-cell scaling to 10,000, `log1p`, 2,000
highly variable genes, gene scaling and 30 PCs, the setting Lun (2018) and
Townes et al. (2019) describe. Depth-confound runs
(`results/analyses/depth_confound_20261107`, default tree, 300 no-group
datasets each) give 2, 50, 238 and 259 false partitions at within-group
log-depth sd 0.25, 0.42, 0.61 and 1.0. Equal-depth downsampling gives 5, 11,
18 and 12; regressing PCs on log depth gives 2, 16, 11 and 231; analytic
Pearson residuals give 1, 1, 7 and 33, with signal ARI 0.655 against 0.859.
The rise with depth spread, the removal by equal depth and the failure of a
linear depth correction at sd 1.0 are consistent with the log-step artefact
described by Lun (2018); this is not a proof of that cause.
