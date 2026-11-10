# The ceiling rule against the literature (lookups of 2026-11-09)

Lookups made for the comparison of the single-population split bound (ceiling rule) and its results with
published work. Each record was confirmed through the Crossref API (`api.crossref.org/works/<DOI>`) unless
marked otherwise. "Read" says how much of each source was seen.

## Records

1. Kimes, Liu, Hayes, Marron 2017, Statistical Significance for Hierarchical Clustering, Biometrics
   73(3):811-821, doi:10.1111/biom.12647. Read: author manuscript, PMC5708128
   (https://pmc.ncbi.nlm.nih.gov/articles/PMC5708128), via a summarising fetch.
   - Null at each node: the node's cells come from a single Gaussian; covariance from a factor model
     (hard, soft or sample eigenvalues).
   - Statistic: the 2-means cluster index CI = (SS1 + SS2) / TSS, or the linkage value.
   - p-value by Monte Carlo ("a closed form of the distribution of CIs under the null is unavailable").
   - Family-wise error by Meinshausen's procedure, cutoff alpha (N_j - 1)/(N - 1), testing down from the root.
   - Stated limits: anti-conservative for heavy tails (t3), "particularly pronounced for large sample sizes";
     conservative for Gaussian data; theory and simulations mainly for Ward's linkage. No closed form or 2/pi
     value for a single Gaussian is given.
2. Liu, Hayes, Nobel, Marron 2008, SigClust, J. Am. Stat. Assoc. 103(483):1281-1293,
   doi:10.1198/016214508000000454 (Crossref record only, 2026-11-08). Described in record 1 as the single-node
   test that SHC applies down the tree.
3. Gao, Bien, Witten 2024 (online 2022), Selective Inference for Hierarchical Clustering, J. Am. Stat. Assoc.
   119(545):332-342, doi:10.1080/01621459.2022.2116331; arXiv:2012.02936v3, pages 1-5 read.
   - Model X ~ MN(mu, I_n, sigma^2 I_q) with sigma known (unknown sigma in their Section 4.3).
   - Null: equal means of two estimated clusters. The Wald test with ||mean difference|| ~
     sigma sqrt(1/|C1| + 1/|C2|) chi_q is "extremely anti-conservative" after clustering; sample splitting
     does not fix it.
   - Their p-value conditions on the event that the two clusters are produced by the clustering; exact for
     agglomerative clustering with common linkages.
4. Keefe, Marron 2023, Powerful Significance Testing for Unbalanced Clusters, arXiv:2308.13079 (abstract page
   read; no journal reference listed). SigClust "is severely underpowered" when candidate clusters have
   unbalanced sizes; they propose a k-means generalisation.
5. Engelman, Hartigan 1969, Percentage Points of a Test for Clusters, J. Am. Stat. Assoc. 64:1647-1648,
   doi:10.1080/01621459.1969.10501084 (Crossref search). Not read. Secondary descriptions (search results):
   best split of the observations into two clusters by the ratio of between- to within-cluster sum of
   squares, percentage points under a single normal population, univariate. Table values unverified.
6. Hartigan 1978, Asymptotic Distributions for Clustering Criteria, Ann. Statist. 6(1):117-131,
   doi:10.1214/aos/1176344071. Not read. Search result: asymptotic distributions of the optimal split criterion,
   used for significance tests of clusters and bimodality; a later paper (arXiv:1205.2123) reproduces the
   N(0,1) case with the split at the median.
7. Hartigan, Hartigan 1985, The Dip Test of Unimodality, Ann. Statist. 13(1):70-84, doi:10.1214/aos/1176346577.
   Abstract via search: the uniform is the asymptotically least favourable unimodal distribution for the dip.
8. Schilling, Watkins, Watkins 2002, Is Human Height Bimodal?, Am. Stat. 56(3):223-229, doi:10.1198/00031300265.
   Abstract via search: an equal-weight mixture of two normals with common sigma is bimodal if and only if the
   means differ by more than 2 sigma.
9. Cover 1965, doi:10.1109/PGEC.1965.264137 (Crossref, 2026-11-08). Not read. Wikipedia "Cover's theorem":
   C(N, d) = 2 sum_{k=0}^{d-1} binom(N-1, k) homogeneously separable dichotomies of N points in general
   position in d dimensions.
10. Sant, Mucke, Corces 2025, CHOIR improves significance-based detection of cell types and states from
    single-cell data, Nat. Genet. 57(5):1309-1319, doi:10.1038/s41588-025-02148-8. Not read beyond search
    results: random-forest classifiers and permutation tests down a hierarchical clustering tree; benchmark
    of 15 methods on 230 simulated and 5 real datasets. Its cell-line results were not looked up.

## Not found

No looked-up source gives g(pi) = r/(1 - r), the value 2/pi for the share of variance between the halves of a
median-cut Gaussian, a ceiling on rho for cuts of a unimodal population, or the combination of a flat-cut count
with a noncentral chi-square tail. Kimes et al. state that no closed form for the single-Gaussian cluster index
is available to them.
