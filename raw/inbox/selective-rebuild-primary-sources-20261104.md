# Primary sources checked for the replacement comparison

Retrieved 2026-11-04. Paraphrases and applicability notes; no production change.

## Gao, Bien and Witten

Selective Inference for Hierarchical Clustering, JASA, DOI
https://doi.org/10.1080/01621459.2022.2116331 . Author manuscript:
https://arxiv.org/html/2012.02936 . Read Section 4.1, equations 22–23,
and the Gaussian conditional construction.

The paper explicitly identifies poor ordinary Monte Carlo sampling when the
observed contrast lies far in the chi tail. It proposes normal proposals near
the observed contrast and density-ratio weighting. The general event can be
checked by rerunning the clustering algorithm. Its exact polynomial event
calculations for standard linkages do not transfer to fitted diffusion.
Our additional chi/normal mixture and both event definitions are diagnostic
choices, not the paper's exact implementation.

## Kimes, Liu, Hayes and Marron

Statistical Significance for Hierarchical Clustering, Biometrics (2017), DOI
https://doi.org/10.1111/biom.12647 . Primary author preprint:
https://arxiv.org/html/1411.5259 . Read Sections 2.2 and 3.1.

SHC simulates a one-Gaussian-population reference, repeats hierarchical
clustering and compares separation scores. Its cluster index is within-group
sum of squares divided by total sum of squares. Covariance estimation is part
of the procedure and depends on assumptions. The present native-diffusion root
experiment uses the complementary score BSS/TSS; it is inspired by this idea,
not an SHC reproduction and not a transfer of its hierarchical error theorem.

## Chen and Gao, new preprint

Selective Inference for Clustering Under Unknown Covariance, submitted
2026-10-12: https://arxiv.org/html/2610.17458 . Read Sections 2.2, 3.1, 4.1–4.3
and discussion.

The preprint derives a conditional polar-factor test for unknown common
Gaussian covariance. The null requires all observations in the tested pair
to share a mean, stronger than equality of their average means. Its full-rank
polar construction requires r-1 >= p; the discussion leaves the weaker null,
dependent rows and non-Gaussian counts open. It is not a ready cancer-data fix.

The paper also uses a reversible-chain rank construction attributed to Besag
and Clifford (1989): insert the observation at a uniformly random path position,
simulate both directions, and rank it including ties. Our radial known-scale
control uses this construction on a different state space; it does not implement
the unknown-covariance test. The paper distinguishes exact chain ranks from
importance estimates and reports low effective sample sizes in selected tails.
