# TBS multiplicity audit: primary sources

Accessed 2026-11-04. Source-grounded notes; the implementation comparison is our inference.

- [Benjamini and Hochberg (1995)](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.2517-6161.1995.tb02031.x)
  proves FDR control for its procedure under independent test statistics.
  Applying its arithmetic does not make invalid input p-values valid.
- [Bogomolov, Peterson, Benjamini and Sabatti (2021), Hypotheses on a tree](https://pmc.ncbi.nlm.nih.gov/articles/PMC9945647/)
  defines level-specific selective error measures. Section 4's TreeBH family
  threshold includes ancestor rejection proportions. Its results state
  dependence, parent-combination and hypothesis-structure conditions; Theorem 1
  uses PRDS, intersection hypotheses and Simes combination. Other results have
  their own conditions. The TBS mean-contrast tree is not automatically such an
  intersection hierarchy, and its directly computed parent-edge p-values are
  not Simes combinations of descendant p-values. This comparison does not prove
  that no other guarantee is possible; it identifies missing premises for
  transferring the cited theorem to the implemented TBS pipeline.
