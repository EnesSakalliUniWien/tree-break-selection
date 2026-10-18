# Lower-dimensional reference-law sources checked on 2026-09-28

These are paraphrased source notes, not a production calibration proposal.

- Richard Lockhart, STAT 350 General Theory, quadratic-form lecture notes:
  https://www.sfu.ca/~lockhart/richard/350/08_2/lectures/Theory/slides.pdf
  A fixed symmetric quadratic operator acting on a standard Gaussian vector
  gives a weighted sum of independent squared standard normals. Idempotence
  yields a chi-square law, with degrees of freedom equal to the operator rank.
  The fixed-operator requirement matters for data-selected projections.
- Ryan Tibshirani, Recent Advances in Selective Inference (2016):
  https://www.stat.cmu.edu/~ryantibs/talks/ssc-2016.pdf
  The forward-stepwise example illustrates how selecting a single variable
  before applying an ordinary one-degree-of-freedom test invalidates that
  reference. The presentation develops conditioning on the selection event.
- Lucy L. Gao, Jacob Bien and Daniela Witten, Selective Inference for
  Hierarchical Clustering, accepted version (2022):
  https://arxiv.org/abs/2012.02936v3
  The abstract distinguishes groups defined in advance from groups defined
  by clustering. It presents inference accounting for selected hypotheses,
  and warns that separate training/testing datasets alone do not universally
  resolve the clustering inference problem. This is not a theorem validating
  this repository's diffusion, adaptive projection or empirical inflation.

Access: the arXiv abstract and both authors' university-hosted lecture PDFs
were read. The PMC article page returned a browser challenge. The attempted
1999 American Statistician article endpoint returned HTTP 403; no full-text
claim relies on that inaccessible page.
