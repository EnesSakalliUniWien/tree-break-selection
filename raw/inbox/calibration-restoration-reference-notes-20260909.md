# Calibration restoration reference notes — 2026-09-09

Primary sources consulted for the restoration validation protocol; these are
paraphrased retrieval notes, not full-text copies.

- Lucy L. Gao, Jacob Bien and Daniela Witten, *Selective Inference for
  Hierarchical Clustering*, JASA 119 (2024), 332–342; accepted manuscript
  https://arxiv.org/abs/2012.02936, accessed 2026-09-09. The abstract describes
  Type I error inflation when testing data-selected clusters and selective
  conditioning for hierarchical clustering. Its Gaussian cluster-mean result
  does not establish a law for this repository's adaptive Hamming diffusion,
  NNLS branch lengths, estimated projections or empirical calibration weights.
- Belinda Phipson and Gordon K. Smyth, *Permutation P-values Should Never Be
  Zero*, Statistical Applications in Genetics and Molecular Biology 9 (2010),
  article 39; author manuscript
  https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf, corrected 2011,
  accessed 2026-09-09. Randomly drawn permutation/Monte Carlo reference samples
  require valid discrete p-value calculations; replacing them with a raw
  exceedance fraction can understate p-values, especially under multiplicity.
  Exchangeability/null-model assumptions must hold before applying this result.

The project-specific protocol and local empirical evidence are in
`reports/calibration_restoration_20260909/validation_protocol.md`.

## 2026-09-10 — Population truth follow-up

Revisited Gao, Bien and Witten's accepted manuscript, equations (2)–(3), at
https://arxiv.org/html/2012.02936v3#S1. The selected-set population mean is an
average of underlying observation means, distinguished from its empirical mean.
The null compares these population means for the data-selected sets. This
informs the oracle truth definition in
`reports/calibration_restoration_20260909/false_split_definition.md`; it does not
transfer the paper's Gaussian inference guarantee to diffusion NNLS.
