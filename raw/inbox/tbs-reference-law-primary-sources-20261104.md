# Primary references for the TBS reference-law audit

Accessed 2026-11-04. These references concern restricted statistical laws.
The current implementation and the controls establish their scope in TBS.

- Fithian, Sun and Taylor, *Optimal Inference After Model Selection*,
  version 4 (2017). The primary paper defines selective Type I error
  conditional on selection and treats data splitting and conditioning.
  A TBS edge split applied after building the tree does not provide
  the independence assumed by inference on data untouched by selection.
  [Author preprint](https://arxiv.org/abs/1410.2597)
- Hotelling, *The Generalization of Student's Ratio* (1931), Annals of
  Mathematical Statistics 2, 360–378. Classical Gaussian sample-mean and
  inverse-covariance framework for the finite-sample F control. The original
  journal URL did not open through the browsing tool. The result used here
  is also derived directly by the Wishart Schur complement in the audit.
  [Original publication record](https://cir.nii.ac.jp/crid/1360011144703680000)
- Muirhead and Chikuse, *Approximations for the Distributions of the Extreme
  Latent Roots of Three Matrices* (1975). Their central Wishart ratio setup
  explicitly uses independent matrices with the same covariance. The audit
  derives the parent-total scaled-Beta and adaptive-energy Gaussian controls
  separately from row rotational invariance and the gamma/Dirichlet identity.
  [Journal-hosted primary article](https://www.ism.ac.jp/editsec/aism/pdf/027_3_0473.pdf)
