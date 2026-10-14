# Vandermonde decomposition source notes

Accessed 2026-09-28. Primary source: Zai Yang and Lihua Xie,
*Frequency-Selective Vandermonde Decomposition of Toeplitz Matrices with
Applications*, https://arxiv.org/abs/1605.02431,
https://arxiv.org/pdf/1605.02431, Theorem 1 and its construction, pp. 4–6.

Paraphrased mathematical contract: a Hermitian positive semidefinite Toeplitz
matrix of rank r admits T = A diag(p) A^H, with positive weights and
sinusoidal columns a(f)[j] = exp(2 pi i j f). The decomposition is unique
when r is smaller than the matrix dimension. In the full-rank case an
arbitrary frequency can be included with weight 1/(a^H T^{-1}a); subtracting
that atom produces a rank-deficient Toeplitz matrix. Frequencies can then be
obtained from the shift relation between adjacent rows of a covariance
factor. Frequencies use the circle convention f modulo 1.

Implementation interpretation: numerical rank, structural residuals,
positive weights and reconstruction must be checked explicitly. Applying
this theorem to an unrestricted covariance matrix would require additional
structure or a separately justified fitting step. No such fitting step is
authorized by the theorem itself.
