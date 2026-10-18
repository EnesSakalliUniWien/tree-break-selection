# Eigensystem API reference notes — 2026-09-28

These are short paraphrased source notes, not copies of the reference pages.
The current official SciPy pages displayed version 1.18.0 when read; the local
numerical reproductions use installed SciPy 1.16.3. Version-specific 1.16.3
page URLs were unavailable through the browsing tool. Numerical claims are
therefore established by the saved local probes, not inferred from a newer API.

- [SciPy eigh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html):
  solves real symmetric or complex Hermitian eigenproblems, returning normalized
  eigenvectors and eigenvalues repeated according to multiplicity. The function
  does not check the symmetry assumption itself.
- [SciPy eig](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html):
  distinguishes right eigenvectors satisfying `A r = lambda r` and left vectors
  satisfying `A.H l = conjugate(lambda) l`. Both families can be requested.
  Biorthogonal scaling is an additional mathematical condition, not assumed here.
- [SciPy subspace_angles](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.subspace_angles.html):
  computes principal angles between the full column spaces of two input matrices.
  The output is in radians, in descending order. The API supports different
  column counts, which is relevant to the unequal-rank diagnostic counterexample.

The Markov similarity formulas, projector identities and Gaussian quadratic-form
consequences in the review are direct algebraic arguments by the reviewing agents.
They are not presented as quotations or conclusions from these API pages.
