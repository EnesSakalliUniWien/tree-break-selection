# Vandermonde applicability and application protocol

Frozen before execution, 2026-09-28. This companion study checks the covariance
matrices used in the known-branch experiment. It does not alter its saved data,
tree inference methods, covariance estimation, PCA/MP rules or production code.

1. Reconstruct all 10,000 original datasets from their declared random streams
   and verify their saved SHA256 hashes.
2. At each fixed root inspect sample pooled covariance (denominator N-1),
   known-covariance-whitened PCA covariance (denominator N), and sample-
   covariance-whitened PCA covariance (denominator N). Record positive
   definiteness and Toeplitz residual relative to the first-row/column Toeplitz
   matrix. Use structure tolerance 1e-10, with no Toeplitz fitting or averaging.
3. Apply Vandermonde factorization to both exact population covariances with
   anchors 0 and .137 cycles per feature index. Record their nonuniqueness.
   Also try the three sample matrices for the first global-null replicate in
   each profile, retaining unsupported matrices and reasons.
4. Include known rank-two complex frequencies {.12,.37}, weights {.7,1.3},
   and a real conjugate-pair example {.125,.875}, weights {.5,.5} as recovery
   controls. These are separate algebraic controls, not new calibration data.
5. Record frequencies, positive weights, reconstruction errors, atom condition,
   normalized overlaps, principal angles and projector differences between
   the atom span and positive covariance eigenspace. Orthonormalize before
   comparing subspaces. Full-rank span agreement is automatic and is not
   evidence of a preferred lower-dimensional projection.

Algorithm: guarded numerical implementation of the rank-deficient shift
relation, with one-atom subtraction for full rank. Rank tolerance is 1e-10
relative to the largest absolute eigenvalue; reconstruction tolerance is
1e-7 relative Frobenius error. Reject matrices failing Hermitian, Toeplitz,
PSD, positive-weight or reconstruction checks. No exact-rank certification or
near-colliding-frequency guarantee is claimed.

Store source hashes before and after execution. The theorem source and
interpretation are captured in `raw/inbox/vandermonde-reference-notes-20260928.md`.
