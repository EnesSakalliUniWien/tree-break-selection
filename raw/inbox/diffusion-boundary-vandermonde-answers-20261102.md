# Diffusion boundary recovery and Vandermonde applicability — 2026-11-02

Read-only production research. No production changes, dependencies, services,
commits, or repaired inference claim. Q3 and Q4 are answered below; parent agent
owns wiki synthesis. This note distinguishes a procedure that would certify an
answer from the currently available floating-point diagnostic implementation.

## Q3: Can Newton–Raphson recover all selection boundaries?

**No. Newton–Raphson can refine a locally isolated smooth boundary; it does not
supply complete boundary recovery.** Even a converged scalar solve does not
prove that no other boundaries exist. SciPy's documentation explicitly separates
Newton's local root search and its step-size stopping rule from bracketed
continuous sign-changing solves:
[Newton](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html),
[Brent](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html).
Brent can refine a valid sign-changing bracket, but a grid of brackets misses
even-multiplicity roots, close pairs of roots between grid points, and jumps.

For the known-sigma Gaussian radial path X(r)=A+rB, selection means replaying the
specified event at every r>=0. For the full observed average-linkage merge
sequence on diffusion Euclidean distances, at every observed merge t and every
other currently active pair (C,D), require

    g[t,C,D](r) = mean_{i in A_t,j in B_t} d_diff(i,j;r)
                 - mean_{i in C,j in D} d_diff(i,j;r) <= 0,

with the actual tie-breaking policy handled at equality. This describes a
sequence event; a root-partition event can be a union over sequences and must
not silently be replaced by one sequence. If the event also includes edge
openings or traversal outcomes, include their test/rank/BH/stopping inequalities;
linkage inequalities alone are not that larger event.

Unlike the already implemented direct squared-Euclidean Gaussian diagnostic,
g is not generally a quadratic polynomial. The actual continuous path is:

1. Euclidean nearest-neighbor search, requested k=15.
2. pydiffmap variable-bandwidth KDE (`-1/(d+2)`), including an internal BGH
   epsilon-grid maximizer and rounded dimension estimate.
3. Repository outer epsilon = median positive scaled squared distances.
4. Directed sparse Gaussian weights; remove diagonal and symmetrize by maximum.
5. S(r)=D(r)^(-1/2) W(r) D(r)^(-1/2).
6. Retain target min(50,n-2,p) eigenvectors, remove leading mode only when
   isolated, extend through numerical ties, clip negative eigenvalues, time=3.
7. Euclidean diffusion-coordinate distances and average linkage.

Source: `tree_break_selection/space_separation/diffusion.py`, functions
`feature_family_diffusion_geometry`, `adaptive_diffusion_geometry`,
`compute_diffusion_coordinates`, `_leading_diffusion_spectrum`; installed
`.venv/lib/python3.12/site-packages/pydiffmap/kernel.py`, `Kernel`, `NNKDE`,
`choose_optimal_epsilon_BGH`. This is the current continuous preset; discrete
Hamming paths have a different support/domain problem.

### Where differentiability helps, and where it stops

On a branch where neighbor identities, KDE discrete choices, median order,
positive degrees, and retained spectral group are fixed and smooth, a simple
normalized eigenpair has

    lambda_i' = u_i.T S' u_i,
    u_i' = sum_{j != i} u_j (u_j.T S' u_i)/(lambda_i-lambda_j),

using a locally normalized eigenvector sign convention. These are local
simple-eigenvalue formulas, not formulas across a zero spectral gap. The primary
perturbation reference is Greenbaum, Li and Overton,
[First-order Perturbation Theory for Eigenvalues and Eigenvectors](https://arxiv.org/abs/1903.00785).

At ties, individual eigenvectors are not identifiable; use the projector of an
entire separated eigenspace group when the downstream quantity is invariant to
within-group rotations. Distances derived from all equally weighted modes in an
exact repeated group are invariant; differentiating a selected arbitrary vector
is not. Numerical grouping still has its own switching boundaries: the live
`eigenvalue_groups.py` uses relative tolerance 1e-10 and absolute tolerance
1e-12 times spectral scale, with transitive adjacent ties. Entering/leaving that
numerical group can change retained dimension even without an exact crossing.

Other required branch boundaries are:

- Neighbor-distance equalities, including support changes. Along A+rB, raw
  squared-distance differences are quadratic, but this does not make later
  scaled-distance or linkage comparisons quadratic.
- BGH maximizer changes, rounded KDE dimension changes, local-bandwidth
  singularities/duplicate rows, outer median order changes.
- Maximum-symmetrization branch changes and degree-zero safeguards.
- Eigenvalue crossings, clipping at zero, numerical tie-group changes, leading
  mode policy changes, and any downstream MP or contrast-energy rank changes
  if that rank participates in the requested event.
- Final linkage ties, including tangencies and an identically zero comparison.

Some are continuous but nonsmooth; others can jump. A derivative through a
single execution trace does not describe the other side of these branches.
Unstable eigenvector signs are a coordinate issue, not independently a distance
boundary when the represented spectral subspace and weights are unchanged.

### What a complete/certified procedure would require

The following is a requirements list, **not an implemented complete algorithm**:

1. Declare the exact selection event, tie policy, radial null law and domain.
   Make every replay branch/unsupported state observable. Freezing a graph or
   KDE state can simplify derivatives but defines a stronger conditioning event
   if it is imposed on the statistical law; it is not automatically equivalent
   to conditioning only on the resulting tree.
2. Cover a bounded interval by cells with validated enclosures for every
   relevant branch function. Isolate algebraic raw-neighbor transitions using
   polynomial root isolation where justified. For later functions, interval or
   ball enclosures must include KDE exponentials, bandwidths, normalization,
   spectral gaps/projectors, coordinate distances and linkage comparisons.
3. Classify a cell only when enclosure signs and branch membership are certified.
   Otherwise subdivide or retain it as unresolved. Use validated interval Newton
   to certify a unique simple zero when its hypotheses hold. Ordinary Newton
   may propose candidates; it must not exclude cells. Tangencies need explicit
   multiple-root/stationary-point treatment or unresolved enclosures; an
   identically zero comparison needs the actual tie policy. One-sided replay
   and enclosure checks are needed at nonsmooth boundaries or jumps.
4. Cover the whole domain, including r=0, every unresolved transition, and
   r=infinity. A proof that choices/signs stabilize beyond R would close the
   tail. Mapping r=tan(theta), theta in [0,pi/2), merely compactifies the domain;
   it does not certify the endpoint or polynomialize the diffusion pipeline.
   Without a tail certificate, all roots have not been recovered.
5. Integrate only certified selection components, retaining upper bounds for
   unresolved interval/tail null mass. A chi_d survival probability beyond R
   bounds unconditional omitted mass, not automatically conditional error:
   it must be compared with a positive lower bound on selection probability.
   For computed numerator N and denominator D>0 and unresolved null mass delta,
   safe coarse bounds are N/(D+delta) <= p <= min(1,(N+delta)/D).

Verified all-root box methods are an established tool class, but this repository
has not implemented validated enclosures for this whole pipeline. See the
software author's [INTLAB global roots demonstration](https://www.tuhh.de/ti3/rump/intlab/demos/html/dglobal.html).
No termination/completeness guarantee has been established here for arbitrary
inputs, degenerate roots or an unbounded diffusion path. A finite successful
grid is evidence of sampled replay only.

## Q4: Does trigonometric/Vandermonde structure solve this problem?

**Not for the current arbitrary-data diffusion operator.** The relevant theorem
requires a Hermitian positive-semidefinite Toeplitz matrix T. It represents
T as A diag(w) A*, where A[j,k]=exp(2 pi i j f_k), positive weights w and distinct
frequencies f_k. Rank below matrix dimension gives uniqueness; full rank does
not. These frequencies and atoms are not generally matrix eigenvalues and
orthogonal eigenvectors. Source: Yang and Xie,
[Frequency-Selective Vandermonde Decomposition of Toeplitz Matrices with Applications](https://arxiv.org/html/1605.02431),
Theorem 1. A circular shift-invariant/circulant operator is a useful separate
special case where Fourier modes diagonalize it; no such shift structure is
established for arbitrarily indexed observed samples here.

The runtime probe below checks the actual continuous n80,p24,seed2026110203
geometry and reproduces its coordinate Gram matrix through the existing public
kernel/operator APIs. It does not monkeypatch production. Toeplitz residual is
the existing helper's Frobenius residual against the matrix fixed by first row
and first column, not a best Toeplitz fit.

| Matrix | Symmetry residual | Toeplitz residual | Minimum eigenvalue of symmetric part |
| --- | ---: | ---: | ---: |
| Directed sparse raw kernel | 0.8689988865 | 1.095146985 | -2.589272804 |
| Symmetric zero-diagonal W | 0 | 1.165773822 | -4.843196459 |
| Normalized S | 0 | 1.262411810 | -0.3174827773 |
| Diffusion-coordinate Gram | 0 | 1.397353908 | -1.6694e-18 |
| Sample feature covariance | 0 | 0.6594936689 | 0.2256126790 |

The directed kernel fails Hermitian structure. W and S fail Toeplitz and PSD;
the coordinate Gram and feature covariance are PSD up to numerical roundoff but
fail Toeplitz. There is also a general elementary obstruction for this exact
operator construction: any nonzero symmetric zero-diagonal S has trace zero,
so it cannot be PSD, since otherwise all eigenvalues would be nonnegative and
their zero sum would force S=0. This argument does not depend on the one seed.

Retained evidence agrees: in
`reports/independent_projection_vandermonde_20260928/vandermonde/structure_summary.csv`,
all 10,000 sample-pooled and all 10,000 oracle-whitened PCA covariances failed
Toeplitz, while sample-whitened PCA covariances passed because they were scalar
identities. Structured population identity/AR and rank-two controls were
separate successful applications. That report concerns covariance controls,
not current diffusion operators, so it does not prove their structure.

The existing guarded `benchmarks/diagnostics/spectral/vandermonde.py` appropriately
rejects non-Hermitian/non-Toeplitz/non-PSD inputs and never fits a Toeplitz matrix.
Forcing a Toeplitz/circulant approximation, reordering samples to seek stationarity,
adding a PSD shift, or retaining a Fourier truncation changes the operator or
geometry. Those would need a declared model, approximation error and renewed
selection-law validation; they are not consequences of the theorem. Expressing
an arbitrary matrix in a Fourier basis is always possible, but its transformed
matrix is generally dense, so the selected-boundary problem remains.

## Reproducible structural probe

Environment: Darwin/zsh, Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, existing
pydiffmap. Run the following with `uv run --no-sync python` from repository root.

```python
import json
import numpy as np
from pydiffmap import kernel
from benchmarks.diagnostics.spectral.vandermonde import (
    toeplitz_residual, decompose_psd_toeplitz,
)
from tree_break_selection.space_separation.diffusion import (
    feature_family_diffusion_geometry, compute_diffusion_coordinates,
    resolve_adaptive_epsilon,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

x = np.random.default_rng(2026110203).normal(size=(80, 24))
space = continuous_feature_space_from_columns([f'F{i}' for i in range(24)])
k = kernel.Kernel(epsilon=1., k=15, metric='euclidean',
                  neighbor_params={'n_jobs': -1}, bandwidth_type='-1/(d+2)')
k.fit(x)
epsilon, _ = resolve_adaptive_epsilon(k, 'median', metric='euclidean')
k.epsilon_fitted = epsilon
K = k.compute().toarray()
W = K.copy()
np.fill_diagonal(W, 0)
W = np.maximum(W, W.T)
d = W.sum(1)
d[d == 0] = 1e-10
S = W / np.sqrt(d[:, None] * d[None, :])
y = compute_diffusion_coordinates(K, diffusion_time=3, n_components=24)
native = feature_family_diffusion_geometry(x, space)
assert np.allclose(y @ y.T, native.coordinates @ native.coordinates.T,
                   rtol=1e-10, atol=1e-12)
print('seed=2026110203 shape80x24 epsilon=', epsilon, 'components=', y.shape[1])
for name, A in [('raw_kernel', K), ('symmetrized_zero_diagonal_W', W),
                ('symmetric_normalized_S', S), ('diffusion_coordinate_Gram', y@y.T),
                ('sample_feature_covariance', np.cov(x, rowvar=False))]:
    sym = np.linalg.norm(A-A.T) / np.linalg.norm(A)
    eig = np.linalg.eigvalsh((A+A.T)/2)
    try:
        decompose_psd_toeplitz(A)
        result = 'applicable'
    except ValueError as e:
        result = str(e)
    print(json.dumps(dict(matrix=name, dimension=len(A),
          symmetry_residual=float(sym), toeplitz_residual=toeplitz_residual(A),
          min_symmetric_part_eigenvalue=float(eig[0]),
          max_symmetric_part_eigenvalue=float(eig[-1]), guard=result)))
assert np.trace(S) == 0 and np.linalg.eigvalsh(S)[0] < 0
```

Captured output:

```text
seed=2026110203 shape80x24 epsilon= 0.8588652016711084 components= 24
{"matrix": "raw_kernel", "dimension": 80, "symmetry_residual": 0.8689988864816581, "toeplitz_residual": 1.0951469846395296, "min_symmetric_part_eigenvalue": -2.5892728036050308, "max_symmetric_part_eigenvalue": 13.853880913699847, "guard": "not_hermitian"}
{"matrix": "symmetrized_zero_diagonal_W", "dimension": 80, "symmetry_residual": 0.0, "toeplitz_residual": 1.1657738216939653, "min_symmetric_part_eigenvalue": -4.843196458542079, "max_symmetric_part_eigenvalue": 18.12582949777284, "guard": "not_toeplitz"}
{"matrix": "symmetric_normalized_S", "dimension": 80, "symmetry_residual": 0.0, "toeplitz_residual": 1.262411809568348, "min_symmetric_part_eigenvalue": -0.3174827772853742, "max_symmetric_part_eigenvalue": 0.9999999999999998, "guard": "not_toeplitz"}
{"matrix": "diffusion_coordinate_Gram", "dimension": 80, "symmetry_residual": 0.0, "toeplitz_residual": 1.3973539076262556, "min_symmetric_part_eigenvalue": -1.6694281186869647e-18, "max_symmetric_part_eigenvalue": 0.0073253418553650505, "guard": "not_toeplitz"}
{"matrix": "sample_feature_covariance", "dimension": 24, "symmetry_residual": 0.0, "toeplitz_residual": 0.6594936688973454, "min_symmetric_part_eigenvalue": 0.22561267899048568, "max_symmetric_part_eigenvalue": 2.120358487496096, "guard": "not_toeplitz"}
```

Checks: executed probe assertions passed; no production files edited. No
certified interval implementation, complete boundary recovery, or automatic
Vandermonde applicability beyond the tested/structurally excluded matrices was
claimed.
