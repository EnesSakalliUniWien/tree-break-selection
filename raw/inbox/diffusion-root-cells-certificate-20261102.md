# Restricted diffusion selection-cell certificates — 2026-11-02

This records the subsequent diagnostic implementation, extending the earlier
read-only `diffusion-boundary-vandermonde-answers-20261102.md` investigation.
Production code is unchanged. This is partial validated coverage of an
exact-real mathematical constructor, not bitwise certification of arbitrary
NumPy/SciPy executions or complete boundary recovery.

## Owning implementation and contract

- `benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_cells.py`
- `tests/statistics/62_test_diffusion_root_cells.py`
- Upstream `diffusion_root_kernel.py`: outward kernel/operator enclosures.
- Downstream `diffusion_root_mass.py`: outward chi mass and conditional p bounds.
- `diffusion_root_replay.py`: actual floating-point point replay, never a cell
  certificate by itself.

The scanner `analyze_diffusion_root_cells` accepts only the path plus partition
controls. It calls the owning kernel encloser itself; it does not accept sampled
labels, arbitrary caller certificates, estimated Lipschitz constants, or a
caller-supplied `validated=True` flag.

Original observations are treated as exact dyadic numbers. The Gaussian path is
re-derived from those observations and rational root contrast eta. Its exact
nuisance offset is `X-eta(eta.T X)/||eta||²`; its slope contains the exact positive
square root of `||eta||² ||eta.T X||²`. The kernel helper encloses this path.
The scanner independently bounds the exact observed chi radius using integer
square-root arithmetic. Approximate replay offset/slope are not proof inputs.

## Shared spectral residual proof

For each operator box `[S_lower,S_upper]`, a floating-point eigensolver proposes
ordered eigenvalues lambda and vectors Q at its midpoint. Convert these floats
to exact dyadic rationals. Compute

    epsilon >= ||Q.T Q-I||_F,
    R(r) = S(r)-Q diag(lambda) Q.T.

If epsilon<0.1, the orthogonal polar factor U of Q exists and satisfies
`||Q-U||_op <= epsilon`. Therefore the exact symmetric reference
`A*=U diag(lambda) U.T` has the proposed eigenvalues, and

    delta = uniform_Frobenius_bound(R)
            + max|lambda| epsilon(2+epsilon)

bounds `||S(r)-A*||_op` on the whole cell. All scalar algebra here is exact
Fraction arithmetic; square-root upper bounds use integer `isqrt` with an
explicit positive rounding increment.

The upstream enclosure follows the native dense route through n=2000; this is an
applicability limit, not a scaling/performance guarantee. The currently supported
spectral route requires an isolated leading eigenvalue,
an isolated retained cutoff, and strictly positive retained eigenvalues after
allowing delta. It also proves that native numerical eigenvalue-group tolerances
cannot merge either cutoff, using the exact dyadic values of the native Python
`1e-10` and `1e-12` literals. Arithmetic failures leave a cell unresolved. Ties within the retained span are harmless to the
spectral projector, but a tie at a cutoff remains unresolved.

For each prefix cut with reference gap g, use the conservative projector bound
`2 delta/(g-2 delta)` with g>2 delta. The selected projector is the difference of
two prefix projectors (remove the isolated leading mode); add their errors and
cap at one, since differences of orthogonal projectors have operator norm at
most one. The applicable primary reference is the operator-norm form of
Theorem 1 in Yu, Wang and Samworth,
[A useful variant of the Davis–Kahan theorem for statisticians](https://www.statslab.cam.ac.uk/~rjs57/Biometrika-2015-Yu-315-23.pdf).
Weyl's eigenvalue bound supplies the required spectral separation.

The initial implementation used the following Gram bound. For native time t=3
and positive selected eigenvalues, the coordinate Gram matrix
is `G=S^(2t)P`. With L=max|lambda| and M=L+delta, telescoping matrix powers gives

    ||G-Q diag(selected lambda^(2t)) Q.T||_op
      <= 2t M^(2t-1) delta
         + L^(2t) projector_error
         + L^(2t) epsilon(2+epsilon).

A pair squared distance is `(e_i-e_j).T G (e_i-e_j)`. Its error is at most twice
the Gram operator error. Integer-root outward bounds then enclose every
Euclidean pair distance, including distances whose lower bound reaches zero.
Averages over descendant pairs enclose each average-linkage comparison.

A cell is `selected` only if every observed chosen merge uniformly beats every
competitor beyond the specified tie-tolerance bound. It is `nonselected` if a
competitor is uniformly strictly closer at any step of that observed sequence.
Any remaining ambiguity is `unresolved`. A failed neighbor/KDE/spectral proof
is also unresolved, never nonselection.

## Exhaustive bookkeeping and current positive result

Default partitioning includes a small window around the observed radius, all
intervening gaps, a finite upper bound R, and `[R,infinity]`. Optional
`anchor_ratios` add windows; they influence partition locations only, never
labels. Subdivision similarly requires a new whole-cell proof for every label.
The unbounded tail is currently unresolved and retains its probability mass.

`selection_certified=True` requires that a proved selected cell contain the
entire outward interval for the exact observed radius. Otherwise the root gate
must remain unresolved. All source hashes used by replay, kernel, cells and
mass helpers are recorded for stale-certificate checks.

For original observations `default_rng(20261105).normal(size=(9,1))`, anchors
`(1,1.1,1.2)` and half-width `1e-12`, the scanner proves positive-width windows:

| Radius ratio | Uniform certificate |
| --- | --- |
| 1 | selected |
| 1.1 | nonselected |
| 1.2 | selected |

The first selected window has a lower linkage margin around `4.78028985e-6`
and an upper Gram error around `3.19082801e-11`. The second selected window's
lower margin is around `1.69726849e-5`. These are interval certificates of
full-history reentry, not point or grid labels. Each window has actual float
endpoint width about `2.0001778e-12`; detailed proof bounds are reported as
outward Decimal strings, while internal algebra remains rational.

The remaining radial regions and tail are unresolved. Thus selected probability
has a positive lower bound, but conditional p bounds remain too wide to decide
the root at alpha 0.01. The diagnostic gate returns unresolved and one cluster.
This establishes useful local certification and bookkeeping, not all roots or
an informative final selected p-value.

## Verification

The focused cell tests cover exact integer square-root bounds, a nonzero
symmetric-operator box, competing full histories, spectral-gap ambiguity,
linkage ties, actual adaptive-kernel reentry windows, entire-domain coverage,
observed exact-radius membership, positive selected probability lower bound,
and unresolved root behavior. Tests comparing exact-real to floating-point
values are consistency checks; certificates are established by the enclosure
calculations, not those comparisons.


## Tightened factor certificate and bounded growth

The current implementation replaces the initial squared-distance Gram route with
`H=S^t P`, a symmetric diffusion factor satisfying `G=H H.T`. Its exact rational
reference is `Hhat=Q diag(selected lambda^t) Q.T`. With the same delta, L, M and
projector bounds as above,

    ||H-Hhat|| <= t M^(t-1) delta
                 + L^t projector_error
                 + L^t epsilon(2+epsilon).

The true pair distance is `||H(e_i-e_j)||`, and the reference norm is computed
from all entries of `Hhat[:,i]-Hhat[:,j]` with exact rational arithmetic and
outward integer square roots. The reverse triangle inequality gives distance
error at most `sqrt(2)||H-Hhat||`. This bound is linear in operator error even
when the reference pair distance is near zero; it avoids the square-root loss
from perturbing a squared distance by an absolute Gram error.

`expansion_steps` independently bounds window growth at each requested radial
anchor. Each doubled window must prove the same selected/nonselected label on
its entire enlarged domain. A failed proof stops that expansion. Windows are
clipped at adjacent-anchor midpoints, so the final partition has no overlap.
All gaps and the unbounded tail remain represented. Existing `subdivision_depth`
can subdivide unresolved finite cells; a child label still requires its own
complete proof. Candidate locations and expansion are numerical work controls,
not statistical thresholds or sampled classifications.

In `raw/assets/diffusion-region-tightening-20261102/spectral_same_kernel_comparison.json`,
the old and new spectral methods receive the SAME outward operator box at each
width, isolating this spectral change from upstream kernel improvements. For
seed20261105, n9,p1, the old Gram method certifies half-width1e-12 with selected
chi mass lower bound about7.02870404e-14; it cannot certify half-width1e-10 or1e-8.
The factor method certifies half-width1e-8 with selected mass lower bound about
7.02807920e-10. Both methods remain unresolved at half-width1e-7 in that recorded
comparison. These are genuine positive-width interval proofs; they still cover
only a small part of the selected region, and do not make conditional p bounds
informative by themselves. Source hashes were unchanged during the comparison.

Additional focused regressions require the actual1e-8 window to become certified
(the test failed against the old engine), prove growth increases selected mass,
and verify that expansion preserves complete domain and unresolved-tail
bookkeeping. No production threshold was changed.
