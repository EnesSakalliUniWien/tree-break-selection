# Diffusion factor-bound tightening — 2026-11-02

This is a subsequent implementation record. It preserves the earlier
`diffusion-root-cells-certificate-20261102.md` proof history. Production changes
and native/exact-real reconciliation belong to the other owning records; this
note concerns the interval classifier in
`benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_cells.py`
and `tests/statistics/62_test_diffusion_root_cells.py`.

## Distance factor rather than squared-distance perturbation

Let S be a real symmetric normalized diffusion operator, P its certified retained
spectral projector, and t=3. The supported certificate establishes positive
retained eigenvalues, an isolated leading mode and an isolated retention cutoff,
including the exact dyadic native numerical tie thresholds. Diffusion-coordinate
Gram G=S^(2t)P can be factored as H H.T, where H=S^t P. Consequently

    d(i,j) = ||H(e_i-e_j)||.

The candidate eigenpairs Q,lambda are floating-point proposals. Exact dyadic
orthogonality/residual calculations establish an exact orthogonal polar factor
U, exact symmetric reference A=U diag(lambda) U.T and uniform bound
`||S-A||<=delta`. Let P* denote the matching reference projector, and
`H*=A^t P*`. Native adaptivity is enclosed upstream for the entire cell; the
candidate eigendecomposition does not certify any branch by itself.

The exact rational reference factor is `Hhat=Q diag(selected lambda^t) Q.T`.
If `epsilon>=||Q.T Q-I||F<0.1`, its polar correction is at most

    ||H*-Hhat|| <= mA^t epsilon(2+epsilon),

where `mA=max(abs(lambda[selected]))`. For exact reference factor difference
vectors, integer square-root bounds enclose `||Hhat(e_i-e_j)||`. The triangle and
reverse triangle inequalities then enclose actual pair distances with additive
error `sqrt(2)||H-Hhat||`. This error is linear in the operator perturbation even
at near-zero reference distances; the former absolute Gram perturbation followed
by square roots lost that property.

## Restrict the power bounds to the retained spectrum

A full-spectrum power norm unnecessarily includes the dropped leading mode near
one. The following exact identity avoids it:

    H-H* = H(I-P*) - (I-P)H* + P(S^t-A^t)P*.

It follows by expansion, using PS=SP and P*A=AP*. With
`E_P>=||P-P*||`, `mS=mA+delta`, Weyl's bound gives

    ||H(I-P*)|| <= mS^t E_P,
    ||(I-P)H*|| <= mA^t E_P.

The projected power telescope is

    P(S^t-A^t)P*
      = sum_{k=0}^{t-1} P S^(t-1-k) (S-A) A^k P*.

Each factor adjacent to a projector is restricted to the retained spectrum, so

    ||H-Hhat|| <= (mS^t+mA^t) E_P
                 + delta sum_{k=0}^{t-1} mS^(t-1-k) mA^k
                 + mA^t epsilon(2+epsilon).

The implementation takes the smaller of this valid restricted bound and the
previous full-spectrum factor bound. The projector error continues to use the
conservative Davis–Kahan bounds at both prefix cuts. The applicable
operator-norm theorem is Theorem 1 in Yu, Wang and Samworth,
[A useful variant of the Davis–Kahan theorem for statisticians](https://www.statslab.cam.ac.uk/~rjs57/Biometrika-2015-Yu-315-23.pdf).
No statistical alpha, rank policy, or native tie tolerance was tuned.

## Bounded expansion and unresolved-region coverage

`expansion_steps` controls independent doubling of the interval around each
requested anchor. An expanded interval is accepted only when the whole interval
proves the SAME selected/nonselected label. Failed proof stops that expansion.
Adjacent-anchor midpoints prevent overlap. Cached accepted certificates are
reused when writing the final partition; every gap and `[R,infinity]` is retained.
`subdivision_depth` separately permits bounded subdivision of unresolved finite
gaps. It never converts samples to labels.

Each anchor's validated reference eigenpairs may be reused for its larger trial
intervals. This saves candidate/residual setup work. Every larger kernel box is
still recomputed and the exact residual, orthogonality, spectral gap and
linkage inequalities are revalidated against that box. Reusing a reference
matrix does not freeze S, neighbors, KDE, epsilon, eigenvalues or event history.

The observed-event validity flag still requires a certified selected interval
containing the entire outward enclosure of the exact Gaussian observed radius.
Unknown regions, arithmetic failures, ambiguous branches, unresolved ties and
the unbounded tail retain probability mass. The gate remains unresolved unless
its p bounds determine a decision.

## Matched evidence and limits

`raw/assets/diffusion-region-tightening-20261102/compare_spectral_bounds.py` is the
captured reproducible probe. It supplies the SAME native kernel/operator box to
three methods: archived prior Gram engine, archived full-spectrum factor engine,
and current restricted-spectrum factor engine. The archived source snapshots
are `gram_engine_before.py.gz` and `factor_global_before.py.gz`; current owning
source files and hashes are recorded in `spectral_postkernel_comparison.json`.

For seed20261105, n9,p1, the observed selected window is certified at tested
half-widths:

- Gram: through1e-12; unresolved at1e-10 and above.
- Full-spectrum factor: through1e-7; unresolved at1e-6 and above.
- Restricted-spectrum factor: through1e-5; unresolved at1e-4 and above.

These tested limits are not claimed to be optimal boundaries. The comparisons
are interval proofs using one matched kernel box, not Monte Carlo estimates.
Focused regressions failed against the earlier methods at1e-8 and1e-6 and pass
with the corresponding tighter methods.

A fresh post-reconciliation native n80,p24,seed2026110203 path additionally
certifies the observed interval at half-width1e-8 using the restricted factor;
the previous full-spectrum factor left that same-width interval unresolved.
The certified selected chi mass lower bound is about2.442396209e-14. The
complete-domain p bounds remain unresolved and return one cluster. See
`actual80_restricted_factor.json`. This is an actual positive selection-cell
certificate, not complete selected-region coverage or calibrated production
traversal.

All claims concern the declared exact-real constructor with outward enclosures,
not proof of identical floating-point execution at every radius. Broad unknown
regions and the tail remain the material limitation even after substantial
local improvements.

Final three-anchor post-kernel growth capture (`window_growth_postkernel.json`)
uses the same seed and ratios1,1.1,1.2. Its total certified selected chi mass lower
bound grows from about8.807706536e-14 at initial half-width1e-12, to
8.806923596e-7 at initial half-width1e-5, and to1.058576795e-6 with six permitted
expansion attempts. The accepted observed selected window has width2e-5;
the reentry selected and intervening nonselected windows each reach width4e-5.
Further expansion fails a proof and stops. All three gates remain unresolved.
Source snapshots are unchanged during the final matched spectral comparison,
all three final growth runs, and the actual80 restricted-factor run.
