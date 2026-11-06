# Self-distance reconciliation and tighter diffusion selection bounds

Date: 2026-11-02. Subsequent user-confirmed implementation with three Sol agents.
This record supersedes the earlier unresolved self-distance implementation step;
it preserves the earlier raw diagnostic captures and their mathematical limits.

## Production correction

The repository-owned `adaptive_diffusion_geometry` now sets the identity
diagonal of `kernel_object.scaled_dists` to zero after fitting and before outer
epsilon selection and kernel computation. The change retains all off-diagonal
entries, including small positive distances and stored zero distances between
distinct observations. It does not call `eliminate_zeros`, which would remove
unit affinity for stored zero off-diagonal distances. Metadata identifies the rule.

The deterministic regression injects positive diagonal artifacts, a positive
off-diagonal distance of 1e-12 and a stored zero off-diagonal distance. It failed
before the correction and passes after it, including preservation of unit affinity.
No vendor, KDE algorithm, rank policy, inflation rule or statistical threshold changed.
Smooth native floating-point effects remain; exact-real interval certificates are
not bitwise floating-execution certificates.

## Genuine matched topology comparison

Before the production edit, the comparison agent captured original arrays,
source snapshot and native histories for 55 datasets: the previous 50 Gaussian
80-by-24 replicates, three small Gaussian seeds and two continuous signal cases.
All 50 Gaussian baseline histories match the earlier factorial hashes. The
fixed-source run used those captured arrays, with original and floating radial
reconstruction arms, giving 110 observations per version.

Final versioned evidence is under
`raw/assets/diffusion-self-distance-20261102/postfix_operator_v2/`.
The genuine baseline and its input archive are in the parent directory.

| Original-input comparison | 50 Gaussian replicates | Three small controls | Two signal cases |
| --- | ---: | ---: | ---: |
| Changed complete merge sequence | 8 | 0 | 0 |
| Changed rooted clade set | 0 | 0 | 0 |
| Changed root partition | 0 | 0 | 0 |

All 55 original/reconstructed pairs have matching complete histories within
each version. Known-identity root statistics are exactly unchanged. The
descriptive, unselected chi-square tail still rejects all 50 selected Gaussian
roots at .01. Thus numerical reconciliation does not repair selection bias.

The original Gaussian arm has 1,194 positive self entries before correction and
zero afterward. The maximum absolute difference between the actual native
eigensolver input and a floating zero-self reference operator is
5.551115123125783e-17 across the 110 corrected rows. This is floating consistency
evidence, not an interval proof. The baseline operator comparison was not
captured, and no legacy monkeypatch was used to invent it afterward.

One initial post-fix output was replaced while improving direct operator
observation. The versioned manifest discloses that replacement. The genuine
pre-fix baseline was preserved; the final CLI refuses raw-output overwrites.

## Kernel interval tightening

The exact Gaussian path gives pair squared distances

    q_ij(r) = c + b r + a r^2,
    c = ||A_i-A_j||^2,
    b = 2 sigma (eta_i-eta_j) (A_i-A_j).v / sqrt(H V),
    a = sigma^2 (eta_i-eta_j)^2/H,
    H = eta.T eta, V = v.T v, v = X.T eta.

Here a and c are exact rational values; b receives a directed enclosure. For
r>=0, q is monotone in b. The lower quadratic uses b_lower and its vertex when
inside the interval; the upper quadratic uses b_upper and its endpoint maximum.
Equal within-child translation cancels before interval evaluation, so within-child
distances retain the same bounds even across r in [0,10]. The old implementation
enclosed each coordinate first and lost this cancellation.

Neighbor certification now proves only the seven-neighbor pilot set and the
k-1 nonself kernel support set. Ordering within either selected set, and among
irrelevant distant neighbors, does not affect the exact-real sums. A boundary
ambiguity remains unresolved; it is not assigned a convenient neighbor set.

On the fixed scalar fixture at half-width 1e-6, maximum operator-entry interval
width decreases from about 5.99674656e-7 to 2.73876040e-7. Before/after evidence
and kernel hashes are `kernel_bounds_before.json` and `kernel_bounds_after.json`
under `raw/assets/diffusion-self-distance-20261102/`.

Wolfram checked the quadratic vertex substitution:

```wolfram
FullSimplify[
 c + (2 s w d/Sqrt[h v2]) r + (s^2 w^2/h) r^2 /.
 r -> -(s w d/Sqrt[h v2])/(s^2 w^2/h),
 Assumptions -> h>0 && v2>0 && s>0 && w!=0]
```

Returned `c - d^2/v2`. This checks the algebraic minimum formula, not the full
selection event or numerical implementation.

## Spectral tightening and final integrated bounds

The new factor bound uses H=S^t P, so distances are ||H(e_i-e_j)||. It projects
the matrix-power difference onto the retained spectra before bounding it, rather
than including the removed leading eigenvalue in every power norm. The full
derivation is in `diffusion-factor-bound-tightening-20261102.md`.

On identical upstream operator boxes, the scalar fixture's largest tested
certified half-widths are 1e-12 for the old Gram bound, 1e-7 for the global
factor bound and 1e-5 for the restricted factor bound. These are tested widths,
not estimated selection boundaries. Interval expansion accepts a larger window
only after a whole-window proof; finite subdivision and the infinite tail retain
all unresolved mass.

Final integrated reports are under
`raw/assets/diffusion-region-tightening-20261102/final-small-controls/` and
`final-attribution-first-root/`. Small controls use half-width 1e-5, four expansion
attempts and finite subdivision depth 4. The 80-row case uses half-width 1e-8,
two expansion attempts and subdivision depth 1. All use full features, known
identity covariance, zero branch lengths, inflation factor 1 and no
contrast-energy rank selection.

For seed 20261105, certified selected probability has lower bound approximately
1.058576795e-6; p lies approximately in [7.16085e-7, 0.999972]. The first
80-row observation now belongs to a certified selected cell, with selected-mass
lower bound approximately 4.884791985e-14. All four gates remain unresolved.
The third small control still lacks certified observed-event membership.
Neither local certification nor one-cluster unresolved outcomes establishes
selected-null calibration or descendant traversal validity.

## Verification

The complete ordered suite passes 977 tests. Repository Ruff lint, dependency
and dead-code audits, generated-file checks and 228-page wiki lint pass.
Formatting checks pass for all 16 checked Python files, including the new
reproduction runner. The repository-wide format check still flags the unchanged earlier
`raw/assets/diffusion-selection-event-20261102/reentry_probe.py`; that raw
capture was preserved. Final panel source hashes match current files. Of 207
previously captured production/shared source files, only the intended
`tree_break_selection/space_separation/diffusion.py` changed. All unrelated dirty
working changes were preserved.
No dependency installation, commit, push, release or persistent service occurred.
