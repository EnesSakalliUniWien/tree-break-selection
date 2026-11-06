# Diffusion selected-root probability bounds and Wolfram check

Date: 2026-11-02. This is a diagnostic implementation record, not a production calibration claim.

## Model and exact path

For independent Gaussian rows with known identity covariance, fix the root contrast eta
specified by the observed complete merge history M. Under eta.T E[X]=0, put
v=X.T eta and H=eta.T eta. The conditional radial path is

    A = X - eta v.T/H
    B = eta v.T/sqrt(H * v.T v)
    X(r) = A + r B
    r_obs = sqrt(v.T v/H).

Original binary floating observations are treated as exact dyadic numbers. A is
computed using rational arithmetic; B and r_obs receive outward square-root bounds.
The floating replay's offset, slope and observed radius do not establish this exact path.

Conditional on the nuisance A and contrast direction, R has chi_p density, with p
equal to the full input feature count. Further conditioning on the complete merge
history restricts R to S_M. It does not change p to the retained diffusion rank.

    p_M = integral_{S_M intersect [r_obs,infinity)} f_chi_p(r) dr
          / integral_{S_M} f_chi_p(r) dr.

The event in the current certificate is explicitly the exact-real interpretation of
the adaptive diffusion formulas. It is not a bitwise NumPy/SciPy selection event.
The native positive-self-distance discrepancy documented separately is material.

## Unknown-region probability bounds

Let A be certified selected mass above the observed radius, B selected mass below,
U unresolved mass above, and V unresolved mass below. With positive certified
selected mass, monotonicity gives

    A/(A+B+V) <= p_M <= (A+U)/(A+U+B).

The implementation combines outward endpoints conservatively, includes every missing
finite region and the infinite tail, and uses [0,1] when no positive selected-mass
lower bound is available. An enclosing observed-radius interval is handled using
tail monotonicity: its upper endpoint supplies the lower p bound and its lower
endpoint supplies the upper p bound.

The root-only diagnostic rejects only when the upper p bound is at most alpha;
it reports nonrejection only when the lower bound exceeds alpha. All straddling,
unsupported or uncertified cases keep one cluster without passthrough. A rejection
would open two children; it would not authorize calibrated descendant traversal.

## Wolfram evaluator capture

Request through the connected Wolfram Language evaluator:

```wolfram
FullSimplify[
 {D[a/(a+b+v),a], D[a/(a+b+v),b],
  D[(a+u)/(a+u+b),u], D[(a+u)/(a+u+b),b]},
 Assumptions -> a>0 && b>0 && u>=0 && v>=0]
```

Returned:

```text
{(b + v)/(a + b + v)^2,
 -(a/(a + b + v)^2),
 b/(a + b + u)^2,
 -((a + u)/(a + b + u)^2)}
```

This checks the stated monotonic directions. It does not prove the selection
region, the numerical implementation, or native floating-point calibration.

## Arithmetic implementation

`diffusion_root_mass.py` uses directed Decimal arithmetic with 80 significant
digits. Correctly rounded exp and sqrt results are widened by adjacent Decimal
values. A positive incomplete-gamma series has an explicit geometric remainder;
even-dimensional survival probabilities use the finite exponential-polynomial
formula. Odd-dimensional far tails use alternating bounds derived from the erfc
integral and a half-integer recurrence. Pi is enclosed using the exact-rational
Machin identity. No quadrature error estimate is treated as a rigorous bound.

Tests include tiny lower-tail mass, far tails, disconnected selected unions,
missing/unsupported regions, observed-radius intervals, source-provenance checks,
and root rejection/nonrejection/unresolved decisions as mathematical controls.
The independently inspected primary spectral reference is Yu, Wang and Samworth,
Theorem 1 and its operator-norm remark:
https://www.statslab.cam.ac.uk/~rjs57/Biometrika-2015-Yu-315-23.pdf
Its application and exact residual derivation are in
`raw/inbox/diffusion-root-cells-certificate-20261102.md`.
