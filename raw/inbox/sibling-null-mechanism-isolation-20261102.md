# Sibling null mechanism isolation — 2026-11-02

Diagnostic counterexamples, not production calibration or a proposed repair.
Executed from repository root with `uv run --no-sync python`; Darwin, zsh,
Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3. No dependencies installed.

The projection experiment calls the production shared projected-Wald kernel.
Each replicate is an independent 40-by-8 iid standard Gaussian matrix with a
fixed first-20/last-20 split. Covariance is known, so the mean contrast has
component variance 0.1. Rank two is fixed except in the explicit contrast-energy
control. Same-data PCA uses all 40 centered leaf rows, without whitening by its
eigenvalues. The branch control divides z by sqrt(2), holding that PCA fixed.
These controls omit topology selection, edge selection, MP rank, estimated
covariance, empirical calibration and traversal. They do not replicate a full
production sibling record. All decisions use p < 0.05.

The fitted-scale experiment uses independent focal and calibration chi-square
draws and true scale three. D is total calibration degrees of freedom. The
plug-in control includes the production no-deflation floor one. The exact F
control uses a **separate independent calibration draw**, not the plug-in draw;
it is a marginal validity control, not a paired comparison. For each D the RNG
consumes focal, plug-in calibration, then F calibration draws; D runs in order
4, 20, 100 without resetting the scale RNG. The projection RNG is separate.

## Exact executable capture

```python
import numpy as np
from scipy.stats import chi2, f
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import run_projected_wald_kernel

# Analytic and implementation self-checks: squared norm and exact df=2 tails.
check = run_projected_wald_kernel(
    np.array([1., 2.]), spectral_k=2,
    pca_projection=np.eye(2), pca_eigenvalues=np.array([2., 1.]),
)
assert check.statistic == 5. and check.degrees_of_freedom == 2.
assert np.isclose(check.p_value, np.exp(-2.5))
assert np.isclose(f.sf(3., 2, 4), (1. + 2. * 3. / 4.) ** -2.)

rng = np.random.default_rng(2026110201)
N = 3000
n, d, k = 40, 8, 2
counts = dict(fixed=0, same_pca=0, energy_fixed_basis=0, same_pca_branch2=0)
for _ in range(N):
    x = rng.normal(size=(n, d))
    z = (x[:20].mean(0) - x[20:].mean(0)) / np.sqrt(.1)
    vals, vec = np.linalg.eigh(np.cov(x, rowvar=False))
    vals, rows = vals[::-1], vec[:, ::-1].T
    for key, basis, eigs, frac, zarg in [
        ('fixed', np.eye(d), np.arange(d, 0, -1), None, z),
        ('same_pca', rows, vals, None, z),
        ('energy_fixed_basis', np.eye(d), np.arange(d, 0, -1), .9, z),
        ('same_pca_branch2', rows, vals, None, z / np.sqrt(2)),
    ]:
        result = run_projected_wald_kernel(
            zarg, spectral_k=k, pca_projection=basis, pca_eigenvalues=eigs,
            adaptive_dimension_energy_fraction=frac,
        )
        counts[key] += result.p_value < .05
print('projection seed=2026110201 independent datasets=', N,
      'n=', n, 'd=', d, 'k=', k, 'tail counts', counts)

rng = np.random.default_rng(2026110202)
N, nu, c = 200000, 2, 3
for D in [4, 20, 100]:
    focal = c * rng.chisquare(nu, N)
    estimate = np.maximum(c * rng.chisquare(D, N) / D, 1)
    print('scale seed=2026110202 N=', N, 'df=', nu, 'true_c=', c, 'D=', D,
          'known=', np.count_nonzero(chi2.sf(focal / c, nu) < .05),
          'plugin_floor1=', np.count_nonzero(chi2.sf(focal / estimate, nu) < .05),
          'F_without_floor=', np.count_nonzero(
              f.sf(focal / (nu * (c * rng.chisquare(D, N) / D)), nu, D) < .05))
```

## Captured output

```text
projection seed=2026110201 independent datasets= 3000 n= 40 d= 8 k= 2 tail counts {'fixed': 132, 'same_pca': 467, 'energy_fixed_basis': 225, 'same_pca_branch2': 59}
scale seed=2026110202 N= 200000 df= 2 true_c= 3 D= 4 known= 9878 plugin_floor1= 26477 F_without_floor= 9896
scale seed=2026110202 N= 200000 df= 2 true_c= 3 D= 20 known= 10017 plugin_floor1= 14448 F_without_floor= 9931
scale seed=2026110202 N= 200000 df= 2 true_c= 3 D= 100 known= 10039 plugin_floor1= 10906 F_without_floor= 9976
```

## Interpretation and limits

Same-data all-leaf PCA and contrast-selected rank each produce excess lower-tail
p-values in this known-covariance fixed-split Gaussian control. Adding an arbitrary
branch variance factor can mask or reverse the excess. These counts do not
establish which factor dominates a production case or establish signal power.

The fitted-scale example demonstrates finite calibration uncertainty even with
independent observations and a valid common-scale family. Ordinary F handles
that restricted family; production selected weighted, overlapping stopping-event
records need a different justified reference law. No production repair follows
from substituting F here. Counts across methods are paired within each projection
dataset; do not treat their differences as independent binomial counts.
