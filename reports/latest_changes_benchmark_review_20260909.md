# Latest changes and benchmark review

Reviewed on 2026-09-09: commit `6441b2e6` against
`641c0e7a9145be3b4d18f147af95042a7c5e7fc5`, the three uncommitted
`selected_gaussian_*` modules and their tests/documentation, and the two latest
2026-08-11 benchmark directories. No production code was changed or full
benchmark rerun. The four findings below reproduce in the current checkout.

## Findings

### P2: Hierarchy interval certification accepts an incorrect merge sequence

Source: `tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_hierarchy.py:211`,
with polynomial replay at lines 97-98 and 139-142.

For `e=1e-8`, `a=1-0.2*e`, `X=[[a],[-a],[-1],[-1+e]]`, contrast
`[1,-1,0,0]`, covariance `[[1]]`, and dimension 1, the observed first merge is
leaves 1/2. With `t=r/sqrt(2)`, their distance is `abs(1-t)`, leaves 1/3 have
distance `abs(1-e-t)`, and leaves 2/3 have distance `e`. Comparing these gives
the full observed merge-sequence interval `sqrt(2)*[1-e/2,1+e]`; the later
merge order adds no boundary here.

The constructor returns `[1.4142135555282114,1.4142135834465193]` instead of
approximately `[1.4142135553020274,1.4142135765152306]`. It gives selected
p-value **0.8561353965**, versus **0.8000000164** from integrating the chi density
on the analytic interval. At radius `1.41421358`, within the returned interval,
SciPy linkage on reconstructed squared distances first merges leaves 2/3;
the polynomial replay reports leaves 1/2.

The certification reuses the same unstable expanded-polynomial evaluation.
Resolve coefficient/root cancellation with sufficient precision or fail on
unresolved boundaries, and verify reconstructed distances independently.
Existing narrow-interval tests pass despite this counterexample. This affects
the internal prototype, which is not connected to production gates.

### P2: Production annotation bypasses support-threshold decisions

Source: `tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflated_projected_wald_annotation/pipeline.py:191`.

The unresolved-law branch bypasses `decide_empirical_null_calibration` for
positive-dimensional focal records. On the existing
`tests/statistics/50_test_sibling_skip_annotation.py` fixture with
`enforce_support_thresholds=True`, the direct decision returns
`undefined_sparse_context` with six threshold failures. The pipeline instead
stamps `undefined_unvalidated_reference_law` regardless of that setting.

It discards the fitted model on return: no fitted scale, dependency-group IDs
or threshold-failure diagnostics survive. Null calibration parent N3 also has
no raw statistic/p-value in its output row. This conflicts with the diagnostic
preservation claim in `wiki/analyses/empirical-null-calibration-reference-law-contract.md`.
Retain the diagnostic decision before withholding calibrated p-values. Both
current outcomes fail closed; this is a diagnostic/enforcement issue, not an
accepted invalid tail.

### P2: Tuple-valued dependency groups crash the decision API

Source: `tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py:117`.

`CalibrationSample` accepts hashable IDs, but converting groups `('tree',0)`
and `('tree',1)` with `np.asarray(...,dtype=object)` creates a two-dimensional
array. Two valid calibration records fit successfully; their decision then
raises `IndexError: too many indices for array` when masking the weight vector.
Construct a one-dimensional mask with scalar object comparisons.

### P2: Chi-tail underflow rejects a positive selected event

Source: `tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_radial.py:133`.

With observed radius `40.01`, dimension 1 and selection interval `[40,inf)`,
absolute chi probabilities underflow and trigger the positive-probability
validation error. The conditional p-value is well defined:
`exp(log_ndtr(-40.01)-log_ndtr(-40)) = 0.6701192098`.
Use scaled/log-domain probabilities to distinguish underflow from impossibility.

## Benchmark assessment

Primary inputs are the `full_benchmark_comparison.csv` files under:

- `benchmarks/results/run_20260811_123017Z_tbs_full_quality/`
- `benchmarks/results/run_20260811_selected_law_variants_full/`

The variants CSV contains **2,790 unique case/run cells, 25 methods, 122 cases**:
1,596 `ok`, 1,139 `unsupported`, 55 `skip`. Twenty-two methods have 122 rows;
bandwidth-context and neighbor-joining have 52 each, and IQ-TREE has two.
No invocation manifest establishes whether those shorter panels are complete.

**Its saved summary describes only a subset.** The performance report cites
the full CSV but summarizes only 732 rows and six NNLS methods. Its summary CSV
and metric/status grids also omit the other 19 methods. The actual best fully
covered method is consequently missing from that ranking. The accompanying
`latest_changes_benchmark_review_20260909_methods.csv` summarizes all 25 methods
without rewriting the original artifacts.

| Method | OK / rows | Unsupported | Skip | Mean ARI on OK | Exact K on OK |
| --- | ---: | ---: | ---: | ---: | ---: |
| `tbs` | 11 / 122 | 111 | 0 | 0.3636 | 4 / 11 |
| `tbs_nnls` | 13 / 122 | 109 | 0 | 0.3077 | 4 / 13 |
| `tbs_diffusion` | 2 / 122 | 107 | 13 | 1.0000 | 2 / 2 |
| `tbs_diffusion_adaptive_nnls` | 5 / 122 | 96 | 21 | 0.2000 | 1 / 5 |
| `tbs_fixed_coordinate_by` | 122 / 122 | 0 | 0 | 0.7123 | 30 / 122 |
| `tbs_fixed_coordinate_by_nnls` | 122 / 122 | 0 | 0 | 0.6990 | 28 / 122 |
| `tbs_fixed_coordinate_bh` | 122 / 122 | 0 | 0 | 0.6732 | 24 / 122 |
| `tbs_fixed_coordinate_bh_nnls` | 122 / 122 | 0 | 0 | 0.6627 | 24 / 122 |

The standalone canonical run agrees with the variants run's canonical status
and quality. Of 111 unsupported cases, 85 lack a validated reference law and
26 lack internal support. **All 11 successful canonical cases return one
cluster**: four match one-cluster truth and seven under-split. Mean ARI 0.3636
therefore does not describe full-suite performance. Diffusion's ARI 1.0 comes
from only two successful cases.

Fixed-coordinate/block methods are recorded as `diagnostic`. They use
fixed-coordinate/block chi-square tails on the selected hierarchy; complete
execution does not validate the selected-tail law. BY over-splits 74 cases,
under-splits 18 and recovers exact K in 30; mean absolute cluster-count error
is 7.2049. It returns three clusters on the one-cluster `gauss_null_large`
case. BH returns four there and two on `gauss_null_small`. These are individual
null outcomes, not Type I error estimates.

### Matched NNLS comparisons

For each case, recorded seeds, sample/feature counts, true K, noise and feature
representations match across methods. Case seeds vary; each case has repeat
index 0 only. All rows use sibling alpha 0.01 and edge alpha 0.001.

| Method | Matched cases | NNLS ARI wins / ties / losses | Mean ARI change |
| --- | ---: | ---: | ---: |
| Coordinate BY | 122 | 10 / 106 / 6 | -0.013298 |
| Coordinate BH | 122 | 14 / 101 / 7 | -0.010570 |
| Coordinate Bonferroni | 122 | 14 / 101 / 7 | -0.010595 |
| Block BH | 122 | 15 / 104 / 3 | +0.005883 |

NNLS does not uniformly improve these variants. Multiple realizations and an
alpha sweep are absent. The result schema has no diffusion-construction or
branch-length-optimization timings, so it cannot establish total NNLS overhead.
These results do not evaluate the new selected-Gaussian prototype.

### Other evidence limits

- No duplicate case/run cells, non-finite ARI on `ok` rows, or finite ARI on
  unsupported/skipped rows were found.
- Hamming/native-continuous incompatibility and duplicate-heavy adaptive
  diffusion skips are distinct from calibration failures.
- The latest variants run contains no non-TBS comparator methods.
- Neither latest directory contains an exact code-revision manifest.
- Both failure reports show missing audit logs for all seven listed cases;
  their requested run-local audit directories are absent.
- Holm and Bonferroni have identical outputs: the minimum Holm-adjusted
  p-value equals the Bonferroni minimum in the inspected implementation.
  Coordinate BH and block-Simes BH also coincide for singleton blocks.
  These rows are not independent confirmations of performance.

## Verification

Inspected manifests, Make targets, test configuration, fixtures and relevant
callers. Verified Darwin arm64, zsh 5.9, uv 0.9.5, repository Python 3.11.12,
pytest 8.4.2, NumPy 2.3.4, SciPy 1.16.3, pandas 2.3.3,
scikit-learn 1.7.2, pydiffmap 0.2.0.1 and Ruff 0.15.14.

- 99 targeted calibration, prototype and unsupported-outcome tests passed in
  33.49 seconds; 35 related integration/report/wiki tests passed in 7.85 seconds.
- Ruff, baseline 190-page and final 191-page wiki lint, and `git diff --check`
  passed. The saved method summary agrees with all 2,790 source rows.
- The four counterexamples above reproduced separately despite passing tests.
- Full `make check` and full benchmark reruns were not performed.

## Numerical reproduction

Run from the repository root using `.venv/bin/python`:

```python
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from scipy.special import log_ndtr
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_radial import (
    build_gaussian_selected_projection_path, compute_selected_gaussian_tail,
    RadialSelectionInterval,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_hierarchy import (
    construct_squared_euclidean_average_linkage_region,
    replay_squared_euclidean_average_linkage,
)
e = 1e-8
a = 1 - 0.2 * e
path = build_gaussian_selected_projection_path(
    np.array([[a], [-a], [-1], [-1 + e]]),
    np.array([1., -1., 0., 0.]), np.eye(1), projection_dimension=1,
)
region = construct_squared_euclidean_average_linkage_region(path)
print(region.selection_intervals)
print(compute_selected_gaussian_tail(
    observed_radius=path.observed_radius, projection_dimension=1,
    selection_intervals=region.selection_intervals,
).p_value)
probe = 1.41421358
print(replay_squared_euclidean_average_linkage(path, radius=probe)[0])
print(linkage(pdist(path.reconstruct(probe), metric='sqeuclidean'),
              method='average')[0, :2])
print(np.exp(log_ndtr(-40.01) - log_ndtr(-40)))
compute_selected_gaussian_tail(
    observed_radius=40.01, projection_dimension=1,
    selection_intervals=(RadialSelectionInterval(40, np.inf),),
)  # Raises despite the well-defined conditional p-value.
```

Tuple-group reproduction:

```python
import runpy
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    fit_empirical_null_inflation_model, decide_empirical_null_calibration,
)
make = runpy.run_path('tests/statistics/35_test_empirical_null_inflation_estimation.py')['_make_record']
records = [make(f'n{i}', stat=4., degrees_of_freedom=2.,
                calibration_dependency_group=('tree', i)) for i in range(2)]
model = fit_empirical_null_inflation_model(records)
focal = make('f', stat=5., degrees_of_freedom=2.,
             is_null_like=False, sibling_null_weight=0.)
decide_empirical_null_calibration(model, focal)  # Raises IndexError.
```
