# Native diffusion positive-self-distance artifact versus the exact-real diagnostic

Date: 2026-11-02. Diagnostic investigation only; no production change.

## Finding

The `n<=32` kernel enclosure cap was computational, not mathematical. It has been replaced by the native dense-spectral-route limit (`DENSE_DIFFUSION_EIGENSOLVER_MAX_SAMPLES`, currently2000). The actual first80-by24 Gaussian attribution pair has a successful exact-real operator enclosure on the standardized-radius cell `r_obs ± 1e-12`. One timed execution took4.997seconds, with BGH index37, internal epsilon1/8, rounded intrinsic dimension2, target diffusion rank24 and maximum operator-entry enclosure width approximately2.296e-14.

However, native floating-point midpoint agreement exposed a material difference from the exact-real formulas. Scikit-learn's brute-force Euclidean nearest-neighbor query can produce tiny positive self distances. The native outer epsilon uses the median of ALL strictly positive sparse scaled-distance squares, so these self artifacts enter epsilon even though the coordinate operator later zeroes the weight diagonal.

For the original first attribution matrix and its reconstructed radial observation:

| Input passed to native geometry | Positive self entries | Zero self entries | All positive sparse entries | Native positive-entry median | Off-diagonal positive-entry median |
|---|---:|---:|---:|---:|---:|
| Original dyadic selected matrix |29|51|1149|0.8588652016711084|0.8639359788816151|
| Floating reconstruction at observed radius |28|52|1148|0.859336163375187|0.8639359788816152|

The path's observed standardized radius is8.765613678802739. The native complete merge signature is unchanged by this particular reconstruction (the path constructor verifies it), but the positive-self mask and epsilon change. Thus even tiny reconstruction roundoff can have a discrete effect on the positive-entry median.

The exact-real formulas have exactly zero self distances, so1120 off-diagonal positive entries for this directed15-neighbor graph (14 nonself edges per row). Reconstructing the zero-self exact-real operator from native row-dependent bandwidths agrees with the rigorous enclosure within the separately stated native-rounding tolerance1e-14. The literal native operator differs materially; this is not repaired by inflating that tolerance. Tests record the literal mismatch separately.

A floating eigendecomposition/linkage proposal using the zero-self operator changes the complete history beginning at zero-based merge index5 compared with the native observed history. This proposal is not a proof of the exact-real history. The scanner's interval-based `observed_event_validated` check must establish observed membership before any conditional p-value decision. If the observed native history is not selected by the declared exact-real constructor at its exact observed radius, the diagnostic must remain undefined/retain one cluster.

## What was checked

- The reduced density/bandwidth calculation matches the native final bandwidths up to one shared normalization factor across all80 rows: approximately3.348131185123405 (row variation at floating rounding scale).
- The seven-neighbor pilot bandwidths agree at approximately9e-16 absolute error.
- Internal BGH epsilon/dimension agree: native epsilon0.12500000000000003 versus exact-real1/8, dimension2.
- The discrepancy is the outer median's inclusion of tiny positive self artifacts, rather than a missing KDE row factor or a neighbor-count formula error.
- `DiffusionReplayPoint.state` now records positive/zero self counts, total positive-entry count, actual positive-square median and off-diagonal positive-square median. These fields are point observations, not certificates.
- The rigorous kernel specifies exact-zero self distances and rational bandwidth exponent `-1/(d+2)`. It does not claim bitwise NumPy/NumExpr agreement.

## Reproduction

Run from the repository root using the existing environment:

```bash
uv run --no-sync python - <<'PY'
from benchmarks.diagnostics.calibration.diffusion_attribution import gaussian_tree_pairs
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path, _geometry,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_kernel import (
    enclose_diffusion_operator,
)
pair = next(gaussian_tree_pairs(replicates=1))
path = make_diffusion_root_path(pair.selected_data)
for label, state in [
    ('original', _geometry(path.values, path.labels)[1]),
    ('reconstructed', path.evaluate(path.observed_radius).state),
]:
    print(label, {key: state[key] for key in [
        'positive_self_scaled_distance_count',
        'zero_self_scaled_distance_count',
        'positive_scaled_distance_count',
        'actual_positive_scaled_sq_median',
        'off_diagonal_positive_scaled_sq_median',
    ]})
r = path.observed_radius
operator = enclose_diffusion_operator(path, r-1e-12, r+1e-12)
print('BGH', operator.proof_details['kde_bgh_index'])
print('dimension', operator.proof_details['kde_dimension'])
print('target', operator.retained_count_target)
print('max entry width', max(float(upper-lower)
      for lower_row, upper_row in zip(operator.lower, operator.upper)
      for lower, upper in zip(lower_row, upper_row)))
PY
```

The80-row regression in `tests/statistics/64_test_diffusion_root_kernel.py` independently extracts graph self counts/medians, checks the exact-real zero-self operator against the enclosure, and records a separate literal-native discrepancy when positive self artifacts are present. It does not change production epsilon or claim literal floating certification. Counts may vary across versions, numerical kernels and reconstructed coordinates; package/source provenance must accompany the result.

## Limits

No production fix was authorized or applied. General bitwise floating-point selective inference is not established. The exact-real conditional diagnostic and the native production map are materially different on this fixture; this limitation must remain visible in any report or panel. Kernel enclosure support at80 rows does not establish a certified selected event or an informative conditional p bound for that native observed history.

## Subsequent final interval evidence

The final source-verified panel in
`raw/assets/diffusion-root-diagnostic-20261102/final-attribution-first-root/`
strengthens the earlier floating proposal. The interval
`[8.765613678801738, 8.765613678803739]` contains the entire rational enclosure
of the exact observed radius. Its certified label is `nonselected` for the
native observed complete history under the declared exact-real constructor.
At zero-based merge step 5, the chosen pair `[9,72]` loses to competitor
`[22,66]`: the upper bound on competitor-minus-chosen distance is
`-5.551285682872848044498961817444603630489e-7`.

This is a whole-cell operator/spectral/linkage certificate excluding the
native history at the exact observation, not just the earlier point proposal.
It does not identify the entire alternative history. The root gate consequently
does not use this different event for a selected p-value; it returns [0,1],
unresolved, and one cluster. All recorded source hashes matched at final audit.
