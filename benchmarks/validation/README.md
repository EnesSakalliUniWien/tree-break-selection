# Validation Evidence Scaffolds

This directory contains strict scaffolds for tracking and generating validation
evidence that is not yet complete in the manuscript.

Programs are grouped by responsibility:

- `contracts/` defines report, provenance, and evidence-manifest contracts.
- `statistics/` validates covariance, selected tests, and error control; its
  `calibration_support.py` module owns the synthetic continuous covariance
  profiles.
- `sweeps/` explores alpha and method-family parameter grids.
- `manifests/` stores validation manifests; `contracts/method_constants_manifest.py`
  reads `method_constant_validation_manifest.json`, and the other two are records
  of retired runners (see below).

The method-constant manifest is deliberately strict:

- it enumerates every method constant currently requiring validation;
- it defines required output fields for future validation runs;
- it can record existing benchmark output paths as candidate sources;
- it does not read, summarize, or infer metrics from those paths;
- it marks evidence as `missing` until an explicit validation artifact supplies
  the required fields.

## Method Constants Covered

- `edge_alpha`
- `sibling_alpha`
- `mp_upper_edge_threshold`
- `min_spectral_dimension`
- `sibling_projection_dimension_rule`
- `empirical_null_weight_rule`
- `context_bandwidth_rule`
- `pass_through_traversal`

## Retired Runners

`statistics/feature_covariance_calibration.py` and
`statistics/selected_pca_projected_wald_calibration.py` were retired on
2026-10-17 (last present at commit 5c554bc5). The missing-evidence manifests
`manifests/feature_covariance_validation_manifest.json` and
`manifests/selected_pca_projected_wald_validation_manifest.json` remain as
records of the targets those runners covered.

## Usage

Run an alpha-grid benchmark diagnostic over the active TBS path:

```bash
uv run python -m benchmarks.validation.sweeps.alpha_grid_search \
  --suite full \
  --edge-alphas 0.0001,0.0003,0.001,0.003,0.01 \
  --sibling-alphas 0.001,0.003,0.01,0.03,0.1 \
  --output-dir benchmarks/results/alpha_grid_full
```

This diagnostic compares benchmark behavior across edge and sibling alpha
constants. It is not a proof of selected-tree Type-I error control.

Run a selected-edge Type-I geometry smoke over binary null cases:

```bash
uv run python -m benchmarks.validation.statistics.selected_edge_type1_geometry run \
  --suite binary \
  --case-names binary_2clusters \
  --modes fixed_tree,selected_tree \
  --edge-alphas 0.0001,0.001 \
  --sibling-alpha 0.01 \
  --replicates 5 \
  --base-seed 20260604 \
  --output-dir benchmarks/results/selected_edge_type1_geometry_smoke
```

This diagnostic separates fixed-tree edge calibration from same-data
selected-tree edge behavior. It is evidence for the selected-tree Type-I
problem and does not change production alpha defaults or add a calibration
fallback.

Run a traversal-aligned sibling-FDR smoke:

```bash
uv run python -m benchmarks.validation.statistics.traversal_sibling_fdr_null \
  --layers synthetic_valid_p \
  --case-names synthetic_balanced_binary_tree \
  --replicates 200 \
  --alpha 0.01 \
  --base-seed 20260604 \
  --output-dir benchmarks/results/traversal_sibling_fdr_smoke/synthetic
```

Run the binary fixed-tree, selected-tree, and inflated layers:

```bash
uv run python -m benchmarks.validation.statistics.traversal_sibling_fdr_null \
  --layers fixed_tree_wald,selected_tree_wald,selected_tree_inflated \
  --suite binary \
  --case-names binary_2clusters \
  --replicates 20 \
  --alpha 0.01 \
  --edge-alpha 0.001 \
  --base-seed 20260604 \
  --output-dir benchmarks/results/traversal_sibling_fdr_smoke/binary
```

This diagnostic separates algorithmic FDR behavior from fixed-tree
projected-Wald calibration, same-data selected-tree effects, and
empirical-inflation support failures. It does not add a production fallback or
change alpha defaults.

Create a manifest skeleton:

```bash
uv run python -m benchmarks.validation.contracts.method_constants_manifest create \
  benchmarks/results/04_generic_benchmark_runs/run_20260325_134017Z \
  --output benchmarks/validation/manifests/method_constant_validation_manifest.json
```

Validate a manifest:

```bash
uv run python -m benchmarks.validation.contracts.method_constants_manifest validate \
  benchmarks/validation/manifests/method_constant_validation_manifest.json
```

List constants:

```bash
uv run python -m benchmarks.validation.contracts.method_constants_manifest constants
```

## Evidence Policy

Existing benchmark outputs are recorded under `source_paths` only. They are not
used as evidence by this scaffold, because historical files may not contain the
locked seeds, configuration manifest, ablation grid, confidence intervals, and
constant-specific endpoints required for manuscript validation.

For each constant, the generated manifest contains:

- `evidence_status: "missing"`;
- `evidence.metrics: {}`;
- `evidence.source_path: null`;
- `evidence.missing_required_fields` listing every required field.

Future validation jobs should either extend this manifest with complete evidence
records or generate a separate results artifact that can be validated against the
same required field contract.
