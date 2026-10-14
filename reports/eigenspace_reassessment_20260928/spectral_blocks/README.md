# Matched spectral-block reassessment — 2026-09-28

The complete same 11-case `method_proof` panel was rerun after the numerical
eigenvalue-group rank repair: 47 blocks, eight tree inference methods, both raw
and adaptive-diffusion geometries, and 752 attempts. The historical results in
`reports/eigenvalue_block_comparison_20260926/` were preserved.

## Matched outcomes

- Status counts before: `{'ok': 525, 'failed_tree': 170, 'unsupported': 57}`; after:
  `{'ok': 525, 'failed_tree': 170, 'unsupported': 57}`. There are 0 status changes.
- 525 runs returned partitions in both versions;
  33 of these partitions changed, including
  32 changes in cluster count.
- 582 trees are available in both versions;
  18 rooted clade sets changed.
- 22 changed partitions occur on
  unchanged rooted trees; 11 occur
  on changed trees. The former changes arise downstream of tree inference.
- Single-cluster outputs: 436 before and
  459 after. Among jointly completed comparisons,
  435 are single-cluster outputs in both
  versions; 90 have at least one
  multi-cluster output. Partition agreement in the former group is not evidence
  of useful cluster recovery.
- Against generator labels, ARI increased in
  15 and decreased in
  15 jointly completed runs. These dependent
  block/method counts are descriptive; no pooled quality score is reported.

The one-population `mp_spike_below_bbp_continuous` example had
10/36 completed outputs with multiple clusters
before and 0/36 after. This is one generated
dataset evaluated repeatedly, not an estimated null error rate or proof of
calibration. The phylogenetic Brownian example's generator labels have 16 groups
and are not recoded as a one-population partition truth.

### Explicit outcome transitions

| status_before | status_after | count |
| --- | --- | --- |
| failed_tree | failed_tree | 170 |
| ok | ok | 525 |
| unsupported | unsupported | 57 |

### Per-case and geometry counts

| case | geometry | ok_before | ok_after | changed_trees | changed_partitions | single_cluster_before | single_cluster_after |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bar_binary_balanced_4c | adaptive_diffusion | 18 | 18 | 0 | 0 | 17 | 17 |
| bar_binary_balanced_4c | raw | 20 | 20 | 0 | 0 | 15 | 15 |
| bar_binary_unbalanced_95_5 | adaptive_diffusion | 18 | 18 | 0 | 0 | 13 | 13 |
| bar_binary_unbalanced_95_5 | raw | 20 | 20 | 0 | 0 | 16 | 16 |
| cat_dirichlet_overdispersed_20cat | adaptive_diffusion | 29 | 29 | 0 | 0 | 26 | 26 |
| cat_dirichlet_overdispersed_20cat | raw | 32 | 32 | 0 | 0 | 27 | 27 |
| cat_simplex_face_rare_20cat | adaptive_diffusion | 48 | 48 | 0 | 0 | 40 | 40 |
| cat_simplex_face_rare_20cat | raw | 50 | 50 | 0 | 0 | 45 | 45 |
| cont_lowrank_pggn_shrinkage | adaptive_diffusion | 15 | 15 | 0 | 0 | 15 | 15 |
| cont_lowrank_pggn_shrinkage | raw | 17 | 17 | 0 | 0 | 17 | 17 |
| mp_spike_above_bbp_continuous | adaptive_diffusion | 18 | 18 | 0 | 3 | 15 | 18 |
| mp_spike_above_bbp_continuous | raw | 20 | 20 | 0 | 2 | 18 | 20 |
| mp_spike_below_bbp_continuous | adaptive_diffusion | 18 | 18 | 0 | 5 | 13 | 18 |
| mp_spike_below_bbp_continuous | raw | 18 | 18 | 0 | 5 | 13 | 18 |
| phylo_brownian_null_16taxa | adaptive_diffusion | 36 | 36 | 6 | 3 | 33 | 35 |
| phylo_brownian_null_16taxa | raw | 38 | 38 | 0 | 6 | 32 | 38 |
| support_selected_only_fail_closed | adaptive_diffusion | 17 | 17 | 6 | 3 | 13 | 12 |
| support_selected_only_fail_closed | raw | 19 | 19 | 0 | 0 | 13 | 13 |
| traversal_deep_branch_recovery_stress | adaptive_diffusion | 10 | 10 | 6 | 6 | 4 | 5 |
| traversal_deep_branch_recovery_stress | raw | 14 | 14 | 0 | 0 | 7 | 7 |
| traversal_deep_signal_under_same_parent | adaptive_diffusion | 24 | 24 | 0 | 0 | 22 | 22 |
| traversal_deep_signal_under_same_parent | raw | 26 | 26 | 0 | 0 | 22 | 22 |

## Geometry and provenance

Input case definitions, seeds, data hashes, retained spectra and all block records
match the historical output exactly. Package versions and all benchmark flags
also match. Direct replay compared all 94 geometries, using historical cosine and
diffusion source snapshots whose SHA-256 hashes exactly match the historical
manifest. All 47 raw geometries are exact; 44 of 47 adaptive-diffusion geometries
are exact. The three changed adaptive geometries are recorded in
`geometry_comparison.csv`. No component-count changes occurred in this panel.

The three changed geometries retain complete repeated leading groups: four
connected components for `support_selected_only_fail_closed` modes 2–6, two for
`phylo_brownian_null_16taxa` modes 2–5, and six for
`traversal_deep_branch_recovery_stress` modes 2–6. Their old diffusion selection
used modes 2–31; the repaired selection uses modes 1–30. The complete repeated
leading group is retained. See `changed_geometry_leading_modes.json`.

The source comparison records 232 paths; all remain unchanged during
the inference run. Differences from the historical source manifest are:

- `tree_break_selection/eigenvalue_groups.py`
- `tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/projection.py`
- `tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py`
- `tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_projection_basis.py`
- `tree_break_selection/hierarchy_analysis/statistics/projection/projection_dimension_estimation/projection_dimension_estimators.py`
- `tree_break_selection/space_separation/adaptive_cosine.py`
- `tree_break_selection/space_separation/diffusion.py`

The historical diffusion snapshot was recovered from its recorded Git commit.
The historical cosine module contained an uncommitted tie fix: the comparison
helper reverses the later shared-helper extraction, then requires an exact
historical SHA-256 match before executing that reconstructed snapshot. It never
replaces working-tree source files.

## Reproduction and verification

The actual invocation, environment overrides, start/end times and exit code are
in `invocation.json`; all runner output is retained in `run.log`. Inference wall
time was 915.9 seconds. Concurrent agent work
means this is not a controlled speed comparison. The existing benchmark command
creates its output directory and refuses to overwrite an earlier run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg uv run --no-sync python -m benchmarks.diagnostics.spectral.adaptive_cosine.compare_spectral_blocks --suite method_proof --max-rank 40 --methods average complete weighted single centroid median ward neighbor_joining --geometries raw adaptive_diffusion --output-dir <new-output-directory>
```

Recompute these comparisons against the preserved baseline with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg .venv/bin/python reports/eigenspace_reassessment_20260928/spectral_blocks/compare_before_after.py
```

The helper checks identical settings, input/spectrum/block records, unique
matched run keys, successful assignment coverage, truth labels, cluster counts,
recomputed ARI, paired tree availability and source stability. Failed and
unsupported runs retain their statuses and are excluded from partition metrics;
they are never turned into single-cluster assignments. The shared CLI calls
diffusion failures `failed_geometry`; none occurred in either version of this
panel. The local helper was linted successfully.

## Scope

The target rank remains 40, cosine weighting is raw values, minimum segment
length four and maximum segments eight. Diffusion uses K=15, time three,
30 target components, Euclidean metric, bandwidth `-1/(d+2)` and median epsilon.
The tree methods are average, complete, weighted, single, centroid, median,
Ward and neighbor joining with MAD rooting. Edge alpha is 0.001 and sibling
alpha is 0.01, with support-threshold enforcement false, branch-time policy
`none`, and contrast-energy adaptive Wald rank disabled in both runs.

This is a matched behavioral comparison, not held-out calibration evidence.
Failures, unsupported cases, one-cluster agreement and generator-label ARI must
remain separate. No inference about a method's general superiority or error
control follows from these dependent one-seed-per-case comparisons. IQ-TREE,
NNLS/branch-time variants and the full repository benchmark suite are outside
this panel.

Detailed changes: `run_comparison.csv`, `changed_runs.csv`,
`summary_by_case_geometry.csv`, `source_comparison.csv` and
`comparison_summary.json`. Original-format outputs (`runs.csv`, `trees.csv`,
`assignments.csv`, spectra, blocks and pair tables) remain alongside them.
