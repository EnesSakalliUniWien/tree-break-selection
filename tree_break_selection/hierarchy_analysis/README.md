# hierarchy_analysis/

Decomposition entrypoints, gate orchestration, and cluster-assignment helpers for the
clustering pipeline.

## Current Layout

| Path | Purpose |
| ---- | ------- |
| `tree_decomposition.py` | Public decomposition engine. Prepares or accepts explicit gate annotation bundles, evaluates gates, and emits cluster assignments. |
| `cluster_assignments.py` | Builds cluster-root and per-sample assignment tables from explicit final tree boundaries. |
| `bootstrap_consensus.py` | Bootstrap-based stability helpers layered on top of decomposition results. |
| `decomposition/gates/orchestrator.py` | Runs the edge-divergence gate then the configured sibling-divergence gate and returns a node-indexed annotation bundle. |
| `decomposition/gates/profiles.py` | Named sibling-gate profile registry and profile-resolution contract. |
| `decomposition/gates/guards.py` | Optional root-stability and selected-permutation guard layers used by auditable profiles. |
| `decomposition/gates/gate_evaluator.py` | Encapsulates the binary-structure prerequisite, edge-divergence gate, sibling-divergence gate, and pass-through traversal decision as one action. |
| `decomposition/gates/column_contracts.py` | Shared checks for the gate-column contract carried through the explicit annotation bundle. |

## Decomposition Flow

1. `TreeDecomposition` receives an annotation DataFrame or an explicit `GateAnnotationBundle`.
2. `run_gate_annotation_pipeline()` populates edge-divergence and sibling-divergence columns when they are missing.
3. `GateEvaluator.decision()` applies the runtime gates and returns a boundary, split, or pass-through action.
4. `TreeDecomposition.decompose_tree()` advances the worklist, optionally honoring pass-through mode for deeper descendant splits.
5. `cluster_assignments.py` converts the final boundary nodes into cluster-level and sample-level outputs.

## Notes

- `TreeDecomposition` defaults to `trace_level="compact"`: results contain
  cluster assignments, independence metadata, and live traversal counters.
  Row-level `traversal_trace` and `full_edge_traversal_trace` payloads are
  produced only when diagnostics opt into `trace_level="full"`.
- The default sibling gate remains `projected_wald_inflation`. Opt-in
  `fixed_global_chi_square`, `fixed_coordinate_bh`, and `fixed_block_bh` gates
  avoid learning parent PCA projection rows and projection dimension from the
  same sibling contrast. `fixed_global_chi_square` is the full fixed-subspace
  chi-square reference; the BH variants aggregate predeclared coordinates or
  feature blocks. `sibling_gate_alpha_penalty` can be used to apply the
  selected-topology alpha penalty validated in diagnostics while keeping the
  public `sibling_alpha` value recorded in metadata.
  `root_stability_guard_threshold` and its subsampling settings expose the
  default-off selected-root stability guard; an unstable open root is closed
  fail-closed and recorded in `Root_Stability_*` diagnostic columns.
  `root_selective_permutation_guard_replicates` exposes a default-off
  validation diagnostic for fixed-subspace sibling gates. It is not part of the
  TBS runtime method: when explicitly enabled for a validation run, it records
  `Root_Selective_Permutation_*` diagnostics and can close an open root whose
  selected-root p-value is above the diagnostic guard alpha. The guard scope
  defaults to `root`; validation profiles can also request selected subtree
  checks for `open_internal` nodes or for the narrower `passthrough_descendant`
  contexts that are reachable only through an ordinary closed sibling ancestor.
  The `global_sibling_min_passthrough_descendant` scope evaluates those
  pass-through candidates against a whole selected-family minimum sibling
  p-value under feature-block permutations.
- `sibling_gate_profile="fixed_coordinate_guarded_v1"` packages the same-data
  fixed-coordinate repair as a named diagnostic candidate. The profile id is
  recorded in gate annotation metadata and participates in cache reuse checks.
- `sibling_gate_profile="fixed_coordinate_global_passthrough_refined_v1"` adds
  a higher-resolution replay for global pass-through families that land on the
  base Monte Carlo p-value floor. It is the current diagnostic follow-up after
  support-level binary runs exposed boundary false splits in the unrefined
  profile.
- Pure Bernoulli and pure categorical fixed-coordinate BH paths use exact
  vectorized Wald whitening before the existing coordinate-wise BH aggregation;
  larger support runs are still limited by repeated selected-family
  permutations.
- The current runtime centers on the gate annotation pipeline under `decomposition/gates/`.
- Statistical test details live in `statistics/README.md`.
