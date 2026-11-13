# The peel pass under the opt-in in production code

Implemented and run on 2026-11-13, Europe/Vienna. Local folder: `results/analyses/peel_pass_production_20261113`
(git-ignored); protocol, deviation, scripts and outputs are in `peel-pass-production-scripts-20261113.py.txt`.
Follow-up of `peel-pass-opt-in-20261113.md`. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5,
ruff 0.15.14, pytest 8.4.2. Nothing was committed.

## The production change

- `tree_break_selection/tree/topology.py`: `peel_continuations(tree, small_peel_leaf_limit=None)` replaces
  `singleton_peel_continuations(tree)`. Without a limit it returns the singleton peels as before. With a limit, a
  binary node whose smaller child has fewer leaves than the limit and whose larger child has at least the limit
  is a small peel and maps to its larger child.
- The limit is passed explicitly to the four users of the mapping: Tree-BH (`apply_tree_bh_correction`, through
  `annotate_child_parent_divergence`), the sibling calibration records (`collect_sibling_pair_records`), the
  sibling frontier (`apply_traversal_aligned_sibling_bh_results`) and `GateEvaluator` (attribute `peels`,
  formerly `singleton_peels`).
- `run_gate_annotation_pipeline` sets the limit to the reference's `MIN_CELLS` (15) when
  `sibling_best_flat_cut_reference=True` and to `None` otherwise, and records it as
  `GateAnnotationConfigMetadata.small_peel_leaf_limit`. `TreeDecomposition` reads it from the gate annotation
  bundle; with a bare annotation table it uses singleton peels only.
- `TreeDecomposition`: the trace reason of a small peel is `small_peel_continue`; a peeled small group whose
  edges are not significant joins the nearest cluster below the peel, as a peeled leaf does; a peeled small
  group that was split itself keeps its clusters.
- Reason for tying the limit to `MIN_CELLS`: the reference returns p = 1 for fewer than 15 cells, so a group
  below 15 cells can never be tested, like a single leaf. Reason for tying it to the opt-in: today's default
  sibling test over-splits one Gaussian population more when small peels pass (825 of 860 fits against 703).
- Tests added: `tests/tree/test_peel_continuations.py` (3), one each in
  `tests/statistics/41_test_multiple_testing_contracts.py`,
  `tests/statistics/43_test_traversal_aligned_sibling_fdr.py` and
  `tests/hierarchy_analysis/decomposition/gates/test_gate_evaluator.py`, two in
  `tests/hierarchy_analysis/test_tree_decomposition.py`, and assertions on the config in
  `tests/hierarchy_analysis/decomposition/gates/test_orchestrator.py` and
  `tests/integration/66_test_best_flat_cut_reference_runner.py`. They were written first and failed before the
  change.
- Documentation: `CHANGELOG.md`, `tree_break_selection/tree/README.md`,
  `tree_break_selection/hierarchy_analysis/README.md`,
  `tree_break_selection/hierarchy_analysis/statistics/README.md` and the module text of
  `best_flat_cut_reference.py`.

## Verification against the run-time trial (production code, nothing patched)

| Check | Result |
| --- | --- |
| CHOIR cells, default tree, opt-in | 432 clusters, ARI 0.6210, 42 lines recovered; the cluster of every cell equals the trial fit in 48,879 of 48,879 cells |
| One Gaussian population, 880 fits, opt-in | found clusters equal the trial's arm with the peel pass in 880 of 880 fits; 1 fit with more than 1 cluster |
| Benchmark, 14 cases, both registered methods, opt-in | found clusters and ARI equal the opt-in without the peel pass in 26 of 26 recorded rows; null case `gauss_null_256n_30f_continuous` at 1 cluster under both methods |

CHOIR cells, Markov/Ward tree, opt-in (current rule: table, k up to 15); there was no trial fit of this cell:

| Peel pass | Found clusters | ARI | NMI | Lines recovered | Split only | Merged only | Both | Clusters below 20 cells | Calls | Tested | Tested at the smallest p (0.0005) | Accepted by BH |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| no (11 November) | 230 | 0.8219 | 0.9358 | 116 | 21 | 38 | 15 | 25 | 1,838 | 1,475 | 140 | 134 |
| yes (production) | 223 | 0.8214 | 0.9334 | 110 | 20 | 43 | 17 | 23 | 2,210 | 1,689 | 144 | 130 |

Of the calls at the smallest p-value the table can give, 133 of 140 are accepted without the peel pass and 129
of 144 with it.

## Predictions

- Q1 held: the default-tree fit equals the trial fit cell by cell.
- Q2 held: 880 of 880 one-population fits equal the trial; 1 fit split.
- Q3 held: 26 of 26 benchmark rows equal; the null case at 1 cluster under both methods.
- Q4 not held: on the Markov/Ward tree the ARI is 0.8214 against 0.8219 and 110 lines are recovered against
  116. With the third version of the draft k the peel pass raised the ARI on this tree (0.8349 to 0.8411
  became 0.8642 to 0.8706).

## Reading

- Production code does what the run-time trial did.
- For the opt-in's current rule the peel pass is a large gain on the default tree (ARI 0.0275 to 0.6210) and a
  small loss on the Markov/Ward tree (6 fewer lines recovered, 5 more merged only).
- The loss on the Markov/Ward tree goes with the table's smallest p-value: the peel pass adds 372 calls, and
  more calls at 0.0005 are refused by Benjamini-Hochberg (15 of 144 against 7 of 140). The third version, whose
  smallest p-value is 0.00002, did not lose there. That the smallest p-value is the cause is a reading; no run
  with a finer table was made.

## Deviation

The copied one-population script imported the trial module, which imports the removed name
`singleton_peel_continuations`; the five processes stopped with ImportError before any fit. The import and the
trial's arm were removed from the copy and the fits were run again. The scripts of the earlier trial folders
that patch `singleton_peel_continuations` no longer run against the current code.

## Follow-up on the same day: the canonical full benchmark

`uv run --no-sync --offline python -m benchmarks.full.run` with the default settings, output in
`benchmarks/results/run_20261113_peel_pass_full/` (git-ignored). 127 cases and the 11 default methods, 1,397
rows (1,353 ok, 44 skipped as unsupported), exit code 0. Compared with the full run of 10 November
(`benchmarks/results/run_20261110_122947Z_full/`, 126 cases):

- All 1,386 paired rows have the same status, the same found clusters and the same ARI (within 1e-9), for every
  method. Mean ARI of the TBS methods: `tbs` 0.7847, `tbs_family_diffusion` 0.8279,
  `tbs_diffusion_adaptive_nnls` 0.8231 (101 supported cases), `tbs_diffusion` 0.5380 (109 supported cases).
- The one new case is `gauss_null_256n_30f_continuous`: `tbs` 1 cluster, `tbs_family_diffusion` 10 clusters
  (ARI 0.0), as recorded when the case was added.
- The default method set does not use `sibling_best_flat_cut_reference`, so this run checks only that the
  defaults are unchanged by the work since 10 November.

## Follow-up on the same day: where the loss on the Markov/Ward tree comes from (`markov_ward_calls.py`)

The two fits of the Markov/Ward tree with the opt-in's current rule, call by call:

- 1,819 calls are in both fits, and every one has the same p-value and the same reason in both.
- 7 nodes that are accepted without the peel pass are refused with it. All 7 have p = 0.0005, the smallest
  p-value the table can give; they have 84, 246, 454, 517, 540, 1,004 and 1,022 cells. No node goes the other
  way.
- The peel pass adds 391 calls and removes 19 (none of the 19 was accepted). Of the 391, 380 have p = 1 (170
  without invariant directions), 160 have fewer than 30 cells, 4 are at 0.0005 and 3 are accepted.
- Accepted: 134 - 7 + 3 = 130.

Reading: the statistic does not change; the multiple-testing step does. Benjamini-Hochberg at alpha 0.01 can
reject a p-value of 0.0005 only in a family where at least 0.0005 / 0.01 = 5% of the calls are at or below it,
because the table cannot rank calls below 0.0005. The newly reached nodes are small and nearly all have p = 1;
they enlarge the families, and in some of them the calls at 0.0005 fall below that share. The families
themselves were not reconstructed, so the share per family is not measured.

## Not checked

- Row orders do not enter the current rule; no other real dataset; limits other than 15 cells.
- The opt-in with guards or with fixed-subspace sibling methods (they exclude the opt-in).
- `collect_sibling_pair_records` with a limit has no test of its own; it is covered by the fits above.
- Populations that are not Gaussian. A table with more draws.
