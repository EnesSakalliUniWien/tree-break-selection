# Passthrough traversal benchmark

This reruns the existing fixed-tree Gaussian diagnostic on 1,000 matched datasets for each of five patterns and two covariance profiles (identity, AR(1) 0.6). The same dataset and exact known-covariance sibling tests feed both `passthrough=False` and `True`; all edge prerequisites are deliberately open. Each decision is from production traversal-aligned sibling BH and `TreeDecomposition`.

Here, `pass_through` bypasses a closed sibling gate only. Binary and edge prerequisites must be open, a descendant split must be available, and any optional support guard must be open. It is recorded separately from `split` and never counted as a false split. For example, at R, A’s child populations `+mu` and `-mu` average to zero, as do B’s two zero populations. R is null, while A is a signal; pass-through can continue from R to test A without counting R as a split.

| Pattern / profile | Signal node split power, disabled → enabled | Any reached false split, disabled → enabled | Mean reached false-split FDP, disabled → enabled |
|---|---:|---:|---:|
| Global null / identity | — | 12/1000 → 20/1000 | 0.012 → 0.020 |
| Global null / AR(1) | — | 15/1000 → 23/1000 | 0.015 → 0.023 |
| A deep only / identity | A: 11/1000 → 1000/1000 | 11/1000 → 26/1000 | 0.0055 → 0.0130 |
| A deep only / AR(1) | A: 12/1000 → 1000/1000 | 12/1000 → 19/1000 | 0.0060 → 0.0095 |
| Separate deep branches / identity | A,B: 9/1000 → 1000/1000 each | 9/1000 → 9/1000 | 0.0030 → 0.0030 |
| Separate deep branches / AR(1) | A,B: 13/1000 → 1000/1000 each | 13/1000 → 13/1000 | 0.0043 → 0.0043 |

The A-only and separate-branch patterns show why passthrough helps: a null root otherwise blocks strong descendant signals. In the separate-branch cases, the final partition count rises from mostly 1 to 4 for every enabled run; in A-only it rises from mostly 1 to mostly 3. The global-null increase is 8 additional datasets with a false split per 1,000 in each profile. These rates have per-cell Wilson 95% intervals in the summary CSVs; paired upward/downward false-split changes are also recorded.

This is a bounded fixed-tree, known-covariance Gaussian study with deliberately open edge prerequisites. It does not establish full selected-tree calibration or production-data behavior. The run exactly reproduces the earlier seed-20260928 output, so the rerun confirms reproducibility and adds no independent simulation evidence. The baseline diagnostic tests passed (4 targeted tests; root also reports 22 tests across its broader baseline). Pre-run source hashes, configuration, command, and environment are in [`pre_run_manifest.json`](pre_run_manifest.json); run outputs and status counts are in [`evaluation/`](evaluation/). Root’s independent traversal reconstruction and prior-run equality check passed in [`independent_verification.json`](independent_verification.json).

Full per-node reach/split power, pass-through counts, and final partition frequencies are in:

- [`passthrough_signal_node_summary.csv`](evaluation/passthrough_signal_node_summary.csv)
- [`passthrough_pattern_summary.csv`](evaluation/passthrough_pattern_summary.csv)
- [`passthrough_partition_summary.csv`](evaluation/passthrough_partition_summary.csv)

The output statuses were `ok` for all 120,000 local-test rows, 10,000 selected-root rows, and 20,000 traversal rows; there were no unavailable or runtime-error outcomes.
