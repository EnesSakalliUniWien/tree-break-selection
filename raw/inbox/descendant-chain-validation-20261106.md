# Descendant chain family test: implementation and validation

Work on 2026-11-06 on `dev` at `3606bfd0` plus the uncommitted descendant test.

## Test

`selected_diffusion_family_chain_tests` in
`benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py`
tests a pre-set family: the root and each internal child of the root. Each
node's mean contrast between its two children uses the random-position
reversible chi chain on that contrast's radius; only the node's rows move
(`_node_path`, adapted from `descendant_path` in
`conditional-small-benchmarks-20261104.py.txt`). A proposal is in the event when
every split from the root down to the node recurs (`_ancestor_path_pairs`),
which is what places the node in the family. At the root this is the validated
root-split event, and the family's root position reproduces
`selected_diffusion_root_chain_test` exactly for the same seed. Decisions use a
fixed Bonferroni level `alpha / 3`; position `j` uses chain seed `seed + j`.
Model: independent Gaussian rows with known covariance `sigma**2 * I`.

The chain runs joblib sequentially: merge signatures are identical, and each
replay drops from about 80 ms to 8 to 20 ms for 30 to 200 rows, because
low-dimensional neighbor searches otherwise start workers on every call.

Six tests in `tests/statistics/69_test_diffusion_descendant_chain.py` failed
before the implementation and pass after it; the full ordered suite passed
1,051 tests. A first strong-signal test at 16 rows failed because the diffusion
tree split each child 2/6 across the true subgroups; the test now uses 32 rows,
where the children split along them.

## Pre-set validation

`descendant-chain-validation-protocol-20261106.md` was fixed before any
panel result, with one implementation addendum (seed offsets) also added before
running. 390 families, sigma 1, 999 steps, alpha 0.05, seven processes, 665 s.
Per-dataset rows: `descendant-chain-validation-results-20261106.jsonl`
(SHA-256 56e8654608af801c6e320265faf80d5cc758216d5ebe4784499010ae252b91fd);
summary: `descendant-chain-validation-summary-20261106.csv`.

| Panel | Families | False rejections at null nodes | Rejections at non-null nodes |
| --- | ---: | --- | --- |
| A global null, 30 rows | 200 | root 4/200, children 9/400; any position 13/200 (upper 95% bound 0.1014) | — |
| B root signal, null children, 30 rows | 100 | children 1/141 | root 70/100; mixed children 2/59 |
| C 2 x 2 hierarchy, 40 rows | 50 | none | root 49/50; children 5/100 |
| D nulls and 2 x 2 hierarchies, 120 and 200 rows | 40 | 0/60 null nodes | roots 20/20; children 6/40 |

- No family or node was unresolved; one chain stalled.
- Across 801 null nodes, the share of p-values at most 0.05, 0.0167 and 0.01 is
  0.0537, 0.0175 and 0.0112.
- The family-wise rate 13/200 is consistent with 0.05 but not established below
  it: the one-sided 95% upper bound is 0.1014.
- In every one of the 140 hierarchical children (C and D) the observed child
  split mixed the true subgroups; its smaller side had a median of 4 of 20
  rows. Low child power reflects these peel-shaped child splits, not the test.

## Benchmark cases

`descendant-chain-benchmark-script.py.txt` ran the family on the six registered
cases with independent Gaussian rows and known spherical noise and unchanged
raw inputs: the four `blobs_continuous` cases and the two no-hub `cancer_pc`
surrogates. The other 120 cases are binary or categorical, median-thresholded,
non-spherical, outlier mixtures, phylogenetically dependent, graphs, or hub
surrogates with a different hub-row noise scale. Results:
`descendant-chain-benchmark-results-20261106.jsonl`.

- All six single-population nodes are kept (p 0.47 to 0.89), including all
  three positions of `gauss_null_large_continuous`.
- Eleven of twelve multi-population nodes are rejected. The miss is the root of
  `gauss_dense_signal_highd_continuous` (20,000 features; p 0.03, 14 chain
  moves).
- Runtime per case was 13 to 153 s; the 496-row surrogates took about 90 s.
- Production `tbs_family_diffusion` recovers the true K on five of these cases
  and ARI 0.697 on `cancer_pc_nohub_connected`; the family test agrees and
  finds no structure production missed.

## Limits

Known-noise Gaussian data only; root and immediate children only; no
production traversal, pass-through, counts, fitted PCA or unknown covariance.
