# The 13 benchmark cases with the best-flat-cut reference extended to 15 invariant directions

Run on 2026-11-12, Europe/Vienna. Local folder:
`results/analyses/best_flat_cut_reference_benchmark_k15_20261112` (git-ignored); protocol, deviations, the three
scripts, the comparison output and the 52 result rows are in
`best-flat-cut-reference-benchmark-k15-scripts-20261112.py.txt`. The opt-in `sibling_best_flat_cut_reference`
and its extension (`MAX_DIRECTIONS` = 15, `MIN_CELLS_ABOVE_10_DIRECTIONS` = 50, table with k = 1 to 15) were
uncommitted working-tree code on top of `bc044e07`. No production code was changed by this run.

## Setup

`run.py` and `reasons.py` are byte copies of the 11 November scripts
(`best-flat-cut-reference-benchmark-scripts-20261111.py.txt`): the 9 cases of
`get_test_cases_by_suite("continuous")` and the 4 cases of `get_test_cases_by_suite("cancer_pc")` with
`tbs_family_diffusion` and `tbs_family_diffusion_markov_ward`, registered cell ("today") and the same cell with
`sibling_best_flat_cut_reference=True` ("reference"), 52 runs (25 s), then a replay of the 26 reference runs
that records every call of `best_flat_cut_reference`. `compare.py` compares both outputs with the 11 November
`results.csv` and `reasons.json`. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

The protocol was written before the run and predicted no change: on 11 November the 151 recorded calls had
k = 0 (58 calls) or k = 1 to 5 (93 calls), and the table cells for k up to 10 are bit-identical.

## Result

No change.

- 52 of 52 runs finished ok. Found clusters, ARI, NMI and status equal the 11 November run in 52 of 52 rows
  (largest ARI difference 0.0).
- 151 recorded calls, as on 11 November: 93 `tested`, 58 `no_invariant_directions`; k = 0 in 58 calls, 1 in
  35, 2 in 38, 3 in 10, 4 in 8, 5 in 2. No call has more than 5 invariant directions, so the new cells for
  k = 11 to 15 are never read. Subtree size, left size, reason, k and p-value are equal in 151 of 151 calls.
- `rho` differs between the two runs in 80 calls by at most 1.0e-12 relative (floating point).
- The replay's found clusters and ARI equal the first run in 26 of 26.

Mean ARI (unchanged): continuous, default tree, today 0.8802, reference 0.6950; continuous, Markov/Ward, today
0.8525, reference 0.6950; cancer_pc, default tree, today 0.8218, reference 0.7221; cancer_pc, Markov/Ward,
today 0.3033, reference 0.7221.

## Reading

The extension changes only subtrees with 11 to 15 invariant directions. These benchmark cases have none, so
the losses reported on 11 November stay as they were: subtrees with no invariant directions (30 cells in 40
features, 30 in 30, 40 in 20,000, dim_diffuse_6c_136f) and the two connected cancer_pc cases (4 clusters
instead of 6). The run is a regression check of the extension, not new evidence on the method.

## Deviation

`compare.py` was changed after its first output, which counted the `rho` round-off as differing calls and read
the CSV with a parser that does not round-trip floats. `run.py`, `reasons.py` and the module were not changed.

## Not checked

No benchmark case with 11 to 15 invariant directions exists in these two suites, so the panel says nothing
about the new cells. The binary and categorical cases were not run (the reference is derived for Gaussian
coordinates). One run.
