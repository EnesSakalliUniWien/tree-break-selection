# Full benchmark after the singleton-peel merge and the L-BFGS-B NNLS switch

Run on 2026-11-06 with the default command and settings
(`uv run --no-sync --offline python -m benchmarks.full.run`, plots on) on merge
`50455316` plus uncommitted working-tree changes, including the column-scaled
L-BFGS-B NNLS fit with its 5,000-iteration budget. Output:
`benchmarks/results/run_20261106_lbfgsb_full/` (local, ignored); its
`full_benchmark_comparison.csv` has SHA-256
4d257dbc330bdb320bbc72b5741b85b8be7797f947eeba7da54ef0794697b512.
The run completed in 608 s with exit code 0.

## Scope

The run has 126 cases and 1,386 rows: the 122 earlier cases plus the four
uncommitted `cancer_pc` cases. `full-benchmark-lbfgsb-compare-script.py.txt`
matches the 1,342 common rows against the "after" results of
`reports/matched_full_benchmark_20261103/` (self-distance fix and inflation cap,
before the singleton-peel merge and the solver switch). It is an integrated
comparison of those two changes, not an isolated estimate of either.

## Results on the 122 common cases

| Method | Mean ARI | Exact K | Improved / degraded | Null splits |
| --- | --- | --- | --- | --- |
| `tbs` | 0.805953 → 0.810390 | 93 → 94 | 2 / 1 | 0 → 0 |
| `tbs_diffusion` | 0.537087 → 0.537957 | 41 → 40 | 1 / 2 | 0 → 0 |
| `tbs_diffusion_adaptive_nnls` | 0.823085 → 0.823085 | 84 → 84 | 0 / 0 | 0 → 0 |
| `tbs_family_diffusion` | 0.825843 → 0.828107 | 98 → 98 | 1 / 1 | 0 → 0 |

The seven non-TBS methods are identical, and no row changed status: every NNLS
fit passed the new duality-gap check and was applied. The eight changed rows
are `phylo_conserved_4taxa`, `cat_clear_3cat_4c` and `gauss_outlier_cluster_4c`
(the cases the singleton-peel change reports, including its false 38 + 2
split, ARI 1.0 → 0.984) and three `phylo_*` cases under `tbs_diffusion`, which
fits no NNLS lengths. `tbs_diffusion_adaptive_nnls`, an NNLS method, changes on
no case. On these cases the solver switch changes no result by itself.

## New cancer surrogate cases

`tbs_family_diffusion` returns 6, 5, 6 and 6 clusters (ARI 0.998, 0.593, 0.998,
0.697) on the hub/disconnected, hub/connected, no-hub/disconnected and
no-hub/connected cases; canonical `tbs` returns one cluster on all four. Both
Hamming diffusion methods skip the continuous inputs. K-means and spectral
clustering receive the true K.

## Files

- `full-benchmark-lbfgsb-summary-20261106.csv`: per-method matched summary.
- `full-benchmark-lbfgsb-changed-20261106.csv`: the eight changed rows.
