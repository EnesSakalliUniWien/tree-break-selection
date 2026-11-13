# Rebuilding the best-flat-cut reference table with 19,999 draws per cell

Run on 2026-11-13, Europe/Vienna. Local folder: `results/analyses/table_rebuild_20261113` (git-ignored); protocol,
deviations, scripts and outputs are in `table-rebuild-scripts-20261113.py.txt`. Follow-up of
`peel-pass-production-20261113.md`. Environment here: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5; on
lisc: Python 3.12, NumPy 2.3.4, SciPy 1.16.3 from the locked requirements.

## Why

With the peel pass the opt-in's current rule refused 7 nodes on the Markov/Ward tree of the CHOIR cells, all at
the table's smallest p-value 1/2,000. Benjamini-Hochberg at 0.01 rejects that value only where 5% of a family is
at it.

## Step 1: best case before building (`floor_check.py`)

Every call at 1/2,000 was given 1/20,000 at run time. Markov/Ward tree: 243 clusters, ARI 0.8648, 124 lines
recovered (144 of 1,689 calls at the smallest value). Default tree unchanged: 432 clusters, ARI 0.6210 (188 of
676 calls at the smallest value). Predictions T1 and T2 held; the table was rebuilt.

## Step 2: the build

`scripts/build_best_flat_cut_reference.py --datasets 19999` on lisc, job 7980350, 64 hardware threads on 32
physical cores of an AMD EPYC 9554 node, 9,930 s (2 h 45 min), commit 10d971cd. The same build stopped here after
25 min at 350 of 32,000 blocks. The new file has 2,000 stored values per cell, 0.94 MB.

## Step 3: checks

| Check | Result |
| --- | --- |
| S1 99% points within 5% of the old table | not held: 5 of 160 cells differ by more than 5%, largest 10.8% at n = 15, k = 7, all at n = 15 or 20 |
| Same draws | the first 1,999 draws of four cells, recomputed here, equal the old table and are in the new one |
| Tests 82, 83, 66 with the new table | pass; full ordered suite passes |
| S2 Markov/Ward tree | 244 clusters, ARI 0.8648, NMI 0.9410, 124 lines recovered, 146 accepted (old table 223, 0.8214, 110, 130): held |
| S3 one Gaussian population, 880 fits | 2 fits with more than 1 cluster (old table 1): held (at most 2) |
| S4 benchmark, 28 rows | found clusters and ARI equal in 28 of 28: held |
| S5 default tree | 439 clusters, ARI 0.6255, 42 lines recovered (old table 432, 0.6210): held |

S1 failed because of my bound, not the table: the old 99% point is the 20th largest of 1,999 draws, whose tail
probability itself varies (Beta(20, 1980)), so it ranks about 116 to 284 among 19,999 draws; observed 124 to 312
(2.5% to 97.5%). At k = 1 the new 99% points match the exact computation of 12 November (5.402 against 5.343 at
n = 15, 2.529 against 2.531 at 100, 2.353 against 2.356 at 150, 2.091 against 2.091 at 400); the old table had
4.888 at n = 15.

## Reading

- The loss of the peel pass on the Markov/Ward tree came from the table's smallest p-value. With 1/20,000 the
  opt-in's current rule recovers 124 lines, more than without the peel pass (116), and the default tree keeps its
  gain.
- The new table is in production (`best_flat_cut_reference_table.npz`), with the builder's default at 19,999.

## Not checked

- One Gaussian population: 2 of 880 fits split; the second split fit was not examined.
- Row orders do not enter the current rule; one real dataset; no table with 99,999 draws.
