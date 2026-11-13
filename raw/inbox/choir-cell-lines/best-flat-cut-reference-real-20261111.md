# Gaussian best-flat-cut reference on the CHOIR cell lines

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/best_flat_cut_reference_real_20261111`
(git-ignored); protocol (written before the run), runner, summary script and summary output are in
`best-flat-cut-reference-real-scripts-20261111.py.txt`. The opt-in `sibling_best_flat_cut_reference` was
uncommitted working-tree code on top of `bc044e07`. No production change; nothing was promoted.

## Setup

Setup of the 9 November ceiling-3 run: 48,879 cells, 30 PCs, the saved diffusion coordinates of
`choir_cell_lines_tbs_lbfgsb_20261106`, registry parameters, edge alpha 1e-3, sibling alpha 1e-2, both
registered methods. Arms through the production opt-in paths: `rule_out`
(`sibling_split_count_rule_out=True`, replay check) and `reference` (`sibling_best_flat_cut_reference=True`,
as implemented: n from 15, k up to 10, n rounded down to the table grid with largest value 800, lowest
possible p-value 0.0005). All four fits finished with status ok (110 to 183 s each).

## Replay check

`rule_out` gives the earlier result exactly: 69 clusters (ARI 0.0182) on the default tree and 664 clusters
(ARI 0.7652) on Markov/Ward. The same 90 and 1,838 nodes were tested as in the ceiling-3 run, with equal k in
all of them and equal `rho` in all 22 and 1,285 tested splits.

## Results (190 cell lines)

| Method, arm | Clusters | ARI | NMI | Lines 90/90 | ARI vs CHOIR_1 | Clusters >= 90% pure | Clusters < 20 cells | Largest cluster | Lines split only / merged only / both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Markov/Ward, today (8 Nov) | 954 | 0.549 | 0.863 | 41 | 0.445 | 761 | 374 | 1,451 | not recorded |
| Markov/Ward, rule_out | 664 | 0.765 | 0.901 | 80 | 0.514 | 494 | 253 | 1,990 | 72 / 11 / 27 |
| Markov/Ward, ceiling3 (9 Nov) | 217 | 0.833 | 0.940 | 119 | 0.538 | 164 | 25 | 1,990 | 16 / 33 / 22 |
| Markov/Ward, reference | 138 | 0.297 | 0.851 | 76 | 0.163 | 98 | 10 | 7,678 | 9 / 91 / 14 |
| Default, today (8 Nov) | 75 | 0.028 | 0.548 | 27 | 0.014 | 48 | 20 | 28,611 | not recorded |
| Default, rule_out | 69 | 0.018 | 0.497 | 29 | 0.009 | 49 | 12 | 32,680 | 4 / 145 / 12 |
| Default, ceiling3 (9 Nov) | 64 | 0.028 | 0.549 | 29 | 0.014 | 46 | 9 | 28,611 | 4 / 145 / 12 |
| Default, reference | 53 | 0.018 | 0.488 | 25 | 0.009 | 37 | 7 | 32,890 | 3 / 150 / 12 |

## Tested splits of the reference arm

| Tree | Subtree cells | Reason | Splits | p <= 0.01 | Accepted by BH | Median k (largest) | Median `rho` |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: |
| Markov/Ward | at most 14 | no_invariant_directions | 78 | 0 | 0 | 0 | |
| Markov/Ward | 15 to 31 | no_invariant_directions | 285 | 0 | 0 | 0 | |
| Markov/Ward | 32 to 159 | tested | 867 | 15 | 5 | 6 (10) | 0.97 |
| Markov/Ward | 160 to 999 | tested | 418 | 42 | 26 | 9 (10) | 0.92 |
| Markov/Ward | 160 to 999 | too_many_invariant_directions | 87 | 0 | 0 | 11 (13) | |
| Markov/Ward | 1,000 and more | too_many_invariant_directions | 103 | 0 | 0 | 14 (15) | |
| Default | 15 to 31 | no_invariant_directions | 8 | 0 | 0 | 0 | |
| Default | 32 to 159 | tested | 12 | 3 | 3 | 6 (8) | 2.49 |
| Default | 160 to 999 | tested | 10 | 3 | 3 | 10 (10) | 0.49 |
| Default | 160 to 999 | too_many_invariant_directions | 11 | 0 | 0 | 12 (13) | |
| Default | 1,000 and more | too_many_invariant_directions | 49 | 0 | 0 | 14 (15) | |

Totals: Markov/Ward 1,838 splits, 57 with p <= 0.01, 31 accepted; default tree 90, 6, 6. p-values of the
1,285 tested Markov/Ward splits: 48 at 0.0005, 9 above it up to 0.01, 5 up to 0.05, 1,223 at 1.0.

## Why the reference refuses what ceiling 3 accepted (Markov/Ward; added after the endpoints)

Of the 120 splits accepted under ceiling 3, the reference accepts 25 and refuses 95; it accepts 6 that
ceiling 3 refused.

- 77 of the 95 are refused as `too_many_invariant_directions` (k of 11 to 15, median subtree 1,022 cells).
  Their `rho` in the ceiling-3 run has median 7.75 and minimum 3.7; the table's 99% point at n = 800 is 2.48
  for k = 10. All 190 splits with more than 10 invariant directions are refused, including every split with
  1,000 cells or more.
- 13 of the 95 were tested and have the lowest possible p-value 0.0005 (`rho` 6.9 to 55.6) and are still not
  accepted by the BH step. Overall 17 of the 48 splits at 0.0005 are not accepted; every accepted split has
  p = 0.0005. Whether the cause is the floor 1 / (B + 1) alone or also refused ancestors was not separated.
- 5 of the 95 were tested with p-values of 0.001 to 0.021.

## Reading

- As implemented, the reference is worse than both earlier rules on the real data: Markov/Ward ARI 0.297
  against 0.833 (ceiling 3) and 0.765 (rule-out), lines recovered 76 against 119 and 80, lines merged only 91
  against 33 and 11, largest cluster 7,678 cells against 1,990.
- The cause is the implemented limits, not the statistic: the limit k up to 10 removes every split at the top
  of the tree, and the lowest possible p-value 0.0005 is too large for the BH step over 1,838 splits. Both
  were stated in the protocol (the first) or in the implementation's limits (the second) before the run.
- The limit k up to 10 comes from the search: for k = 15 and 20 the measured 99% points are lower bounds only
  (`gaussian-best-flat-cut-reference-20261111.md`), so extending the table with the same search would give a
  reference that is too low.
- The default tree does not change materially; its giant cluster comes from the topology, as before.
- Expectation check: the refusal of all splits with more than 10 invariant directions held (190 of 190). The
  expectation that a similar share of the splits with k of 1 to 10 would be kept held before BH (57 with
  p <= 0.01; ceiling 3 kept 44 of the same 1285 nodes) but not after it (31 accepted against 43).

## Limits

- One run. These data have no known no-group parts, so error control cannot be read from them.
- Not run: a larger reference sample (lower floor), a table beyond k = 10, or a restriction to the 10 leading
  invariant directions.
- ARI against CHOIR_1 is agreement with a published clustering, not accuracy.
