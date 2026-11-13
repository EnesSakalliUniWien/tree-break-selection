---
title: Gaussian Best-Flat-Cut Reference on the Continuous and cancer_pc Benchmark Suites
type: source
status: draft
updated: 2026-11-12
sources:
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-20261111.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-scripts-20261111.py.txt
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-k15-20261112.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-k15-scripts-20261112.py.txt
tags:
  - benchmark
  - calibration
  - sibling-gate
---

# Gaussian Best-Flat-Cut Reference on the Continuous and cancer_pc Benchmark Suites

## Summary

On the 9 continuous and 4 `cancer_pc` benchmark cases the opt-in
`sibling_best_flat_cut_reference` removes over-splitting but loses real groups
wherever a subtree has no invariant directions, and in the connected
`cancer_pc` cases. Mean ARI on the continuous suite drops from 0.85-0.88
(today) to 0.695, the same value as the split-count rule-out; on `cancer_pc`
the Markov/Ward tree rises from 0.303 to 0.722 and the default tree drops from
0.822 to 0.722.

## Key Points

- 52 of 52 runs finished ok; the "today" cells reproduce the 8 and 10 November
  runs.
- Continuous suite: found clusters and ARI equal the split-count rule-out and
  the ceiling-3 rule in every case. All three use `invariant_directions`.
- Gains: dim_consolidated to the exact 4 clusters (ARI 1.000); Markov/Ward on
  `cancer_pc` from 10 to 16 clusters to 6, 4, 6, 4.
- Losses from `no_invariant_directions` (58 of 151 calls): 30 cells in 40
  features, 30 in 30, 40 in 20,000; dim_diffuse_6c_136f goes to 1 cluster
  (root tested, k = 3, `rho` = 1.61, p = 1.0).
- Losses in the connected `cancer_pc` cases: 4 clusters instead of 6 (ARI
  0.446; today's default tree 0.593 and 0.697). The refused splits are tested
  with k = 1 and `rho` of 0.15 to 1.03.
- With the reference both trees give the same found clusters and ARI in all 13
  cases. The null case stays at one cluster.
- Of 93 tested p-values, 56 were 0.0005, 36 were 1.0 and one was 0.056;
  `too_few_cells` and `too_many_invariant_directions` did not occur.
- True-group check on `cancer_pc`: the refused splits in the connected cases
  are impure (42 to 52 of the 64 cells of group 5 plus 4 to 9 cells of group
  0), and the exact true-group split would be refused too: `rho` 0.81 to 0.88
  in the invariant direction, 0.92 to 0.98 in all 30 features, against a 99%
  point of 2.151.
- Cause: a group of 64 next to 256 at offset 3 has population `rho` = 1.44,
  below the 1.752 of the median cut of one Gaussian population; the 99% point
  needs an offset of 3.67 at these sizes.
- A group whose complement holds another distant group is refused even at
  offset 20 (`rho` = 1.51); the pass-through traversal recovers it after the
  next split passes.
- Rerun on 12 November with the reference extended to 15 invariant directions
  ([[best-flat-cut-reference-k15]]): no change. Found clusters, ARI and NMI
  are equal in 52 of 52 rows; the 151 calls have at most 5 invariant
  directions, so the cells for k = 11 to 15 are never read; `rho` differs in
  80 calls by at most 1.0e-12 relative.
- One run per date, synthetic cases only, no real data; the true-group check
  compares one group against the rest only; no case has 11 to 15 invariant
  directions.

## Evidence

`raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-20261111.md` holds
the setup, the per-case table and the per-split reasons; protocol, runner,
result rows and all 151 recorded calls are in
`raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-scripts-20261111.py.txt`.
The rerun of 12 November is in
`raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-k15-20261112.md`
with protocol, scripts and comparison output in
`raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-k15-scripts-20261112.py.txt`.

## Links

- [[gaussian-best-flat-cut-reference]]
- [[cancer-pc-over-splits]]
- [[best-flat-cut-reference-k15]]
- [[gap-options-power]]
- [[rule-out-continuous-benchmark]]
- [[ceiling-rule-real-and-benchmark]]
- [[split-count-exclusion]]
