---
title: Plan for a Non Heuristical Statistic, Separated Spaces, Topology Alternatives and the Neighborhood
type: question
status: draft
updated: 2026-11-13
sources:
  - raw/inbox/choir-cell-lines/plan-non-heuristical-statistic-20261110.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-k15-20261111.md
  - raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-20261111.md
  - raw/inbox/choir-cell-lines/spr-alternatives-best-flat-cut-20261111.md
  - raw/inbox/choir-cell-lines/spr-subspaces-equidistances-20261111.md
  - raw/inbox/choir-cell-lines/spr-eigenvalue-blocks-20261111.md
  - raw/inbox/choir-cell-lines/spr-refit-20261111.md
  - raw/inbox/choir-cell-lines/cancer-pc-edge-test-20261112.md
  - raw/inbox/choir-cell-lines/cancer-pc-edge-k-20261112.md
  - raw/inbox/choir-cell-lines/selected-tree-null-over-splitting-20261112.md
  - raw/inbox/choir-cell-lines/peel-pass-20261112.md
  - raw/inbox/choir-cell-lines/peel-pass-row-orders-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-default-test-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-opt-in-20261113.md
  - raw/inbox/choir-cell-lines/peel-pass-production-20261113.md
  - tree_break_selection/tree/topology.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/best_flat_cut_reference.py
  - tree_break_selection/tree/spr.py
tags:
  - plan
  - sibling-gate
  - edge-gate
  - topology
  - calibration
---

# Plan for a Non Heuristical Statistic, Separated Spaces, Topology Alternatives and the Neighborhood

## Question

The plan of 10 November 2026 asks for one opt-in profile in which the sibling
statistic is derived from a stated population model, the separated spaces of
each subtree identify topology alternatives, and all contrasts in a
neighborhood use one common Sigma. Which of its 13 steps are done, which
changed, and which are open, with the evidence up to 12 November 2026?

## Current State

Of the 13 steps, 2 are done in a changed form, 2 exist as measurements only,
2 are partly done, 1 has a prototype outside production and 6 are not started. Defaults and the legacy configuration
are unchanged, as the plan requires.

| Step | Plan | State on 12 November |
| --- | --- | --- |
| A1 | Correct the single-population bound; proof or counterexample | Changed. No proved closed-form bound; replaced by the simulated Gaussian null of the best flat cut's `rho` per (n, k) |
| A2 | Derive the noise band that defines k | Prototype outside production ([[k-margin-prototype]]): simulated margin for an independent sample, also for more features than cells; it recovers three benchmark cases and fails on the CHOIR cells (ARI 0.8219 to 0.1285), because one cell line is far from spherical. Production k still comes from `invariant_directions`; a common whitening, with or without labels, does not help because each line has its own covariance ([[line-non-sphericity]]); a draft with k = 1 from a basis half and the test on the other half ([[k-split-direction-draft]]): with the direction from the tree's split it gives continuous benchmark mean ARI 0.8302 (0.6950) and on the CHOIR cells an ARI within the spread of the row order around the current rule's 0.8219, nulls at 1 cluster; the level of that variant is measured, not derived: p <= 0.01 at 0.35% of 31,482 tested nodes of one Gaussian population, more than 1 cluster in 12 of 1,960 fits. An independent check found an uncaught error for a constant feature and a dependence on the order of the rows (CHOIR ARI 0.8043 to 0.8397 over five orders); both are fixed in a third version (20 random divisions, twice the median p-value, exact null of 99,999 draws for one direction): CHOIR ARI 0.8349 to 0.8411 over five orders, 0 of 800 one-population fits split; it is outside production |
| A3 | `rho` against the bound with a count of splits | Changed. Opt-in `sibling_best_flat_cut_reference`: `rho` against the table, no count, k from 1 to 15; it replaces the sibling p-value and, since 13 November, lets a small peel (smaller child below 15 leaves) pass the edge gate like a singleton peel |
| B1 | Move `spectral_subspaces`, `spectral_geometry`, `compare_subspaces` into `tree_break_selection/space_separation` | Not started; they are in `benchmarks/diagnostics` |
| B2 | Margins per space, equidistant and contradicting cells per split | Measured in an analysis (tree, pcs, invariant, fuzzy), not in the package |
| B3 | Block equality test; enumerate separations the tree cannot hold where m >= 3 | Measured in an analysis; no numerical tie on the real data, so nothing to enumerate; sprs took the place of the enumeration; under the third version of the draft k the sprs do not improve the CHOIR fit either (ARI 0.7691 for all sprs of 11 November, 0.8266 and 0.8257 for sprs at selected refused nodes, against 0.8359 unchanged, [[spr-refit]]) |
| B4 | Return primary and alternative partitions through the gate metadata | Not started |
| C1 | One common Sigma for a node, its parent, children and siblings; drop the edge half-split, the global branch mean and the covariance fallback at 8 leaves | Not started as planned; one defect of the edge gate under the opt-in is removed: a peel of fewer than 15 cells no longer closes the descent ([[k-split-direction-draft]]) |
| C2 | Move `compare_surviving_paths` and the neighborhood records into reporting | Not started |
| D1 | Fixtures | Partly: tests of the reference; benchmark null case `gauss_null_256n_30f_continuous` |
| D2 | Frozen 120,000-dataset and 3,200-dataset panels | Not started |
| D3 | Acceptance: paired AMI improvement, runtime at most 5x | Not evaluated |
| D4 | Record each task as supported, rejected or inconclusive | Done for the tasks that were run |

Three changes of direction:

1. Statistic. The reference is a simulated null, not a proved bound, and the
   population class is one Gaussian population, not "unimodal, ceiling 3". The
   count of splits is gone because the best flat cut is the best of all
   splits. On the CHOIR cells the Markov/Ward fit reaches ARI 0.8219 with the
   reference at k up to 15; ceiling 3 gave 0.8329.
2. Topology alternatives. Sprs inside a node are bounded by its best flat cut:
   0 of 1,007 refused splits pass after sprs. With the reference at k up to 15
   the unchanged tree (ARI 0.8219) is better than the tree after sprs
   (0.7376). The plan's default "a stronger alternative is reported beside the
   tree's split, not swapped in" stands. Set based node identity (`Partition`,
   `PartitionSet`) and `spr.py` are in the package and were not in the plan.
3. Neighborhood. The edge test and the sibling test hold only for a tree that
   does not depend on the tested cells. One Gaussian population in 30 features
   is split by today's complete method in 46 of 60 runs and in 0 of 80 with
   the opt-in. The three things C1 drops are measured as not protecting: the
   half-split leaves 15.2% to 26.1% rejections at 0.05 for a tree built from
   the data against 1.2% to 2.8% for an independent tree; the branch
   multiplier is 1.00 to 1.03 at the median; the fallback at 8 rows together
   with the eigenvalue tie rule sets the edge test's k to 30, 30 minus test
   rows, or 2.

What the evidence implies for the open steps (reading, not a decision):

- A2 is the bottleneck of part A. Subtrees with no invariant directions are
  refused and lose real groups (continuous suite mean ARI 0.8802 and 0.8525
  to 0.6950), the earlier limit of 10 directions caused the poor real-data
  result, and k chosen from the same data shows a small excess at 100 cells
  (152 of 12,000).
- B3 and B4 have little reason left in their planned form.
- C1 is now the main open defect, because the edge gate stays as it is when
  the opt-in is on. On the default tree of the CHOIR cells the sibling
  statistic changes nothing (64 and 63 clusters, ARI 0.0275 under the
  current rule and the third version of the draft k): the chain topology
  and the gate before the sibling test are the limit there
  ([[k-split-direction-draft]]). The gate closes at a peel of 4 cells above
  28,606 (edge raw p 0.044) and Tree-BH then leaves 101 nodes of at least
  1,000 cells untested. A trial at run time that lets peels of fewer than
  15 cells pass the gate like a peel of one cell raises the default tree
  from ARI 0.0275 to 0.6172 (402 clusters, largest 1,990 cells) and the
  Markov/Ward tree from 0.8349 to 0.8411 to 0.8642 to 0.8706 over five row
  orders (`raw/inbox/choir-cell-lines/peel-pass-row-orders-20261113.md`), with the benchmark unchanged in
  28 of 28 rows and one Gaussian population at 1 cluster in 880 of 880 fits
  under the third version (`raw/inbox/choir-cell-lines/peel-pass-20261112.md`).
  The edge gate itself is unchanged by the trial. With today's default
  sibling test the same peel pass splits one Gaussian population in 30
  features in 825 of 860 fits against 703 without it
  (`raw/inbox/choir-cell-lines/peel-pass-default-test-20261113.md`), so it
  cannot enter the gate without a sibling rule that holds its level. The
  opt-in's current rule does: with and without the peel pass 1 of 880 fits
  is split (`raw/inbox/choir-cell-lines/peel-pass-opt-in-20261113.md`). The
  peel pass is in production under the opt-in since 13 November
  (`raw/inbox/choir-cell-lines/peel-pass-production-20261113.md`): default
  tree ARI 0.6210, Markov/Ward tree 0.8214 against 0.8219 with 110 against
  116 lines recovered.
- D2 and D3 remain the release condition; the opt-in is not a default.

The plan's third default, continuous features only, is unchanged. Its closing
question, whether the count still covers the alternatives once k is chosen from
the same data, no longer applies to a count; the same-data choice of k is open
under A2.

## Evidence

- `raw/inbox/choir-cell-lines/plan-non-heuristical-statistic-20261110.md`
  holds the plan text verbatim.
- The reference and its extension:
  `tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/best_flat_cut_reference.py`,
  `raw/inbox/choir-cell-lines/best-flat-cut-reference-k15-20261111.md` and
  `raw/inbox/choir-cell-lines/best-flat-cut-reference-benchmark-20261111.md`.
- Sprs, margins per space and eigenvalue blocks:
  `raw/inbox/choir-cell-lines/spr-alternatives-best-flat-cut-20261111.md`,
  `raw/inbox/choir-cell-lines/spr-subspaces-equidistances-20261111.md`,
  `raw/inbox/choir-cell-lines/spr-eigenvalue-blocks-20261111.md`,
  `raw/inbox/choir-cell-lines/spr-refit-20261111.md` and
  `tree_break_selection/tree/spr.py`.
- Edge test, its k and the root cause:
  `raw/inbox/choir-cell-lines/cancer-pc-edge-test-20261112.md`,
  `raw/inbox/choir-cell-lines/cancer-pc-edge-k-20261112.md` and
  `raw/inbox/choir-cell-lines/selected-tree-null-over-splitting-20261112.md`.

## Links

- [[k-split-direction-draft]]
- [[line-non-sphericity]]
- [[k-margin-prototype]]
- [[marchenko-pastur-multiplicity-literature]]
- [[gaussian-best-flat-cut-reference]]
- [[best-flat-cut-reference-k15]]
- [[best-flat-cut-reference-benchmark]]
- [[single-population-split-bound]]
- [[flat-cut-sup-count]]
- [[ceiling-rule-real-and-benchmark]]
- [[spr-alternatives-best-flat-cut]]
- [[spr-subspaces-equidistances]]
- [[spr-eigenvalue-blocks]]
- [[spr-refit]]
- [[cancer-pc-edge-test]]
- [[selected-tree-null-over-splitting]]
- [[subtree-descriptions-and-statistical-tests]]
- [[choir-cell-lines-recovery-plan]]
