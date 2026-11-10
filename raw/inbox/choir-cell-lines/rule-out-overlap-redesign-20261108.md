# Redesigning the zero-overlap rule of the split-count rule-out

Analysis on 2026-11-08. Local folder: `results/analyses/rule_out_overlap_20261108` (git-ignored); protocol (with
one amendment), deviations, code and outputs in `rule-out-overlap-redesign-scripts-20261108.py.txt`.

## Redesign

- New overlap rule: the Fisher discriminant of the split in the invariant directions (pooled covariance,
  observed proportions) may misplace at most the Bayes overlap of the weakest separation the count of all splits
  can accept, pi(1 - pi) delta^2 = 2 ln 2 (balanced: delta 2.36, 11.9%). Self-check: 2.2% misplaced at delta = 4
  (theory 2.3%).
- Amendment after E1: the old 5% equidistant-cell rule in the full whitened space also tightened with n and p,
  so a tie rule on the discriminant coordinate was added, allowing the share expected at the same weakest
  separation (balanced 4.7%; simulated 4.72%).

## E1 grid (8,640 runs, same seeds as `rule-out-gradient-20261108.md`)

With both new rules the kept rate rises with n at fixed separation (p = 10, delta = 3: 0.25 at 20 cells, 0.85 at
640; p = 30, delta = 4: 0.78 to 1.00; delta = 6: 0.98-1.00 everywhere above the n <= p + 1 edge); random splits of
one population kept 0 of 1,080. With only the new overlap rule the old tie rule bound at delta 3-4 and declined
with n (p = 30, delta = 4: 0.98 at 40 cells, 0.08 at 640).

## E2 end to end (935 datasets, both trees; today from the exclusion run)

| Markov/Ward | today | old rule (zero overlap) | new overlap | new overlap and ties |
| --- | --- | --- | --- | --- |
| Large no-group splits passing (sd 0.25 / 0.42) | 71.8% / 73.1% | 0.0% / 0.9% | 0.4% / 1.6% | 0.6% / 1.9% |
| Small no-group splits passing (sd 0.25 / 0.42) | 1.0% / 3.1% | 0.7% / 1.6% | 0.7% / 2.0% | 0.7% / 1.6% |
| False partitions of 300 (sd 0.25 / 0.42 / 1.0) | 298 / 296 / 300 | 2 / 25 / 15 | 13 / 42 / 36 | 17 / 46 / 141 |
| ARI balanced / imbalanced / singleton / density | 0.32 / 0.08 / 0.09 / 0.03 | 1.00 / 1.00 / 0.88 / 0.87 | 0.99 / 0.82 / 0.73 / 0.53 | 0.96 / 0.64 / 0.76 / 0.26 |
| Median clusters, density variation (true 4) | 81 | 5 | 15 | 17 |

Default tree: false partitions at sd 1.0 0 (old), 9, 47; real-group ARI unchanged in every arm. All three arms
meet the frozen rule (per-split passes <= 2% at sd 0.25 and 0.42, default-tree ARI within 0.02); the rule does not
check false partitions per dataset or Markov/Ward ARI, where the new rules are worse.

## Why

Tested splits with an invariant direction, Fisher misplaced fraction (median) and share at or below the allowed
overlap (no-group data, Markov/Ward): 32-79 cells 0.147 / 40%, 80-159 0.200 / 15%, 160-399 0.187 / 29%, 400+ 0.130
/ 41%; fully separable 0-12%. Splits the tree picked from one population are much cleaner than random cuts and
overlap about as much as real groups at the weakest separation the count accepts, so no overlap fraction tied
to the count can tell them apart. The zero-overlap rule worked because it is far stricter than tree-chosen noise
gets; the price is that its required separation grows with the number of cells.

## Status

Production keeps the zero-overlap rule; the protocol's replacement clause was not applied (deviations.md). Not
run: E3 (real CHOIR) and E4 (continuous benchmark) with the new rules.

## Correction (2026-11-08, later)

`allowed_overlap` had its two Bayes-error terms swapped: right for balanced splits (11.9%), too large otherwise
(7.3% instead of 1.3% at a 10/90 split). Rerun with the formula fixed in
`single-population-split-bound-20261108.md`: Markov/Ward false partitions 13 / 17 / 28 instead of 13 / 42 / 36;
the conclusion (an overlap fraction does not separate tree-picked noise splits from weak groups) stands.
