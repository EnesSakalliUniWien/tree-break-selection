<!-- Plan text of 2026-11-10 13:37 UTC from the working session on the CHOIR cell lines, saved verbatim on 2026-11-12.
     It was not in the repository before. Nothing below this comment was edited. -->

# Plan: non heuristical statistic, separated spaces, topology alternatives, neighborhood

**First step is A1: correct the single-population bound. Every later step uses it.**

## Objective

One opt-in profile in which:
- the sibling statistic is derived from a stated population model;
- the separated spaces of each subtree identify topology alternatives;
- all contrasts in a neighborhood use one common Sigma.

Defaults and the legacy configuration stay unchanged.

## A — Statistic

1. **A1. Correct the bound.** Fix the incorrect uniform-boundary variance in `g(pi)`. Turn the recorded 40/1,000 rejections at nominal 0.01 into a regression fixture. Result: a proof for a stated population class, or a counterexample.
2. **A2. Define k.** Replace the median-eigenvalue convention in the noise band by a derived band for selected subtrees. No-group Markov/Ward subtrees currently show 3–7 directions above the band.
3. **A3. Implement the statistic.** Add `rho = pi(1-pi) D^2` against the bound, with the count `sum_{i<=k} binom(n-1, i)`, next to `split_count_rule_out.py`. Reuse `invariant_directions` and `is_separable`. In this profile it replaces adaptive rank, energy truncation, the branch multiplier, the sibling dimension rule, the empirical inflation and the rule-out constants 0.05 and 0.1.

## B — Separated spaces and topology alternatives

1. **B1. Move the code.** Move `spectral_subspaces`, `spectral_geometry` and `compare_subspaces` from benchmarks into `tree_break_selection/space_separation`. Benchmarks import them from there.
2. **B2. Margins per space.** Compute the cell margin in the tree space, the test space, the invariant directions and the fuzzy directions. Record equidistant cells and cells on contradicting sides for every split.
3. **B3. Multiplicity.** Add the block equality test beside the numerical ties in `eigenvalue_groups.py`. Where m ≥ 3, enumerate the separations the tree cannot hold and test each with A3.
4. **B4. Output.** Return the primary partition plus the alternative partitions, each with its space, through the existing gate metadata.

## C — Neighborhood

1. **C1. One common Sigma.** Compute the contrasts of a node, its parent, children and siblings on one Sigma, so the telescoping and covariance identities hold. The profile drops the edge half-split, the global branch mean and the covariance fallback at 8 leaves.
2. **C2. Move the records.** Move `compare_surviving_paths` and the neighborhood records into the package's reporting. Pass-through is off in this profile until it has a derivation.

## D — Validation and release

1. Fixtures: 42 versus 4.5, the 16-leaf example, the four node cases, crossed partitions, and the 47/50 Gaussian-null rejections.
2. The frozen 120,000-dataset and 3,200-dataset panels, run only after a source-and-parameter lock.
3. Acceptance unchanged: positive lower confidence bound for paired AMI improvement, runtime at most 5× baseline, every failed or unresolved attempt retained.
4. Record each task as supported, rejected or inconclusive in the wiki.

## What changes against the old plan

- The conditional sampler becomes an independent check of A1, not the production reference. That removes the conflict with the 5× cap, because A3 is closed form.
- Alternatives are tested with the same statistic, not attached as an optional report.

## Defaults I chose; say so if any is wrong

1. **Population class: unimodal, ceiling 3.** It measured 0/0/0 false partitions and ARI 1.00 on all four group types. The Gaussian ceiling 1.75 over-splits skewed groups. Long flat tails reach 8 and are out of scope.
2. **A stronger alternative is reported beside the tree's split, not swapped in.** Crossing partitions cannot both be splits of one tree.
3. **Continuous features only.** The rule-out and the ceiling were never measured on Bernoulli or categorical data.

Order: A1 → A2 → A3; B1 can start alongside A1; B2–B4 and C1 need A3; C2 needs B1.

Not checked: whether the count in A3 still covers the alternatives once k is chosen from the same data; that is part of A1's proof. I gave no durations because A1 is a proof-or-counterexample task with no basis for an estimate.
