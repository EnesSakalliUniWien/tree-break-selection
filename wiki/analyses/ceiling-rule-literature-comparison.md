---
title: The Ceiling Rule Against the Literature
type: analysis
status: draft
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/ceiling-rule-literature-20261109.md
  - raw/inbox/choir-cell-lines/single-population-split-bound-20261108.md
  - raw/inbox/choir-cell-lines/ceiling-rule-real-and-benchmark-20261109.md
  - raw/inbox/choir-cell-lines/three-tests-audit-20261109.md
tags:
  - calibration
  - literature
  - sibling-gate
---

# The Ceiling Rule Against the Literature

## Summary

The ceiling rule shares SigClust and SHC's one-population testing aim; their
Gaussian null and the proposed unimodal null are different. In one dimension the 2-means
cluster index is 1/(1 + rho). Where SHC simulates the null and SigClust-type tests
are anti-conservative for heavy tails at large n, the ceiling rule uses a closed-form
ceiling plus a count of flat cuts. The Gaussian ceiling repeats SHC's weakness on
skewed groups. Ceiling 3 is the value of the uniform, the least favourable unimodal
shape of the dip test. The closed form g(pi) and the finite-n count were not found in
the sources checked.

## Details

### Formula by formula

| Ours | Literature | Agreement |
| --- | --- | --- |
| rho = pi(1-pi) D^2, T = n rho | Gao et al. 2024 eq. (4): Wald mean difference, scaled chi_q, known sigma^2 I | Same statistic family; ours uses the pooled within-sides covariance |
| Null: one population cut in two | SHC and SigClust: the node is one Gaussian | Same null; ours also tries any unimodal shape (ceiling 3) |
| g(pi), 1.75 at 50/50 | SHC 2-means CI; CI = 1/(1+rho) in 1D, so 1.75 is CI 0.363 = 1 - 2/pi | Same quantity; SHC simulates it ("no closed form"); 2/pi derived here, not found in a source |
| Ceiling 3 | Dip test: the uniform is the least favourable unimodal law | Same shape; the dip is a different statistic, so the 8 of long flat tails is no contradiction |
| delta_0 = 2.65 (Gaussian), 3.46 (ceiling 3), 50/50 | Equal mixture is bimodal only above 2 sigma (Schilling et al. 2002) | Rule asks for more than bimodality; at delta 2, rho = 1.0 < 1.75 |
| M_k(n) = sum_{i<=k} binom(n-1, i) | Cover: C(N, d) = 2 sum_{k<d} binom(N-1, k); affine in k dims = C(n, k+1)/2 | Same count (checked against Wikipedia, paper not read) |
| Error control: count of flat cuts, then BH | Gao: exact p-value conditioned on the clustering; SHC: Monte Carlo + Meinshausen FWER | No general guarantee for the ceiling rule: the later [[three-tests-construction-audit]] disproves its broad null and finds boundary miscalibration |
| Best split of one normal sample, between/within | Engelman and Hartigan 1969 percentage points | Same statistic in 1D; their tables unverified |

### Results against the literature

- The later [[three-tests-construction-audit]] checks Cover's original count,
  distinguishes the mode-cut ceiling from arbitrary unimodal cuts, and confirms
  that the dip statistic's least-favourable shape cannot establish a reference
  law for this different statistic.
- Gaussian ceiling over-splits skewed groups (density-variation ARI 0.57): matches SHC
  being anti-conservative for heavy tails, more so for large samples.
- Ceiling 3 merges more cell lines (11 to 33 merged only): consistent with SigClust's
  low power for unbalanced splits (Keefe and Marron 2023), since a constant ceiling
  asks small groups for a larger gap. The merge check found this for 12 of 39 merged
  line pairs; the other 27 pass ceiling 3 on their true labels, so the tree is the cause
  there ([[ceiling-rule-real-and-benchmark]]).
- Real data, ARI against the published CHOIR clustering 0.538 (ceiling 3): agreement,
  not accuracy. CHOIR's own cell-line results were not looked up.

## Evidence

Lookups and what was read of each source are in
`raw/inbox/choir-cell-lines/ceiling-rule-literature-20261109.md`. Derivation and
simulations: [[single-population-split-bound]]; real data and benchmark:
[[ceiling-rule-real-and-benchmark]].

## Links

- [[single-population-split-bound]]
- [[ceiling-rule-real-and-benchmark]]
- [[split-count-exclusion]]
- [[three-tests-construction-audit]]

## Open Questions

- Do Engelman and Hartigan's 1969 percentage points match g(0.5) at large n?
- Does conditioning on the tree (Gao et al.) give more power than the flat-cut count?
- Why does the Markov/Ward tree not offer the line boundary for 27 merged pairs whose
  true split passes ceiling 3?
