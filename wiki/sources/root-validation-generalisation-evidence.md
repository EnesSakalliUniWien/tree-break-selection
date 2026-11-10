---
title: Root Validation and Clustering Generalisation Evidence
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/root-validation-followup-20261109.md
  - raw/inbox/choir-cell-lines/root-validation-followup-20261109-evidence.json
  - raw/inbox/choir-cell-lines/root-validation-followup-20261109-sources.py.txt
  - raw/inbox/choir-cell-lines/root-validation-calibration-20261109.md
  - raw/inbox/choir-cell-lines/root-validation-calibration-20261109-sources.py.txt
  - raw/inbox/choir-cell-lines/root-validation-metrics-20261109.md
  - raw/inbox/choir-cell-lines/root-validation-metrics-20261109-evidence.json
  - raw/inbox/choir-cell-lines/root-validation-metrics-20261109-scripts.py.txt
  - raw/inbox/choir-cell-lines/root-validation-mathematics-20261109.md
  - raw/inbox/choir-cell-lines/root-validation-wolfram-20261109.json
  - raw/inbox/choir-cell-lines/root-validation-claude-context-20261109.md
tags:
  - calibration
  - selective-inference
  - clustering-metrics
  - generalisation
---

# Root Validation and Clustering Generalisation Evidence

## Summary

Three Sol agents and a coordinated count panel extend the original root validation
framework. The known-noise chain is conservative on 1,200 small Gaussian nulls at
a fixed 119-step budget. B's CHOIR Markov/Ward recovery gains persist across
several metrics, but small groups, weaker separation and anisotropy expose losses.
Neither these gains nor zero failures on finite count panels establish general
calibration or final-partition error control.

## Key Points

- All 24 historical root chains reproduce exactly, and seven sampled descendant
  summaries match. Across 1,664 fresh families, five chains stall and one observed
  tree is unsupported; every attempt remains in the denominator.
- At n=9,p=3, root .05 rejects 38/1,200 and family .05 rejects 19/1,200; one-sided
  exact upper limits are .04130 and .02315. The n=12,p=8 and n=24,p=5 family
  upper limits are .06334 and .10677, so those panels do not substantiate .05.
- On 20 paired balanced signals, 999 versus 119 steps increases root .05
  rejection from 6 to 12, with median family runtime 1.42 to 11.65 seconds.
  Poor child splits remain a separate limitation.
- CHOIR Markov/Ward A/B/C/rule-out macro F1 is .617/.800/.746/.741; among the
  48 smallest lines it is .701/.715/.557/.804. Scores use one-to-one F1 matching,
  zero for unmatched lines, and equal line weights.
- All 448 fresh Gaussian-mixture fits and 224 jitter comparisons were rescored.
  In the imbalanced Ward panel B/C improve ARI .332 to .616 while reducing
  macro F1 .717 to .510 and losing the 12-cell component on every seed. At
  weaker separation their single-cluster result has truth ARI 0 and stability 1.
- All 800 fresh count fits succeed after refitted preprocessing; B/C give zero
  false partitions in each 30-dataset depth stratum and recover the four strong
  profiles. Zero of 30 still has a one-sided .09503 upper limit. All 200 matched
  dataset/tree geometries and 1,298 B/C split records agree.
- In the separate all-coordinate Gaussian n=32,k=20 cut diagnostic, B rejects
  818/1,000 searched thresholds at .01; C clips to 1 throughout the attainable
  statistic range. These isolate reference behaviour and are not full TBS fits.
- Four Wolfram checks and Gaussian numerical controls support the restricted
  exact root geometry in [[gaussian-root-calibration-and-cut-geometry]]. A
  selection-aware reference can retain the sibling score in principle. It was
  unimplemented in this source; the later restricted implementation and its
  limits are recorded in [[subtree-topology-reference-evidence]].

## Evidence

The combined raw report links all protocols, source archives, numerical summaries,
Wolfram inputs/results and selected Claude source context. Detailed chains,
assignments and PNG/PDF/SVG comparisons remain in the ignored local run directory
`results/analyses/root_validation_followup_20261109/`. There were 1,248 new full
fits; saved CHOIR partitions were rescored without repeating the full-cell fit.
No production code, defaults, dependencies, commit or push changed.

## Links

- [[root-conditional-sampler]]
- [[descendant-chain-family-test]]
- [[three-tests-construction-audit]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[flat-cut-sup-count]]
- [[rule-out-gradient]]
