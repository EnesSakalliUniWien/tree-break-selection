---
title: Root Conditional Sampler Implementation
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/root-validation-calibration-20261109.md
  - raw/inbox/descendant-chain-validation-20261106.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py
  - benchmarks/diagnostics/calibration/README.md
  - tests/statistics/68_test_diffusion_root_conditional_chain.py
  - reports/root_conditional_sampler_20261104/README.md
  - reports/root_conditional_sampler_20261104/protocol.md
  - reports/root_conditional_sampler_20261104/manifest.json
  - reports/root_conditional_sampler_20261104/results.csv
  - reports/root_conditional_sampler_20261104/summary.csv
  - reports/root_conditional_sampler_20261104/verification.json
tags:
  - selective-inference
  - diffusion
  - implementation
---

# Root Conditional Sampler Implementation

## Summary

`selected_diffusion_root_chain_test` is a maintained opt-in diagnostic for one
native diffusion/average root. It implements the random-position reversible
chi chain from [[selective-rebuild-comparison]], preserving the observed root
children. It does not replace default production TBS gates.

## Key Points

- Known spherical Gaussian noise scale is mandatory. Independent rows and
  raw features satisfying that model remain caller assumptions; input checks
  do not establish the model for fitted-PC cancer data.
- Canonical labels, observed merge signature, selected root children, source
  hashes, chain states and observation position are retained. The root event
  permits other native merges to differ in simulation proposals.
- `decision` is `reject`, `do_not_reject` or `unresolved`. Geometry failures
  give no p-value. A stalled chain retains p=1 and an unresolved decision;
  insufficient rank resolution also prevents a decision. Failure is never
  silently treated as nonselection.
- The proposal scale is fixed at one. The budget and seed must be fixed before
  inspecting the result; there is no convergence stopping or repeated-seed
  search. Complete computation and movement do not certify exact-real event
  membership. The shared geometry replay is serial within a process.
- Existing rejection sampling remains available. The root API deliberately
  omits empirical inflation, selected projection rank and branch scaling;
  their production reference-law problems are not presumed solved.

## Evidence

The test-first implementation has 28 new tests covering the random anchor,
both chain directions, upper-tail ties, stalled paths, geometry failure,
argument checks, label ordering and a saved strong-signal regression. The
44-test focused replay/sampler batch passes, and the full ordered suite passes
1,031 tests with one unrelated scikit-learn MDS FutureWarning.

The fixed validation protocol completed 86 runs: 24 exact chain/rank parity
replays, 60 fresh cases and two registered raw Gaussian benchmarks. At 0.05
the fresh panel rejects 3/40 nulls, 4/10 moderate signals and 7/10 strong
signals; at 0.01 the counts are 0, 1 and 5. The registered clear and moderate
roots give p=0.005 and p=0.001 with their known noise scales. Verification
checks 258 native chain memberships and 149 source hashes. No run was
unresolved. These small panels do not establish cancer-data error control.

The same module now also tests the root's children conditional on their
ancestor paths; [[descendant-chain-family-test]] records that extension and
its validation.

The 9 November follow-up reproduces all 24 historical chains and runs 1,664
fresh families at a fixed 119-step budget. The primary n=9,p=3 null has
38/1,200 root rejections at .05 (one-sided exact upper limit .04130).
On 20 preselected balanced signals, 999 versus 119 steps raises root .05
rejection from 6 to 12, at higher cost. Larger/wider panels, unresolved cases,
model boundaries and descendant power are in
[[root-validation-generalisation-evidence]].

## Links

- [[root-validation-generalisation-evidence]]
- [[gaussian-root-calibration-and-cut-geometry]]
- [[descendant-chain-family-test]]
- [[selective-rebuild-comparison]]
- [[conditional-diffusion-monte-carlo-evidence]]
- [[statistical-construction-audit]]
- [[choir-cell-lines-recovery-plan]]
