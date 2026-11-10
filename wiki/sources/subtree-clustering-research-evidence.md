---
title: Subtree Clustering Research Implementation Evidence
type: source
status: reviewed
updated: 2026-11-09
sources:
  - raw/inbox/choir-cell-lines/subtree-clustering-research-20261109.md
  - raw/inbox/choir-cell-lines/subtree-clustering-research-20261109-manifest.json
  - raw/inbox/choir-cell-lines/subtree-clustering-research-20261109-sources.py.txt
  - reports/subtree_clustering_research_20261109/protocol.md
  - reports/subtree_clustering_research_20261109/promotion_gate.json
  - reports/subtree_clustering_research_20261109/input_integrity.json
  - reports/subtree_clustering_research_20261109/verification.json
  - reports/subtree_clustering_research_20261109/calibration/wolfram_queries.md
tags:
  - subtree
  - calibration
  - benchmarks
  - implementation
---

# Subtree Clustering Research Implementation Evidence

## Summary

The approved research plan was implemented with three Sol agents. Shared metrics
and statistical profiles, native score attribution, bounded all-pair topology
search, a common-event Ward family reference and evidence validators are
available. Development results reject promotion of the topology-search candidate.
Production algorithm sources and default presets were preserved.

## Key Points

- Completed 2,400 native factorial fits and 1,600 recovery fits on 850 distinct
  matrices. The 100 signal datasets in the factorial are reused recovery inputs.
  All fits returned software status `ok`; this does not establish calibration.
- Both symmetric/average and Markov/Ward NNLS with .90 energy truncation reject
  the selected Gaussian root in 47/50 datasets at sibling alpha .01. Independent
  dispatcher replay and branch/rank/correction reconstruction verify the effect.
- Markov/Ward minus symmetric/average mean AMI is -0.05550, with paired interval
  [-0.06360,-0.04740], across the 16-cell development grid. Nine cells have
  Holm-adjusted evidence of AMI loss. This compares complete configurations.
- Matched group F1 exposes fragmentation hidden by majority-remapped macro F1.
  Statistical profiles retain group size, covariance spectrum, separation,
  overlap and explicitly supplied subtree splits. Singular covariance leaves
  Mahalanobis distance unavailable; reference labels never enter fitting.
- Exact compiled NNLS agrees with dense fits for all 945 rooted six-leaf trees
  in both target-distance arms. Larger candidates use whole eigenspaces and
  leading spans, with incomplete coverage explicitly reported up to 256 leaves.
- The six-cell isolated timing panel includes 36 complete fits. Candidate time
  is 14.10 times baseline at n=128,p=72 and 24.88 times at n=240,p=72, exceeding
  the agreed fivefold cap. The new search lacks a production selection reference.
- The 300-record Gaussian/attribution study retains numerical corrections.
  All 50 seven-hypothesis mixed-root families are unresolved at 99 steps; at
  999 steps there is one false family rejection and no true family rejections.
  This restricted family construction has not demonstrated useful mixed-root power.
- A further 10,000 perturbations compare five laws. These frequencies do not
  count all possible topologies or constitute topology confidence probabilities.
  CHOIR and earlier synthetic profiles are descriptive reuse, not fresh validation.

## Evidence

The primary report, protocol, source archive and manifest retain source/data
hashes, all attempts, full replicate-zero native traces, exact Wolfram inputs
and outputs, corrections, figures and reproducible summary commands. Regenerating
all development inputs matches every stored data/truth hash. The 110 recorded
production source hashes are unchanged from the pre-edit snapshot.

The broad check passed 214 tests with two cases deselected. The final validator
rerun passed 18 tests; their union is 216 distinct passing cases. Final style,
wiki and input-integrity checks are recorded in the verification artifact.
The 120,000 frozen calibration datasets and 3,200 frozen recovery datasets were
not run because prerequisites failed. No production-improvement claim follows.

## Links

- [[subtree-clustering-production-gates]]
- [[subtree-distance-topology-families]]
- [[subtree-topology-reference-evidence]]
- [[root-validation-generalisation-evidence]]
- [[diffusion-topology-branch-rank-isolation]]
