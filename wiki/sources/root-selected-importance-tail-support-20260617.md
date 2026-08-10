---
title: Root Selected Importance Tail Support 2026-06-17
type: source
status: draft
updated: 2026-08-10
sources:
  - benchmarks/diagnostics/calibration/root/tie_rank/root_tie_rank_null_proposal_frontier.py
  - benchmarks/diagnostics/calibration/root/selected/root_selected_spectral_tail_law_panel.py
tags:
  - source
  - diagnostics
  - root
  - spectral
  - importance-sampling
---

# Root Selected Importance Tail Support 2026-06-17

## Summary

The retained selected-root spectral-tail path supports importance-weighted
external-null rows. Proposal families may be labeled as external null support
only when they carry the likelihood-ratio metadata needed to compare the
matched iid Bernoulli null (P_0) with the proposal law (Q):

\[
\log w(X)=\log P_0(X)-\log Q(X).
\]

The tail panel uses those weights only inside the same reduced
\((T,A,E,H_u)\) stratum. It does not treat tilted proposal rows as ordinary
empirical null observations.

## Key Points

- Importance rows require explicit external-null calibration roles and finite
  likelihood-ratio metadata.
- The tail estimate uses an importance effective sample size and remains
  conservative when support is sparse.
- Selected spectral excess is the tail variable, not part of the conditioning
  key.
- The retired generated-neighborhood replay and topology-join path is no
  longer part of this contract.

## Evidence

- `root_tie_rank_null_proposal_frontier.py` computes proposal and target-null
  log probabilities for supported tilted Bernoulli families.
- `root_selected_spectral_tail_law_panel.py` separates ordinary calibration
  support from importance-weighted external support before computing a tail
  estimate.
- No non-empty retained result capture or focused current test supports the
  former numerical smoke-run claims, so this page remains draft.

## Links

- [[root-selected-spectral-tail-law-with-legacy-overlay-20260617]]
- [[root-tie-rank-null-proposal-frontier-20260616]]
- [[selected-neighborhood-measurability-law]]
