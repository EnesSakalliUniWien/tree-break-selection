# Descendant chain-sampler validation protocol

Fixed on 2026-11-06 before any result of the descendant test was computed.

## Test

`selected_diffusion_family_chain_tests` in
`benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py`
tests a pre-set family: the root and each internal child of the root. Each
node's mean contrast between its two children is tested by the random-position
reversible chi chain on that contrast's radius (`descendant_path`: only the
node's rows move), conditional on the node's ancestor-path event: every split
from the root down to the node must recur. At the root this is the validated
root-split event. Model: independent Gaussian rows with known covariance
`sigma**2 * I`, here `sigma = 1`. Steps 999, proposal scale 1, family level
`alpha / 3` with `alpha = 0.05` regardless of how many children are internal.
A node is rejected only when its chain completes, moves and has
`p <= alpha / 3`. Seeds are fixed per dataset and family position.

## Panels

Data seeds are `20261106 * 1000 + replicate` per panel offset; chain seeds are
`2026110600 + 10 * replicate + position` (root 0, children 1 and 2).

| Panel | n | p | Truth | Datasets |
| --- | ---: | ---: | --- | ---: |
| A global null | 30 | 5 | one population | 200 |
| B root signal, null children | 30 | 5 | 15 vs 15, shift 4 in feature 1 | 100 |
| C hierarchical signal | 40 | 5 | 2 x 2 groups of 10: shift 6 in feature 1, then 3 in feature 2 within each side | 50 |
| D larger nulls and signals | 120, 200 | 5 | A and C designs scaled (groups of n/2 and n/4) | 10 each |

A tested node is null when all its rows come from one true population; its
rejection is then a false rejection. A family has a false rejection when any
null node in it is rejected.

## Endpoints

- Per-position false-rejection counts and the family-wise false-rejection count
  with one-sided 95% exact binomial upper bounds; the family-wise target is 0.05.
- Rejection rates for non-null nodes (root in B, C, D; children in C, D).
- Unresolved, stalled and non-internal positions are reported and counted as
  non-rejections; no dataset is rerun, reseeded or dropped.

## Limits fixed in advance

These are known-noise Gaussian panels for the validated replay geometry only.
They do not validate counts, fitted PCA, unknown covariance, nodes below the
root's children, or the production traversal with pass-through.

## Implementation detail (added before any panel was run)

Panel data-seed offsets: A 0, B 10,000, C 20,000, D-null-120 30,000,
D-signal-120 40,000, D-null-200 50,000, D-signal-200 60,000; the data seed is
`20261106 * 1000 + offset + replicate`. Signal shifts are split symmetrically
(for example -2 and +2 for a shift of 4). Runs use one BLAS thread per process
and seven worker processes.
