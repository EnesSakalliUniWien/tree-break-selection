# Adaptive alphas on the tree: protocol (fixed before any result)

Written 2026-11-06 on `dev` at `69360b99`, before any panel below was run.

## Question

Can the alphas of the conditional chain tests adapt to the tree that opens,
while still controlling error? Three procedures are compared with fixed
Bonferroni, and the assumptions they rest on (the "tree assumption") are tested
on null nodes.

## Model and family

- Independent Gaussian rows, 5 features, known covariance `I` (sigma 1).
- Tree: the native diffusion/average tree of `make_diffusion_root_path`.
- Family: every internal node of depth at most 2 (root, its internal children,
  their internal children), at most 7 positions, in breadth-first signature
  order. Position `j` uses chain seed `seed + j`.
- Node p-value: `_node_chain_test` (random-position reversible chi chain on
  the node's contrast radius, event = the node's split and every ancestor split
  recur), 999 steps. A stalled chain keeps p = 1; an unresolved node has no
  p-value and is never rejected.
- A node is null when all its rows come from one population; otherwise mixed.

## Panels (fresh seeds, base 2026110700)

| Panel | Rows | Design | Families |
| --- | ---: | --- | ---: |
| N30 | 30 | global null | 300 |
| R30 | 30 | root signal: feature 0 shifted -2 / +2 by halves; children null | 150 |
| H40 | 40 | 2 x 2: feature 0 +/-3 by halves, feature 1 +/-1.5 by quarters | 100 |
| H80 | 80 | 2 x 2 x 2: features 0, 1, 2 shifted +/-4, +/-3, +/-2 by halves, quarters, eighths | 100 |

Even replicates are the training half, odd replicates the test half.

## Procedures (alpha = 0.05; tree constraint for all: a node can be rejected only if its parent was)

- **B7** fixed Bonferroni: every position at `alpha / 7`.
- **B3** reference, the current family test: root and depth-1 nodes at
  `alpha / 3`; depth-2 nodes not tested.
- **AI** alpha-investing (Foster-Stine form as tabulated by Ramdas et al.
  2017): wealth `W0 = alpha`, reward constant `B0 = alpha`. Eligible nodes are
  tested in breadth-first order; the root is eligible, and a rejected node's
  internal children in the family become eligible. Penalty
  `phi_t = W(t-1) / 2`, level `alpha_t = phi_t / (1 + phi_t)`, wealth
  `W(t) = W(t-1) - phi_t + R_t (phi_t + B0)`. Target: mFDR.
- **CA** chain-conditioned alphas: `alpha_v = alpha * c_d * n_v / n`,
  `c = (1/3, 1/3, 1/3)` for depths 0, 1, 2. `n_v` and depth are fixed by the
  node's event, and nodes at one depth are disjoint, so the alphas of any tree
  sum to at most `alpha` (FWER by the union bound, given each node's
  conditional validity).
- **LA** learned from node features: `alpha_v = alpha * c_d * (lambda_d * n_v / n
  + (1 - lambda_d) * 2^-d)`. Both allocations sum to at most one per depth (at
  most `2^d` nodes at depth d), so any mixture keeps FWER at most `alpha`.
  `c` on the simplex grid of step 0.05 and `lambda_d` in {0, 0.25, 0.5, 0.75, 1}
  are chosen on the training half to maximise the total number of rejected
  mixed nodes; ties go to the earliest grid point. Applied unchanged to the
  test half.
- **LM** exploratory manifold weighting (no validity argument): within each
  depth, `c_d` (CA's) is split across the depth's nodes in proportion to their
  merge height in the diffusion tree, a feature computed from the full data and
  not fixed by any node's event.

AdaPT and STAR are not run: their FDP estimates `(1 + A_t) / R_t` and
`(h(1) + sum h(p_i)) / (1 + |R|)` (with `h(1) >= 1`) cannot fall to 0.05 with
fewer than 20 and 19 rejections, and a family has at most 7 nodes.

## Primary endpoints (test half)

Per procedure and panel: FWER (share of families with any null node rejected),
mean FDP, mFDR (total false over total rejections), share of mixed nodes
rejected overall and by depth. One-sided 95% exact binomial upper bounds for
FWER.

## Tree-assumption tests (all families, null nodes)

- **T1 (alpha-investing, SuperCoAD):** for null nodes at depth 1 or 2, compare
  `P(p <= 0.05)` and `P(p <= 0.01)` given parent `p <= 0.05` versus parent
  `p > 0.05` (Fisher exact test, one-sided binomial test of the
  parent-rejected rate against 0.05), and the Spearman correlation between
  parent and child p-values; separately for null parents and mixed parents.
- **T2 (STAR/AdaPT independence):** Spearman correlation between the p-values
  of two null siblings; the null p-value histogram in ten bins.
- **T3 (features):** Spearman correlation between null p-values and depth,
  `n_v / n`, child balance `min(child) / n_v` (all fixed by the event), and the
  merge height relative to the root (not fixed); null rejection rate at 0.05 by
  feature tertile.

An assumption counts as contradicted when its test gives p <= 0.05; small
counts are reported with their bounds rather than as a pass.

## Limits fixed in advance

Known-noise Gaussian rows only; depth at most 2; no production traversal,
counts, fitted PCA or unknown covariance. One run, no reruns with new seeds.
