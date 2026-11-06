# Adaptive alphas on the tree: three procedures and the tree assumption

Work on 2026-11-06 on `dev` at `69360b99`. Protocol fixed before any result:
`adaptive-alpha-tree-protocol-20261106.md`. Library check:
`adaptive-alpha-tree-libraries-20261106.md`.

## Design

650 known-noise Gaussian families (sigma 1, 5 features, 999 chain steps, alpha
0.05), each the internal nodes of depth at most 2 of the native diffusion/average
tree (up to 7 nodes), tested with `_node_chain_test` conditional on the ancestor
path. Panels: N30 global null (300), R30 root signal with null children (150),
H40 2 x 2 hierarchy (100), H80 2 x 2 x 2 hierarchy (100). Even replicates train
the learned allocation; odd replicates are the test half for every procedure.
All procedures reject a node only after its parent.

- B7 fixed `alpha / 7`; B3 the current family test (`alpha / 3`, depth at most 1).
- AI alpha-investing: `W0 = B0 = alpha`, `phi = W / 2`, level `phi / (1 + phi)`,
  breadth-first. Target mFDR.
- CA alphas fixed by each node's event: `alpha * c_d * n_v / n`, `c = 1/3` each.
  FWER at most alpha by the union bound.
- LA learned from node features (depth, size): `alpha * c_d * (lambda_d * n_v / n +
  (1 - lambda_d) * 2^-d)`; training chose `c = (0.35, 0.45, 0.20)`,
  `lambda = (1, 0.25, 1)`. FWER at most alpha by the same bound.
- LM exploratory: depth share split by diffusion merge height (not fixed by any
  event).

Run: 650/650 families, 0 unresolved or stalled nodes, 0 family errors; median 23 s
(30 rows) to 37 s (80 rows) per family, 42 minutes on seven processes. The log
holds only scikit-learn's notice that loky runs single-job inside
multiprocessing. Results SHA-256
000ed184303d9951f562edb2a4c9fb8b32f1755241113fb5285105b389a28f0f.

## Primary endpoints (test half)

Families with any false rejection, and mixed nodes rejected:

| Procedure | N30 false | R30 false | H40/H80 false | R30 mixed | H40 mixed | H80 mixed | All mixed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B7 | 1/150 | 0/75 | 0/100 | 43/141 | 51/292 | 135/347 | 229/780 |
| B3 (current) | 2/150 | 1/75 | 0/100 | 49/141 | 54/292 | 137/347 | 240/780 |
| AI | 3/150 | 1/75 | 0/100 | 55/141 | 63/292 | 195/347 | 313/780 |
| CA | 2/150 | 0/75 | 0/100 | 49/141 | 53/292 | 141/347 | 243/780 |
| LA | 2/150 | 1/75 | 0/100 | 51/141 | 53/292 | 139/347 | 243/780 |
| LM | 2/150 | 0/75 | 0/100 | 49/141 | 54/292 | 137/347 | 240/780 |

- Alpha-investing rejects 30% more mixed nodes than the current test (313 versus
  240). Most of the gain is at depth 2 of H80 (52/197 versus 0 for B3) and depth 1
  of H40 (11/100 versus 4/100). Its global-null family-wise rate is 3/150 (upper
  95% bound 0.051); its guarantee is mFDR, not FWER.
- CA and LA keep provable FWER and gain only 3 mixed nodes over B3; learning the
  allocation from features adds nothing over the fixed one.
- LM gains nothing and has no validity argument.
- H40 depth-1 power stays at 3 to 11 of 100 for every procedure: the observed
  child splits are peels across the true subgroups, which no alpha can repair.

## Tree-assumption tests (all families, null nodes)

Null calibration: p at most 0.05 in 11/300 roots, 37/815 depth-1 and 66/1701
depth-2 null nodes.

- T1, alpha-investing's condition (null p-values super-uniform given earlier
  decisions). Null children rejected at 0.05 after their parent was rejected:
  3/92 under null parents (upper bound 0.082; 78/2067 when the parent was not
  rejected; Fisher p 1.0) and 6/190 under mixed parents (upper bound 0.061;
  16/167 when not rejected; Fisher p 0.015). Super-uniformity given a parent
  rejection is not contradicted. Independence is: parent and child p-values are
  negatively correlated (Spearman -0.14, p 1e-10, null parents; -0.12, p 0.02,
  mixed parents), so a rejected parent makes its null children less likely to
  reject. By the protocol's rule these are contradictions; their direction is
  conservative for top-down testing.
- Null children of an unrejected mixed parent reject at 16/167 (9.6%). Every
  procedure here leaves them untested, but a procedure that tests children
  without the parent gate would be anticonservative there.
- T2, STAR/AdaPT independence: two null siblings, Spearman 0.045 (p 0.13, 1,115
  pairs); null p-value histogram in ten bins 255, 313, 286, 282, 264, 279, 281,
  268, 279, 309.
- T3, features: null p-values are unrelated to depth (rho 0.005, p 0.80), size
  share (0.004, p 0.85) and child balance (0.002, p 0.93), the event-fixed
  features CA and LA use. Diffusion merge height, which no event fixes, is
  related (rho -0.061, p 0.001): taller merges give smaller null p-values, so
  manifold weighting tends to give more alpha to null nodes that already look
  significant. Its 0.05 rejection rates by tertile were 36, 41 and 37 of about
  939.

## Limits

Known-noise Gaussian rows, 30 to 80 rows, depth at most 2; no production
traversal, counts, fitted PCA or unknown covariance. One run with fixed seeds.
AI's mFDR target differs from the FWER of B3, CA and LA.
