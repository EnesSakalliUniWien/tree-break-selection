---
title: Cancer Cell-Line TBS Recovery Plan
type: question
status: draft
updated: 2026-11-06
sources:
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-20261106.md
  - raw/inbox/choir-cell-lines/tbs-merged-null-rerun-20261106.md
  - raw/inbox/choir-cell-lines/tbs-merged-cell-lines-pvalues-20261106.md
  - raw/inbox/choir-cell-lines/tbs-topology-comparison-20261106.md
  - raw/inbox/choir-cell-lines/tbs-topology-comparison-script.py.txt
  - raw/inbox/choir-cell-lines/tbs-markov-ward-full-run-20261106.md
  - raw/inbox/choir-cell-lines/tbs-markov-ward-5000-summary-20261106.json
  - raw/inbox/choir-cell-lines/tbs-nnls-lbfgsb-switch-20261106.md
  - raw/inbox/choir-cell-lines/tbs-nnls-solver-benchmark-20261106.json
  - raw/inbox/choir-cell-lines/tbs-lbfgsb-null-summary-20261106.csv
  - reports/root_conditional_sampler_20261104/README.md
  - reports/selective_rebuild_comparison_20261104/README.md
  - raw/inbox/choir-cell-lines/tbs-remedy-proposal.md
  - raw/inbox/choir-cell-lines/tbs-remedy-proposal-evidence.json
  - raw/inbox/choir-cell-lines/diffusion-maps-coifman-lafon-2006.pdf
  - raw/inbox/choir-cell-lines/tbs-diagnosis.md
  - tree_break_selection/tree/construction/defaults.py
  - tree_break_selection/space_separation/README.md
  - reports/calibration_restoration_20260909/validation_protocol.md
  - raw/inbox/choir-cell-lines/tbs-recovery-comparison.md
  - raw/inbox/choir-cell-lines/tbs-recovery-protocol.json
  - raw/inbox/choir-cell-lines/tbs-recovery-comparison.csv
  - raw/inbox/choir-cell-lines/tbs-recovery-validation-null.csv
  - raw/inbox/choir-cell-lines/tbs-recovery-validation-signal.csv
  - raw/inbox/choir-cell-lines/tbs-recovery-verification.json
  - raw/inbox/choir-cell-lines/tbs-nnls-convergence-diagnosis.md
  - raw/inbox/choir-cell-lines/tbs-nnls-convergence-evidence.json
  - raw/inbox/choir-cell-lines/tbs-nnls-convergence-comparison.csv
  - raw/inbox/choir-cell-lines/tbs-root-groups.md
  - raw/inbox/choir-cell-lines/tbs-root-groups-evidence.json
  - raw/inbox/choir-cell-lines/tbs-root-groups-composition.csv
  - raw/inbox/choir-cell-lines/tbs-component-statistics.md
  - raw/inbox/choir-cell-lines/tbs-component-statistics-evidence.json
  - raw/inbox/choir-cell-lines/tbs-diffusion-eigenvector-trace.md
  - raw/inbox/choir-cell-lines/tbs-diffusion-eigenvector-trace-evidence.json
  - raw/inbox/choir-cell-lines/tbs-diffusion-eigenvector-trace.csv
  - raw/inbox/choir-cell-lines/tbs-statistical-logic-trace.md
  - raw/inbox/choir-cell-lines/tbs-statistical-logic-trace-evidence.json
  - raw/inbox/conditional-diffusion-statistics-20261104.md
  - raw/inbox/conditional-small-benchmarks-review-20261104.md
  - raw/inbox/cancer-pc-surrogate-review-20261104.md
  - reports/cancer_pc_methods_20261104/README.md
  - raw/inbox/diffusion-attribution-20261102.md
tags:
  - choir
  - diffusion
  - validation
  - recovery
---

# Cancer Cell-Line TBS Recovery Plan

## Question

Can a controlled coordinate/linkage change avoid the singleton root failure and
improve recovery while meeting a defined, independently evaluated error target?

## Current State

The executed first experiment is a two-by-two comparison: current symmetric
versus Markov right-eigenvector coordinates, each with average versus Ward
linkage. The current symmetric/average case is the existing one-cluster baseline.
All cells, PCA, kernel and retained eigenpairs were reused. The remaining TBS
settings were fixed and NNLS was refitted separately for each topology. No
cluster count was forced, cell removed, alpha relaxed or Tree-BH bypassed.

For positive kernel degrees, the candidate is Phi_M = D^(-1/2) Phi_S. A read-only
check gives maximum Markov eigen-relation residual 4.51e-15 and reduces ASPC1
coordinate-norm CV from 0.190099 to 0.050017. This is geometric evidence only.
Ward is already supported and must use Euclidean coordinate distances; it does
not guarantee balanced roots or improved results.

All three new full-data NNLS fits failed to converge under the frozen
1,000-iteration budget; the native runner refused to apply them. They have no
valid TBS partitions or gate outcomes. Topology diagnostics remove the singleton
root: symmetric/Ward gives 250/48,629, Markov/average 151/48,728 and Markov/Ward
267/48,612. These topology changes do not establish clustering recovery.

The smaller root children are pairwise disjoint and each contains one cell
line. Symmetric/average separates one ASPC1 cell from its 88 same-line graph
neighbors. Symmetric/Ward isolates the entire 250-cell A204 graph component.
Markov/average separates 151/161 OAW28 cells, exactly CHOIR cluster 103, and
Markov/Ward separates 267/270 NCIH727 cells. The last two cuts cross two and
nine graph edges, respectively. The larger groups each contain 189 or 190
cell lines. The frozen NNLS sample contains 1, 567, 326 and 574 crossing pairs,
respectively. Both root-edge columns equal the crossing-pair indicator in
each design, so the predictions use the combined root length and the linkage
ratio determines its division between the two edges.

Independent full-pipeline validation completed 600 homogeneous count nulls and
35 signal mixtures, evaluated with every variant (2,540 successful method
evaluations). False final partitions under each 300-dataset null are:

| Variant | Low library-depth variation | High library-depth variation |
| --- | --- | --- |
| Symmetric / average | 0/300 | 262/300 |
| Symmetric / Ward | 295/300 | 228/300 |
| Markov / average | 34/300 | 298/300 |
| Markov / Ward | 297/300 | 300/300 |

The endpoint was P(final K > 1), with frozen target 0.01 and separate one-sided
95% exact binomial bounds. No variant meets the target across both null models.
Symmetric/Ward improves balanced-mixture mean ARI from 0.700000 to 0.988412,
while its null failures prevent adoption. The Markov variants often
overpartition signal mixtures. No new variant was adopted or production default
changed. Nulls have 400 cells and signals 1,000–3,000; this study does not certify
48,879-cell behavior, selected p-values or descendant FDR.

The follow-up NNLS diagnosis reproduced all three 1,000-iteration failures.
Changing only the inner LSMR tolerance to 0.001, with the initial solve and
outer settings fixed, reduced the stops to 311–388 iterations and substantially
lowered cost. All stopped on relative cost change. The accepted original
baseline has scaled first-order residual 2.8945, above the separate `1e-5`
criterion. All 2,540 earlier simulations also stopped on cost change with
reported optimality above that criterion.

An independent bounded solver lowered the objective on each of the four
unchanged sampled systems and produced numerical relative duality gaps below
`2e-5`. Those systems have 56,414 equations for 97,756 edge lengths. Feasible
perturbations outside the root pair preserve all sampled predictions even near
the independent solutions. Accurate optimization and identifiable individual
lengths are separate requirements. A production acceptance rule, stable branch
lengths and robustness of the gates to equivalent fits remain open.

The frozen-tree diagnostic replay applied the independent lengths for each
alternative and ran the existing gates. All three split at the root and each
returned the entire A204 and IM95 components as separate clusters. Every final
alternative cluster stayed inside one graph component. The baseline stopped at
its ASPC1 singleton root and mixed all three components in one cluster. A204
and IM95 were intact subtrees in all four topologies. Direct PCA sibling Wald
statistics for the three component pairs were 403,636, 23,981 and 50,158 on 30
degrees of freedom. The unequal-covariance check also showed strong observed
separation. The p-values use nominal references on selected groups.

The target is one hierarchy. The saved diffusion spectrum has three stationary
modes, one per graph component. All 27 retained nonstationary modes come from
the main component. A204 and IM95 have local second eigenvalues 0.936762 and
0.863919, below the global cutoff 0.999052, so their saved Markov coordinates
have only roundoff-level variation within each component. Both left and right
eigenvectors have maximum relative residual below `2e-14`. The stationary
subspace accounts for 99.9997% of the squared centroid separation at the
symmetric/Ward root. Mode 17 accounts for 95.02% at the Markov/Ward root.

A single-tree candidate can allocate nonstationary modes within components,
preserve each as a clade, and define upper merges using common PCA geometry.
Between-component contrasts need independent measurements with frozen
membership or null calibration that accounts for the graph selection. The
candidate tree construction and its error control remain to be tested.

The full-data and null results expose distinct problems. The alternatives
originally stopped before testing at the NNLS guard. The baseline stopped at
a within-line singleton split, even with branch variance removed. The same
gate code separated the disconnected components in the alternative replays
with different root groups and diagnostic lengths. Those replays jointly
change the topology and fitted lengths.

The instrumented symmetric/average run computed 97,756 edge statistics, but
Tree-BH reached only the two root edges. Both failed, so the sibling gate was
skipped despite a raw root record. On the Markov/average root with diagnostic
lengths, both edges passed. With its raw sibling Wald statistic and empirical
model fixed, the committed local inflation rule gave `c_hat=746.50` and
`p=0.999989`. The working-tree family fallback gave `c_hat=20.34` and
`p=1.50e-28`. The root exceeds the largest calibration parent and local
effective sample size is 1.239. The comparison isolates the nominal root
sibling rejection. The committed rule was evaluated at the root only.

In the fixed-tree Gaussian control, the nominal root sibling test with known
identity covariance rejected 50/50 selected-data nulls and 0/50 independent
nulls. Guarded covariance gave 8/50 rejections even on the independent data.
Selection and covariance reference mismatch are separate calibration problems
in that control. Any later change allowing exploration below a failed immediate
contrast requires a defined hypothesis family and calibration argument.

The conditional-inference prototype now refits native diffusion for each
Gaussian radial proposal while conditioning on one complete observed merge
history. The known-identity nine-row panel had zero conditional rejections
in 60 nulls and 40 signals. Every strong signal had zero accepted proposals
among 999, making p=1 through limited Monte Carlo resolution. Wolfram checked
the density and tail formulas. This prototype has separate observation-model,
event-computation and descendant-traversal gaps before a cancer-data application.
[[conditional-diffusion-monte-carlo-evidence]] records its assumptions and results.

The later [[selective-rebuild-comparison]] separates rare-event sampling from
the choice of event. On six of those strong examples, a random-position
reversible-chain rank detects three under full-history conditioning and five
under root-split conditioning at 0.05; ordinary proposals still accept none.
Known-I whole-procedure rebuilding detects six, while its fitted-covariance
variant detects none. This is a root-only development comparison, not cancer
validation or a simultaneous rule for descendants of the observed tree.

The root-split sampler is now available as an opt-in maintained diagnostic;
[[root-conditional-sampler]] records its API, explicit unresolved decisions
and validation. All 24 earlier chain traces match exactly. A fresh 60-case
panel gives 3/40 null, 4/10 moderate and 7/10 strong rejections at 0.05, and
the two registered raw Gaussian signal roots reject with known scales. The
cancer-data model, descendant rules and exact-real membership remain open.

The subsequent registered continuous benchmark check gives numerical root
p=0.001 for both the clear and moderate signal cases with their known noise
scales. Exact-real certification remains [0,1] because their graphs have four
and three components and repeated leading eigenvalues. Immediate-descendant
checks use full global-tree replays and a fixed three-node Bonferroni family.
The moderate mixed child still has zero accepted proposals. The certificate
code itself reconstructs only root contrasts. Count-model preprocessing and
nonspherical or dependent-row benchmark inputs require their own model treatment.

### Compact cancer-inspired benchmark

Four registered 496-by-30 cases now reproduce root failure with known
generating populations and paired random draws. The development case has
three components, a dense singleton root, covariance fallback and 988
ancestor-blocked edges. Restoring ordinary noise to the designated observation
changes one cluster to six. All five additional seeds retain this cluster-count
contrast; three reproduce the singleton root, while two stop at the sibling
calibration gate after significant edge tests. The latter provide a compact
reference for investigating the inflation decision separately.

The known heterogeneous row covariance prevents applying the current
common-sigma conditional prototype directly. This PC-stage model does not
reproduce counts, fitted cancer covariance, rare-component proportions,
full-scale NNLS iteration limits or a certified selected law.
[[cancer-pc-surrogate-evidence]] records the reproducible scope and controls.

The integrated `cancer_pc` suite has also been evaluated with all eleven default
methods across six seeds. Nine methods run; the two Hamming diffusion variants
skip continuous inputs. Across the five additional seeds, spectral clustering
with supplied K=6 retains ARI near 0.992 in both disconnected row controls,
while family-diffusion TBS retains the 0 versus 0.9922 contrast. The report and
plots preserve method-specific geometry and known-K qualifications. This
comparison does not validate a selected reference law or identify a cancer fix.

### Singleton-peel merge rerun

After merge `50455316`, the original symmetric/average baseline was rerun on
unchanged inputs. The linkage is byte-identical, with the same 1/48,878 root,
but traversal now passes through the singleton peel: 70 clusters, 148
significant edges, ARI 0.027593 and NMI 0.548219 against cell lines. One
cluster still holds 28,611 cells from 144 lines; 27 lines are recovered at
90%/90%. Rerunning the 600 count nulls, the baseline splits 2/300 low-depth
nulls (previously 0/300), both through a singleton root, so it no longer meets
the frozen 0.01 target; high-depth splits go from 262 to 264.
Most remaining clusters stop at small peels rather than singletons: 65% of
cells sit in clusters whose top split has a child of at most 10 cells, including
the giant cluster's 4 versus 28,606 split. [[choir-singleton-peel-rerun]]
records the comparisons, cell-line breakdown and p-value checks.

### Topology limit and the Markov/Ward method

The symmetric/average tree contains only 42 of 190 lines as clean (90%/90%)
subtrees, so no gate can recover the rest; TBS finds 27 of the 42. With the
same kernel, Markov coordinates give 134 clean lines under either linkage, and
Ward removes most small peels: a 190-cluster cut of the Markov/Ward tree scores
ARI 0.731. Markov coordinates are now a production option and
`tbs_family_diffusion_markov_ward` an opt-in method, not a default. Its earlier
count-null failures (297/300, 300/300) stand, and it oversplits the four
`cancer_pc` surrogates (10 to 17 clusters for six populations).

On all cancer cells the method stops at the registered 1,000-iteration NNLS
budget. With a one-off 5,000 budget the fit stops at 1,017 iterations and TBS
returns 1,015 clusters: ARI 0.5337 and NMI 0.8602 against cell lines, 45 lines
recovered, 83.7% of cells in clusters at least 90% one line, but a median of
four clusters per line. The NNLS solver (`lsq_linear` TRF with automatic LSMR
tolerance) stops about three times above an independent L-BFGS-B cost, and
the sampled design leaves edge lengths non-identified.

Production now fits NNLS by column-scaled SciPy L-BFGS-B and accepts a fit
only with a certified relative duality gap of at most the solver tolerance;
the registered budget is 5,000 iterations. CHOIR fits reach the optimum
(cost 535.35 instead of 1856.37 for the default method) in 8 to 18 s. Results
barely move: `tbs_family_diffusion` gives 75 clusters (ARI 0.0275, 27 lines),
Markov/Ward 954 (ARI 0.549, 41 lines), and the baseline count-null splits are
2/300 low-depth and 259/300 high-depth.

## Evidence

- `raw/inbox/choir-cell-lines/tbs-remedy-proposal.md` specifies the comparison,
  fixed settings, metrics, validation and limits.
- `raw/inbox/choir-cell-lines/tbs-remedy-proposal-evidence.json` captures the
  numerical coordinate check and primary references.
- `raw/inbox/choir-cell-lines/tbs-diagnosis.md` records the observed failure and
  why removing branch variance alone does not fix it.
- `raw/inbox/choir-cell-lines/tbs-recovery-comparison.md` records actual outcomes,
  independent count models, verification and limits; the protocol, CSV tables
  and verification JSON preserve the frozen settings and measured results.
- `raw/inbox/choir-cell-lines/tbs-nnls-convergence-diagnosis.md` records native
  replays, the single-parameter tolerance controls, independent objective checks
  and nonunique edge lengths; its JSON and CSV preserve numeric evidence.
- `raw/inbox/choir-cell-lines/tbs-root-groups.md` records the root children,
  their composition and graph boundaries, PCA distances and sampled NNLS fits.
- `raw/inbox/choir-cell-lines/tbs-component-statistics.md` records gate replays,
  direct component statistics and the proposed component policy.
- `raw/inbox/choir-cell-lines/tbs-diffusion-eigenvector-trace.md` records the
  paired left and right modes, local eigenvalues and their root contributions.
- `raw/inbox/choir-cell-lines/tbs-statistical-logic-trace.md` follows the
  root through edge tests, sibling calibration and traversal. The evidence
  JSON preserves both inflation decisions on the same statistic.
- `raw/inbox/choir-cell-lines/tbs-merged-rerun-20261106.md` records the
  singleton-peel merge rerun, held-fixed inputs, hashes and recovery counts.
- `raw/inbox/choir-cell-lines/tbs-merged-null-rerun-20261106.md` records the
  merged-code count-null rerun and its paired comparison with the original nulls.
- `raw/inbox/choir-cell-lines/tbs-topology-comparison-20261106.md` compares the
  four saved topologies by peel structure and clean cell-line subtrees, and
  records the opt-in Markov/Ward method with its surrogate smoke run.
- `raw/inbox/choir-cell-lines/tbs-markov-ward-full-run-20261106.md` records the
  full-data Markov/Ward runs at 1,000 and 5,000 NNLS iterations and the solver
  analysis.

## Links

- [[choir-singleton-peel-rerun]]
- [[cancer-pc-surrogate-evidence]]
- [[choir-cell-lines-data]]
- [[diffusion-nnls-calibration-open-points]]
- [[adaptive-diffusion-nnls-method-library-audit]]
