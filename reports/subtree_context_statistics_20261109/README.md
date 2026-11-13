# Subtree alternatives, eigenvalue magnitudes and native statistics

Date: 2026-11-09. This study implements the confirmed comparison of the four
node-context cases, relates them to the previous spectral and calibration
results, and traces the actual native tests. Production algorithms, thresholds
and defaults were not changed.

## What the four cases mean statistically

A node C is its exact set of labeled leaves. Its original parent P and unordered
child partition A|B are recorded separately. Linkage numbers are not identities.

| Structural case | What happens to the original statistical comparison |
| --- | --- |
| 1. C disappears | There is no native test for that clade in the alternative tree. Record absence, not p=1 or non-rejection. |
| 2. Both original children disappear | If C survives, its new sibling test compares a different partition. If C also disappears, there is no corresponding test. |
| 3. C survives but its original parent disappears | The incoming C-versus-parent comparison changes. If A|B survives, the outgoing sibling contrast can remain the same, although fitted quantities and eligibility can change. |
| 4. C survives but its original parent and both children disappear | Both the incoming comparison and the outgoing sibling partition change; only C's data rows remain fixed. |

These cases overlap: case 4 is in cases 2 and 3; case 1 can overlap case 2.
The root has no incoming parent comparison. The trace also retains each child's
presence and immediate adjacency: a parent clade can survive higher in the tree
while the immediate parent changes.

The six-leaf witness is
`((((0,1),(2,3)),4),5)` versus `((((0,2),(1,3)),5),4)`.
Clade `0123` survives while its parent changes from `01234` to `01235` and its
children change from `01|23` to `02|13`. Distinct whole covariance eigendirections
can generate these Ward trees. The [structural report](structure/report.md)
gives the exact coordinates and shows why six leaves are necessary for case 4.

## What was measured

The [protocol](protocol.md) fixes 16 reused matrices, four from each previous
regime, with 64 observations and six features. Each receives two native reference
trees and eight covariance-subspace alternatives: six individual whole blocks,
a leading rank-at-least-two span, and full coordinates. All 160 fits completed.

Subspaces propose topologies. The native gate replay uses the original feature
data in every topology. All branches are refitted by native all-pair NNLS against
the same original mean-squared standardized Euclidean target. The source-space
distances are also saved; they are not substituted for that common target.
This design separates changes of tree from directly changing the tested data.

The 16,128 comparisons contain case counts 15,316 / 5,271 / 401 / 102 for cases
1 / 2 / 3 / 4. Of the case-2 rows, 399 have a surviving focal clade, including
200 roots. These are dependent node comparisons, not independent experiments
or estimates of population topology probabilities. The [native summary](native_summary.md)
stratifies node size so small terminal clades are not confused with factor recovery.

There are 19 nonroot comparisons retaining the same immediate parent and child
partition. In every one, the raw sibling statistic and branch factor change,
while degrees of freedom stay fixed. Thus identical local groups do not imply
identical native numerical evidence when the rest of the tree and its branch
fit change. For 296 nonroot comparisons retaining the same children, the median
absolute difference after multiplying raw T by its branch factor g is
2.53e-11. The 23 tolerance exceptions range from 7.68e-10 to 3.50e-7;
this is a diagnostic of approximate scaling, not an exact identity with ridge.
For 260 comparisons with changed children, T*g changes in 257: those scores
also compare different hypotheses.

Of 812 present alternative clades, 479 have a finite collected raw p-value but
no available active BH decision; 331 record active rejection and two active
non-rejection. Traversal visits 523 and leaves 289 unvisited. The trace retains
absence, BH availability, visits and final boundaries independently.
All 32 native reference fits reproduce the previous data hashes, root splits,
decisions and partitions. Raw covariance eigenvalues agree for the same clade;
the independent recomputation's maximum absolute error is 3.38e-14.

## Three spectra with different roles

| Spectrum | Where its magnitude enters | What it does not establish |
| --- | --- | --- |
| Diffusion operator | Retained coordinates have weights max(lambda,0)^t; squared-distance coefficients are max(lambda,0)^(2t), here lambda^6 for positive modes. Vector coordinates and degrees also determine actual distances. | A coefficient share is not a share of pair-distance energy, a Wald variance, or a probability of a topology. |
| Raw or guarded covariance | Covariance used for whitening determines how large the contrast is relative to estimated variation. Raw trace shares describe retained variation in declared feature units. | Large retained variance does not guarantee preservation of a target grouping. |
| Native null-whitened PCA context | Eigenvalues order directions, delimit numerical tie blocks, affect positive-rank availability and MP counts. | They do not directly weight or divide the final projected Wald sum once its basis and retained rank are fixed. |

For the continuous sibling path, with d = mean(A)-mean(B), the useful schematic is

```text
V0 = (1/n_A + 1/n_B) Sigma_hat
g  = 1 + (b_A + b_B)/(2 mean_positive_branch)
V  = g V0 + ridge I
z  = whiten(d, V)
T  = ||Q z||^2,     nominal df = retained rank
```

Native covariance blocks and their jitter policies refine this schematic.
The adaptive .90 rule uses cumulative squared contrast coordinates
`(q_i^T z)^2`, not PCA eigenvalues as weights. A boundary tie extends the retained
block. With Q=I and z=(3,1), changing PCA eigenvalues from (4,1) to (100,1)
leaves T=9, df=1. Tying them at (4,4) instead yields T=10, df=2. The exact
[Wolfram checks](math/wolfram_evidence.json) and
[native kernel controls](math/kernel_controls.json) support the algebra and
these conditional examples; they are not a calibration proof.

The actual incoming-edge basis uses the parent's basis halves. The sibling
test uses the all-leaf context returned by the edge pipeline. In the first
crossed fixture, the root's six all-leaf eigenvalues are about .984375, while
the edge-half spectrum includes 3.4315 and five values about .9375. The runner
captures both actual calls. Half-row basis construction does not make the whole
test independent: covariance fallback, full-node whitening and tree selection
remain relevant. The [code and mathematical audit](math/README.md) specifies
those conditions and the subsequent empirical inflation, BH and traversal steps.

## Connection to earlier findings

1. **Weaker directions can preserve an omitted partition.** In the 20 crossed
   datasets, PC2 alone recovers the omitted factor B with two-cluster AMI 1 in
   every dataset. Its median raw variance share is .2885. A leading span with
   median share .9300 and the full tree both omit B's factor split, while the
   full tree still contains the four joint groups. The best B-covering
   shortlisted tree has median common-target Euclidean RMSE 1.3962 higher than
   the full tree. This earlier Euclidean fitting objective is distinct from
   the native standardized target used in the current replay.
2. **Projection loss and incompatible partitions are different.** In the
   compatible regime, the leading span retains median variance .8022 but all
   four joint clades only in 6/20 datasets, versus 20/20 for the full tree.
   Crossing factor splits cannot both appear in one tree; compatible splits can
   be lost merely by truncating the coordinates. Details are in the
   [saved-data magnitude reanalysis](relationships/README.md).
3. **Inherited geometry is not rebuilt geometry.** Restricting parent scores,
   inducing the parent graph, and refitting a child graph are separate operations.
   The earlier induced/refitted operators differ in 316/320 comparisons. The
   current covariance comparisons use common feature coordinates and distinguish
   rank mismatch from equal-rank span difference. Review found three local
   rebuilt bases slicing numerical tie blocks; corrected span-only sidecar data
   retain whole positive blocks. Original native traces and the original runner
   remain archived, and no native fit was rerun for this correction.
4. **Geometry does not replace a statistical reference.** The original A/B/C
   audit's fitted-covariance, contrast-selected-rank, ceiling-3 and best-cut-tail
   limitations remain. Neither an eigenvalue nor a count of observed alternative
   topologies supplies the missing probability bound. The new null replay records
   root rejection on all four common-mean null matrices in each of ten methods;
   these are four independent null datasets reused across methods, not 40.

The earlier exact Gaussian root Beta reference fixes both the tested union and
label contrast, with a shared Gram-only subspace and rank. The deterministic
all-label correction covers label selection on that fixed union. Its critical
ratio tends to 3 as sample size grows at fixed rank; this is not a population
ceiling. Neither result applies automatically to a selected descendant or native
contrast-adaptive rank. A change of parent or child partition must enter the
chosen hypothesis family and selection event.

The [Chen–Gao polar construction](https://arxiv.org/html/2610.17458v1) supports
conditioning under independent Gaussian rows with a common covariance and the
strong common-mean null on the tested union. It does not make arbitrary selected
TBS gates calibrated. The earlier implementation now includes an opt-in native
raw-score replay under a declared Ward event and refits fixed-topology NNLS at
accepted states. Calibration of the jointly selected full diffusion, subspace
search, branch fit, eligibility and traversal procedure, and final-partition
control, remain unestablished. Bootstrap frequencies depend on the resampling law and
can reach alternatives outside a fixed geometric family.

## Consequences for production and benchmarks

The implemented diagnostic identifies which labeled partitions each subtree
preserves and whether a changed number belongs to the same contrast. Candidate
evaluation should retain the subspace and distance provenance, exact clades,
parent/child context, fitted branches, covariance/rank, raw score, correction,
gate status and traversal. This is information needed to assess a change;
it is not a new production selection rule.

The next production decision remains gated by separate evidence: preservation
of declared factor and joint partitions; recovery measured by AMI, matched F1
and pair split/merge errors; fit on a common declared distance target; runtime;
and a reference covering the actual selection procedure. A lower residual or
larger eigenvalue cannot stand in for all of them. The earlier failed calibration
and runtime prerequisites remain in force, so no preset was promoted and no
frozen validation panel was started.

## Reproduction and evidence

Run the native module with `--output` pointing to a new directory; existing
manifests are protected against overwrite. `summarize_native.py` reuses saved
traces, regenerates input hashes and checks parity; `relationships/run.py`
reanalyses the earlier 5,120 saved alternative rows without fitting new trees.
The original 160 fits, 16,128 comparisons, source hashes and all attempt statuses
are in `native/`. `source_during_native_run.py.txt` preserves the exact runner
used before the diagnostic tie correction. Final check commands, counts and
source/archive verification are recorded in `verification.json`.
