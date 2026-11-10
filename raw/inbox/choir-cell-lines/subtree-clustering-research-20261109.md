# Subtree clustering research: implementation and measured results

The benchmark and diagnostic implementation is complete for the development
stage. **No new clustering configuration is promoted.** The bounded topology
candidate fails the agreed runtime gate, its full production selection reference
is not established, and the conditional Ward family study has insufficient
demonstrated power. Production algorithm files and default presets are unchanged.

## Implemented capabilities

- Shared benchmark metrics now include one-to-one matched macro/small-group F1,
  pair-splitting and pair-merging shares, and an AMI performance grid. Existing
  macro-F1 semantics remain unchanged. Matched scores penalize pure fragmentation
  that majority-remapped scores can conceal.
- Statistical profiles describe known reference groups and explicitly supplied
  subtree splits separately: sample size/balance, effective sample size,
  covariance spectrum and entropy effective rank, Euclidean and pooled
  Mahalanobis separation, descriptive centroid overlap and supplied geometry
  disagreement. Singular covariance does not produce a misleading pseudoinverse
  distance. Reference labels never enter clustering.
- The candidate selector retains the current topology and adds one canonical
  Ward topology per complete eigenspace and leading span. It supports 2–256
  leaves, refits all leaf pairs against one declared target, preserves real Ward
  heights, and retains the current topology on a residual tie. A compressed
  split-Gram NNLS objective agrees with the dense reference for all 945 rooted
  six-leaf trees in each distance arm. Large-subtree coverage is incomplete.
- The native score attribution module separates branch factor, covariance,
  projection/rank and displayed empirical correction. The conditional family
  module tests node homogeneity under one common Ward subtree/ancestry event and
  uses Holm after deriving conditional validity for that fixed labeled family.
- The study runner records complete native decisions, raw records, fitted
  branches, covariance/rank profiles, annotations, traversal, source hashes,
  data/truth hashes, status and complete-pipeline timing. Frozen evidence
  validators require exact coverage and explicit phase-2 provenance; missing,
  failed or unresolved attempts cannot produce a promotion pass.

## Completed experiments

| Evidence | Completed work | Interpretation |
| --- | --- | --- |
| Native factorial | 2,400 fits: 50 datasets in each of a Gaussian null, balanced signal and imbalanced signal setting; 16 operator/linkage/branch/energy combinations | Dataset-paired attribution |
| Recovery grid | 1,600 fits: two registered methods on 800 datasets, 16 statistical-profile cells, 50 seeds each | Development comparison, not frozen validation |
| Unique inputs in those panels | 850 matrices; the 100 signal datasets are shared across factorial and recovery analyses | Do not count reused inputs as independent replications |
| Conditional/attribution development | 300 final dataset records and four matched native-orbit illustrations | Restricted Gaussian reference and power evidence |
| Geometry | Eight larger all-pair fits; exhaustive small-tree parity checks | Distance fit and bounded candidate coverage |
| Perturbation laws | 10,000 six-leaf draws, 2,000 under each declared law | Sensitivity, not confidence probabilities for a true topology |
| Runtime | 36 isolated complete-pipeline runs on six inputs; initial 36 fits preserved separately | Three runs per method/input; repeated timings are not independent recovery datasets |
| Saved-data profiles | 1,520 CHOIR and 128 synthetic per-group recovery rows; 209 selected CHOIR subtree profiles | Descriptive reanalysis, not new clustering fits |

Every native factorial and recovery fit completed with software status `ok`.
That status does not establish a calibrated statistical reference.

## Calibration and mechanisms

On the 50 independent n=64,p=12 single-population Gaussian datasets, both native
symmetric/average NNLS and Markov/Ward NNLS with .90 energy truncation rejected
the selected root in **47/50** cases at nominal sibling alpha .01. Each returned
more than one cluster in 48/50 cases. The exact 95% interval for 47/50 is
[.835,.987]. The first fixture reproduces through the canonical dispatcher.
This is a substantial calibration failure in the tested setting, not a claim
about every input distribution.

Within fixed Ward topology construction and energy truncation, replacing linkage
lengths with NNLS added 41 root rejections in the symmetric arm and 44 in the
Markov arm; these reversals entered from an unavailable sibling p-value. The
traces show changed branch factors and edge eligibility. A separate same-root
comparison kept raw T, degrees of freedom and edge p-values unchanged while the
empirical correction changed from about 4 to 1, turning a boundary into a split.
These effects explain why examining only the raw sibling statistic is inadequate.

Linkage heights are a diagnostic control, not established stochastic variance.
The no-energy arm retains the native edge/MP projection; it is not a fixed-rank
Gaussian Beta test. See `calibration/factorial_cell_0_audit.md` for independent
hash, dispatcher and statistic reconstructions.

Global positive rescaling of branch lengths preserves
g=1+(bL+bR)/(2 mean_positive_branch). The ideal unregularized score scales as 1/g
with fixed covariance and projection. Production adds a fixed ridge after
variance scaling, and adaptive rank can also change; the exact g=1 control must
recompute the score. The initial overly strict ridge-free check and corrected
same-input rerun are retained in the calibration report.

The common-event Holm argument requires independent Gaussian rows with a common
mean and covariance in each true-null node, independent of its outside rows.
It covers false rejections of those homogeneous nodes, including partial-null
families with heterogeneous outside means. It does
not establish equal-average-only inference, final-partition error control, or
validity for the new diffusion/minimum-residual selection and production gates.
At 999 steps the mixed-root development panel had one false family rejection and
**no true family rejections in 50 datasets**. All 50 seven-hypothesis mixed-root
families were unresolved at 99 steps: the smallest attainable p-value, .01,
exceeds the first Holm cutoff, .05/7. At least 139 steps are required even to
reach that cutoff. These limitations prevent promotion.

## Recovery and statistical profiles

On the 16-cell development grid, Markov/Ward minus symmetric/average has mean
AMI difference **−0.05550**, with paired 95% Student interval
**[−0.06360,−0.04740]**. All 800 dataset pairs are complete. Nine cells show an
AMI decrease after Holm adjustment of the one-sided cell-loss tests. This is
an effect of the complete registered configurations, not an isolated operator
or linkage effect and not a general ranking on all benchmark distributions.

| Equal-weight grid mean | Symmetric/average | Markov/Ward |
| --- | ---: | ---: |
| AMI | .46986 | .41436 |
| ARI | .41434 | .31762 |
| Matched macro F1 | .59199 | .50048 |
| Matched small-group F1 | .49877 | .41475 |
| Pair split share | .41914 | .48448 |
| Pair merge share | .31989 | .34386 |

In the earlier four-seed imbalance data, the smallest planted group has nearest
centroid separation 5.108 on average, larger than the other groups, yet B/C leave
it unmatched on every Markov/Ward fit. Its sample covariance has rank 11 in 12
features. Size, separation, rank and information are distinct quantities.

CHOIR profiles use the actual archived 30-PC matrix. The 190 cell lines have
effective covariance ranks 5.714–20.784. B's mean matched F1 rises from .484 to
.969 across nearest-neighbor Mahalanobis separation quartiles and falls from
.947 to .510 across overlap quartiles. These are selected, confounded descriptive
strata. Cell-line labels are an external reference partition; they do not imply
that each line is one homogeneous Gaussian population. See `profiles/README.md`.

## Runtime and promotion decision

| Leaves | Features | Candidate / baseline median time | 5× gate |
| ---: | ---: | ---: | --- |
| 32 | 12 | 1.21× | Pass |
| 32 | 72 | 2.05× | Pass |
| 128 | 12 | 2.29× | Pass |
| 128 | 72 | 14.10× | Fail |
| 240 | 12 | 4.65× | Pass |
| 240 | 72 | 24.88× | Fail |

The isolated panel has stable source hashes and 36 successful runs. The widest
candidate evaluates 143 trees over 28,680 leaf pairs; exact compiled NNLS is the
dominant search cost. Its existing empirical reference is not repaired by the
topology search. On the single 240-row,72-feature preflight input, AMI drops from
.8633 to .5087; repeated timings do not turn that example into a power study.

`promotion_gate.json` records the rejected prerequisites. The 120,000-dataset
frozen Gaussian panel and 3,200-dataset frozen recovery panel were **not run**.
The full candidate development recovery grid was also not started after the
preflight failure, as specified before running the panels. This is not a claim
that those validation tasks passed or that the approved production objective was
achieved. No defaults were changed to force an apparent improvement.

## Task disposition

T0–T7: implemented protocols, shared metrics/profiles, native traces, geometry,
fits and restricted equation controls with tests and development evidence.
T8/T12: implemented and independently reviewed common-event Ward references;
full diffusion/search/gate/traversal inference remains unsupported and useful
mixed-root power is unestablished. T9: all 2,400 factorial fits completed.
T10: all 10,000 declared perturbations completed. T11: bounded all-pair search
implemented and verified, failed two runtime cases. T13: prerequisite gate
evaluated; frozen validation and production promotion were not entered.

The next mathematical work should address the jointly selected reference and
its usable information before further threshold changes. A faster candidate
policy would need a prespecified or proved pruning rule, its coverage measured
against the small exact oracle, and the complete search included in calibration.
Neither a smaller observed topology count nor a better distance residual is a
valid p-value correction by itself.

## Verification and provenance

The broad check passed 214 tests with two slow/optional cases deselected. A final
18-test validation rerun covers phase-provenance hardening; the union contains
216 distinct passing cases. Changed-file Ruff checks pass. Final verification
and wiki lint are recorded in `verification.json`.

The initial timing run saved all 36 fits before its summary writer encountered a
NumPy boolean serialization error. The fix only converts gate-result booleans to
Python bools. `source_change_review.json` retains the exact old source and AST
comparison. Isolated timings were then run with no other research panel active.
The calibration report separately preserves its ridge-check correction and
observed-numerical-failure guard. No original evidence was overwritten.

Reproduce new panels with `python -m benchmarks.diagnostics.calibration.subtree_research`
using `preflight`, `factorial` or `recovery`, a fresh `--output`, and one BLAS
thread. Existing outputs are not overwritten. `summarize.py` verifies counts and
paired hashes; `figures/plot.py` regenerates the PNG/PDF/SVG. Commands and seeds
are documented in `protocol.md` and each run's manifest. Frozen summary test
fixtures are deterministic software tests, not simulated calibration evidence.


Primary report location: reports/subtree_clustering_research_20261109/README.md.
Report-relative paths above resolve there. The report directory retains protocols,
all attempts, numerical corrections, query provenance and complete figure files.
