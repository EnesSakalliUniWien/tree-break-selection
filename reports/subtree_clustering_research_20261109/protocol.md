# Subtree clustering research protocol

Approved plan implementation, 2026-11-09, Europe/Vienna. This protocol is written
before the new development panels. Existing 9 November studies remain historical
evidence; their source archives and outputs must not be overwritten.

## Questions and constraints

Separate topology construction, branch-length fitting, covariance/rank selection,
empirical correction, conditional reference and traversal. Alternatives are
evaluated inside one declared subtree. Outside observations remain fixed in a
conditional replay; quantities depending on changed rows must be recomputed.

AMI is the primary recovery measure. Report matched group F1, pair splitting and
merging, ARI and statistical profiles. Small-group recovery is supplementary;
there is no 0.02 F1 acceptance margin. Truth labels never enter fitting.

Production candidates require a justified reference, successful frozen checks,
and complete-pipeline runtime at most five times the current method. Source and
data hashes, every attempt, failure, unsupported/unresolved case, and timing are
retained. A failed prerequisite stops promotion; it is not a successful test.

## Ownership and outputs

- T0/T3/T9/T13 coordinator: protocol, runner, traces, factorial panels, gate report,
  source snapshots and wiki integration.
- T1/T2 metrics Sol: shared metrics, aligned group/subtree profiles, reporting.
- T4/T5/T10/T11 geometry Sol: exact bounded oracle, scalable candidates, common
  target branch fits, invariance and resampling evidence.
- T6/T7/T8/T12 calibration Sol: variance/rank attribution, conditional score and
  family reference, mathematical assumptions and restricted calibration evidence.

Agents own disjoint source/test files. Existing working changes are preserved.
No installs, commits, pushes, services or default-method changes are authorized.
The existing `.venv` and project commands are used with BLAS threads fixed to one
for timed panels. Independent jobs use separate processes, never concurrent
instrumentation of native replay inside one process.

## Development panels

Data seeds use `SeedSequence([20261109, phase, panel, cell, replicate])`; phase 1
is development and phase 2 validation. Sampler and perturbation streams have
separate panel keys. No seed is chosen in response to a p-value.

Recovery: 240 rows, four square-vertex means, 16 Cartesian cells with p in
{12,72}, group sizes balanced or {132,60,36,12}, neighbor mean gap in {2.5,5},
and covariance identity or the existing 45-degree rotated diag(4,1/4,1,...).
There are 50 development datasets per cell. Initial methods are native
symmetric/average and Markov/Ward; the bounded-search candidate is added only
after correctness and runtime preflight. Profile and recovery methods see the
same rows. A new frozen candidate is not chosen using validation outcomes.

Factorial attribution: 50 paired datasets in each of three declared settings:
Gaussian null n=64,p=12; balanced square n=240,p=12,gap=5; imbalanced square
n=240,p=12,gap=5. Cross coordinate operator {symmetric,markov}, linkage
{average,ward}, branch source {linkage,NNLS}, rank {fixed/native contrast energy}.
Linkage and NNLS both retain their actual nonnegative branches; an explicit
zero-branch control belongs to the separate variance study. All changes refit
the same empirical-correction policy. Trace all internal nodes for replicate 0.

Geometry examples include square, simplex, rectangle, repeated eigenspaces,
zero-edge contraction and root-edge identifiability. Exhaustive topology and cut
claims remain limited to the existing six- and ten-leaf bounds. Large-subtree
candidate coverage is explicitly incomplete. Refit candidates against a common
target and all leaf pairs for supported n<=256; do not compare residuals across
target metrics or silently substitute a sampled pair set.

Runtime preflight uses n=32,128,240 and p=12,72, balanced gap5, one fixed data
instance per cell and three complete runs per method. Every method includes its
own geometry, selection, NNLS and traversal. A candidate exceeding the 5x limit,
failing a numerical/representation invariant, or lacking the required selection
argument remains a diagnostic. Its expensive frozen panel is not started.
Development evidence for rejected candidates is retained, with unrun cells
explicitly marked as not run because the prerequisite failed.

## Frozen validation, only after a candidate lock

The lock identifies source hashes, candidate parameters, supported input model,
selection event, statistic, family correction and fixed sampling budget.
No validation command may run without a complete passing prerequisite record.

- Recovery: 200 fresh phase-2 datasets in each of the same 16 cells.
- Gaussian calibration: 24 cells, n in {32,128}, p in {4,40}, covariance in
  {identity,AR(0.8),diagonal geometric spectrum from 1 to 100}, and root versus
  prespecified canonical first-child path. Use 5,000 independent datasets per
  cell, thresholds .01 and .05. The data, not the node or chain position, is the
  replication unit.
- Report exact binomial intervals and Holm-adjusted one-sided binomial tests
  for excess rejection across the 48 declared cell/threshold comparisons.
  Passing a finite panel does not prove a probability bound.
- Report unresolved and failed attempts separately and include all attempts in
  rate denominators. Never describe a failed call as a completed non-rejection.
- Use paired dataset differences, equal cell weights, confidence intervals and
  multiplicity-adjusted profile comparisons for recovery. AMI improvement must
  have a positive confidence bound; a detected profile deterioration blocks
  promotion. Report missing outcomes and paired-support coverage.

Count, heavy-tail, unequal-covariance and cancer-inspired inputs are robustness
panels, not evidence extending a Gaussian theorem. CHOIR is an already-examined
application dataset, never a fresh confirmatory validation set.

## Completion and publication

Every task reports supported, rejected or inconclusive. A theorem or a verified
counterexample is an acceptable mathematical research result. A new production
configuration is added only if all prerequisite and frozen gates pass. Otherwise
the production change is explicitly not promoted, with the failed gates and
remaining work recorded. Benchmark metrics and verified diagnostic capabilities
can be integrated independently of a new clustering preset.

Final evidence goes into this report directory and cited raw source captures;
durable synthesis uses only the existing wiki, its index and append-only log.
Run focused tests, relevant pipeline/regression checks, changed-file Ruff checks
and `make wiki-lint`. Do not claim unrun experiments or general partition error
control from local Gaussian homogeneity tests.

## Terminology clarification after the development panels

The factorial's `energy=False` arm disables contrast-energy truncation. It still
uses the native edge/MP projection dimensions, which can be data dependent. The
earlier shorthand "fixed" does not mean a prespecified Gram rank in this arm.
The independent fixed-Gram-rank control belongs to T7. No settings, datasets or
decisions changed with this clarification; figures use "No energy truncation".
