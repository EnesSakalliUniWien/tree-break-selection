# Branch-length internal distributions and variance comparison

All 48,600 runs completed on 900 matched Gaussian dataset pairs: 32,400 current
single-matrix runs and 16,200 independent-equation runs. The effects depend on
how internal distributions are constructed. Enabling empirical internal rows
can increase null false partitions considerably; branch-length internal states
are different and do not produce a uniform improvement in signal recovery.
No production defaults or calculations were changed.

## Current diffusion NNLS

At sibling alpha .01 and edge alpha .001, with passthrough enabled:

| Internal spectral context | Edge variance policy | Null false partitions: identity2 / ar1_8 / identity40 | Signal mean ARI: identity2 / ar1_8 / identity40 |
|---|---|---|---|
| leaf_only | none | 16 / 0 / 0 | .02701 / .03773 / .00000 |
| empirical_barycenter | none | 13 / 23 / 37 | .02701 / .06814 / .06989 |
| branch_length_state | none | 15 / 9 / 0 | .02701 / .06370 / .00000 |
| leaf_only | normalized_branch_length | 6 / 0 / 0 | .00000 / .00000 / .00000 |
| empirical_barycenter | normalized_branch_length | 3 / 4 / 25 | .00000 / .04451 / .04992 |
| branch_length_state | normalized_branch_length | 4 / 0 / 0 | .00000 / .02427 / .00000 |

Every null count has denominator 250; every signal mean uses 50 datasets.
False partition means more than one cluster under the global null. The
empirical-barycenter high-dimensional result is 25/250 with edge variance,
compared with 0/250 using leaf-only or branch-length-state spectral rows.
The additional signal recovery is accompanied by additional false partitions;
these results do not establish a calibrated improvement.

With edge variance enabled, branch-length states change 4/900 final partitions
relative to leaf-only; empirical barycenters change 46/900. Holding branch-length
states fixed, enabling edge variance changes 27/900 partitions. These are matched
changes in full partitions, not just cluster counts.

Without NNLS, branch-length states give null false-partition counts 1/0/1
without edge variance and 0/0/0 with it; signal mean ARI remains zero in all
three profiles. Raw linkage branch variance uses an explicit diagnostic
override and is not promoted as a production setting.

## Independent equations with discovery-only internal states

These rows use the experimental report adapter described in
[compatibility.md](compatibility.md), with known covariance and no extra
branch-length variance. The adapter changes discovery PCA only; inference data
remain independent. It preserves the leaf-only candidate-rank cap.

| Branch source | Internal spectral context | Null false partitions: identity2 / ar1_8 / identity40 | Signal mean ARI: identity2 / ar1_8 / identity40 |
|---|---|---|---|
| diffusion | leaf_only | 0 / 0 / 0 | .64529 / .25120 / .15531 |
| diffusion | empirical_barycenter | 0 / 0 / 0 | .64529 / .26448 / .16442 |
| diffusion | branch_length_state | 0 / 0 / 0 | .64529 / .25707 / .16442 |
| diffusion_nnls | leaf_only | 0 / 0 / 0 | .64529 / .25120 / .15531 |
| diffusion_nnls | empirical_barycenter | 0 / 0 / 0 | .64529 / .26448 / .16442 |
| diffusion_nnls | branch_length_state | 1 / 0 / 0 | .64529 / .25707 / .14621 |

Branch-length states make NNLS relevant to the independent equation through
the discovery basis: NNLS changes 3/900 final partitions relative to diffusion
with branch-length states. With leaf-only or empirical-barycenter discovery
rows, NNLS changes none. Within NNLS, branch-length states also change 3/900
partitions relative to leaf-only. These counts refer to different comparisons.
Across 27,900 matched nodes, NNLS with branch-length states changes 6,887
p-values and 1,988 selected ranks. Thus few changed final partitions can hide
substantial node-level changes. The other two independent contexts have no
NNLS p-value or rank changes.

At root p-value cutoff .01, branch-length states yield 3/2/5 rejections out of
250 null datasets in the three profiles for either branch source, versus 3/3/6
for leaf-only. None of the prespecified 72 independent root-tail cells rejects
excess-tail behavior after Holm correction (minimum adjusted p=1). This is
finite simulation evidence with limited tail precision, not proof of calibration
or whole-tree error control. Signal changes are descriptive means over 50 cases.

## Where the p-values change

- Internal distributions change the spectrum and projection. Internal rows
  remain dependent combinations of leaves and do not increase independent row
  counts. The current method's same-data covariance/MP selection remains a
  different statistical construction from discovery-only known-covariance PCA.
- `edge_branch_length_variance_policy` controls child-parent edges. Sibling
  branch variance remains active even with the edge policy set to `none`.
  Changing this flag can change edge support, empirical inflation, which sibling
  tests are available, and final partitions without changing raw sibling tests.
- Under the independent generator, adding a branch variance multiplier m
  yields statistic T/m. Feeding it to the original survival function is a
  scale mismatch. The correct reference gives S(m * T/m) = S(T), preserving
  the original p-value. The recorded uncorrected tail is a negative control,
  not a production p-value or a replacement degrees-of-freedom formula.

For NNLS branch-length states, that negative control reduces root rejections at
.01 from 3/2/5 to 0/0/3. The smaller counts reflect the changed scale; they do
not demonstrate improved calibration. [The derivation](compatibility.md)
also explains why the scalar correction applies to each child-parent edge.

For current NNLS with branch-length states, toggling edge variance leaves all
27,900 raw sibling statistics and p-values unchanged. It changes gate-p-value
availability at 374 nodes and changes all 28 p-values available in both runs.
Reconstructing the empirical inflation fit gives 8,171 changed factors across
all records; many records are support-only or skipped and do not supply an
active sibling gate. The separate counts avoid treating those factors as
8,171 changed active tests.

Alpha sensitivity is material: with NNLS branch-length states, identity2 null
false partitions at sibling alpha .005/.01/.05 are 7/15/78 without edge
variance and 2/4/18 with it (each out of 250). The corresponding independent
variant has 1/1/1, but edge gates can limit final partitions. These counts do
not replace raw-tail checks.

## Coverage, provenance and remaining limitations

The full crossing uses two branch sources, three internal contexts, two current
edge-variance policies and sibling alpha .005/.01/.05. The independent adapter
uses three internal contexts for both branch sources at those alphas. Edge alpha
is .001 and passthrough is true throughout. Covariance profiles are identity in
2 and 40 dimensions, and AR1(.6) in 8 dimensions. All cases have 32 leaves.
Only average linkage and Euclidean diffusion are covered here.

The adaptive diffusion dispatcher currently does not forward the internal-state
flags. This benchmark calls the existing distance runner with the actual adaptive
diffusion distances and verifies the resolved settings; it does not silently
assume extra dispatcher parameters were honored. That dispatcher limitation
remains unfixed. The independent internal-row adapter is report-local, not a
new public API option.

The evaluation reuses the earlier 900 pairs and verifies data hashes. It is an
extension of that study, not an independent confirmation. No tuning followed
evaluation. Current final-gate root tails are conditional on availability;
missing/skipped sibling p-values stay unavailable. Raw root tails are separately
recorded. Do not pool dependent nodes as independent null replicates or treat
unsupported tests as p=1. Full alpha sensitivity, matched effects and Wilson
intervals are in the tables below.

The first pilot exposed a report-export `KeyError('L0')`: spectral metadata
included leaves that had no internal-row count. Its original records are kept
in `pilot/` and are superseded by `pilot_reviewed/`, whose 486 runs and 5,022 node
records passed verification. The exporter now restricts node diagnostics to
internal nodes. Production code was unaffected. All 49 focused regression
tests pass; [test output](verification-tests.txt) is saved.

Evaluation verification covers 502,200 node records, 5,400 frozen independent
contexts, and 7,200 reproduced baseline run rows from the previous benchmark.
There are no execution errors or captured warnings. Correcting the negative
control's reference scale restores p-values with maximum absolute difference
1.72e-15. Source hashes remained unchanged throughout evaluation.

## Reproduction and evidence

```sh
TBS_N_JOBS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python -m reports.branch_distribution_variance_20260928.run_benchmark \
  --phase evaluation --workers 6 --output-dir /tmp/branch-distribution-reproduction
.venv/bin/python -m reports.branch_distribution_variance_20260928.analyze_results \
  /tmp/branch-distribution-reproduction
```

The output directory must not exist. Baseline verification uses the previous
evaluation saved in this repository. The script defaults to 250 null and 50
signal replicates per profile. Pilot seeds differ from evaluation seeds.

- [Frozen protocol](protocol.md), [mathematical compatibility](compatibility.md).
- [Summary and intervals](evaluation/summary.csv), [root tails](evaluation/root_tails.csv).
- [Matched partition effects](evaluation/run_effects.csv), [matched node effects](evaluation/node_effects.csv).
- [Run records](evaluation/runs.csv), [node records](evaluation/nodes.jsonl.gz),
  [actual configurations](evaluation/contexts.jsonl), [manifest](evaluation/manifest.json).
- [Verification](evaluation/verification.json), including baseline reproduction,
  matching, warnings, frozen contexts and reference-scale restoration.
