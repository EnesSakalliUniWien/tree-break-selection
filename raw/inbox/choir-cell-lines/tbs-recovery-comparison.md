# Cancer cell-line TBS coordinate/linkage experiment

The controlled comparison found no usable replacement under the frozen settings.
All three new full-data runs were refused by the native runner because their
fixed-topology NNLS fits did not converge. Independent count simulations also
showed false partitions far above the prespecified target under at least one
null model for every variant. No new variant was adopted.

## Frozen experiment

The protocol was written before examining the new clustering results. The
existing symmetric/average run was reused as the baseline. All 48,879 cells,
2,000 variable genes, 30 PCs, the kernel, retained eigenpairs, diffusion time 3,
NNLS target/pair-sampling seed and tolerances, covariance policy, projection,
edge alpha 0.001, sibling alpha 0.01, Tree-BH and traversal were fixed. NNLS was
refitted on the original PCA target for each topology, with a 1,000-iteration
budget and `apply_nonconverged=False`.

The Markov candidate is Phi_M = D^(-1/2) Phi_S, where Phi_S = U Lambda^3 and
S = D^(-1/2) W D^(-1/2). Ward and average both use Euclidean distances from their
explicit coordinates. No cells were removed, no cluster count was forced, and
no threshold was adjusted. Known labels were available only for diagnostic
scoring, not as inputs to the fit.

## Full cancer-data outcomes

| Variant | Root child sizes | Native TBS result | Top-100 singleton splits |
| --- | --- | --- | --- |
| symmetric_average | [1, 48878] | 1 cluster; ARI/NMI 0 | 12 |
| symmetric_ward | [250, 48629] | NNLS did not converge; no partition | 0 |
| markov_average | [151, 48728] | NNLS did not converge; no partition | 24 |
| markov_ward | [267, 48612] | NNLS did not converge; no partition | 2 |

Each alternative removed the singleton root. Those topology diagnostics do not
establish that TBS would recover the cell lines: the native runner refused each
NNLS fit before gate annotation or traversal. Their cluster counts, ARI/NMI,
fragmentation, mixed-cluster composition and gate decisions are unavailable,
not zero. The original baseline reached two root-edge tests, rejected neither
and ancestor-blocked the remaining 97,754 non-root nodes.

Detailed NNLS diagnostics were captured for Markov/Ward: 1,000 iterations,
maximum-iteration termination, cost 1494.830081, optimality 1.436924 and RMSE
0.230206. The fit was not applied to the tree. For the first two alternatives,
the original exception preserved the nonconvergence status but not detailed
solver iteration counts. Their topologies were subsequently rebuilt with the
same native `build_tree` and saved as topology-only diagnostics. Original
failure records were retained. No nonconverged fit was forced into testing.

## Independent validation

The frozen validation contains 635 independent count datasets: 300 replicates
for each of two homogeneous negative-binomial nulls, 10 balanced four-group
mixtures, 10 imbalanced four-group mixtures, 10 mixtures with an additional
singleton, and five density-variation mixtures. Every dataset was evaluated
with all four variants, giving 2,540 method evaluations; all completed
successfully. These small-data successes do not resolve the real-data solver
failures.

For each dataset, counts were generated afresh and converted to CPM, followed
by log1p(CPM/100), Seurat HVG selection, scaling, PCA, kernel/eigenpair fitting,
linkage, NNLS, covariance/projection estimation, both gates and traversal.
Paired variants used identical counts within a replicate. Fingerprints confirm
635 distinct count matrices and identical inputs across each four-way pairing.
Gene means and cell library factors were newly drawn per dataset; the two nulls
have library log-standard-deviations 0.25 and 1.0. All null cells share the same
biological expression profile and noise law. Full generator details, sample
sizes, seeds and limitations are in `protocol.json`.

### False final partitions under a one-population null

The endpoint is P(final K > 1), with prespecified target 0.01. It is distinct
from nodewise p-value validity and final-method FDR. Each null has 400 cells.

| Variant | Low depth variation: false partitions / datasets | High depth variation: false partitions / datasets |
| --- | --- | --- |
| symmetric_average | 0/300 | 262/300 |
| symmetric_ward | 295/300 | 228/300 |
| markov_average | 34/300 | 298/300 |
| markov_ward | 297/300 | 300/300 |

One-sided 95% exact binomial bounds are in `validation_null_summary.csv`; each
bound is marginal, with no simultaneous calibration claim. Symmetric/average
on the low-variation null has upper bound 0.009936. Its high-variation null has
262/300 false partitions and separate one-sided bounds 0.837299 and 0.903722.
Every variant's high-variation lower bound exceeds 0.01. No variant meets the
target across both models. No simulation outcome was silently dropped or
replaced with a conservative partition.

### Signal recovery

Mean ARI by independently generated signal model:

| Signal model | Symmetric / average | Symmetric / Ward | Markov / average | Markov / Ward |
| --- | --- | --- | --- | --- |
| signal_balanced | 0.700000 | 0.988412 | 0.740493 | 0.305492 |
| signal_imbalanced | 1.000000 | 0.917122 | 0.809758 | 0.074621 |
| signal_singleton | 0.879876 | 0.917483 | 0.877034 | 0.085487 |
| signal_density_variation | 0.582483 | 0.106331 | 0.069223 | 0.023421 |

Symmetric/Ward improves mean ARI on the balanced mixtures, while its null
failures prevent adoption. The Markov variants often overpartition the signal
mixtures, especially the density-variation cases. Per-dataset cluster counts,
purity, dominant-cluster recall, split/merge errors and solver diagnostics are
retained in `validation_dataset_results.csv` and the primary JSONL capture.
The singleton model's ARI is not an isolated-singleton recovery guarantee.

## Verification, artifacts and limits

- The 45 relevant existing eigenspace/tree/branch-length tests passed before
  the experiment. Assert-based checks verified the Markov relation, preservation
  of input coordinates, zero-degree rejection, seed reproducibility and an
  end-to-end 400-cell count simulation. All four experiment scripts compile.
- Aggregation asserts exactly 635 independent datasets and 2,540 unique
  dataset/variant records, complete frozen seed/model coverage and matching
  paired fingerprints. No simulation exceptions or nonconverged fits occurred.
- Diffusion, native TBS runner and method-registry hashes remained unchanged
  throughout the experiment. Execution-only sharding and failure capture are
  recorded separately; no numerical setting was changed.
- The nulls have 400 cells; signal datasets have 1,000, 1,001 or 3,000. They do
  not certify behavior on arbitrary distributions or 48,879-cell nulls. No
  selected-hierarchy p-value theorem or descendant-FDR guarantee was established.
- The three alternatives have no valid full-data clustering export. The
  original one-cluster AnnData remains available. Increasing the solver budget
  alone would leave the observed false-partition failures unresolved.
- Production defaults and dependencies were not changed. Temporary pairwise
  distance files were removed and experiment workers ended. No commit or push.

Local outputs are in `results/analyses/choir_tbs_recovery_comparison_20261103_174622/`.
`compare.py` runs real-data alternatives; `validate.py` runs frozen simulations;
`aggregate.py` checks completeness and builds tables. Original local PCA/kernel
checkpoints and bulk expression files are required; they are absent from fresh
clones. Primary compact captures and script snapshots are in
`raw/inbox/choir-cell-lines/tbs-recovery-*`.
