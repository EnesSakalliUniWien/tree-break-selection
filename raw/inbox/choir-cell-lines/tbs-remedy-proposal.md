# Proposed recovery experiment for cancer cell-line TBS

Proposed on 2026-11-03 after the user asked for a concrete solution. This is a
testable proposal, not an implemented or validated repair.

## First experiment: isolate coordinates from linkage

Reuse the existing full-data PCA, kernel and retained eigenpairs. Compare:

| Coordinates | Linkage | Status |
| --- | --- | --- |
| Current symmetric coordinates | Average | Existing one-cluster baseline |
| Current symmetric coordinates | Ward | Proposed control |
| Markov right-eigenvector coordinates | Average | Proposed coordinate control |
| Markov right-eigenvector coordinates | Ward | Proposed combined candidate |

Keep all 48,879 cells, the 2,000 variable genes, 30 PCs, kernel parameters,
diffusion time 3, complete eigenvalue groups, NNLS target/seed/sample/tolerances,
covariance policy, edge alpha 0.001, sibling alpha 0.01, Tree-BH and traversal
fixed. Fit NNLS separately on each resulting topology using the same original
PCA target. Do not copy fitted edge lengths across different trees. Do not
remove cells, force 190/367 clusters or tune thresholds using the published labels.

For the saved symmetric kernel W and positive degree matrix D, the current
coordinates are Phi_S = U Lambda^3 for S = D^(-1/2) W D^(-1/2). The candidate
Phi_M = D^(-1/2) Phi_S uses right eigenvectors of P = D^(-1) W. This changes the
coordinate geometry; it must be a named variant rather than a silent change
to the current symmetric-coordinate definition. Retain all three leading
component modes and the existing complete-tie rule for this comparison.

A read-only calculation verifies the candidate eigenvector relation with
maximum residual 4.51e-15. Within the 290 ASPC1 cells, coordinate-norm coefficient
of variation (population standard deviation / mean) changes from 0.190099 to
0.050017. This motivates the comparison but does not establish improved tree
structure, cluster recovery, biological interpretation or error control.

Ward linkage is already supported by the project and SciPy. It minimizes
within-cluster variance and requires Euclidean pairwise distances. Compute
those distances from each explicit coordinate matrix. It does not guarantee
balanced splits or successful TBS gates. No new dependency is needed.

## Evaluation and adoption

Record root and upper-tree child sizes, covariance fallbacks, reached versus
ancestor-blocked tests, gate decisions, final K, ARI/NMI, cell-line fragmentation
and mixed-cluster composition. Use the known cell-line and published CHOIR labels
only for diagnostic scoring, with no target-K constraint. A higher cancer-data
ARI after selecting a variant is a development result, not unbiased validation.

Before adopting a candidate, freeze it and evaluate the complete procedure on
independent null/signal datasets and seeds, including large imbalanced groups,
singletons and density variation. Define and measure false final partitions
under a homogeneous null as well as power on mixtures. State the null model
and report uncertainty; a fitted or mismatched simulation is not an exact
reference. Rebuild preprocessing, diffusion, topology, NNLS and testing wherever
the evaluation's selection mechanism requires it. Fixed-tree label permutations
or an ordinary row split do not automatically remove selection bias.

If these geometry controls still stop at singleton roots, a later methodological
change must separate subtree exploration from accepting an immediate split.
That needs a defined hypothesis family and validated selection-aware calibration.
Do not simply bypass Tree-BH, relax alpha or label additional clusters as a fix.
The existing same-data p-values already lack a general selected-hierarchy error
guarantee; unchanged gates do not confer one on a new tree.

## Sources

- `tbs-diagnosis.md`, `tbs-root-variance-diagnostic.json`, and
  `tbs-structure-diagnostic.json`: observed root failure and frozen-tree control.
- `tbs-remedy-proposal-evidence.json`: coordinate calculation, frozen settings
  and source metadata; no new clustering or production changes.
- `diffusion-maps-coifman-lafon-2006.pdf`, Appendix A: Markov/symmetric spectral
  relation. Original: https://doi.org/10.1016/j.acha.2006.04.006.
- SciPy linkage documentation:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html.
- `tree_break_selection/tree/construction/defaults.py` and
  `tree_break_selection/space_separation/README.md`: current method support and
  the distinction between symmetric and Markov coordinates.
- `reports/calibration_restoration_20260909/validation_protocol.md`: full-pipeline
  selection, held-out evaluation and uncertainty requirements.
