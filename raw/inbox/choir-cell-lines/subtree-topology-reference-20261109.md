# Subtree distances, alternative topologies and a Gaussian reference

9 November 2026, Europe/Vienna. This continues the root validation framework and assesses the user's claim that distances, equidistances and conflicting distances across eigenspaces identify alternative relationships and constrain their possible topologies. Three Sol agents checked the mathematics, geometry and calibration/bootstrapping; the coordinator connected the diagnostic to native sibling statistics and reran the integrated tests. This is an opt-in diagnostic. Production code and defaults were preserved.

## What the geometry establishes

For a fixed subtree S with m labeled leaves, a distance matrix, a reconstruction rule and a rule for equidistances define a family of possible outputs. Denote it by T_b(D_b; A, tau), where b identifies a subspace, A is the reconstruction rule and tau is its declared tie tolerance. The implemented family uses Ward's increase in within-cluster squared error at each merge and retains every permitted minimizing merge. Different merge orders that yield the same rooted topology are counted once. The count is conditional on these choices, not a universal number inferred from eigenspace dimension alone.

For orthogonal spectral projectors P_b, projected squared distances are `(x_i-x_j)' P_b (x_i-x_j)`. Within a repeated eigenvalue block the whole projector is invariant to an arbitrary change of eigenvector basis. An individual eigenvector inside that block is not a uniquely defined source of relationships. With specified nonnegative metric weights, the total squared distance is the weighted sum of block squared distances. Disagreement between blocks can therefore reflect separate measured structure. It does not require that one block be noise.

The union of the families counts topologies admitted by at least one space; their intersection identifies common topologies. These questions differ from counting all mathematically possible trees or all bipartitions. On m leaves there are `2^(m-1)-1` unordered nonempty cuts and `(2m-3)!!` rooted labeled binary topologies. Rank one still permits every rooted topology across suitably spaced point clouds. For one fixed general-position line cloud there are only m-1 strict affine cuts. Whole-tree topology, an edge split and a flat cut are distinct objects.

| Four-leaf geometry | Ward topology count | Interpretation |
| --- | ---: | --- |
| Generic distinct line positions | 1 | A unique sequence up to commuting independent merges in this example |
| First rectangle axis | 1 | AB\|CD |
| Second rectangle axis | 1 | AC\|BD; union with the first has 2 trees, intersection 0 |
| Square | 2 | Four tied first-merge histories reduce to two rooted trees |
| Equidistant simplex | 15 | Every rooted binary topology; 18 ranked merge histories |
| Near-square in the recorded example | 1 at zero tolerance, 2 at 1e-10 | Near-tie family depends on the declared tolerance |

The rectangle axes reverse four of the 15 orders between pair distances; eight other comparisons are tied in only one axis. The two rooted trees have Robinson–Foulds distance 4. The full rectangle and first-axis Ward trees share a topology while their cophenetic distances differ, with RMSE 1.155. Thus topology equality does not imply equal branch lengths or equal distance fit.

Given each candidate topology, the diagnostic fits nonnegative additive edge lengths to the same observed leaf distances. The residuals distinguish a shared branching structure from its geometric fit:

| Input distances | RMSE for AB\|CD | RMSE for AC\|BD |
| --- | ---: | ---: |
| First rectangle axis | 0 | 1.8856 |
| Second rectangle axis | .9428 | 0 to numerical precision |
| Full rectangle | .1927 | 1.0716 |
| Full square | .3382 | .3382 |

All 15 rooted simplex topologies fit to numerical precision, but their internal fitted edges collapse to zero: the distances support an unresolved star, not 15 distinguishable resolved trees. Rerootings of an identical unrooted topology have identical additive fits. The fitter therefore combines complementary root edges instead of claiming that pair distances identify their separate lengths or the root position. Five input matrices and all 75 fitted matrices satisfy the triangle checks. These RMSEs concern fitting unrooted additive distances; they are distinct from Ward's merge criterion and from the cophenetic comparison above.

The counts are exhaustive only inside the implemented bounds: all cuts through ten leaves and full rooted topology/tied-Ward enumeration through six. Relative tolerance is part of the defined numerical problem. A sequence of allowed near-tie choices need not correspond to any one perturbation of the entire distance matrix within a statistical uncertainty region.

## When a topology can be called correct

Distances can identify a topology under structural assumptions. Exact additive distances with positive internal edges identify the reduced unrooted tree. Exact ultrametric distances identify a rooted hierarchy; simultaneous merge heights leave binary refinements unresolved. In a quartet, the additive four-point sums have two equal largest values. The remaining gap is twice the internal edge length. An error bound on distances can protect that ordering when it is smaller than the gap allows. Euclidean projected distances need not be additive or ultrametric. A small fitting residual is evidence of compatibility with a chosen tree model, not proof of a biological population tree. Conflicting spaces may encode crossing groupings that one tree cannot represent. The mathematical report gives the precise distinction and primary references. [Buneman](https://doi.org/10.1016/0095-8956(74)90047-1), [Carlsson and Mémoli](https://jmlr.csail.mit.edu/papers/volume11/carlsson10a/carlsson10a.pdf)

## Gaussian mathematical construction and implemented reference

On a prespecified subtree of m independent Gaussian rows with one common mean and arbitrary unknown shared covariance, condition on its mean and centered feature Gram matrix. For a fixed normalized label contrast u and a rank-k subspace chosen from that Gram only,

`R² = u' Proj(col(YV)) u ~ Beta(k/2, (m-1-k)/2)`, for `0 < k < m-1`.

For the unbiased pooled within covariance ratio rho, `R²=m rho/(m-2+m rho)`. The all-label Bonferroni factor is `2^(m-1)-1`; its critical rho tends to 3 with fixed k as m increases. This is a conservative Gaussian all-label threshold limit, not a general population unimodality bound. Replacing all labels by the observed candidate count is unjustified with the same fixed-label tail. Contrast-selected rank does not satisfy the Gram-only hypothesis. Full centered row rank makes R² identically one, so that Beta reference degenerates.

The implemented conditional reference uses a different argument for arbitrary scores: rotate centered subtree rows on their conditional Gaussian orbit, preserve the mean and Gram, freeze outside rows, and reject proposals outside a declared selection event. A symmetric orthogonal proposal includes reflections. The observed state is placed at an independently uniform position in a fixed-length reversible path, and the p-value is its weak upper rank over all path positions. Numerical failures return unresolved results. Canonical row-ID order fixes the proposal basis, including when input indices arrived in feature-sorted order. The mathematical rank argument is exact in the model; floating-point event membership is not formally certified. [Chen and Gao](https://arxiv.org/html/2610.17458v1), [Besag and Clifford](https://doi.org/10.1093/biomet/76.4.633)

The null requires every row in the tested subtree to have the same mean. It is stronger than equality of two aggregate means while individual rows remain heterogeneous. It does not assume known covariance. Selecting a subtree requires replaying its selection; the raw Beta helper does not supply that conditioning.

## All three testing objectives

1. Individual cuts: the fixed-subtree Beta helper and conditional rank are implemented. A selected native sibling is tested by preserving a prespecified Ward node path and its observed split, while refitting the native raw sibling score on every accepted orientation.
2. Topology comparisons: rooted RF, edge-split comparisons, cophenetic distances and additive distance residuals describe the relationship between trees and each space. A paired conditional reference replays the absolute difference between two prespecified summed sibling scores on the same orientation. Its null is common Gaussian homogeneity, not equal reconstruction accuracy of arbitrary methods.
3. Joint search: one reference regenerates all strict affine cuts in the full spectral span and their maximum raw full-space sibling score. Another regenerates the union of the tied-Ward families across whole spectral blocks and cumulative spans, then takes the maximum summed sibling score. Every reported candidate receives a tail against the same maximum path. These give simultaneous control under the complete common-mean subtree null, not strong FWER across heterogeneous subtrees or final-partition FDR.

Including the full cumulative span makes the union of all smaller-span affine cuts equal to the full-span cut family. It does not make the corresponding Ward topology families identical. The topology-search score is an explicitly declared diagnostic sum, not a likelihood or an estimator of the probability that a topology is correct.

The native adapter refits Ward distances, fixed-topology NNLS branch lengths, the edge spectral context, guarded covariance and the .90 adaptive contrast rank. It retains `record.stat`, the raw sibling statistic. The displayed `Sibling_Test_Statistic` has already been divided by empirical inflation. Neither the empirical-inflation p-value nor eligibility and traversal are promoted to an exact reference. The alternative-cut control uses the established full-space raw sibling statistic with zero branch lengths and adaptive rank disabled, and is labeled separately.

## Rerun evidence and power limits

The core reference panels used fixed protocols and independent data/sampler seeds. When selecting the best of four Gaussian candidate subtrees, the naive chain rejected 219/1,200 at nominal .05; conditioning on that selection reduced this to 60/1,200, with exact 95% interval [.0384,.0639]. The exact selected truncated-Beta reference rejected 59/1,200. The full all-label maximum rejected 13/320 on fixed eight-row subtrees, and a median-selected descendant rejected 8/160 after its membership event was replayed. These last intervals remain broad. The full-rank orbit supports an arbitrary-score reference while the full-rank Beta helper explicitly reports its degeneracy.

The integrated experiment retained all 138 attempts. At nominal .05, cut-search nulls rejected 0/50 (exact interval [0,.0711]), topology-search nulls 1/40 ([.00063,.1316]), and prespecified topology pairs 2/40 ([.00611,.1692]). No attempt was unresolved or stalled. All 530 cut-adjusted and 86 topology-adjusted tails were recomputed independently and were at least their joint tail. All six native observed scores refitted exactly, and a complete native null-root chain reproduced exactly.

Power is unresolved for broad use. Both planted native roots had the correct 8+8 split but did not reject with 99 transitions (p=.99,.07). A disclosed exploratory follow-up retained both inputs: 999 transitions gave p=.504,.140; doubling the planted shift gave .034,.017. The two small six-row cut-search examples also did not reject (.09,.48). These results establish working integration and show sensitivity to sample size, selection and signal strength; they do not establish general effectiveness.

## What bootstrapping approximates

The experiment preserved leaf identities and compared 1,000 feature resamples, 1,000 small Gaussian perturbations, 1,000 larger perturbations and 1,000 row permutations for each of three configurations: 12,000 replicates. On the square, the exact geometry family had two topologies. Low measurement noise stayed within those two in all 1,000 replicates. Feature resampling produced seven topologies, with 109 draws outside the fixed family. Larger measurement noise produced 14, with 280 outside. The planted six-leaf hierarchy was recovered in 934/1,000 feature draws, 1,000/1,000 small-noise draws and 977/1,000 larger-noise draws.

Thus geometry can describe the alternatives around which a specified resampling law moves. It does not alone determine their probabilities, and finite perturbations can reach other topologies. Feature resampling assumes appropriate resampling units; this demonstration deliberately includes dependent coordinates, so those frequencies are descriptive. Canonical row permutation preserved labeled distances and the tree in every control; fitting in shuffled order merely switched between the square's two existing ties. Null permutations, conditional Gaussian rotations and bootstrap perturbations answer different questions. None is automatically the probability law of the geometric count. [Felsenstein](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x)

## Evidence and reproducibility

- `subtree-topology-mathematics-20261109.md` and its Wolfram JSON: five symbolic checks, proofs, counts, assumptions and primary references.
- `subtree-topology-geometry-20261109.md`: exact bounded families, distance conflicts, tie conventions and topology-conditioned fits.
- `subtree-topology-calibration-20261109.md`: 2,900 primary Gaussian datasets and 480 all-label/descendant datasets, including the preserved canonicalization correction.
- `subtree-topology-bootstrap-20261109.md`: all bootstrap frequencies, RF distributions and assumptions.
- `subtree-topology-reference-20261109-evidence.json`: combined summaries, verification and source provenance.
- `subtree-topology-reference-20261109-sources.py.txt`: runnable scripts and frozen protocols; maintained implementations are under `benchmarks/diagnostics/calibration/selected/hierarchy/subtree_*.py`.
- Local detailed runs and PNG/PDF/SVG figure: `results/analyses/subtree_topology_reference_20261109/`.

Final verification: 123 focused reference, topology, native-score, eigenspace and existing root/descendant tests passed. Ruff checks and formatting passed on all six new maintained Python files. Wiki lint passed on 264 pages. Production source hashes remained unchanged during integration and at completion; recorded geometry additions preserve the earlier definitions and their provenance. The new diagnostics leave production defaults untouched. No dependencies, commits, pushes or services were introduced. Numerical near-ties, larger-subtree computation, non-Gaussian preprocessing, general native-score power, full TBS selection and final-partition control remain open.
