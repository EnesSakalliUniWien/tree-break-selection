# Actual diffusion selection event: answers to Q1 and Q2

Date: 2026-11-02. Scope: research and deterministic diagnostic replay; no production changes, threshold changes or selected p-value implementation.

## Q1 — What preserves the observed complete merge sequence?

### Answer

Fix sample labels/order, software configuration, and the complete observed merge sequence M. Let L,R be its root children, eta_i=1/n_L on L and -1/n_R on R, e=eta/||eta||. For X with independent Gaussian rows, known common feature covariance sigma² I_p and null eta^T E[X]=0, write

z=X^T e/sigma, r_obs=||z||, u=z/||z||,
A=(I-ee^T)X,
X(r)=A+sigma e r u^T, r>=0.

The sigma=1 convention is X(r)=A+e r u^T. Condition on A and u, treating the observed sequence and hence e as fixed within its selection event. For a fixed e, Gaussian orthogonal decomposition makes X^T e independent of A; under the contrast null z~N_p(0,I), its radius is independent of its direction and has chi_p law. Selecting e by the observed sequence is accounted for by conditioning on that sequence; it is not an independent-topology assertion.

Let Phi(r) be the ACTUAL refitted feature_family_diffusion_geometry(X(r), continuous_feature_space), delta_ij(r)=||Phi_i(r)-Phi_j(r)||_2, and

D_AB(r)=(|A||B|)^(-1) sum_{i in A,j in B} delta_ij(r).

At observed merge step t, use the observed active clusters C_t and selected pair (A_t,B_t). The exact full-sequence event is

S_M={r>=0: the geometry procedure is valid at X(r), and for every t and every competing active pair (C,D), the algorithm selects (A_t,B_t) ahead of (C,D)}.

Away from merge-score ties this is precisely the intersection of D_AtBt(r)-D_CD(r)<=0. At ties the actual configured implementation's deterministic ordering must be included; the repository delegates linkage tie behavior to SciPy rather than declaring a symbolic global tie order. This event is sufficient for the observed root split and stronger than conditioning only on the root partition or the unordered tree clades. The final merge has no competitor. Active clusters are fixed by the observed history when writing these constraints, rather than refitted arbitrarily in each inequality.

For positive conditional radial mass, the exact model-level selected tail is

p_sel = integral_{S_M intersect [r_obs,infinity)} f_chi_p(r) dr / integral_{S_M} f_chi_p(r) dr.

This law is valid independently of how complicated S_M is. It does NOT supply its boundaries, certify a grid, justify an estimated covariance, or establish a production adaptive-PCA/sibling/traversal law. Here p is the tested raw full-feature dimension, not the retained diffusion-coordinate count. The global iid common-mean Gaussian null suffices; known sigma and isotropic feature covariance are assumptions, not estimated production quantities.

There is no reason to replace these inequalities with the squared-Euclidean quadratic hierarchy engine: actual linkage averages EUCLIDEAN diffusion-coordinate distances. Its diffusion map is refitted at each r. Original pairwise squared distances are quadratic in r, but the refitted final diffusion distances are not.

Conditioning additionally on discrete intermediate-state cells (neighbor membership/order, BGH winning index, rounded dimension, retained spectral groups) is stronger than S_M. This can still define a positive-mass truncated chi law if the extra selection event is explicitly included. Freezing exact continuous bandwidths, epsilon, kernel entries or eigenvectors is substantially different: these values generally constrain or determine r, potentially leaving zero radial mass. They cannot silently be treated as independent nuisance values. Holding a graph fixed without conditioning on its observed support is also a changed replay procedure.

### Evidence: actual implementation

Primary inspected sources:

- `tree_break_selection/space_separation/diffusion.py`: `feature_family_diffusion_geometry`, `adaptive_diffusion_geometry`, `resolve_adaptive_epsilon`, `compute_diffusion_coordinates`, `_leading_diffusion_spectrum`, `_sparse_leading_eigenpairs`.
- `tree_break_selection/tree/construction/build.py`: canonical sample order, average-linkage invocation and input distance reordering.
- `tree_break_selection/eigenvalue_groups.py`: numerical spectral groups.
- `.venv/lib/python3.12/site-packages/pydiffmap/kernel.py`: installed `Kernel`, `NNKDE`, `choose_optimal_epsilon_BGH`, `_scale_by_bw`.
- `wiki/analyses/root-selected-region-model.md`: earlier conditional radial construction and its diffusion compatibility boundary.

Continuous preset: Euclidean input metric, requested k=15, diffusion time=3, target coordinate count min(50,n-2,p), bandwidth_type='-1/(d+2)', outer epsilon='median'. Duplicate-aware neighbor guard applies. In the nine-row fixture the resolved k is8 and target count1.

Kernel fitting uses scikit-learn NearestNeighbors. NNKDE(k=8) constructs local pilot bandwidths from seven nearest positive squared distances. Its internal BGH calculation scans 81 epsilon candidates 2^-40,...,2^40, takes the first argmax of adjacent log-kernel-sum slopes, and rounds twice that maximum slope to an estimated intrinsic dimension. KDE recomputes local query bandwidths/densities; the final bandwidth is rho_i=q_i^(-1/(d+2)). Thus specifying outer epsilon='median' does NOT remove BGH adaptation: its internal epsilon and rounded dimension still control KDE bandwidths.

The kernel's scaled distance is Euclidean distance divided by sqrt(rho_i rho_j); outer epsilon is the median of its positive finite squared entries. The directed sparse Gaussian weights are exp(-scaled_distance²/(4 epsilon)). The repository uses pydiffmap Kernel, not pydiffmap DiffusionMap's density-alpha normalization. It zeroes the diagonal, takes elementwise maximum with the transpose, computes degrees, and diagonalizes D^-1/2 W D^-1/2. For n<=2000 it uses dense scipy.linalg.eigh; larger sparse operators have a deterministic-start ARPACK route, with dense fallback for disconnected support or unresolved retained groups. It drops the leading mode only when isolated, extends the retained cutoff through numerical ties, clips negative retained eigenvalues to0, and multiplies vectors by eigenvalue^3. Diffusion distances are ordinary Euclidean distances of these coordinates.

Neighbor support/ranking, BGH argmax, rounded dimension, epsilon order statistics, spectral groups and linkage ordering can all change along r. Eigenvector signs do not change distances; rotations inside an exactly repeated group retained whole also preserve distances. The numerical group tolerance is adjacent np.isclose(rtol=1e-10, atol=1e-12 times full spectral scale), with transitive grouping. This is not simply a fixed-rank polynomial map.

### Unresolved

No exact boundary-enumeration or certified integration routine exists for this complete diffusion event. Interior replay and a finite grid do not exclude missed narrow components. A tractable conditional model must preserve the actual refitting operations or explicitly state a stronger conditioning event with its own positive-mass law. Floating-point library tie behavior and the gap between an exact arithmetic map and numerical implementations remain relevant to a certification claim.

## Q2 — Are multiple boundaries/disconnected admissible sets possible, and how do n,p,rank,geometry matter?

### Answer

Yes for the COMPLETE-MERGE-SEQUENCE event of the actual continuous feature-family diffusion procedure. A deterministic n=9,p=1 standard-normal fixture exhibits return to the observed sequence after an intervening nonmatching sequence. This occurs with retained diffusion dimension1. At the four localized transition brackets, neighbor lists/order, internal BGH epsilon, rounded intrinsic dimension and retained rank agree across the bracket endpoints. Therefore graph or rank switches are not necessary for multiple selection boundaries; diffusion merge-score crossings can suffice.

This fixture does NOT prove a disconnected ROOT-PARTITION-ONLY event: the root partition at normalized radii1,1.1,1.2 is the same. The gap swaps the temporal order of two disjoint intermediate merges. Root-only selection unions over compatible histories and can bridge the gap.

Neither n nor p nor retained rank alone determines interval count. The number of observed candidate-pair comparisons is sum_{k=2}^n (choose(k,2)-1), but each comparison is a nonlinear diffusion-score function, and this combinatorial count is not a bound on its number of radial roots. Larger n adds comparisons, potential graph order changes and spectral modes; p changes original distance coefficients and the preset target dimension; retained dimension and spectral gaps change geometry. None establishes monotonicity or a general interval-count formula. Even fixed n,p,rank can yield different observed radial event geometry. No monotone empirical dependence is established by this small search.

### Evidence: reproducible bounded probe

Artifacts:

- `raw/assets/diffusion-selection-event-20261102/reentry_probe.py`: independent executable probe with reconstruction/nuisance assertions and actual-map replay.
- `raw/assets/diffusion-selection-event-20261102/reentry.json`: exact input matrix, contrast/path/direction, root leaf sets, observed/grid signatures, nearest-neighbor lists, epsilon/dimension metadata, transition brackets, versions and source hashes.
- `raw/assets/diffusion-selection-event-20261102/initial_panel.json`: preceding 17-case finite panel, including twelve n12,p2 cases, four n16,p4 cases, and the n80,p24 attribution seed2026110203.

Run:

```bash
uv run --no-sync python raw/assets/diffusion-selection-event-20261102/reentry_probe.py
```

Fixture uses numpy.default_rng(20261105).normal(size=(9,1)), sigma1. Observed standardized radius is2.4990448450502516; root leaf sets are {3,4,5} and {0,1,2,6,7,8}. The observed complete signature (zero-based SciPy IDs) is [(2,6),(3,5),(7,8),(4,10),(1,9),(0,13),(11,14),(12,15)]. At radius ratio1.1 it is [(2,6),(3,5),(7,8),(1,9),(4,10),(0,12),(11,14),(13,15)], and at1.2 the observed signature returns.

A grid over [0,4] (121 points) plus20 logarithmic points from4.1 to100 locates two sampled matching runs. Thirty bisection iterations per observed transition yield the following normalized-radius brackets; endpoints have opposite membership:

| Transition | Lower ratio | Upper ratio |
|---|---:|---:|
| enter |0.9479473985421161|0.9479473985731601|
| exit |1.070727343143275|1.070727343174319|
| reenter |1.1195953217335042|1.1195953217645485|
| exit |1.3991668048004307|1.399166804831475|

At every bracket, both endpoint replays use NNKDE dimension1, internal epsilon0.12500000000000003, identical neighbor lists/order and retained diffusion dimension1. Outer median epsilon changes continuously across these tiny brackets. This identifies localized changes in merge selection without those discrete state switches; it is not a global continuity proof or unique-root certificate.

The initial17-case scan used81 equally spaced ratios in[0,3] plus the observed ratio1 and found only one sampled matching run per case. The subsequent n9,p1 bounded search found reentry on its fourth seed (20261105) and stopped. Failure to find another component in other fixtures does not establish connectedness; no frequency estimate or dimensional trend is claimed.

### Unresolved

There may be further components between samples, beyond ratio100, or inside a bisection bracket. The output brackets are localized numerical changes, not exact algebraic endpoints. No selected-tail integral is computed from this grid. The example concerns complete history, not necessarily root partition, and does not show that n9,p1 behavior generalizes to the n80,p24 experiment. Establishing exact component counts or a dependence law on n,p,rank and diffusion settings requires a certified treatment of all piecewise kernel/spectral/merge changes.

## Verification performed

The reentry probe ran successfully in the existing uv --no-sync environment, reconstructed X at its observed radius, verified e^T A=0, and asserted matching complete signatures at radius ratios1 and1.2 with a nonmatching replay at1.1. It records package versions/source hashes. No tests of production calibration were added or claimed. No production, wiki, dependency, threshold, commit or service change was made by this investigation.

An additional nine independent replays at ratios1±1e-5,1.1±1e-5,1.2±1e-5 (including the centers) reproduced matching neighborhoods around1 and1.2 separated by a nonmatching neighborhood around1.1. These local checks strengthen reproducibility without certifying every point in an interval.
