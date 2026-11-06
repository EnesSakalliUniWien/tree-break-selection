# Tree-Break Selection: statistical construction audit

2026-11-04. Three Sol agents independently audited Marchenko–Pastur and rank,
raw reference laws, and accumulated method knowledge. The coordinator checked
calibration, multiplicity and traversal, integrated the findings, and corrected
four current wiki explanations. Production inference and historical evidence
were preserved. This is a source-level and controlled numerical diagnosis,
not a validated replacement statistic.

## Main conclusion

**The stages do not share one established null law.** Marchenko–Pastur (MP)
provides a nominal spectral heuristic. Estimated covariance, estimated PCA,
contrast-dependent rank, fitted branch lengths, selected membership and
empirical calibration can each alter the distribution tested afterward.
Reporting the final projection rank correctly does not establish that the
quadratic form has a chi-square law with that many degrees of freedom.

The audit separates four quantities that previous descriptions conflated:

1. Raw MP count: eigenvalues strictly above the unit-spectrum upper edge.
2. Nominal projection dimension: floor, positive-rank cap and tie completion;
   siblings additionally derive a nominal dimension from child contexts.
3. Retained projection dimension: the actual rows used after optional
   contrast-energy selection and tie completion.
4. Reference degrees of freedom: the kernel uses the retained row count as a
   chi-square parameter. That parameter is justified only under additional
   assumptions about the standardized contrast conditional on the projection.

## Step-by-step verdicts

| Stage | What is established | Defect, assumption mismatch or open condition |
| --- | --- | --- |
| Feature representation and hypothesis | Node summaries are leaf-count weighted means/frequencies. | Equal means do not imply homogeneous descendants. Bernoulli/categorical tests need their discrete or asymptotic laws; separate blocks omit cross-block covariance. |
| Tree construction and row ownership | Edge means use alternating test rows within already selected children; sibling means use all leaves. | Fitting the tree used both halves. A split after selection is not independent discovery and testing. Diffusion and full-scale NNLS were not rerun in this audit. |
| Nested/disjoint contrast variance | The familiar `1/n_child−1/n_parent` and `1/n_left+1/n_right` factors are correct for fixed groups, independent rows and common covariance. | Heterogeneous row covariances require the weighted covariance sum. Selected membership changes the contrast law. |
| Continuous covariance | Guarded pooling and fallback are implemented; counts, residual df, storage width and memory limits choose full versus diagonal matrices. | Parent-total and pooled covariance have different dependence and finite-sample laws. The edge denominator still uses full-data covariance. |
| Whitening and ridge | Cholesky/diagonal solves implement the stated block model. | Fitted covariance does not produce an exactly standard Gaussian vector. Ridge and omitted correlations can give a weighted chi-square even under Gaussian data. |
| Branch-variance factor | The code applies its configured scalar factors. | Fitted NNLS or linkage lengths do not by themselves justify a sampling-variance model. The edge switch and projected sibling behavior differ. |
| Spectral matrix and MP edge | Primal and dual positive spectra share the actual row divisor. Unit MP edge is an asymptotic reference. | Supplied leaf counts are not measured independent sample sizes. Centering, same-row whitening, optional internal rows and tree selection change the ensemble. |
| PCA and rank | Available directions, numerical ties and adaptive rank are traced. | In energy mode a positive nominal MP request is replaced using the tested contrast. No invented orthonormal complement is appended. |
| Raw statistic and df | `Q = norm(P z)^2`, scale one, df equal to retained rows. | Correct rank bookkeeping does not establish chi-square calibration. Covariance, basis and rank need compatible conditional laws. |
| Empirical-null calibration | The weighted scale objective, local kernel, family fallback and support diagnostics are reproducible. | Stopped-edge roles can contain signal. Dependency-group ESS does not turn the record-level scale fit into an independent chi-square denominator. |
| Multiple testing | Tree-BH family scaling and per-depth sibling BH arithmetic pass focused checks. | Local p-value validity, dependence and hierarchy assumptions are unestablished. The mean-contrast hierarchy is not generally an intersection hierarchy. |
| Traversal | Binary structure, an open edge, sibling flags and optional supported pass-through produce deterministic boundaries. | An edge-reachable sibling frontier differs from live traversal. A closed edge blocks descendants; pass-through may cross a closed sibling gate. A final-partition error theorem is separate. |

## Marchenko–Pastur findings

At energy fraction 0.90, the same 30-direction basis returned rank 1 for a
contrast along its first direction and rank 30 for a contrast along its last,
for each supplied nominal rank 1, 2 and 30. Thus the MP count is not the final
rank controller in this mode. Rank zero exits before adaptive selection.

For the same fixed 60-by-30 Gaussian matrix, fitted full covariance whitening
made every eigenvalue approximately `59/60`. The raw MP count was zero, but
complete numerical tie retention expanded the nominal rank to 30. Fixed known
identity and fitted diagonal whitening each yielded nominal rank two. This is
a change in the spectral ensemble, not evidence that 30 signal dimensions were
detected. The current two-row fitted-diagonal case has rank one, eigenvalue
approximately 15 and zero MP spikes.

With fixed known whitening, centered Gaussian scatter has matrix law
`Wishart_d(m−1,I)/m`. The implemented `(1+sqrt(d/m))²` edge and the centered
scale-matched edge are different asymptotic conventions. Neither is a chosen
finite-sample significance quantile. The primary literature and 3,000-draw
controls for each of three matrix shapes are documented in the
[MP audit](mp_audit.md) and its [source notes](../../raw/inbox/tbs-mp-primary-sources-20261104.md).

## Degrees of freedom and reference-law controls

All counts below use 20,000 fixed-group Gaussian draws and nominal alpha 0.05.
These are separate mechanism controls, not a comparison of complete clustering
methods. Exact laws omit ridge and extra branch scaling and require the stated
independence and projection conditions.

| Control | Nominal chi-square rejections | Appropriate restricted-reference rejections |
| --- | ---: | ---: |
| Known covariance, fixed rank-three projection | 989 | 989 |
| Same-row parent-total covariance, fixed rank three, total residual df 11 | 305 | 972 with scaled Beta |
| Full means, default-guard-compatible pooled covariance, 8+8 rows and 6 features | 6,540 | 1,027 with Hotelling F |
| Known covariance, fixed eight-axis basis, contrast-energy fraction 0.90 | 1,562 | No replacement law implemented |
| Known covariance, rank-one direction fitted to the tested 30-feature pair | 20,000 | 1,001 with chi-square(30) |

The parent-total control has `Q/nu ~ Beta(k/2,(nu−k)/2)` under the report's
conditional-Haar assumptions. The full-mean pooled-covariance control has
`Q ~ nu*p/(nu−p+1) * F(p,nu−p+1)` under the independent Wishart model. These
are different distributions, not different guesses for chi-square df.

The edge's alternating-row mean shares residual variation with its full-data
within-child denominator. In the 8+8 control, using the full-mean F law for
that edge gives 241 rejections, versus 4,060 under nominal chi-square. Neither
count establishes a calibrated edge replacement. Native covariance resolution
and both test wrappers matched the analytic quadratic forms in 25 datasets to
maximum absolute error `7.42e−10`. The independent Sol review caught an earlier
6+6 mathematical control below the production guard; it remains explicitly
labeled, and the 8+8 native replay resolves that implementation-comparison gap.

See the [reference-law audit](reference_law_audit.md) for the derivations,
discrete-family conditions, adaptive-rank mixture and covariance dependencies.

## Calibration, multiplicity and traversal checks

The saved cancer surrogate shows 11 mean-different records among 494 supported
calibration records at a failed root. The independent downstream probe also
shows that adding records within one dependency group can change the fitted
family scale from 5.5 to 1.81818 while group ESS remains two. This is the
implemented record-level objective, not a new bug in its arithmetic. Support
diagnostics do not make those records independent or null.

The Tree-BH ancestor-proportion recurrence matches the local implementation.
Transferring the published error guarantees requires their additional
conditions. In particular, an ancestor mean contrast can be zero while both
descendant contrasts are nonzero; its null is not the intersection of those
descendant nulls. Parent-edge p-values are also directly computed rather than
formed by Simes combination. This identifies missing premises, not a proof
that no alternative guarantee could hold.
[TreeBH primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC9945647/)

A deterministic frontier replay tests a descendant sibling despite a parent
sibling p-value of 0.9 when edge prerequisites remain open. The live root is a
boundary with pass-through disabled and passes through with suitable deeper
evidence when enabled. Closing its child edges blocks both configurations.
The probe also demonstrates that applying correction to invalid `U²` p-values
does not repair their distribution. These checks concern algorithm boundaries;
they do not estimate complete-pipeline FDR.

## Knowledge corrections and unresolved implementation gaps

The [knowledge audit](knowledge_audit.md) records exact inspected locations and
claim verdicts. Corrections applied to current wiki synthesis:

- Removed the general independence explanation for leaf-only PCA; retained
  all historical results with their original settings and dates.
- Distinguished edge basis halves, full-leaf sibling context, optional sibling
  internal rows, physical spectral row count and MP leaf count.
- Distinguished MP count, nominal dimension, contrast-selected retained rank
  and nominal reference df; documented the missing-edge-basis sentinel.
- Qualified deformed MP as a candidate for a specified ensemble, corrected
  the exact-K oracle scope, and identified the active empirical plug-in tail.

The manuscript still has identified stale descriptions, including an exact-F
API removed on 2026-10-16, old PCA construction, and rank/reference conflation.
They are listed for a coordinated manuscript revision, rather than silently
changing a mathematical specification during this diagnostic.

Two software-level gaps are recorded separately from the statistical failures:
an augmented zero-rank branch reports physical rather than leaf threshold rows
although no threshold was evaluated, and direct low-level kernel callers can
supply nonorthonormal rows without a local orthonormality check. No malformed
production eigensolver output was observed. Neither gap explains the principal
reference-law failures; neither was changed in this audit.

## Reproduction and verification

From the repository root, with existing dependencies:

```bash
uv run --no-sync --offline python reports/statistical_construction_audit_20261104/mp_probe.py
uv run --no-sync --offline python reports/statistical_construction_audit_20261104/reference_probe.py
PYTHONPATH=. uv run --no-sync --offline python reports/statistical_construction_audit_20261104/downstream_probe.py
```

Agent checks: 37 MP/tie tests and 92 reference-law tests passed; these sets
overlap. Coordinator checks: 33 calibration-metadata/multiplicity/frontier tests
passed. Passing implementation tests does not establish a null law. Probe
assertions, native replays and an independent mathematical review support the
specific findings above. The [verification record](verification.json) checks
114 source hashes, including 109 production files unchanged since the preceding
cancer-PC diagnosis. Compact [MP](mp_results.json),
[reference-law](reference_results.json) and [downstream](downstream_results.json)
results are retained; their generated originals are under `generated/` and are
classified as local generated artifacts.

No dependencies were installed, production statistic changed, commit made or
service started. Full cancer data, full-scale solver behavior, general
selected-event certification and final-partition error control remain open.
