# Known-branch skew: controlled diagnostic protocol

Frozen design before the main evaluation, 2026-09-28. This implements a bounded
Gaussian experiment from Q10, Q17, Q19 and Q24–Q25 of
`raw/inbox/recursive-method-program-20260604.txt`. It locates statistical
distortion in controlled constructions; it does not tune or change production.

## Known tree and population contrasts

The fixed tree is `R -> (A, B)`, `A -> (A0, A1)`, `B -> (B0, B1)`.
Each terminal group contains 32 independent observations in eight dimensions.
Its population covariance is either identity or AR(1), `Sigma[i,j]=0.6**abs(i-j)`.
All groups in a dataset share that covariance. Set `mu=1.25*e1`, `nu=e2`.

| Pattern | A0 | A1 | B0 | B1 | True signal nodes |
|---|---|---|---|---|---|
| global_null | 0 | 0 | 0 | 0 | none |
| root_only | mu | mu | -mu | -mu | R |
| A_deep_only | mu | -mu | 0 | 0 | A |
| separate_deep_branches | mu | -mu | nu | -nu | A, B |
| strong_root_weak_A | mu+0.25nu | mu-0.25nu | -mu | -mu | R, A |

Population-null truth is equality of the children's population mean vectors.
Mixed descendants may therefore form a null sibling contrast. Preserve the
actual vectors, child sample sizes and truth at each fixed node. At selected
nodes, derive truth from exact group counts and the declared mean coefficients,
not sample means, majority labels, edge decisions or floating-point equality of
separately accumulated averages.

## Stage comparisons

1. **Covariance:** fixed R/A/B membership and full whitened coordinates;
   compare known Sigma with sample pooled covariance, including between-child
   mean variation. This latter estimator is not within-child residual covariance.
2. **Projection/MP:** at the same memberships, cross each covariance estimator
   with full coordinates versus same-data PCA/MP. Use every raw descendant
   observation and the same whitened coordinate frame for PCA and the contrast.
   Use production eigendecomposition, MP dimension selection with diagnostic
   minimum dimension 1, complete numerical tie groups and the Wald kernel.
   Record actual rank and raw MP signal count. Covariance/projection interactions
   are reported explicitly rather than attributed to one factor.
3. **Selection:** construct an average-linkage Euclidean tree from the same
   128 observations and test its selected root with known Sigma/full coordinates.
   The selected root is a different membership estimand. Report it separately;
   global-null datasets permit a clean null comparison, but signal-selected
   memberships cannot be treated as the original fixed node.
4. **FDR/traversal:** inject fixed-tree known-Sigma full-coordinate p-values
   into production sibling BH and TreeDecomposition. Deliberately open every
   edge prerequisite. Evaluate pass-through disabled and enabled. Count actual
   reached `split` decisions, retain `pass_through` separately, and distinguish
   unreached nodes from non-rejections. This isolates algorithm behavior without
   empirical inflation or estimated edge selection.

Production Cholesky covariance whitening is used consistently for the contrast
and PCA input. Ridge is zero in this positive-definite Gaussian design: adding
ridge before versus after the contrast variance factor would change the two
coordinate frames and perturb the otherwise tied sample-whitened spectrum.
Singular covariance configurations are unavailable rather than silently
regularized. This is not a claim that the entire construction equals the
production pipeline's covariance policies or sibling dimension policy.

## Exact control

For fixed groups, independent Gaussian observations with known covariance give
`z = solve(cholesky(Sigma), mean_L-mean_R) / sqrt(1/n_L+1/n_R)` and
`T=z.T@z`. Under the declared mean-null, `T` has a chi-square law with eight
degrees of freedom; signal nodes have a noncentral law. Never compare signal
raw-statistic skew to a central-null reference.

The balanced fixed-tree R/A/B contrasts are independent Gaussian contrasts.
With all edges open and all nodes null, root and depth-one BH are independent
families. Any reached false-split probability is alpha without pass-through,
and `1-(1-alpha)**2` with pass-through: 0.01 and 0.0199 respectively.
These references apply only to this controlled three-test experiment.

## Evaluation and reporting

- Main design: 1,000 independent datasets per pattern/covariance profile,
  10,000 datasets total, seed 20260928, alpha 0.01. A small pilot checks runtime
  and implementation only; it is excluded from evaluation. If runtime requires
  a smaller evaluation, record the amendment before the main run.
- Reuse each dataset across its stage variants. Record the seed, input hash,
  command/configuration, source hashes and package versions. Preserve pilot and
  evaluation in separate output directories; refuse to overwrite prior results.
- Report null rejection counts/rates and pointwise Wilson 95% intervals per
  named node and variant. Report signal power separately. Each such count has
  at most one observation per independent dataset. Do not pool overlapping
  nodes as independent trials.
- Compare raw-statistic skew to the appropriate central chi-square df-mixture
  reference only on null rows. Report p-value quantiles/tails and actual rank.
- Save matched per-dataset changes between covariance/projection variants.
  For traversal, report any false split, split FDP, planted-node recovery,
  pass-through behavior and final terminal-group partition count.
- Retain unsupported and failed outcomes with reasons and explicit denominators.
  Intervals and stage differences are descriptive; do not declare general
  calibration or causal sufficiency from this finite panel.

Categorical/Bernoulli covariance, feature dependence beyond Gaussian AR(1),
branch-length estimation, diffusion, other tree inference methods, empirical
inflation and full production error control remain outside this bounded run.
