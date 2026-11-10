# The largest split one population allows: a rule in place of the allowed overlap

Analysis on 2026-11-08. Local folder: `results/analyses/gaussian_cut_rule_20261108` (git-ignored); protocol,
deviations, code and outputs in `single-population-split-bound-scripts-20261108.py.txt`.

## Question

Infer the mathematical rule for the allowed overlap against noise splits the tree picks, instead of an
empirical calibration, keeping the rule from tightening as the number of cells grows.

## Derivation

Split of n cells into shares pi, 1 - pi; D = Mahalanobis distance of the side means under the pooled
within-sides covariance; rho = pi (1 - pi) D^2 = between / within variance along the split's Fisher direction;
Hotelling T = n rho.

1. Among all splits with a given share, the one with the largest between/within ratio along a direction u is the
   threshold cut on u (it maximises the gap of the side means; total variance along u is fixed). D^2 is the
   maximum over u. So the largest rho over all 2^(n-1) - 1 splits is a flat cut. Check A2: in 20 sets of 12 points
   (all 2,047 splits each) the best split, overall and per smaller-side size, was linearly separable 140/140.
2. One Gaussian population, any covariance: a flat cut at share pi has rho = g(pi) = r / (1 - r),
   r = phi(Phi^-1(pi))^2 / (pi (1 - pi)). Balanced g = 1.752 (r = 2/pi), pi = 0.2: 0.960, pi = 0.05: 0.289
   (2 x 10^6 draws: 1.753, 0.961, 0.288). No split of one Gaussian population goes above g(pi).
3. Other one-population shapes go higher: uniform 3.00 at the middle, exponential 1.84 at an 80/20 cut, logistic
   1.41, lognormal(1) 1.21, t5 1.18, Laplace 1.00. For cuts at the mode of any unimodal population rho <= 3
   (Khintchine representation X = mode + U Y; Cauchy-Schwarz). A thin, long, flat tail cut a third of the way out
   reaches 8 (3 (1 - a)(1 + a)^2 / ((1 - a)^3 + 4 a^3) at a = 1/3; computed 7.998); a random search over 20,000
   unimodal mixtures found 7.84. So "unimodal" alone gives no useful ceiling; the ceiling is a statement about the
   shape of one population.
4. Two Gaussian groups at separation delta > 0: a flat cut through the gap goes above g(pi), so any delta > 0 is
   reached as n grows. If the split's labels are the true groups, rho = pi (1 - pi) delta^2, above g(pi) only for
   delta > delta_0(pi) = sqrt(g / (pi (1 - pi))): balanced 2.647, overlap (Bayes error) 9.3%; pi = 0.2: 2.45,
   8.0%; pi = 0.1: 2.40, 5.4%.
5. Overlap does not decide: a flat cut of one population has zero overlap and is separable. This is why the
   zero-overlap rule and the overlap-fraction rule both judge the wrong quantity.
6. Finite n, on the k invariant directions of the subtree: at the boundary (a fixed flat cut of one Gaussian
   population) T is about noncentral chi2_k with noncentrality n g(pi). Check A4 (fixed median cut, n = 1,000,
   4,000 replicates): share p <= 0.01 0.0095; mean T 1,751 vs 1,753; variance 6,727 vs 7,010 (noncentral chi2)
   and 6,582 (delta method), so the tail is slightly conservative. The tree can pick any split; by (1) the
   worst picks are flat cuts, M_k(n) = sum_{i<=k} binom(n-1, i) of them (Cover 1965). Rule:

       p = min{1, M_k(n) P(chi2_k(n bound(pi_hat)) >= n rho_hat)},  p = 1 when k = 0,

   with bound = g(pi) (one Gaussian population) or a constant (3 = the at-the-mode ceiling).

## Allowed overlap that follows (k = 1, p <= 0.01, true-label groups)

| share | 100 cells | 1,000 | 10,000 | 48,879 | limit |
| --- | --- | --- | --- | --- | --- |
| 0.5, Gaussian ceiling | 4.5% | 7.2% | 8.5% | 8.9% | 9.3% |
| 0.2, Gaussian ceiling | 3.4% | 6.1% | 7.3% | 7.7% | 8.0% |
| 0.1, Gaussian ceiling | 1.8% | 3.8% | 4.8% | 5.1% | 5.4% |
| 0.5, ceiling 3 | 1.8% | 3.1% | 3.8% | 4.0% | 4.2% |

The allowance rises with cells toward its limit (the required excess of rho over the ceiling shrinks like
sqrt(bound log n / n)); the zero-overlap rule's requirement grew instead (4.7 SD at 100 cells, 8.2 at 48,879).

## Tree-picked splits on simulated populations (Ward's top split, 30 coordinates, A5)

- Gaussian gradient (variance 9 on one axis): kept 3/200 at 64 cells, 0/200 at 128-4,096; Ward's splits have
  median rho 1.15-1.29, below g (1.50-1.61). Ceiling 3: 0 everywhere.
- Flat (uniform) gradient: kept under the Gaussian ceiling rises with n, 10.5% at 64 cells to 63.5% at 4,096
  (predicted: a flat gradient is not Gaussian); ceiling 3: 0 everywhere.
- Two groups at delta 4: kept rises with n (balanced 30% at 64 cells to 100% from 2,048; share 0.2: 23% to 100%);
  delta 3 at share 0.2: 1% to 13%; delta 3 balanced and delta 2: 0 (Ward's top split in 30 coordinates is too
  blurred, median rho 1.25-1.43). Ceiling 3: delta <= 2 stays out (one-mode mixtures); delta 4 balanced 2-14%.
- Old rule for contrast: Ward's splits are separable in 0-13% of runs at 64 cells and 0% from 256 cells on.

## End to end (935 simulated datasets, both trees, production fits; harness of `rule_out_overlap_20261108`)

Replay of the old rule: 1,870/1,870 datasets with the same cluster count; cell-count check 11,488/11,488.

| Markov/Ward | old rule (zero overlap) | overlap rule, formula fixed | Gaussian ceiling | ceiling 3 |
| --- | --- | --- | --- | --- |
| False partitions of 300 (sd 0.25 / 0.42 / 1.0) | 2 / 25 / 15 | 13 / 17 / 28 | 0 / 0 / 5 | 0 / 0 / 0 |
| Large no-group splits passing (sd 0.25 / 0.42 / 1.0) | 0.0 / 0.9 / 0.4% | 0.4 / 0.6 / 0.9% | 0 / 0 / 0.2% | 0 / 0 / 0 |
| ARI balanced / imbalanced / singleton / density | 1.00 / 1.00 / 0.88 / 0.87 | 0.99 / 0.87 / 0.73 / 0.55 | 1.00 / 1.00 / 1.00 / 0.57 | 1.00 / 1.00 / 1.00 / 1.00 |
| Median clusters, density variation (true 4) | 5 | 12 | 7 | 4 |

Default tree: false partitions 0 in every arm except the fixed overlap rule (6 at sd 1.0); real-group ARI
unchanged (density variation 0.667; 0.685 with the Gaussian ceiling).

Tree-picked splits with an invariant direction (no-group data, Markov/Ward): median rho / g 0.23-0.50 by size;
3-14% above g; 0-0.3% above 3. Real-group splits of 400+ cells: median rho / g 170-214. Under the Gaussian
ceiling, splits inside a density-variation group pass: such groups are skewed, and a skewed population exceeds
g(pi) at unbalanced cuts (as the exponential does), but not 3.

## Reading

- The rule is not an overlap fraction. A split is evidence of two groups when its between/within ratio along its
  own discriminant exceeds the most that a cut of one population can produce; the overlap allowance follows from
  that ceiling and grows with cells.
- The ceiling depends on what "one population" may look like: 1.75 (balanced) for Gaussian, 3 at the mode for
  any unimodal shape, up to 8 for long flat tails. On these simulations the constant 3 removed every false
  partition and recovered all four simulated group types with Markov/Ward; the Gaussian ceiling over-splits
  skewed groups.
- Not measured: real CHOIR data, the continuous benchmark, and populations with long flat tails, where 3 is not a
  ceiling. No production change in this task.

## Correction to the overlap redesign

`allowed_overlap` in `rule_out_overlap_20261108/overlap_rule.py` had its two Bayes-error terms swapped. It is
right for balanced splits (11.9%) but too large otherwise: 7.3% instead of 1.3% at a 10/90 split (2 x 10^6
simulated cells: 1.32%). With the formula fixed, Markov/Ward false partitions are 13 / 17 / 28 (was 13 / 42 / 36)
and real-group ARI 0.99 / 0.87 / 0.73 / 0.55, still worse than the zero-overlap rule.

## Literature (looked up)

- Liu, Hayes, Nobel, Marron 2008, Statistical Significance of Clustering for High-Dimension, Low-Sample Size
  Data, J. Am. Stat. Assoc. 103(483):1281-1293, doi:10.1198/016214508000000454 (Crossref record). SigClust also
  tests a split against a single Gaussian population (2-means cluster index); not used in the code here.
- Cover 1965, Geometrical and Statistical Properties of Systems of Linear Inequalities with Applications in
  Pattern Recognition, IEEE Trans. Electron. Comput. EC-14(3):326-334, doi:10.1109/PGEC.1965.264137 (Crossref
  record). The count formula used here was not checked against the paper (unverified).
- Ion, Klaassen, van den Heuvel 2022, Sharp Inequalities of Bienayme-Chebyshev and Gauss Type for Possibly
  Asymmetric Intervals around the Mean, arXiv:2208.08813 (title page read; keyword "Khintchine representation").
  The form Z = M + U Y was seen in a search snippet of its appendix; Khintchine's 1938 paper was not retrieved
  (unverified).
