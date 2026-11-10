# A tighter count than all flat cuts: the tail of one population's best cut

Analysis on 2026-11-09 (user choice: derive a tighter count of the splits a Ward tree can pick). Local folder:
`results/analyses/flat_cut_sup_count_20261109` (git-ignored); protocol, deviations, code and outputs in
`flat-cut-sup-count-scripts-20261109.py.txt`.

## Derivation

1. Whatever rule the tree uses, its split's ratio is at most the best flat cut's. So p = P(best flat cut of one
   population >= rho_obs) is valid for Ward without counting cuts. The old rule bounded this by a union over all
   M_k(n) ~ n^k flat cuts, which ignores that neighbouring cuts share almost all cells.
2. In coordinates whitened by the total covariance, for the cut through the centre along v,
   R(v) = mean|v'x| / sd(v'x) and rho = R^2/(1 - R^2) (times (n - 2)/n). Over all sign splits, the best |mean of
   signed cells| equals the best direction's mean absolute deviation.
3. To first order the best cut sits on the population's best cuts (envelope argument: Hartigan 1978, Pollard
   1982); thresholds and side shares add O(1/n). Per direction sqrt(n)(R(v) - R*) is the mean of
   f = |t| - c0 (t^2 + 1)/2: variance 1 - 3/pi (Gaussian), 1/40 (uniform).
4. Gaussian population: every direction is a best cut, so the best cut is the maximum of an even Gaussian field
   on the sphere S^{k-1}; correlation (2/pi)(sqrt(1 - r^2) + r arcsin r - 1 - r^2/2) / (1 - 3/pi), second
   spectral moment lambda = (1 - 2/pi)/(1 - 3/pi) = 8.06. Ceiling 3: the uniform law on the sphere of radius
   sqrt(3) in k = 3 has a uniform marginal in every direction; its lambda = 13.0 is used for every k
   (assumption: the ceiling alone does not fix lambda).
5. Tail of the maximum: expected Euler characteristic with the Lipschitz-Killing curvatures of the sphere of
   radius sqrt(lambda) (Adler-Taylor notes, eq. 3.6.4; EC densities (2 pi)^{-(j+1)/2} H_{j-1}(u) e^{-u^2/2}),
   halved for the even field; on the R scale (amended, below); as a non-increasing envelope (amended, below).
   Self-check: with lambda = 1 the formula equals P(chi_k >= u) (Adler-Taylor notes, eq. 4.5.15).

## Checks

- Frozen rule (u on the rho scale) failed the null check: 2-22% of single populations at p <= 0.01, worst for
  k = 3 at 100 cells. Amended to the R scale (rho is convex in R). Fresh seeds, 1,000 replicates, k = 1-3,
  100-1,600 cells, Gaussian and uniform-marginal nulls: 0-1.3% at p <= 0.01, 1.0-7.4% at p <= 0.05 (old count:
  about 0%). One-dimensional Gaussian, 4,000 replicates: 1.3-1.5% at 1%, shrinking with n.
- For k >= 4 the EC expansion swings sign at moderate u; clipping gave p = 0 at the ceiling. Fixed with a
  non-increasing envelope before any result was used.
- k = 6 and 9 (best cut found by search, a lower bound): 0-0.3% at 1%. Null 99th percentile of the best cut vs
  the rho needed for p = 0.01:

| Null, k, cells | 99th percentile | sup rule needs | flat-cut count needs |
| --- | ---: | ---: | ---: |
| Gaussian, 6, 400 | 2.58 | 2.66 | 2.96 |
| Gaussian, 9, 1,600 | 2.18 | 2.23 | 2.52 |
| uniform product, 6, 400 | 3.86 | 4.74 | 4.53 |
| uniform product, 9, 1,600 | 3.38 | 3.90 | 3.98 |

  For the Gaussian null the sup rule is close to exact and tighter than the count; for ceiling 3 with
  lambda = 13 it is no tighter at 400 cells.

## On the saved Markov/Ward ceiling3 tree (descriptive, no refit)

- The 27 hidden-boundary pairs: none of their top splits passes either rule. Needed rho, old 4.3-7.4, new 4.3
  to infinite (36 and 57 cells with 6-9 directions cannot pass: R is at most 1).
- Needed rho with ceiling 3, k = 6 / 9 / 12: 100 cells 9.5 / 16.2 / 37.3; 300 cells 5.2 / 6.0 / 7.0; 1,000
  cells 3.9 / 4.2 / 4.5; 5,000 cells 3.4 / 3.5 / 3.5.
- 1,838 tested splits kept at p <= 0.01 (old -> new): 32-159 cells 14 -> 1, 160-999 65 -> 62, 1,000+ 42 -> 43.

## Reading

- The count of all flat cuts can be replaced by the tail of one population's best cut, an even Gaussian field
  on the sphere of directions. For a Gaussian population this is near exact and tighter than the count; its
  margin shrinks like 1/sqrt(n) instead of sqrt(k log n / n).
- For ceiling 3 the field's smoothness lambda depends on the shape of the least favourable population, which the
  ceiling does not fix; with the k = 3 value (13.0) the rule is no tighter than the count below about 1,000
  cells and only tighter above.
- At the sizes of the merged real clusters (100-500 cells, 8-12 directions) the best cut of one population
  itself reaches rho 3.9-4.1 (uniform-product null, 400 cells, 99th percentile), so even an exact count needs
  rho near 4 or more. The tree's splits there reach 1.6-4.5; two clean-boundary splits (4.38 at 258 cells,
  4.53 at 360) are the only ones near that line. The binding problem is the tree's blurred splits, not the count.
- Not done: end-to-end runs with the new rule, a least-favourable law for ceiling 3 in k > 3, a Ward-specific
  count (selective inference conditions on the tree instead, Gao et al. 2024). No production change.

## Literature (looked up)

- Pollard 1982, A Central Limit Theorem for k-Means Clustering, Ann. Probab. 10(4):919-926,
  doi:10.1214/aop/1176993713 (Crossref; abstract via search results: within-groups sum of squares, empirical
  process CLT). Full text not read.
- Hartigan 1978, Asymptotic Distributions for Clustering Criteria, Ann. Statist. 6(1):117-131,
  doi:10.1214/aos/1176344071 (Crossref; not read).
- Adler, Taylor, Applications of Random Fields and Geometry (draft notes, web.stanford.edu/class/stats317/hrf.pdf),
  read: eq. 3.6.4 (L_j of the sphere of radius lambda: 2 binom(N-1, j) s_N / s_{N-j} lambda^j when N - 1 - j is
  even), eq. 4.5.15 (P(chi_k >= u) = sum_j L_j(S^{k-1}) rho_j(u)), EC densities (2 pi)^{-(j+1)/2} H_{j-1}(u)
  e^{-u^2/2}. Book: Adler, Taylor 2007, Random Fields and Geometry, Springer, doi:10.1007/978-0-387-48116-6
  (Crossref).
- Taylor, Takemura, Adler 2005, Validity of the expected Euler characteristic heuristic, Ann. Probab.
  33(4):1362-1396, doi:10.1214/009117905000000099 (Crossref; not read).
