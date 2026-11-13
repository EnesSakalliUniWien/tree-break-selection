# Marchenko-Pastur, geometrical and algebraic multiplicity, and the definition of k

Lookups and run on 2026-11-12, Europe/Vienna. Local folder:
`results/analyses/marchenko_pastur_multiplicity_20261112` (git-ignored); protocol, script and output of the four
numerical checks are in `marchenko-pastur-multiplicity-scripts-20261112.py.txt`. No production code was changed.
Environment: Python 3.12.11, NumPy 2.3.4, uv 0.9.5.

## Subject

k is the number of eigenvalues of a node's covariance above
`edge = median eigenvalue / MP median(p/n) x (1 + sqrt(p/n))^2` (`invariant_directions` in
`tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/split_count_rule_out.py`), and 0 when
n <= p + 1. The question is what the literature on the Marchenko-Pastur law and on eigenvalue multiplicity says
about this definition.

## Looked-up sources

Every entry was retrieved on 2026-11-12; "confirmed" names what the retrieved text states. Nothing is taken
from memory.

1. V. A. Marcenko, L. A. Pastur, "Distribution of eigenvalues for some sets of random matrices", Mathematics
   of the USSR-Sbornik 1(4), 457-483, 1967. https://doi.org/10.1070/SM1967v001n04ABEH001994
   Confirmed: bibliographic record only (Crossref). The paper's content was not retrieved; the law is used
   below as restated in sources 4 and 5.
2. I. M. Johnstone, "On the distribution of the largest eigenvalue in principal components analysis", The
   Annals of Statistics 29(2), 295-327, 2001. https://doi.org/10.1214/aos/1009210544
   Confirmed (abstract in the Crossref record): for a Gaussian matrix with both dimensions large, the largest
   eigenvalue, centred and scaled, approaches the Tracy-Widom law of order 1; simulations suggest the
   approximation is usable for dimensions as small as 5.
3. J. Baik, G. Ben Arous, S. Peche, "Phase transition of the largest eigenvalue for nonnull complex sample
   covariance matrices", The Annals of Probability 33(5), 1643-1697, 2005.
   https://doi.org/10.1214/009117905000000233 ; https://arxiv.org/abs/math/0403022
   Confirmed (arXiv abstract): complex Gaussian samples; when all but finitely many population eigenvalues are
   equal, the limiting law of the largest sample eigenvalue depends on those eigenvalues and a phase
   transition is observed. The abstract states no critical value.
4. J. Baik, J. W. Silverstein, "Eigenvalues of large sample covariance matrices of spiked population models",
   Journal of Multivariate Analysis 97(6), 1382-1408, 2006. https://doi.org/10.1016/j.jmva.2005.08.003 ;
   https://arxiv.org/abs/math/0408165
   Confirmed (arXiv text, pages 1 to 4 read): with p/n -> c, the sample eigenvalues of i.i.d. samples follow
   the Marchenko-Pastur density on [(1 - sqrt(c))^2, (1 + sqrt(c))^2] and the largest converges to
   (1 + sqrt(c))^2. In the spiked model (all population eigenvalues 1 except finitely many alpha_j with
   multiplicities k_j; entries with finite fourth moment, not necessarily Gaussian), Theorem 1.1 (c < 1): for
   alpha_j > 1 + sqrt(c) all k_j associated sample eigenvalues converge almost surely to
   alpha_j + c alpha_j / (alpha_j - 1); a population eigenvalue at or below 1 + sqrt(c) leaves the largest
   sample eigenvalue at (1 + sqrt(c))^2. Theorem 1.2 states the same limit for c > 1.
5. Z. Bai, J.-F. Yao, "Central limit theorems for eigenvalues in a spiked population model", Annales de
   l'Institut Henri Poincare - Probabilites et Statistiques 44(3), 447-474, 2008.
   https://doi.org/10.1214/07-AIHP118 ; https://arxiv.org/abs/0806.2503
   Confirmed (arXiv text, pages 1 to 3 read): the limiting distribution of a sample eigenvalue of a spike is
   Gaussian only when the population spike is simple; for a spike of multiplicity k, the k consecutive sample
   eigenvalues converge jointly to the eigenvalue distribution of a k x k symmetric Gaussian random matrix,
   so a single one is in general not Gaussian. The critical values are 1 + sqrt(y) and 1 - sqrt(y).
6. M. Gavish, D. L. Donoho, "The Optimal Hard Threshold for Singular Values is 4/sqrt(3)", IEEE Transactions on
   Information Theory 60(8), 5040-5053, 2014. https://doi.org/10.1109/TIT.2014.2323359 ;
   https://arxiv.org/abs/1305.5870
   Confirmed (arXiv abstract): low-rank matrix in noise, hard threshold on singular values; with unknown noise
   level the recommended threshold is about 2.858 times the median empirical singular value.
7. M. Okamoto, "Distinctness of the Eigenvalues of a Quadratic form in a Multivariate Sample", The Annals of
   Statistics 1(4), 763-765, 1973. https://doi.org/10.1214/aos/1176342472
   Confirmed (abstract in the Crossref record): under a rank condition on the defining matrix and a sample
   distribution absolutely continuous with respect to Lebesgue measure, the nonzero eigenvalues of the
   quadratic form are distinct with probability one.
8. D. Passemier, J.-F. Yao, "Estimation of the number of spikes, possibly equal, in the high-dimensional case",
   Journal of Multivariate Analysis 127, 173-183, 2014. https://doi.org/10.1016/j.jmva.2014.02.017 ;
   https://arxiv.org/abs/1110.1364
   Confirmed (arXiv abstract): an estimator of the number of spikes in high dimension, based on random matrix
   theory, extended to the case where some spikes are equal.
9. S. Kritchman, B. Nadler, "Non-Parametric Detection of the Number of Signals: Hypothesis Testing and Random
   Matrix Theory", IEEE Transactions on Signal Processing 57(10), 3930-3941, 2009.
   https://doi.org/10.1109/TSP.2009.2022897
   Confirmed (institutional record of the abstract): the number of signals is detected by a sequence of
   hypothesis tests built on random matrix theory. Unverified: that the thresholds are Tracy-Widom quantiles;
   this appears only as a keyword in a search listing, not in a retrieved text.
10. T. W. Anderson, "Asymptotic Theory for Principal Component Analysis", The Annals of Mathematical Statistics
    34(1), 122-148, 1963. https://doi.org/10.1214/aoms/1177704248
    Confirmed: bibliographic record only. Content not retrieved; no statement is taken from it.
11. Definitions: https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors and
    https://en.wikipedia.org/wiki/Symmetric_matrix
    Confirmed: algebraic multiplicity is the multiplicity of the eigenvalue as a root of the characteristic
    polynomial; geometric multiplicity is the dimension of its eigenspace and is never larger; a matrix is
    diagonalizable when the two agree for every eigenvalue and defective otherwise; every real symmetric
    matrix is diagonalizable by an orthogonal matrix with real eigenvalues.

## Geometrical and algebraic multiplicity in this project

- A covariance matrix is real symmetric, so for each of its eigenvalues the geometrical and the algebraic
  multiplicity are equal (source 11). The distinction cannot occur for the spectrum behind k.
- The Markov operator D^-1 K equals D^-1/2 (D^-1/2 K D^-1/2) D^1/2, a similarity transform of a symmetric
  matrix (one-line identity; the wiki records Phi_M = D^-1/2 Phi_S), so it is diagonalizable too. A defective
  eigenvalue could only come from a kernel that is not symmetric; whether every kernel path symmetrizes was not
  checked.
- Three different things are called multiplicity in the project:
  1. the multiplicity of a population eigenvalue (a model property; four balanced groups give 3);
  2. numerical ties in a computed spectrum (`tied_eigenvalue_neighbors`, relative tolerance 1e-10);
  3. an estimated multiplicity (the block equality test of 8 November).
- By source 7, ties of kind 2 have probability zero for data with a density; by sources 4 and 5 a population
  eigenvalue of multiplicity m shows as m sample eigenvalues near one limit, not as a tie. Measured before:
  no numerical tie in 9,250 blocks on the CHOIR data. The exact ties found on 12 November (29 tied eigenvalues;
  30 minus test rows) are produced by whitening rows with a covariance computed from the same rows, not by the
  population.

## Four numerical checks (protocol written before the run)

| Check | Setting | Result |
| --- | --- | --- |
| C1 | One Gaussian population, independent sample, 2,000 draws; share with k >= 1 | (64, 30) 23.2%; (128, 30) 18.2%; (256, 30) 17.7%; (800, 30) 16.0%; (200, 80) 21.6%; largest k 2 |
| C2 | 64 and 256 cells, 30 features, sqrt(p/n) = 0.306; share with k >= 1 by population `rho` | 0.15: 19.7%; 0.306: 31.5%; 0.6: 89.9%; 1.44: 100% |
| C3 | 30 cells in 40 features; current rule against a prototype rule on the nonzero eigenvalues | current rule 0 of 1,000 in all four settings; prototype: one population 21.5%; `rho` 2.3: 75.4%; `rho` 5 and 50: 100% |
| C4 | Four balanced groups (population multiplicity 3), 30 features, 500 draws | numerical ties 0 of 500 at 400 cells and 0 of 500 at 4,000; relative spread of the three eigenvalues 0.115 and 0.038 at the median; k = 3 in 499 and 498 of 500 |

All four predictions of the protocol hold.

## What this says about the definition of k

1. The edge is a limit, not a quantile (source 4), and the largest eigenvalue fluctuates around its centring
   (source 2). With no margin, an independent sample of one Gaussian population gets a direction in 16% to 23%
   of draws (C1), before any selection by the tree.
2. The median convention is the unknown-noise device of source 6, which is stated for a low-rank matrix in
   noise. It is not derived for a subtree chosen by the tree.
3. Own derivation, supported by C2: two groups with shares pi and 1 - pi at distance d add the population
   eigenvalue 1 + pi(1 - pi) d^2, so by source 4 a direction exists in the limit only when the population
   `rho` = pi(1 - pi) d^2 exceeds sqrt(p/n). This states when "no invariant directions" is the correct answer.
4. Source 4 has the same threshold for more features than cells (Theorem 1.2). The rule "k = 0 when
   n <= p + 1" is therefore not required by the theory; a prototype on the nonzero eigenvalues finds the
   direction in 75% to 100% of draws where the current rule returns none (C3). This is the layout of the
   benchmark losses (30 cells in 40 features, 30 in 30, 40 in 20,000).
5. Multiplicity does not enter k through ties: in data the eigenvalues are distinct (source 7, C4). A
   population multiplicity can only be estimated (sources 5 and 8).

## Not checked

Sources 1 and 10 were verified as records only. Source 9's use of Tracy-Widom quantiles is unverified. No
source on spectra of subsets selected by clustering was looked up. C1 to C4 use identity noise and Gaussian
data, one seed; the prototype rule of C3 is not a proposal for production and its null rate (21.5%) has the
same excess as C1. Whether every kernel path is symmetric was not checked.
