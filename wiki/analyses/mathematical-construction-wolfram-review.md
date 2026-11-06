---
title: Mathematical Construction Wolfram Review
type: analysis
status: reviewed
updated: 2026-11-02
sources:
  - raw/inbox/wolfram-mathematical-construction-checks-20261102.md
  - README.md
  - manuscript/sections/method/representation.tex
  - manuscript/sections/method/edge_test.tex
  - manuscript/sections/method/sibling_test.tex
  - manuscript/sections/method/decomposition.tex
  - tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py
  - tree_break_selection/hierarchy_analysis/decomposition/backends/eigen/operators.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py
  - tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py
  - tree_break_selection/hierarchy_analysis/statistics/sibling_divergence/inflation_correction/empirical_null_inflation_estimation.py
  - tree_break_selection/hierarchy_analysis/decomposition/gates/orchestrator.py
  - tree_break_selection/hierarchy_analysis/decomposition/gates/gate_evaluator.py
  - tests/statistics/49_test_projected_wald_reference_distribution.py
  - tests/statistics/42_test_contrast_covariance.py
  - tests/hierarchy_analysis/decomposition/gates/test_gate_evaluator.py
tags:
  - mathematics
  - review
  - wolfram
  - calibration
---

# Mathematical Construction Wolfram Review

## Summary

The barycentric and fixed-membership covariance algebra is consistent. The
implemented projected quadratic is well-defined, but its chi-square law needs
isotropy and a fixed or appropriately independent projection. Same-data
selection, estimated covariance, and empirical scale fitting require further
statistical justification. This review records those boundaries and three
current manuscript/code differences; it does not modify the method.

## Details

### Construction map

\[
(X,\mathcal F,\mathcal T)
\longrightarrow \{L(u),n_u,\theta_u\}
\longrightarrow \{\Delta,\Sigma_0,z\}
\longrightarrow \{P,k,T=\|Pz\|^2\}
\longrightarrow \{p_{\mathrm{edge}},\hat c,p_{\mathrm{sib}}\}
\longrightarrow \{\text{Tree-BH, sibling BH, traversal}\}
\longrightarrow \mathcal C_{\mathrm{TBS}}.
\]

The middle reference-law step is conditional on statistical assumptions. The
last step is an algorithm reading annotations. Neither the existence of a
partition nor a passing software test establishes its inferential error control.
[[mathematical-glossary]] defines each object without renaming existing concepts.

### Verified algebra

Wolfram evaluations 1 and 3 in the raw capture establish:

1. **Barycenter contrasts.** For positive child counts \(a,b\) and
   \(\mu=(ax+by)/(a+b)\),
   \(x-\mu=b(x-y)/(a+b)\) and \(y-\mu=-a(x-y)/(a+b)\).
2. **Nested covariance cancellation.**
   \([b/(a+b)]^2(1/a+1/b)=1/a-1/(a+b)\).
   Combined with independent leaf means of common covariance, this derives
   the nested child-parent covariance rather than treating child and parent
   as independent samples.
3. **Standardized direction identity.** With a common positive variance,
   no branch multipliers and no ridge, the standardized left edge equals the
   sibling contrast; the right edge is its negative. This is a coordinate
   check; compatible common matrix whitening gives the vector form.
4. **Split action.**
   \(a(x-\mu)^2+b(y-\mu)^2=ab(x-y)^2/(a+b)\).
   Coordinate summation gives the Euclidean vector identity.
5. **Categorical interior covariance.** For a three-category block, the
   two-coordinate covariance determinant is
   \(q_1q_2(1-q_1-q_2)>0\) in the simplex interior. This verifies the
   two-dimensional example; the boundary remains singular without regularization.
6. **Branch-time qualification.** If both edge and sibling branch factors
   are applied with the same covariance baseline, the squared ratio of
   standardized left edge to sibling contrast is
   \[
   \frac{1+(t_L+t_R)/(2\bar t)}{1+t_L/\bar t}.
   \]
   Thus the sampling-only identity does not hold unchanged for arbitrary
   unequal branch times. Under the default edge policy `none`, the denominator
   branch factor is absent.

### Reference law: two separate requirements

Follow-up inspection of the annotation caller established a distinction missing
from the initial overview: the edge basis uses alternate rows within each
selected child, while edge contrasts use the remaining rows. The sibling basis
uses all leaves. The kernel-level selection warning below remains relevant,
but direct same-row PCA reuse is not the current default edge path. Splitting
after hierarchy selection does not establish independent membership, and the
continuous covariance resolver can use all leaves. See
[[analytic-selected-root-laws]] for the subsequent restricted analytic model.

**Covariance requirement.** If \(z\sim N(0,R)\) and fixed \(P\) has
orthonormal rows, the projected covariance is \(PRP^\top\). For eigenvalues
\(\eta_j\) of that covariance,
\[
\|Pz\|^2\ \stackrel d=\ \sum_{j=1}^k\eta_j Z_j^2,
\qquad Z_j\ \text{independent }N(0,1).
\]
The ordinary chi-square law requires \(PRP^\top=I_k\). Coordinatewise unit
variance is insufficient. Wolfram verified that
\[
R=\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix},\qquad
P=(1,1)/\sqrt2
\quad\Longrightarrow\quad PRP^\top=1+\rho.
\]
Hence \(T\sim(1+\rho)\chi_1^2\) for this fixed projection when
\(-1<\rho<1\). With full two-dimensional projection the weights are
\(1-\rho,1+\rho\), yielding mean 2 and variance \(4(1+\rho^2)\), rather
than variance 4 unless \(\rho=0\). The evaluator differentiated the weighted
chi-square MGF to check these moments. Declared block whitening must therefore
be distinguished from verified whitening of the actual joint contrast.

**Selection requirement.** Even with \(z\sim N(0,I_d)\), choosing
\(P(z)=z^\top/\|z\|\) gives an orthonormal one-row projection with
\(T=\|z\|^2\sim\chi_d^2\), rather than \(\chi_1^2\) when \(d>1\).
Wolfram checked the two-coordinate algebra. The theorem can hold for a
random independent projection conditionally, but not from orthonormality alone.
The current kernel may choose both PCA directions and contrast-energy rank
from data also involved in the test; a fixed-projection proof does not cover this.

### Selection and empirical scale do not supply a tail theorem

For the illustrative selection event \(W>c>0\), with \(W\sim\chi_2^2\),
Wolfram verified
\[
\Pr(W>w\mid W>c)=e^{-(w-c)/2},\qquad w\ge c.
\]
The naive unselected p-value \(e^{-W/2}\) has conditional mean
\(e^{-c/2}/2\), rather than \(1/2\). This is a counterexample, not a derived
law for TBS's actual hierarchy event.

Wolfram also verified that \(\chi_2^2\) and \(2\chi_1^2\) share mean 2 but
have variances 4 and 8. Matching an empirical inflation scale cannot establish
the full selected-null tail. The current code explicitly describes its
plug-in chi-square tail as historical and its same-selected-hierarchy reference
law as unvalidated. `internal_admissible` reports support-contract state; it
does not establish a mathematical calibration theorem.

### Current manuscript/code differences

| Item | Manuscript | Current implementation / consequence |
| --- | --- | --- |
| Edge branch-time policy | `edge_test.tex` presents the branch-inflated variance as its main formula. | `orchestrator.py` defaults to `none`; branch inflation is opt-in for edges. The adaptive sibling test still applies its own branch factor. A derivation should name the policy. |
| Parent spectral matrix | `edge_test.tex` describes a Pearson-correlation matrix from unit-variance standardized rows. | The backend centers already null-whitened tangent rows and computes covariance without another sample-variance standardization. The glossary uses the current covariance definition. |
| Pass-through descendant condition | `decomposition.tex` defines a significant sibling at any proper descendant where tested. | `GateEvaluator` requires a descendant full split reachable through open intermediate prerequisites and any pass-through support. A blocked intermediate node hides its subtree; an existing test checks this case. |

The manuscript also mixes parent-to-child edge subscripts with a child-to-parent
label for the standardized contrast in `representation.tex`. The numerator
\(\theta_c-\theta_u\) fixes the meaning; the glossary records the difference.

The pre-existing dirty implementation of
`empirical_null_inflation_estimation.py` includes a family-baseline cap when
predicting above every calibration parent size, if the local estimate exceeds
the family baseline. This review describes the live working tree; it does not
attribute that rule to a clean commit or verify its statistical validity.

### Verification and limits

- Existing focused baseline: 48 tests passed in 1.31 seconds across wiki lookup,
  projected-reference helper, contrast covariance and gate evaluator.
- Baseline `make wiki-lint`: 214 pages passed before documentation additions.
- Environment: macOS, zsh task shell; existing project `.venv`, Python 3.12.11,
  NumPy 2.3.4, SciPy 1.16.3, pytest 8.4.2, uv 0.9.5; no installation or service
  launch was needed. The existing `uv run --no-sync` workflow was used.
- Wolfram's context tool failed internally. Its evaluator succeeded; the raw
  capture preserves all three calls, including correction of the second call's
  unevaluated scaled-distribution expression.
- No full benchmark, simulation calibration, manuscript compilation, selected
  conditional-law derivation, or method-code changes were performed. Tests
  exercise the formulas and decisions, not the validity of the selected law.

## Evidence

- `raw/inbox/wolfram-mathematical-construction-checks-20261102.md` contains
  exact inputs and returned outputs. Only evaluation results, not assigned
  illustrative expressions in evaluation 2, are treated as computation evidence.
- The cited live implementation was traced from contrasts through the kernel,
  empirical inflation, multiplicity and gate evaluation. The manuscript was
  checked separately rather than assumed to match the code.
- `README.md` explicitly retains selected-hierarchy calibration as an open
  research question. Earlier synthesis in [[method-proof-web]] is consistent
  with that boundary but is not used as current runtime proof.

## Links

- [[mathematical-glossary]]
- [[wolfram-mathematical-construction-checks]]
- [[method-proof-web]]
- [[projected-wald-statistic]]
- [[empirical-null-calibration-reference-law-contract]]
- [[open-mathematical-questions]]

## Open Questions

- Establish actual joint covariance and the role of estimated covariance,
  ridge and small discrete samples before claiming isotropic Gaussian contrasts.
- Derive or validate the law conditional on hierarchy, projection/rank, edge
  opening, focal context and empirical-null fitting.
- State applicable multiplicity and final-partition guarantees for that law.
- Reconcile the three manuscript/code differences above in a separate scoped
  change; this glossary/review records them without altering source behavior.
