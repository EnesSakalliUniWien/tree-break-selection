# Conditional selective inference for diffusion root splits

## Model and statistic

This experiment uses independent Gaussian rows with known covariance
\(\sigma^2 I_d\). For the observed root children define
\(v_i=1/n_L\) on the left and \(v_i=-1/n_R\) on the right, and
\(e=v/\|v\|\). The tested null is \(e^T\mu=0\). Set

\[
r_{obs}=\|X^Te\|/\sigma,\qquad
u=X^Te/\|X^Te\|,\qquad A=(I-ee^T)X.
\]

Under the null, \(X^Te/\sigma\) is standard multivariate Gaussian and
independent of \(A\). Its direction and radius are independent. Conditioning
on \(A,u\) leaves the path \(X(r)=A+\sigma eru^T\), with \(r\sim\chi_d\).
The scalar statistic is \(r^2\), using the full input feature dimension.

Let \(M(X)\) be the complete native symmetric-diffusion average-linkage merge
history. The diagnostic conditions on
\(S=\{r>0:M(X(r))=M(X)\}\). Every radius refits the neighbors, bandwidths,
epsilon, kernel, eigenspectrum, diffusion coordinates and linkage. This
preserves one observed hierarchy in the selection event. Conditioning on all
merges is stronger than conditioning only on its root partition and can lose
power. NNLS branch lengths, fitted empirical inflation and adaptive projected
Wald ranks are outside this prototype statistic.

The conditional null tail is

\[
p_{sel}=\frac{\int_{S\cap[r_{obs},\infty)}f_{\chi_d}(r)\,dr}
                  {\int_S f_{\chi_d}(r)\,dr},\qquad
f_{\chi_d}(r)=\frac{2^{1-d/2}}{\Gamma(d/2)}r^{d-1}e^{-r^2/2}.
\]

The set can have disconnected intervals. The existing native replay already
demonstrates full-history selection, exclusion and reentry at increasing radii.
A single lower truncation threshold therefore needs additional justification.

## Conditional Monte Carlo computation

The new diagnostic draws a fixed number \(N\) of independent chi radii and
retains the \(K\) whose reconstructed diffusion history matches the observation.
It returns

\[
p_{MC}=\frac{1+\sum_{j=1}^K 1\{r_j\ge r_{obs}\}}{K+1}.
\]

Conditional on \(A,u,M,K\), the observed radius and accepted radii are
exchangeable under the null. With continuous radii, the rank is uniform on
\(1,\ldots,K+1\). Consequently rejection probability conditional on \(K\)
is \(\lfloor\alpha(K+1)\rfloor/(K+1)\le\alpha\). The proposal count
and random seed are set independently of the observation. This is a finite
Monte Carlo test, and its rank resolution depends on the accepted count.
Zero acceptances give p=1. At alpha 0.01, at least 99 acceptances are required
for rejection to be possible.

An unclassified proposal stops the diagnostic with an unresolved result. A
later near-tie can coexist with an earlier strict violation of a required
merge inequality. When the native history differs and such a strict violation
exists, the proposal is excluded. The initial experiment treated every near-tie
as unresolved. Its output is retained, and the paired rerun changes only that
exclusion rule. The regression test failed before the rule was added and passed
afterward.

The rank argument assumes correct event evaluation. The implementation uses
native floating arithmetic and retains the existing replay tolerance. It has
no exact-real membership certificate. This differs from the existing interval
diagnostic, whose certified mass bounds remain wide. Numerical point replay
and symbolic density checks alone cannot certify all diffusion event boundaries.

## Wolfram and independent probability control

`conditional-diffusion-wolfram-20261104.json` contains the submitted Wolfram
Language and returned results. Wolfram integrated the chi density to one and
verified the derivative of its incomplete-gamma survival function. For df2,
\(S=[1,2]\cup[3,\infty)\) and observed radius 3.2, it returned conditional
p=0.01239056336490902, compared with nominal p=0.005976022895005943.
SciPy independently reproduced this value and the binomial mixtures of the
finite-rank rejection probabilities.

A separate exact-region control used 20,000 independent observations from this
truncated chi distribution and 199 proposals each. It rejected 59, 863 and
1,925 times at alpha 0.01, 0.05 and 0.10. Each Clopper-Pearson interval contains
the corresponding analytical finite-rank rejection probability. The region
probability is 0.48230437301426304. Reproduce with NumPy seed 2026110444,
\(a=e^{-1/2},b=e^{-2},c=e^{-9/2},q=a-b+c\). Draw 20,000 uniforms \(u\) on
\([0,q]\), set survival \(a-u\) when \(u<a-b\) and \(q-u\) otherwise, and
take radius \(\sqrt{-2\log(\mathrm{survival})}\). Then draw a 20,000-by-199
chi-square-2 array, take its square roots, retain entries in S and form the ranks.
The counts and intervals are in
`raw/assets/conditional-diffusion-20261104/model-control.json`.

Wolfram verifies the stated symbolic calculations. The exchangeability
derivation above supplies the statistical argument. The diffusion event
implementation and its numerical behavior are tested separately.

## Native diffusion experiment

The fixed panel contains 60 homogeneous null datasets and 20 paired datasets
at each of two signal levels. Each has nine rows and three features with known
identity covariance. Four rows have first-coordinate mean -delta and five have
mean +delta, for delta 0, 1 and 2. Data seeds are 2026110401 plus replicate.
Proposal seeds are 2026120401 plus replicate. Every run requests 999 proposals.
Alpha levels are 0.01 and 0.05. All cases, including unresolved cases, stay in
the panel denominator. These are root detection measurements. They supply no
descendant traversal or final cluster-recovery validation.

The executable protocol is
`raw/inbox/conditional-diffusion-experiment-20261104.py.txt`.
Run it with `uv run --no-sync --offline python`, followed by the script path and
one of `null`, `moderate` or `strong`. It creates one output directory per
scenario and refuses to overwrite its existing replay capture. Results retain
accepted radii, seeds, source hashes, root statistics and package versions.

All 100 final cases completed. At alpha 0.01 the observed outcomes were

| Scenario | Datasets | Nominal chi rejections | Conditional MC rejections | Median accepted radii | Zero-acceptance cases |
| --- | ---: | ---: | ---: | ---: | ---: |
| Homogeneous Gaussian null | 60 | 20 | 0 | 19.5 | 3 |
| Signal, delta 1 | 20 | 17 | 0 | 2 | 5 |
| Signal, delta 2 | 20 | 20 | 0 | 0 | 20 |

At alpha 0.05 the nominal counts were 36/60, 20/20 and 20/20, while the
conditional counts remained zero. Only six null cases, and no signal cases,
accepted the 99 radii needed to attain p<=0.01 with this rank statistic.
Every strong-signal p=1 arose from zero acceptances. These values reflect
insufficient conditional Monte Carlo resolution. They give no evidence of
equal group means. The finite panel therefore exposes an ineffective sampling
strategy for these strongly selected histories. It does not measure the power
of an accurately evaluated continuous conditional tail.

The 95% Clopper-Pearson interval for zero of 60 null rejections is
[0, 0.0596295]. This sample does not establish a 0.01 false-rejection bound.
The exact-region probability control and exchangeability argument establish
different facts from the finite native-diffusion panel.

Median elapsed time per nine-row root was 4.44 to 4.70 seconds with 999
proposals. The full panel evaluated 99,900 reconstructed diffusion trees.
The initial replay rule left 29 null, 10 moderate and 7 strong cases unresolved.
The strict-exclusion rerun resolved all of them. Every initially complete
case retained exactly the same rank p-value with the paired data and seeds.
`raw/assets/conditional-diffusion-20261104/summary.json` records the
independently recomputed ranks, paired checks, source hashes and intervals.

The current implementation has eleven focused tests. The strict-exclusion
case includes a native-data regression as well as a controlled margin case.

## Scope for the cancer cell-line data

The actual cell-line input is transformed deposited CPM with variable-gene
selection and fitted PCA. The known independent Gaussian covariance model used
here has not been established for those measurements. A valid application
needs a justified observation and covariance model, treatment of upstream
preprocessing, and a defined selection event for every tested descendant and
traversal decision. This experiment tests the root law and its computation.
Conditioning only on the observed root partition would use a larger selection
event. Computing that event and evaluating its power remain separate work.

## Sources

- Gao, Bien and Witten, *Selective Inference for Hierarchical Clustering*,
  https://doi.org/10.1080/01621459.2022.2116331.
- `benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_replay.py`
  owns the native replay reused here.
- `benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_conditional_mc.py`
  implements the conditional rank statistic.
- `tests/statistics/67_test_diffusion_root_conditional_mc.py` checks rank
  accounting, empty acceptance, ambiguous membership, strict exclusion,
  reproducibility and input preservation.
- `raw/inbox/analytic-selected-root-derivation-20261102.md` and
  `raw/inbox/diffusion-self-distance-and-region-tightening-20261102.md` contain
  the existing restricted derivations and interval certification evidence.
