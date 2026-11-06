# Analytic selected-root derivation — 2026-11-02

User chose analytic selected-region inference before production integration.
Three Sol agents investigated selection, covariance, and sibling mechanisms.

## Primary references inspected

- Gao, Bien and Witten, *Selective Inference for Hierarchical Clustering*,
  DOI 10.1080/01621459.2022.2116331; abstract inspected at
  https://pubmed.ncbi.nlm.nih.gov/38660582/ . It establishes the need to account
  for clustering selection and develops selective cluster-mean tests.
- Yun and Barber, *Selective inference for clustering with unknown variance*,
  https://arxiv.org/html/2301.12999v3 . Inspected the known-variance conditional
  radius construction and the unknown-variance model. Our scalar global-null
  construction below is a restricted implementation, not a new theorem claim.

## Independently derived model

For known sigma, independent Gaussian rows have covariance sigma² I. For the
observed root contrast v, e=v/||v||, z=X' e/sigma, u=z/||z|| and
A=(I-ee')X. Condition on A, u and the complete observed average-linkage merge
sequence. Along X(r)=A+sigma e r u', every squared-distance average-linkage
comparison is a quadratic inequality in r. Their intersection is a finite
interval union S. The conditional chi_d tail over S gives the root p-value.
The contrast null is e'E[X]=0. A common mean is sufficient but not necessary.

For unknown sigma and scalar rows, require the stronger global iid N(mu,sigma²)
null. Let a=e'x, w=x-mean(x)1-ea, nu=n-2, eta=w/||w||, s=sign(a),
R=||x-mean(x)1||, h=|a|/||w||. Condition on mean, R, eta, s and the merge
sequence. The path is mean1+R(eta+s e h)/sqrt(1+h²). Common squared-distance
factor R²/(1+h²) cancels from every linkage comparison; use the same quadratic
interval solver. nu h² has F(1,nu) law truncated to that region. Nonconstant
within-root mean structure would invalidate this central F law.

## Wolfram check

Input:
```wolfram
ClearAll[h,nu]; FullSimplify[2*nu*h*PDF[FRatioDistribution[1,nu],nu*h^2],h>0&&nu>0]
```
Returned:
```text
(2*((1 + h^2)^(-1))^((1 + nu)/2))/Beta[1/2, nu/2]
```

## Independent executable verification

Run from the repository with existing `uv run --no-sync python`:

```python
import numpy as np
from scipy.stats import kstest
from benchmarks.diagnostics.calibration.selected.hierarchy.gaussian_root_region import selected_gaussian_root_test
from benchmarks.diagnostics.calibration.selected.hierarchy.gaussian_root_unknown_variance import selected_scalar_gaussian_root_test
rng = np.random.default_rng(20261103)
for name, fn in [('known_sigma', lambda x: selected_gaussian_root_test(x[:, None], sigma=1)),
                 ('unknown_sigma', selected_scalar_gaussian_root_test)]:
    ps = np.array([fn(rng.normal(size=4)).p_value for _ in range(300)])
    print(name, len(ps), int((ps < .05).sum()), kstest(ps, 'uniform').statistic)
x = np.array([-100.07, -99.94, 99.91, 100.13])
r = selected_scalar_gaussian_root_test(x)
assert r.residual_degrees_of_freedom == 2 and r.p_value < .001
print(r.p_value, r.statistic, r.intervals)
```

Observed output:
```text
known_sigma 300 20 0.044752619581589714
unknown_sigma 300 13 0.03266094299608624
3.1601693694910593e-06 2450842.304747382 ((1.7719463435936105, inf),)
```

P-values use deterministic analytic interval construction and numerical tail
integration. Simulation is verification only, never a calibration fallback.
These small panels cannot certify all parameter settings or floating boundaries.
The symmetric earlier linspace examples have equal-distance merge ties and are
outside the new diagnostic's declared no-ties domain.
