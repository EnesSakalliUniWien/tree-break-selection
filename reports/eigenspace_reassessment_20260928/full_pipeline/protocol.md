# Frozen complete-pipeline reassessment protocol

Frozen before the evaluation-seed results, 2026-09-28. No tuning or production
change is permitted from these outcomes within this study. A separate pilot
namespace checks execution and audit extraction; its outcomes are excluded.

## Data and independent units

Six families, with one-population null and two-population signal datasets:

| Family | Samples and features | Population laws |
|---|---|---|
| Bernoulli iid | n=64, p=24 | Null p=.5; signal .2 vs .8 in first 12 features, remaining .5; equal population sizes |
| Bernoulli dependent | n=64, p=24 | Eight independent Bernoulli features copied three times each; signal first four underlying probabilities .2 vs .8, remaining .5; population sizes 48/16 |
| Categorical | n=64, eight 3-category blocks | Null uniform; signal (.65,.20,.15) vs (.15,.20,.65) in every block; equal sizes |
| Gaussian identity | n=64, p=8 | Common identity covariance; null mean zero; signal means -1 vs +1 in first two coordinates, others zero; equal sizes |
| Gaussian AR1 | n=64, p=8 | Covariance 0.6^abs(i-j); signal means as above; sizes 48/16 |
| Gaussian high dimensional | n=40, p=64 | Identity covariance; signal -1 vs +1 in first eight coordinates; equal sizes |

There are 500 null and 200 signal datasets per family: 4,200 independent
datasets. All methods and alpha settings share each exact input. Each dataset
seed is generated from SeedSequence([20260928, phase, family_index, signal,
replicate]), where evaluation phase=17 and pilot phase=91. Disjoint seeds
separate evaluation from pilot and prior rank experiments. A seeded random
permutation jointly shuffles each matrix and its population labels before any
method, so input-order tie breaking is not aligned with population membership.
Population profiles are fixed independently of tree inference. Input hashes
allow matched-case verification. Dependent Bernoulli deliberately stresses the
model's independent-feature covariance; it is not an iid null-model guarantee.

## Methods and fitting

Use the current registered `tbs`, `tbs_complete`, `tbs_neighbor_joining`,
`tbs_diffusion` and `tbs_diffusion_adaptive_nnls` parameters. For native Gaussian
inputs, direct distance-tree methods use the established standardized-Euclidean
metric. The two registered Hamming diffusion methods remain explicitly skipped
on continuous inputs. Do not replace them with unlabeled Euclidean variants.

Every method refits geometry, tree, branch lengths where configured, covariance,
PCA, complete-group rank, edge selection, internal calibration, sibling BH and
traversal. Adaptive NNLS uses original features for its branch-length target.
The actual executed gate configuration is saved per run. The primary sibling
alpha is .01 and edge alpha is .001. Prespecified sibling-alpha sensitivity
.005/.05 applies to iid Bernoulli null and signal data for `tbs` and adaptive
NNLS, using the same independent datasets. No edge-alpha sensitivity is claimed:
the diffusion dispatcher currently does not forward that argument.

The generic methods use `adaptive_projection_dimension_energy_fraction=None`.
The separate isolated Wald experiment assesses the .9 adaptive rule. This
full-pipeline study is not a test of scRNA application defaults or every method
variant.

## Truth, endpoints and uncertainty

The primary global-null endpoint is Pr(final K>1), compared diagnostically with
the supplied sibling alpha. This is a whole-pipeline criterion, not a claimed
theorem of the existing per-depth procedure. Report coverage and a conditional
Wilson 95% interval among completed datasets. Successful method outputs remain
in this primary denominator even if subsequent diagnostic extraction fails;
split/root/support metrics have a separate audited denominator. Also report the observed event
count per attempted dataset and bounds that leave unavailable results unknown.
Never count unsupported/error results as correct one-cluster answers.

For a sibling contrast, H0 is equality of generator population-mixture
parameters, determined by exact integer cross-products of population counts
and child sizes. It does not use observed means, purity cutoffs or edge tests.
Actual false rejections are reached traversal decisions equal to `split` at
generator-null nodes. Count `pass_through` separately. Signal outputs report
ARI, K, true/false split counts and dataset-level false-discovery proportions.
These signal FDPs describe generator-null contrast rejections; extra final
clusters can also arise through pass-through and are a different endpoint.

Record all binary internal nodes, raw and adjusted p-values, degrees of freedom,
projection dimension, support status, calibration base weight, generator truth
and traversal use. Count nonnull calibration records and their base weight;
base weight is not the focal-context kernel weight. Root-null tails use at most
one root per independent dataset, conditional on a completed run and a finite
root-null p-value; eligibility counts remain explicit. Thresholds are
.001/.005/.01/.05. Do not pool dependent nodes into binomial sample sizes.

One-sided exact binomial tests diagnose excess global-null partitions and
root-tail mass. Apply Holm adjustment separately across all 34 prespecified
global-null partition cells and all 272 root-tail cells, including unavailable
cells with p=1 for multiplicity while retaining NA in displayed results.
Wilson intervals are pointwise descriptive intervals, not simultaneous bounds.
No favorable nonrejection certifies calibration. At 500 replicates, Monte Carlo
SE near p=.01 is about .00445; rare eligible contexts have worse precision.

## Execution and reporting

Use the repository virtual environment, six independent worker processes,
`TBS_N_JOBS=1` and one BLAS thread per process. Keep source hashes fixed throughout evaluation.
Retain errors with reasons and preserve the prior reports. Replay traversal
using the already fitted gate bundle only to expose its full trace, and assert
that replay labels match dispatcher labels. No parameter selection, correction,
commit, push or dependency installation is part of this study.
