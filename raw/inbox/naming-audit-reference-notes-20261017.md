# Naming audit reference notes

Captured on 2026-10-17 for the naming-consistency audit. These are short
paraphrased reference notes, not copies of external documents.

## Multiple-testing terminology

Source: [R stats: Adjust P-values for Multiple Comparisons](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).

The official documentation distinguishes family-wise error control, including
Holm and Bonferroni, from false discovery rate (FDR) control by
Benjamini-Hochberg (BH) and Benjamini-Yekutieli (BY). FDR is the expected
proportion of false discoveries among rejected hypotheses. The page cites
Holm (1979), Benjamini and Hochberg (1995), and Benjamini and Yekutieli (2001).
The original papers were not retrieved in full for this audit.

Repository implication: a dispatch accepting all four methods needs the
generic term `adjustment`, rather than `fdr`, in its common API. This naming
change does not establish the validity of the pipeline's selected p-values.

## Standardized Euclidean distance

Source: [SciPy: scipy.spatial.distance.seuclidean](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.seuclidean.html).

The documented squared standardized Euclidean distance is the sum over
coordinates of squared differences divided by coordinate variances. The
repository's NNLS target further divides this sum by the total input feature
count. Its name therefore needs `mean_squared_standardized_euclidean`.
Constant columns remain in that divisor under the existing implementation.

The source reviewed was the current online manual. Local numerical checks
used the existing SciPy 1.16.3 environment; no package was installed or
updated. The naming change leaves the target calculation unchanged.
