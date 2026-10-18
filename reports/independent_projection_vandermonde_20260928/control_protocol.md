# Independent-data projection control: frozen design

This companion to `reports/known_branch_skew_20260928/protocol.md` changes no
production method. Freeze this protocol, configuration and imported project
source hashes before the main run. Main: 1,000 replicates in each of the same
five patterns and two covariance profiles, 32 samples per leaf group, eight
dimensions, known covariance, zero ridge, fixed R/A/B, alpha 0.01. Population
means, Gaussian sampling and node truth are imported from the original diagnostic.

Evaluation draws exactly reproduce the original profile/pattern/replicate order
and `SeedSequence(20260928).spawn(10000)` children. Basis training and rank
training are independent draws from the same population, with distinct spawn
keys `(17021983,1,evaluation_key)` and `(17021983,2,evaluation_key)`. Record all
three input hashes and keys. These are simulation controls with extra data,
not a recommendation to halve the production sample or a trained calibration.

Cross basis source {evaluation PCA, independent basis-training PCA, prespecified
coordinate prefix} with rank source {evaluation MP, independent rank-training
MP, prespecified 1, 2, 8}. Use production covariance whitening, eigendecomposition,
MP dimension estimator, tie-completing PCA projection and projected-Wald kernel.
The fixed coordinate prefix has no PCA eigenspectrum or eigenvalue ties, and
uses the identical unweighted squared norm with its exact coordinate rank.
Rank-training and basis-training use different samples: both are independent
of the evaluation contrast and of one another. Known covariance prevents
covariance estimation from becoming another factor.

Record MP raw count, requested rank and actual rank after complete numerical
tie groups. Under the declared Gaussian mean-null, z has law N(0,I), including
fixed mixtures whose child population averages agree. Independent basis and
rank yield a central chi-square law conditional on the entire training data,
with actual rank as degrees of freedom. Fixed coordinate/prespecified-rank
controls also do. Evaluation-dependent basis or rank generally does not have
this conditional guarantee. Prespecified full rank is basis invariant. A
prespecified partial PCA rank may expand through ties; its actual rank is used.

Compare basis sources at held-fixed rank source, and rank sources at held-fixed
basis source. Preserve matched outcomes and their interaction; changed selected
subspaces are not proof that one factor alone caused every observed difference.
Report canonical MP rank variability explicitly. A floor of one can make the
all-null canonical MP dimension constant even when its raw signal count varies;
then the canonical panel does not identify a dimension-selection effect there.

Prespecified supplementary stress: global_null only, both covariance profiles,
same evaluation/training data and all three basis sources, crossed with the two
MP rank sources, but set `mp_threshold_rows=64*actual_rows` and retain the actual
row count in `effective_independent_rows`. This deliberately lowers the threshold
to exercise variable ranks. It is a changed-threshold diagnostic, not canonical
MP, not evidence that production uses these rows, and not a proposed fix.
Never pool this stress with canonical MP or tune its multiplier after outcomes.

Retain error statuses/reasons and denominators. Summarize per named node/pattern/
profile/variant, at most one trial per independent dataset: null rejection counts
and pointwise Wilson 95% intervals, signal power separately, rank distributions,
statistic/df means, null-only skew and p-value quantiles. The central chi-square
df-mixture skew is a reference comparator, not an asserted law for
evaluation-selected projections. Save each matched pair's
rejections/statistic/ranks, discordant counts, mean paired rejection difference
and descriptive normal intervals (which may degenerate with no discordance).
Do not pool dependent R/A/B as independent trials. This finite Gaussian panel
does not establish production calibration, tree selection validity or causality.

A ten-replicate/cell pilot with base seed 20260929 checks execution only and is
excluded from evaluation. Main source/config hashes must match the freeze before
execution. Existing original diagnostic source and outputs remain unchanged.

Before the main freeze, source verification was strengthened to recheck all
frozen hashes after computation; completion records any changed sources and
raises if present. The conditional-law flag for the prespecified full rank
also requires actual rank to equal feature dimension. These checks do not
alter default pilot statistics or the prespecified main design.
