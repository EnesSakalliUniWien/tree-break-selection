# Wolfram checks for fixed-partition Gaussian attribution

Executed 2026-11-02 using the explicitly requested Wolfram plugin. The context lookup did not supply a useful selected-tree derivation; the Wolfram Language evaluator performed the checks below. These checks verify stated algebra and distribution calculations, not the diffusion selection event or production calibration.

## Model and scope

For iid Gaussian rows with covariance I and a partition independent of testing data, let s=1/nL+1/nR, z=(meanL-meanR)/sqrt(s), and Q=z'z. Conditional on the entire independently selected tree, z is N(0,I_p), hence Q is chi-square_p. In the selected-data arm the partition depends on the testing observations, so that reference is not justified by this calculation.

The production covariance audit found three regimes for n=80,p=24: dense pooled covariance when both children have at least 25 leaves; diagonal pooled covariance when the size guard passes but a child has fewer than 25 leaves; dense parent-total fallback when the size guard fails. The following finite-sample laws apply to fixed/independent partitions and ideal unregularized matrices. Production uses tiny numerical ridge terms; measured differences are retained separately.

## Check 1: contrast equivalence, known covariance and pooled reference

For a binary root the full-row child-parent statistic equals the sibling statistic under the same covariance. This identity does not identify the split-row production edge statistic with the full-row diagnostic.

For ordinary dense within-child pooled S with nu=n-2, T2=z' S^-1 z obeys ((nu-p+1)/(nu*p))*T2 ~ F(p,nu-p+1) under the independent-data Gaussian null. This law is not asserted for selected partitions, diagonal estimates or fallback covariance.

```wolfram
Clear[l,r,n,d,nu,t2,q]; n=l+r;
edgeEquality=FullSimplify[((r/n)^2)/(1/l-1/n)==1/(1/l+1/r),Assumptions->l>0&&r>0];
mgfCheck=FullSimplify[MomentGeneratingFunction[ChiSquareDistribution[24],t]==(1-2t)^(-12),Assumptions->t<1/2];
totalMap=FullSimplify[(n-1)*t2/(n-2+t2)];
threshold=InverseCDF[ChiSquareDistribution[24],99/100];
pooledFalsePositive=N[SurvivalFunction[FRatioDistribution[24,55],threshold*55/(78*24)],12];
{edgeEquality,mgfCheck,totalMap,N[threshold,12],pooledFalsePositive,78*24/53}
```

Evaluator response:

```text
Out[1]= {True, True, ((-1 + l + r)*t2)/(-2 + l + r + t2), 42.97982013935163628768655608158866243684`12., 0.23427010839488665222928250107107680099`12., 1872/53}
```

At n=80,p=24, comparing the ideal dense-pooled statistic to chi-square 24's .01 cutoff yields rejection probability .23427010839488665 under the independent null. Its mean is 1872/53=35.3207547 rather than 24. This is an analytic explanation of the covariance/reference mismatch, not an empirical error-rate estimate or production repair.

## Check 2: parent-total covariance and diagonal pooled moments

Write S_total=(nu*S_pool+z*z')/(nu+1). The quadratic statistic is Q_total=(nu+1)*T2/(nu+T2). The evaluator verifies its rank-one inverse identity symbolically on a general symmetric 2x 2 matrix; the same general identity follows from the Sherman-Morrison formula. For p 24,nu 78 it verifies the change of variables from F24,55 to Q_total/79 ~ Beta(12,55/2). This shows why parent-total and pooled covariance cannot use the same finite-sample reference.

The independent, identity-covariance diagonal-pooled case has Q_diag=sum of 24 independent F(1,78) variables. Its mean and variance are checked below. It is not the dense Hotelling F law.

```wolfram
Clear[a,b,c,x,y,nu,v]; mat={{a,b},{b,c}}; z={x,y}; q=z.Inverse[mat].z; total=(nu mat+Outer[Times,z,z])/(nu+1); identity=FullSimplify[z.Inverse[total].z-(nu+1)q/(nu+q)];
betaChange=FullSimplify[PDF[FRatioDistribution[24,55],55v/(24(1-v))]*55/(24(1-v)^2)==PDF[BetaDistribution[12,55/2],v],Assumptions->0<v<1];
cut=InverseCDF[ChiSquareDistribution[24],99/100];
{identity,betaChange,N[SurvivalFunction[BetaDistribution[12,55/2],cut/79],12],24 Mean[FRatioDistribution[1,78]],24 Variance[FRatioDistribution[1,78]]}
```

Evaluator response:

```text
Out[1]= {0, True, 0.00105855085946097943476419614048236605`12., 468/19, 702702/13357}
```

The ideal parent-total rejection probability at the same chi-square 24 cutoff is .0010585508594609794. The diagonal-pooled mean is 468/19 and variance 702702/13357. No distribution check here calibrates data-dependent projection, empirical-null fitting, or final traversal.

