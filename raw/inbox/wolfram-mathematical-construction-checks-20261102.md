# Wolfram mathematical construction checks — 2026-11-02

Source: connected Wolfram plugin, Wolfram Language evaluator.
Scope: symbolic identities and counterexamples extracted from current local sources.
These evaluations do not execute Tree-Break Selection or validate statistical calibration.

The initial WolframContext call failed with an internal error. The evaluator calls below succeeded.
The second call contained one unevaluated distribution expression; the third call replaces that
check with TransformedDistribution and differentiates the MGF for the moment check.

## Evaluation 1

```wolfram
ClearAll[a,b,x,y,v,rho,tL,tR,tbar,w,c];
mu=(a*x+b*y)/(a+b);
edgeL=1/a-1/(a+b); edgeR=1/b-1/(a+b); sib=1/a+1/b;
ass=a>0&&b>0&&v>0&&Element[{x,y},Reals];
out=<|
"barycenterContrasts"->FullSimplify[{x-mu-b/(a+b)*(x-y),y-mu+a/(a+b)*(x-y)},ass],
"nestedCovariance"->FullSimplify[(b/(a+b))^2*(1/a+1/b)-(1/a-1/(a+b)),ass],
"standardizedLeftIdentity"->FullSimplify[(x-mu)/Sqrt[v*edgeL]-(x-y)/Sqrt[v*sib],ass],
"standardizedRightIdentity"->FullSimplify[(y-mu)/Sqrt[v*edgeR]+(x-y)/Sqrt[v*sib],ass],
"splitAction"->FullSimplify[a*(x-mu)^2+b*(y-mu)^2-a*b/(a+b)*(x-y)^2,ass],
"correlatedUnitVarianceProjectionVariance"->FullSimplify[({1,1}/Sqrt[2]).{{1,rho},{rho,1}}.({1,1}/Sqrt[2])],
"correlatedUnitVarianceEigenvalues"->Eigenvalues[{{1,rho},{rho,1}}],
"fixedChiSquareMGF"->FullSimplify[MomentGeneratingFunction[ChiSquareDistribution[3],t]],
"selectedTailExample"->FullSimplify[SurvivalFunction[ChiSquareDistribution[2],w]/SurvivalFunction[ChiSquareDistribution[2],c],w>=c&&c>0],
"naivePConditionalMeanExample"->Integrate[Exp[-w/2]*(Exp[-w/2]/2)/Exp[-c/2],{w,c,Infinity},Assumptions->c>0]
|>; out
```

Returned text:

```text
Out[1]= <|"barycenterContrasts" -> {0, 0}, "nestedCovariance" -> 0, "standardizedLeftIdentity" -> 0, "standardizedRightIdentity" -> 0, "splitAction" -> 0, "correlatedUnitVarianceProjectionVariance" -> 1 + rho, "correlatedUnitVarianceEigenvalues" -> {1 - rho, 1 + rho}, "fixedChiSquareMGF" -> (1 - 2*t)^(-3/2), "selectedTailExample" -> E^((c - w)/2), "naivePConditionalMeanExample" -> 1/(2*E^(c/2))|>
```

## Evaluation 2

```wolfram
ClearAll[q1,q2,a,b,x,t,rho,k];
p={q1,q2}; cat=DiagonalMatrix[p]-Outer[Times,p,p];
ass=q1>0&&q2>0&&q1+q2<1;
out=<|
"categoricalDropLastDeterminant"->Factor[Det[cat]],
"categoricalPositiveDefiniteInterior"->FullSimplify[{cat[[1,1]]>0,Det[cat]>0},ass],
"branchTimeLeftRatioSquared"->FullSimplify[(1+(x+y)/(2*t))/(1+x/t),t>0&&x>=0&&y>=0],
"correlatedQuadraticMGF"->FullSimplify[(1-2*(1-rho)*t)^(-1/2)*(1-2*(1+rho)*t)^(-1/2)],
"correlatedQuadraticMeanVariance"->{2,4*(1+rho^2)},
"adaptiveUnitRowCounterexample"->FullSimplify[({a,b}/Sqrt[a^2+b^2]).{a,b},Element[{a,b},Reals]&&a^2+b^2>0],
"equalMeanDifferentVariance"->{Mean[ChiSquareDistribution[2]],Mean[2*ChiSquareDistribution[1]]}
|>; out
```

Returned text:

```text
Out[1]= <|"categoricalDropLastDeterminant" -> -(q1*q2*(-1 + q1 + q2)), "categoricalPositiveDefiniteInterior" -> {True, True}, "branchTimeLeftRatioSquared" -> (2*t + x + y)/(2*(t + x)), "correlatedQuadraticMGF" -> 1/(Sqrt[1 + 2*(-1 + rho)*t]*Sqrt[1 - 2*(1 + rho)*t]), "correlatedQuadraticMeanVariance" -> {2, 4*(1 + rho^2)}, "adaptiveUnitRowCounterexample" -> Sqrt[a^2 + b^2], "equalMeanDifferentVariance" -> {2, Mean[2*ChiSquareDistribution[1]]}|>
```

## Evaluation 3

```wolfram
ClearAll[q1,q2,x,y,t,rho,s,a,b];
p={q1,q2}; cat=DiagonalMatrix[p]-Outer[Times,p,p];
mgf=(1-2*(1-rho)*s)^(-1/2)*(1-2*(1+rho)*s)^(-1/2);
scaled=TransformedDistribution[2*x,x\[Distributed]ChiSquareDistribution[1]];
<|
"categoricalDropLastDeterminant"->Factor[Det[cat]],
"categoricalPositiveDefiniteInterior"->FullSimplify[{cat[[1,1]]>0,Det[cat]>0},q1>0&&q2>0&&q1+q2<1],
"branchTimeLeftRatioSquared"->FullSimplify[(1+(x+y)/(2*t))/(1+x/t),t>0&&x>=0&&y>=0],
"correlatedQuadraticMeanVarianceFromMGF"->FullSimplify[{D[mgf,s]/.s->0,(D[mgf,{s,2}]/.s->0)-(D[mgf,s]/.s->0)^2}],
"adaptiveUnitRowCounterexample"->FullSimplify[({a,b}/Sqrt[a^2+b^2]).{a,b},Element[{a,b},Reals]&&a^2+b^2>0],
"equalMeanDifferentVariance"->{{Mean[ChiSquareDistribution[2]],Variance[ChiSquareDistribution[2]]},{Mean[scaled],Variance[scaled]}}
|>
```

Returned text:

```text
Out[1]= <|"categoricalDropLastDeterminant" -> -(q1*q2*(-1 + q1 + q2)), "categoricalPositiveDefiniteInterior" -> {True, True}, "branchTimeLeftRatioSquared" -> (2*t + x + y)/(2*(t + x)), "correlatedQuadraticMeanVarianceFromMGF" -> {2, 4*(1 + rho^2)}, "adaptiveUnitRowCounterexample" -> Sqrt[a^2 + b^2], "equalMeanDifferentVariance" -> {{2, 4}, {2, 8}}|>
```

