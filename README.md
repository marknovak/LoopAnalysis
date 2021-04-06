# Loop Analysis
Loop analysis using Mathematica and R

Loop analysis is an approach for calculating how press perturbations to a community might affect each species’ equilibrium abundance.  In simplified terms, it calculates the net effect of all the direct and indirect interactions between species given that only the presence and sign (+/-) of each pairwise interaction is known.

The following Mathematica commands calculate the classical adjoint (adjA, i.e. the net effects, a.k.a. the adjugate matrix, or the transpose of the cofactor matrix), as well as [Dambacher et al.’s (2002)](https://doi.org/10.1890/0012-9658(2002)083[1372:ROCSIA]2.0.CO;2) absolute feedback (T), and weighted-predictions (W) matrices given a so-called qualitatively-specified Community matrix.  (Mathematica seems to use computer resources a bit more efficiently than Maple — for which Dambacher et al. provided code — especially for larger networks.)

The following was used in the analyses of [Novak, Wootton, Doak, Emmerson, Estes & Tinker (2011). Predicting community responses to perturbations in the face of imperfect knowledge and network complexity. Ecology, 92, 836-846](https://doi.org/10.1890/10-1354.1).

Download the file `LoopAnalysis.m` which contains the following:

(*Specify a community matrix*)
```
A = {{-1, 1, 1, 1},{-1, -1, 1, 1},{-1, -1, -1, 1},{-1, -1, -1, -1}}
```

(*Or, import community matrix from file in working directory*)
```
Directory[]
A=Import[“A.csv”];
A//MatrixForm
```

(*Create functions to compute matrix minors and matrix permanent*)
```
SetAttributes[ZD, Listable]; ZD[x_, y_]:= If[y == 0, 1, x/y];
Minor[m_List?MatrixQ,{i_Integer, j_Integer}]:=Abs[Drop[Transpose[Drop[Transpose[A],{j}]], {i}]]
Permanent[m_List]:= With[{v = Array[x, Length[m]]},Coefficient[Times @@ (m.v), Times @@ v]]
```

(*Calculate adjoint, absolute feedback, and weighted-predictions matrices *)
```
n = Length[A];
adjA = Inverse[-A]*Det[-A];
adjA//MatrixForm
T = Outer[Permanent[Minor[Abs[A], {##}]] &, Sequence @@ Range /@ Dimensions[A], 1];
T//MatrixForm
W = N[ZD[Abs[adjA], T]];
W//MatrixForm
Export["Adjoint.csv",adjA,"CSV"];
Export["AbsFeed.csv", T,"CSV"];
Export["WPred.csv",W,"CSV"];
```

To source these Mathematica commands from within R, place the Mathematica Package (`LoopAnalysis.m`) into a working directory.  In _R_, export your qualitatively-specified community matrix (*A*) into the same working directly...
```
fout1<- file(paste('A.csv'),"w");
write.table(A,file=fout1,sep=',',row.names=F,col.names=F,quote=F)
close(fout1)
```
...and run one of the following lines (depending on your platform)...

```
shell('math < LoopAnalysis.m',mustWork=T) # in Windows
system('Mathkernel < LoopAnalysis.m',wait=T) # on Mac
system('math < LoopAnalysis.m',wait=T) # on server
```
...then import the results back into R and continue with your analyses.
