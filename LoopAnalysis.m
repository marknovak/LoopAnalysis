(* ::Package:: *)


(*  Paste in the community matrix here.  This should be the output from converting the database-exported matrix in Mathematica format using the R-script *)


SetAttributes[ZD,Listable];
ZD[x_,y_]:=If[y==0,1,x/y];
Minor[m_List?MatrixQ,{i_Integer,j_Integer}]:=Abs[Drop[Transpose[Drop[Transpose[A],{j}]],{i}]]
Permanent[m_List]:=With[{v=Array[x,Length[m]]},Coefficient[Times@@(m.v),Times@@v]]


A=Import["A.csv"]


n=Length[A];



							(* Calculate the adjoint matrix - adjA *) 
Export["MathematicaOutComm.csv",A,"CSV"];
adjA=Inverse[-A]*Det[-A];
Export["Adj.csv",adjA,"CSV"];


							(* Calculate Absolute Feedback Matrix - T *)


										(* Method 2 *)


T=Outer[Permanent[Minor[Abs[A],{##}]]&,Sequence@@Range/@Dimensions[A],1];
Export["AbsFeed.csv",T,"CSV"];


					(* Calculate the weighted predictions matrix -W *)


W=N[ZD[Abs[adjA],T]];
Export["WPred.csv",W,"CSV"];
