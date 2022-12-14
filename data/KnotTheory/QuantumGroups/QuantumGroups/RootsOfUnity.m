(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



BeginPackage["QuantumGroups`RootsOfUnity`",{"QuantumGroups`","QuantumGroups`RootSystems`"}];


AlcoveDefiningRoot;WeightInAlcoveQ;AlcoveWeights;AlcoveWeightsInLattice;AlcoveRoots;LevelFromRoot;RootFromLevel;


Begin["`Private`"];


AlcoveDefiningRoot[\[CapitalGamma]_,l_]:=AlcoveDefiningRoot[\[CapitalGamma],l]=With[{lp=If[EvenQ[l],l/2,l]},
If[Divisible[lp,LacingNumber[\[CapitalGamma]]],
LongDominantRoots[\[CapitalGamma]][[1]],
ShortDominantRoots[\[CapitalGamma]][[1]]]
]


WeightInAlcoveQ[\[CapitalGamma]_,l_Integer][\[Lambda]:{___Integer}]:=(And@@(NonNegative/@\[Lambda]))\[And](KillingForm[\[CapitalGamma]][AlcoveDefiningRoot[\[CapitalGamma],l],\[Lambda]+\[Rho][\[CapitalGamma]]]<If[EvenQ[l],l/2,l])


AlcoveWeights[\[CapitalGamma]_,l_]:=AlcoveWeights[\[CapitalGamma],l]=Module[{ar=AlcoveDefiningRoot[\[CapitalGamma],l],\[Lambda]=ZeroVector[Rank[\[CapitalGamma]]],p},
Reap[
While[
Sow[\[Lambda]];
\[Lambda]+=UnitVector[Rank[\[CapitalGamma]],1];
While[\[Not]ZeroVectorQ[\[Lambda]]\[And]\[Not]WeightInAlcoveQ[\[CapitalGamma],l][\[Lambda]],
(* odometer *)
p=Position[\[Lambda],x_/;x>0][[1,1]];
\[Lambda][[p]]=0;
If[p<Rank[\[CapitalGamma]],++\[Lambda][[p+1]]];
];
\[Not]ZeroVectorQ[\[Lambda]]
]
][[2,1]]
]


AlcoveWeightsInLattice[\[CapitalGamma]_,l_,All]:=AlcoveWeights[\[CapitalGamma],l]


AlcoveWeightsInLattice[\[CapitalGamma]_,l_,lattice_]:=AlcoveWeightsInLattice[\[CapitalGamma],l,lattice]=Cases[AlcoveWeights[\[CapitalGamma],l],\[Lambda]_/;WeightInLatticeQ[\[CapitalGamma],\[Lambda],lattice]]


AlcoveRoots[\[CapitalGamma]_,l_]:=AlcoveRoots[\[CapitalGamma],l]=Cases[AlcoveWeights[\[CapitalGamma],l],\[Lambda]_/;RootWeightQ[\[CapitalGamma],\[Lambda]]]


LevelFromRoot[\[CapitalGamma]_,l_]:=l/(2LacingNumber[\[CapitalGamma]])-DualCoxeterNumber[\[CapitalGamma]]


RootFromLevel[\[CapitalGamma]_,k_]:=2LacingNumber[\[CapitalGamma]](k+DualCoxeterNumber[\[CapitalGamma]])


AlcoveRoots[\[CapitalGamma]_,l_]:=AlcoveRoots[\[CapitalGamma],l]=Cases[AlcoveWeights[\[CapitalGamma],l],\[Lambda]_/;RootWeightQ[\[CapitalGamma],\[Lambda]]]


End[];


EndPackage[];
