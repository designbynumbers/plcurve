#!/usr/local/bin/WolframScript -script

(* Classify KnotTheory format pd-codes from a file *)

file  = ToString[$ScriptCommandLine[[2]]];

PrependTo[$Path, "../data"];
<< KnotTheory`
KnotsWithInvariants = ReadList["../data/invariantlist_old.txt"];

ClassifyKnot[PD_] := Module[{Invariants, PossibleKnots},
  Invariants = {Expand[Kauffman[PD][q, y]], KnotSignature[PD]};
  PossibleKnots = 
   Select[KnotsWithInvariants, Invariants == #[[3]] &];
   Table[PossibleKnots[[i]][[1]], {i, 1, Length[PossibleKnots]}]
  ];

data = 
  ToExpression[
   StringSplit[
     Import[$ScriptCommandLine[[2]]], "\n"][[2 ;; ;; 2]]];
  
dataRaw = ClassifyKnot /@ data;

FixSpace[s_] := 
  StringReplace[s, "," ~~ x : DigitCharacter :>  ", " <> ToString[x]];
ConnectSum[s_] := 
  If[TrueQ[Head[s[[1]]] == List], StringJoin[Riffle[#, "#"]] & /@ s, s];

ClassifyPDList[L_] := 
  Length /@ GroupBy[FixSpace /@ ConnectSum /@ L, First];
  
$ScriptCommandLine[[3]] = ClassifyPDList[dataRaw];

Save[$ScriptCommandLine[[4]],ScriptCommandLine[[3]]];


