(* ::Package:: *)

(* ::Title:: *)
(*Whitten Operations in KnotTheory*)


(* ::Subtitle:: *)
(*This package implements the action of the various Whitten group elements on KnotTheory's PD codes, and includes code for automatically determining some facts about the Whitten symmetry groups of knots and links.*)


BeginPackage["whitten`"]
Needs["KnotTheory`"]
Needs["Combinatorica`"]
Needs["PlotLegends`"]


(* ::Section:: *)
(*Utility Functions for Managing PD Codes*)


(* ::Text:: *)
(*This section of the library deals with managing PD codes. We often need to figure out (for instance) what the orientation of an arc is at a crossing, who is the next arc, and so forth. These are some helper functions to get those jobs done. Our operating principle is that we're going to be able to deduce a lot of information about a standard PD code from context, but it's convenient to expand this data into a new type called EPD codes in order to do the framing code.*)
(**)
(*In an ExpandedPD, each arc is a 3-tuple in the form*)
(**)
(*{x,y,{tc,hc}}*)
(**)
(*where x is the Primary index (inheirited from original PD), y is the Subscript index (used to subdivide arcs), tc is the Tail crossing, and hc is the Head crossing.*)
(**)
(*Each of these functions has a _PD version and an _EPD version (if the same code won't work for both):*)
(**)
(*"CrossingsForArc" is guaranteed to order the crossings so that the arc is going OUT of the first crossing and IN to the second crossing. Since we will later decorate our arc numbers by turning them into tuples, we take some care to use the Order function (which correctly orders lists) instead of just "<" when comparing arc numbers. We can also get to these with the convenience functions "ArcHead" and "ArcTail".*)
(**)
(*It is easy to get a list of all the arc in a PD code, and this is done by AllArcs[ThisPD]. To get the arcs by component number requires a scheme for deciding which loop of arcs corresponds to which component, and that depends currently on the way Skeleton partitions a PD code into components. Unfortunately, Skeleton is rather weird in its choices sometimes, and we need the component numbers to be a stable and predictable function of the arc numbers. So our rewrite of AllArcs[component_,ThisPD_] orders the components so that if i > j, then all arcs of component i are > all arcs of component j.*)
(**)
(*"OrientationAtCrossing" returns the strings "In" or "Out", giving the orientation of the arc at the crossing.*)
(**)
(*"OrientationsForArc" is a helper function for debugging, giving the two orientations of the arc at each incident crossing.*)
(**)
(*"ArcToLeft" returns the number (or tuple) for the arc to the left of the given arc at the given crossing.*)
(**)
(*"PDtoList" and "ListToPD" fix the heads in PD objects (or nested lists of depth two) to go back and forth to crossings*)
(**)
(*"ListToPD" turns a list of (ordered) quadruples and returns the corresponding KnotTheory PD code.*)


PDtoEPD[ThisPD_PD] := Module[{rules,arcs,i,newEPD},
	arcs = AllArcs[ThisPD];
	rules = Table[arcs[[i]] -> {arcs[[i]],1,CrossingsForArc[arcs[[i]],ThisPD]},{i,1,Length[arcs]}];
	newEPD = ThisPD /. rules;
	newEPD[[0]] = EPD;
	newEPD
];

EPDtoPD[ThisEPD_EPD] := Module[{arcs,rules,allarcs,newPD},
	arcs = Sort[AllArcs[ThisEPD]]; 
	rules = Table[arcs[[i]] -> i,{i,1,Length[arcs]}]; 
	newPD = ThisEPD /. rules;
	newPD[[0]] = PD;
	newPD
];

ToList[X_] := Module[{Y}, Y = X; Y[[0]] = List; Y]

PDtoList[ThisPD_] := Module[{L},
	L = ThisPD;
	L[[0]] = List;
	ToList /@ L
];

ListToPD[L_] := Module[{MyQuads,MyPD},
	MyPD="PD[";
	MyQuads = "X["<>StringDrop[StringDrop[ToString[#],1],-1]<>"],"& /@ L;
	ToExpression[StringDrop[StringJoin["PD[",MyQuads],-1]<>"]"]
]

AllArcs[ThisPD_] := DeleteDuplicates[Flatten[PDtoList[ThisPD],1]];
AllComponents[ThisPD_] := Module[{compRule},
	compRule = (compList_ :> DeleteDuplicates[Append[compList,NextArc[compList[[-1]],ThisPD]]]);
	Sort[DeleteDuplicates[Sort /@ ((# //. compRule)& /@ Partition[AllArcs[ThisPD],1])],Min[#1] < Min[#2]&]
];
AllArcs[component_,ThisPD_] := AllComponents[ThisPD][[component]];
NumComponents[ThisPD_] := Length[AllComponents[ThisPD]];

(* Here's a cute example of overloading: ComponentNumber can either be called *)
(* with a PD or with the list of components already worked out (for speed) *)
ComponentNumber[arc_,ThisPD_PD] := ComponentNumber[arc,AllComponents[ThisPD]];
ComponentNumber[ThisPD_PD,arc_Integer] := ComponentNumber[arc,ThisPD];
ComponentNumber[arc_,components_List] := If[MemberQ[components,arc,2],Position[components,arc][[1,1]],
	Print["CompNum:arc ",arc," is not in list of components ",components]];


ArcCounts[ThisPD_] := Module[{pdl}, pdl = Flatten[PDtoList[ThisPD],1]; {#,Count[pdl,#]}& /@ DeleteDuplicates[pdl]];

OrientationAtCrossing[arc_,crossing_,ThisPD_PD] := Module[{pos,otherArc,otherAdjacent},

	pos = Position[ThisPD[[crossing]],arc];

	If[!MemberQ[{{{1}},{{2}},{{3}},{{4}}},pos],
		Print["Arc ",arc," not part of crossing ",crossing," (",ThisPD[[crossing]],")"];Return[];
	];

	Which[
		pos == {{1}},
		Return["In"],
		pos == {{3}},
		Return["Out"]
	];

	otherArc = If[pos == {{2}},ThisPD[[crossing,4]],ThisPD[[crossing,2]]];
	(*Print["otherArc = ",otherArc];*)
	otherAdjacent = TrueQ[Abs[arc-otherArc]==1];
	(*Print[If[otherAdjacent,"adjacent","not adjacent"]];*)
	If[Xor[otherAdjacent,Order[arc,otherArc]==1],"Out","In"]

];

ArcOverAtCrossingQ[arc_,crossing_,ThisPD_PD] := If[MemberQ[ThisPD[[crossing]],arc], 
	Or[Position[ThisPD[[crossing]],arc] == {{2}}, Position[ThisPD[[crossing]],arc] == {{4}}],
	Print["ArcOverAtCrossingQ:arc",arc," not part of crossing ",crossing," (",ThisPD[[crossing]],")"];];

ArcUnderAtCrossingQ[arc_,crossing_,ThisPD_PD] := !ArcOverAtCrossingQ[arc,crossing,ThisPD];

ArcOverAtHeadQ[arc_,ThisPD_PD] := ArcOverAtCrossingQ[arc,ArcHead[arc,ThisPD],ThisPD];
ArcUnderAtHeadQ[arc_,ThisPD_PD] := !ArcOverAtHeadQ[arc,ThisPD];

ArcOverAtTailQ[arc_,ThisPD_PD] := ArcOverAtCrossingQ[arc,ArcTail[arc,ThisPD],ThisPD];
ArcUnderAtTailQ[arc_,ThisPD_PD] := !ArcOverAtTailQ[arc,ThisPD];

ArcAlternatingQ[arc_,ThisPD_PD] := Xor[ArcOverAtHeadQ[arc,ThisPD],ArcOverAtTailQ[arc,ThisPD]];
(* This is True when exactly one of the head and tail go over and under (or under and over) *) 

OrientationAtCrossing[arc_,crossing_,ThisEPD_EPD] := 
	Which[arc[[3,1]] == crossing,"Out",arc[[3,2]] == crossing,"In",
	True,Print["OrientationAtCrossing: arc ",arc," doesn't hit crossing ",crossing]];

CrossingsForArc[arc_,ThisPD_PD] := Module[{crossings,pdlist},
	pdlist = PDtoList[ThisPD];
	crossings = Position[pdlist,arc];
	crossings = Transpose[crossings][[1]];
	If[OrientationAtCrossing[arc,crossings[[1]],ThisPD] == "In",Reverse[crossings],crossings]
];

CrossingsForArc[arc_,ThisEPD_EPD] := arc[[3]];

ArcHead[arc_,ThisPD_] := CrossingsForArc[arc,ThisPD][[2]];
ArcTail[arc_,ThisPD_] := CrossingsForArc[arc,ThisPD][[1]];

OrientationsForArc[arc_,ThisPD_] := OrientationAtCrossing[arc,#,ThisPD] & /@ CrossingsForArc[arc,ThisPD];

ArcToLeft[arc_,crossing_,ThisPD_] := Module[{pos,or,leftpos},
	pos = Position[ThisPD[[crossing]],arc];
	(*Print["pos:",pos];*)
	If[Length[pos]==0,
	Print["Arc ",arc," not part of crossing ",crossing," (",ThisPD[[crossing]],")"];Return[];];

	or = OrientationAtCrossing[arc,crossing,ThisPD];
	(*Print["or:",or];*)
	leftpos = 
		Which[
			pos == {{1}},4,
			pos == {{2}},If[or == "In",1,3],
			pos == {{3}},4,
			pos == {{4}}, If[or == "In",3,1]
		];
	ThisPD[[crossing,leftpos]]
];

NextArc[arc_,ThisPD_] := Module[{headCrossing,pos},

	headCrossing = ArcHead[arc,ThisPD];
	pos = Position[ThisPD[[headCrossing]],arc][[1,1]];
	RotateLeft[ThisPD[[headCrossing]],2][[pos]]

];


(* ::Section:: *)
(*BlackboardDouble*)


(* ::Text:: *)
(*This part of the package allows us to double a given component of the link using the blackboard framing. *)
(*This process has several steps:*)
(**)
(*1) Find the arcs which belong to the given component, and pick the first. Compute the "self" and "other" crossings of the given component in the original PD code. We will add 3*self + other crossings to the diagram before we're through. *)
(**)
(*2) Convert the PD code into an "Expanded" PD code which maintains orientation information for each arc. *)
(**)
(*3) Follow the arcs of the given component, adding new arcs as we go and splitting each left arc as we find it, fixing the crossings we encounter at the same time. *)
(**)
(*4) Each time we process an arc, we add a new arc, split an old arc, and add a new crossing. We put a 'T' in the PD called a "tie" which is a symbol specifying that we will later tie an arc into that position in the crossing. (We use a replacement rule to perform the tie). *)
(**)
(*5) After adding the correct number of crossings, we go ahead and call TieLooseArc to find the first arc we created (and did NOT tie in, since no "T" existed at that point)*)
(**)
(*6) Finally, we call EPDtoPD to relabel the arcs and convert back to a standard PD format.*)
(**)
(*The functions in this block are all intended to by internal functions for the BlackboardDouble function. The new component of the resulting PD is guaranteed to be the LAST component. *)


CrossingCounts[ThisPD_PD] := Module[{comps,compList,compRules,i,modPd,CC},
	compList = AllComponents[ThisPD];	
	comps = Length[compList];
	compRules = Table[i -> Position[compList,i][[1,1]],{i,1,Length[Flatten[compList]]}];
	modPd = DeleteDuplicates /@ (PDtoList[ThisPD] /. compRules);
	CC = Outer[List,Range[comps],Range[comps]];
	CC /. {i_Integer,j_Integer} :> If[i!=j,Count[modPd,{i,j}|{j,i}],Count[modPd,{i}]]
];

NewCrossingsToBBD[comp_,ThisPD_PD] := Module[{CC},
	CC = CrossingCounts[ThisPD];
	Total[CC[[comp,;;]]] + 2*CC[[comp,comp]]
]; 

NewPrimaryIndex[ThisEPD_EPD] := Module[{aa,last},
	aa = AllArcs[ThisEPD];
	last = Sort[aa][[-1]];
	last[[1]]+1
];

TieCrossing[ThisEPD_EPD] := If[MemberQ[ThisEPD,"T",2],Position[ThisEPD,"T"][[1,1]],"None"];
NewCrossing[ThisEPD_EPD] := Length[ThisEPD] + 1;
ArcsWithPrimary[primary_,ThisEPD_] := Cases[AllArcs[ThisEPD],{arc[[1]],_,_}];

SplitArc[arc_,ThisEPD_EPD] := Module[{headCr,newEPD,splitOut,splitIn},
    (* Shift the subscripts of all arcs with the same primary as us and larger subscripts up by one to make room *)
	newEPD = ThisEPD /. ({arc[[1]],k_,cr_} /; k > arc[[2]]) :> {arc[[1]],k+1,cr};
    (* the new outgoing head-half and tail-half arcs are missing a crossing at the split and out's subscript is altered *)
	splitOut = arc; splitOut[[2]] += 1; splitOut[[3,1]] = "S"; (* This arc has no tail crossing (yet) *)
	splitIn = arc; splitIn[[3,2]] = "S"; (* This arc has no head crossing (yet) *)
	(* We now update the original EPD with these changes *) 
	newEPD[[ArcHead[arc,ThisEPD]]] = newEPD[[ArcHead[arc,ThisEPD]]] /. arc -> splitOut;
	If[ArcTail[arc,ThisEPD] != "None",
		newEPD[[ArcTail[arc,ThisEPD]]] = newEPD[[ArcTail[arc,ThisEPD]]] /. arc -> splitIn;];
	(* And return the new arcs and EPD *)
	{splitIn,splitOut,newEPD}
];

AddLeftArcAndCrossing[arc_,ThisEPD_EPD] := 
	Module[{newEPD,newArc,leftArc,newCross,tieCross,splitOut,splitIn,newX,pos},
		(* First, we generate the new arc, and tie its tail into the Tie symbol if it exists *)
		tieCross = TieCrossing[ThisEPD];
		newCross = NewCrossing[ThisEPD];
		newArc = {NewPrimaryIndex[ThisEPD],1,{tieCross,newCross}};
		newEPD = ThisEPD /. "T" -> newArc;
		(* Now we figure out what arc the head of the new arc will impact at its head crossing and split it *)
		leftArc = ArcToLeft[arc,ArcHead[arc,newEPD],newEPD];
		{splitIn,splitOut,newEPD} = SplitArc[leftArc,newEPD];
        (* We now construct the new crossing *)
		pos = Position[ThisEPD[[ArcHead[arc,ThisEPD]]],arc];
		newX = Which[
			   pos == {{1}},
			   If[OrientationAtCrossing[leftArc,ArcHead[arc,ThisEPD],ThisEPD] == "In",X[newArc,splitOut,"T",splitIn],X[newArc,splitIn,"T",splitOut]],
			   pos == {{2}},
			   X[splitIn,newArc,splitOut,"T"],
			   pos == {{4}},
			   X[splitIn,"T",splitOut,newArc]];
		newEPD = AppendTo[newEPD,newX];
		(* Finally, all the "split" ends come together at the new crossing *)
		newEPD /. "S" -> newCross
];

TieLooseArc[ThisEPD_EPD] := Module[{loosePos,looseArc,pdlist,arcs},
	(* The loose arc has a tail crossing of "None" *)
	If[MemberQ[ThisEPD,"None"],
		Print["TieLooseArc: No loose arc to tie in."];
		ThisEPD,
		loosePos = Position[ThisEPD,"None"];
	    (*Print["loosePos:",loosePos];*)
		looseArc = ThisEPD[[loosePos[[1,1]],loosePos[[1,2]]]];	
		(*Print["TieLooseArc: Tying loose arc ",looseArc," into crossing ",TieCrossing[ThisEPD]];*)		
		looseArc[[3,1]] = TieCrossing[ThisEPD];
		ThisEPD /. {{_,_,{"None",_}} -> looseArc, "T" -> looseArc}]
];	

BlackboardDouble[component_,ThisPD_PD] := Module[{newEPD,arc,startArc,startPrimary,i,nc},

	Print["BBDouble:Building blackboard double of ",Length[AllArcs[component,ThisPD]]," arc component ",component," of ",Length[ThisPD]," crossing PD."];
	startPrimary = AllArcs[component,ThisPD][[1]];

	newEPD = PDtoEPD[ThisPD];
	(*Print["Orig EPD:",newEPD];*)
	startArc = {startPrimary,1,CrossingsForArc[startPrimary,ThisPD]};
	(*Print["startArc:",startArc];*)
	newEPD = AddLeftArcAndCrossing[startArc,newEPD];
	(*Print["After adding along startArc:",newEPD];*)
	nc = NewCrossingsToBBD[component,ThisPD];
	Print["BBDouble:Will add ",nc," crossings to double this component."];

	For[arc = NextArc[startArc,newEPD];i=1;,i < nc,arc = NextArc[arc,newEPD];i++,
		(*Print["arc now:",arc];*)
	    newEPD = AddLeftArcAndCrossing[arc,newEPD];
	];

	EPDtoPD[TieLooseArc[newEPD]]

];	


(* ::Section:: *)
(*Faces of PD diagrams*)


(* ::Text:: *)
(*This function makes a list of faces. This is often useful for various relatively simpleminded knot invariants, such as making counts of the number of edges in various faces of a knot diagram.*)


DebugPDFaces = False;

CrossingSlot[arc_,crossing_,ThisPD_] := 
	If[MemberQ[ThisPD[[crossing]],arc],Position[ThisPD[[crossing]],arc][[1,1]],
	Print["CrossingSlot: Arc ",arc," is not part of crossing ",\
	   crossing," (",ThisPD[[crossing]],") of PD."];];

OrientedArcHead[OArc_,ThisPD_] := If[OArc[[2]] == "f",ArcHead[OArc[[1]],ThisPD],ArcTail[OArc[[1]],ThisPD]];

OrientedArcToLeft[OArc_,ThisPD_] := Module[{leftArc,xRec,cr,slot},
	(*Print["OArc:",OArc]; Print["Oarc head:",OrientedArcHead[OArc,ThisPD]];*)
	cr = OrientedArcHead[OArc,ThisPD];
	xRec = ThisPD[[cr]];
	slot = CrossingSlot[OArc[[1]],cr,ThisPD];
	(*Print["OA2L: arc:",OArc," cr:",cr," xRec:",xRec," slot:",slot];*)
	(* Now that we've assembled information, we can compute *)
	leftArc = RotateRight[xRec][[slot]];
	(*Print["leftArc:",leftArc];*)
	{leftArc,If[OrientationAtCrossing[leftArc,cr,ThisPD]=="Out","f","r"]}
];

DeleteOrientedArc[OArc_,Arcs_] := 
	If[MemberQ[Arcs,OArc],(*Print["Deleting ",OArc," from ",Arcs];*) Delete[Arcs,Position[Arcs,OArc]],Arcs];

OrientedPDFaces[ThisPD_PD] := Module[{arcs,faces,newFace,thisArc,i},
	(* Make a list of pairs in the form {n,"f"} and {n,"r"} for arcs. *)
	arcs = Partition[Append[Riffle[AllArcs[ThisPD],"f"],"f"],2];
	arcs = Sort[Join[arcs,arcs /. "f" -> "r"]];
	(*Print[arcs];*)
	faces = {}; i = 0;
	(* Now we start to walk *)
	While[Length[arcs] > 0,
			newFace = {arcs[[1]]}; 
			thisArc = OrientedArcToLeft[newFace[[1]],ThisPD]; 
			arcs=DeleteOrientedArc[arcs[[1]],arcs];
			(*Print["First arc:",newFace[[1]]]; i++; *)
			(*Print["Moved to arc:",thisArc," (should stop at ",newFace[[1]],")"];*)
			While[thisArc != newFace[[1]],
				AppendTo[newFace,thisArc]; 
				arcs = DeleteOrientedArc[thisArc,arcs];
				thisArc = OrientedArcToLeft[thisArc,ThisPD];
				(*Print["Moved to arc:",thisArc," (should stop at ",newFace[[1]],")"];*)
			    i++;
			];
			(*Print["Created new face: ",newFace, " ",Length[arcs]," arcs left"];*)
			AppendTo[faces,newFace];
			(*Print["arcs now: ",arcs];*)
	];
	faces
];

PDFaces[ThisPD_PD] := OrientedPDFaces[ThisPD] /. {i_Integer,"f"|"r"} :> i



(* ::Section:: *)
(*NiceDrawPD*)


(* ::Text:: *)
(*DrawPD is a really clever piece of code, but it has some weaknesses. Among them, you can't display orientations, have much control over how different components come out, or label arcs or crossings. NiceDrawPD is a wrapper which solves these problems by post-processing the output of DrawPD. We need to make a small modification to the DrawPD code in order to correctly parse the results.*)


 (* Begin source file src/DrawPD.m*)

(* BeginPackage["KnotTheory`"] *)

(* PD; X; OuterFace; Gap; Colour; StrandColour *)

WhittenDrawPD::usage = "
  WhittenDrawPD[pd] takes the planar diagram description pd and creates a
  graphics object containing a picture of the knot.
  WhittenDrawPD[pd,options], where options is a list of rules, allows the user
  to control some of the parameters.  OuterFace->n sets the face at
  infinity to the face numbered n.  OuterFace->{e_1,e_2,...,e_n} sets
  the face at infinity to a face which has edges e_1, e_2, ..., e_n in
  the planar diagram description.  Gap->g sets the size of the gap
  around a crossing to length g.
"

WhittenDrawPD::about = "
  WhittenDrawPD was written by Emily Redelmeier at the University of Toronto in
  the summers of 2003 and 2004. Modified in 2011 by Jason Cantarella at University of Georgia.
"

Begin["WhittenDrawPD`"] (* was `WhittenDrawPD` *)

(* Representation Manipulation *)

(* Positions of the various fields *)
neighbours=1;
type=2;
r=3;
centre=4;
graphicsObjs=5;

FieldValues[triangulation_,field_]:=
  Table[triangulation[[i,field]],{i,
      Length[triangulation]}]

AddField[triangulation_,field_,values_]:=
  Table[ReplacePart[
      If[Length[triangulation[[i]]]<
          field,PadRight[
          triangulation[[i]],field],
        triangulation[[i]]],
      values[[i]],field],{i,
      Length[triangulation]}]

ChangeField[triangulation_,field_,f_]:=
  MapAt[f,triangulation,Table[{i,field},{i,Length[triangulation]}]]

DeriveField[triangulation_,field_,f_]:=
  AddField[triangulation,field,Map[f,triangulation]]

(* PD Graph Manipulation *)

OtherVertex[pd_,coordinates_]:=
  Complement[
      Position[pd,
        Extract[pd,
          coordinates]],{coordinates}][[1\
]]

Faces[pd_]:= Module[{i,j},
  Select[Flatten[
      Table[NestWhileList[
          Function[{coordinates},{OtherVertex[pd,
                  coordinates][[1]],
              Mod[OtherVertex[pd,
                      coordinates][[2]]-1,Length[pd[[OtherVertex[pd,\
                        coordinates][[1]]]]],1]}],{i,j},Unequal,All,Infinity,-1],\
		  {i,Length[pd]},{j,Length[pd[[i]]]}],1],
    Function[{face},
      face[[1]]==
        Sort[face][[1]]]]];

Triangulate[pd_]:=Module[{i},
	If[DebugDrawPD,Print["madeTcall"]];
	(facelist=Faces[pd];
    Join[Table[{Flatten[
            Table[{Length[pd]+
                  pd[[vertex,i]],
                Length[pd]+Length[Union[Flatten[pd,1,X]]]+
                  Position[facelist,face_/;MemberQ[face,{vertex,i}],
                      1][[1,1]]},{i,
                Length[pd[[vertex]]]}]],
          "X"},{vertex,Length[pd]}],
      Table[{Flatten[
            Table[{Length[pd]+Length[Union[Flatten[pd,1,X]]]+
                  Position[facelist,
                      face_/;MemberQ[face,
                          Position[pd,edge,2][[
                            i]]],2][[1,
                    1]],
                Position[pd,edge,2][[i,
                  1]]},{i,Length[Position[pd,edge,2]]}]],
          "e"},{edge,Length[Union[Flatten[pd,1,X]]]}],
      Table[{Flatten[
            Table[{facelist[[face,i,1]],
                Length[pd]+
                  Extract[pd,
                    facelist[[face,
                      i]]]},{i,
                Length[facelist[[
                    face]]]}]],"f"},{face,
          Length[facelist]}]])]

GetOuterFace[triangulation_,edges_]:=
  Select[Range[Length[triangulation]],
      Function[v,
        triangulation[[v,type]]==
            "f"&&Union[
              Map[triangulation[[#]]&,
                triangulation[[v,
                  neighbours]]],
              Select[triangulation,#[[type\
]]=="e"&][[edges\
]]]==Union[
              Map[triangulation[[#]]&,
                triangulation[[v,
                  neighbours]]]]]][[1\
]]

NthOrderNeighbours[triangulation_,v_,0]:={v}
NthOrderNeighbours[triangulation_,v_,n_/;n>0]:=
  Apply[Union,
    Map[triangulation[[#,neighbours]]&,
      NthOrderNeighbours[triangulation,v,n-1]]]

DefaultOuterFace[triangulation_]:=
  Sort[Select[Range[Length[triangulation]],
        triangulation[[#,
              type]]=="f"&],
      Length[triangulation[[#1,
              neighbours]]]>=Length[
            triangulation[[#2,
              neighbours]]]&][[1\
]]

(* generalize to accept other types of vertices (v, etc.) *)

(* Circle Packing: Radii *)

CircleAngle[r_,r1_,r2_]:=
  ArcCos[((r1+r)^2+(r2+r)^2-(r1+r2)^2)/(2(r1+r)(r2+r))]

FlowerAngle[triangulation_,v_,radii_]:=
  Plus@@(CircleAngle[radii[[v]],#1,#2]&@@@
        Map[radii[[#]]&,
          Transpose[{triangulation[[v,
                neighbours]],
              RotateRight[
                triangulation[[v,
                  neighbours]]]}],{2}])

AdjustRadius[triangulation_,v_,targetAngle_,radii_]:=
  N[radii[[v]]((1-
              Cos[FlowerAngle[triangulation,v,radii]/
                  Length[triangulation[[v,
                      neighbours]]]]+
              Sqrt[2-2Cos[
                      FlowerAngle[triangulation,v,radii]/
                        Length[
                          triangulation[[v,
                            neighbours]]]]])/(1+
              Cos[FlowerAngle[triangulation,v,radii]/
                  Length[triangulation[[v,
                      neighbours]]]]))(Sqrt[
            2/(1-Cos[
                    targetAngle/
                      Length[triangulation[[v,
                          neighbours]]]])]-1)]

PackingStep[triangulation_,targetAngles_]:=(
    radii=Table[Unique[radius],{Length[triangulation]}];
    Compile[Evaluate[radii],
      Evaluate[Table[
          If[targetAngles[[v]]==0,
            radii[[v]],
            AdjustRadius[triangulation,v,
              targetAngles[[v]],
              radii]],{v,Length[triangulation]}]]]
    )

GetRadii[triangulation_,targetAngles_,radii_]:=(
    EvaluatedPackingStep=PackingStep[triangulation,targetAngles];
    NestWhile[EvaluatedPackingStep@@#&,radii,Unequal,2]
    )

DefaultDirichlet[triangulation_]:= Module[{},
	If[DebugDrawPD,Print["Made DD call on tri:",triangulation]];
  AddField[triangulation,r,
    GetRadii[triangulation,
      ReplacePart[Table[2Pi,{Length[triangulation]}],
        0,{{1},{triangulation[[1,neighbours,
              1]]},{triangulation[[1,
              neighbours,2]]}}],
      Table[1,{Length[triangulation]}]]]];

(* Circle Packing: Positions *)

PlaceFlower[triangulation_,v_,neighbour_]:=Module[{i},
  For[w=triangulation[[v,neighbours,
        neighbour]];
    theta=Arg[(z[w]-
              z[v])/(triangulation[[w,
                r]]+
              triangulation[[v,r]])];
    lastr=triangulation[[w,r]];i=1,
    i<Length[triangulation[[v,
          neighbours]]],i++,
    w=triangulation[[v,neighbours,
        Mod[neighbour+i,
          Length[triangulation[[v,
              neighbours]]],1]]];
    currentr=triangulation[[w,r]];
    theta+=CircleAngle[
        triangulation[[v,r]],lastr,
        currentr];lastr=currentr;
    z[w]=z[v]+(triangulation[[v,r]]+
              currentr)Exp[I*theta];placed=Union[placed,{w}]]];

PackCircles[triangulation_]:=(Clear[z];placed={};surrounded={};v1=1;z[v1]=0;
    placed=Union[placed,{v1}];
    v2=triangulation[[v1,neighbours,1]];
    z[v2]=triangulation[[v1,r]]+
        triangulation[[v2,r]];
    placed=Union[placed,{v2}];
    While[Length[placed]!=Length[triangulation],
      v=Complement[placed,
            surrounded][[1]];
      PlaceFlower[triangulation,v,
        Position[
            triangulation[[v,
              neighbours]],
            w_/;MemberQ[placed,w]][[1,
          1]]];surrounded=Union[surrounded,{v}]];
    Table[z[i],{i,Length[triangulation]}])

PackCirclesBound[triangulation_]:=(Clear[z];placed={};surrounded={};
    v1=Select[Range[Length[triangulation]],
          FlowerAngle[triangulation,#,
                FieldValues[triangulation,
                  r]]==2Pi&][[1\
]];z[v1]=0;placed=Union[placed,{v1}];
    v2=triangulation[[v1,neighbours,1]];
    z[v2]=triangulation[[v1,r]]+
        triangulation[[v2,r]];
    placed=Union[placed,{v2}];
    While[Length[placed]!=Length[triangulation],
      v=Select[Complement[placed,surrounded],
            FlowerAngle[triangulation,#,
                  FieldValues[triangulation,
                    r]]==2Pi&][[1\
]];
      PlaceFlower[triangulation,v,
        Position[
            triangulation[[v,
              neighbours]],
            w_/;MemberQ[placed,w]][[1,
          1]]];surrounded=Union[surrounded,{v}]];
    Table[z[i],{i,Length[triangulation]}])

AddPositions[triangulation_]:=
  AddField[triangulation,centre,PackCircles[triangulation]]

AddPositionsBound[triangulation_]:=
  AddField[triangulation,centre,PackCirclesBound[triangulation]]

(* Fractional Linear Transformations *)

NewRadius[z_,radius_,{{a_,b_},{c_,d_}}]:=
  radius*Abs[a*d-b*c]/(Abs[c*z+d]^2-Abs[c]^2*radius^2)

NewPosition[z_,
    radius_,{{a_,b_},{c_,d_}}]:=((a*z+b)*Conjugate[c*z+d]-
        a*Conjugate[c]*radius^2)/((c*z+d)*Conjugate[c*z+d]-
        c*Conjugate[c]*radius^2)

ApplyFLMap[
    triangulation_,{{a_,b_},{c_,d_}}]:=(newRadii=
      Table[NewRadius[
          triangulation[[v,centre]],
          triangulation[[v,
            r]],{{a,b},{c,d}}],{v,Length[triangulation]}];
    newPositions=
      Table[NewPosition[
          triangulation[[v,centre]],
          triangulation[[v,
            r]],{{a,b},{c,d}}],{v,Length[triangulation]}];
    AddField[AddField[triangulation,r,newRadii],centre,newPositions])

Moebius[a_]:={{1,-a},{-Conjugate[a],1}}

ComposeMoebius[a_,b_]:=(a+b)/(1+a*Conjugate[b])

(* Inversion *)

PutInside[triangulation_,outerFace_]:=
  ApplyFLMap[
    triangulation,{{0,
        triangulation[[outerFace,
          r]]},{1,-triangulation[[
            outerFace,centre]]}}]

(* Balancing *)

BalanceStep[triangulation_,moebiusConst_]:=
  ComposeMoebius[
    Plus@@FieldValues[ApplyFLMap[triangulation,Moebius[moebiusConst]],centre]/
      Length[triangulation],moebiusConst]

BalanceMoebius[triangulation_]:=FixedPoint[BalanceStep[triangulation,#]&,0]

Balance[triangulation_]:=
  ApplyFLMap[triangulation,Moebius[BalanceMoebius[triangulation]]]

(* Graphics *)

xyCoords[z_]:={Re[z],Im[z]}

(* Graphics: Circle Packing *)

PackingGraphics[triangulation_]:=Module[{tr},
	If[DebugDrawPD,Print["Triangulation:",triangulation]];
	Graphics[Join[
      Map[
        Circle[xyCoords[#[[centre]]],
            Abs[#[[r]]]]&,triangulation],
      Table[Text[i,
          xyCoords[
            triangulation[[i,
              centre]]]],{i,Length[triangulation]}]],
    AspectRatio->1]
];
(* Knot Manipulation *)

ConnectedNeighbours[triangulation_,v_]:={{1,5},{3,7}}/;
    triangulation[[v,type]]=="X"
ConnectedNeighbours[triangulation_,v_]:={{2,4}}/;
    triangulation[[v,type]]=="e"
ConnectedNeighbours[triangulation_,v_]:={}/;
    triangulation[[v,type]]=="f"

Xunder=1;
Xover=2;

OtherEnd[triangulation_,v_,neighbour_]:=
  Select[Select[
          Map[triangulation[[v,
                neighbours,#]]&,
            ConnectedNeighbours[triangulation,v]],
          MemberQ[#,
              neighbour]&][[1]],# !=neighbour&][[1]]

AdjacentComponents[triangulation_,{v_,n_}]:=
  Select[Map[Take[#,2]&,
      Position[Table[
          Map[triangulation[[w,
                neighbours,#]]&,
            ConnectedNeighbours[triangulation,w]],{w,Length[triangulation]}],
        v]],Function[component,
      MemberQ[Map[
          triangulation[[v,
              neighbours,#]]&,
          ConnectedNeighbours[triangulation,v][[
            n]]],
        component[[1]]]]]

GetStrand[triangulation_,{v_,n_}]:=
  FixedPoint[
    Apply[Union,
        Append[Map[
            Function[component,
              AdjacentComponents[triangulation,component]],#],#]]&,{{v,n}}]

ListStrands[triangulation_]:=
  Union[Map[GetStrand[triangulation,#]&,
      Flatten[Table[{v,n},{v,Length[triangulation]},{n,
            Length[ConnectedNeighbours[triangulation,v]]}],1]]]

(* Graphics: Graphs *)

gapParam=1;

ArcCentre[z_,radius_,{z1_,z2_}]:=2*radius/Conjugate[Sign[z1-z]+Sign[z2-z]]

ArcRadius[z_,radius_,{z1_,z2_}]:=
  Sqrt[Abs[ArcCentre[z,radius,{z1,z2}]]^2-radius^2]

LinearProject[{v1_,v2_},w_]:=Re[v2-v1]*Cos[Arg[w]]+Im[v2-v1]*Sin[Arg[w]]

ArcProject[{v1_,v2_},o_]:=Mod[Arg[v2-o]-Arg[v1-o],2*Pi,-\[Pi]]

ArcOrientation[z_,radius_,{z1_,z2_}]:=
  Sign[ArcProject[{radius*Sign[z1-z],radius*Sign[z2-z]},
      ArcCentre[z,radius,{z1,z2}]]]

(*find a way to consolidate the crossing functions*)

ArcCrossing[
    radius_,{arc1_,
      arc2_}]:=(Re[arc1]*Im[arc2]-Im[arc1]*Re[arc2]-
          Sign[Re[arc1]*Im[arc2]-Im[arc1]*Re[arc2]]*
            Sqrt[(Re[arc1]*Im[arc2]-Im[arc1]*Re[arc2])^2-
                Abs[arc1-arc2]^2*radius^2])/
      Abs[arc1-arc2]^2*(arc1-arc2)*I

ArcLineCrossing[
    radius_,{arc_,
      line_}]:=((Re[arc]*Re[line]+Im[arc]*Im[line]-
            Sign[Re[arc]*Re[line]+Im[arc]*Im[line]]*
              Sqrt[(Re[arc]*Re[line]+Im[arc]*Im[line])^2-
                  Abs[line]^2*radius^2])/Abs[line]^2)*line

Crossing[z_,radius_,{{z1_,z2_},{z3_,z4_}}]:=
  If[Sign[z1-z]+Sign[z2-z]==0,
    If[Sign[z3-z]+Sign[z4-z]==0,0,
      ArcLineCrossing[radius,{ArcCentre[z,radius,{z3,z4}],z2-z1}]],
    If[Sign[z3-z]+Sign[z4-z]==0,
      ArcLineCrossing[radius,{ArcCentre[z,radius,{z1,z2}],z4-z3}],
      ArcCrossing[
        radius,{ArcCentre[z,radius,{z1,z2}],ArcCentre[z,radius,{z3,z4}]}]]]

arcConst=10^(-5);

ArcDistance[z_,radius_,{z1_,z2_},{w1_,w2_}]:=
  Which[Sign[z1-z]+Sign[z2-z]==0,LinearProject[{w1,w2},z2-z1],
    Abs[ArcProject[{w1,w2},ArcCentre[z,radius,{z1,z2}]]](*==0*)(**)<
      arcConst(**),LinearProject[{w1,w2},
      ArcOrientation[z,
          radius,{z1,z2}]*I*((w1+w2)/2-
            ArcCentre[z,radius,{z1,z2}])],True,
    ArcOrientation[z,radius,{z1,z2}]*ArcRadius[z,radius,{z1,z2}]*
      ArcProject[{w1,w2},ArcCentre[z,radius,{z1,z2}]]]

(* The circles in the drawing are produced with this call. *)
(* THIS IS THE ONLY CHANGE IN SOURCE (EXCEPT PRINTS) *)

GetArc[z_,radius_,{z1_,z2_},{w1_,w2_},{l1_,l2_}]:=
  Which[ArcDistance[z,radius,{z1,z2},{w1,w2}]<
      l1+l2,{"Arc Deleted"},(*Mod[
          Arg[w1-ArcCentre[z,radius,{z1,z2}]]+
            ArcOrientation[z,radius,{z1,z2}]*l1/ArcRadius[z,radius,{z1,z2}],
          2*Pi,Arg[w1-ArcCentre[z,radius,{z1,z2}]]]==
        Mod[Arg[z2-ArcCentre[z,radius,{z1,z2}]]-
            ArcOrientation[z,radius,{z1,z2}]*l2/ArcRadius[z,radius,{z1,z2}],
          2*Pi,Arg[w1-ArcCentre[z,radius,{z1,z2}]]]*)(**)
      Abs[ArcProject[{w1,w2},ArcCentre[z,radius,{z1,z2}]]-
          l1/ArcRadius[z,radius,{z1,z2}]-l2/ArcRadius[z,radius,{z1,z2}]]<
      arcConst(**),{Line[{xyCoords[
            ArcCentre[z,radius,{z1,z2}]+
              ArcRadius[z,radius,{z1,z2}]*
                Exp[I*(Arg[w1-ArcCentre[z,radius,{z1,z2}]]+
                        ArcOrientation[z,radius,{z1,z2}]*
                          l1/ArcRadius[z,radius,{z1,z2}])]],
          xyCoords[
            ArcCentre[z,radius,{z1,z2}]+
              ArcRadius[z,radius,{z1,z2}]*
                Exp[I*(Arg[w2-ArcCentre[z,radius,{z1,z2}]]-
                        ArcOrientation[z,radius,{z1,z2}]*
                          l2/ArcRadius[z,radius,{z1,z2}])]]}]},
    True,{Circle[xyCoords[z+ArcCentre[z,radius,{z1,z2}]],
        ArcRadius[z,radius,{z1,z2}],
        Sort[{Mod[
              Arg[w1-ArcCentre[z,radius,{z1,z2}]]+
                ArcOrientation[z,radius,{z1,z2}]*
                  l1/ArcRadius[z,radius,{z1,z2}],2*Pi,
              Which[ArcOrientation[z,radius,{z1,z2}]>0,
                Arg[radius*Sign[z1]-ArcCentre[z,radius,{z1,z2}]]-Pi/2,
                ArcOrientation[z,radius,{z1,z2}]<0,
                Arg[radius*Sign[z2]-ArcCentre[z,radius,{z1,z2}]]-Pi/2]],
            Mod[Arg[w2-ArcCentre[z,radius,{z1,z2}]]-
                ArcOrientation[z,radius,{z1,z2}]*
                  l2/ArcRadius[z,radius,{z1,z2}],2*Pi,
              Which[ArcOrientation[z,radius,{z1,z2}]>0,
                Arg[radius*Sign[z1]-ArcCentre[z,radius,{z1,z2}]]-Pi/2,
                ArcOrientation[z,radius,{z1,z2}]<0,
                Arg[radius*Sign[z2]-ArcCentre[z,radius,{z1,z2}]]-Pi/2]]},
          Less]]}]

CircleParams[triangulation_,v_,
    n_]:={triangulation[[v,centre]],
    triangulation[[v,r]],
    Map[triangulation[[triangulation[[v,
            neighbours,#]],centre]]&,
      Map[ConnectedNeighbours[triangulation,
              v][[#]]&,n,{-1}],{-1}]}

ExtraGap[triangulation_,v_,neighbour_,gap_]:=
  If[MemberQ[
        Map[triangulation[[v,
              neighbours,#]]&,
          ConnectedNeighbours[triangulation,v][[
            Xunder]]],neighbour],
      Max[0,gap-
          Abs[Apply[ArcDistance,
              Join[CircleParams[triangulation,v,
                  Xunder],{{Apply[Crossing,
                      CircleParams[triangulation,v,{Xunder,Xover}]],
                    triangulation[[v,r]]*
                      Sign[triangulation[[neighbour,
                            centre]]-
                          triangulation[[v,
                            centre]]]}}]]]],0]/;
    triangulation[[v,type]]=="X"
ExtraGap[triangulation_,v_,neighbour_,gap_]:=
  If[MemberQ[
        Map[triangulation[[v,
              neighbours,#]]&,
          ConnectedNeighbours[triangulation,
              v][[1]]],neighbour],
      Max[0,ExtraGap[triangulation,OtherEnd[triangulation,v,neighbour],v,gap]-
          Apply[ArcDistance,
            Join[CircleParams[triangulation,v,
                1],{Map[
                  triangulation[[v,r]]*
                      Sign[triangulation[[#,
                            centre]]-
                          triangulation[[v,
                            centre]]]&,
                  Map[triangulation[[v,
                        neighbours,#]]&,
                    ConnectedNeighbours[triangulation,
                        v][[1]]]]}]]],0]/;
    triangulation[[v,type]]=="e"

DefaultGap[triangulation_]:=
  
  Min[Table[
      If[triangulation[[v,
            type]]=="X",
        Map[Abs[Apply[ArcDistance,
                Join[CircleParams[triangulation,v,
                    Xunder],{{triangulation[[v,
                          r]]*
                        Sign[triangulation[[#,
                              centre]]-
                            triangulation[[v,
                              r]]],
                      Apply[Crossing,
                        CircleParams[triangulation,v,{Xunder,Xover}]]}}]]]&,
          Map[triangulation[[v,
                neighbours,#]]&,
            ConnectedNeighbours[triangulation,v][[
              Xunder]]]],Infinity],{v,
        Length[triangulation]}]]

(* Draw a crossing *)

GetGraphicsObjs[triangulation_,v_,
    graphicsParams_]:={Join[
        Apply[GetArc,
          Join[CircleParams[triangulation,v,
              Xunder],{{triangulation[[v,
                    r]]*
                  Sign[triangulation[[
                        triangulation[[v,neighbours,
                          ConnectedNeighbours[triangulation,\
                              v][[Xunder,\
                            1]]]],
                        centre]]-
                      triangulation[[v,
                        centre]]],
                Apply[Crossing,
                  CircleParams[triangulation,v,{Xunder,Xover}]]},{ExtraGap[
                  triangulation,
                  triangulation[[v,neighbours,
                    ConnectedNeighbours[triangulation,
                        v][[Xunder,
                      1]]]],v,
                  graphicsParams[[gapParam]]],
                graphicsParams[[gapParam]]}}]],
        Apply[GetArc,
          Join[CircleParams[triangulation,v,
              Xunder],{{Apply[Crossing,
                  CircleParams[triangulation,v,{Xunder,Xover}]],
                triangulation[[v,r]]*
                  Sign[triangulation[[triangulation[[v,neighbours,ConnectedNeighbours[triangulation,v][[Xunder,2]]]],
                        centre]]-
                      triangulation[[v,
                        centre]]]},{graphicsParams[[gapParam]],
                ExtraGap[triangulation,
                  triangulation[[v,neighbours,
                    ConnectedNeighbours[triangulation,
                        v][[Xunder,
                      1]]]],v,
                  graphicsParams[[gapParam\
]]]}}]]],
      Apply[GetArc,
        Join[CircleParams[triangulation,v,
            Xover],{Map[
              triangulation[[v,r]]*
                  Sign[
                    triangulation[[#,
                        centre]]-
                      triangulation[[v,
                        centre]]]&,
              Map[triangulation[[v,
                    neighbours,#]]&,
                ConnectedNeighbours[triangulation,v][[
                  Xover]]]],
            Map[ExtraGap[triangulation,#,v,
                  graphicsParams[[
                    gapParam]]]&,
              Map[triangulation[[v,
                    neighbours,#]]&,
                ConnectedNeighbours[triangulation,v][[
                  Xover]]]]}]]}/;
    triangulation[[v,type]]=="X"

(* Draw an edge *)

GetGraphicsObjs[triangulation_,v_,
    graphicsParams_]:={Apply[GetArc,
        Join[CircleParams[triangulation,v,
            1],{Map[triangulation[[v,r]]*
                  Sign[triangulation[[#,
                        centre]]-
                      triangulation[[v,
                        centre]]]&,
              Map[triangulation[[v,
                    neighbours,#]]&,
                ConnectedNeighbours[triangulation,
                    v][[1]]]],
            Map[ExtraGap[triangulation,#,v,
                  graphicsParams[[1]]]&,
              Map[triangulation[[v,
                    neighbours,#]]&,
                ConnectedNeighbours[triangulation,
                    v][[1]]]]}]]}/;
    triangulation[[v,type]]=="e"

(* Draw a face -- actually no graphics are added here *)

GetGraphicsObjs[triangulation_,v_,graphicsParams_]:={}/;
    triangulation[[v,type]]=="f"

AddGraphicsObjs[triangulation_,graphicsParams_]:=
  AddField[triangulation,graphicsObjs,
    Table[GetGraphicsObjs[triangulation,v,graphicsParams],{v,
        Length[triangulation]}]]

Draw[triangulation_]:=
  Graphics[Flatten[
      FieldValues[triangulation,graphicsObjs]],{AspectRatio->1}]

(* Colours *)

colourList={RGBColor[0,0,0],RGBColor[1,0,0],RGBColor[0,1,0],RGBColor[1,1,0],
      RGBColor[0,0,1],RGBColor[0.5,0.25,0],RGBColor[1,0,1],
      RGBColor[1,0.5,0.5],RGBColor[1,0.5,0],RGBColor[0.5,0.5,0.5]};

AddColour[triangulation_,components_,colour_]:=
  Insert[Insert[triangulation,colour,
      Table[{components[[i,1]],
          graphicsObjs,components[[i,2]],
          1},{i,Length[components]}]],RGBColor[0,0,0],
    Table[{components[[i,1]],
        graphicsObjs,
        components[[i,2]],-1},{i,
        Length[components]}]]

ColourStrands[triangulation_,colouredStrands_]:=
  Fold[AddColour[#1,#2[[1]],#2[[2]]]&,triangulation,\
    Transpose[{ListStrands[triangulation],
        Take[Fold[
            Insert[#1,#2[[2]],#2[[1]]]&,\
            Select[colourList,
              FreeQ[If[Length[colouredStrands]==0,{},
                    Transpose[colouredStrands][[2]]],#]&],colouredStrands],\
          Length[ListStrands[triangulation]]]}]]

OuterFace="OuterFace";
(* Commented out by Dror: Gap="Gap"; *)
Colour="Colour";
StrandColour="StrandColour";
DebugDrawPD=False;

(* Dror: Add line and pd -> pd_PD *)
WhittenDrawPD[L_] := WhittenDrawPD[PD[L]]
WhittenDrawPD[L_,options_] := WhittenDrawPD[PD[L],options]
WhittenDrawPD[pd_PD]:=Module[{t},
    CreditMessage["WhittenDrawPD was written by Emily Redelmeier at the University of Toronto in the summers of 2003 and 2004."];
    t=AddPositionsBound[DefaultDirichlet[Triangulate[pd]]];
	If[DebugDrawPD,Print["madeit"]];
    t=PutInside[t,DefaultOuterFace[t]];t=Balance[t];
	If[DebugDrawPD,Print["beforegrapx:",t]];
    t=AddGraphicsObjs[t,{DefaultGap[t]}];t=ColourStrands[t,{}];Draw[t]];

WhittenDrawPD[pd_PD,options_]:=Module[{t},(optionsList=Map[Apply[List,#]&,options];
    t=AddPositionsBound[DefaultDirichlet[Triangulate[pd]]];
    t=PutInside[t,
        Which[Length[
              Select[optionsList,#[[1]]\
==OuterFace&]]==0,DefaultOuterFace[t],
          Depth[Select[
                  optionsList,#[[1]]\
==OuterFace&][[1,2]]]==1,
          Select[optionsList,#[[1]]\
==OuterFace&][[1,2]],True,
          GetOuterFace[t,
            Select[optionsList,#[[1]]\
==OuterFace&][[1,2]]]]];
    t=Balance[t];
    graphicsParams={If[
          Length[Select[
                optionsList,#[[1]]\
==Gap&]]==0,DefaultGap[t],
          Select[optionsList,#[[1]]\
==Gap&][[1,2]]]};
    t=AddGraphicsObjs[t,graphicsParams];
    t=If[Length[
            Select[optionsList,#[[1]]\
==Colour&]]==0,t,
        ColourStrands[t,
          If[Length[
                Select[optionsList,#[[1]]\
==StrandColour&]]==0,{},
            MapAt[Position[
                    ListStrands[t],{#+Length[pd],1}][[1,
                  1]]&,
              Select[optionsList,#[[1]]\
==StrandColour&][[1,2]],
              Table[{1,i},{i,
                  Length[Select[
                        optionsList,#[[1\
]]==StrandColour&][[1,\
                      2]]]}]]]]];Draw[t])
];

End[];

(* End source file src/WhittenDrawPD.m*)


(* ::Text:: *)
(*Some examination of the source code of DrawPD shows that the output is structured in the following way. *)
(**)
(*1) All outputs are circle arcs.*)
(*2) The first 3*crossings circles draw the crossings; they are arranged as a list of two circles (with color) for the undercrossing followed by a single circle (with color) for the overcrossing. These appear in the order that the crossings appear in the original PD.*)
(*3) These are followed by # arcs circles which draw the arcs. These appear in sorted order.*)
(**)
(*We can use this to write a wrapper which parses and decorates the DrawPD output, giving the arc numbers, orientations, component numbers and crossing numbers. We design a new primitive, OrientedCircle, which adds a simple triangular arrowhead denoting the orientation of the circle arc.*)


(* Build a triangle on the line segment AB *)
TriangleHead[A_,B_,C_] := Module[{ax,l}, 
ax = B - A;
l = (1/Sqrt[3])*{-ax[[2]],ax[[1]]};
Polygon[{{A+l,B,A-l,C}}]
];

(* Add an arrowhead at the head of a circle arc. *)
OrientedCircle[C_,r_,Theta_] :=Module[{ptA,ptB,ptC,size,rads,lpt,rpt,angles,orient},
size = 0.01; (* Try to draw the orientation triangle at this (absolute) size *)
rads = size / r; (* Which is this many radians *)
If[rads > (2/3)*Abs[Theta[[2]]-Theta[[1]]],Return[{{Circle[C,r,Theta]}}]]; (* arc too small 4 arrow*)
orient = Sign[Theta[[2]]-Theta[[1]]];
angles = {Theta[[2]] - orient*rads,Theta[[2]]+0.5*orient*rads,Theta[[2]] - orient*(2/3)*rads};
(* Now build the arrowhead and return the graphics as a list *)
{ptA,ptB,ptC} = C + r*{Cos[#],Sin[#]}&  /@ angles; (* Build the start, end of arrow *)
{Circle[C,r,Theta],TriangleHead[ptA,ptB,ptC]}
];

(* Add an arrow at the head of a line *)
OrientedLine[{A_,B_}] :=Module[{ptA,ptB,ptC,size,svals,rel},
size = 0.01; (* Try to draw the orientation triangle at this (absolute) size *)
rel = size/Norm[B-A];
If[rel > (2/3),Return[{{Line[{A,B}]}}]]; (* segment too small 4 arrow*)
svals = {1 - rel,1+0.5*rel,1 - (2/3)*rel};
(* Now build the arrowhead and return the graphics as a list *)
{ptA,ptB,ptC} = #*B + (1-#)*A&  /@ svals; (* Build the start, end of arrow *)
{Line[{A,B}],TriangleHead[ptA,ptB,ptC]}
];

LabelledCircle[C_,r_,Theta_,Lab_] :=Module[{pt},
{pt} = C + r*{Cos[#],Sin[#]}&  /@ {Mean[Theta]}; 
{Circle[C,r,Theta],Inset[Style[Lab,FontFamily->"Helvetica",FontWeight->Bold],pt,Background->White]}
];

LabelledLine[{A_,B_},Lab_] := Module[{pt},
{pt} = (1/2)*A + (1/2)*B;
{Line[{A,B}],Inset[Style[Lab,FontFamily->"Helvetica",FontWeight->Bold],pt,Background->White]}
];

ArcEndpoints[C_,r_,Theta_] := C + r*{Cos[#],Sin[#]} & /@ Theta;

OrientedCrossingGraphic[crossingArg_,cr_,styles_] := Module[{inOrient,ep,epd,C,Chead},
	
	C = crossingArg;
	Chead = If [Or[ToString[C[[1]]] == "Arc Deleted",ToString[C[[2]]] == "Arc Deleted"],
	   C[[1]],
	   ep = C[[1;;2]] /. {Circle[C_,r_,T_]:> ArcEndpoints[C,r,T],Line[{A_,B_}]:>{A,B}};
	   epd = Min[Norm[# - ep[[2,1]]],Norm[# - ep[[2,2]]]] & /@ ep[[1]];
	   (*Print["arcs:",C[[1;;2]]," ep:",ep," epd:",epd];*)
       (* epd is now a two-element list of distances from endpoints of arc crossing[[1]] to arc crossing[[2]]*)
	   If[TrueQ[Head[C[[1]]] == Circle], (* Sometimes the arc is a line *)
	       If [epd[[2]] < epd[[1]], 
		        OrientedCircle @@ C[[1]], (* the arc is correctly oriented *)
		        OrientedCircle @@ {C[[1,1]],C[[1,2]],Reverse[C[[1,3]]]}
	       ],
           If [epd[[2]] < epd[[1]],
	             OrientedLine @@ C[[1]], (* correct orientation *)
                 OrientedLine @@ Reverse[C[[1]]] 
           ]
        ]
	]; 

	(*Print[Chead];*)
	Graphics[{Flatten[Join[styles[[1]],{Chead}]],Join[styles[[2]],{C[[2]]}],Join[styles[[3]],{C[[3]]}]}] /. {{Style_,Circle[C_,r_,T_]} :> Tooltip[{Style,Circle[C,r,T]},ToString[cr]],Join[Style_,"Arc Deleted"] -> Null,{Style_,"Arc Deleted"} -> Null}];

CrossingGraphic[crossing_,cr_,styles_] := Module[{i},
	Graphics[Table[Tooltip[Join[styles[[i]],{crossing[[i]]}],ToString[cr]],{i,1,3}]] /. "Arc Deleted"->Null];

LabelledEdgeGraphic[edge_,ed_,style_] := 
	edge /. {Line[{A_,B_}] :> Graphics[Join[style,LabelledLine[{A,B},ToString[ed]]]],Circle[C_,r_,Th_] :> Graphics[Join[style,LabelledCircle[C,r,Th,ToString[ed]]]],"Arc Deleted"->Graphics[]};

EdgeGraphic[edge_,ed_,style_] := 
	edge /. {Line[{A_,B_}] :> Graphics[{style[[1]],Line[{A,B}]}],Circle[C_,r_,Th_] :> Graphics[{style[[1]],Circle[C,r,Th]}],"Arc Deleted"-> Graphics[]};


FaceScore[face_,cmps_] := Module[{facecmps},
	facecmps = Length[DeleteDuplicates[ComponentNumber[#,cmps]& /@ face]];
	4*((Length[cmps])/(facecmps)) + Length[face]
];

ComponentStyles = ComponentStyles;
Orientation = Orientation;
ArcLabels = ArcLabels;
CrossingLabels = CrossingLabels;
OuterFace = OuterFace;

Options[NiceDrawPD] = {ComponentStyles -> {{AbsoluteThickness[1]},{AbsoluteThickness[2]},{AbsoluteThickness[3]},{AbsoluteThickness[4]}},Gap -> 0.02,
	Orientation->True,EdgeLabels->True,OuterFace->Automatic};

NiceDrawPD[ThisLink_Link,opts:OptionsPattern[]] := NiceDrawPD[PD[ThisLink],opts];

NiceDrawPD[ThisPD_PD,OptionsPattern[]] := Module[

	{faces,outerFace,cmps,dpdout,crossings,edges,nCr,nE,edgeList,\
	 crgraphics,edgegraphics,styles,comps,crossingarcs,crossingstyles,edgestyles},

	If[TrueQ[OptionValue[OuterFace] == Automatic],
		faces = PDFaces[ThisPD];
		cmps = AllComponents[ThisPD];
		outerFace = Sort[faces,FaceScore[#1,cmps] > FaceScore[#2,cmps] &][[1]];
		Print["Setting outerFace automatically to ",outerFace];
	    ,
	    Print["Setting outerFace to ",OptionValue[OuterFace]];
		outerFace = OptionValue[OuterFace]];
	
	(* Sometimes the outerface call can cause DrawPD to crash. We time out after 20 seconds to avoid this. *)

	If[TrueQ[OptionValue[OuterFace] == "DrawPD"],
	    Print["Drawing with DrawPD default outer face."];
		dpdout = TimeConstrained[WhittenDrawPD[ThisPD,{Gap->OptionValue[Gap]}],40,Return["DrawPD can't draw this link"]];
	,
		dpdout = TimeConstrained[List @@ TimeConstrained[WhittenDrawPD[ThisPD,{Gap->OptionValue[Gap],OuterFace->outerFace}],40,Print["Timed out on outerface"]; outerFace = Reverse[outerFace]; WhittenDrawPD[ThisPD,{Gap->OptionValue[Gap],OuterFace->outerFace}]],40,Return["DrawPD can't draw this link"]];
	];

	nCr = Length[ThisPD];
	edgeList = Union[Flatten[List @@ (List @@@ ThisPD)]];
	nE = Length[edgeList];

	Print["DrawPD constructed drawing with ",Length[dpdout[[1]]]," arcs."];
	Print["PD has ",nCr," crossings and ",nE," edges."];

	styles = OptionValue[ComponentStyles];
	styles = PadRight[styles,NumComponents[ThisPD],styles]; (* Cyclically repeat styles as needed *)

    comps = AllComponents[ThisPD];

	crossings = Partition[Take[dpdout[[1]],3 nCr],3];
	edges = Take[dpdout[[1]],{3 nCr + 1,3 nCr + nE}];

	(* Now we need to build a list of styles for each arc by component number *)
	crossingarcs = List @@ Map[ComponentNumber[#,comps] &,((List @@ #[[{1,3,2}]]) & /@ ThisPD),{2}];
	crossingstyles = Map[styles[[#]] &,crossingarcs,1];
	(*Print[crossingarcs];
    Print[crossingstyles];*)

	(* We'll now apply styles for the edges*)
	edgestyles = styles[[#]] & /@ (ComponentNumber[#,comps] & /@ Range[nE]);
	
	crgraphics = If[OptionValue[Orientation],
		MapThread[OrientedCrossingGraphic,{crossings,List @@ ThisPD,crossingstyles}],
        MapThread[CrossingGraphic,{crossings,List @@ ThisPD,crossingstyles}]];
	
	edgegraphics = If[OptionValue[EdgeLabels],MapThread[LabelledEdgeGraphic,{edges,edgeList,edgestyles}],
	MapThread[EdgeGraphic,{edges,edgeList,edgestyles}]];

	(*Print[crgraphics,edgegraphics];*)

	Show[Join[crgraphics,edgegraphics]]

];

NiceDrawPD::usage = "
  NiceDrawPD[pd] takes the planar diagram description pd and creates a
  graphics object containing a picture of the knot.

  NiceDrawPD[pd,options], where options is a list of rules, allows the user
  to control some of the parameters.  

	Gap->g sets the gap around a crossing to g
	Orientation->True/False controls the display of orientation arrows
	EdgeLabels->True/False controls the display of labels for the edges
	ComponentStyles->{{a},{b},{c},...} sets the style for the components
		to {a}, {b}, {c} (repeated cyclically if there are more components
		than styles. For instance, {{Red},{Green},{Blue}} displays component 1 
		in Red, component 2 in Green and component 3 in Blue. The components
		are numbered in the natural order (different from the order in Skeleton),
		where the arcs in component 1 are numbered {1,...,k}, those in comp 2
		are numbered {k+1, ...., l} and so forth. The default sets the components
		to have increasingly dark line weights.
	OuterFace->{e_1,e_2,...,e_n} sets the face at infinity to a face which 
		has edges e_1, e_2, ..., e_n in the planar diagram description.  
"

NiceDrawPD::about = "
  NiceDrawPD is a wrapper which post-processes the output of the DrawPD function written by 
  Emily Redelmeier at the University of Toronto in the summers of 2003 and 2004. 
  NiceDrawPD was written in 2011 by Jason Cantarella at the University of Georgia.
"




(* ::InheritFromParent:: *)
(*"\n  NiceDrawPD[pd] takes the planar diagram description pd and creates a\n  graphics object containing a picture of the knot.\n\n  NiceDrawPD[pd,options], where options is a list of rules, allows the user\n  to control some of the parameters.  \n\n\tGap->g sets the gap around a crossing to g\n\tOrientation->True/False controls the display of orientation arrows\n\tEdgeLabels->True/False controls the display of labels for the edges\n\tComponentStyles->{{a},{b},{c},...} sets the style for the components\n\t\tto {a}, {b}, {c} (repeated cyclically if there are more components\n\t\tthan styles. For instance, {{Red},{Green},{Blue}} displays component 1 \n\t\tin Red, component 2 in Green and component 3 in Blue. The components\n\t\tare numbered in the natural order (different from the order in Skeleton),\n\t\twhere the arcs in component 1 are numbered {1,...,k}, those in comp 2\n\t\tare numbered {k+1, ...., l} and so forth.\n\tOuterFace->{e_1,e_2,...,e_n} sets the face at infinity to a face which \n\t\thas edges e_1, e_2, ..., e_n in the planar diagram description.  \n"*)


(* ::Section:: *)
(*SpliceTwists, FramedDouble, and ZeroFrame*)


(* ::Text:: *)
(*After constructing the Blackboard Double of a given component in a link, the next natural goal is to construct doubles of various framings (that is to say, we want to prescribe the linking number of the two components). To do so, there are several substeps:*)
(**)
(*FindFourGonSides[a,b,ThisPD] looks for a four-gon among the faces of the given PD code so that *)
(*	*)
(*	1. Two opposite sides, arcs SideA, SideB, belong to components a and b.*)
(*	2. The crossings differ from under to over on the two sides of the crossing. *)
(*	2. The orientations of these sides agree.*)
(*	3. One side is part of a face with as many sides as possible (hence likely to be external in NiceDrawPD). *)
(*	    The total number of arcs in faces incident to a given arc is given by ArcScore[arc,ThisPD].*)
(*	4. Returns {SideA,SideB}.*)
(**)
(*Assuming that such a four-gon can be found between the given component and the new component in the Blackboard framing, we then have *)
(**)
(*SpliceTwists[n,{SideA,SideB}, n, ThisPD] which converts ThisPD to an extended PD, then twists SideA and SideB together with n (positive or negative) twists, converting back to PD and returning a PD at the end of the process.*)
(**)
(*FramedDouble organizes the work, calculating the linking number of the new component in the BlackboardDouble, then using the functions above to splice in twists as require to get the desired linking number.*)
(**)
(*ZeroFramedDouble is a convenience function which calls FramedDouble to frame a component with linking number zero.*)


ArcScore[arc_,ThisPD_] := Length[Flatten[Select[PDFaces[ThisPD],MemberQ[#,arc]&]]];

FindFourGonSides[a_,b_,ThisPD_PD] := Module[{faces,arcpairs,cpts},
	faces = Select[OrientedPDFaces[ThisPD],Length[#] == 4 &]; (* fourgons *)
	faces = Join[faces,RotateLeft /@ faces]; (* now we need only consider 1,3 pairs *)
	faces = Select[faces,#[[1,2]] != #[[3,2]] &]; (* orientations agree *)
	(* We can now throw out orientation string and the 2nd and 4th entry *)
	arcpairs = (faces[[;;,{1,3}]]) /. {n_,_String} :> n;
	(*Print["All pairs of edges opposite in a fourgon:",arcpairs];*)
	(* decorate with component numbers *)
	cpts = AllComponents[ThisPD];	
	(*Print[cpts];*)
	arcpairs = arcpairs /. n_Integer :> {n,ComponentNumber[n,cpts]};
	(* search for pairs with the right component numbers *)	
	(*Print["With components:",arcpairs];*)
	arcpairs = Join[Cases[arcpairs,{{c_,a},{d_,b}} :> {c,d}],
	                Cases[arcpairs,{{c_,b},{d_,a}} :> {d,c}]];
	(*Print["On components ",a," and ",b,":",arcpairs];*)
	(* search for pairs with over/under correct. We'd like both arc to pass over *)
	(* and under the next strands together. *)
	arcpairs = Select[arcpairs,\
		ArcAlternatingQ[#[[1]],ThisPD] && ArcAlternatingQ[#[[2]],ThisPD] && \
	    !Xor[ArcOverAtHeadQ[#[[1]],ThisPD],ArcOverAtHeadQ[#[[2]],ThisPD]] &];
	(* sort by ArcScore *)
	(*Print[arcpairs];*)
	Sort[arcpairs,Total[ArcScore /@ #1] > Total[ArcScore /@ #2]&][[1]]
];

TwistArc[Primary_,Index_,nc_] := {Primary,Index,{nc+Index-1,nc+Index}}; (* An intermediate arc *)

TwistCrossing[SideArcs_,Index_,Sign_,nc_] := Module[{minusVersion,cr},
	
	minusVersion = X[TwistArc[SideArcs[[2]],Index,nc],TwistArc[SideArcs[[1]],Index,nc],\
		            TwistArc[SideArcs[[2]],Index+1,nc],TwistArc[SideArcs[[1]],Index+1,nc]];

	minusVersion = If[OddQ[Index],minusVersion,minusVersion[[{2,1,4,3}]]];  

	cr = If[Sign>0,RotateLeft[minusVersion],minusVersion];
	(*Print["TwistCrossing:",Index," (",Sign,")"];*)
	cr
];

SpliceTwists[n_,SideArcs_,ThisPD_PD] := Module[{ThisEPD,i,nc,TCs,heads,tails},
	heads = ArcHead[#,ThisPD]& /@ SideArcs; (*Print["heads:",heads];*)
	tails = ArcTail[#,ThisPD]& /@ SideArcs; (*Print["tails:",tails];*)
	(* We begin by adding 2*Abs[n] twist crossings of the correct sign *)
	ThisEPD = PDtoEPD[ThisPD]; nc = Length[ThisPD];
	TCs = Table[TwistCrossing[SideArcs,i,n,nc],{i,1,2*Abs[n]}]; TCs[[0]] = EPD;
	(* The first and last of the twist crossings need to be fixed manually *)
	TCs[[1]] = TCs[[1]] //. {{SideArcs[[1]],1,{t_,h_}} :> {SideArcs[[1]],1,{tails[[1]],h}},
        {SideArcs[[2]],1,{t_,h_}} :> {SideArcs[[2]],1,{tails[[2]],h}}};
	TCs[[Length[TCs]]] = TCs[[Length[TCs]]] //. {{SideArcs[[1]],2*Abs[n]+1,{t_,h_}} :> {SideArcs[[1]],2*Abs[n]+1,{t,heads[[1]]}},
		{SideArcs[[2]],2*Abs[n]+1,{t_,h_}} :> {SideArcs[[2]],2*Abs[n]+1,{t,heads[[2]]}}};
	(* Now we're ready to put these into the EPD *)
	ThisEPD = Join[ThisEPD,TCs];
	(* Now the first and last of the new arcs refer (incorrectly) to crossings nc and nc + 2 Abs[n] + 1 *)
	ThisEPD[[tails[[1]]]] = ThisEPD[[tails[[1]]]] /. {SideArcs[[1]],1,{t_,h_}} :> {SideArcs[[1]],1,{t,nc+1}};
	ThisEPD[[tails[[2]]]] = ThisEPD[[tails[[2]]]] /. {SideArcs[[2]],1,{t_,h_}} :> {SideArcs[[2]],1,{t,nc+1}};
	(* We also need to alter the old heads of the SideArcs to reflect the subdivisions *)
	ThisEPD[[heads[[1]]]] = ThisEPD[[heads[[1]]]] /. {SideArcs[[1]],1,{t_,h_}} :> {SideArcs[[1]],2*Abs[n]+1,{nc+2*Abs[n],h}};
	ThisEPD[[heads[[2]]]] = ThisEPD[[heads[[2]]]] /. {SideArcs[[2]],1,{t_,h_}} :> {SideArcs[[2]],2*Abs[n]+1,{nc+2*Abs[n],h}};
	(*Print[ThisEPD];*)
	EPDtoPD[ThisEPD]
];

ComponentwiseLinkingNumber[a_,b_,ThisPD_PD] := Module[{cmps,crs},
	cmps = AllComponents[ThisPD];
	crs = Position[Map[ComponentNumber[#,cmps]&,PDtoList[ThisPD],{2}],{a,b,a,b}|{b,a,b,a}];
	Total[CrossingSign[ThisPD[[#[[1]]]],ThisPD]& /@ crs]/2
];

FramedDouble[DesiredLK_,Component_,ThisPD_] := Module[{newPD,currentLK,sideArcs},
	newPD = BlackboardDouble[Component,ThisPD];
	currentLK = ComponentwiseLinkingNumber[Component,NumComponents[newPD],newPD];
	Print["FramedDouble:CurrentLk:",currentLK," DesiredLk:",DesiredLK];
	If[currentLK == DesiredLK,
		newPD,
		sideArcs = FindFourGonSides[Component,NumComponents[newPD],newPD];
		If[Length[sideArcs] > 0,
			SpliceTwists[DesiredLK - currentLK,sideArcs,newPD],
	        Print["FramedDouble: Component ",Component," doesn't have an acceptable fourgon with component," NumComponents[newPD]," in double."];
		]
	]
];

ZeroFramedDouble[Component_,ThisPD_] := FramedDouble[0,Component,ThisPD];



(* ::Section:: *)
(*Mirror Symmetry*)


(* ::Text:: *)
(*Now we try to come up with a transformation on PD codes. Suppose we have a crossing code. We know that it's in the form *)
(**)
(*X[k,m,k+1,m+1] or X[k,m+1,k+1,m]*)
(**)
(*or in the form*)
(**)
(*X[k,m,k+1,1] or X[k,1,k+1,m]*)
(**)
(*if the crossing connects the last arc to the first.*)
(**)
(*If we switch the crossing, we should make a cyclic permutation to put the m term first (cyclic because we aren't changing orientation in the plane, and putting the m term first because it's now the undercrossing). This works out to the following Mathematica code.*)


SwitchCrossing[x_,PD_] := Module[{NextGuy,Pos,thisLoop},

If[x[[4]] == NextArc[x[[2]],PD],RotateLeft[x,1],RotateRight[x,1]]

];


(* ::Text:: *)
(*We can now write the method that performs the mirror operation on a PD code. This is a one-liner:*)


MirrorPD[sign_,ThisPD_PD] := If [sign == 1,ThisPD,PD @@ (SwitchCrossing[#,ThisPD] & /@ (List @@ ThisPD))];


(* ::Section:: *)
(*Reversing the Orientation on a Component*)


(* ::Text:: *)
(*To reverse the orientation on a component in a PD code, we need to first make a transformation which relabels the arcs of the component in reverse order. We then have the problem that the arcs in the PD code might be in the wrong order if the newly reversed component is an undercrossing. To fix this, we write the following code:*)


FixCrossing[x_,PD_] := Module[{},

If[x[[1]] == NextArc[x[[3]],PD],RotateLeft[x,2],x]

];


(* ::Text:: *)
(*And now the code for reversing a component. To make later code make sense, we think of this as applying an orientation (+1 or -1) to the component. *)


OrientComponent[ThisPD_PD,{sign_,n_}] := Module[{comp,revcomp,rules,i},
(*Print["RunningPD"];*)
If [sign == 1, Return[ThisPD],
	comp = AllComponents[ThisPD][[n]];
	revcomp = Reverse[comp];
	rules = Table[comp[[i]] -> revcomp[[i]],{i,1,Length[comp]}];
	(*Print["comp:",comp,"revcomp:",revcomp,"rules:",rules];*)
	newPD = ThisPD /. rules;
	(*Print[newPD];*)
	FixCrossing[#,ThisPD]& /@ newPD ]
];

OrientComponent[L_Link,{sign_,n_}] := OrientComponent[PD[L],{sign,n}];


(* ::Section:: *)
(*Permuting Components. *)


(* ::Text:: *)
(*Using transformation rules, permuting the components of a PD is as simple as figuring out how we want the relabelling to go. Remember that component Perm[[i]] of the original PD will be component i of the permuted PD. That means that we should start by taking AllComponents, which gives us a list of the arcs in the PD, broken into sublists for each component, then permute the sublists by Perm. If we flatten the resulting list of lists, we have the order in which these arcs should appear in the NEW PD. This means that we just have to make a list of numbers {1, ...,  numArcs*)


PermutePD[ThisPD_PD,Perm_List] := Module[{AllComps,PermComps,Rules},

	AllComps = AllComponents[ThisPD];
	If[Length[Perm] != Length[AllComps],Print["Permutation ",Perm," has the wrong number of entries for ",Length[AllComps]," component PD."]; Return[]];
	PermComps = Flatten[AllComps[[Perm]]]; (* Permuted list of components *)
	Rules = Dispatch[Table[PermComps[[i]] -> i,{i,1,Length[Flatten[AllComps]]}]];
	
	ThisPD /. Rules

];	


(* ::Section:: *)
(*Putting it Together: Applying a Whitten Group Element to a Given Link*)


(* ::Text:: *)
(*We now want to apply a given Whitten group element to a given link. We decide that a Whitten group element will look like this in Mathematica:*)
(**)
(*{+1,{+-1, ..., +-1}, p}*)
(**)
(*where p = {p1, ..., pn} is a list of the n components of the link. The permutation will be in the form*)
(**)
(*1 -> p1, 2 -> p2, 3 -> p3 .... n -> pn. We follow the Whitten group rule of "permute first, then reverse". Remember that we apply "epsilon[i]" (the guy in position [i] of the list of reversals) NOT to element i in the skeleton, but to the permuted index element p[i]. The only clever bit of code here is the Fold, which applies OrientComponent sequentially to newL for each sign in the list Rev and returns the overall result. *)
(**)
(*The result applies all three procedures we just wrote: *)
(**)
(*	Permute (on the inside, with PermutePD)*)
(*	Reverse component orientations (at the middle level, by Folding OrientComponent onto the permuted PD)*)
(*	Mirror (on the outside, using MirrorPD)*)


ApplyWhitten[W_,L_] := Module[{Mir,Rev,P},
	{Mir, Rev, P} = W;
	MirrorPD[Mir,Fold[OrientComponent,PermutePD[L,P],Partition[Riffle[Rev,Range[Length[Rev]]],2]]]
];


(* ::Section::Closed:: *)
(*Listing Whitten Group Elements*)


(* ::Text:: *)
(*The joy of the Whitten group is that while the group operation is weird, the group structure is basically trivial. So we can list all of the various Whitten group elements without too much work.*)


WhittenGroup[n_] := Module[{Revs,Perms,G,Orientations}, 

Revs = If[n ==1 , {{1},{-1}},

If [n == 2, Flatten[Outer[List,{1,-1},{1,-1}],1],

If [n == 3, Flatten[Outer[List,{1,-1},{1,-1},{1,-1}],2],

If [n == 4, Flatten[Outer[List,{1,-1},{1,-1},{1,-1},{1,-1}],3], 

If [n == 5, Flatten[Outer[List,{1,-1},{1,-1},{1,-1},{1,-1},{1,-1}],4],

If [n == 6, Flatten[Outer[List,{1,-1},{1,-1},{1,-1},{1,-1},{1,-1},{1,-1}],5],

{"error"}]]]]]];

Perms =Permutations[Range[n]];

Flatten[Outer[List,{-1,1},Revs,Perms,1],2]

];


WhittenIdentity[n_] := Module[{i},{1,Table[1,{i,1,n}],Range[n]}];


Needs["Combinatorica`"]

StandardOrderCycle[C_] := Module[{Least,Pos},
	
	Least = Sort[C][[1]];
	Pos = Position[C,Least][[1]];
	RotateLeft[C,Pos-1]

];
	
WhittenPrettyPrint[W_] :=  Module[{Mir,Rev,P,CycleForm,i,StringCycles,StringPerm,n},

{Mir, Rev, P} = W;
n = Length[Rev];
CycleForm = Select[ToCycles[P],Length[#1] > 1 &];
CycleForm = StandardOrderCycle /@ CycleForm;

If[Length[CycleForm] == 0,StringPerm="e",

	StringCycles = StringJoin[Table[ToString[StringForm["``",#1[[i]]]],{i,1,Length[#1]}]] & /@ CycleForm;
	StringCycles = StringJoin["(",#1,")"] & /@ StringCycles;
	StringPerm = StringJoin[StringCycles];
];

(* Now we replace e, Pure Mirror, Pure Invert, Pure Exchange, and Mirror Exchange with abbreviations *)

If[Mir == 1 && Rev == Table[1,{i,1,n}] && StringPerm=="e","e",
	If[Mir == 1 && Rev == Table[-1,{i,1,n}] && StringPerm=="e","PI",
		If[Mir == -1 && Rev == Table[1,{i,1,n}] && StringPerm=="e","PM",
			If[Mir == 1 && Rev == Table[1,{i,1,n}],"PE " <> StringPerm, 
				If[Mir == -1 && Rev == Table[1,{i,1,n}],"ME "<> StringPerm,
					"(" <> ToString[Mir] <> ",(" <> StringJoin[Riffle[ToString /@ Rev,","]] <> ")," <> StringPerm <> ")" 
			]]]]]
];

WhittenPrettyPrint::usage = "
	WhittenPrettyPrint formats Whitten group elements and replaces certain common elements with 
	abbreviated names. The abbreviations include

		e - the identity element {1,{1,...,1},{1,2,...,n}}
		PI - the \" pure inversion \" {1,{-1,...,-1},{1,2,...,n}}
		PM - the \" pure mirror \" {-1,{1,...,1},{1,2,...,n}}
		PE (x) - the \" pure exchange \" {1,{1,...,1},{x}}
		ME (x) - the \" mirror exchange \" {-1,{1,....,1},{x}}

	";

WhittenTablePrint[L_,n_] := Module[{LStrings,Rows,Ands,i},

	LStrings = Sort[WhittenPrettyPrint /@ L]; (* Apply WPP to everything *)
	Print[LStrings];
	Rows = Partition[LStrings,n];  (* split into sublists of length n *)

	StringJoin[Table[StringJoin[Riffle[Rows[[i]]," & "],"\\\\ \n"],{i,1,Length[Rows]}]]

];


(* ::Subtitle:: *)
(*Additional Knot Invariants*)


(* ::Text:: *)
(*We now add a few knot invariants to KnotTheory's basic collection.*)


(* ::Section:: *)
(*Linking Numbers and the Linking Matrix*)


(* ::Text:: *)
(*We now come up with some code to figure out the sign of a crossing. *)


CrossingSign[x_,PD_] := If[x[[2]] == NextArc[x[[4]],PD],+1,-1];


(* ::Text:: *)
(*We also want to figure out which components are involved in a crossing. *)


CrossingComponents[x_,PD_PD] := {ComponentNumber[x[[1]],PD],ComponentNumber[x[[2]],PD]};
CrossingComponents[x_,AllComps_List] := {ComponentNumber[x[[1]],AllComps],ComponentNumber[x[[2]],AllComps]};


(* ::Text:: *)
(*We are now in a position to work out the linking matrix by computing the sign of each component and adding the appropriate term to each entry of the linking matrix. We note that the linking matrix should be symmetric, but we're only adding crossingsign in one of the two places it ought to go. *)
(**)
(*Luckily, we can solve the problem by adding the transpose of the matrix to the matrix. We don't keep track of the main diagonal, because it is only an invariant for alternating knots. The result is a matrix giving the pairwise linking numbers of the components of L.*)


LinkingMatrix[L_PD] := Module[{LinkM,i,j,k,Cross,Sign,AllComps},

LinkM = Table[0,{NumComponents[L]},{NumComponents[L]}];
AllComps = AllComponents[L];

Do[  
  Cross = L[[k]];
  {i,j} = CrossingComponents[Cross,AllComps];
  LinkM[[i,j]] = LinkM[[i,j]] + CrossingSign[Cross,L];
,
{k,1,Length[PD[L]]}];

LinkM = LinkM + Transpose[LinkM];
LinkM = LinkM - DiagonalMatrix[Diagonal[LinkM]];
LinkM = (1/2)*LinkM;

LinkM
];


(* ::Section:: *)
(*The (Corrected) Multivariable Alexander Polynomial*)


(* ::Text:: *)
(*KnotTheory gives some code for determining the Multivariable Alexander Polynomial of a knot or link. The ordering of the variables in this polynomial is given by the ordering in Skeleton. Unfortunately, as we've seen above, this ordering is somewhat unpredictable, and requires the application of the correction permutation in order to be consistent. We now give a wrapper for the MVA which applies the correction. *)
(**)
(*In the end, while we leave this code in, we don't use it. *)


CorrectedMultivariableAlexander[{L_,P_},tt_] := 
  MultivariableAlexander[L][tt] /. Table[tt[i] -> tt[P[[i]]],{i,1,NumComponents[L]}];


(* ::Section:: *)
(*Knot types of different components*)


(* ::Text::Italic:: *)
(*One of the more effective tools for detecting exchange symmetries should be determining when the two components to be exchanged have different knot types. Here is some code for that. The component numbering here is in Natural order.*)


WhittenSubLink[pd_PD, js_List] := Module[
  {k, t0, t, t1, t2, S, P},
  t0 = Flatten[List @@@ AllComponents[pd][[js]]];
  t = pd /. x_X :> Select[x, MemberQ[t0, #] &];
  t = DeleteCases[t, X[]];
  k = 1;
  While[
   k <= Length[t],
   If[ Length[t[[k]]] < 4, 
     t = Delete[t, k] /. (Rule @@ t[[k]]), ++k];
   ];
  t1 = List @@ Union @@ t;
  t2 = Thread[(t1) -> Range[Length[t1]]];
  S = t /. t2;
  P = If[S != PD[] && Length[S] >= 3, S, PD[Knot[0, 1]], S]
  ];
WhittenSubLink[pd_PD, j_] := WhittenSubLink[pd, {j}];
WhittenSubLink[L_, js_] := WhittenSubLink[PD[L], js];

ComponentsOk[L_,M_] := Module[{LComps,MComps,QQ,AllMatch},

	LComps = Jones[#1][QQ]& /@ Table[WhittenSubLink[L,i],{i,1,NumComponents[L]}];
	MComps = Jones[#1][QQ]& /@ Table[WhittenSubLink[M,i],{i,1,NumComponents[M]}];

	AllMatch = And @@ Table[TrueQ[LComps[[i]] == MComps[[i]]],{i,1,NumComponents[L]}];
	AllMatch
];
 
	


(* ::Subsection:: *)
(**)


(* ::Section:: *)
(*The Satellite Lemma.*)


(* ::Text:: *)
(*One way to determine whether two components of a link can be exchanged is to double the components (preserving framing) and then compare two versions of the link with different components doubled. *)
(**)
(*The idea is that we call FrameComponents, which applies a framing to each of the components in the link which are involved in the permutation to tell them apart. FrameComponents does not frame a component when the associated index is zero, but otherwise frames the components with linking number given by the index.*)
(**)


FrameComponents[ThisPD_PD,Indices_] := Module[{workingPD,i},
	Print["Trying to frame with",Indices];
	For[workingPD = ThisPD; i=1,i<=Length[Indices],i++,
		workingPD = If[Indices[[i]] != 0,FramedDouble[Indices[[i]],i,workingPD],workingPD]];
	workingPD
];

SatelliteLemmaTestQ[a_,b_,ThisPD_PD] := Module[{AFrame,AJones,BFrame,BJones,QQ},

	AFrame = FramedDouble[1,a,ThisPD];
	BFrame = FramedDouble[1,b,ThisPD];
	
	AJones = Jones[AFrame][QQ];
	BJones = Jones[BFrame][QQ];
	
	Print["AJones:",AJones /. QQ -> "z"];
	Print["BJones:",BJones /. QQ -> "z"];

	TrueQ[AJones == BJones]
];

SatelliteTest[L_Link,M_Link] := SatelliteTest[PD[L],PD[M]];

SatelliteTest[L_PD,M_PD] := Module[{LFramed,MFramed,LJones,MJones,nc,LHOMFLYPT,MHOMFLYPT,AA,ZZ},

	(*Print["Got to sat test"];*)

	nc = NumComponents[L];
	LFramed = FrameComponents[L,Range[nc]-1];
	MFramed = FrameComponents[M,Range[nc]-1];

		(*Print[LFramed,MFramed];*)

	LJones = Jones[LFramed][QQ];
	MJones = Jones[MFramed][QQ];
	
	Print["L Jones:",LJones /. QQ -> "z"];
	Print["M Jones:",MJones /. QQ -> "z"];

	TrueQ[LJones == MJones]

	(*If[!TrueQ[LJones == MJones],TrueQ[LJones == MJones],

	  Print["Checking HOMFLYPT polynomial..."];
	
	  LHOMFLYPT = HOMFLYPT[LFramed][AA,ZZ];
	  MHOMFLYPT = HOMFLYPT[MFramed][AA,ZZ];

	  Print["L HOMFLYPT:",LHOMFLYPT /. {ZZ -> "z",AA->"a"}];
	  Print["M HOMFLYPT:",MHOMFLYPT /. {ZZ -> "z",AA->"a"}];

	  TrueQ[LHOMFLYPT == MHOMFLYPT]
	] *)
];		

SatellitePureExchanges[ThisLink_Link] := SatellitePureExchanges[PD[ThisLink]];

SatellitePureExchanges[ThisPD_PD] := Module[{cases,nc},
	cases = Subsets[Range[NumComponents[ThisPD]],{2}];
	Print["SatellitePureExchanges:Testing ",Length[cases]," exchanges"];
	{ToString[#[[1]]] <> "<->" <> ToString[#[[2]]],SatelliteLemmaTestQ[#[[1]],#[[2]],ThisPD]}& /@ cases
];
	


(* ::Section:: *)
(*Comparing links *)


(* ::Text:: *)
(*We now give some code for determining when two links are equivalent using the various invariants provided by KnotTheory, along with the two new invariants (LinkingMatrix and CorrectedMultivariableAlexander) that we defined in the last section. This code sets the global "WhyDifferent".*)


DifferentLinkQ[L_,M_] := Module[{QQ,a,z,TT},
 
	If [!TrueQ[LinkingMatrix[L] ==  LinkingMatrix[M]], WhyDifferent = "LinkingMatrix";True,
		If [!TrueQ[Jones[L][QQ] == Jones[M][QQ]],  WhyDifferent = "Jones"; True,
			If[!TrueQ[HOMFLYPT[L][a,z] == HOMFLYPT[M][a,z]],WhyDifferent = "HOMFLYPT";  True,
				If[!ComponentsOk[L,M],
					WhyDifferent = "Components Have Different Types"; True, 
				  If[!SatelliteTest[L,M], WhyDifferent = "Satellite Lemma"; True,False]]]]]

];


(* ::Section:: *)
(*Computing Whitten symmetry subgroups.*)


(* ::Text:: *)
(*We now do our best to compute Whitten symmetry subgroups for a given link by applying the various elements of the Whitten group and checking whether the resulting links are different using the DifferentLinkQ test given above.*)


SigmaGroup[L_Link] := SigmaGroup[PD[L]];

SigmaGroup[L_PD] := Module[{WhittenG,nComp,Possibles,Diff,BaseLink,NewLink,i,StringF},

nComp = NumComponents[L];
WhittenG = WhittenGroup[nComp];
Possibles = {};	

RuledOut = {}; (* Clear the global RuledOut and the individual lists of ruled out guys*)
RuledOutByLinkingMatrix = {}; 
RuledOutByJones = {};
RuledOutByHOMFLYPT = {};
RuledOutByComponents = {};
RuledOutBySatelliteLemma = {};

Print["Trying ",Length[WhittenG], " elements from ",nComp," component Whitten group.\n"];

BaseLink = ApplyWhitten[WhittenIdentity[nComp],L];
Print[BaseLink];

Do[
	StringF = IdentifyWhittenElement[WhittenG[[i]]];
	If [Or[StringF == "e",StringF == "PI"], (* These can't be ruled out *)
	   PrependTo[Possibles,WhittenG[[i]]]; Print[WhittenG[[i]], " possible "],

	NewLink = ApplyWhitten[WhittenG[[i]],L];
    Diff =  DifferentLinkQ[BaseLink,NewLink];

    If [Diff == True, 
		PrependTo[RuledOut,{WhittenG[[i]],WhyDifferent}]; Print[WhittenG[[i]]," ruled out ", WhyDifferent]; , 
		PrependTo[Possibles,WhittenG[[i]]]; Print[WhittenG[[i]], " possible "]; 
	];

	If[WhyDifferent == "LinkingMatrix",PrependTo[RuledOutByLinkingMatrix,WhittenG[[i]]]];
	If[WhyDifferent == "Jones",PrependTo[RuledOutByJones,WhittenG[[i]]]];
	If[WhyDifferent == "HOMFLYPT",PrependTo[RuledOutByHOMFLYPT,WhittenG[[i]]]];
	If[WhyDifferent == "Components Have Different Types",PrependTo[RuledOutByComponents,WhittenG[[i]]]];
	If[WhyDifferent == "Satellite Lemma",PrependTo[RuledOutBySatelliteLemma,WhittenG[[i]]]];
	];

,{i,1,Length[WhittenG]}];

Print["Ruled out: ",RuledOut // ColumnForm];
Print["Sigma appears to be a ",Length[Possibles]," element subgroup:"];

Possibles

];


(* ::Section:: *)
(*Processing Whitten groups and elements*)


(* ::Text:: *)
(*This section contains some helper functions for manipulating Whitten group elements. First, we define a natural order on Whitten group elements: we go in dictionary order on the mirror and reversals, counting +1 as < -1. *)


WhittenEltCompare[A_,B_] := Module[{mrA,mrB},
	mrA = StringJoin[Join[{A[[1]]},A[[2]]] /. {1 -> "a",-1 -> "b"}];
	mrB = StringJoin[Join[{B[[1]]},B[[2]]] /. {1 -> "a",-1 -> "b"}];
	If[Order[mrA,mrB] != 0, Order[mrA,mrB] == 1,
		OrderedQ[{A[[3]],B[[3]]}]]
];		


gamma1Rules = Dispatch[{
 {{1,{1},{1}}} -> "No Symmetry",
 {{1,{-1},{1}},{1,{1},{1}}} -> "Reversible",
 {{-1,{1},{1}},{1,{1},{1}}} -> "+Amphichiral",
 {{-1,{-1},{1}},{1,{1},{1}}} -> "-Amphichiral",
 {{-1,{-1},{1}},{-1,{1},{1}},{1,{-1},{1}},{1,{1},{1}}} -> "Full Symmetry"}];

gamma2Rules = Dispatch[{	
 {{1, {1, 1}, {1, 2}}} -> "No Symmetry", 
 {{1, {-1, -1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,1}", 
 {{-1, {1, 1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,2}", 
 {{-1, {-1, -1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,3}", 
 {{1, {-1, 1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,4}", 
 {{1, {1, -1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,4}", 
 {{1, {1, 1}, {1, 2}}, {1, {1, 1}, {2, 1}}} -> "Sigma_{2,5}", 
 {{1, {-1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,5}", 
 {{-1, {1, 1}, {2, 1}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,6}", 
 {{-1, {-1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,6}", 
 {{-1, {-1, 1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,7}", 
 {{-1, {1, -1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{2,7}", 
 {{1, {-1, -1}, {1, 2}}, {1, {-1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}, 
   {1, {1, 1}, {2, 1}}} -> "Sigma_{4,1}", 
 {{1, {-1, -1}, {1, 2}}, {1, {-1, 1}, {1, 2}}, {1, {1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,2}", 
 {{-1, {-1, 1}, {1, 2}}, {-1, {1, -1}, {1, 2}}, {1, {-1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,3}", 
 {{1, {-1, -1}, {1, 2}}, {1, {-1, 1}, {2, 1}}, {1, {1, -1}, {2, 1}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,4}", 
 {{-1, {-1, 1}, {2, 1}}, {-1, {1, -1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,5}", 
 {{-1, {-1, 1}, {1, 2}}, {-1, {1, 1}, {1, 2}}, {1, {-1, 1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,6}", 
 {{-1, {1, -1}, {1, 2}}, {-1, {1, 1}, {1, 2}}, {1, {1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,6}", 
 {{-1, {1, 1}, {1, 2}}, {-1, {1, 1}, {2, 1}}, {1, {1, 1}, {1, 2}}, 
   {1, {1, 1}, {2, 1}}} -> "Sigma_{4,8}", 
 {{-1, {-1, -1}, {2, 1}}, {-1, {1, 1}, {1, 2}}, {1, {-1, -1}, {2, 1}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,8}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {1, 1}, {1, 2}}, {1, {-1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,9}", 
 {{-1, {-1, -1}, {2, 1}}, {-1, {1, 1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,11}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {-1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}, 
   {1, {1, 1}, {2, 1}}} -> "Sigma_{4,12}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {1, 1}, {2, 1}}, {1, {-1, -1}, {2, 1}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,12}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {-1, 1}, {1, 2}}, {1, {1, -1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,14}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {1, -1}, {1, 2}}, {1, {-1, 1}, {1, 2}}, 
   {1, {1, 1}, {1, 2}}} -> "Sigma_{4,14}", 
 {{1, {-1, -1}, {1, 2}}, {1, {-1, -1}, {2, 1}}, {1, {-1, 1}, {1, 2}}, 
   {1, {-1, 1}, {2, 1}}, {1, {1, -1}, {1, 2}}, {1, {1, -1}, {2, 1}}, 
   {1, {1, 1}, {1, 2}}, {1, {1, 1}, {2, 1}}} -> "Sigma_{8,1}", 
 {{-1, {-1, 1}, {1, 2}}, {-1, {-1, 1}, {2, 1}}, {-1, {1, -1}, {1, 2}}, 
   {-1, {1, -1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, {1, {-1, -1}, {2, 1}}, 
   {1, {1, 1}, {1, 2}}, {1, {1, 1}, {2, 1}}} -> "Sigma_{8,2}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {-1, 1}, {1, 2}}, {-1, {1, -1}, {1, 2}}, 
   {-1, {1, 1}, {1, 2}}, {1, {-1, -1}, {1, 2}}, {1, {-1, 1}, {1, 2}}, 
   {1, {1, -1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{8,3}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {-1, 1}, {2, 1}}, {-1, {1, -1}, {2, 1}}, 
   {-1, {1, 1}, {1, 2}}, {1, {-1, -1}, {1, 2}}, {1, {-1, 1}, {2, 1}}, 
   {1, {1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{8,4}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {-1, -1}, {2, 1}}, {-1, {1, 1}, {1, 2}}, 
   {-1, {1, 1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, {1, {-1, -1}, {2, 1}}, 
   {1, {1, 1}, {1, 2}}, {1, {1, 1}, {2, 1}}} -> "Sigma_{8,5}", 
 {{-1, {-1, -1}, {2, 1}}, {-1, {-1, 1}, {1, 2}}, {-1, {1, -1}, {1, 2}}, 
   {-1, {1, 1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, {1, {-1, 1}, {2, 1}}, 
   {1, {1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{8,6}", 
 {{-1, {-1, -1}, {2, 1}}, {-1, {-1, 1}, {2, 1}}, {-1, {1, -1}, {2, 1}}, 
   {-1, {1, 1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, {1, {-1, 1}, {1, 2}}, 
   {1, {1, -1}, {1, 2}}, {1, {1, 1}, {1, 2}}} -> "Sigma_{8,7}", 
 {{-1, {-1, -1}, {1, 2}}, {-1, {-1, -1}, {2, 1}}, {-1, {-1, 1}, {1, 2}}, 
   {-1, {-1, 1}, {2, 1}}, {-1, {1, -1}, {1, 2}}, {-1, {1, -1}, {2, 1}}, 
   {-1, {1, 1}, {1, 2}}, {-1, {1, 1}, {2, 1}}, {1, {-1, -1}, {1, 2}}, 
   {1, {-1, -1}, {2, 1}}, {1, {-1, 1}, {1, 2}}, {1, {-1, 1}, {2, 1}}, 
   {1, {1, -1}, {1, 2}}, {1, {1, -1}, {2, 1}}, {1, {1, 1}, {1, 2}}, 
   {1, {1, 1}, {2, 1}}} -> "Full Symmetry"}];

IdentifyWhittenGroup[G_] := Module[{SortedG},

SortedG = Sort[G];

SortedG = ((SortedG /. gamma2Rules) /. gamma1Rules);

If [TrueQ[Head[SortedG] == String], SortedG,

	Clear[Sigma631];
	Clear[Sigma834];
	Clear[Sigma833];
	Clear[Sigma731];
	Clear[Identity3];
	Clear[PI3];
	Clear[PM3];

	Sigma631=Sort[{{1,{1,1,1},{1,2,3}},{1,{1,1,1},{2,1,3}},{1,{1,1,1},{3,2,1}},
{1,{1,1,1},{1,3,2}},{1,{1,1,1},{2,3,1}},{1,{1,1,1},{3,1,2}},
{1,{-1,-1,-1},{1,2,3}},{1,{-1,-1,-1},{2,1,3}},{1,{-1,-1,-1},{3,2,1}},
{1,{-1,-1,-1},{1,3,2}},{1,{-1,-1,-1},{2,3,1}},{1,{-1,-1,-1},{3,1,2}}}];

	Sigma834 = Sort[{{1,{1,1,1},{1,2,3}} , {-1,{1,1,1},{2,1,3}}, {1,{-1,-1,-1},{1,2,3}} , {-1,{-1,-1,-1},{2,1,3}}}];
	Sigma833 = Sort[{{1,{1,1,1},{1,2,3}} , {1,{1,1,1},{2,1,3}} , {1,{1,1,1},{3,2,1}} ,   {1,{1,1,1},{1,3,2}} , {1,{1,1,1},{2,3,1}} , {1,{1,1,1},{3,1,2}} ,   {1,{-1,-1,-1},{1,2,3}} , {1,{-1,-1,-1},{2,1,3}},{1,{-1,-1,-1},{3,2,1}},   {1,{-1,-1,-1},{1,3,2}},{1,{-1,-1,-1},{2,3,1}},{1,{-1,-1,-1},{3,1,2}}}];
	Sigma731 = Sort[{{1,{1,1,1},{1,2,3}}, {1,{1,-1,-1},{1,2,3}}, {1,{1,1,1},{2,3,1}},  {1,{-1,-1,-1},{2,1,3}}  ,{1,{1,1,1},{2,1,3}}, {1,{1,1,-1},{1,2,3}} , {1,{-1,-1,-1},{1,3,2}}, {1,{-1,-1,-1},{3,1,2}},  {1,{-1,-1,1},{1,3,2}} , {1,{-1,1,-1},{3,1,2}}, {1,{1,-1,1},{1,2,3}}, {1,{-1,1,-1},{3 2 1}} ,  {1,{-1,1,1},{1,3,2}}, {1,{1,1,1},{3 2 1}}, {1,{-1,-1,1},{2,3,1}}, {1,{1,-1,-1},{2,3,1}},  {1,{-1,1,-1},{1,2,3}}, {1,{1,-1,-1},{3 2 1}} , {1,{-1,1,1},{1,2,3}}, {1,{1,-1,1},{3 2 1}},  {1,{-1,-1,-1},{2,3,1}}, {1,{-1,1,-1},{1,3,2}}, {1,{1,1,-1},{2,1,3}}, {1,{1,-1,1},{3,1,2}} ,  {1,{-1,-1,-1},{1,2,3}}, {1,{1,1,-1},{1,3,2}}, {1,{-1,1,1},{3,1,2}} , {1,{-1,-1,1},{1,2,3}},  {1,{1,1,-1},{3,1,2}}, {1,{1,1,1},{3,1,2}} , {1,{-1,1,1},{3 2 1}}, {1,{1,1,-1},{2,3,1}},  {1,{1,-1,-1},{2,1,3}}, {1,{1,-1,1},{1,3,2}}, {1,{1,-1,-1},{1,3,2}}, {1,{1,-1,1},{2,3,1}} ,  {1,{-1,-1,-1},{3 2 1}}, {1,{1,-1,-1},{3,1,2}}, {1,{-1,-1,1},{3 2 1}}, {1,{1,1,-1},{3 2 1}},  {1,{-1,1,-1},{2,1,3}}, {1,{-1,1,-1},{2,3,1}} , {1,{-1,1,1},{2,3,1}}, {1,{-1,-1,1},{2,1,3}},  {1,{-1,-1,1},{3,1,2}}, {1,{-1,1,1},{2,1,3}}, {1,{1,1,1},{1,3,2}}, {1,{1,-1,1},{2,1,3}}  }];
	Identity3 = {{1,{1,1,1},{1,2,3}}};
	PI3 = Sort[{{1,{1,1,1},{1,2,3}},{1,{-1,-1,-1},{1,2,3}}}];
	PM3 = Sort[{{1,{1,1,1},{1,2,3}},{-1,{1,1,1},{1,2,3}}}];

	Which[
	SortedG == Sigma631,
	"Sigma631",
	SortedG == Sigma731,
	"Sigma731",
	SortedG == Sigma833,
	"Sigma833",
	SortedG == Sigma834,
	"Sigma834",
	SortedG == Identity3,
	"No Symmetry",
	SortedG == PI3,
	"<e,PI>",
	SortedG == PM3,
	"<e,PM>",
	True,
	ToString[Length[G]]<>" elt subgroup " ]
]
];



IdentifyWhittenElement[gamma_] := 
Which[
    gamma == {1,{1},{1}},
	"e",
	gamma == {1,{-1},{1}},
	"PI",
	gamma == {-1,{1},{1}},
	"PM",
	gamma == {-1,{-1},{1}},
	"MI",
	gamma == {1,{1,1},{1,2}},
	"e",
	gamma == {1,{-1,1},{1,2}},
	"I1",
	gamma == {1,{1,-1},{1,2}},
	"I2",
	gamma == {1,{-1,-1},{1,2}},
	"PI",
	gamma == {1,{1,1},{2,1}},
	"PE",
	gamma == {-1,{1,1},{1,2}},
	"PM",
	gamma == {1,{-1,1},{2,1}},
	"I1E",
	gamma == {1,{1,-1},{2,1}},
	"I2E",
	gamma == {1,{-1,-1},{2,1}},
	"IE",
	gamma == {1,{1,1,1},{1,2,3}},
	"e",
	True,
	ToString[gamma]
];
	


CompileSigmaGroups[LinkList_] := Module[{SigmaGroups,ThisGroup},

SigmaGroups = {};


Do[
	Print["Computing Symmetry Group for ",LinkList[[i]]];
	ThisGroup = SigmaGroup[PD[LinkList[[i]]]];
	SigmaGroups  = Append[SigmaGroups,{Evaluate[LinkList[[i]]],ThisGroup}];
	Print[ThisGroup],
{i,1,Length[LinkList]}];

SigmaGroups

];



(* ::Section:: *)
(*Converting from Hoste-Thistlethwaite Link Tables to Rolfsen's Tables (and back).*)


(* ::Text:: *)
(*We now give a list of transformation rules (thanks to Jacob Rooney and Ellie Dannenberg!) which translate from KnotTheory's numbering scheme for links to the standard Rolfsen notation.*)


RolfsenNumber[L_] := Module[{Name},
Name = L /. 
{Link[2,Alternating,1]->"2^2_1",
Link[4,Alternating,1]->"4^2_1", 
Link[5,Alternating,1]->"5^2_1", 
Link[6,Alternating,1]->"6^2_3", 
Link[6,Alternating,2]->"6^2_2", 
Link[6,Alternating,3]->"6^2_1", 
Link[6,Alternating,4]->"6^3_2", 
Link[6,Alternating,5]->"6^3_1", 
Link[7,Alternating,1]->"7^2_6", 
Link[7,Alternating,2]->"7^2_5", 
Link[7,Alternating,3]->"7^2_4", 
Link[7,Alternating,4]->"7^2_3", 
Link[7,Alternating,5]->"7^2_2", 
Link[7,Alternating,6]->"7^2_1", 
Link[7,Alternating,7]->"7^3_1", 
Link[8,Alternating,1]->"8^2_13", 
Link[8,Alternating,2]->"8^2_10", 
Link[8,Alternating,3]->"8^2_9", 
Link[8,Alternating,4]->"8^2_12", 
Link[8,Alternating,5]->"8^2_11", 
Link[8,Alternating,6]->"8^2_6", 
Link[8,Alternating,7]->"8^2_14", 
Link[8,Alternating,8]->"8^2_7", 
Link[8,Alternating,9]->"8^2_8", 
Link[8,Alternating,10]->"8^2_5", 
Link[8,Alternating,11]->"8^2_3", 
Link[8,Alternating,12]->"8^2_2", 
Link[8,Alternating,13]->"8^2_4", 
Link[8,Alternating,14]->"8^2_1", 
Link[8,Alternating,15]->"8^3_3", 
Link[8,Alternating,16]->"8^3_5", 
Link[8,Alternating,17]->"8^3_2", 
Link[8,Alternating,18]->"8^3_1", 
Link[8,Alternating,19]->"8^3_6", 
Link[8,Alternating,20]->"8^3_4", 
Link[8,Alternating,21]->"8^4_1", 
Link[9,Alternating,1]->"9^2_32", 
Link[9,Alternating,2]->"9^2_31", 
Link[9,Alternating,3]->"9^2_33", 
Link[9,Alternating,4]->"9^2_18", 
Link[9,Alternating,5]->"9^2_30", 
Link[9,Alternating,6]->"9^2_29", 
Link[9,Alternating,7]->"9^2_17", 
Link[9,Alternating,8]->"9^2_25", 
Link[9,Alternating,9]->"9^2_37", 
Link[9,Alternating,10]->"9^2_36", 
Link[9,Alternating,11]->"9^2_26", 
Link[9,Alternating,12]->"9^2_14", 
Link[9,Alternating,13]->"9^2_16", 
Link[9,Alternating,14]->"9^2_13", 
Link[9,Alternating,15]->"9^2_15", 
Link[9,Alternating,16]->"9^2_28", 
Link[9,Alternating,17]->"9^2_27", 
Link[9,Alternating,18]->"9^2_10", 
Link[9,Alternating,19]->"9^2_38", 
Link[9,Alternating,20]->"9^2_42", 
Link[9,Alternating,21]->"9^2_34", 
Link[9,Alternating,22]->"9^2_35", 
Link[9,Alternating,23]->"9^2_22", 
Link[9,Alternating,24]->"9^2_21", 
Link[9,Alternating,25]->"9^2_8", 
Link[9,Alternating,26]->"9^2_11", 
Link[9,Alternating,27]->"9^2_12", 
Link[9,Alternating,28]->"9^2_20", 
Link[9,Alternating,29]->"9^2_19", 
Link[9,Alternating,30]->"9^2_3", 
Link[9,Alternating,31]->"9^2_39", 
Link[9,Alternating,32]->"9^2_40", 
Link[9,Alternating,33]->"9^2_24", 
Link[9,Alternating,34]->"9^2_6", 
Link[9,Alternating,35]->"9^2_9", 
Link[9,Alternating,36]->"9^2_1", 
Link[9,Alternating,37]->"9^2_7", 
Link[9,Alternating,38]->"9^2_5", 
Link[9,Alternating,39]->"9^2_2", 
Link[9,Alternating,40]->"9^2_4", 
Link[9,Alternating,41]->"9^2_23", 
Link[9,Alternating,42]->"9^2_41", 
Link[9,Alternating,43]->"9^3_4", 
Link[9,Alternating,44]->"9^3_3", 
Link[9,Alternating,45]->"9^3_7", 
Link[9,Alternating,46]->"9^3_10", 
Link[9,Alternating,47]->"9^3_2", 
Link[9,Alternating,48]->"9^3_5", 
Link[9,Alternating,49]->"9^3_6", 
Link[9,Alternating,50]->"9^3_1", 
Link[9,Alternating,51]->"9^3_11", 
Link[9,Alternating,52]->"9^3_8", 
Link[9,Alternating,53]->"9^3_12", 
Link[9,Alternating,54]->"9^3_9", 
Link[9,Alternating,55]->"9^4_1",
Link[6,NonAlternating,1] -> "6^3_3",
Link[7,NonAlternating,1] -> "7^2_7",
Link[7,NonAlternating,2] -> "7^2_8",
Link[8,NonAlternating,1] -> "8^2_16",
Link[8,NonAlternating,2] -> "8^2_15",
Link[8,NonAlternating,3] -> "8^3_7",
Link[8,NonAlternating,4] -> "8^3_8",
Link[8,NonAlternating,5] -> "8^3_9",
Link[8,NonAlternating,6] -> "8^3_10",
Link[8,NonAlternating,7] -> "8^4_2",
Link[8,NonAlternating,8] -> "8^4_3",
Link[9,NonAlternating,28]->"9^3_20",
Link[9,NonAlternating,22]->"9^3_15",
Link[9,NonAlternating,9]->"9^2_60",
Link[9,NonAlternating,12]->"9^2_59",
Link[9,NonAlternating,11]->"9^2_57",
Link[9,NonAlternating,8]->"9^2_56",
Link[9,NonAlternating,6]->"9^2_55",
Link[9,NonAlternating,7]->"9^2_48",
Link[9,NonAlternating,3]->"9^2_47",
Link[9,NonAlternating,2]->"9^2_46",
Link[9,NonAlternating,5]->"9^2_44",
Link[9,NonAlternating,16]->"9^2_51",
Link[9,NonAlternating,17]->"9^2_52",
Link[9,NonAlternating,20]->"9^3_16",
Link[9,NonAlternating,24]->"9^3_14",
Link[9,NonAlternating,10]->"9^2_58",
Link[9,NonAlternating,13]->"9^2_54",
Link[9,NonAlternating,19]->"9^2_61",
Link[9,NonAlternating,18]->"9^2_53",
Link[9,NonAlternating,15]->"9^2_49",
Link[9,NonAlternating,4]->"9^2_43",
Link[9,NonAlternating,23]->"9^3_13",
Link[9,NonAlternating,21]->"9^3_17",
Link[9,NonAlternating,25]->"9^3_18",
Link[9,NonAlternating,26]->"9^3_19",
Link[9,NonAlternating,27]->"9^3_21",
Link[9,NonAlternating,1]->"9^2_45",
Link[9,NonAlternating,14]->"9^2_50"
};

If [Name[[0]] == String, Return[Name]];
ToString[Name[[1]]]<>"^"<>ToString[NumComponents[Name]]<>"_?"

];

ThistlethwaiteNumber[L_] := Module[{Name,Cr,Comp,Ind},
Name = L /. {
"2^2_1"->Link[2,Alternating,1],
"4^2_1"->Link[4,Alternating,1],
"5^2_1"->Link[5,Alternating,1],
"6^2_3"->Link[6,Alternating,1],
"6^2_2"->Link[6,Alternating,2],
"6^2_1"->Link[6,Alternating,3],
"6^3_2"->Link[6,Alternating,4],
"6^3_1"->Link[6,Alternating,5],
"7^2_6"->Link[7,Alternating,1],
"7^2_5"->Link[7,Alternating,2],
"7^2_4"->Link[7,Alternating,3],
"7^2_3"->Link[7,Alternating,4],
"7^2_2"->Link[7,Alternating,5],
"7^2_1"->Link[7,Alternating,6],
"7^3_1"->Link[7,Alternating,7],
"8^2_13"->Link[8,Alternating,1],
"8^2_10"->Link[8,Alternating,2],
"8^2_9"->Link[8,Alternating,3],
"8^2_12"->Link[8,Alternating,4],
"8^2_11"->Link[8,Alternating,5],
"8^2_6"->Link[8,Alternating,6],
"8^2_14"->Link[8,Alternating,7],
"8^2_7"->Link[8,Alternating,8],
"8^2_8"->Link[8,Alternating,9],
"8^2_5"->Link[8,Alternating,10],
"8^2_3"->Link[8,Alternating,11],
"8^2_2"->Link[8,Alternating,12],
"8^2_4"->Link[8,Alternating,13],
"8^2_1"->Link[8,Alternating,14],
"8^3_3"->Link[8,Alternating,15],
"8^3_5"->Link[8,Alternating,16],
"8^3_2"->Link[8,Alternating,17],
"8^3_1"->Link[8,Alternating,18],
"8^3_6"->Link[8,Alternating,19],
"8^3_4"->Link[8,Alternating,20],
"8^4_1"->Link[8,Alternating,21],
"9^2_32"->Link[9,Alternating,1],
"9^2_31"->Link[9,Alternating,2],
"9^2_33"->Link[9,Alternating,3],
"9^2_18"->Link[9,Alternating,4],
"9^2_30"->Link[9,Alternating,5],
"9^2_29"->Link[9,Alternating,6],
"9^2_17"->Link[9,Alternating,7],
"9^2_25"->Link[9,Alternating,8],
"9^2_37"->Link[9,Alternating,9],
"9^2_36"->Link[9,Alternating,10],
"9^2_26"->Link[9,Alternating,11],
"9^2_14"->Link[9,Alternating,12],
"9^2_16"->Link[9,Alternating,13],
"9^2_13"->Link[9,Alternating,14],
"9^2_15"->Link[9,Alternating,15],
"9^2_28"->Link[9,Alternating,16],
"9^2_27"->Link[9,Alternating,17],
"9^2_10"->Link[9,Alternating,18],
"9^2_38"->Link[9,Alternating,19],
"9^2_42"->Link[9,Alternating,20],
"9^2_34"->Link[9,Alternating,21],
"9^2_35"->Link[9,Alternating,22],
"9^2_22"->Link[9,Alternating,23],
"9^2_21"->Link[9,Alternating,24],
"9^2_8"->Link[9,Alternating,25],
"9^2_11"->Link[9,Alternating,26],
"9^2_12"->Link[9,Alternating,27],
"9^2_20"->Link[9,Alternating,28],
"9^2_19"->Link[9,Alternating,29],
"9^2_3"->Link[9,Alternating,30],
"9^2_39"->Link[9,Alternating,31],
"9^2_40"->Link[9,Alternating,32],
"9^2_24"->Link[9,Alternating,33],
"9^2_6"->Link[9,Alternating,34],
"9^2_9"->Link[9,Alternating,35],
"9^2_1"->Link[9,Alternating,36],
"9^2_7"->Link[9,Alternating,37],
"9^2_5"->Link[9,Alternating,38],
"9^2_2"->Link[9,Alternating,39],
"9^2_4"->Link[9,Alternating,40],
"9^2_23"->Link[9,Alternating,41],
"9^2_41"->Link[9,Alternating,42],
"9^3_4"->Link[9,Alternating,43],
"9^3_3"->Link[9,Alternating,44],
"9^3_7"->Link[9,Alternating,45],
"9^3_10"->Link[9,Alternating,46],
"9^3_2"->Link[9,Alternating,47],
"9^3_5"->Link[9,Alternating,48],
"9^3_6"->Link[9,Alternating,49],
"9^3_1"->Link[9,Alternating,50],
"9^3_11"->Link[9,Alternating,51],
"9^3_8"->Link[9,Alternating,52],
"9^3_12"->Link[9,Alternating,53],
"9^3_9"->Link[9,Alternating,54],
"9^4_1"->Link[9,Alternating,55],
"6^3_3"->Link[6,NonAlternating,1],
"7^2_7"->Link[7,NonAlternating,1],
"7^2_8"->Link[7,NonAlternating,2],
"8^2_16"->Link[8,NonAlternating,1],
"8^2_15"->Link[8,NonAlternating,2],
"8^3_7"->Link[8,NonAlternating,3],
"8^3_8"->Link[8,NonAlternating,4],
"8^3_9"->Link[8,NonAlternating,5],
"8^3_10"->Link[8,NonAlternating,6],
"8^4_2"->Link[8,NonAlternating,7],
"8^4_3"->Link[8,NonAlternating,8],
"9^3_20"->Link[9,NonAlternating,28],
"9^3_15"->Link[9,NonAlternating,22],
"9^2_60"->Link[9,NonAlternating,9],
"9^2_59"->Link[9,NonAlternating,12],
"9^2_57"->Link[9,NonAlternating,11],
"9^2_56"->Link[9,NonAlternating,8],
"9^2_55"->Link[9,NonAlternating,6],
"9^2_48"->Link[9,NonAlternating,7],
"9^2_47"->Link[9,NonAlternating,3],
"9^2_46"->Link[9,NonAlternating,2],
"9^2_44"->Link[9,NonAlternating,5],
"9^2_51"->Link[9,NonAlternating,16],
"9^2_52"->Link[9,NonAlternating,17],
"9^3_16"->Link[9,NonAlternating,20],
"9^3_14"->Link[9,NonAlternating,24],
"9^2_58"->Link[9,NonAlternating,10],
"9^2_54"->Link[9,NonAlternating,13],
"9^2_61"->Link[9,NonAlternating,19],
"9^2_53"->Link[9,NonAlternating,18],
"9^2_49"->Link[9,NonAlternating,15],
"9^2_43"->Link[9,NonAlternating,4],
"9^3_13"->Link[9,NonAlternating,23],
"9^3_17"->Link[9,NonAlternating,21],
"9^3_18"->Link[9,NonAlternating,25],
"9^3_19"->Link[9,NonAlternating,26],
"9^3_21"->Link[9,NonAlternating,27],
"9^2_45"->Link[9,NonAlternating,1],
"9^2_50"->Link[9,NonAlternating,14]};

If [Name[[0]] == Link, Return[Name]];
{Cr,Comp,Ind} = StringSplit[Name,{"^","_"}];

"Link["<>ToString[Cr]<>",?,?]"

];



(* ::Section:: *)
(*Display of Whitten Groups *)


SortRolfsen[A_,B_] := Module[{CrA,CompA,IndA,CrB,CompB,IndB},

	{CrA,CompA,IndA} = ToExpression[StringSplit[A[[1]],"^" | "_"]];
	{CrB,CompB,IndB} = ToExpression[StringSplit[B[[1]],"^" | "_"]];
	
	If[CrA < CrB, True, (* Dictionary Sort on Cr, Comp, Ind pairs *)
	 If[CrA == CrB,
	   If[CompA < CompB, True,
         If[CompA == CompB,
            If[IndA < IndB,True,False],
         False],
        False],
      False],
    False]
];


WhittenGroupTable[T_]:= Module[{},
Print["hi"];
Sort[
	Table[{RolfsenNumber[T[[i,1]]], 
		   IdentifyWhittenGroup[T[[i,2]]], 
		   StringJoin[Riffle[WhittenPrettyPrint /@ Sort[T[[i,2]],WhittenEltCompare],", "]]},
    {i,1,Length[T]}],
SortRolfsen]
]



WhittenGroupGraphic[ThisLink_Link,G_] := WhittenGroupGraphic[PD[ThisLink],G];
WhittenGroupGraphic[ThisPD_PD,G_] := {NiceDrawPD[ThisPD,{EdgeLabels->False}],IdentifyWhittenGroup[G],StringJoin[Riffle[WhittenPrettyPrint /@ Sort[G,WhittenEltCompare],", "]]};
WhittenGroupGraphic[Pic_Graphics|Pic_Rotate,G_] := {Pic,IdentifyWhittenGroup[G],StringJoin[Riffle[WhittenPrettyPrint /@ Sort[G,WhittenEltCompare],", "]]};


WhittenGroupTableWithGraphics[T_] := Column[Join[WhittenGroupGraphic[#[[1]],#[[2]]],{RolfsenNumber[#[[1]]]}]] & /@ T;


(* ::Section:: *)
(*Digrammatic Symmetries*)


(* ::Text:: *)
(*This section contains some code for detecting symmetries in PD codes.*)


VertexGraph[ThisPD_PD] := Module[{verts,edges},

verts = Union[List @@ ThisPD,AllArcs[ThisPD]];
edges =Flatten[{# <->#[[1]], # <-> #[[2]], # <-> #[[3]], # <-> #[[4]]} & /@ (List @@ ThisPD),1];
Graph[verts,edges]

];

FaceGraph[ThisPD_PD] := Module[{edges,faces,fgEdges,i},

faces = PDFaces[ThisPD];
edges = AllArcs[ThisPD];
fgEdges = Table[Select[faces,MemberQ[#,edges[[i]]]&],{i,1,Length[edges]}];
Graph[faces, #[[1]] <-> #[[2]] & /@ fgEdges]

];

OrderFace[F_List] := RotateLeft[F,Position[F,Min[F]][[1,1]]-1];

FaceMap[ThisPD_PD,EdgeMap_List] := Module[{Faces,NewFaces,i,P,R},
	
	Faces = OrderFace /@ PDFaces[ThisPD];
	NewFaces = OrderFace /@ (Faces /. EdgeMap);
	
	(*Print["Original faces (ordered):",Faces];*)
	(*Print["New faces (ordered):",NewFaces];*)

	P =Quiet[PermutationList[FindPermutation[NewFaces,Faces],Length[Faces]]] /. _PermutationList -> {};
	R = Quiet[
		PermutationList[FindPermutation[
		OrderFace/@ (Reverse /@ NewFaces),Faces],Length[Faces]
		]] /. _PermutationList -> {};
	Which[
	Length[P] != 0,
	(*Print["Rotational Permutation (O->N):",P];*){"F",P},
	Length[R] != 0,
	(*Print["Reflection Permutation (O->N):",R];*){"R",R},
	True,
	(*Print["No map found"];*){}
  ]

];

RotateSort[X_] := RotateLeft[X,Position[X,Min[List @@ X]][[1]]-1];

CrossingMap[ThisPD_PD,EdgeMap_List] :=Module[{F,R,SignedP,Signs},
F = Quiet[PermutationList[FindPermutation[RotateSort /@ (ThisPD /. EdgeMap),RotateSort /@ ThisPD],Length[ThisPD]]] /. {{} -> Range[Length[ThisPD]],_PermutationList -> {}};
R = Quiet[PermutationList[FindPermutation[RotateSort /@ Reverse /@ (ThisPD /. EdgeMap),RotateSort /@ ThisPD],Length[ThisPD]]] /. {{} -> Range[Length[ThisPD]],_PermutationList -> {}};
Which[F != {},
SignedP = Table[If[MemberQ[{1,3},Position[ThisPD[[i]],ThisPD[[F[[i]],1]] /. EdgeMap][[1,1]]],+F[[i]],-F[[i]]],{i,1,Length[F]}];
(* Now check whether all the signs are consistent *)
Signs = Union[Sign /@ SignedP];
If [Length[Signs] == 1,{"F",SignedP},{}],
R != {},
SignedP = Table[If[MemberQ[{1,3},Position[ThisPD[[i]],ThisPD[[R[[i]],1]] /. EdgeMap][[1,1]]],-R[[i]],+R[[i]]],{i,1,Length[R]}];
Signs = Union[Sign /@ SignedP];
If [Length[Signs] == 1,{"R",SignedP},{}],
True,
{}
]
];

PossibleComponentPermutations[ThisPD_] := Module[{Components,Clist,TPerms,StdOrder,PermRules,i,j},

Components = AllComponents[ThisPD];
(*Print[Components];*)
Clist = Gather[Range[Length[Components]],Length[Components[[#1]]] == Length[Components[[#2]]] &];
StdOrder = Join @@ Clist;
TPerms = Apply[Join,Tuples[Permutations /@ Clist],{1}];
(*Print[StdOrder];*)
(*Print[TPerms];*)
PermRules = Table[StdOrder[[i]] -> TPerms[[j,i]],{j,1,Length[TPerms]},{i,1,Length[StdOrder]}];
Table[Range[Length[StdOrder]] /. PermRules[[j]],{j,1,Length[TPerms]}]

];

ComponentwisePerms[{Op_,A_,B_},Components_] := Module[{i,j,Bop,Code},
Bop = Op[Components[[B]]];
Code = If[TrueQ[Op == Identity],"F","R"];
Table[Table[Components[[A,j]] -> RotateRight[Bop,i][[j]],{j,1,Length[Bop]}],{i,0,Length[Bop]-1}]
];

EdgeMaps[ThisPD_] := Module[{CPerms,Signs,SignPerms,Components,Needed,PermRules},

Components = AllComponents[ThisPD];
Signs = Tuples[{Identity,Reverse},NumComponents[ThisPD]];
SignPerms = Map[Flatten,Map[MapIndexed[List,#] &,Transpose /@ Tuples[{Signs,PossibleComponentPermutations[ThisPD]}]],{2}];
Needed = Union[Flatten[SignPerms,1]];
PermRules = # -> ComponentwisePerms[#,Components] & /@ Needed;
(*Print[PermRules];*)
(*Print[SignPerms];*)
SignPerms = Tuples[{{Transpose[#[[1]]]},Apply[Join,Tuples[#[[2]]],1]}] & /@ (Map[Transpose,Map[{{#[[1]],#[[2]]},#} & ,SignPerms,{2}],{1}] /. PermRules);
(*Print[SignPerms];*)
Flatten[SignPerms,1]
];

RawDiagrammaticSymmetries[ThisPD_PD] := Select[{#[[2]],#[[1]], CrossingMap[ThisPD,#[[2]]],FaceMap[ThisPD,#[[2]]]} & /@Select[Select[EdgeMaps[ThisPD],CrossingMap[ThisPD,#[[2]]] != {} &],FaceMap[ThisPD,#[[2]]] != {} &],#[[3,1]] == #[[4,1]] &];

WhittenForm[Sym_] := {Sign[Sym[[3,2,1]]],Sym[[2,1]] /. {Identity->+1,Reverse->-1},Sym[[2,2]]}

DiagrammaticSymmetries[ThisPD_PD] := Module[{RawSyms},
 RawSyms = RawDiagrammaticSymmetries[ThisPD];
{#[[1,1]],#[[2]]} & /@ Map[Union,Map[Transpose,Gather[{WhittenForm[#],#[[1]]} & /@ RawSyms,#1[[1]] == #2[[1]] &]],{2}]
];

DiagrammaticSymmetryGroup[ThisPD_PD] := IdentifyWhittenGroup[Transpose[DiagrammaticSymmetries[ThisPD]][[1]]];



(* ::Section:: *)
(*Classify PD*)


(* ::Text:: *)
(*This section is dedicated to a method for classifying PD-codes. The module needs the auxiliary file "invariantlist.txt" to be in the same directory as your working notebook. The module BuildInvariantList can be used to generate this file.*)
(**)
(*For now, the only invariants being tested are the Kauffman Polynomial and Knot Signature. These two invariants classify all oriented knots (prime and composite) through (with the usual caveat of additivity of crossing number under connected sums) except for the list below.*)
(**)
(*{"Knot[8, 17]", "Knot[8, 17]r"}*)
(*{"Knot[9, 32]", "Knot[9, 32]r"}*)
(*{"Knot[9, 32]m", "Knot[9, 32]mr"}*)
(*{"Knot[9, 33]", "Knot[9, 33]r"}*)
(*{"Knot[9, 33]m", "Knot[9, 33]mr"}*)
(**)
(*Note that the argument is a KnotTheory-style PD-code and the return value is a string.*)


KnotsWithInvariants=ReadList[NotebookDirectory[]<>"invariantlist_old.txt"];

ClassifyKnot[PD_]:=Module[{Invariants,PossibleKnots},
Invariants={Expand[Kauffman[PD][q,y]],KnotSignature[PD]};
PossibleKnots=Select[KnotsWithInvariants,Invariants==#[[3]]&];
(*Print[Invariants];
Print[PossibleKnots];*)
Table[PossibleKnots[[i]][[1]],{i,1,Length[PossibleKnots]}]
]

(*Ignore the following code for now!*)
ClassifyKnotNew[PD_,InvariantList_,InvariantNames_]:=Module[{Invariants,PossibleKnots,CheckKnot,GetInvariants},
Invariants=Table[InvariantList[[i]][PD],{i,1,Length[InvariantList]}];
(*Print[Invariants]*)
CheckKnot[CurrentInvariants_]:=If[Invariants==CurrentInvariants,True,False];
GetInvariants[CurrentKnot_]:=InvariantNames/.CurrentKnot[[3]];
PossibleKnots=Select[KnotsWithInvariants,CheckKnot[GetInvariants[#]]&];
Print[PossibleKnots]
]


(* ::Text:: *)
(*The following code will build the auxiliary file "invariantslist.txt".*)


BuildInvariantList[KnotList_,KnotNames_,InvariantList_,InvariantNames_]:=Module[{ComputedInvariants,KnotsWithInvariants},
ComputedInvariants=Table[Table[InvariantNames[[i]] -> InvariantList[[i]][KnotList[[j]]],{i,1,Length[InvariantList]}],{j,1,Length[KnotList]}];
(*Print[ComputedInvariants]*)
KnotsWithInvariants=Table[{KnotNames[[i]],KnotList[[i]],ComputedInvariants[[i]]},{i,1,Length[KnotList]}];
(*Print[KnotsWithInvariants]*)
Export["invariantlist.txt",KnotsWithInvariants];
]


EndPackage[]
