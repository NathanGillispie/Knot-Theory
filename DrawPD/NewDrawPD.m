(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Nathan's edits*)
BeginPackage["KnotTheory`"];

PD; X; OuterFace; Gap; Colour; StrandColour

NewDrawPD::loading="Requires KnotTheory"

NewDrawPD::usage = "
  NewDrawPD[pd_] is a faster and more reliable DrawPD
"

NewDrawPD::about = "
  These are some words about NewDrawPD
"

ShowCirclePacking::usage="ShowCirclePacking[pd_]"
ShowCirclePacking::about=""

ShowTriangulationGraph::usage="ShowTriangulationGraph[pd_]"
ShowTriangulationGraph::about=""

Begin["`DrawPD`"]

ShowCirclePacking[L_]:= ShowCirclePacking[PD[L]]
ShowCirclePacking[pd_PD]:= (
t=AddPositionsBound[DefaultDirichlet[Triangulate[pd]]];
    t=PutInside[t,DefaultOuterFace[t]];
t=ApplyFLMap[t,Moebius[Nest[BalanceStep[t,#]&,0,100]]];
Graphics[Join[Table[{Circle[xyCoords[t[[i,centre]]],Abs@t[[i,r]]]},{i,Length@t}],
Table[Text[i,
          xyCoords[t[[i,centre]]]],{i,Length[t]}]]]
)


ShowCirclePacking[pd_PD,options_]:= Module[{optionsList= Map[Apply[List,#]&,options],newt},
t=AddPositionsBound[DefaultDirichlet[Triangulate[pd]]];
    t=PutInside[t,DefaultOuterFace[t]];
t=ApplyFLMap[t,Moebius[Nest[BalanceStep[t,#]&,0,100]]];
t= AddGraphicsObjs[t,{DefaultGap[t]}];
Switch[Length@optionsList,
3,
If[optionsList[[1,1]]=="X"&&optionsList[[2,1]]=="e"&&optionsList[[3,1]]=="f",
Graphics[Join[Table[{Switch[t[[i,type]],optionsList[[1,1]],optionsList[[1,2]],optionsList[[2,1]],optionsList[[2,2]],optionsList[[3,1]],optionsList[[3,2]], _, Black],Circle[xyCoords[t[[i,centre]]],Abs@t[[i,r]]]},{i,Length@t}],{Black},
Table[Text[i,
          xyCoords[t[[i,centre]]]],{i,Length[t]}]]]
],
4,
If[optionsList[[1,1]]=="X"&&optionsList[[2,1]]=="e"&&optionsList[[3,1]]=="f"&&optionsList[[4,1]]=="ShowDiagram"&&optionsList[[4,2]]==True,
newt=Flatten[FieldValues[t,graphicsObjs]];
(*for composite knots*)
Do[While[newt[[i,3,2]]-newt[[i,3,1]]>Pi,newt[[i,3,2]]-=2Pi];
While[newt[[i,3,1]]-newt[[i,3,2]]>Pi,newt[[i,3,2]]-=2Pi];
,{i,Length@newt}];
Graphics[Join[newt,
Table[{Switch[t[[i,type]],optionsList[[1,1]],optionsList[[1,2]],optionsList[[2,1]],optionsList[[2,2]],optionsList[[3,1]],optionsList[[3,2]], _, Black],Circle[xyCoords[t[[i,centre]]],Abs@t[[i,r]]]},{i,Length@t}],{Black},
Table[Text[i,
          xyCoords[t[[i,centre]]]],{i,Length[t]}]]]
]
]
]

(*shows the graph too*)
ShowTriangulationGraph[pd_PD]:= Module[{t,defaultGap,g,defaultOuterFace},
t=AddPositionsBound[DefaultDirichlet[Triangulate[pd]]];
defaultOuterFace=DefaultOuterFace[t];
    t=PutInside[t,defaultOuterFace];
t=ApplyFLMap[t,Moebius[Nest[BalanceStep[t,#]&,0,100]]];
defaultGap=DefaultGap[t]/2.0;
g=Join[Table[If[defaultOuterFace!=i,
Disk[xyCoords[t[[i,centre]]],defaultGap]],{i,Length@t}],
Table[If [(t[[i,neighbours,j]]!= defaultOuterFace)&&(i!=defaultOuterFace),
Line[{xyCoords[t[[i,centre]]],xyCoords[t[[t[[i,neighbours,j]],centre]]]}]]
,{i,Length@t},{j,Length[t[[i,neighbours]]]}]];
t=AddGraphicsObjs[t,{DefaultGap[t]}];
t=Flatten[FieldValues[t,graphicsObjs]];
(*for composite knots*)
Do[While[t[[i,3,2]]-t[[i,3,1]]>Pi,t[[i,3,2]]-=2Pi];
While[t[[i,3,1]]-t[[i,3,2]]>Pi,t[[i,3,2]]-=2Pi];
,{i,Length@t}];
Graphics[Join[{Red,g},{Black,t}],AspectRatio->1]
]

NewDrawPD[L_] := NewDrawPD[PD[L]]
NewDrawPD[L_,options_] := NewDrawPD[PD[L],options]
NewDrawPD[pd_PD]:=( 
    CreditMessage["DrawPD was written by Emily Redelmeier at the University of Toronto in the summers of 2003 and 2004."];
    t=AddPositionsBound[DefaultDirichlet[Triangulate[pd]]];
    t=PutInside[t,DefaultOuterFace[t]];
(*same as drawPD but replaces the Balance function*)
t=ApplyFLMap[t,Moebius[Nest[BalanceStep[t,#]&,0,100]]];
    t=AddGraphicsObjs[t,{DefaultGap[t]}];
(*t=ColourStrands[t,{}];*)
t=Flatten[FieldValues[t,graphicsObjs]];
(*for composite knots*)
(*Do[While[t\[LeftDoubleBracket]i,3,2\[RightDoubleBracket]-t\[LeftDoubleBracket]i,3,1\[RightDoubleBracket]>Pi,t\[LeftDoubleBracket]i,3,2\[RightDoubleBracket]-=2Pi];
While[t\[LeftDoubleBracket]i,3,1\[RightDoubleBracket]-t\[LeftDoubleBracket]i,3,2\[RightDoubleBracket]>Pi,t\[LeftDoubleBracket]i,3,2\[RightDoubleBracket]-=2Pi];
,{i,Length@t}];*)
 Graphics[t,AspectRatio->1]
)

NewDrawPD[pd_PD,options_]:=(optionsList=Map[Apply[List,#]&,options];
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
   t=ApplyFLMap[t,Moebius[Nest[BalanceStep[t,#]&,0,100]]];
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
]]==StrandColour&][[1,
                      2]]]}]]]]];
t=Flatten[FieldValues[t,graphicsObjs]];
(*for composite knots*)
Do[While[t[[i,3,2]]-t[[i,3,1]]>Pi,t[[i,3,2]]-=2Pi];
While[t[[i,3,1]]-t[[i,3,2]]>Pi,t[[i,3,2]]-=2Pi];
,{i,Length@t}];
 Graphics[t,AspectRatio->1])
 
 Message[NewDrawPD::loading]
 End[]; EndPackage[]
