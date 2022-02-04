(* ::Package:: *)

(* Mathematica Version: 5.0; 5.1; 5.2  *)

(* Name: `LinKnots` *)

(* Title: LinKnots *)

(* Authors: M.Ochiai and N.Imafuji are the authors of Knot 2000 (K2K), 
S. Jablan and R. Sazdanovic of LinKnot (K2KC), and Dror Bar-Natan of \
KnotTables *)

(* Copyright: Copyright 2006 *)

(* History: current Version 1.1 at 24.03.2006 *)

(* Summary: The complete package K2K.m works with two packages: 
knotbycomp.m by M.Ochiai and N.Imafuji, and LinKnot.m  
by S. Jablan and R. Sazdanovic. The package LinKnot.m, 
contained in K2K.m, works with knots and links in Conway 
notation, with basic polyhedra with N <= 20 vertices, 
producing their pdata, Dowker codes and Millet codes, 
and calculating linking numbers and estimated unknotting 
and unlinking numbers. The package KnotTheory by Dror Bar-Naran
computes new KL invariants: Vassiliev and Khovanov polynomials *)
  
(* Context: `LinKnots` *)

(* Keywords: Knot Theory, Knots, Links, Conway Notation, Basic Polyhedra, \
Invariants *)

(* Requirements: DiscreteMath`Combinatorica,
"LinearAlgebra`MatrixManipulation`",
"DiscreteMath`ComputationalGeometry`"*)

Off[General::"newpkg"];
Off[General::"obspkg"];
Off[General::"poly"];
Off[PolynomialGCD::"argm"];
Off[Solve::"svars"];
Off[AppendTo::"rvalue"];
Off[Divide::"infy"];
Off[Syntax::"com"];
Off[General::"compat"];



BeginPackage["LinKnots`", "Combinatorica`","ComputationalGeometry`"
,"PolyBase`", "NalgBase`","KnotLinkBase`","LinksGouldExplorer`","KnotTheory`","PolyhedronOperations`"];


(* BeginPackage["LinKnots`", "DiscreteMath`Combinatorica`",
"DiscreteMath`ComputationalGeometry`",
"LinearAlgebra`MatrixManipulation`","PolyBase`",
"NalgBase`","KnotLinkBase`" 
,"LinksGouldExplorer`","KnotTheory`","PolyhedronOperations`"]; *)


 Get["knotbycomp.m"];
 Get["RecBase.txt"];
 Get["Classic.txt"];
 Get["lllll.txt"]; 
 Get["knotlists.txt"];
 Get["doubleedges.txt"];
 Get["pl.txt"];
 Get["uu.txt"]; 
 Get["uuu.txt"];
 Get["uuuu.txt"]; 
 Get["polyhedra.txt"]; 
 Get["ttt.txt"];  
 Get["morlinks.txt"]; 
 Get["hypvoltable.txt"]; 
 Get["virtknotab.m"];
 Get["virttab8.txt"]; 
 Get["fordetect01.txt"]; 
 Get["fordetect02.txt"];
 Get["tableallknotspoly.txt"];



 
 
(* USAGE MESSAGES *)

fFaces::usage="" 
fKLTab::usage=""  
            
fCreatePData::usage = 
"fCreatePData[Conway_String] calculates pdata for 
given knot or link given in Conway notation."

Dow::usage = 
"Dow[Conway_String] produces Dowker code of a knot or link K.
The first data of result are lengths of components, and the second Dowker 
code with the signs of crossing points."

fGaussExt::usage="fFormDL[Ul_] produces extended Gauss code with signs 
of knot or link given by Conway symbol or Dowker code." 

LinkingNo::usage="LinkingNo[Ul_] computes the linking number of a 
link given by Conway symbol or Dowker code."

fGenSign::usage =
"fGenSign[Conway_String] returns the signlist of crossing points 
of a KL given by Conway symbol."

fDowfromPD::usage = 
"fDowfromPD[pdata_List] computes Dowker code from pdata."

fDowfromPData::usage = 
"Correction for fDowfromPD"


fMillett::usage="fMillett[Conway_String] writes in the file 
izlaz11.txt the Millet code of a knot or link and calculates
polynomial invariants, including Homfly polynomial."

UnKnotLink::usage="UnKnotLink[Ul_] returns an (estimated)
unknotting or unlinking number of a given knot or link given 
by its Conway symbol, Dowker code or pdata. For the 
reduction of knots or links it uses ReductionKnotLink."

RatReduce::usage="RatReduce[Conway_String] reduces 
Conway symbol of a rational KL."

UnR::usage="UnR[Conway_String] returns an estimated unknotting 
or unlinking number of a given rational knot or link. It is
computed by using continued fractions."

fCreateGraphics::usage="fCreateGraphics[Conway_String] writes in 
the file c:\\Program Files\\KnotPlot\\graphics.txt the coordinates 
of a knot or link. The file obtained can be loaded in the program 
KnotPlot by writing load graphics.txt in the KnotPlot command line."
        
GetKnotLink::usage="GetKnotLink[Str_String, No_Integer].
Str is the name of a list of alternating or non-alternating 
knots and links with a specified number of crossings (for example, 
a10 or n10, where a stands for alternating and n for non-alternating 
knots and links) and No is the number of a particullar desired knot ot 
link from that list. As a result, GetKnotLink returns its Conway symbol."

NumberOfKL::usage="NumberOfKL[Str_String] Str is the name of a list
of knots and links with a specified number of crossings 
(for example, a10 or n10, where a stands for alternating and n 
for non-alternating knots and links, written as a string). The 
function NumberOfKL gives the number of alternating or 
non-alternating knots and links with a specified number of 
crossings."

Dowker::usage ="Interior function"
fConvert::usage="Interior function"
fConvertPoly::usage="Interior function"
fBlock::usage="Interior function"
fComma::usage="Interior function"
fPlus::usage="Interior function"
fInv::usage="Interior function"
fSpace::usage="Interior function"
fRef::usage="Interior function"
fGen::usage="Interior function"
fMoja::usage="Interior function"
fDowkerCode::usage="Interior function"

fDToD::usage="fDtoD[Con_String] calculates Dowker code 
of a KL given in Conway notation."

fProjections::usage="fProjections[Con_String] calculates 
Conway symbols of all projections of a KL given by its Conway symbol."

fPr::usage=""

fProject::usage=""

fFixP::usage=""
fChangePart::usage=""

fNewGraph::usage=""

fBal::usage=""

UnKnotLinkNo1::usage=""

parr::usage=""

fDiffProjectionsAltKL::usage="fDiffProjectionsAltKL[Con_String] calculates 
Conway symbols and minimal Dowker codes of all non-isomorphic 
projections of a KL given by its Conway symbol."

MinDowProjAltKL ::usage="MinDowProjAltKL[Ul_] calculates minimal 
Dowker code without signs (in Knotscape form) for a given 
alternating KL projection. An input is its Conway symbol, Dowker code, 
or pdata."

MinDowAltKL::usage="MinDowAltKL[Con_String] calculates 
minimal Dowker code of any alternating KL given by its Conway symbol."

SameAltProjKL::usage="SameAltProjKL[Con1_String,Con2_String] compares two 
alternating KL projections given by Conway symbols, 
Dowker codes, or pdata. The result is 1 
for equal, and 0 for non-equal projections."

SameAltConKL::usage="SameAltConKL[Con1_String,Con2_String] compares 
two alternating KLs given by Conway symbols. 
The result is 1 for equal, and 0 for non-equal KLs."

fFindCon::usage="fFindCon[Ul_] finds Conway symbol of 
any alternating KL with at most 12 vertices 
given by its Dowker code with signs, or by pdata."

fAlexPoly::usage="fAlexPoly[Conway_String] calculates multi-variable 
Alexander polynomial of a KL from its Conway symbol 
by using Writinger presentation (generators of KL)."

fSeifert::usage="fSeifert[Ul_] calculates Seifert matrix of a 
KL given by its Conway symbol, Dowker code, 
or pdata. Its basic function SeifertMatrix  
is written by Stepan Orevkov."

SeifertMatrix::usage="Interior function"

fSignature::usage=""
fSeifertJ::usage=""
fSignatureJ::usage=""
fSeifertM::usage=""
fSignatureM::usage=""
fDeterminantM::usage=""

fSignat::usage="fSignat[Ul_] calculates signature of a KL 
given by its Conway symbol, Dowker code, or pdata. 
Its basic function ssmW is written by Stepan Orevkov."

fOrientedLink::usage="fOrientedLink[Ul_] calculates Gauss 
codes for a link projection given by Conway symbol, 
Dowker code, or pdata. The Gauss codes obtained correspond 
to different orientations of components (the first part of data 
obtained). The second part gives the orientations of components 
(where + is denoted by 1, and - by 0), and the third part is the 
(signed) linking number of the corrsponding oriented link."

fGaussExtSigns::usage="fGaussExtSigns[Ulaz_] calculates 
Gauss code with signs for a link projection 
given by its Conway symbol, Dowker code, or pdata."

fGraphInc::usage="fGraphInc[Ul_] calculates from a 
Conway symbol, Dowker code, or pdata of a KL, the 
corresponding graph of KL given by edges (unordered pairs), 
and by the list of vertex signs."

fPrimeGraph::usage="fPrimeGraph[GG_List] tests that a graph given 
by a list of unordered pairs, after transforming it 
into the corresponding alternating KL, will result in a  
prime or composite alternating KL. The output is 1 for a prime, 
and 0 for a composite KL."

fPlanarEmbKL::usage="fPlanarEmb[Ul_] calculates the planar 
embedding of a prime KL given by Conway symbol, pdata or 
Dowker code. An output is the list that consists of the graph of 
input KL, its planar embedding, and faces of planar 
embedded graph. As the basis of 
this program it is used the external program 
planarity.exe written by J.M.Boyer."

fPlanarEmbGraph::usage="fPlanarEmbGraph[LP_] gives the planar embedding 
of a 3-connected planar graph given by a list of unordered pairs. 
An output is the list that consists of the graph 
given by a list of unordered pairs, its
planar embedding given by vertex cycles, 
and faces of planar embedded graph. As the basis of 
this program it is used the external program planarity.exe written by \
J.M.Boyer."

DrawPlanarEmbKL::usage="DrawPlanarEmbKL[Ul_] draws the planar embedding 
of a prime KL given by Conway symbol, Dowker code, or pdata, without 
showing digonal faces. As the basis of this program it is used the 
program 3-Dimensional Convex Drawings of 3-Connected Planar Graphs 
by M.Ochiai, N.Imafuji and N.Morimura."

DrawPlanarEmbGraph::usage="DrawPlanarEmbGraph[Ul_] draws the planar 
embedding of a 3-connected graph given by the list of unordered pairs.  
As the basis of this program it is used the program 3-Dimensional 
Convex Drawings of 3-Connected Planar Graphs by M.Ochiai, N.Imafuji 
and N.Morimura."

DrawPlanarEmbGraphNew::usage=""

fMidEdgeGraph::usage="fMidEdgeGraph[L_List] gives the graph 
defined by mid-edge points of a given polyhedral graph G 
given by a list of unordered pairs of vertices."

fKLfromGraph::usage="fKLfromGraph[L_List] gives the KL
defined by a graph G given as a list of unordered pairs."

fKLinGraph::usage="fKLinGraph[UOPair_List] gives the graph 
defined by mid-edge points of a polyhedral graph G 
given by a list of unordered pairs of vertices."

fAddDig::usage="fAddDig[GL_List] produces from a given graph 
all 4-regular non-isomorphic graphs by replacing single edges by  
double (digonal) edges."

fSignsKL::usage="fSignsKL[Dow_List] calculates Dowker code with signs 
of a KL given by its  Dowker code in Knotscape form. 
For an alternating knot, an input is Dowker code without signs, 
and for a nonalternatng KL, an input is Dowker code containing 
only signs of points with signs changed with regard to the 
corresponding alternating KL."

fGraphKL::usage="fGraphKL[Ulaz_] calculates and draws graph  
of a KL given by Conway symbol, Dowker code, or pdata.  
The result is the graph given by the list of unordered pairs."

fForGraphKL::usage=""

fGraphKLNew::usage=""
fDualGraphKLNew::usage=""

fGenerators::usage="fGenerators[Ul_] calculates 
generators of a KL and the list of the signs corresponding 
to crossing points for a KL given by Conway symbol,  
Dowker code, or pdata. The first generator incommes to  
a vertex, the seccond is a the outgoing generator, and the  
third is the incomming generator (IOP=Incomming-Outgoing-Passing).  
The result is divided according to the components of KL."

fGeneratorsMirr::usage="Generators of the mirror image of a KL"

fDet::usage=""

fColTest::usage="fColTest[Ul_,cn_Integer] calculates from a given 
Conway symbol, Dowker code,  pdata, and from a number of colors 
(greater then 2), a coloring of KL. The result is the list of 
generators, list of their labelling, and the
list of generator colors."

fColTestMirr::usage="Coltest of the mirror image of a KL"


fColNo::usage=""

fColNoMirr::usage=""

fSchubertBridges::usage=""

fPDataFromDow::usage="fPDataFromDow[Dow_List] calculates pdata 
from Dowker code in Knotscape format (without signs for an 
alternating prime KL, or from Dowker code with signs of changed 
crossings for nonalternating prime KL)."

fPDataFromDowker::usage="fPDataFromDowker[Dow_List] calculates pdata 
from Dowker code with signs."

fKnotscapeDow::usage="fKnotscapeDow[Con_String] calculates from a 
Conway symbol of a KL its Dowker code in a the Knotscape format: 
Dowker code without signs for alternating KL, or Dowker code with 
signs of changed crossings for nonalternating KL."

fPrimeKL::usage="fPrimeKL[Ul_] checks that KL given by Dowker code, 
or pdata is a prime or composite KL, i.e. a direct product of some 
prime KLs. The result is 1 for a prime KL, and 0 for a composite KL."

fTorusKL::usage="fTorusKL[a_Integer,b_Integer] calculates for a 
torus KL [a,b] its braid word, min. number of crossings, unknotting 
number or number of components, bridge number, Alexander polynomial 
and (Murasugi) signature (see: Murasugi K., Knot Theory and its 
Applications, Birkhauser, Boston, Basel, Berlin, 1996)."

ShowTorusKL::usage=""

fComponentNo::usage="fComponentNo[Ul_] calculates the number of 
components of KL given by Conway symbol, Dowker code, or pdata."

fBreakComp::usage="fBreakComp[Ul_,k_Integer] in a link given by 
Conway symbol, Dowker code, or pdata, and by ordering number 
k of a cutted component calculates the pdata of a link with 
k-th component cutted."

BreakCoAll::usage="BreakCoAll[Ul_] in a link given by its Conway 
symbol, Dowker code, or pdata, cutts all components. The results 
are all different KLs obtained by cutting all components, where 
for unknot (unlink) the result is {0}."

CuttNo::usage="CuttNo[Ul_] calculates the cutting number of a  
link given by its Conway symbol, Dowker code, or pdata."

NoSelfCrossNo::usage="U-infinity number is the minimum number of 
cuts in self-crossing points of a KL necessary to obtain unknot, 
or a link without self-crossings. The function NoSelfCrossNo[Ul_] 
calculates the U-infinity number (see: Jablan S.: Unknotting number 
and Infty-Unknotting Number of a Knot, Filomat (Nis), 12:1 (1998), 
113-120) of a  KL given by its Conway symbol, Dowker code, or pdata."

fCuttRealKL::usage="fCuttRealKL[UL_] calculates the number of 
real cuttings of a given projection of KL given by its Conway 
symbol, Dowker code, or pdata, i.e., the number of cuttings 
with a different cutting point in the projection. The result is 
the number of different real cutting classes with preserved 
signs, or with preserved or reversed signs."

SplittNo::usage="SplittNo[PData_] calculates splitting number 
of a link given by its Conway symbol, Dowker code, or pdata, i.e., 
the minimum number of crossing changes necessary in order to obtain a \
splitted 
link. It is illustrated by the example from Adams. (see: Adams C: 
Splitting Versus Unlinking, Journal of Knot Theory and Ramifications, 
Vol.5, No. 3 (1996) 295-299)."

AmphiProjAltKL::usage="AmphiProAltjKL[Ul_] tests the amphicheirality of 
a given projection of an alternating KL given by its Conway symbol, 
Dowker code, or pdata."

AmphiQ::usage="AmphiQ[Ul_List] tests the amphicheirality of 
KL given by its pdata."

PeriodProjAltKL::usage="PeriodProjAltKL[Ul_] calculates the period 
of a given projection of an alternating knot given by its 
Conway symbol, Dowker code, or pdata."

PeriodAltKL::usage="PeriodAltKL[Conway_String] calculates the 
period of a given alternating KL given by its Conway symbol."

fAmphiAltKL::usage="fAmphiAltKL[Conway_String] tests the 
amphicheirality of an alternating KL diagram given by its Conway symbol."

AmphiAltKL::usage="AmphiAltKL[Conway_String] tests the 
amphicheirality of an alternating KL given by its Conway symbol."

Symm::usage="Symm[Ul_] calculates all automorphisms of an 
alternating KL given by its Conway symbol, Dowker code, 
or pdata. The first part of the result is the list of 
automorphisms given by permutations, and the other the 
list of the corresponding cycles."

MaxSymmProjAltKL::usage="MaxSymmProjAltKL[Con_String] finds maximum 
symmetrical projection of an alternating KL given by its Conway symbol."

fBasicPoly::usage="fBasicPoly[Ul_List] finds the basic polyhedron 
for a KL given by Dowker code, or pdata."

fCompositePoly::usage="fCompositePoly[Conway_String] checks that basic 
polyhedron given by its Conway symbol is elementary or composite.  
The result is 1 for composite, and 0 for elementary basic polyhedron."

fPolyFlype::usage="fPolyFlype[Conway_String] checks that basic 
polyhedron given by its Conway symbol permits flypes or not.  
The result is 1 for basic polyhedra permiting flypes, and 0 for the others."

fAllStatesProj::usage="fAllStatesProj[UL_] calculates all states of 
a given alternating KL, i.e., all different variations of signs in 
a given KL projection. Projection can be given by its Conway symbol, 
Dowker code, or pdata. The result is the list of reduced KL obtained 
from it by all sign-changes. The first datum is  of crossing changes, 
and the second the list of all KL obtained."

fAllStPrk::usage=""

RationalKL::usage="RationalKL[n_Integer] calculates the number and 
Conway symbols of all rational KLs for a given number of crossings n."

RationalAmphiK::usage="RationalAmphiK[n_Integer] calculates the number 
and Conway symbols of all rational amphicheiral knots for a given 
number of crossings n."

RationalAmphiL::usage="RationalAmphiL[n_Integer] calculates the number 
and Conway symbols of all rational amphicheiral links for a given 
number of crossings n."

RatGenSourKL::usage="RatGenSourKL[n_Integer, m_Integer] calculates the 
number and Conway symbols of all rational generating KLs (for m=3) 
and rational source KLs (for m=2) with n crossings."

RatSourceKLNo::usage="RatSourceKLNo[n_Integer] calculates the 
number of rational source KLs (for k=2) with n crossings 
according to the general recursion formula: b[0]=1, b[1]=1, 
b[2n-2]+b[2n-1]=b[2n], b[2n]+b[2n-1]-f[n-1]=b[2n+1], where 
f is the Fibonacci sequence given by recursion 
f[0]=1, f[1]=1, f[n-2]+f[n-1]=f[n]."

RatLinkU1::usage="RatLinkU1[n_Integer] calculates the number 
and Conway symbols of all rational links with the unlinking 
number 1 with n crossings."

RatLinkU0::usage="RatLinkU0[n_Integer] calculates the number 
and Conway symbols of all rational unlinks (U=0) with n crossings."

RatKnotGenU0::usage="RatKnotGenU0[n_Integer] gives Conway symbols 
of all rational links with the unlinking number 0 with n crossings."

RatKnotGenU1::usage="RatKnotGenU1[n_Integer] gives Conway symbols 
of all rational knots with the unknotting number 1 with n crossings."

MSigRat::usage="MSigRat[Conway_String] calculates the fraction 
and Murasugi signature of a rational KL given by its Conway symbol."

AllStatesRational::usage="AllStatesRational[Conway_String] calculates 
all states of a given rational KL given by its Conway symbol, i.e., 
all different variations of signs in a given projection of a rational 
KL. The first datum in the result is the unknotting (unlinking) 
number obtained from the fixed projection of a given rational KL,  
followed by  the list of minimum numbers of crossing changes necessary 
to obtain the corresponding rational KL given by its Conway symbol, 
i.e., the list of the KL distances from a given rational KL projection."

fAlexPolyOne::usage="Computes one-variable AlexanderPolynomial of a link"

AllStPR::usage=""

AllStatesRatFast::usage="Interior function"

fDirectProduct::usage="Interior function"

fDToDDirect::usage="fDToDDirect[Conway_String] computes Dowker code 
of a direct product of KLs given in Conway notation."

fDToDDirectPD::usage="fDToDDirectPD[Conway_String] computes p-data 
of a direct product of KLs given in Conway notation."

fGenSignDirProd::usage="fGenSignDirProd[Conway_String] computes 
signs of a direct product of KLs given in Conway notation."

fConvertDirect::usage="Interior function"

fClassicToCon::usage="fClassicToCon[Ulaz_String] gives Conway 
symbol of a KL given in classical notation (according to D. Rolfsen)."

fConToOther::usage="fConToOther[Con_String] gives for a Conway symbol of
a KL its classical (Rolfsen) and Dowker-Thistlethwaite notation."

fAllCodes::usage=""

fDowThistToCon::usage=""

fConNotation::usage="fConNotation[Ul_String] gives Conway symbol 
of a KL given in clasical (Rolfsen) or Dowker-Thistlethwaite notation."

fGetBraidRepresent::usage="fGetBraidRepresent[Ul_] calls the 
external program Braids-9.0 by A. Bartholomew. It returns a B-word of 
a braid representation for the 
given Conway symbol or P-data."

fGaussExtSignsBraid::usage=""

BraidWReduced::usage="BraidWReduced[PD_List] gives a 
reduced braid word of a KL given by its pdata."

fEdmonds::usage="fEdmonds[UlL_List] gives list of polygons, 
characteristic and genus for a graph given as a list of 
unordered pairs."

fStellarBasic::usage="fStellarBasic[n_Integer] produces the list 
of basic stellar KLs with n crossings."

fStellar::usage="fStellar[n_Integer] produces the list 
of stellar KLs with n crossings."

fCompl::usage="Interior function"
fStelRot::usage="Interior function"
fGener::usage="Interior function"
fMakeStel::usage="Interior function"

fPlanarEmbNew::usage=""

fPlanarEmb::usage="fPlanarEmb[neu_List] calculates 
the planar embedding of a prime KL given by 
Conway symbol, pdata or  Dowker code. An output is 
the list that consists of the graph of  input KL, 
its planar embedding given  by vertex cycles, and 
faces of planar embedded graph. As the basis of  
this program it is used the external MS DOS 
program planarity.exe written by J.M.Boyer."

fStellarPlus::usage="fStellarPlus[n_Integer] produces the list 
of distinct stellar KLs with n crossings and with pluses."

fStellarNalt::usage="fStellarNalt[n_Integer] produces the list 
of nonalternating stellar KLs with n crossings."

NMoveRat::usage="NMoveRat[n_Integer,Conway_String] calculates 
for a given integer n and Conway symbol of a KL the minimum 
number of n-moves neccessary to unknott a minimal projection 
of a rational KL."

fJablanPoly::usage="fJablanPoly[Conway_String] calculates 
multivariable Jablan polynomial for a projection of KL 
given by its Conway symbol."

fGener::usage="Interior function"

fGapRat::usage="fGapRat[Conway_String] finds rational KLs with 
an unlinking gap. The result is the unlinking gap, followed by 
the unlinking number of KL and the unlinking number of its 
fixed minimal projection."

fUnKLFixed::usage="fUnKLFixed[UL_,n_Integer] checks that a 
fixed projection of KL can be unknotted (unlinked) by n crossing 
changes. The result is {} for a KLs that cannot be unknotted 
(unlinked) by n crossing changes, and n for the others."

fGap::usage="fGap[Conway_String] finds KLs with the unlinking gap. 
The first part of result is 1 for KLs with an unlinking gap, 
and 0 otherwise, and second part is unlinking number."

fUnRFixProj::usage="fUnRFixProj[Conway_String,k_Integer] 
calculates the minimum number of crossing changes 
necessary to unlink a given fixed KL projection."

UnRFixPr::usage=""

ShowKnotfromPdataNew::usage="ShowKnotfromPdataNew[Ul_List, s_:7, opts___] \
shows KL using
Mathematica 6.0 or Mathematica 7.0 3D graphic."
ShowKnot3DNew::usage=""
ShowKnot3DNewTop::usage=""


fAutoSigInp::usage="fAutoSigInp[Gr_List] calculates stable 
states of an automaton, given by a list of outgoing edges, 
signs of vertices and inputs in vertices. In the vertices 
with a sign 1 is used the operation NOR, and in the 
vertices -1 the operation OR. The result is the list of 
edge colorings corresponding to stable states and the 
list of stable states according to vertices (see Kauffman). 
If the signlist is empty, it is computed as {1,...,1}, 
i.e. with NOR in all vertices."

fAutoKL::usage="fAutoKL[Ul_List] calculates the stable 
states of an automaton obtained from a KL given in Conway 
notation, followed by a list of signs of vertices, and 
inputs in vertices. In the vertices with a sign 1 is 
used the operation NOR, and in the vertices -1 the 
operation OR. The result is the list of edge colorings 
corresponding to stable states and the list of stable states 
according to vertices (see Kauffmann). If the signlist 
is empty, it is computed as the signlist of the given KL."

LiangPoly::usage="LiangPoly[Conway_String] calculates Liang 
polynomial of a KL projection given by its Conway symbol 
and distinguishes amphicheiral 
KL projections with a Liang polynomial satisfying the 
relationship L(t)=L(1/t)."

LiangPoly1::usage=""

fBoolean::usage="fBoolean[LL_] calculates the square-free 
polynomial for a given logical term."

fDiffSeq::usage="fDiffSeq[n_Integer,p_Integer] calculates 
all different periodic sequences with n variables, 
of period p. As the result, every sequence is given 
by its code, Kauffman code, and square-free polynomials 
corresponding to it."

fBalanced::usage="fBalanced[Diff_List,redBr_Integer] 
calculates stable (balanced) states for 
a given periodic sequence."

fDowCodes::usage="fDowCodes[n_Integer] 
generates all different knot and link projections
with n crossings, including non-prime KLs. It uses
external program plantri.exe by B.McKay and G.Brinkmann."

fGenKL::usage="fGenKL[b_Integer, pp_String] generates 
all different alternating KLs with n crossings 
from a given source KL."

fMinDowKnotPr::usage=""

fWrGraph::usage=""

fPlanEmbedding::usage=""

fForknotFind::usage=""

fKnotFind::usage=""

fKnotFindOld::usage=""

fKnotFindML::usage=""

fKnotFindDirProd::usage=""

fGapKnotsc::usage=""

fVarNewP::usage=""

RedKauffmanPolynomial::usage=""

RedKauffman::usage="RedKauffman[Ulaz_] computes normalized two-variable 
Kauffman polynomial (up to factor) of a KL given by P-data or Conway symbol."

RedAlex::usage="RedAlex[Ulaz_] computes normalized Alexander 
polynomial (up to factor) of a KL given by P-data or Conway symbol."

RedAlexander::usage=""

RedCon::usage="RedCon[Ulaz_] computes normalized Conway 
polynomial (up to factor) of a KL given by P-data or Conway symbol."

RedJones::usage="RedJones[Ulaz_] computes normalized Jones 
polynomial (up to factor) of a KL given by P-data or Conway symbol."

RedHom::usage="RedHom[Ulaz_] computes normalized HOMFLYPT 
polynomial (up to factor) of a KL given by P-data or Conway symbol."

fJon::usage=""

fHom::usage=""

fKauff::usage=""

fKh::usage=""

fCheckGap::usage=""

fCheckSign::usage=""

fSpecial::usage=""

fProdTangles::usage="fProdTangles[s1_String, s2_String, tt_String] 
gives product of non-algebraic tangles m5*,m7*,m81*-m82*,m91*-m96*,
m101*-m1011*,m111*-m1138* with an algebraic tangle tt placed in it."

fSumTangles::usage="fSumTangles[s1_String, s2_String, tt_String] 
gives sum of non-algebraic tangles m5*,m7*,m81*-m82*,m91*-m96*,
m101*-m1011*,m111*-m1138* with an algebraic tangle tt placed in it."

ListOfOneFactors::usage:="ListOfOneFactors[n] gives the list of all 
possible one-factors in a cubic graph with a given Hamilton cycle 
on n vertices. Some combinations can be isomorphic."

fDiffViae::usage="fDiffViae[n_Integer] gives all different
viae with n mirrors"

fViaToKL::usage="fViaToKL[LL_List] from given via 
computes Dowker code of corresponding KL."

ShowChord::usage="ShowChord[Ul_] shows chord diagram"

AdmissibleEdge::usage=""

ConnectedComponents::usage =""

DFS::usage =""

fDualGraphKL::usage ="fDualGraphKL[Ul_] gives the dual graph of a KL graph
for a KL given by its Conway symbol, Dowker code or P-data."

fConwayToPD::usage =""
fConToPD::usage =""
fConwayToPDDirProd::usage =""
fKnotscapeDowToPD::usage =""
fDowkerToPD::usage =""
fPdataToPD::usage =""

fPDfromBW::usage = "
fPDfromBW[Ul_String] computes for the braid word given as a 
string its PD (planar diagram)"

BraidWord::usage = "
  BraidWord[L] is an alphabetic string which is a braid word 
  presentation for 
  the named link L. That is, it is a presentation of a braid whose closure 
  (within the miminum-possible braid group) is L."


PlanarDiagramFromSignedIntegerBraidRepresentative::usage = "
  PlanarDiagramFromSignedIntegerBraidRepresentative[SIBR[n, sib]] 
  is a planar diagram constructed from a signed integer braid."

StringIndexOfBraidWord::usage = "
  StringIndexOfBraidWord[bw] is the string index of the given braid word. \
This 
  is not the braid index of the braid (of the associated link) formed by its 
  closure."

SoberBraidWordQ::usage = "
  SoberBraidWordQ[bw] is True if bw is a well-formed braid word. 
  It tolerates 
  the unused \"\"."

StringIndexOfBraidWord::usage = "
  StringIndexOfBraidWord[bw] is the string index of the given braid word. \
This 
  is not the braid index of the braid (of the associated link) formed by its 
  closure."


SignedIntegerBraidWidth::usage = "
  SignedIntegerBraidWidth[bw]."


BraidWordFromSignedIntegerBraid::usage = "
  BraidWordFromSignedIntegerBraid[sib]."


SignedIntegerBraidFromBraidWord::usage = "
  SignedIntegerBraidFromBraidWord[bw]."
  
LinksGould::usage="LinksGould[m,n,L] computes 
Links-Gould invariant for a given braid L, and 
the parameters m,n."

fBRtoBW::usage=""

fBraidW::usage="fBraidW[Ul_] gives braid word of
KL given by its Conway symbol, Dowker code with
signs, or pdata."

LinksGouldInv::usage="LinksGouldInv[m,n,L] computes 
Links-Gould invariant for a KL given by its 
Conway symbol, Dowker code with signs, or pdata, and 
the parameters m,n."

fKauffAlg::usage="The function fKauffAlg[n_Integer, p_Integer] 
produces all different periodic sequences of the 
period p with n variables."

fViaToKL::usage="The function fViaToKL[LL_List]for every curve 
given by its (uncolored) chord diagram finds its basic prime KL."

ListOfOneFactors::usage="ListOfOneFactors[n] gives the list 
of all possible one-factors in a cubic graph with a 
given Hamilton cycle on n vertices. Some combinations 
can be isomorphic. This function is written by Tamara 
Bertok and corrected by S.Jablan."
  
AdmissibleEdge::usage=""

fDifSeq::usage=""

ShowTorusKLNew::usage=""
 
ZerosAlex::usage="ZerosAlex[Ul_String] gives sum of absolute values 
and list plot of the zeros of Alexander polynomial of a KL 
given by its Conway symbol."

ZerosJones::usage="ZerosJones[Ul_String] gives sum of absolute values 
and list plot of the zeros of Jones polynomial of a KL 
given by its Conway symbol." 

fDiffViae::usage="fDiffViae[n_Integer] for a given number n 
derives all different self-avoiding curves
with n mirrors that can be obtained from prime 
KLs with n crossings."

fStellarNinv::usage="fStellarNinv[n_Integer] 
derives non-invertible stellar knots with n crossings"

fTangleType::usage="fTangleType[Ul_], where Ul is string
or list, computes type of R-tangle, giving as the result 
0 for atangle of type [0], 1 for tangle of type [1], 
and 2 for tangle of type [infty]"

fMakeType::usage="For given n fMakeType[n_Integer,k_Integer] 
gives all R-tangles of type k=0,1,2, where 2 corresponds 
to R-tangles of type [infty]"

fStelTypeTest::usage="fStelTypeTest[Ul_String] checks 
type condition for non-invertible stellar knots. 
Complete check of non-ivertibility gives symmetry test"

fNinvStellar::usage="fNinvStellar[n_Integer] derives all 
non-invertible stellar knots with n crossings."

fNinvStellarCorr::usage=""

fMulTan::usage="fMultTan[ll_List, ll1_List] multiplies two chord diagrams."

fMulTanTab::usage="fMultTanTab[Ul_List] multiplies diagrams
from the input set and gives the multiplication table
of the output set."

fGenSet::usage="fGenSet[Ul_List] checks that a given set of 
diagrams is the generator set of all n-diagrams, gives the number
of the elements of the set obtained by 
multiplying generators, and if it is a generator set 
gives 1 as the third part of the result. Otherwise, it gives 0."



fTutte::usage ="fTutte[ll_List] computes Tutte polynomial
of a graph given by the list of unordered pairs of vertices."

fCorr::usage=""

fMakeTan::usage =""

fMinTangle::usage =""

fMakeTangle::usage =""

fTanglefromCode::usage =""

fMakeEdges::usage =""

fRegions::usage =""

fRegList::usage =""

fAllClosures::usage ="fAllClosures[n_Integer] gives list of 
all possible closures of n-tangle."

fClosedTangle::usage =""

fBasicTan::usage ="fBasicTan[n_Integer,k_Integer] gives the list 
of all basic tangles derived from n-tangle with k crossings
and the ordering numbers of their closures. "

fInConnect::usage =""

fGrEdgGen::usage =""

fDowfromBP::usage =""

fDowfromCode::usage =""

fBasicPolyTan::usage ="fBasicPolyTan[k_Integer, n_Integer] gives the list 
of basic polyhedra as k-tangles for given k 
and given number of crossings n"

fMakeBP::usage ="fMakeBP[LL_List, 
n_Integer, s_Integer, k_Integer, l_Integer] 
from given list of tangles LL makes all basic polyhedra with 1,2 
sequence of the length k, tangle from the position s, 
and 1,2 or 2,1 sequence 
of length l"

fMakeAllnsBP::usage ="fMakeAllnsBP[LL_List,n_Integer,s_Integer,c_Integer] 
from s-sequence from the list of minimal tangle codes computes
all basic polyhedra with c crossings derived from the 
corresponding n-tangle.
"
fChangeCrossing::usage =""

fVarP3::usage =""

fVarP::usage =""

fKnotscapeDowfromDow::usage =""

fPolyNorm::usage ="fPolyNorm[Ul_] normalizes different polynomials."

VirtKnotTab::usage ="VirtKnotTab[n_Integer] gives n-th virtual
knot from the table of virtual knots which contains 2999 virtual knots
derived from knots with at most 8 crossings."
VirtLinkTab::usage ="VirtLinkTab[n_Integer] gives n-th virtual
knot from the table of prime virtual links which contains 3687 virtual links
derived from links with at most 8 crossings."
fVirtPD::usage ="fVirtPD[Con_String] computes PD (planar diagram)
for a virtual KL given by its Conway symbol."
fOddWrithe::usage ="fOddWrithe[Ul_String] computes odd writhe of a virtual \
knot 
given by its Conway symbol."
fFindSignsPart::usage =""
fFindSigns::usage =""
fSawollek::usage ="fSawollek[Ul_String] computes Sawollek polynomial 
of a virtual KL given by its Conway symbol."
fSawollekNorm::usage ="fSawollek[Ul_String] computes normalized Sawollek \
polynomial 
of a virtual KL given by its Conway symbol."
ShowVirtKL::usage ="ShowVirtKL[Ul_String, s_:7, opts___] shows virtual KL as \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
a 3D
graphic using Mathematica 6.0 or Mathematica 7.0."
fGaussExtSignsDirPrKnot::usage =""
fGaussVirtKnot::usage ="fGaussVirtKnot[Ul_String] computes Gauss code 
of a virtual knot given by its Conway symbol."
fAlexVirt::usage ="fAlexVirt[Ul_String] computes normalized Alexander \
polynomial
of a virtual KL given by its Conway symbol."
fJonesVirt::usage ="fJonesVirt[Ul_String] computes normalized Jones \
polynomial
of a virtual KL given by its Conway symbol."
fCabledJonesVirt::usage ="fCabledJonesVirt[Ul_String,n_Integer] computes 
normalized n-cabled Jones polynomial of a virtual knot given by its Conway \
symbol."

fGaussVirtPD::usage=""
fCabledJonesVirtPD::usage=""

fFindConVirt::usage="fFindConVirt[Ul] gives standard (according to our \
tables) Conway symbol 
of a virtual knot given by its Conway symbol or PD."

(* fKauffmanPD::usage="" *)


fKamadaCode::usage="fKamadaCode[Ul_String] computes Kamada code of a 
virtual knot given by its Conway symbol."
fMiyazawaPoly::usage="fMiazawaPolu[Ul_String] computes Miazawa polynomials
of a virtual knot given by its Conway symbol."
fBothPoly::usage="fBothPoly[pp_] computes both Jones polynomials of a \
(virtual) KL 
and its mirror image."
fGaussExtSignsDirPrKnot::usage=""
fFindSignsPart::usage=""
fFindSigns::usage=""
fMakeStrPart::usage=""
fLabeledImmersionCode::usage="fLabeledImmersionCode[Ul_String] computes
labeled immersion code of a virtual knot given by its Conway symbol."
fBraidVirt::usage="fBraidVirt[Ul_String] computes virtual braid corresponding
to a virtual knot given by its Conway symbol."
fSawollekBraid::usage="fSawollekBraid[Ul_String] computes Sawollek polynomial
of a virtual knot given in Conway notation by using its virtual braid."
fRealCrossingChange::usage=""
fRealUnknottingNo::usage="fRealUnknottingNo[Ul_String] computes real \
unknotting
number of a virtual knot with at most 12 crossings given by its Conway \
symbol."
fVirtualCrossingChange::usage=""
fVirtualUnknottingNo::usage="fVirtualUnknottingNo[Ul_String] computes virtual \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
unknotting
number of a virtual knot with at most 12 crossings given by its Conway \
symbol."
fMixedCrossingChange::usage=""
fMixedUnknottingNo::usage="fMixedUnknottingNo[Ul_String] computes mixed \
unknotting
number of a virtual knot with at most 12 crossings given by its Conway \
symbol."
fVirtCrossChangesFix::usage="fVirtCrossChangesFix[Ul_String] computes minimal \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

number of virtual crossing changes for unknotting fixed virtual knot diagram
given by its Conway symbol."
fRealCrossChangesFix::usage="fRealCrossChangesFix[Ul_String] computes minimal \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

number of real crossing changes for unknotting fixed virtual knot diagram
given by its Conway symbol."
fUnknottingVirt::usage="fUnknottingVirt[Ul_String] computes virtual fixed, 
virtual, real fixed, and real unknotting number of a virtual 
knot given by its Conway symbol."
framing2::usage=""
expandpower::usage=""
ShowBraidNew::usage="ShowBraidNew[Ul_,n_:4,s_:5] shows braid as a 3D graphic
using Mathematica 6.0 or Mathematica 7.0."
framed::usage=""
fBasVirt::usage=""
fMakePolyhedral::usage=""
fMakeAllNonalt::usage=""
fMakeConVirt::usage="fMakeConVirt[Ul_String] derives all virtual knots from
a knot given by its Conway symbol."
Cable::usage="Cable[n_][L_PD] computes n-cable for a knot given by its PD."
fAdequateMixed::usage=""
fAdequateSignedMixed::usage="fAdequateSignedMixed[Ul_,s_List] computes state \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
circles of a KL given 
in Conway notation with a sequence of markers."
fAllStates::usage="fAllStates[Ul_] computes all state graphs of a KL given in \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
Conway notation."
fGenerate::usage=""
fGenSign::usage=""
fGenSignsPD::usage=""
fGenSignsPDNew::usage=""
fGenSignfromPD::usage=""
fGenSignfromPDNew::usage=""
fGenSignVirtfromPD::usage=""
fReduceReal::usage=""
fClearDouble::usage=""
fReduceGraphsVirt::usage=""
fDelVirtLoops::usage=""
fFinalGraphsVirt::usage=""
fKauffmanVirt::usage="fKauffmanVirt[Ul_String] computes Kauffman
extended bracket of a virtual KL given in Conway notation."
fGenSignfromPDAdditional::usage=""
fKauffmanVirtfromPD::usage="fKauffmanVirtfromPD[Ul_,vv_List] computes \
Kauffman extended bracket
for a KL given by its PD and the list of virtual crossings."
fAllCycles::usage=""
fAdequate::usage=""
fCycleNo::usage=""
fAdequateSigned::usage=""
fAdequateTest::usage="fAdequateTest[Ul_] checks adequacy of a KL diagram
given by its Conway symbol."
fReducedStateGraphs::usage=""
fReducedAdequateStateGraphs::usage=""
fCriticalLines::usage="fCriticalLines[Ul_] computes critical lines of
Kauffman polynomial of a KL given by its Conway symbol."
fAdequateMixed::usage=""
fAdequateSignedMixed::usage="fAdequateSignedMixed[Ul_]"
fAdequateTestMixed::usage="fAdequateTestMixed[Ul_,s_List] checks adequacy of \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
state graph."
fAdequateMixedFast::usage=""
fAdequateSignedMixedFast::usage=""
fMixedStateGraph::usage=""
fAllAdequateStatesGraphs::usage=""
fStateGraph::usage="fStateGraph[ll_List] computes + and - state graphs of a \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
KL given in Conway notation."
fStateGraphReducedOneDigon::usage=""
fAllAdequateStatesFast::usage="fAllAdequateStatesFast[Ul] computes all \
adequate state 
graphs of a KL given in Conway notation."
fAllStates::usage=""
fEdgeReduction::usage=""
fBasicGr::usage=""
fReducedGraph::usage=""
fFindGen::usage=""
graphHom::usage=""
fDataforPari::usage=""
fMakeKn::usage=""
fMakeMat::usage=""
fTorsionList::usage=""
fGraphTorsion1::usage=""
fGraphTorsion2::usage=""
fKauffPolyAdqChrTor1::usage="fKauffPolyAdqChrTor1[Con_,m_Integer] computes \
adequacy 
polynomial of a KL given in Conway notation using first homology."
fMakeMap::usage=""
graphHom1::usage=""
fC2Gen::usage=""
fMakeMap1::usage=""
graphHom2::usage=""
fKauffPolyAdqChrTor2::usage="fKauffPolyAdqChrTor1[Con_,m_Integer] computes \
adequacy 
polynomial of a KL given in Conway notation using second homology."
fMakeAllPoly::usage=""
UnFix::usage=""
UnKLFix::usage=""
fAllMoves::usage=""
fExpandCon::usage=""
fMakePolyhedral::usage=""
fMakeConV::usage=""
fStandardPoly::usage=""
fDyePoly::usage="fDyePoly[Ul_String] computes k-cabled Jones polynomials
of a virtual KL given by its Conway symbol."
fMakeConVirtLink::usage="fMakeConVirtLink[Ul_String] gives all virtual links 
derived from a link given by its Conway symbol."
fMakeNonalt::usage=""
fAllNonaltKnotsfromAlt::usage="fAllNonaltKnotsfromAlt[Ul_Strring] derive
all non-alternating knots from an alternating knot given by its
Conway symbol."
fAllNonaltLinksfromAlt::usage="fAllNonaltLinksfromAlt[Ul_Strring] derive
all non-alternating links from an alternating link given by its
Conway symbol."
fMakeCycle::usage=""

fFindConway::usage="fFindConway[Ul_List,type_Integer] gives Conway symbol
of a KL with at most 12 crossings given by its (unreduced) Conway symbol \
(type=1), P-data (type=2), 
Knotscape DT code (type=3), PD (type=4), or Doeker code with signs (type=5)."

fFindDirProd::usage=""

fRecCon::usage=""
fRecPdata::usage=""
fRecKnotsc::usage=""
fRecPD::usage=""
fRecDow::usage=""

fForFindConway::usage=""
fFindKSketch::usage=""
fVirtLinkSketch::usage=""



fPDToPData::usage=""
fDrawingPD::usage=""
fDrawingPDNew::usage=""

ShowVirtKLPD::usage=""

fCheckGap::usage=""

fCheckSign::usage=""

fSpecial::usage=""

fKhoHoHom::usage="fKhoHoHom[Ul_String] computes reduced Khovanov homology 
and thickness of a KL given in Conway notation or by P-data by using the \
program 
KhoHo by A.Shumakovitch."
fKhoHoOddHom::usage="fKhoHoOddHom[Ul_String] computes odd reduced Khovanov \
homology 
and thickness of a KL given in Conway notation or by P-data by using the \
program 
KhoHo by A.Shumakovitch."
fQuasialtTest::usage=""
fOddHomologyTor::usage=""
fQuasiAlt::usage="fQuasiAlt[Ul_] checks that a KL (diagram) given by its \
Conway symbol 
is quasi-alternating or not."


X::usage=""
Y::usage=""
del::usage=""
led::usage=""
B::usage=""
F::usage=""
fKauffmanExtendedBracket::usage=""
fKauffmanBracket::usage=""
fKauffmanBracket1::usage=""
fKauffmanArrow::usage="fKauffmanArrow[Ul_String] 
computes Kauffman arrow polynomial of a virtual KL given by its Conway \
symbol."
RawBracket::usage=""
fReplacement::usage=""
fReplacement1::usage=""
fKauffmanPD::usage=""
J::usage=""
K1::usage=""
K2::usage=""
K3::usage=""
fMakeConVirtLinkKauff::usage="fMakeConVirtLinkKauff[Ul_String] produces a \
list
of virtual KLs derived from a KL given by its Conway symbol 
by using Kauffman extended bracket."



fDrawMirrorCurveCol::usage="fDrawMirrorCurveCol[Ul_List] draws a colored mirror curve given by its code, where positive crossings 1 are blue, and negative -1 red." 

fDrawMirrorCurve::usage="fDrawMirrorCurveCol[Ul_List] draws a black-white mirror curve given by its code."

fAnaMirror::usage="" 
fAnaMirrorVirtPD::usage="" 
fGaussVirtMirrCurve::usage=""
fCabledJonesVirtMirrCurve::usage=""
fShowVirtKLMirrCurve::usage=""
fVirtPDfromPD::usage=""
fVirtLSketch::usage=""
fMinLinkSketcher::usage=""
fShowMinLinkSketcher::usage=""
fDrawFromLinkSketcher::usage=""

            
fVarP4::usage=""


fMirrSquare::usage="fMirrSquare[n_Integer] for a given n generates all codes of mirror curves in a square with a side n."

fMirrRect::usage="fMirrRect[p_Integer,q_Integer] for a given p and q  generates all codes of mirror curves in a rectangle wirh sides p and q."

fCompNo::usage="fCompNo[Ul_List] for a mirror curve without crossings, containing only two-sided mirrors 2 or -2 computes number of components."

fPolyNorm::usage="fPolyNorm[Ul_] normalizes different polynomials."

fBracketAlt::usage="fBracketAlt[Ul_List] for an alternating mirror curve computes its Kauffman bracket polynomial exressed in terms of x^4."

fBracket::usage="fBracket[Ul_List] for an arbitrary mirror curve computes its Kauffman bracket polynomial expressed in terms of x^4."

fMirrProduct::usage="fMirrProduct[ll_List,ll1_List] for two mirror-curves ll and ll1 of the same dimensions computes their product."

fRepresent::usage="fRepresent[Ul_List] represents every mirror-curve given by its binary code as the product of two unlinks (the 
first and the second part of the outoput), i.e., mirror curves containing in their codes only 2 and -2, and the third part is the code of the input 
mirror-curve, originally given by its binary code."

fDecompose::usage="fDecompose[Ul_List] gives the binary code of input mirror-curve and decomposes it into its projections: mirror curves 
expresed in terms of 2 and -2. In the code the first two numbers are dimensions of the rectangle p, q, and the other two the binaru coding of the mirror-curve."

fDecomposeFast::usage="fDecompose[Ul_List] decomposes every mirror curve into its projections: mirror curves expresed in terms of 2 and -2."

fProjfromCode::usage="fProjfromCode[Ul_List] computes the mirror-curve from its binary code."

fMinimize::usage="fMinimize[Ul_List] gives the binary code of a mirror-curve. The first two data are dimensions of the rectangle p and q, 
and the other two binary codes of the mirror-curve."


fDowMorwen::usage="" 
fMorwenDow::usage=""
fKnotscapeDT::usage=""
fAllToPD::usage=""
fAllToDT::usage=""
fAllToDTLet::usage=""
fCrNo::usage=""

fDTNotationToDTLet::usage=""
fClassicalNotationToDTLet::usage=""
fSnapPyInp::usage=""
fSnapPyInput::usage=""
fSnapPyUl::usage=""
fForSnapPyDT::usage=""
fForSnapPyBR::usage=""
fSnapPyVol::usage=""
fForSnapPyAmphi::usage=""
fSnapPyAmphi::usage=""
fForSnapPyInv::usage=""
fSnapPyInv::usage=""
fForSnapPyFundGr::usage=""
fSnapPyFundGr::usage=""
fForSnapPyChSim::usage=""
fSnapPyChSim::usage=""
fForSnapPyIsom::usage=""
fSnapPyIsom::usage=""
fPDfromLS::usage=""
fDTfromSnapPy::usage=""
fPdatafromSnapPy::usage=""
fForSnapPySymmGr::usage=""
fSnapPySymmGr::usage=""


fCrossingChange::usage=""
fCrossingChangeExt::usage=""
fDowfromMor::usage=""
fUnKnotFast::usage=""
fNextStep::usage=""
fUnknotting::usage=""
fUnknotting1::usage=""
fForUnknottingDraw::usage=""
fUnknottingDraw::usage=""
fUnRat::usage=""
fUnRatDraw::usage=""
fUnSignat::usage=""

fUnLink::usage=""
fUnLink1::usage=""
fPomDrawUnLink::usage=""
fDrawUnLink::usage=""
fPomForDrawUnLink::usage=""
fForDrawUnLink::usage=""
fLowerBound::usage=""
fKnotsc::usage=""


fHyp001::usage=""
fHyp002::usage=""


fPyramid::usage=""
fPrism::usage="" 
fBipyramid::usage="" 
fAntiPrism::usage=""
fAntiBipyramid::usage=""
fCupola::usage=""
fCupolaDual::usage=""
fElongatedCupola::usage=""
fGyroElongatedCupola::usage=""
fElongatedRotund::usage=""
fGyroElongatedRotund::usage=""
fPrismDrum::usage=""
fGyroDrum::usage=""
fOrthoBicupola::usage=""
fGyroBicupola::usage=""
fPrismAntiprism::usage=""
fElongatedPyramid::usage=""
fElongatedPrism::usage=""
fElongatedBipyramid::usage=""
fElongatedOrthoBicupola::usage=""
fElongatedGyroBicupola::usage=""
fBicupolaAntiprism::usage=""
fTruncation::usage=""
fCyc::usage=""
fPolyPlot::usage=""
fPolyPlotDoubleEdges::usage=""
fPolyPlotKL::usage=""
fPolyFunctions1::usage=""
fPolyFunctions2::usage=""
fPolyFunctions4::usage=""
fPolyFunctions41::usage=""
fMidEdgeDT1::usage=""
fMidEdgeKL1::usage=""
fMidEdgePlot1::usage=""
fMidEdgeGraphPlot1::usage=""
fMidEdgeDT2::usage=""
fMidEdgeKL2::usage=""
fMidEdgePlot2::usage=""
fMidEdgeGraphPlot2::usage=""
fTruncDT1::usage=""
fTruncKL1::usage=""
fTruncPlot1::usage=""
fTruncGraphPlot1::usage=""
fTruncDT2::usage=""
fTruncKL2::usage=""
fTruncPlot2::usage=""
fTruncGraphPlot2::usage=""
fMidEdgeDT0::usage=""
fMidEdgeKL0::usage=""
fMidEdgePlot0::usage=""
fMidEdgeGraphPlot0::usage=""
fTruncDT0::usage=""
fTruncKL0::usage=""
fTruncPlot0::usage=""
fTruncGraphPlot0::usage=""
fPolyTransf::usage=""
fPolyTransfJaeger::usage=""
fPolyTransfInput::usage=""
fMidEdgeTransfDT0::usage=""
fMidEdgeTransfKL0::usage=""
fMidEdgeTransfPlot0::usage=""
fMidEdgeTransfGraphPlot0::usage=""
fTruncTransfDT0::usage=""
fTruncTransfKL0::usage=""
fTruncTransfPlot0::usage=""
fTruncTransfGraphPlot0::usage=""
fTruncDT::usage=""
fMidEdgeDTBasic::usage=""
fMidEdgeKLBasic::usage=""
fMidEdgePlotBasic::usage=""
fMidEdgeGraphPlotBasic::usage=""
fTruncDTBasic::usage=""
fTruncKLBasic::usage=""
fTruncPlotBasic::usage=""
fTruncGraphPlotBasic::usage=""

fJaeger::usage=""
fPolyTransfMidInput::usage=""
fPolyTransfTruncInput::usage=""
fPolyTransfJaegerInput::usage=""
fJaegerDT1::usage=""
fJaegerKL1::usage=""
fJaegerPlot1::usage=""
fJaegerGraphPlot1::usage=""
fJaegerDT2::usage=""
fJaegerKL2::usage=""
fJaegerPlot2::usage=""
fJaegerGraphPlot2::usage=""
fJaegerTransfDT0::usage=""
fJaegerTransfKL0::usage=""
fJaegerTransfPlot0::usage=""
fJaegerTransfGraphPlot0::usage=""
fJaegerDTBasic::usage=""
fJaegerKLBasic::usage=""
fJaegerPlotBasic::usage=""
fJaegerGraphPlotBasic::usage=""
fMidEdgeIteration::usage=""
fMidEdgeIteration1::usage=""
fMidEdgeIteration2::usage=""
fMidEdgeIterationBasic::usage=""
fPlantriPolytopes::usage=""
fPlantriTri::usage=""
fDoubleFull::usage=""
fDoubleFull1::usage=""
fDoubleEdgesFull::usage=""
fFullList::usage=""

fOneParPoly::usage="" 
fTwoParPoly::usage=""
fAllPolyMid::usage=""
fAllPolyCross::usage=""
fAllPolyJaeger::usage=""
fBasPoly::usage=""
fGraph3DKL::usage=""
fFourVal1::usage=""
fFourVal2::usage=""
fIterativeMid1::usage=""
fIterativeMid2::usage=""
fIterativeMid3::usage=""
fIterativeBasic::usage=""
fSurfaceKL::usage=""
fOneParPolyMidSurfaces::usage=""
fOneParPolyCrossSurfaces::usage=""
fTwoParPolyMidSurfaces::usage=""
fTwoParPolyCrossSurfaces::usage=""
fAllPolyMidSurfaces::usage=""
fAllPolyMidTruncSurfaces::usage=""
fAllPolyCrossSurfaces::usage=""
fPlantriMid::usage=""
fPlantriCross::usage=""
fGeneralFull::usage=""
fGeneralFullDouble::usage=""
fFullgen::usage=""
fVirtPlanarLink::usage=""
fVirtPlanar::usage=""
fVirtLinkGraph::usage=""
fVirtGraph::usage=""
IsNotKnot::usage=""


Begin["`Private`"]



(*## ## ## ## # Implementation of +,",","_",- #### ## ## ## ## #*)

fRef[{L1_List,{{a_,b_},{a1_,b1_}},L3_List}]:={L1,{{a,a1},{b,b1}},-L3}
fBlock[x_,n_Integer]:={{},{{4n-3,4n},{4n-2,4n-1}},{1}}
fInv[{L1_List,{{a_,b_},{a1_,b1_}},L3_List}]:={L1,{{a,b},{a1,b1}},-L3}
fPlus[L_List]:=
  Module[{s,i,l,p,d,p1,d1},l=Length[L];p=L[[1]][[2]][[1]];d=L[[1]][[2]][[2]];
    p1=L[[2]][[2]][[1]];d1=L[[2]][[2]][[2]];
    If[l>2,s=L[[1]];  Do[s=fPlus[{s,L[[i]]}],{i,2,l}],
      s={Join[Join[L[[1]][[1]],
              L[[2]][[1]]],{{p[[2]],p1[[1]]},{d[[2]],d1[[1]]}}],{{p[[1]],
              p1[[2]]},{d[[1]],d1[[2]]}}, Join[L[[1]][[3]],L[[2]][[3]]]}];
              s
    ] 
fComma[L_List]:=fPlus[Map[fRef,L]]
fSpace[L_List]:=Module[{s,i,l},l=Length[L];
    If[l>2,s=L[[1]];  Do[s=fSpace[{s,L[[i]]}],{i,2,l}],
      s=fPlus[{fRef[L[[1]]],L[[2]]}]];
    s
    ] 
    
 (* ## ## ## ## Conversion from Conway Notation ## ## ## ## ## ## *) 
   
fConvert[ sInput_String, opts___Rule ] :=
  Module[ {sOut,sCh,iI,iJ,iLevel,iPos,iPrev,iNext,lComm} ,
     sCh = Characters[ sInput ] ;
    lComm = {};     
    For[ iI=1, iI<=Length[sCh], iI++,  
      If[ sCh[[iI]]//DigitQ ,        
        While[ iI+1<=Length[sCh] && DigitQ[sCh[[iI+1]]] ,
           sCh[[iI]] = sCh[[iI]] <> sCh[[iI+1]];
           sCh       = Drop[ sCh , {iI+1} ];
        ];
        AppendTo[ lComm , 4 ];  
      , AppendTo[ lComm ,
          Switch[ sCh[[iI]],
               "(" , 2  ,          
               ")" , -2 ,          
               "+" , 5  ,         
               " " , 6  ,         
               "," , 7  ,         
               "-" , 8  ,          
               _   , 0            
          ]
        ];
      ];
    ];
    
    If[ Position[ lComm, 0 ] =!= {} , Return[False] ]; 
    iLevel = 0;               
    Do[        
      If[ Abs[lComm[[iI]]]===2 , iLevel += lComm[[iI]] ];
    ,{iI,Length[sCh]}];
    If[ iLevel =!= 0 , Return[ False ] ];  
   
    For[ iPos=iI=1, iI<=Length[sCh], iI++,  
      If[ lComm[[iI]] === 4 ,              
        sCh[[iI]] = "fBlock[" <> sCh[[iI]] <> "," <> ToString[iPos] <> "]"; 
        iPos++;                                  
        lComm[[iI]] = 1 ;       
      ];
    ];
   
    If[ (Print /. {opts})===On,
      Print["blocks: ",Thread[{sCh,lComm}]//MatrixForm];
    ];
     
    For[ iI=Length[sCh]-1, iI>0, iI--, 
      If[ lComm[[iI]]===8 ,
         If[ lComm[[iI+1]] == 1 ,
            (* single block *)
            sCh[[iI]]   = "fInv[" <> sCh[[iI+1]] <> "]";
            lComm[[iI]] = 1;
            sCh   = Drop[ sCh   , {iI+1} ];
            lComm = Drop[ lComm , {iI+1} ];
         , sCh[[iI]]   = "fInv[" <> sCh[[iI+1]];
            lComm[[iI]] = 2;
            sCh    = Drop[ sCh   , {iI+1} ];
            lComm  = Drop[ lComm , {iI+1} ];
            iLevel = 0;
            iJ     = iI+1;
            While[ iLevel>-2 ,   
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            sCh[[iJ-1]]   = sCh[[iJ-1]] <> "]";
         ];
      ];
    ];
  
    If[ (Print /. {opts})===On,
      Print["operation '-': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
 
    For[ iI=2, iI<Length[sCh], iI++, 
      If[ lComm[[iI]]===6,
        iPrev = iI-1;              
        If[ lComm[[iI-1]]===-2 , 
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "fSpace[{", iPrev ];
        lComm =  Insert[ lComm, 2         , iPrev ]; 
        iNext = iI+1;  
       
        While[ iNext<=Length[sCh] && lComm[[iNext]]===6,
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 9  ;
          iNext++;              
          If[ lComm[[iNext]]===2 ,    
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ-1;
          ];
          iNext++;   
        ];
        sCh   =  Insert[ sCh  , "}]" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ]; 
      ];
    ];
    
    If[ (Print /. {opts})===On,
      Print["operation ' ': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
    
    For[ iI=2, iI<Length[sCh], iI++,  
      If[ lComm[[iI]]===7,
      
        iPrev = iI-1;              
        If[ lComm[[iI-1]]===-2 ,   
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "fComma[{", iPrev ];
        lComm =  Insert[ lComm, 2         , iPrev ]; 
         iNext = iI+1;  
       
        While[ iNext<=Length[sCh] && lComm[[iNext]]===7,
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 9  ;
          iNext++;              
          If[ lComm[[iNext]]===2 ,  
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ-1;
          ];
          iNext++;   
        ];
        sCh   =  Insert[ sCh  , "}]" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ]; 
      ];
    ];
   
    If[ (Print /. {opts})===On,
      Print["operation ',': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
    
    For[ iI=2, iI<=Length[sCh], iI++,
      If[ lComm[[iI]]===5,
         iPrev = iI-1;             
        If[ lComm[[iI-1]]===-2 ,   
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "fPlus[{", iPrev ];
        lComm =  Insert[ lComm, 2       , iPrev ];
        iNext = iI+1;   
        
        sCh[[iNext]]   = ",";
        lComm[[iNext]] = 9  ; 
        iNext++;            
        If[ lComm[[iNext]]===2 ,   
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ-1;
        ];
        iNext++;    
        sCh   =  Insert[ sCh  , "}]" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ];
      ];       
    ];
   
    If[ (Print /. {opts})===On,
      Print["operation '+': ",Thread[{sCh,lComm}]//MatrixForm];
    ];

    sOut = StringJoin[ sCh ];
    If[ SyntaxQ[sOut]===False, Return[ False ]];    
    ToExpression[ sOut ]                            
     ]
     
 (*## # Conversion from Conway Symbols for Polyhedral Knots and Links ## ##*)
     
fFormVList[Out_List]:= Module[{ll = {}, i},
    Do[Switch[Out[[i]],
               -1,ll=Append[ll,{4i-2,4i-1,4i,4i-3}],
                1,ll=Append[ll,{4i-3,4i,4i-1,4i-2}]]
         , {i, 1, Length[Out]}];
    ll
        ]
        
fMoja[{x_?NonNegative,y_?NonNegative}, br_]:={x+br,y+br}
fMoja[{x_?NonNegative,y_?Negative}, br_]:={x+br,y}
fMoja[{x_?Negative,y_?NonNegative}, br_]:={x,y+br}
fMoja[{x_,y_}, br_]:={x,y}

fConvertPoly[Conway_String]:=
  Module[{OutG=0,str,Con=Conway,Vert,br=0,SignList={},el,min,p,ll,pos,x6,x61,
      x8},pos=Flatten[StringPosition[Con,"*"]];
    If[pos=={},
      If[StringTake[Con,1]==".",OutG="x61";Con=StringDrop[Con,1],
        OutG="x6"],OutG="x"<>StringTake[Con,{1,pos[[1]]-1}];
      If[StringLength[Con]==pos[[1]],Con="",
        Con=StringTake[Con,{pos[[1]]+1,StringLength[Con]}]]];
    OutG=ToExpression[OutG];
    If[Head[OutG]==List,Vert=fFormVList[OutG[[3]]];
      If[Con!="",
        If[StringTake[Con,1]==".",Con=StringInsert[Con,"1",1]]];
      While[StringPosition[Con,".."]!={},
        Con=StringReplace[Con,".."->".1."]];
      If[Con=="",Do[Con=Con<>"1.",{i,1,Length[OutG[[3]]]}],
        pos=Union[Flatten[StringPosition[Con,"."]]];
        Do[Con=Con<>".1",{i,1,Length[OutG[[3]]]-1-Length[pos]}];
        Con=Con<>".";
        If[StringTake[Con,1]==".",Con=StringDrop[Con,1]];];
      pos=Union[Flatten[StringPosition[Con,"."]]];
      Do[str=StringTake[Con,{1,pos[[i]]-1}];
        Con=StringDrop[Con,{1,pos[[i]]}];
        pos=pos-pos[[i]];
        (*Print["str ",str];*)If[
          StringTake[str,-1]=="0"&&StringPosition[str," \
"]!={},
          el=fRef[fConvert[StringDrop[str,-2]]],el=fConvert[str]];
        el[[1]]+=br*4;el[[2]]+=br*4;
        br=br+Length[el[[3]]];
        SignList=Append[SignList,el[[3]]*OutG[[3,i]]];
        min=Min[Vert[[i]]]+4;
        p=4*Length[el[[3]]]-4;
        If[p!=0,ll=OutG[[2]];ll=Map[(#-min)&,ll];
          ll=Map[fMoja[#,p]&,ll];
          ll=Map[(#+min)&,ll];OutG[[2]]=ll];
        ll={Vert[[i]][[1]]->el[[2]][[1]][[1]],
            Vert[[i]][[2]]->el[[2]][[1]][[2]],
            Vert[[i]][[3]]->el[[2]][[2]][[2]],
            Vert[[i]][[4]]->el[[2]][[2]][[1]]};
        OutG[[2]]=ReplaceAll[OutG[[2]],ll];
        Do[Vert[[j]]+=p,{j,i+1,Length[Vert]}];
        OutG[[2]]=Join[OutG[[2]],el[[1]]],{i,1,Length[OutG[[3]]]}];
      OutG[[2]]=Sort[OutG[[2]]];
      OutG[[3]]=Flatten[SignList]];
    OutG]

(*## ## ## ## ## ## ## ## ## ## ## #*)
 
 (* ## ## ## ## ## ## ## ## ## ## ## # *)
 
fFun[k_] := Module[{l},
    If[Mod[k, 4] == 1 \[Or] Mod[k, 4] == 2, l = IntegerPart[k/4] + 1, 
      l = -(IntegerPart[(k - 1)/4] + 1)];
    l]

fGraphHa[ConSym_, OutLinks_List, n_] :=
  	Module[{InLinks = {}, pom, p, i, ModSym, KLGraph, HL},
    ModSym = Map[Mod[#, 2] &, ConSym];
    Do[If[ModSym[[i]] == 1, 
        InLinks =  Union[Flatten[Map[Append[InLinks, #] &,
         {{4i - 3, 4i - 1}, {4i - 2, 4i}}],1]], 
        InLinks =  Union[Flatten[ Map[Append[InLinks, #] &,
         {{4i - 3, 4i},{4i - 2, 4i - 1}}],1]]]
     , {i, 1, Length[ModSym]}];
    KLGraph = FromUnorderedPairs[Join[InLinks, OutLinks]];
    HL = ConnectedComponents[KLGraph];
    If[Length[HL] == 1,(*KNOT*)HL = Map[fFun[#] &, HL[[1]]]; 
      HL = Take[HL, {1, Length[HL], 2}],
       (*LINK*)pom = {};
        Do[p = Map[fFun[#] &, HL[[i]]];
          	p = Take[p, {1, Length[p], 2}];pom = Append[pom, p]
       , {i, 1, Length[HL]}];
       HL = pom];
    {ModSym, HL} ]

fFormElem[p_, k_, q_] := Module[{i, LL = {}},
    If[p > 0, Do[LL = Append[LL, i], {i, k + 1, k + q}], 
      Do[LL = Append[LL, k + q + 1 - i], {i, 1, q}]];
    LL]
(* ## ## ## ## ## ## ## fConSymPoly[Conway_String]## ## ## ## ## #*)
(*nalazi "pun" Konvejev simbol poliedra*)
(*ulaz Konvej*)
(*pomocna za DOWKERA *)

fConSymPoly[Conway_String]:=Module[{Conpom=Conway,pos,ConSym, i},
    If[Union[Flatten[StringPosition[Conpom,"*"]],
          Flatten[StringPosition[Conpom,"."]]]=={},
      (*nije poliedarski*)
      ConSym=ReadList[
          StringToStream[
            StringReplace[
              Conway,{"#"->" ",","->" ","+"->" ","-"->" ",
                "("->" ",")"->" "}]],Number]
      ,
      pos=Flatten[StringPosition[Conpom,"*"]];
      If[pos=={},
          If[StringTake[Conpom,1]==".",Conpom=StringDrop[Conpom,1]]; 
        If[Conpom!="",
                      
          If[StringTake[Conpom,1]==".",
            Conpom=StringInsert[Conpom,"1",1] ]] ;
               pol=6
           (*ako je . na prvom mestu ovo je ".1" koji je x6 " *),  
            If[pos[[1]]==2,
                pol=ToExpression[StringTake[Conpom,1]],  
                 pol=ToExpression[StringTake[Conpom,{1,2}]] 
              ];
                (*sada nam je u pol oznaka poliedra tj. 
              broj njegovih presecnih tacaka*)
            
        If[StringLength[Conpom]==1,
                      Conpom="", 
            Conpom=StringTake[Conpom,{pos[[1]]+1,StringLength[Conpom]}]];
                   (*uzimamo ostatak posle * 
            ili prve . (ako ga ima) *)
            If[Conpom!="",
                      
          If[StringTake[Conpom,1]==".",
            Conpom=StringInsert[Conpom,"1",1] ]]
        ];
      
      (*ubacujemo tackice*)
               
      While[StringPosition[Conpom,".."]!={},
        Conpom=StringReplace[Conpom,".."->".1."]];
                    (*ubacujemo 1 u preskocenim temenima *)
               
      If[Conpom=="",
                      
        Do[Conpom=Conpom<>"1.",{i,1,Length[fConvertPoly[Conway][[3]]]}], 
        	     pos=Union[Flatten[StringPosition[Conpom,"."]]];
              Do[Conpom=Conpom<>".1",{i,1,pol-1-Length[pos]}];
                Conpom=Conpom<>".";
                
        If[StringTake[Conpom,1]==".",Conpom=StringDrop[Conpom,1]]];
                ConSym=ReadList[StringToStream[StringReplace[Conpom,
              {","->" ","+"->" ","-"->" ","("->" ",
                ")"->" ","."->" "," 0"->" "}]],Number]];
    ConSym
    ]
(* ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## #*)
  (*## ## ## ## ## ## # DIRECT PRODUCT ## ## ## ## ## ## ## #*)
   (* 27.8.2003*)  
    (*## ## ## ## ## ## ## DIRECT PRODUCT ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## #*)

fRot[{U1_List,{{a_,b_},{a1_,b1_}},Z1_List}]:={U1,{{a,a1},{b,b1}},Z1}
(*NOVA 17.8.2003 *)

fDirectProduct[L_List]:=Module[{pomList,A1,B1,ind,pom,p,p1,d,d1},
    A1=L[[1]];B1=L[[2]]; 
    pomList=Rest[L];
    While[Not[SameQ[pomList,{}]],
                  ind=0;
                 If[Length[A1[[2]]]>2,ind=1];
                 If[Length[B1[[2]]]>2,If[ind==1,ind=3,ind=2]];
      If[ind==2, 
      (*zameni mesta prvom i drugom*)
       pom=A1;A1=B1;B1=pom];
       (*ako je samo 1 poly on je sada u A1*)
       (*ind oznacava da li ima poliedarskih:
       0-nema;1,2 ima jedan;3 oba poly*)
      p=A1[[2,1]];  d=A1[[2,2]];
      p1=B1[[2,1]];   d1=B1[[2,2]];
      If[ind==0, (*ne poly*)
        A1={Join[A1[[1]],B1[[1]],{{p[[2]],p1[[1]]},{d[[2]],d1[[1]]}}],
                                 {{p[[1]],d[[1]]},{p1[[2]],d1[[2]]}},
                                 Join[A1[[3]],B1[[3]]]},
                      (*jedan ili vise poliedarskih*)
                        \
 A1={Join[A1[[1]],B1[[1]]],
                                 Sort[Join[Rest[A1[[2]]],Rest[B1[[2]]],
                {Sort[{p[[1]],p1[[2]]}],Sort[{p[[2]],p1[[1]]}]}]],
                                  Join[A1[[3]],B1[[3]]]}];
                  (*u A1 je trenutni rezultat*)
                   
      pomList=Rest[pomList];
                   If[Not[SameQ[pomList,{}]],B1=First[pomList]]
      ] (*end of While *);
    A1
    ]
(*## ## ## ## ## ## ## ## ## ## ## ## ##*)
fConvertDirect[Conway_String]:=Module[{pos,pom,ind,llCon={},i,j=0},
    pos=Flatten[Map[Take[#,1]&,StringPosition[Conway,"#"]]];
    ind=SameQ[
        Union[StringPosition[Conway,"*"],StringPosition[Conway,"."]],{}];
    (* Dowker obezbedjuje da ovo nije prazna lista*)
    Do[If[i==1,
                pom=StringTake[Conway,{1,pos[[1]]-1}],
                 If[i==Length[pos]+1,  
                     
          pom=StringTake[Conway,{Last[pos]+1,StringLength[Conway]}],
                    pom=StringTake[Conway,{pos[[i-1]]+1,pos[[i]]-1}]
                 ]];
      If[Union[StringPosition[pom,"*"],StringPosition[pom,"."]]=={},
        If[ind,pom=fRot[fConvert[pom]],
          pom=fConvert[pom]],
        pom=fConvertPoly[pom]];
      pom={pom[[1]]+j,pom[[2]]+j,pom[[3]]};
         j=j+4*Length[pom[[3]]];
           llCon=Append[llCon,pom]
      ,{i,1,Length[pos]+1}];
    llCon=fDirectProduct[llCon];
    llCon
    ]
(*## ## ## ## ## ## ## ## #  fConSymDirect ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## # *)
(*pravi "pun" Konvejev simbol za direktne proizvode *)
(*pomocna za Dowkera*)

fConSymDirect[Conway_String]:=Module[{pom,pos,llCon={},i},
    pos=Flatten[Map[Take[#,1]&,StringPosition[Conway,"#"]]];
    Do[If[i==1,
                pom=StringTake[Conway,{1,pos[[1]]-1}],
                 If[i==Length[pos]+1,  
                pom=StringTake[Conway,{Last[pos]+1,StringLength[Conway]}],
                    pom=StringTake[Conway,{pos[[i-1]]+1,pos[[i]]-1}]
                 ]];
          pom=fConSymPoly[pom];
          llCon=Append[llCon,pom]
      ,{i,1,Length[pos]+1}];
    Flatten[llCon]
    ]
 (*## ## ## ## ## ## ## ## ## ## ## # fGenSignDirProd ## 
#*)
 
fGenSignDirProd[Conway_String] := Module[{pp},
    pp = Map[Sign, fDToDDirect[Conway][[2]]];
    pp
    ]

(*## ## ## ## ## ## ## ## fDToDDirect ## ## ## ## ## ## ## ## ## ## #*)

fDToDDirect[Conway_String] := 
  Module[{pos, llCon = {}, i, vv = {}, vv1 = {}, ll, pom, pom1, pom2, tt, 
      tt1}, pos = Flatten[Map[Take[#, 1] &, StringPosition[Conway, "#"]]];
    Do[If[i == 1, pom = StringTake[Conway, {1, pos[[1]] - 1}], 
        If[i == Length[pos] + 1, 
          pom = StringTake[Conway, {Last[pos] + 1, StringLength[Conway]}], 
          pom = StringTake[Conway, {pos[[i - 1]] + 1, pos[[i]] - 1}]]];
      pom = Dow[pom];
      pom1 = pom[[1]];
      pom2 = pom[[2]];
      llCon = Append[llCon, pom];
      vv = Append[vv, pom1];
      vv1 = Append[vv1, pom], {i, 1, Length[pos] + 1}];
     tt = Table[vv1[[i, 1]], {i, Length[vv1]}];
    tt1 = Table[tt[[i, -1]] + tt[[i + 1, 1]], {i, Length[tt] - 1}];
    tt1 = 
      Join[Table[
          ReplacePart[tt[[i]], tt1[[i]], -1], {i, Length[tt1]}], {tt[[-1]]}];
    tt1 = 
      Flatten[Join[{{tt1[[1, 1]]}}, 
          Table[Drop[tt1[[i]], 1], {i, Length[tt1]}]]];
    ll = 2Table[Length[vv1[[i, 2]]], {i, Length[vv1]}];
    vv = ll[[1]];
    ll = Join[{ll[[1]]}, Table[vv = vv + ll[[i]], {i, Length[ll] - 1}]] - 
        ll[[1]];
    ll = Flatten[
        Table[Table[
            vv1[[j, 2, i]] + Sign[vv1[[j, 2, i]]]*ll[[j]], {i, 
              Length[vv1[[j, 2]]]}], {j, Length[vv1]}]];
    If[Apply[Plus, tt1] < Length[ll], 
      tt1 = ReplacePart[tt1, (tt1[[-1]] + vv1[[-1, 1]])[[1]], -1] , tt1];
    ll = {tt1, ll};
    ll]
    
fDToDDirectPD[Conway_String] := Module[{pp, pp1, pp2, sl, ll, i},
    pp = fDToDDirect[StringReplace[Conway, "-" -> ""]][[2]];
    pp2 = fDToDDirect[Conway];
    pp1 = pp2[[2]];
    sl = Map[Sign, pp1];
    ll = Map[Sign, pp]*sl;
    pp = Sort[
        Table[If[
            SameQ[ll[[i]], 1], {2i - 1, Abs[pp[[i]]], sl[[i]]}, \
{Abs[pp[[i]]],
               2i - 1, sl[[i]]}], {i, Length[pp]}]];
    pp = {pp2[[1]], Table[pp[[i, 2]]*pp[[i, 3]], {i, Length[pp]}]};
    pp
    ]
    

    
    
    
 (* ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## #*)
 
 (* ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## #*)
 
 
    (*## ## ## ## ## ## ## ## ## DOWKER ## ## ## ## ## ## ## ## ## ## ##*)
   (* 16.8.2003*)
(*## ## ## ## ## ## ## ## # DOWKER ## ## ## ## ## ## ## ## ## ## ## ## # *)

Dowker[Conway_String]:=
  Module[{ConList,Conpom,LinInd=0,L,h,n,ConSym=0,DL={},LDL={},LenList,LL={},
      PomList,lparova={},pol,LenDowList={},Ham,l,pos,k=0,pom=0,ll={},ll1={},
      ind,r,ListforSigns,q, i,j},
    Conpom=Conway;
    If[StringPosition[Conpom,"#"]!={},
            (*ne poledarski DIRECT PRODUCT*)
              
      ConList=fConvertDirect[Conpom];  ConSym=fConSymDirect[Conway]
      (*ConSym nije odgovarajuci za poliedre vidi crvemno dole *)
      ,
      If[Flatten[
            Union[Map[
                Flatten[StringPosition[Conpom,#]]&,{"*","."}]]]=={},
                (*ne poledarski Prime*)
                  
        ConList=fConvert[Conway];ConSym=fConSymPoly[Conway],
                 (*poliedri*)
                   ConList=fConvertPoly[Conway];
        ConSym=fConSymPoly[Conway]
          ]
      ];(*kraj prvog if-a od Direct*)
    (*Print["sym ",ConSym];
      Print["Convert ",ConList];*)
    (*ConSym={3,3};
      ConList={{{3,6},{4,5}},{{1,2},{7,8}},{1,1}};
      ConList=ReplacePart[ConList,{1,1,1,1,1,1,1,1,1},3];*)
    
    If[Head[ConSym]==List,n=Length[ConSym];
      L=fGraphHa[ConSym,Join[ConList[[1]],ConList[[2]]],n];
      DL={ConSym,L[[1]],L[[2]]};
      LenList=Map[Length[#]&,L[[2]]];
      Do[pom=pom+LenList[[i]];LL=Append[LL,pom],{i,1,Length[LenList]}];
      PomList=Flatten[Abs[L[[2]]]];k=0;
      Do[pom=0;
        Do[pom=pom+ConSym[[PomList[[k+j]]]],{j,1,LenList[[i]]}];
        LenDowList=Append[LenDowList,pom];
        k=k+LenList[[i]],{i,1,Length[L[[2]]]}];
      pom=0;
      Do[pom=pom+LenDowList[[i]];
        LDL=Append[LDL,pom],{i,1,Length[LenDowList]}];
      l=Length[LenList];
      If[Head[L[[2]][[1]]]==List,LinInd=1];
      Ham=Flatten[L[[2]]];
      k=0;
      Do[q=ConSym[[Abs[Ham[[i]]]]];
        ll=Append[ll,fFormElem[Ham[[i]],k,q]];
        k=k+q,{i,1,Length[Ham]}];
      pom=LenList[[1]];
      Ham=Abs[Ham];
      Do[ind=False;
        For[j=i+1,j<=Length[Ham]&&Not[ind],j++,
          ind=(Ham[[i]]==Ham[[j]]);
          If[ind,lparova=Append[lparova,{i,j}]];
          r=Position[LL,j-1];
          If[r!={},r=r[[1,1]]+1];
          
          If[LinInd==1&&IntegerQ[r]&&EvenQ[ll[[i,1]]-ll[[j,1]]]&&ind,
            Do[ll[[k]]=Map[(#-1)&,ll[[k]]],{k,LL[[r-1]]+1,LL[[r]]}];
            Part[ll,LL[[r-1]]+1,1]=LDL[[r]]]],{i,1,Length[Ham]-1}];
      Do[i=lparova[[h]][[1]];
        j=lparova[[h]][[2]];
        Do[
          Switch[Mod[ll[[i]][[t]],2],0,
            ll1=Union[ll1,{{ll[[j]][[t]],ll[[i]][[t]]}}],1,
            ll1=Union[ll1,{{ll[[i]][[t]],ll[[j]][[t]]}}]],{t,1,
            Length[ll[[i]]]}],{h,1,Length[lparova]}];
      ListforSigns=ll;
      ll1=Sort[Partition[Flatten[ll1,2],2]];
      ll1=Flatten[Map[Drop[#,1]&,Sort[ll1]]];PomList={};
      If[LinInd==1,LDL=Map[(#/2)&,LDL];
        Do[
          If[i==1,PomList=Append[PomList,Take[ll1,{1,LDL[[1]]}]],
            PomList=Append[PomList,Take[ll1,{LDL[[i-1]]+1,LDL[[i]]}]]],{i,1,
            Length[LDL]}];
        ll1=PomList];
      DL=Append[DL,ll1];
      DL=Append[DL,ListforSigns]];
    DL](*Vraca {} ako je greska*)

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
(* ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ##*)

fGenSign[Conway_String]:=
  Module[{BlockSign={},i,DL,par,fConSignList,Ham,VertSign,Dow={}},
    DL=Dowker[Conway];
    Dow=If[Not[SameQ[StringPosition[Conway,"#"],{}]],fGenSignDirProd[Conway],
        Dow=If[DL!={},Dow=Flatten[DL[[4]]];
            VertSign=DL[[5]];
            Do[BlockSign=Append[BlockSign,0],{i,1,Length[DL[[1]]]}];
            
            If[Flatten[
                  Union[Map[
                      Flatten[StringPosition[Conway,#]]&,{"*",
                        "."}]]]=={},fConSignList=fConvert[Conway][[3]],
              fConSignList=fConvertPoly[Conway][[3]]];
            Ham=Flatten[DL[[3]]];
            Do[par=Flatten[Position[Abs[Ham],i]];
              
              BlockSign[[i]]:=
                Sign[Ham[[par[[1]]]]Ham[[par[[2]]]]]*fConSignList[[i]];
              VertSign[[par[[1]]]]=BlockSign[[i]]*VertSign[[par[[1]]]];
              
              VertSign[[par[[2]]]]=BlockSign[[i]]*VertSign[[par[[2]]]],{i,1,
                Length[BlockSign]}];VertSign=Flatten[VertSign];
            
            Dow=Map[(Sign[
                      VertSign[[Flatten[Position[Abs[VertSign],#]][[1]]]]])&,
                Dow]];
        Dow]]
(*Vraca {} ako je greska*)
 
fGen[Conway_String] := 
  Module[{Dow,OddDow={},proL = {}, ll={}, ulL = {}, izL = {}, pomL, 
  pomLL1,Gen={}, pom=0, i, k, p, n1, j},
    Dow=Dowker[Conway];
If[Dow!={},
    n1 = Length[Dow[[4]]];
    OddDow=Flatten[Dow[[4]]];
    OddDow=Sort[Map[(#-1)&,OddDow]];
    If[Head[Dow[[4]][[1]]] == List, 
     Do[If[i == 1, pom += Length[Dow[[4]][[i]]];
           ll = Append[ll, Take[OddDow,{1, pom}]],
           ll = Append[ll,Take[OddDow, 
           {pom+1,pom+Length[Dow[[4]][[i]]]}]]; 
        pom += Length[Dow[[4]][[i]]]], {i, 1, n1}];
     pomLL1 = {};
     Do[If[i == 1, 
        pomLL1 =Append[pomLL1,Flatten[{Dow[[4]][[n1]][[1]], 
                   Drop[Dow[[4]][[i]], 1]}]], 
        pomLL1 =Append[pomLL1,Flatten[{Dow[[4]][[i-1]][[1]],
                   Drop[Dow[[4]][[i]],1]}]]]
        ,{i,1,n1}];
       pomLL1 = Flatten[pomLL1]
    ];
    pomL = Flatten[Dow[[4]]];
    If[Head[Dow[[4]][[1]]] == Integer, pomLL1 = pomL];
    n = Length[pomL];
    Do[ proL = Append[proL, pomL[[i]]/2], {i,1,n}];
    Do[ For[j = 1, j <= n && (OddDow[[j]] != pomL[[i]] - 1), j++];
        ulL = Append[ulL, pomL[[j]]/2]  ;
        p = Mod[pomL[[i]] + 1, 2n];
        For[k = 1, k <= n && (OddDow[[k]] != p), k++];
        izL = Append[izL, pomLL1[[k]]/2]
      , {i, 1, n}];
    Gen = {};
    Do[  For[k = 1, k <= n && (proL[[k]] != i), k++];
         Gen = Append[Gen, ulL[[k]]];
         Gen = Append[Gen, izL[[k]]];
         Gen = Append[Gen, proL[[k]]]
      , {i, 1, n}]
      
      ];
       
    Gen
    ](*Vraca {} ako je greska*)
 
 Dow[Conway_String]:=Module[{dd,ll={},str="{"},
   dd=Dowker[Conway];
   If[dd!={},dd=dd[[4]];
   If[Head[dd[[1]]]==List,ll=Map[Length[#]&,dd]];
   If[ll=={}, str=str<>"{"<>ToString[Length[dd]]<>"},"
            <>""<>ToString[dd*fGenSign[Conway]]<>"}"
            ,str=str<>ToString[ll]<>","
            <>ToString[Flatten[dd]*fGenSign[Conway]]<>"}"
            ]
            
            ];
           
  If[dd=={},str="Input data is not correct. You have probably exceeded 
  the number of possible basic polyhedra for the given nuber of crossings."];
            
  str=ToExpression[str];
  str   
    ]
    
    
(* ## ## ## ## ##  Gauss Code - Extended Dowker Code ## ## ## ## ## # *)    
fGaussExt[Ul_] := 
  Module[{Dow, SL, LL, DowPair = {}, DowPair1, DowFl, i, k = 0, DL = {}, 
      pomL = {}},
    If[Head[Ul] == List,
            Dow = Abs[Ul[[2]]]; 
            If[Length[Ul[[1]]] == 1,
               LL = {},
               LL = 2*Ul[[1]]];
           SL = Flatten[Map[Sign[#] &, Ul[[2]]]]];
    If[Head[Ul] == String,
           Dow = Dowker[Ul];
           If[Dow != {},
               Dow = Dow[[4]];
               SL = fGenSign[Ul];
              If[SameQ[Flatten[Dow], Dow],
                  LL = {},
                  LL = Map[2*Length[#] &, Dow]]
        ]];
    If[Dow != {},
       DowFl = Flatten[Dow];
       Table[DowPair = Append[DowPair, {2i - 1, DowFl[[i]]}];
       DowPair = Append[DowPair, {DowFl[[i]], 2i - 1}], {i, Length[DowFl]}];
       DowPair = Union[Map[Sort[#] &, Sort[DowPair]]];
       DowPair1 =Map[Append[#, 
           SL[[Position[DowFl, Select[#, EvenQ][[1]]][[1, 1]]]]*
              Position[DowPair, #][[1, 1]]] &, DowPair];
      DL =Union[Map[Take[#, -2] &, DowPair1], 
          Map[{First[#], Last[#]} &, DowPair1]];
      DL = Flatten[Map[Take[#, -1] &, DL]];
      If[LL != {},
        Do[If[i == 1,
              k = LL[[1]]; pomL = Append[pomL, Take[DL, {1, LL[[1]]}]], 
              pomL = Append[pomL, Take[DL, {k + 1, k + LL[[i]]}]];
              k = k + LL[[i]]],
             {i, 1, Length[LL]}];  
        DL = pomL]];
    If[Flatten[DL] == DL, DL = {DL}];
    DL](*ako je greska vraca {}*)
    
(* ## ## #  Linking No. #### *)   
(* April 3. 2005 *)   

LinkingNo[Ulaz_] := Module[{GE, LinNoList = {}, d, i, j, k, p,str,Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
     If[Not[SameQ[Head[Ul], String]],
     If[Length[Ul] != 1,
          If[SameQ[Select[Ul[[2]], OddQ[#1] &], {}],
             GE = fGaussExtSigns[Ul],
             GE = Ul],
         GE = Ul],
     GE = fGaussExtSigns[Ul]];
     d = Length[GE];
     If[SameQ[GE, Flatten[GE]],
        (*ako je u pitanju cvor GE nema podlista *)
       str="Knot",
       Do[Do[
          p = Map[Sign[#] &, Intersection[GE[[i]], GE[[j]]]];
         LinNoList = Append[LinNoList, Sum[p[[k]], {k, 1, Length[p]}]]
         , {j, i+1, d}], {i, 1, d - 1}];
        str=
        Sum[Abs[LinNoList[[k]]], {k, 1, Length[LinNoList]}]/2];
      str
    ]
    
(* ## ## ## ## ## UnKnot for Rational Knots and Links ## ## ## ## ## ## ## # \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)

IsNotKnot[Conway_List]:=
  Module[{Con=Conway,Con1,Con2},
    If[SameQ[Union[Conway],{0}],(*Fractioni ne prolaze-same nule*)
      
      Con={},Con1=FromContinuedFraction[Con];
      Con2=Con;
      If[SameQ[Con1,ComplexInfinity],
      Con=FromContinuedFraction[Reverse[Con2]],
        Con=Con1];
      If[SameQ[Con,ComplexInfinity]&&
          SameQ[FromContinuedFraction[Reverse[Con2]],ComplexInfinity],Con=1,
        Con];
      Con=ContinuedFraction[Con];
      If[Length[Con]>1,
        If[Abs[First[Con]]==0,Con=Drop[Con,{1,2}],
          If[Abs[First[Con]]==1,Con[[2]]=Con[[1]]+Con[[2]];
            Con=Drop[Con,{1}]]]];
      Con=If[SameQ[Con,{}], {1},Con]
      ];
    (* Print["Is is not knot ",Con]; *) Con] (*25.9.2003*)
(*## ## ## ## ## ## ## RatReduce-samo 1 konveja ## ## ## ##*)

RatReduce[Conway_String]:=
  Module[{Con,Con1},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con=IsNotKnot[Con];
    Con=Last[Sort[{Con,Reverse[Con]}]];
    Con=Abs[If[SameQ[Con,{0}],{1},Con]];
Con=Last[Sort[{Con,Reverse[Con]}]];
 StringReplace[ToString[Con],{","->"","{"->"","}"->""}]]

(*25.9.2003*)
(*## ## ## ## ## ## ## ## #*)


Redu[Conway_List]:=
    Module[{Con,i=1,j,l,ind=False,Con1,Res={}},
      Switch[Head[Conway[[1]]],Integer,l=1,List,l=Length[Conway]];
      While[i<=l&&Not[ind],
        Switch[Head[Conway[[1]]],Integer,Con=Conway,List,Con=Conway[[i]]];
          Con=IsNotKnot[Con];
          If[Con=={}||Con=={0}||Con=={1}||Con=={-1},
                 ind=True;Res={},  
                j=1;
             While[Not[MemberQ[Res,{0}]]&&j<=Length[Con],
                       
            Res=Union[
                Append[Res,ReplacePart[Con,Con[[j]]-2*Sign[Con[[j]]],j]]];
                (* Print[Res,MemberQ[Res,{}]]; *)
                       j++]
          
          ];
        i++];
      Res];


(* April 3. 2005 - puca zabog Japanaca za 1, 2, -2 *)


UnR[Conway_String]:=
  Module[{Con,br=-1},
 If[SameQ[Conway,"1"],{0},   Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con=IsNotKnot[Con];
    If[Not[Con=={}||Con=={0}||Con=={1}||Con=={-1}], 
      	While[Con!={},br++; Con=Redu[Con]]
               ,br=0];
               br=If[fSignat[Conway]!=0,
                Range[IntegerPart[(fSignat[Conway]+1)/2],br],
  Range[1,br]];
    br;
If[SameQ[br,{}],{0},br]]
]

(* Slavik 24.11.2011 *)
    
 (* ## ## ## ## ## ## ## ## ## ## ## ## Millet ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## *)  
  
NoToLett[Ken_, Br1_, Br2_, br_] := Module[{p, Ken11 = Ken},
      If[br == 3, p = "a"];
      If[br == 5, p = "b"]; If[br == 7, p = "c"]; If[br == 9, p = "d"];
      Ken11[[Br1]] = ReplacePart[Ken11[[Br1]], p, Br2];
      Ken11
      ];

fKenSym[Ken_, Br1_Integer, Br2_Integer] := 
  Module[{el = Ken[[Br1, Br2]], Ken11 = Ken, pos, pos1, ind1, ind2  },
    pos = Flatten[Position[Ken11[[Br1]], el]]; 
    If[Length[pos] == 1, 
           pos1 = Flatten[Position[Ken11[[el]], Br1]]; 
            Ken11 = NoToLett[Ken11, Br1, Br2 + 1, pos1[[1]]],
            pos1 = Flatten[Position[Ken11[[el]], Br1]];
            ind1 = pos[[2]] - pos[[1]] == 2;
            ind2 = pos1[[2]] - pos1[[1]] == 2;
            If[Or[ind1 && ind2, Not[ind1] && Not[ind2]],
                    Ken11 = NoToLett[Ken11, Br1, Br2 + 1, pos1[[2]]];     
                    Ken11 = NoToLett[Ken11, Br1, pos[[2]] + 1, pos1[[1]]],
                    Ken11 = NoToLett[Ken11, Br1, Br2 + 1, pos1[[1]]]; 
                    Ken11 = NoToLett[Ken11, Br1, pos[[2]] + 1, pos1[[2]]]
                ] ];
                If[Length[pos]==1,pos=Append[pos,-1]];
   { Ken11,pos[[2]]}
    ]

DowToKen[Conway_String] := 
  Module[{DL, Con,Dow, VertSign, GenList, Ul, Iz,Pro, i, j,
   p, m, n,pom,Ind,ll={}, Ken},
    Switch[Conway, "1", Ken={},
                    "-1", Ken={},
                   "2",
Ken={{1,"-",2,d,2,c,2,b,2,a},{2,"-",1,d,1,c,1,b,1,a}},
                   "-2", 
Ken={{1,"+",2,b,2,a,2,d,2,c},{2,"+",1,b,1,a,1,d,1,c}},
                   "1 1",
Ken={{1,"-",2,d,2,c,2,b,2,a},{2,"-",1,d,1,c,1,b,1,a}},
                    "-1 -1", 
Ken={{1,"+",2,b,2,a,2,d,2,c},{2,"+",1,b,1,a,1,d,1,c}},
                   "1 -1", 
Ken={{1,"-",2,c,2,b,2,a,2,d},{2,"+",1,c,1,b,1,a,1,d}}, 
                    "-1 1", 
Ken={{1,"+",2,c,2,b,2,a,2,d},{2,"-",1,c,1,b,1,a,1,d}}];
  Con=StringReplace[Conway,"-"->""];
  If[Conway!="1"&&Conway!="-2"&&Conway!="2"&&Conway!="-1"
    &&Conway!="1 1"&&Conway!="1 -1"&&Conway!="-1 1"&&Conway!="-1 -1",
    DL = Dowker[Con];
    VertSign = fGenSign[Con];
    If[VertSign[[1]]==-1,VertSign=-VertSign];
    GenList = fGen[Con];
    Dow = Flatten[DL[[4]]];
    Do[ll=Append[ll,{Dow[[i]],VertSign[[i]]}], {i,1,Length[Dow]}];
    ll={};
    Do[ll=Append[ll,Dow[[i]]/2->i], {i,1,Length[Dow]}];
    GenList=ReplaceAll[GenList,ll];
    GenList=Flatten[Sort[Map[Reverse[#]&,Partition[GenList,3]]]];
    Pro= Take[GenList, {1, Length[GenList], 3}];
    Ul = Take[GenList, {2, Length[GenList], 3}];
    Iz = Take[GenList, {3, Length[GenList], 3}];
    Ken = Table[0, {i, Length[Dow]}, {j, 10}];
    
    Do[ Ken[[i, 1]] = Pro[[i]];
       Switch[VertSign[[i]]
       , 1, 
          Ken[[i, 2]] = "+", -1, Ken[[i, 2]] = "-"];
          Ken[[i, 3]] =Position[Iz,i][[1]][[1]];(*izlazni za i*)
          Switch[Ken[[i, 2]], "+", Ken[[i, 5]] = Iz[[i]], "-", 
          Ken[[i, 5]] = Ul[[i]]];
         Ken[[i, 7]] = Position[Ul,i][[1]][[1]];
         Switch[Ken[[i, 2]], "+", Ken[[i, 9]] = Ul[[i]], "-", 
                                  Ken[[i, 9]] = Iz[[i]]];
        , {i, 1, Length[Dow]}];   
    Ind=fGenSign[Conway]*
        fGenSign[StringReplace[Conway,"-"->""]]*
                      Flatten[DL[[4]]];
    Ind=Flatten[Map[Position[Ind,#]&,Select[Ind,#1<0&]]];
     Do[    m=Take[Ken[[Ind[[i]]]],2];
           Switch[m[[2]],"+",m[[2]]="-","-",m[[2]]="+"];
           n=Take[Ken[[Ind[[i]]]],{3,10}];
           pom=Ken[[Ind[[i]],5]];
           Ken[[Ind[[i]],5]]=Ken[[Ind[[i]],9]];
           Ken[[Ind[[i]],9]]=pom;
           If[First[fGenSign[StringReplace[Conway,"-"->""]]]>0,
              If[fGenSign[Conway][[Ind[[i]]]]>0,
                 n=RotateLeft[n,2],n=RotateRight[n,2]],
              If[fGenSign[Conway][[Ind[[i]]]]>0, 
				 n=RotateRight[n,2],n=RotateLeft[n,2]]] ;
           Ken[[Ind[[i]]]]=Join[m,n]
       ,{i,1,Length[Ind]}];
  Do[ m=-1;n=-1;
       For[j = 3, j <= 9, j = j + 2,
       If[j!=m&&j!=n,
              Ken=fKenSym[Ken, i, j][[1]];
              If[m==-1,n=fKenSym[Ken, i, j][[2]],n=fKenSym[Ken, i, j][[2]]]]]
      , {i, 1, Length[Dow]}]];
    Ken
    ]

fMillett[Conway_String] := 
  Module[{Ken, str, i, j, str1, f = "izlaz11.txt"},
       Ken = DowToKen[Conway];
       str = " " <> ToString[Length[Ken]] <> ".  00:\n";
       Do[str1 = "";
             Do[   str1 = str1 <> ToString[Ken[[i, j]]] ,
        {j, 1, Length[Ken[[i]]]}];
      str = str <> str1 <> "\n", {i, 1, Length[Ken]}];
    str = str <> "\n";
    OpenWrite[f];
     WriteString[f, str];
       Write[f];
    Close[f];
    Run["lmp1" <> " izlaz11.txt" <> 
        " < " <> "izlaz11.txt" <> " > " <> 
        "polinomi1.txt"];
    Import["polinomi1.txt"] 
    ]

    
(* ## ## ## ## ## ## Connections with knotbycomp ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## *)

fCreatePData[Conway_String] := Module[{Dow, i, l, ind, pd = {}},
If[Not[SameQ[StringPosition[Conway,"#"],{}]],pd=fDToDDirectPD[Conway],
    Dow = Dowker[Conway];
 If[Dow!={},Dow=Dow[[4]];
    l = fGenSign[Conway];
    ind = l*fGenSign[StringReplace[Conway, "-" -> ""]];
    pd = Append[pd, If[IntegerQ[First[Dow]], 
          {Length[l]}, Map[Length, Dow]]];
    pd = Append[pd, Table[Sort[Table[
               If[Part[ind, i] == 1,{2i - 1, Flatten[Dow][[i]]*l[[i]]},
                                    {Flatten[Dow][[i]], (2i - 1)*l[[i]]}] 
                            ,{i, 1, Length[l]}]][[i]][[2]]
          , {i, 1, Length[l]}]];
    pd={pd[[1]],-pd[[2]]}
    
    ]];
    
    If[Dow=={},pd={}];
    
    pd ]  
    
 fDowfromPD[Ulaz_]:= 
  Module[{pd, pd1, pd2, ll = {}, i, j, z1, p = 0, q, z = {},pdata},
  If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],pdata=Ulaz,
    pdata=ToExpression[Ulaz]];
    pd = Flatten[pdata[[2]]];
pd1 = Complement[Table[i, {i, 1, 2Length[pd]}], Abs[pd]];
    pd1 = Table[{pd1[[i]], pd[[i]]}, {i, 1, Length[pd]}];
pd1 = Join[pd1, Map[Reverse[#] &, pd1]];
pd1 = Sort[Map[(#*Sign[#[[1]]]) &, pd1]];
If[SameQ[pd1,{}],{},
    pd2 = Take[Flatten[pd1], {2, Length[Flatten[pd1]], 2}];
z1 = Select[pd2, (#1 < 0) &];
    z1 = Map[(Abs[#] -> #) &, z1];
    pd1 = Abs[pd1];
Do[p = p + 2pdata[[1, i]]; 
      ll = Append[ll, p], {i, 1, Length[pdata[[1]]]}];
    Do[z = {};
      If[OddQ[Flatten[pd1][[2ll[[i - 1]] + 2]]],
                 p = ll[[i - 1]] + 1; q = ll[[i]];
        	Table[ If[j == q, z = Append[z, q -> p],                                                     
            z = Append[z, j -> j + 1]], {j, p, q}];
        pd1 = ReplaceAll[pd1, z]
         ]
      , {i, 2, Length[ll]}];
    pd1 = Sort[Map[If[EvenQ[#[[1]]], Reverse[#], #] &, pd1]];
    pd1 = ReplaceAll[pd1, z1];
    pd1 = Sort[Map[{Abs[#[[1]]], #[[2]]} &, Union[pd1]]];
    pd1 = Flatten[Map[Take[#, -1] &, pd1]];
    z1=Union[Mod[Abs[pd1],2]];
    If[MemberQ[z1,1],fDowfromPData[Ulaz],{pdata[[1]], pd1}]]
    ] 
    
    
    fDowfromPData[Ul_List]:=Module[{ll,pp,ss,pp1,pp2,i},
    ll=Length[Ul[[2]]];
    pp=Complement[Range[2ll],Abs[Ul[[2]]]];
    ss=Union[Table[pp[[i]]*Sign[Ul[[2,i]]],{i,ll}],Ul[[2]]];
    ss=Map[Last,Sort[Table[{Abs[ss[[i]]],Sign[ss[[i]]]},{i,Length[ss]}]]];
    pp1=Table[{pp[[i]],Abs[Ul[[2,i]]]},{i,ll}];
    pp2=Map[Reverse,pp1];
    pp=Union[pp1,pp2];
    pp1=iteratedTake[Range[2ll],2Ul[[1]]];
    pp2=Mod[Table[Apply[Plus,pp[[i]]],{i,Length[pp]}],2];
    pp2=Flatten[Position[pp2,0]];
    pp2=Flatten[Table[pp[[pp2[[i]]]],{i,Length[pp2]}]];
    pp2=Table[Position[pp1,pp2[[i]]],{i,Length[pp2]}];
    pp2=Partition[Map[First,Flatten[pp2,1]],2];
    pp2=Union[
        Map[Last,
          Union[Table[
              If[pp2[[i,1]]<pp2[[i,2]],pp2[[i]],Reverse[pp2[[i]]]],{i,
                Length[pp2]}]]]];
    pp2=Flatten[
        Table[If[MemberQ[pp2,i],RotateRight[pp1[[i]]],pp1[[i]]],{i,
            Length[Ul[[1]]]}]];
    pp=Map[Last,pp];
    pp=Table[pp2[[pp[[i]]]],{i,Length[pp2]}];
    pp=Sort[Table[{pp2[[i]],pp[[i]]*ss[[i]]},{i,Length[pp]}]];
    pp=Map[Last,Table[pp[[2i-1]],{i,ll}]];
    pp={Ul[[1]],pp};
    pp]
    
    (* KORIGOVAO SLAVIK 4.04.2007 *) 
    
    
    
    
    
 (* ## ## ## ## # Unlinking or unknotting number ## ## ## ## ## ## ## #*)     \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
  
fGenerate[{LK_List, L_List}] := 
  Module[{pom, p, i, pl = L, l, rez = {}, pDa = {}, br = 1},
    Do[   
      If[MemberQ[Abs[pl], i] == False, 
        pDa = Append[pDa, {i*Sign[pl[[br]]], pl[[br]]}]; br++]
          , {i, 1, 2Length[pl]}];  
    Do[pom = pDa; l = {LK};
          pom[[i]] *= -1;
           pom[[i]] = Reverse[pom[[i]]];
          pom = Sort[pom, Abs[#1[[1]]] < Abs[#2[[1]]] &];
          (*p je oblika originalnog linka*)
            p = Flatten[Map[Take[#, -1] &, pom]];
              l = Append[l, p];
              rez = Append[rez, ReductionKnotLink[l]]
       , {i, 1, Length[pDa]}];
       rez;
     Union[rez]
    ]
    
(* calculates UnKLNo od Konveja, pdata i dowkera *)
(* April 3. 2005- puca za 1, 2, -2 zbog Japanaca*)


 
   
UnKnotLink[Ulaz_]:=
  Module[{pD=Ulaz,str="",R,NoGen=1,un0},  
   If[SameQ[Head[Ulaz],String]&&
    SameQ[Union[StringPosition[Ulaz,","],StringPosition[Ulaz,"("],
        StringPosition[Ulaz,"*"],StringPosition[Ulaz,"."],
        StringPosition[Ulaz,"+"]],{}],NoGen=UnR[Ulaz],   
    If[SameQ[Head[Ulaz],String],pD=fCreatePData[Ulaz]];
   (*  If[SameQ[Head[Ulaz],List]&&SameQ[Select[Ulaz[[2]],OddQ[#]&],{}],
       pD=fPDataFromDow[Ulaz]]; *)
If[SameQ[Head[Ulaz],List],pD=Ulaz];
      un0=MemberQ[ReductionKnotLink[pD],{}];
    If[un0==True,
          str="Unknott",
                  R=fGenerate[pD];
     While[SameQ[MemberQ[R,{{},{}}],False]&&
         SameQ[MemberQ[R,{{0},{}}],False]&&SameQ[MemberQ[R,{{0,0},{}}],
False]&&SameQ[MemberQ[R,{{0,0,0},{}}],False],  
        	R=Flatten[Map[fGenerate[#]&,R],1];
        	NoGen++]];
      NoGen=If[fSignat[pD]!=0,Range[IntegerPart[(fSignat[pD]+1)/2],NoGen],
  Range[1,NoGen]] ]; 
  NoGen] (* 20.10.2010 Radi od Con i pd *)
    
  
  
UnKnotLinkNo1[Ulaz_]:=
  Module[{pD=Ulaz,str="",R,NoGen=1,un0,vv},
  If[SameQ[Head[Ulaz],String]&&
    SameQ[Union[StringPosition[Ulaz,","],StringPosition[Ulaz,"("],
        StringPosition[Ulaz,"*"],StringPosition[Ulaz,"."],
        StringPosition[Ulaz,"+"]],{}],NoGen=UnR[Ulaz],
    If[SameQ[Head[Ulaz],String],pD=fCreatePData[Ulaz]];
    If[SameQ[Head[Ulaz],List]&&SameQ[Select[Ulaz[[2]],OddQ[#]&],{}],
       pD=fPDataFromDow[Ulaz]];
      un0=MemberQ[ReductionKnotLink[pD],{}];
     If[un0==True,
          str="Unknott",
                  R=fGenerate[pD];
                vv=If[MemberQ[R,{{},{}}]==True||
          MemberQ[R,{{0},{}}]==True||MemberQ[R,{{0,0},{}}]==True
          ||MemberQ[R,{{0,0,0},{}}]==True
          ||MemberQ[R,{{0,0,0,0},{}}]==True
          ||MemberQ[R,{{0,0,0,0,0},{}}]==True,  
        	1,0]]];
        	vv
    ](*11.11.2006 *)
  
    
   
 (* ## ## ## ## # Graphics for KnotPlot ## ## ## ## ## ## ## ## ## #*)       

fCreateGraphics[Conway_] := 
  Module[{pdata, k, i, j, pom, str = "", mm,
    filestr = "C:\\Program Files\\KnotPlot2010\\graphics.txt"},
    pdata=If[SameQ[Head[Conway],String],    
    fCreatePData[Conway], Conway];
 If[pdata!={},
    Install["DrawKnot"];
    mm = AccountingForm[ N[GetDrawData[pdata[[2]], pdata[[1]]], 6]][[1]];
    Uninstall["DrawKnot"];
    Table[Do[pom = mm[[i, j]];
             Table[str = str <> StringReplace[
          ToString[AccountingForm[pom[[k]]]], {"(" -> "-", ")" -> ""}] 
<>
               " ", {k, 3}];
                        str = str <> "\n"  , {j, 1, Length[mm[[i]]]}];
                   str = str <> "\n" , {i, Length[mm]}];
    filestr = OpenWrite[filestr];
    WriteString[filestr, str];
    Write[filestr];
    Close[filestr],
    Print["Input data is incorrect"]
    ]
    ]
    
 (* ## ## ## ## #  KnotLinkBase ## ## ## ## #*)
 
 GetKnotLink[str_String,nn_Integer]:=Module[{s,rez},
   s=ToExpression[str];
   If[Head[s]==List,
    If[nn>Length[s],
     s="You have exceeded the number of knots and links from the data base "
     <>str<>".",rez=s[[nn]];
     s="Conway symbol: "<> ToString[rez];
      ]];
  If[Head[s]==Symbol,s="Data base "<>str<> " does not exist.";rez={}];
      (*  Print[s]; *)
  rez   
]

NumberOfKL[str_String]:=Module[{s="There are ",n},
    t=str;
   n=Read[StringToStream[StringDrop[t,1]],Number];
   If[Head[ToExpression[str]]==List,
   s=s<>ToString[Length[ToExpression[str]]];
   If[StringTake[str,1]=="a",s=s<>" alternating",s=s<>" non-alternating"];
   If[n==11,s=s<>" knots with ",s=s<>" knots and links with "];
   s=s<>ToString[n]<>" crossings."];
   If[Head[ToExpression[str]]==Symbol,s="Data base "<>str<> " does not \
exist."];
s]
   
 (* ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ##*) 
 (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
 (*## ## ## ## ## ## ## # K2KC1 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
 
(*****************************************************
fProjections[string] (R.S. & M.S.)
30.05.2003. 
******************************************************
OPERATIONS:  (in the order of the execution)
-----------
o Operation "X n" --> "(X,1,1,...,1)" where "n" is a number 
  and there are "n" aces; "X" can be anything
  "(4,3) 2"  -->  "((4,3),1,1)"
  This also works with the negative number:
  "(4,3) -2"  -->  "((4,3),-1,-1)"
o Operation "X+n" --> "X,1,1,...,1" where "n" is a number and
  there are "n" aces; "X" can be anything
  "-(2,3)+3" --> "-(2,3),1,1,1"
  This also works with the negative number:  "3+-2"  -->  "3,-1,-1"
o Other combinations with " ", "+", "-" or "," are left unchanged
o Operation "X,1" (or "X,-1") --> makes multiples (different projections?)
  and "X" can be anything (separated by: a. "(", b. "," or c. beggining
  "2,1"      --> "2,1"  and  "1,2"
  "2,-1"     --> "2,-1" and  "-1,2"
  "2 (-2),1" --> "2 (-2),1"  and "1,2 (-2)"
o Operation "1,X" (or "-1,X") --> same thing ( "1,3" --> "1,3" and "3,1")

******************************************************)

(* ================================================ *)
(* = function fProjections[string]                = *)
(* = (previously fConvert2)                       = *)
(* = Output is the list of strings.               = *)
(* = The function is not checking all             = *)
(* = syntax errors that user may make!            = *)
(* = option: Print->On prints intermediate results= *)  
(* ================================================ *)
Options[ fProjections ] = { Print -> Off };
fProjKL[ sInput_String, opts___Rule ] :=
  Module[ 
    {
       sOut,sCh,lComm,iI,iJ,iK,iLevel,iPos,iPrev,iNext,iNumRead,
       sOne,tmp,ltmp,
       lAllCh,lAllComm        (* <--- we keep multiple outputs here *)
    } 
  ,

    (* we keep here splitted string, and later recognized commands *)
    sCh = Characters[ sInput ] ;

    (* ----------------------------------------------- *)
    (* - This part is recognizing the numbers and    - *)
    (* - letters, and we put the description in      - *)
    (* - the array lComm with symbols:               - *)
    (* - 0  --> not recognized (yet)                 - *)
    (* - 1  --> single block                         - *)
    (* - 2  --> open bracket    "("                  - *)
    (* - -2 --> closed bracket  ")"                  - *)
    (* - 3  --> operation "+"                        - *)
    (* - 4  --> number                               - *)
    (* - 5  --> letter "+"                           - *)
    (* - 6  --> letter " "                           - *)
    (* - 7  --> letter ","                           - *)
    (* - 8  --> letter "-"                           - *)
    (* - 9  --> separator after recognized operations- *)
    (* - -4 --> -number (e.g. "-2")                  - *)
    (* ----------------------------------------------- *)
    lComm = {};     (* array for code of symbols *)
    For[ iI=1, iI<=Length[sCh], iI++,  
    (* For loop since we change the length *)
      If[ sCh[[iI]]//DigitQ ,
        (* we have a digit *)
        While[ iI+1<=Length[sCh] && DigitQ[sCh[[iI+1]]] ,
           sCh[[iI]] = sCh[[iI]] <> sCh[[iI+1]];
           sCh       = Drop[ sCh , {iI+1} ];
        ];
        AppendTo[ lComm , 4 ];   (* we have a number *)
      ,
        (* we have something else *)
        AppendTo[ lComm ,
          Switch[ sCh[[iI]],
               "(" , 2  ,          (* open bracket   *)
               ")" , -2 ,          (* closed bracket *)
               "+" , 5  ,          (* letter "+"     *)
               " " , 6  ,          (* letter " "     *)
               "," , 7  ,          (* letter ","     *)
               "-" , 8  ,          (* letter "-"     *)
               _   , 0             (* something else *)
          ]
        ];
      ];
    ];
    (* ---------------------------------------------- *)

    (* ---------------------------------------------- *)
    (* - Check for unrecognized symbols and         - *)
    (* - missing brackets                           - *)
    (* ---------------------------------------------- *)
    If[ Position[ lComm, 0 ] =!= {} , Return[False] ];  (* other symbols? *)
    iLevel = 0;               
    Do[          (* missing brackets? *)
      If[ Abs[lComm[[iI]]]===2 , iLevel += lComm[[iI]] ];
    ,{iI,Length[sCh]}];
    If[ iLevel =!= 0 , Return[ False ] ];   (* error, return False *)
    (* ---------------------------------------------- *)

    (* ---------------------------------------------- *)
    (* - lets find "-", output: nothing             - *)
    (* - it just rearanges the arrays...            - *)
    (* ---------------------------------------------- *)
    For[ iI=Length[sCh]-1, iI>0, iI--,  
    (* For loop since we change the length *)
      If[ lComm[[iI]]===8 ,
         If[ lComm[[iI+1]] === 4 ,
            lComm[[iI]] = -4;   (* single block, i.e. -number *)
         ,
            lComm[[iI]] = 2;    (* bracket *)
         ];
         sCh[[iI]]   = "-" <> sCh[[iI+1]] ;
         sCh   = Drop[ sCh   , {iI+1} ];
         lComm = Drop[ lComm , {iI+1} ];
      ];
    ];
    (* ---------------------------------------------- *)
    If[ (Print /. {opts})===On,
       Print["blocks: ",Thread[{sCh,lComm}]//MatrixForm];
    ];
(* Other way of printing:
Print[ Thread[{Range[sCh//Length],sCh,lComm}]//Transpose//MatrixForm ];
*)

    (* ---------------------------------------------- *)
    (* - operation " n"  ("n" is number!)           - *)
    (* ---------------------------------------------- *)
    For[ iI=2, iI<Length[sCh], iI++,       (* For loop since we change the \
length *)
      If[ lComm[[iI]]===6 && Abs[ lComm[[iI+1]] ]===4 ,     (* " num" *)
        (* what is the prevoius element? *)
        iPrev = iI-1;               (* this assumes single block *)
        If[ lComm[[iI-1]]===-2 ,    (* or it is a bracket        *)
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        sCh   =  Insert[ sCh  , "("       , iPrev ];
        lComm =  Insert[ lComm, 2         , iPrev ];
        iNext = iI+1;                                
        (* we inserted one element  *)
        (* now we search for other elements *)
        If[ iNext<=Length[sCh] && lComm[[iNext]]===6,     
        (* bilo While, sad je If *)
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 7  ; (* this is ordinary comma now *)
          iNext++;              (* next element          *)
          If[ lComm[[iNext]]===2 ,    (* is it a bracket? *)
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ; (* go to the next element, if any *)
	  ,  (* no, it is a number *)
            iNumRead = ToExpression[ sCh[[iNext]] ];
	    If[ iNumRead < 0 ,
	       iNumRead = -iNumRead;
               sOne     = "-1";
            ,
               sOne     = "1";
            ];
	    sCh[[ iNext ]] = sOne; 
	    iNext++;
	    Do[
              sCh   =  Insert[ sCh  , "," , iNext ];
              lComm =  Insert[ lComm, 7   , iNext ]; (* comma *)
	      iNext++;
              sCh   =  Insert[ sCh  , sOne , iNext ];
              lComm =  Insert[ lComm, 4*ToExpression[sOne] , iNext ]; 
              (* number    *)
	      iNext++;
            ,{iK,iNumRead-1}];
          ];
        ];
        sCh   =  Insert[ sCh  , ")" , iNext    ];
        lComm =  Insert[ lComm, -2   , iNext    ]; 
        (* we consider this as a bracket *)
      ];
    ];
    (* ---------------------------------------------- *)
    If[ (Print /. {opts})===On,
      Print["operation ' ': ",Thread[{sCh,lComm}]//MatrixForm];
    ];

    (* ---------------------------------------------- *)
    (* - operation "+n"  ("n" is number!)           - *)
    (* ---------------------------------------------- *)
    For[ iI=2, iI<Length[sCh], iI++,       (* For loop since we change the \
length *)
      If[ lComm[[iI]]===5 && Abs[ lComm[[iI+1]] ]===4 ,     (* "+num" *)
        (* what is the prevoius element? *)
        iPrev = iI-1;               (* this assumes single block *)
        If[ lComm[[iI-1]]===-2 ,    (* or it is a bracket        *)
            iLevel = 0;
            iJ     = iI-2;
            While[ iLevel<2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ--;
            ];
            iPrev = iJ+1 ;
        ];
        iNext = iI;       
        (* now we search for other elements *)
        If[ iNext<=Length[sCh] && lComm[[iNext]]===5,     
        (* bilo While, sad je If *)
          sCh[[iNext]]   = ",";
          lComm[[iNext]] = 7  ; (* this is comma now *)
          iNext++;              (* next element          *)
          If[ lComm[[iNext]]===2 ,    (* is it a bracket? *)
            iLevel = 0;
            iJ     = iNext+1;
            While[ iLevel>-2 ,
              If[ Abs[lComm[[iJ]]]===2 , iLevel += lComm[[iJ]] ];
              iJ++;
            ];
            iNext = iJ;         (* go to the next element, if any *)
	  ,  (* no, it is a number *)
            iNumRead = ToExpression[ sCh[[iNext]] ];
	    If[ iNumRead < 0 ,
	       iNumRead = -iNumRead;
               sOne     = "-1";
            ,
               sOne     = "1";
            ];
	    sCh[[ iNext ]] = sOne; 
	    iNext++;
	    Do[
              sCh   =  Insert[ sCh  , "," , iNext ];
              lComm =  Insert[ lComm, 7   , iNext ]; (* comma *)
	      iNext++;
              sCh   =  Insert[ sCh  , sOne , iNext ];
              lComm =  Insert[ lComm, 4*ToExpression[sOne] , iNext ]; (* \
number    *)
	      iNext++;
            ,{iK,iNumRead-1}];
          ];
        ];
      ];
    ];
    (* ---------------------------------------------- *)
    If[ (Print /. {opts})===On,
      Print["operation '+': ",Thread[{sCh,lComm}]//MatrixForm];
    ];
    
    lAllCh   = { sCh   } ;   (* mupltiple outputs *)
    lAllComm = { lComm } ;
    (* ---------------------------------------------- *)
    (* - operation ",1"    --> find multiples      -- *)
    (* ---------------------------------------------- *)
      For[ iK=1, iK<=Length[lAllCh], iK++ ,      
      (* loop over all multiples *)
      For[ iI=2, iI<Length[lAllCh[[iK]]] , iI++ ,    (* characters in one \
realization *)
        sCh   = lAllCh[[iK]];
	lComm = lAllComm[[iK]];      (* new realization to be generated *)
        If[ lComm[[iI]]===7 && (sCh[[iI+1]]==="1" || sCh[[iI+1]]==="-1"), 
        (* ",1" *)

          (* what is the prevoius element: a. begginig, b. "," or c. "(" *)
	  For[ iPrev=iI-1, iPrev>0 && lComm[[iPrev]]=!=7 && lComm[[iPrev]]=!=2, \
iPrev--,
             If[ lComm[[iPrev]]===-2 ,    (* or it is a bracket   ")"  *)
               iLevel = 0;
               While[ iLevel<2 ,
                 iPrev--;
                 If[ Abs[lComm[[iPrev]]]===2 , iLevel += lComm[[iPrev]] ];
               ];
             ];
          ];
          iPrev++;
      
          If[ !(iPrev===iI-1 && sCh[[iPrev]]===sCh[[iI+1]]) ,  
          (* no sence to swap "1,1" *)
            (* and swap them *)
	    tmp  = sCh[[iI+1]];
	    ltmp = lComm[[iI+1]];
	    Do[
                sCh[[ iJ ]]   = sCh[[ iJ-2 ]];
                lComm[[ iJ ]] = lComm[[ iJ-2 ]];
            ,{iJ,iI+1,iPrev+2,-1}];
            sCh[[iPrev]]     = tmp;
	    lComm[[iPrev]]   = ltmp;
            sCh[[iPrev+1]]   = ",";
	    lComm[[iPrev+1]] = 7;
 
	    (* do we have it alreadey in the list? *)
	    If[ !(MemberQ[lAllCh,sCh] && MemberQ[lAllComm,lComm]) ,	
	       AppendTo[ lAllCh, sCh ];
	       AppendTo[ lAllComm, lComm ];
            ];	  
          ];
        ];
      ];
    ];
    (* ---------------------------------------------- *)
    (* - operation "1,"    --> find multiples      -- *)
    (* ---------------------------------------------- *)
    For[ iK=1, iK<=Length[lAllCh], iK++ ,    (* loop over all multiples *)
      For[ iI=2, iI<Length[lAllCh[[iK]]], iI++ , (* characters in one \
realization *)
        sCh   = lAllCh[[iK]];
	lComm = lAllComm[[iK]];      (* new realization to be generated *)
        If[ lComm[[iI]]===7 && (sCh[[iI-1]]==="1" || 
        sCh[[iI-1]]==="-1"),     
 (* "1," *)

          (* what is the next element: a. end, b. "," or c. ")" *)
	  For[ iNext=iI+1, 
                 iNext<=Length[sCh] && lComm[[iNext]]=!=7 && \
lComm[[iNext]]=!=-2
          , iNext++,
             If[ lComm[[iNext]]===2 ,    (* or it is a bracket   ")"  *)
               iLevel = 0;
               While[ iLevel>-2 ,
                 iNext++;
                 If[ Abs[lComm[[iNext]]]===2 , iLevel += lComm[[iNext]] ];
               ];
             ];
          ];
          iNext--;
      
          If[ !(iNext===iI+1 && sCh[[iNext]]===sCh[[iI-1]]) ,  
          (* no sence to swap "1,1" *)
            (* and swap them *)
	    tmp  = sCh[[iI-1]];
	    ltmp = lComm[[iI-1]];
	    Do[
                sCh[[ iJ ]]   = sCh[[ iJ+2 ]];
                lComm[[ iJ ]] = lComm[[ iJ+2 ]];
            ,{iJ,iI-1,iNext-2}];
            sCh[[iNext]]     = tmp;
	    lComm[[iNext]]   = ltmp;
            sCh[[iNext-1]]   = ",";
	    lComm[[iNext-1]] = 7;

	    (* do we have it alreadey in the list? *)
	    If[ !(MemberQ[lAllCh,sCh] && MemberQ[lAllComm,lComm]) ,	
	       AppendTo[ lAllCh, sCh ];
	       AppendTo[ lAllComm, lComm ];
            ];	  
          ];
        ];
      ];
    ];
    (* ---------------------------------------------- *)

    sOut =  Map[ StringJoin , lAllCh ];   (* Join the strings! *)
    If[ (Print /. {opts})===On,
       Print["Input: \"" <> sInput <> "\"" ];
       Do[
             Print[ "  " <> ToString[iK] <> " --> \"" 
             <> sOut[[iK]] <> "\"" \

];
       ,{iK,lAllCh//Length}];
    ];
    If[ Head[ tmp=(Print /. {opts}) ] === OutputStream ,
       WriteString[tmp,"Input: \"" <> sInput <> "\"\n" ];
       Do[
         WriteString[tmp,"  " <> ToString[iK] <> " --> \"" 
         <> sOut[[iK]] <> "\
\"\n" ];
       ,{iK,lAllCh//Length}];
    ];

    sOut    (* Output *)
  ];
(* ================================================ *)


fPr[Con_String, k_Integer]:=Module[{pp},
    pp=fProjections[Con][[k]];
    pp
    ]
    
fProject[Con_String, k_Integer]:=Module[{pp},
    pp=fDiffProjectionsAltKL[Con][[k,1]];
    pp
    ]


(**************************************************)
(*## ## ## ## ## ## # kraj fPROJECTIONS ## ## ## ## ## ## ## ##*)


(*## ## ## ## ## ## # fPROJECTIONS ADD ## ## ## ## ## ## ## ##*)

fPartL[LL_List, LL1_List] := Module[{ss = LL, pp = LL1, i},
    For[i = 1, i < Length[pp] + 1, ss = Part[ss, pp[[i]]]; i++];
    ss]

fPartList[LL_List, LL1_List, i_Integer] := Module[{ss = LL, pp = LL1},
    pp = Join[Drop[pp, -1], {i}];
    pp = fPartL[ss, pp];
    pp
    ]

fBlockOne[LL_List] := 
  Module[{dd2 = LL, mm, mm1, mm2, i}, 
    mm = Table[If[SameQ[dd2[[i]], 1], 1, 0], {i, Length[dd2]}];
    mm = Table[If[SameQ[dd2[[i]], 1], 1, 0], {i, Length[dd2]}];
    mm1 = Table[If[SameQ[Head[dd2[[i]]], List], 2, 0], {i, Length[dd2]}];
    mm2 = 
      Table[If[SameQ[NumberQ[dd2[[i]]], True] && dd2[[i]] >= 2, 3, 0], {i, 
          Length[dd2]}];
    mm = mm1 + mm2;
    mm = Split[mm];
    mm = Flatten[
        Table[If[MemberQ[mm[[i]], 2] && Length[mm[[i]]] > 1, 
            Partition[mm[[i]], 1], {mm[[i]]}], {i, Length[mm]}], 1];
    mm1 = 
      Flatten[Position[Table[If[MemberQ[mm[[i]], 2], 1, 0], 
      {i, Length[mm]}], 1]];
    mm2 = mm1 - 1;
    mm = Table[
        If[MemberQ[mm[[mm2[[i]]]], 0], Length[mm[[mm2[[i]]]]], 0], {i, 
          Length[mm2]}];
    mm
    ]


fChangePart[Ulaz_String] := Module[{ss, dep, dd, ff, dd1, dd2, dd3, k, i},
    ss = ToExpression[StringReplace[Ulaz, {"(" -> "{", ")" -> "}"}]];
    dep = Depth[ss];
    Do[
      dd = Depth[ss];
ff = dd;
      dd = Flatten[ss, dd - k];
      dd = Select[dd, SameQ[Head[#], List] &];
      dd1 = Flatten[Table[Position[ss, dd[[i]]], {i, Length[dd]}], 1];
       dd1 = Select[dd1, Length[#] == ff - k + 1 &];
      dd1 = Union[dd1]; 
      dd3 = fPartL[ss, Drop[dd1[[1]], -1]];
      (* ovo je deo na kome radimo *)
      ll = fBlockOne[dd3];
      dd2 = Table[fPartL[ss, dd1[[i]]], {i, Length[dd1]}];
      dd2 = 
        Table[If[OddQ[ll[[i]]], Reverse[dd2[[i]]], dd2[[i]]], {i, 
            Length[dd2]}];
      Do[ss = ReplacePart[ss, dd2[[i]], dd1[[i]]], {i, Length[dd1]}];
      ss, {k, 3, dep}];
    ss = ToString[ss];
    ss = StringReplace[ss, {"{" -> "(", "}" -> ")", " " -> ""}];
    ss
    ]
    
    fFixP[Ul_String] := Module[{pr = Ul, pp, i}, 
    pp = StringPosition[pr, " "];
    pr = StringJoin["{", 
        StringReplace[pr, {"(" -> "{", ")" -> "}", " " -> ","}], "}"];
    pr = ToExpression[pr];
    pr = Table[
        StringReplace[
          ToString[pr[[i]]], {" " -> "", "{" -> "(", "}" -> ")"}], {i, 
          Length[pr]}];
    pr = Map[fChangePart, pr];
    pr = Table[
        If[StringLength[pr[[i]]] == 3, StringDrop[StringDrop[pr[[i]], 1], \
-1],
           pr[[i]]], {i, Length[pr]}];
    pr = StringDrop[
        StringJoin[Table[StringJoin[pr[[i]], ","], {i, Length[pr]}]], -1];
    pr = StringReplacePart[pr, " ", pp];
    pr]
    

   (*## ## ## ## ## ## ## # fPPoly ## ## ## ## ## ## ## ## ## ## ## #*)

(*computes projections of polyhedral KL *)
(* 15. 05. 2005 *)
(*".2 2..2..2"="6*.2 2..2..2"and 
"2 2..2..2"="6*2 2..2..2"*)
(*computes projections of polyhedral KL *)

fPPoly[Conway_String]:=
  Module[{Con=Conway,pos1,pos2,res={},str,pom},
  If[SameQ[StringPosition[Con,"*"],{}]&&
  Not[SameQ[StringPosition[Con,"."],{}]]
&&Not[SameQ[StringPosition[Con,"."][[1]],{1,1}]],Con="6*"<>Con,Con];
    pos1=Flatten[StringPosition[Con,"*"]];
     If[pos1=={},Con=Con<>".",
      If[pos1[[1]]!=StringLength[Con],Con=Con<>"."]];
         pos2=Union[Flatten[StringPosition[Con,"."]]];
     If[pos1!={},str=StringTake[Con,pos1[[1]]];
      Con=StringDrop[Con,pos1[[1]]],str=StringTake[Con,pos2[[1]]];
      Con=StringDrop[Con,pos2[[1]]];
      pos2=Rest[pos2](*izbacili 1. tackicu*)];
        (*izdvojili smo pocetak i sad idemo od tackice do tackice tacnije 1. 
    Prvo je tackica:prepisemo je;
      2. prvo je broj:uzimamo deo do sledece tackice i radimo fProjections*)
      res={{str}};
    If[Con==".",Con=""];
    While[Con!="",
      If[StringTake[Con,1]==".",(*prva tackica*)Con=StringDrop[Con,1];
        res=Map[(#<>".")&,res],(*pom je parce do sledece tackice*)pom=
          StringTake[Con,StringPosition[Con,"."][[1,1]]-1];
        (*brisemo iz stringa parce*)Con=
          StringDrop[Con,StringPosition[Con,"."][[1,1]]];
        If[StringLength[pom]>2,
          If[SameQ[StringTake[pom,-2]," 0"],(*ako je nula na kraju*)
          pom=StringDrop[pom,-2];
            pom=fProjKL[pom];
            
            pom=Map[#<>" 0"&,pom],(*ako nema _ 0 i duzi od 3*)pom=
              fProjKL[pom]],(*ako nema _ 0 i kraci od 3*)pom=
            fProjKL[pom]];(*sve projekcije datog parceta*)
           pom = Flatten[Map[fFixP, pom]];
          Do[
          If[Con=="",str=Map[(res[[i]]<>#)&,pom],
            str=Map[(res[[i]]<>#<>".")&,pom]];
          res=ReplacePart[res,str,i],{i,1,Length[res]}];
        res=Flatten[res]];(*kraj ifa*)If[Con==".",Con=""]];
    res]
    
    
  
    (*## ## ## ## ## ## ## # fPROJECTIONS  ## ## ## ## ## ## ## ## ##*)
    
    fProjections[Con_String] := 
  Module[{res, s, s1, pr},    
   If[Not[SameQ[StringTake[Con,-1],"*"]],   
    If[Union[Flatten[StringPosition[Con, "*"]], 
          Flatten[StringPosition[Con, "."]]] == {}, res = fProjKL[Con], 
      res = fPPoly[Con]];
      If[Union[Flatten[StringPosition[res[[1]], "*"]], 
          Flatten[StringPosition[res[[1]], "."]]] == {},
    If[Length[res] == 1,
      pr = res,
      pr = Flatten[Map[fFixP, res]]];
    pr = Flatten[pr];
    pr, pr=res], pr={Con}];
    pr]
    
    
(* fProjections[Con_String] := 
  Module[{res, s, s1, pr}, 
    If[Union[Flatten[StringPosition[Con, "*"]], 
          Flatten[StringPosition[Con, "."]]] == {}, res = fProjKL[Con], 
      res = fPPoly[Con]];
      Print[res];
    If[Length[res] == 1,
      pr = res,
      pr = Flatten[Map[fFixP, res]]];
    pr = Flatten[pr];
    pr]*)
     

(*## ## ## ## ## ## ## ## F DIFFERENT PROJECTIONS ## ## ## ## ## ##*)
fDiffProjectionsAltKL[Con_String]:=Module[{ind,l,i,j,res={},pom},
    l=fProjections[Con];
    If[Con=="2",
        res={{"2",{{1, 1},{4, 2}}}},
         If[Con=="-2",
        res={"2",{{1, 1},{4, 2}}}, 
    If[Length[l]==1&&Head[l[[1]]]==List,l=l[[1]]];
    (*red samo za 6* i sl*)
    Do[ pom=Abs[MinDowProjAltKL[l[[i]]]];
              ind=False;j=1;
           While[  ind==False&&j<=Length[res],
                           ind=SameQ[pom,res[[j,2]]];j++];
       If[ind==False,res=Append[res,{l[[i]],pom}]]
      ,{i,1,Length[l]}]
    ]];res]
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #*)


(*## ## ## ## ## ## # fTorus KL ## ## ## ## ## ## ## ## ## ## ## ## ##*)

(*fTorusKL calculates for a torus KL[a,b] braid word,min.no.of
 crossings, unknotting no.,no.of components,brodge no.,
 Alexander polynomial and Murasugi signature*)
(**********)
fTorusKL[a_Integer,b_Integer]:=
  Module[{qq,rr,dd,www,bw,pd,pdr,cr,mcr,br,bwr,Alex},
    qq=Max[{a,b}];
    rr=Min[{a,b}];
    dd=GCD[qq,rr];
    www=FromCharacterCode[Range[rr-1]+96];
    bw=StringJoin[Table[www,{i,qq}]];
    pd=KnotFromBraid[bw];
    pd=ReductionKnotLink[pd];
    Print["Braid word: ",bw];
    If[dd>1,
        Print["Crossings: ",Length[pd[[2]]]];
        Print["Number of components: ",fComponentNo[pd]];
        (*Print["Unlinking number: ",UnKnotLink[pd]];*)
        Print["Alexander polynomial: ", RedAlex[pd]];
      	Print["Signature: ",fSignat[pd]],
      cr=Max[qq (rr-1),rr (qq-1)];
      mcr=Min[qq (rr-1),rr (qq-1)];
      Print["Crossings: ",mcr];
      Print["Number of components: ",fComponentNo[pd]];
      If[SameQ[dd,1],Print["Unknotting number: ",(qq-1)(rr-1)/2],
        Print[dd," -component link"]];
      br=Min[qq,rr];
      Print["Bridge number: ",br];
      Alex=PolynomialQuotient[(1-x)((1-x^(qq*rr)/dd)^dd),(1-x^qq)(1-x^rr),x];
      Print["Alexander polynomial: ",Alex];
      Print["Signature: ",fSignat[pd]];
      Print["Murasugi signature: ",
        If[SameQ[a,b],If[SameQ[Mod[a,2],1],(a^2-1)/2,(a^2-2)/2],
          fSigTor[a,b]]]
      ];
      Print["P-data: ",pd];
     ]
      

  ShowTorusKL[a_Integer,b_Integer]:=
  Module[{qq,rr,dd,bw,pd,cr,mcr,br,bwr,Alex,www},
    qq=Max[{a,b}];
    rr=Min[{a,b}];
    dd=GCD[qq,rr];
    www=FromCharacterCode[Range[rr-1]+96];
    bw=StringJoin[Table[www,{i,qq}]];
    pd=KnotFromBraid[bw];
    pd=ReductionKnotLink[pd];
   If[dd>1,
      cr=Max[qq (rr-1),rr (qq-1)];
      mcr=Min[qq (rr-1),rr (qq-1)]];     
      Graphics3DKnotfromPdata[pd]]
    


 ShowTorusKLNew[a_Integer,b_Integer]:=
  Module[{qq,rr,dd,bw,pd,cr,mcr,br,bwr,Alex,www},
    qq=Max[{a,b}];
    rr=Min[{a,b}];
    dd=GCD[qq,rr];
    www=FromCharacterCode[Range[rr-1]+96];
    bw=StringJoin[Table[www,{i,qq}]];
    pd=KnotFromBraid[bw];
    pd=ReductionKnotLink[pd];
   If[dd>1,
      cr=Max[qq (rr-1),rr (qq-1)];
      mcr=Min[qq (rr-1),rr (qq-1)]];     
      ShowKnotfromPdataNew[pd]]  

      
(*## ## ## ## ## ## ## # fAlexPoly ## ## ## ## ## ## ## ## ## ##*)

(*calculates multivariable Alexander polynomial *)

(*## ## ## ## ## ## ## # fALEX ## ## ## ## ## ## ## ## ## ##*)(*calculates \
multivariable \
Alexander polynomial*)
(*April 3, 2005 *)
fAlexPoly[Conway_String]:=
  Module[{lL0,res,lL,t,i,s,sl,j,p,MAT},
    If[Conway=="1",res=0,
      If[Conway=="2",res=1,lL0=fGenerators[Conway][[1]];
        lL=Flatten[lL0];
       t=Table[0,{i,Length[lL0]}];
 Do[
          If[i==1,t[[i]]=Length[lL0[[i]]]/3,
            t[[i]]=t[[i-1]]+Length[lL0[[i]]]/3],{i,1,Length[lL0]}];
        s=fGenerators[Conway][[2]];
        sl=Length[s];(*duzina*)
 MAT=ConstantArray[0,{sl,sl}];
        j=1;p="a";
        For[i=1,i<=sl,i++,If[i>t[[j]],j++;
            p=FromCharacterCode[ToCharacterCode["a"]+j-1]];
            Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+1]]]=1,-1,
            MAT[[i,Part[lL,3(i-1)+1]]]=-p];
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+2]]]=-p,-1,
            MAT[[i,Part[lL,3(i-1)+2]]]=1];
          MAT[[i,Part[lL,3i]]]=p-1;];
        (*   Print[MAT]; *)
          MAT=
          Apply[PolynomialGCD,
            Table[First[Part[Minors[MAT],i]],{i,1,Length[MAT]}]];
       If[SameQ[MAT,0],res=0,
        MAT=If[SameQ[MAT[[1]],-1],-MAT,MAT];
        res=Expand[Divide[MAT,PolynomialGCD[MAT[[1]],MAT]]]]] ];
        res]
        
        
  fAlexPolyOne[Conway_String]:=
  Module[{lL0,res,lL,t,i,s,sl,j,p,MAT,dd,dd1},
    If[Conway=="1",res=0,
      If[Conway=="2",res=1,lL0=fGenerators[Conway][[1]];
        lL=Flatten[lL0];
         t=ConstantArray[0,{1,Length[lL0]}];
       (* t=ZeroMatrix[1,Length[lL0]][[1]]; *)
        Do[
          If[i==1,t[[i]]=Length[lL0[[i]]]/3,
            t[[i]]=t[[i-1]]+Length[lL0[[i]]]/3],{i,1,Length[lL0]}];
        s=fGenerators[Conway][[2]];
        sl=Length[s];(*duzina*) (* MAT=ZeroMatrix[sl]; *)
         MAT=ConstantArray[0,{sl,sl}];
        j=1;p="a";
        For[i=1,i<=sl,i++,If[i>t[[j]],j++;
            p=FromCharacterCode[ToCharacterCode["a"]+j-1]]; 
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+1]]]=1,-1,
            MAT[[i,Part[lL,3(i-1)+1]]]=-p];
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+2]]]=-p,-1,
            MAT[[i,Part[lL,3(i-1)+2]]]=1];
          MAT[[i,Part[lL,3i]]]=p-1;];
        dd=fComponentNo[Conway];
        dd1=Table[FromCharacterCode[96+i],{i,dd}];
        MAT=
          First[Flatten[
              Minors[ReplaceAll[MAT,Table[dd1[[i]]->a,{i,dd}]]]]];
        res=If[SameQ[MAT,0],MAT,fPolyNorm[MAT]];
        (* If[SameQ[MAT,0],res=0,MAT=If[SameQ[MAT[[1]],-1],-MAT,MAT];
            res=Expand[Divide[MAT,PolynomialGCD[MAT[[1]],MAT]]]] *)] 
      ];
    res]
        
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ##*)
 
 (*## GAUSS CODE ## ##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
##*)
(*## ## fGaussExtSigns ## ## ## ## ## ## ##*)

(*# Input List-
      Dowker Code {{Lengths of Components},{code}} 
      ConwaySymbol_ String 
      PData
Output Extended Dowker=
    Gauss Code with signs #*)
    (*Neded for MinDowPr for Knots*)
    
  fGaussExtSigns[Ulaz_]:=
  Module[{DowL,SL,LL,DowPair={},DowPair1,DowFl,i,k=0,DL,
      pomL={},Ul=Ulaz},
       (*Input:Pdata*)
      If[SameQ[Head[Ul],List],
      If[MemberQ[Map[OddQ,Abs[Ul[[2]] ] ] ,True],
      Ul=fDowfromPD[Ul]]];
        (*sad nam je Ul iliKonvej ili dowker *)
      (*Input:Dowker*)
      If[SameQ[Head[Ul],List],DowL=Abs[Ul[[2]]];
      If[Length[Ul[[1]]]==1,LL={},LL=2*Ul[[1]]];
      SL=Flatten[Map[Sign[#]&,Ul[[2]]]]];
    (*Input:Conway*)
    If[Head[Ul]==String,DowL=Dowker[Ul];
      If[DowL!={},DowL=DowL[[4]];
        SL=fGenSign[Ul];
        If[SameQ[Flatten[DowL],DowL],LL={},LL=Map[2*Length[#]&,DowL]]]];
    (*Print[DowL,SL,LL];*)(*Form ext*)DowFl=Flatten[DowL];
    Table[DowPair=Append[DowPair,{2i-1,DowFl[[i]]}];
      DowPair=Append[DowPair,{DowFl[[i]],2i-1}],{i,Length[DowFl]}];
    DowPair=Union[Map[Sort[#]&,Sort[DowPair]]];
    (*Print["Ovde je",DowPair];*)DowPair1=
      Map[Append[#,
            SL[[Position[DowFl,Select[#,EvenQ][[1]]][[1,1]]]]*
              Position[DowPair,#][[1,1]]]&,DowPair];
    (*Print["Ili ovde ",DowPair1];*)DL=
      Union[Map[Take[#,-2]&,DowPair1],Map[{First[#],Last[#]}&,DowPair1]];
    DL=Flatten[
        Map[Take[#,-1]&,DL]];(*Print["Prvi prosireni ",
          DL];*)(*ako je link moramo da ga podelimo*)If[LL!={},
      Do[If[i==1,k=LL[[1]];pomL=Append[pomL,Take[DL,{1,LL[[1]]}]],
          pomL=Append[pomL,Take[DL,{k+1,k+LL[[i]]}]];
          k=k+LL[[i]]],{i,1,Length[LL]}];(*Print["DL ",DL];*)DL=
        pomL];(*Print["DP ",DowPair];*)DL]

(*## ## fGaussExt ## ## ## ## ## ## ##*)
(*# Input List-
      Dowker Code {{Lengths of Components},{code}} ConwaySymbol_ String \
Output Extended Dowker=Gauss Code without signs #*)
fGaussExt[Ul_]:=Module[{r},r=Abs[fGaussExtSigns[Ul]]]

(*## ## ## ## ## ## Dowker from ExtGauss ## ## ## ## ## ## ## ## #*)
fDowfromGaussExt[Ul_List]:=Module[{LL,ExtD={Ul},pom,k=0,i,Dow={}},
    If[SameQ[Length[Ul],Length[Flatten[Ul]]],
      LL={};
      ExtD=ExtD[[1]],
      ExtD=Flatten[ExtD,1];
      LL=Map[(Length[#]/2)&,ExtD];
      ExtD=Flatten[fFixLink[ExtD]]];
     pom=Union[Map[Append[Flatten[Position[ExtD,#]],Sign[#]]&,ExtD]];
    pom=Sort[Map[If[EvenQ[#[[1]]],{#[[2]],#[[1]],#[[3]]},#]&,pom]];
    pom=Map[(#[[2]]*#[[3]])&,pom];
    If[LL=={},LL={Length[pom]};Dow=pom,
      Do[If[i==1,Dow=Append[Dow,Take[pom,{1,LL[[1]]}]],
          Dow=Append[Dow,Take[pom,k+{1,LL[[i]]}]]];
        k=k+LL[[i]],{i,1,Length[LL]}]];
    {LL,Dow} ]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
##*)

(*## MINIMIZATION ## ##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
##*)
(*# fMinComp ## ## ## ## ## #*)(*# Minimizes ExtDow of a knot #*)
(*# pomocna za fMinDowPr #*)
(*# Could be used for minimizing link components #*)

fMinComp[ExtDowComp_List]:=Module[{d1,d2,m,p,pL,i,pom, t},
    d1=Map[Flatten[Position[ExtDowComp,#]]&,ExtDowComp];
    d1=Union[
        Map[If[Length[#]==1,{#[[1]],
                Length[ExtDowComp]},{#[[1]],#[[2]]-#[[1]]}]&,d1],
        Map[If[Length[#]==1,{#[[1]],Length[ExtDowComp]},{#[[1]],
                Length[d1]-#[[1]]+#[[2]]}]&,Map[Reverse[#]&,d1]]];
    d1=Flatten[Map[Take[#,-1]&,d1]];
    d2=ReplaceAll[Length[d1]-d1,0->Length[d1]];
    m=Min[Join[d1,d2]];
    (*m je minimal distance in both lists*)
    p={Flatten[Position[d1,m]],Flatten[Position[d2,m]]};
   pL=
      Join[Map[{Take[RotateLeft[d1,#-1],{1,Length[d1],2}],#,0}&,p[[1]]],
        Map[{Take[Reverse[RotateRight[d2,Length[d2]-#]],{1,Length[d1],2}],#,
              1}&,p[[2]]]];
    (*Element pl je oblika:odgovarajuca lista razliak, 
      pos od koje smo krenuli i 0 ili 1 za smer*)
      pL=Sort[pL];
    pom=Flatten[Union[Map[Take[#,1]&,pL]],1];
    (*Biramo one koje imaju razlicite razlike bez 
    obzira na smer i poziciju*) 
    m={};
    Do[m=Append[m,Select[pL,SameQ[#[[1]],pom[[i]]]&][[1]]],{i,1,
        Length[pom]}];
    m;
    (* naredni deo dopisan -uzima samo minimalne *)
    t=Select[
        Sort[Table[{Take[Sort[m][[i]][[1]],2],Sort[m][[i,2]],
              Sort[m][[i,3]]},{i,Length[Sort[m]]}]],#[[1]]==
            Sort[Table[{Take[Sort[m][[i]][[1]],2],Sort[m][[i,2]],
                    Sort[m][[i,3]]},{i,Length[Sort[m]]}]][[1,1]]&];
       m=Table[
       If[SameQ[Take[m[[i]][[1]],2],t[[i]][[1]]],m[[i]],{}],{i,Length[t]}];
      m
    ]
    
(*## ## fKnittDow ## ## ## ## ## ## ##*)
(*# pomocna za fMinDowPr #*)
fKnittDow[ExtDow_List,UlList_List]:=Module[{pom,Ps},
    If[Head[UlList[[1]]]==Integer,
      If[UlList[[2]]==0,pom=RotateLeft[ExtDow,UlList[[1]]-1],
        pom=Reverse[RotateRight[ExtDow,Length[ExtDow]-UlList[[1]]]]];
      Ps=Map[Sign[#]&,Sort[Union[pom],Abs[#1]<Abs[#2]&]]];
    {pom,Ps}]

(*## ## fFormDow ZAMENJENA ## ## ## ## ## ## ##*)

fFormDow[KKK_List,Sig_List]:=Module[{LL=KKK,q, i,t,t1},
    LL=Union[Map[Flatten[Position[Abs[KKK],Abs[#]]]&,KKK]];
    t=Table[LL[[i]][[1]],{i,Length[KKK]/2}];
    t1=Table[
        Sign[KKK[[t[[i]]]]],{i,Length[Table[LL[[i]][[1]],{i,Length[LL]}]]}];
    q=Table[{LL[[i]],Sign[KKK[[t[[i]]]]]},{i,Length[t]}];
    q=Table[
        Sort[Table[
                If[EvenQ[q[[i,1,1]]],{q[[i,1,2]],
                    q[[i,1,1]]*q[[i]][[2]]},{q[[i,1,1]],
                    q[[i,1,2]]*q[[i]][[2]]}],{i,Length[LL]}]][[i]][[2]],{i,
          Length[LL]}];
    q
    ]
(*## ## ## ## ## ## ## ## ## ## ## ? #### ## ## ## ## ## ## ## #*)

(*## INCIDENCY GRAPH ## ##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ## # fFormGr ## ## ## #*)
(*# radi za cvorove ili komponente linkova #*)
fFormGr[Ul_List]:=Module[{i,rez={}, j},
    If[SameQ[Length[Ul],Length[Flatten[Ul]]],Do[If[i==Length[Ul],
          rez=Append[rez,{Last[Ul],First[Ul]}],
          rez=Append[rez,{Ul[[i]],Ul[[i+1]]}];],{i,1,Length[Ul]}], 
      rez=Flatten[
          Table[Append[
              Table[{Ul[[j]][[i]],Ul[[j]][[i+1]]},{i,
                  Length[Ul[[j]]]-1}],{Last[Ul[[j]]],First[Ul[[j]]]}],{j,
              Length[Ul]}],1]];
    rez=Map[If[Abs[#[[1]]]>Abs[#[[2]]],{#[[2]],#[[1]]},#]&,rez]; 
    rez]
(*## ## # fMakeS ## ## ## #*)
fMakeS[LL_List,SL_List]:=Module[{r},r=Map[(#*SL[[#]])&,LL];
    r]
(*## ## # fSignSum ## ## ## #*)
fSignSum[LL_List]:=Module[{r,i},r=Sum[Sign[LL[[i]]],{i,1,Length[LL]}];
    r]
(*## ## # fGrInc ## ## ## #*)
(*# radi graf incidencije od prosirenog dowkera #*)
(*# Input ExtDow Output List of UnOrderedPairs (Graph) #*)
(*# Needs:#*)
fGrInc[Ul_List]:=Module[{rez={},sl},
    If[SameQ[Length[Ul],Length[Flatten[Ul]]],
      rez=fFormGr[Ul],
      rez=Flatten[Map[fFormGr[#]&,Ul],1]];
      rez=Sort[rez,Abs[#1]<Abs[#2]&];
      sl=Sort[Union[Flatten[rez]],Abs[#1]<Abs[#2]&];
     sl=Map[Sign[#]&,sl];
    {Sort[Abs[rez]],sl}]
(*## ## ## ## ## ## fGraphInc ## ## ## ## ## ## ## ## ## ## ## ## ## #*)
(*# isto sto i fGrinc, samo direktno iz 
Conwaya ili Dowkera ili pdata #*)
fGraphInc[Ulaz_]:=Module[{Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
fGrInc[fGaussExtSigns[Ul]]]
(*## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ORIENTED LINKS ## ##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
#*)

(*## ## ## ## # OldLinkingNo ## ## ## ## ## ## ## ## ##*)
OldLinkingNo[Ul_]:=
  Module[{GE,LinNoList={},d,i,j,k,p,rez={}},
    If[Not[SameQ[Head[Ul],String]],
      If[Length[Ul]!=1,
        If[SameQ[Select[Ul[[2]],OddQ[#1]&],{}],
        GE=fGaussExtSigns[Ul],GE=Ul],GE=Ul],
      GE=fGaussExtSigns[Ul]];
    d=Length[GE];
    If[Length[GE]==1,Print["Knot"],
      Do[Do[p=Map[Sign[#]&,Intersection[GE[[i]],GE[[j]]]];
          
          LinNoList=Append[LinNoList,Sum[p[[k]],{k,1,Length[p]}]],{j,i+1,
            d}],{i,1,d-1}];
      rez=
        Append[rez,{Sum[LinNoList[[k]],{k,1,Length[LinNoList]}],
            Sum[Abs[LinNoList[[k]]],{k,1,Length[LinNoList]}]}]
      (*Print["Linking number of oriented link is:  ",
          Sum[LinNoList[[k]],{k,1,Length[LinNoList]}]];
        Print["Linking number of nonoriented link is:  ",
          Sum[Abs[LinNoList[[k]]],{k,1,Length[LinNoList]}]]*)];
    Flatten[rez,1]/2]
(*radi od prosirenig Dowkera*)

(*## ## ## # Pomocne za fOrientedLink ## ## ## ## ## ## ## #*)
(*## ## # fFixLink ## ## ## #*)
(*# Rotation of Components due to oddity #*)
fFixLink[Ul_List]:=
  Module[{FlD=Flatten[Ul],p,rez={Ul[[1]]},i},
    Do[p=Flatten[Position[FlD,Ul[[i,1]]]];
      If[EvenQ[p[[2]]-p[[1]]],rez=Append[rez,RotateLeft[Ul[[i]]]],
        rez=Append[rez,Ul[[i]]]],{i,2,Length[Ul]}];
    rez]
(*## ## # fVarP ## ## ## #*)
fVarP[n_Integer]:=Module[{r,i},r=Table[IntegerDigits[i,2,n],{i,0,2^n-1}];
    r]
    
     (*## ## # fVarNewP ## ## ## #*)
    
    fVarNewP[n_Integer]:=Module[{r,i},
    r=Drop[Table[IntegerDigits[i,2,n],{i,0,2^n-1}],1];
    r=Sort[Table[{Count[r[[i]],1],r[[i]]},{i,Length[r]}]];
    r=Drop[
        Union[Table[
            If[r[[i,1]]>IntegerPart[(n+1)/2],{},r[[i,2]]],{i,Length[r]}]],1];
    r]
(*## ## # fChangeCommSign ## ## ## #*)
(*# Menja znake onih elemenata ostalih 
komponenti koji pripadaju komponenti sa kojom radimo #*)
fChangeCommSign[LL_List,p_List]:=
  Module[{i,rez={},q},Do[q=Map[If[MemberQ[Abs[p],Abs[#]],-#,#]&,LL[[i]]];
      rez=Append[rez,q],{i,1,Length[LL]}];
    rez]
    
   
(*## ## # fChangeOne ## ## ## #*)
(*Menja znake u komponenti koju trenutno radimo samo onima koje nisu \
samopreseci*)
fChangeOne[DL_List,CL_List,b_Integer]:=
  Module[{p,q},p=Flatten[Complement[DL,{DL[[b]]}]];
    p=Intersection[Abs[p],Abs[CL]];
    q=Map[If[MemberQ[Abs[p],Abs[#]],-#,#]&,CL];
    (*q je komponenta*)
    {Reverse[q],p}]
(*vraca datu komponentu kako treba i p njene zajednicke sa ostalima*)


(*## ## # fMakeOL ## ## ## #*)
fMakeOL[DL_List,Ind_List]:=
  Module[{rez={},i,pom=DL},
    Do[If[Ind[[i]]==0,rez=Append[rez,First[pom]];
        pom=Drop[pom,1],p=fChangeOne[DL,First[pom],i];
        (*Print["el:",p];*)rez=fChangeCommSign[rez,p[[2]]];
        (*Print[rez];*)rez=Append[rez,p[[1]]];
        (*Print[rez];*)pom=fChangeCommSign[Drop[pom,1],p[[2]]];
        (*Print[pom]*)],{i,1,Length[DL]}];
    rez]

(*## ## ## # fOrientedLink ## ## ## ## ## ## ## ## ##*)

fOrientedLink[Ul_]:=
  Module[{GE,VarL,EL,ExtList,rez={},i},GE=fGaussExtSigns[Ul];
    (*GE je prosireni Gauss dowker*)
    If[Not[SameQ[GE,Flatten[GE]]],
      VarL=fVarP[Length[GE]];
      ExtList=Map[fMakeOL[GE,#]&,VarL];
      ExtList=Map[fFixLink[#]&,ExtList];
      EL=Map[OldLinkingNo[#][[1]]&,ExtList];
      Do[rez=Append[rez,{ExtList[[i]],VarL[[i]],EL[[i]]}],    
           {i,1,Length[ExtList]}],
      Print["Error: Knot"]];
    rez]
(*## ## ## ## ## ## ## ## ##   PLANARNI ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
##*)

(* fPlanarEmb[neu_List]: 
    ulaz-Conway ili Dow izlaz- ulazni graf, graf izomorfan sa ulaznim,  
   planar embedding izomorfnog grafa, izomorfizam i faces izomorfnog grafa *)
   
fFindOrientation[L_List,VerLabel_List]:=Module[{p,i},
    p=Drop[TranslateVertices[L,-L[[1]]],1];
    p=Map[(#[[1]]+I*#[[2]])&,p];
    p=Map[Arg[#]&,p]+Pi;
    Do[p=ReplacePart[p,{p[[i]],VerLabel[[i+1]]},i],{i,1,Length[p]}];
    p=Prepend[Map[#[[2]]&,Reverse[Sort[p]]],First[VerLabel]];
    p]
    
fFind4[Graf_List]:=Module[{l},
    l=Union[Flatten[Graf]];
     l=Union[Map[Count[Flatten[ Graf],#]&,l]];
    Select[l,#!=4&]=={}
    ]
(* Vraca False ako nisu sva temena cetvorovalentna*)
(*Znaci True ako je graf \
u redu*)


fOrientGr[GG_Graph,neUredjeni_List]:=
  Module[{l,g,rez={},CRD,sus,j,i,ind=True,GNew},
    l=Flatten[FromAdjacencyMatrix[GG[[1]]][[1]],1];
    Do[sus=Select[l,MemberQ[#,i]&];
          sus=Map[If[SameQ[#[[1]],i],#[[2]],#[[1]]]&,sus];
         CRD=Map[Extract[GG[[2]],#]&,sus];
          
      rez=Append[rez,
          fFindOrientation[Prepend[CRD,GG[[2,i]]],Prepend[sus,i]]],
      {i,1,Length[GG[[1]]]}];
    g=Flatten[
        Table[Table[{rez[[j]][[1]],rez[[j]][[i]]},{i,2,Length[rez[[j]]]}],{j,
            Length[rez]}],1];
    g=Union[Map[Sort,g]];
    g=Isomorphism[FromUnorderedPairs[Union[neUredjeni]],
    FromUnorderedPairs[g],All];
    j=1;
    While[j<=Length[g]&&ind,                         
      GNew=Sort[
          Map[Sort,ReplaceAll[neu,Table[i->g[[j,i]],{i,Length[g[[j]]]}]]]];
           If[fFind4[GNew],ind=False,j++]
      ];
    {GNew,g[[j]],rez}
    ]
    
    
(*GG je novi graf,
  g[[j]] odgovarajuci izomorfizam a rez je planar embeding dat ciklim u \
temenima*)



fOrientGr1[GG_Graph,neUredjeni_List]:=
  Module[{l,g,rez={},CRD,sus,j,i,ind=True,iso,GNew},
    l=Flatten[FromAdjacencyMatrix[GG[[1]]][[1]],1];
    Do[sus=Select[l,MemberQ[#,i]&];
          sus=Map[If[SameQ[#[[1]],i],#[[2]],#[[1]]]&,sus];
         CRD=Map[Extract[GG[[2]],#]&,sus];
          
      rez=Append[rez,
          fFindOrientation[Prepend[CRD,GG[[2,i]]],Prepend[sus,i]]],
      {i,1,Length[GG[[1]]]}];
    g=Flatten[
        Table[Table[{rez[[j]][[1]],rez[[j]][[i]]},{i,2,Length[rez[[j]]]}],{j,
            Length[rez]}],1];
    g=Union[Map[Sort,g]];(*novi graf- dobijen bez duplih ivica*)
    
    g=Isomorphism[FromUnorderedPairs[Union[neUredjeni]],
    FromUnorderedPairs[g], All];
    j=1;
    If[fFind4[neUredjeni],
      While[j<=Length[g]&&ind,
                   
        GNew=Sort[
            Map[Sort,
              ReplaceAll[neUredjeni,Table[i->g[[j,i]],{i,Length[g[[j]]]}]]]];
             If[fFind4[GNew],ind=False,j++];
        iso=g[[j]]
        ],
      (*ako graf nije 4-
          valentan onda daj trivijalan izomorfizam i sam neu kao nov*)
      
      GNew=neUredjeni;iso=Union[Flatten[neUredjeni]]
      ];
    {GNew,iso,rez}
    ]

(*GG je novi graf,
  g[[j]] odgovarajuci izomorfizam a rez je planar embeding dat ciklim u \
temenima*)
(*za PlanarEmbGraf mi treba samo GG i rez*)


    




  fNasCikl[LL_List]:=Module[{i,kk,l,p},kk=LL;
    (*Print["Nas cikl radi sa: ",kk];*)
    l=kk[[1]];
    kk=Rest[kk];
    Do[p=Position[kk,l[[i+1]]][[1]];
      l=Append[l,kk[[p[[1]],ReplaceAll[Mod[p[[2]]+1,2],0->2]]]];
      kk=Drop[kk,{p[[1]]}],{i,1,Length[kk]-1}];
    (*Print["Nas cikl vraca: ",l];*)
    l
    ]
    
(*pomocna za fPljosni *)
fDtEd[LL_List]:=Map[{LL[[1]],#}&,LL[[2]]]
fSve[LL_List,elL_List]:=Union[Map[MemberQ[LL,#]&,elL]]
fMakeVert[LL_List]:=Map[First[#]&,LL]
fMakeEd[FF_List]:=Module[{i,ed={}},
    Do[
      ed=Append[ed,Take[FF,{i,i+1}]],{i,1,Length[FF]-1}];
      ed=Append[ed,{Last[FF],First[FF]}];
    ed=Map[Sort[#]&,ed]
    ]
    
fMakeVert[LL_List]:=Map[First[#]&,LL]
(* radi od Konture spoljne ,
    adj mat- matrica incidencije
      grafa,UnPair -neuredjeni parovi sa dodatim digonima*)

fPljosni[Kontura_List,adjmat_List,graf_,gUnPair_List]:=
  Module[{DT,isogr,p,n,isoLvert={},i,pos,pom,res={},lepi,tro},
    isogr=Union[gUnPair];
    n=Length[Union[Flatten[isogr]]];
    isogr=Union[isogr, Map[Reverse[#]&,isogr]];
    Do[p=  Select[isogr,#[[1]]==i&];
       isoLvert=Append[isoLvert,
          Flatten[Map[ Take[#,{2}]&,p],1]]
       ,{i,1,n}];
    DT=DelaunayTriangulation[graf[[2]]];
    lepi=Flatten[Map[fDtEd[#]&,DT],1];
    pom=Map[fMakeEd[#[[2]]]&,DT];
      Do[ pom=ReplacePart[pom,Map[Prepend[#,i]&,pom[[i]]],i] ,{i,1,
        Length[pom]}];
    pom=Map[fMakeEd[#]&,Flatten[pom,1]];
    res=Select[pom,fSve[isogr,#]=={True}&];
    res=Map[Union[Flatten[#]]&,res];
    res=Union[Map[Sort[#]&,res]];
    pom=Complement[pom,res];
    pom=Union[Map[Sort[#]&,pom]];
    lepi=Union[Map[Sort[#]&,Complement[lepi,isogr]]];
    res=Map[fMakeEd[#]&,res];
    While[lepi!={},
                pos=Position[pom,lepi[[1]]];
        p=Select[Union[   pom[[pos[[1,1]]]],
                    pom[[pos[[2,1]]]]],
                   (#1!=lepi[[1]])&];
     If[fSve[isogr,p]=={True},
                             res=Append[res,p],
                             pom=ReplacePart[pom,p,pos[[1,1]]];
        pom=Drop[pom,{pos[[2,1]]}]           ];    
        lepi=Rest[lepi]
         ];
    res=Map[fNasCikl[#]&,res];
    res
   ]
     
     
fFacOrient[faces_List]:=Module[{pom,ind,t,newfac,ind1,f1,ind0,indold,i},
    pom={};
    ind=Join[{1},Table[0,{i,Length[faces]-1}]];
    While[MemberQ[ind,0],
      indold=ReplaceAll[ind,0->1];
      newfac=Table[
          If[ind[[i]]>=0,faces[[i]],Reverse[faces[[i]]]],{i,
            Length[ind]}];
         ind1=Abs[ind];
      pom=
        Union[Flatten[
            Table[If[SameQ[ind1[[i]],0],{},
                If[SameQ[ind1[[i]],1],
                  Join[Table[{newfac[[i,j]],newfac[[i,j+1]]},{j,
                        Length[newfac[[i]]]-1}],{{Last[newfac[[i]]],
                        First[newfac[[i]]]}}],
                  Map[Reverse,
                    Join[Table[{newfac[[i,j]],newfac[[i,j+1]]},{j,
                          Length[newfac[[i]]]-1}],{{Last[newfac[[i]]],
                          First[newfac[[i]]]}}]]]],{i,Length[ind]}],1]];
      pom=Map[Reverse,pom];
      f1=Table[
          Join[Table[{newfac[[j]][[i]],newfac[[j]][[i+1]]},{i,
                Length[newfac[[j]]]-1}],{{Last[newfac[[j]]],
                First[newfac[[j]]]}}],{j,Length[newfac]}];
      f1=Drop[f1,1];
      ind0=
        Table[If[
            SameQ[Intersection[f1[[i]],pom],{}]&&
              SameQ[Intersection[Map[Reverse,f1[[i]]],pom],{}],0,1],{i,
            Length[f1]}];
      ind1=
        Table[If[
            SameQ[If[
                  SameQ[Intersection[f1[[i]],pom],{}]&&
        SameQ[Intersection[Map[Reverse,f1[[i]]],pom],{}],0,2],2]&&
              Not[SameQ[Intersection[f1[[i]],pom],{}]],1,-1],{i,
            Length[faces]-1}];
      ind=Join[{1},ind0*ind1];
      ind=indold*ind;
      newfac=
        Table[If[ind[[i]]>=0,faces[[i]],Reverse[faces[[i]]]],{i,
            Length[ind]}];
      ];
    newfac
    ]  



(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #*)
 (*## ## ## ## ## ## # fPlanar Emb ## ## ## ## ## ## ## ## ## ## ## ## ## #*)



(* ## ## # CONVEX DRAW by M.Ochiai, N.Imafuji, N.Morimura ## ## ## # *)


makematrix[mtx_, len_,xy_]:=Block[{i,j,k,d,len1,amt,bmt,lumat},
	len1=Length[mtx];
	amt=Table[0,{i,len1-len},{j,len1-len}];
	bmt=Table[0,{i,len1-len},{j,len}];
	i=len+1;k=1;
	While[i<=len1,
		d=Sum[mtx[[i,j]],{j,len1}];
		For[j=1,j<=len1,j++,
			If[mtx[[i,j]]==1,
				If[j<=len,bmt[[k,j]]=1/d,amt[[k,j-len]]=1/d]]];
		k++;
		i++];
    	lumat=LUDecomposition[IdentityMatrix[len1-len]-amt];
    	Return[LUBackSubstitution[lumat,bmt . xy]]
    ]
rearrange[mtx_,cyl_]:=Block[{i,j,len1,len2,cset,tmm1,tmm2},
	len1=Length[mtx];len2=Length[cyl];
	cset=Complement[Table[i,{i,len1}],cyl];
	tmm1=mtx;
	tmm2=Table[0,{i,len1},{j,len1}];
	For[j=0,j<2,j++,
		For[i=1,i<=len2,i++,tmm2[[i]]=tmm1[[cyl[[i]]]]];
		For[i=len2+1,i<=len1,i++,tmm2[[i]]=tmm1[[cset[[i-len2]]]]];
		tmm1=Transpose[tmm2];
	];
	Return[tmm1]
]
circluarvertices[n_]:=Block[{i,s=Pi/2-2 Pi / n,x=N[2 Pi / n]},
	  Return[
      Chop[N[Map[({{Cos[s],-Sin[s]},{Sin[s],Cos[s]}} . #)&,
            Table[{(Cos[x  i]),(Sin[x  i])},{i,n}]]]]]
	]


Next1[v_List]:=Part[v,2]


FindPlotRange[v_List] :=
  Module[{xmin=Min[Map[First,v]], xmax=Max[Map[First,v]],
			ipsmin=Min[Map[Next,v]], ipsmax=Max[Map[Next,v]]},
		{ {xmin - 0.05 Max[1,xmax-xmin], xmax + 0.05 Max[1,xmax-xmin]},
		  {ipsmin - 0.05 Max[1,ipsmax-ipsmin], 
        ipsmax + 0.05 Max[1,ipsmax-ipsmin]}}] 
showgraphics2D[vtx_,prs_]:=Show[Graphics[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]
		],{AspectRatio->1,PlotRange->FindPlotRange[vtx]}]
showgraphics[vtx_,prs_]:=
		Show[Graphics3D[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]],{AspectRatio->1,PlotRange->All}]
innerproduct[v1_,v2_]:=
  Sqrt[(First[v1]-First[v2])^2+(Next1[v1]-Next1[v2])^2+(Last[v1]-Last[v2])^2]
averagepos[vtx_,vw_]:={
	Apply[Plus,Map[(vtx[[#,1]])&,vw]],
  Apply[Plus,Map[(vtx[[#,2]])&,vw]],
  Apply[Plus,Map[(vtx[[#,3]])&,vw]]}/Length[vw];
centerpos[vtx_]:=Block[{sx=0,sy=0,sz=0,len,i},
    len = Length[vtx];
    For[i=1,i<=len,i++,sx += vtx[[i,1]];sy += vtx[[i,2]]; 
      sz += vtx[[i,3]]];
    Return[{sx/len,sy/len,sz/len}]]
drawgraph[mtx_, cyl_]:=
  Block[{len1,len2,xy0, xy,pairs,mx,xyz,minlen,cnt,ccentar},
	len1=Length[mtx];len2=Length[cyl];
         xy0=circluarvertices[len2];
    	If[len1==len2,xy=circluarvertices[len1];mx=mtx,
      			mx=rearrange[mtx,cyl];
      			xy=Join[xy0,makematrix[mx,len2,xy0]]
      	];
    	pairs=Select[Position[mx, _?(Function[n,n != 0])],OrderedQ];
    	xy=Join[xy0,Chop[normalizedmatrix[mx,xy,Length[cyl], xy0]]];
        showgraphics2D[xy,pairs];
        xyz=Map[(Append[#,1/2])&,xy];
        For[i=1,i<=len1,i++, 
            xyz[[i]]  *= 
            xyz[[i,3]]/(2 xyz[[i,1]]^2+2 xyz[[i,2]]^2+2 xyz[[i,3]]^2)];
             ccentar=centerpos[xyz];
        While[N[innerproduct[{0.0,0.0,0.0},ccentar],20]>0.0000001,
          For[i=1,i<=len1,i++,  
              xyz[[i]] -= ccentar;   
        xyz[[i]]/=innerproduct[{0.0,0.0,0.0},xyz[[i]]]];
        ccentar=centerpos[xyz]
      ];
          If[Length[len1]>30,cnt=10,cnt=5];
          For[j=1,j<cnt,j++,
      		For[i=1,i<=len1,i++,
        			
        xyz[[i]]=(xyz[[i]]+
                averagepos[xyz,
                  Flatten[Position[ mx[[i]], _?(Function[n,n != 0]) ]]])/2.0;
        			xyz[[i]]/=innerproduct[{0.0,0.0,0.0},xyz[[i]]]]];		
    	xyz =N[xyz,20];
    	showgraphics[xyz,pairs]; 
    Return[{mx,xyz}]
    ]	
    
    
    
Edgesnn[Graph[e_,_]]:=e
Vert[Graph[_,v_]]:=v
Ve[Graph[e_,_]]:=Length[e]
Me[Graph[g_,_],___]:=Apply[Plus,Map[(Apply[Plus,#])&,g]]/2
Me[Graph[g_,_],Directed]:=Apply[Plus,Map[(Apply[Plus,#])&,g]]
ChVertices[g_Graph,v_List]:=Graph[Edgesnn[g],v]
ChEdges[g_Graph,e_List]:=Graph[e,Vert[g]]
DFS[v_Integer]:=(dfi[[v]]=cnt++;
    AppendTo[visit,v];
    Scan[(If[dfi[[#]]==0,AppendTo[eedgs,{v,#}];DFS[#]])&,e[[v]]])
DepthFirstTrans[g_Graph,start_Integer,flag_:Vertex]:=
  Block[{visit={},e=ToAdjacencyL[g],eedgs={},dfi=Table[0,{Ve[g]}],cnt=1},
    DFS[start];
    If[flag===Edge,eedgs,visit]]
ToAdjacencyL[Graph[g_,_]]:=
  Map[(Flatten[Position[#,_?(Function[n,n!=0])]])&,g]
ComplGraph[0]:=Graph[{},{}]
ComplGraph[1]:=Graph[{{0}},{{0,0}}]
ComplGraph[n_Integer?Positive]:=CircGraph[n,Range[1,Floor[(n+1)/2]]]
CircGraph[n_Integer?Positive,l_List]:=
  Module[{i,r},r=Prepend[MapAt[1&,Table[0,{n-1}],Map[List,Join[l,n-l]]],0];
    Graph[Table[RotateRight[r,i],{i,0,n-1}],CircVertices[n]]]
EmptyGr[n_Integer?Positive]:=
  Module[{i},Graph[Table[0,{n},{n}],Table[{0,i},{i,(1-n)/2,(n-1)/2}]]]
ComplGraph[l__]:=
  Module[{ll=List[l],t,i,x,rroww,stages=Length[List[l]]},
      t=FoldList[Plus,0,ll];
      Graph[
        Apply[Join,
          Table[rroww=
              Join[Table[1,{t[[i-1]]}],Table[0,{t[[i]]-t[[i-1]]}],
                Table[1,{t[[stages+1]]-t[[i]]}]];
            Table[rroww,{ll[[i-1]]}],{i,2,stages+1}]],
        Apply[Join,
          Table[Table[{x,i-1+(1-ll[[x]])/2},{i,ll[[x]]}],{x,stages}]]]]/;
    TrueQ[Apply[And,Map[Positive,List[l]]]]&&(Length[List[l]]>1)
CircVertices[0]:={}
CircVertices[n_Integer]:=
  Module[{i,x=N[2 Pi/n]},Chop[Table[N[{(Cos[x i]),(Sin[x i])}],{i,n}]]]
CircVertices[Graph[g_,_]]:=Graph[g,CircVertices[Length[g]]]
showgraphics2D[vtx_,prs_]:=Show[Graphics[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]
		],{AspectRatio->1,PlotRange->FindPlotRange[vtx]}]
defese[v_Integer]:=(dfi[[v]]=cnt++;
    AppendTo[visit,v];
    Scan[(If[dfi[[#]]==0,AppendTo[eedgs,{v,#}];defese[#]])&,e[[v]]])
depthfirstsearch[mtx_List,start_Integer]:=
  Block[{visit={},e,eedgs={},dfi=Table[0,{Length[mtx]}],cnt=1},
    e=Map[(Flatten[Position[#,_?(Function[n,n!=0])]])&,mtx];
    defese[start];
    Return[eedgs]]
getalledges[adjmat_]:=
  Select[Position[ adjmat, _?(Function[n,n != 0]) ],OrderedQ]
gettreeedges[adjmat_,vtx_,stpt_]:=
  If[Length[vtx]==0,
    DepthFirstTrans[Graph[adjmat,Vert[ComplGraph[Length[adjmat]]]],stpt,
      Edge],
    DepthFirstTrans[Graph[adjmat,vtx],stpt,Edge]]
getleafedges[alledges_,treeedges_]:=
  
  Complement[alledges,
    Union[Select[treeedges, OrderedQ],
      Map[Reverse,Complement[treeedges,Select[treeedges, OrderedQ]]]]]
ShowTree[adjmat_,vtx_,stpt_]:=
  showgraphics2D[vtx,gettreeedges[adjmat,vtx,stpt]]
Dijks[g_Graph,start_Integer]:=First[Dijks[g,{start}]]
Dijks[g_Graph,l_List]:=
  Module[{x,start,e=ToAdjacencyL[g],i,p,pparentt,untraversed},
    p=Edgesnn[PathConditionGraph[g]];
    Table[start=l[[i]];
      pparentt=untraversed=Range[Ve[g]];
      distan=p[[start]]; distan[[start]]=0;
      Scan[(pparentt[[#]]=start)&,e[[start]]];
      While[untraversed!={},x=First[untraversed];
        Scan[(If[distan[[#]]<distan[[x]],x=#])&,untraversed];
        untraversed=Complement[untraversed,{x}];
        Scan[(If[distan[[#]]>distan[[x]]+p[[x,#]],
                distan[[#]]=distan[[x]]+p[[x,#]];
                pparentt[[#]]=x])&,e[[x]]];];
      {pparentt,distan},{i,Length[l]}]]
ShortPath[g_Graph,s_Integer,e_Integer]:=
  Module[{pparentt=First[Dijks[g,s]],i=e,lst={e}},
    While[(i!=s)&&(i!=pparentt[[i]]),
      PrependTo[lst,pparentt[[i]]];
      i=pparentt[[i]]];
    If[i==s,lst,{}]]
    
MakeUndir[Graph[g_,v_]]:=
  Module[{i,j,n=Length[g]},
    Graph[Table[
        If[g[[i,j]]!=0||g[[j,i]]!=0,1,0],{i,n},{j,n}],v]]
FromOrdP[{}]:=Graph[{},{}]
FromOrdP[l_List,v_List]:=
Graph[MapAt[1&,Table[0,{Length[v]},{Length[v]}],l],v]
FromUnordP[l_List]:=MakeUndir[FromOrdP[l]]
FromUnordP[l_List,v_List]:=MakeUndir[FromOrdP[l,v]]
ShortestPathSTree[g_Graph,s_Integer]:=
  Module[{pparentt=First[Dijks[g,s]],i},
    FromUnordP[Map[({#,pparentt[[#]]})&,Complement[Range[Ve[g]],{s}]],
      Vert[g]]]
AllPairsShortestP[g_Graph]:=
  Module[{p=Edgesnn[PathConditionGraph[g]],i,j,k,n=Ve[g]},
      Do[p=Table[Min[p[[i,k]]+p[[k,j]],p[[i,j]]],{i,n},{j,n}],{k,n}];
      p]/;Min[Edgesnn[g]]<0
AllPairsShortestP[g_Graph]:=Map[Last,Dijks[g,Range[Ve[g]]]]
RemoveSelfL[g_Graph]:=Module[{i,e=Edgesnn[g]},Do[e[[i,i]]=0,{i,Ve[g]}];
    Graph[e,Vert[g]]]
PathConditionGraph[Graph[e_,v_]]:=
  RemoveSelfL[Graph[ReplaceAll[e,0->Infinity],v]]
characteristiccyclematrix[adjmat_,vtx_,stpt_]:=
  Block[{i,len,len1,len2,at,al,alle,tttree,leaf,idm,bpm},
    len=Length[adjmat];
    alle=getalledges[adjmat];
    tttree=depthfirstsearch[adjmat,stpt];
    leaf=getleafedges[alle,tttree];
    len1=Length[tttree];
    len2=Length[leaf];
    at=Table[0,{i,1,len},{j,1,len1}];
    al=Table[0,{i,1,len},{j,1,len2}];
    For[i=1,i<=len1,i++,Map[(at[[#,i]]=1)&,tttree[[i]]]];
    For[i=1,i<=len2,i++,Map[(al[[#,i]]=1)&,leaf[[i]]]];
    idm=IdentityMatrix[len2];
    bpm=Mod[Transpose[Inverse[Delete[at,len],Modulus->2] . Delete[al,len]],
        2];
    Return[Map[(Flatten[Append[idm[[#]],bpm[[#]]]])&,Table[i,{i,1,len2}]]]
    ]
showgraphics2D[vtx_,prs_]:=Show[Graphics[
			Join[{PointSize[0.02]},Map[Point,Chop[vtx]],
					Map[(Line[Chop[vtx[[#]]]])&,prs]]
		],{AspectRatio->1,PlotRange->FindPlotRange[vtx]}]
ToOrdPairs[g_Graph]:=Position[Edgesnn[g],_?(Function[n,n!=0])]
TranslateVert[v_List,{x_,y_}]:=Map[(#+{x,y})&,v]
TranslateVert[Graph[g_,v_],{x_,y_}]:=Graph[g,TranslateVert[v,{x,y}]]
DilateVert[v_List,d_]:=(d*v)
DilateVert[Graph[e_,v_],d_]:=Graph[e,DilateVert[v,d]]
NormalizeVert[v_List]:=Module[{v1},v1=TranslateVert[v,{-Min[v],-Min[v]}];
    DilateVert[v1,1/Max[v1,0.01]]]
NormalizeVert[Graph[g_,v_]]:=Graph[g,NormalizeVert[v]]
PointsAndLines[Graph[e_List,v_List]]:=
  Module[{pairs=ToOrdPairs[Graph[e,v]]},
    Join[{PointSize[0.025]},Map[Point,Chop[v]],
      Map[(Line[Chop[v[[#]]]])&,pairs]]]
ShowLabelGraph[g_Graph]:=ShowLabelGraph[g,Range[Ve[g]]]
ShowLabelGraph[g1_Graph,labels_List]:=
  Module[{pairs=ToOrdPairs[g1],g=NormalizeVert[g1],v},v=Vert[g];
    Show[Graphics[
        Join[PointsAndLines[g],Map[(Line[Chop[v[[#]]]])&,pairs],
          GraphLabels[v,labels]]],{AspectRatio->1,
        PlotRange->FindPlotRange[v]}]]
GraphLabels[v_List,l_List]:=
  Module[{i},Table[Text[l[[i]],v[[i]]-{0.03,0.03},{0,1}],{i,Length[v]}]]


    
     (*## ## ## ## ## ## # fPlanar Emb ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
#*)

fFindPre[trazi_Integer, gde_Integer, CC_List] := 
  Module[{l = CC[[gde]], nadjen, pom},
    pom = Position[l, trazi][[1, 1]];
   If[pom == 1, nadjen = Last[l], nadjen = l[[pom - 1]]]; nadjen
    ]
    
fCloseFace[CC_List, Ed_List, PP_List] := 
  Module[{sledeci, izbaci = {}, poc = PP},
    sledeci = fFindPre[poc[[1]], poc[[2]], CC];
    While[sledeci != poc[[1]],
            izbaci = Append[izbaci, {Last[poc], sledeci}];
           poc = Append[poc, sledeci];
           sledeci = fFindPre[poc[[-2]], Last[poc], CC]
            ];
    izbaci = Append[izbaci, {Last[poc], sledeci}];
    {poc, Complement[Ed, izbaci]}
    ]
    
fFaces[CiklL_List] := Module[{Ul = CiklL, t, res, konacna = {}},
    (*Lista svih stranica*)
    t = Sort[
        Flatten[Table[
            Table[{i, Ul[[i, j]]}, {j, Length[Ul[[i]]]}], {i, Length[Ul]}], 
          1]];
        poc = First[t];  
        t = Rest[t];
      While[ t != {},
         res = fCloseFace[CiklL, t, poc];
        konacna = Append[konacna, res[[1]]];
        t = res[[2]];
      If[t != {}, poc = First[t]; t = Rest[t]]
      ];
    konacna
    ]  
    
fWrGraph[g_List, file_] := Module[{g1, edg, v, i, x, y, p},
    g1 = FromUnorderedPairs[g];
    p = Length[Union[Flatten[g]]];
    edg = ToAdjacencyMatrix[g1];
    edg = Table[Flatten[Position[edg[[i]], 1]], {i, Length[edg]}];
     v = Flatten[N[NormalizeVertices[g1[[2]]]], 1]; 
    OpenWrite[file];
    WriteString[file, "N=", ToString[p] <> "\n"];
    Do[WriteString[file, "", ToString[i - 1] <> ":"];
       {x, y} = Chop[v[[i]]];
      Scan[(WriteString[file, " ", ToString[# - 1]]) &, edg[[i]]] ;
      WriteString[file, " -1"];
      Write[file], {i, Length[v]}];
    Close[file];]
    
fPlanEmbedding[g_List] := Module[{vv, vv1, vv2, p,i},
    fWrGraph[g, "graph.txt"]; 
    Run["planarity" <> " graph.txt" <> " graph1.txt"];
    vv = Import["graph1.txt"];
    DeleteFile["graph.txt"];
    DeleteFile["graph1.txt"];
    vv = StringReplace[
        vv, {"\n" -> "", " -1" -> "}", ": " -> "{", " " -> ","}];
    vv1 = StringPosition[vv, "{"];
    vv2 = StringPosition[vv, "}"];
    p = Length[vv1];
    vv = StringDrop[vv, vv1[[1, 1]] - 1];
    For[i = 1, i < p, vv1 = Drop[StringPosition[vv, "{"], 1];
      vv2 = Drop[StringPosition[vv, "}"], -1]; 
      vv = StringReplacePart[vv, ",", {vv2[[i, 1]] + 1, vv1[[i, 1]] - 1}]; 
      i++];
    vv = ToExpression["{" <> vv <> "}"] + 1;
    vv
    ]
    
fPlanarEmb[Ul_] := Module[{vv, kk,i},
    kk = fGraphInc[Ul][[1]];
    vv = fPlanEmbedding[Union[kk]];
    vv = {kk, kk, Union[Flatten[kk]], 
        Table[Join[{i}, vv[[i]]], {i, Length[vv]}], 
        fFaces[vv]};
    vv
    ]


fPlanarEmbNew[Ul_] := Module[{vv, kk,i},
    kk = fForGraphKL[Ul][[1]];
    vv = fPlanEmbedding[Union[kk]];
    vv = {kk, kk, Union[Flatten[kk]], 
        Table[Join[{i}, vv[[i]]], {i, Length[vv]}], 
        fFaces[vv]};
    vv
    ]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## #*)
    
 fPlanarEmbKL[Ulaz_] := Module[{vv,kk,i,Ul},
 If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    kk = fGraphInc[Ul][[1]];
    vv = fPlanEmbedding[Union[kk]];
    vv = {kk,Table[Join[{i}, vv[[i]]], {i, \
Length[vv]}],fFaces[vv]};
    vv
    ]
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## #*)    

(*## ## # fPlanarEmbGraph ## #*)

fPlanarEmbGraph[Ul_List] := Module[{vv,kk,i},
    kk = fPlanEmbedding[Union[Ul]];
    vv = Table[Join[{i}, kk[[i]]], {i, Length[kk]}];
    {Ul,vv,fFaces[kk]}
    ]
    




(*## ## ## ## ## # DrawPLanarEmbKL ## ## ## ## ## ## ## ## ## ## ##*)



 
 Dijks[g_Graph,start_Integer]:=First[Dijks[g,{start}]]

Dijks[g_Graph,l_List]:=
  Module[{x,start,e=ToAdjacencyL[g],i,p,pparentt,untraversed},
    p=Edgesnn[PathConditionGraph[g]];
    Table[start=l[[i]];
      pparentt=untraversed=Range[Ve[g]];
      distan=p[[start]]; distan[[start]]=0;
      Scan[(pparentt[[#]]=start)&,e[[start]]];
      While[untraversed!={},x=First[untraversed];
        Scan[(If[distan[[#]]<distan[[x]],x=#])&,untraversed];
        untraversed=Complement[untraversed,{x}];
        Scan[(If[distan[[#]]>distan[[x]]+p[[x,#]],
                distan[[#]]=distan[[x]]+p[[x,#]];
                pparentt[[#]]=x])&,e[[x]]];];
      {pparentt,distan},{i,Length[l]}]]
    
findfacialcycle[adjmat_,vtx_,stpt_]:=
  Block[{i,j,len,len1,len2,sptreem,ssptr,al,alle,tttree,leaf,PPT,fpath,LPT,
      cccount,adjmat1,vtx0,vtx1},
    len=Length[adjmat];
    alle=getalledges[adjmat];
    tttree=depthfirstsearch[adjmat,stpt];
    leaf=getleafedges[alle,tttree];
    len1=Length[tttree];
    len2=Length[leaf];
    sptreem=Table[0,{i,1,len},{j,1,len}];
    For[i=1,i<=len1,i++,sptreem[[tttree[[i,1]],tttree[[i,2]]]]=1];
    sptreem=sptreem + Transpose[sptreem];
    If[Length[vtx]==0,vtx0=Vert[ComplGraph[len]],vtx0=vtx];
    ssptr=Graph[sptreem,vtx0];
    fpath={};
 For[i=1,i<=len2,i++,
      cccount=0;adjmat1=.;vtx1=.;
      PPT=ShortPath[ssptr,leaf[[i,1]],leaf[[i,2]]];
       LPT=Map[({#})&,PPT];
     For[j=1,j<=len2,j++,
        If[
          i!=j && MemberQ[PPT,leaf[[j,1]]] && 
            MemberQ[PPT,leaf[[j,2]]],cccount++];
        If[cccount==1,Break[]]];
      If[Length[adjmat]==Length[PPT],adjmat1=={};vtx1={},
        adjmat1=Delete[Transpose[Delete[adjmat,
                LPT]],LPT]];
      vtx1=Delete[vtx0,LPT];
     If[adjmat1=={} || vtx1=={} ,cccount++,
        If[Length[DepthFirstTrans[Graph[adjmat1,vtx1],1]]!=Length[adjmat1],
          cccount++]];
      If[cccount==0,If[Length[fpath]<Length[PPT],fpath=PPT]]
      ];
    Return[fpath]
    ]  
    

DrawPlanarEmbKL[Ul_]:=
  Module[{q=0,ind,Ul1,neu,neu1,adjmat,cikl,len1,len2,xy,xy0,mx,pairs,vv},
    Ul1=Ul;
    If[SameQ[Head[Ul],String],
      ind=Length[
            ReadList[
              StringToStream[
                StringReplace[
                  Ul,{"*"->" 6 ","."->" 6 ","("->" 6 ",")"->" 6 ",
                    ","->" 6 "}]],Number]]==1
       ];
If[SameQ[Head[Ul],List],
      If[Length[Ul[[1]]]==1,q=Ul[[2]],q=iteratedTake[Ul[[2]],Ul[[1]]]];
          ind=SameQ[q,Dowker[ToString[Apply[Plus,Ul[[1]]]]][[4]]];
         If[ind,Ul1=ToString[Apply[Plus,Ul[[1]]] ]]
      ];
    If[ind==False,
           neu=fGaussExtSigns[Ul];
           neu=fGrInc[neu][[1]];
            neu1=neu;
           adjmat=ToAdjacencyMatrix[FromUnorderedPairs[Union[neu]]];
           cikl=findfacialcycle[adjmat,NULL,1];
           len1=Length[adjmat];len2=Length[cikl];
          If[len1==len2,xy=circluarvertices[len1];  mx=adjmat,
        			mx=rearrange[adjmat,cikl];
        			xy0=circluarvertices[len2];
        			xy=Join[xy0,makematrix[mx,len2,xy0]]
                         ];			
           pairs=Position[mx, _?(Function[n,n != 0])];
           vv=Graph[mx,xy];
        vv=ShowLabelGraph[Graph[mx,xy]],
    vv=ShowLabeledGraph[FromUnorderedPairs
    [fMakeEd[Range[ToExpression[Ul1]]]]]
      ];
     vv ](*16.8.2003*)
 
 (*## ## ## ## # Draw Planar Embeding of the Graph ## ## ## ## ## ## ##*)
 
 
 fNewGraph[LL_,LL1_]:=Module[{gg,pp,qq},
    pp=Flatten[
        Table[Table[{j,Flatten[Position[LL[[j]],1]][[i]]},{i,
              Length[Flatten[Position[LL[[j]],1]]]}],{j,Length[LL]}],1];
    pp=Union[Map[Sort,pp]];
    pp=Table[{pp[[i]]},{i,Length[pp]}];
    qq=Table[{LL1[[i]]},{i,Length[LL1]}];
    gg=Graph[pp,qq];
    gg
    ]
 
DrawPlanarEmbGraph[Ul_]:=
  Module[{ind,neu,adjmat,cikl,len1,len2,xy,mx,xy0,pairs,vv,dd},
  dd=Union[Select[Ul,Length[Position[Ul,#]]>1&]];





    neu=Union[Ul];
   vv=If[IsomorphicQ[FromUnorderedPairs[neu],
   Cycle[Max[Union[Flatten[neu]]]]],
  vv=ShowLabeledGraph[Highlight[Cycle[Max[Union[Flatten[neu]]]],{dd}]],
    ind=Union[Map[Count[Ul,#]&,Ul]];
    If[SameQ[ind,{4}]||SameQ[ind,{2}],
           (*slucaj 2 li -2 ili CIKL *)        
      ShowLabeledGraph[FromUnorderedPairs[neu]],
        adjmat=ToAdjacencyMatrix[FromUnorderedPairs[neu]];
      cikl=findfacialcycle[adjmat,NULL,1];
      len1=Length[adjmat];len2=Length[cikl];
      If[len1==len2,xy=circluarvertices[len1];
             mx=adjmat,
            mx=rearrange[adjmat,cikl];xy0=circluarvertices[len2];
            xy=Join[xy0,makematrix[mx,len2,xy0]]
        ];
      pairs=Position[mx,_?(Function[n,n!=0])];  
     vv=ShowLabeledGraph[Highlight[fNewGraph[mx,xy],{dd}]]
     ]];
   vv ]
(* 26.06.2006 *)


DrawPlanarEmbGraphNew[Ul_]:=
  Module[{ind,neu,adjmat,cikl,len1,len2,xy,mx,xy0,pairs,vv,dd,ss,ss1,ss2},
  dd=Union[Select[Ul,Length[Position[Ul,#]]>1&]];

ss=fPlanarEmbGraph[Ul][[3]];
ss1=Map[Range,Map[Length,ss]];
ss2=Flatten[Position[Table[Complement[ss1[[i]],ss[[i]]],{i,Length[ss1]}],{}]][[1]];
dd=ReplaceAll[dd,Table[ss[[ss2,i]]->i,{i,Length[ss[[ss2]]]}]];





    neu=Union[Ul];
   vv=If[IsomorphicQ[FromUnorderedPairs[neu],
   Cycle[Max[Union[Flatten[neu]]]]],
  vv=ShowLabeledGraph[Highlight[Cycle[Max[Union[Flatten[neu]]]],{dd}]],
    ind=Union[Map[Count[Ul,#]&,Ul]];
    If[SameQ[ind,{4}]||SameQ[ind,{2}],
           (*slucaj 2 li -2 ili CIKL *)        
      ShowLabeledGraph[FromUnorderedPairs[neu]],
        adjmat=ToAdjacencyMatrix[FromUnorderedPairs[neu]];
      cikl=findfacialcycle[adjmat,NULL,1];
      len1=Length[adjmat];len2=Length[cikl];
      If[len1==len2,xy=circluarvertices[len1];
             mx=adjmat,
            mx=rearrange[adjmat,cikl];xy0=circluarvertices[len2];
            xy=Join[xy0,makematrix[mx,len2,xy0]]
        ];
      pairs=Position[mx,_?(Function[n,n!=0])];  
     vv=ShowLabeledGraph[Highlight[fNewGraph[mx,xy],{dd}]]
     ]];
   vv ]


(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # *)  
(*## ## ## ## ## # fSignsKL ## ## ## ## ## ## ## ## ## ## ##*)
(* fSignsKL racuna znake Alt KL. 
      Ulaz: Dow u smislu Knotscape. 
            Izlaz: Dow sa znacima. Dow za Alt KL unosimo bez znakova, 
  a za Nonalt samo znake tacaka koje su izmenjene u odnosu na Alt (sign \
change) *) 

UnsortedUnion[x_]:=Module[{f},f[y_]:=(f[y]=Sequence[];y);f/@x]

(* fSignsKL[Dow_List]:=
  Module[{DS,DowA,OldGaussExt,Emb,iso,i,j,CC,cik,q=1,QQQ={},uu,zakraj,eee},
    DS=Table[Sign[Dow[[2]][[i]]],{i,Length[Dow[[2]]]}];
    DowA=Abs[Dow];
    If[SameQ[Length[DowA[[1]]],1],OldGaussExt={fGaussExtSigns[DowA]},
      OldGaussExt=fGaussExtSigns[DowA]];
      eee=fPlanarEmb[DowA];
      Emb=eee[[4]];
     iso=eee[[3]];
    Emb=Sort[
        ReplaceAll[Emb,
          Table[iso[[i]]->i,{i,Length[iso]}]],#1[[1]]<#2[[1]]&];
    zakraj=Emb;
    Emb=Map[Rest,Emb];(*izomorfni PlanarEmb*)
    
    CC=Table[Position[OldGaussExt,i],{i,
          Length[Union[Flatten[OldGaussExt]]]}];
    CC=Map[If[OddQ[#[[1,2]]],#,Reverse[#]]&,CC];
    (*izgleda da pravimo novi Dow tj.ovih parova cemo uzimati druge \
koordinate za nove parove kod kojih je prvi neparan i po njima cemo \
sortirati*)

        CC=Table[{First[RotateLeft[OldGaussExt[[CC[[i,1,1]]]],
        CC[[i,1,2]]-2]],
          First[RotateLeft[OldGaussExt[[CC[[i,2,1]]]],CC[[i,2,2]]-2]],
          First[RotateLeft[OldGaussExt[[CC[[i,1,1]]]],CC[[i,1,2]]]],
          First[RotateLeft[OldGaussExt[[CC[[i,2,1]]]],CC[[i,2,2]]]]},{i,
          Length[CC]}];
    (*pravimo malu tabelicu ??? pomocu GaussExt i CC*)
    
    cik=Table[UnsortedUnion[CC[[i]]],{i,Length[CC]}];
    (*cik su orginalni CC bez dvojnih veza*)
    uu=Map[UnsortedUnion,Emb];
    (*uporedjujemo orijentaciju cik i uu.Ako neki 
    deo cik i Emb ima samo dva 
clana,pisemo 0*)
QQQ=Table[If[SameQ[Length[Union[cik[[j]]]],2],
          			0,
          			
          If[MemberQ[Table[RotateLeft[cik[[j]],i],{i,Length[cik[[j]]-1]}],
              uu[[j]]],
            		               1,
                                          -1]],
        		{j,Length[cik]}];
    QQQ=If[SameQ[Union[QQQ],{0}],Table[1,{i,Length[QQQ]}],QQQ];
    (*temena sa nulama cuvaju znak prethodnika*)
    (*Sad jos treba da \
odredimo u sta idu nule !!!!!*)
    (*u Emb trazimo kom lancu digona to teme \
sa 0 pripada i uzmemo taj znak*)
      pp=Flatten[Union[Position[QQQ,0]]];
    Do[ cik=Rest[zakraj[[pp[[i]]]]];
            cik=Select[cik,Not[SameQ[QQQ[[#]],0]]&][[1]];
        QQQ=ReplacePart[QQQ,QQQ[[cik]],pp[[i]]]
       ,{i,1,Length[pp]}];
    uu=Table[Position[OldGaussExt,i],{i,Length[QQQ]}];
    Do[OldGaussExt=
        ReplacePart[OldGaussExt,i*QQQ[[i]],{uu[[i,1]],uu[[i,2]]}],{i,
        Length[uu]}];
    OldGaussExt=fDowfromGaussExt[OldGaussExt];
    OldGaussExt={OldGaussExt[[1]],-Flatten[OldGaussExt[[2]]]*Flatten[DS]}] *)

fSignsKL[Ul_List]:=Module[{pp0,pp,pp1,pp2},
pp0=DTCode@@If[SameQ[Length[Ul[[1]]],1],
iteratedTake[Ul[[2]],Ul[[1]]][[1]],iteratedTake[Ul[[2]],Ul[[1]]]];
pp=PositiveQ/@PD[pp0];
pp1=ReplaceAll[Table[pp[[i]],{i,Length[pp]}],{True->1,False->-1}];
pp2={Ul[[1]],Abs[Ul[[2]]]*pp1};
pp2]

(*4.10.2011*)


(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ##*)
(*************************)
(*pomocna za GraphKL*)
    
fNeighFac[LL_List,ff_List]:=Module[{l,p,dig,Lose,pom,i},
    dig=Union[Flatten[Select[LL,Length[#]==2&],1]];
    p=Position[LL,ff][[1,1]];
    l=Complement[Select[LL,Length[Intersection[#,ff]]>=1&],{ff}];
    Lose=Select[l,Length[#]!=2&];
    Lose=Map[{#,Complement[Intersection[#,ff],dig]}&,Lose];
    Lose=Select[Lose,#[[2]]=={}&];
    If[Lose!={},Lose=Map[Take[#,1]&,Lose]];
    Lose=Flatten[Lose,1];
    l=Complement[l,Lose];
    l=Map[Position[LL,#][[1,1]]&,l];
    Prepend[l,p]
    ]
fMakeCircEd[LL_List]:=Module[{i,res={}},
       Do[  res=Append[res,{LL[[1]],LL[[i]]}],{i,2,Length[LL]}];
    res]
  (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ## ## ## ## ## # fGraphKL ## ## ## ## ## ## ## ## ## ## ##*)

fMakeEd[FF_List]:=Module[{i,ed={}},
    Do[
      ed=Append[ed,Take[FF,{i,i+1}]],{i,1,Length[FF]-1}];
      ed=Append[ed,{Last[FF],First[FF]}];
    ed=Map[Sort[#]&,ed]
    ]
    
(*fGraphKL:
    input Conway ili Dow,output:graf KL*)
   
fGraphKL[Ulaz_]:=
  Module[{Ul,ind,i,pe,faces,NewGr,pom,facesEd,pos,p,mm,nbs},
    If[Head[Ulaz]==List,Ul=Abs[Ulaz];ind=Ul[[1]]];
    If[Head[Ulaz]==String,Ul=Ulaz;ind=fDToD[Ulaz][[1]]];
    If[IsomorphicQ[FromUnorderedPairs[Union[fGraphInc[Ulaz][[1]]]],
    Cycle[Max[Flatten[Union[fGraphInc[Ulaz][[1]]]]]]],
  pom= Union[fGraphInc[Ulaz][[1]]],
    If[SameQ[ind,{1,1}],
      (* slucaj "2","-2" *)
      pom={{1,2},{1,2}},
      (*ostali *)
      pe=fPlanarEmb[Ul];
      NewGr=pe[[2]];
      pom=Union[NewGr];
      faces=pe[[5]];(*Print[faces];*)             
      Do[If[Count[NewGr,pom[[i]]]==2,faces=Append[faces,pom[[i]]]],{i,1,\
          Length[pom]}];
      (*fac su faces+digons*)(*sad treba da nadjemo susedne*)
      facesEd=
        Map[fMakeEd,faces];
      nbs=Map[fNeighFac[facesEd,#]&,facesEd];
      (*cikli susednih oko faca*)pom=
        Map[{#,0,0}&,Range[Length[nbs]]];(*brojevi svih ivica*)pom=
        ReplacePart[pom,{1,0,1},1];(*bela je 0;crna je 1*)(*redni broj,
        nula za boju,nula za da li je do sada obojena*)While[
        Union[Flatten[Map[Take[#,-1]&,pom]]]!={1},
        Do[If[pom[[i,3]]==1,
            pom=Map[If[#[[3]]!=1&&MemberQ[nbs[[i]],#[[1]]],{#[[1]],
                      1-pom[[i,2]],1},#]&,pom]],{i,1,Length[nbs]}]];
      pom=
        Select[pom,#[[2]]==1&];(*pom je lista koj asdarzi crne \
pljosni*)
       pom=Flatten[Map[Take[#,1]&,pom]];
      (*redni brojevi crnih pljosni*)faces=
        Map[Prepend[#,Position[faces,#][[1,1]]]&,faces];
      faces=Select[faces,MemberQ[pom,#[[1]]]&];
      Do[faces=ReplacePart[faces,Flatten[{i,Rest[faces[[i]]]}],i],{i,1,
          Length[faces]}];
      faces=Flatten[Map[fMakeCircEd,faces],1];
      (*lista parova {stranica,teme} a mi cemo sad dobiti {str,str}*)n=
        Length[Union[Take[Flatten[faces],{2,Length[Flatten[faces]],2}]]];
      (*ovo su u sustini smao sva temena*)pom={};
      Do[p=Select[faces,#[[2]]==i&];
        pos=Map[Position[faces,#][[1,1]]&,p];
        pom=Append[pom,{p[[1,1]],p[[2,1]]}];
        faces=Drop[faces,{pos[[1]]}];
        faces=Drop[faces,{pos[[2]]-1}],{i,1,n}];
    (*  ShowLabeledGraph[FromUnorderedPairs[pom]] *)
      ]];
   (*   Print["Graph of KL: ", pom]; *)
      mm=Union[Select[pom,Length[Position[pom,#]]>1 &]];
    (*  Print["Double edges: ", mm]; *)
      pom=Sort[pom];
      pom]
  

fForGraphKL[Con_]:=Module[{uu,pp,ll,pp1,pp2,res,jj,jj1,ttt,ttt1,i,j},
    uu=fConwayToPD[Con];
    pp=Table[Table[uu[[j,i]],{i,4}],{j,Length[uu]}];
    ll=Length[Union[Flatten[pp]]];
    pp1=Map[Sort,Table[Map[First,Position[pp,i]],{i,ll}]];
    pp2=Union[
        Flatten[Table[
            Map[Sort,{{pp1[[pp[[i,1]]]],pp1[[pp[[i,2]]]]},{pp1[[pp[[i,3]]]],
                  pp1[[pp[[i,4]]]]}}],{i,Length[pp]}],1]];
    res={Sort[pp1],pp2};
    res]
 
  fGraphKLNew[Ulaz_]:=
    Module[{Ul,ind,i,pe,faces,NewGr,pom,facesEd,pos,p,mm,ttt,nbs,jj,jj1,ttt1},
    ttt=fForGraphKL[Ulaz];
    If[Head[Ulaz]==String,Ul=Ulaz;ind=fDToD[Ulaz][[1]]];
   
 
    If[IsomorphicQ[FromUnorderedPairs[Union[fGraphInc[Ulaz][[1]]]],
    Cycle[Max[Flatten[Union[fGraphInc[Ulaz][[1]]]]]]],
  pom= Union[ttt[[1]]],
    If[SameQ[ind,{1,1}],
      (* slucaj "2","-2" *)
      pom={{1,2},{1,2}},
      (*ostali *)
      pe=fPlanarEmbNew[Ul];
      NewGr=pe[[2]];
      pom=Union[NewGr];
      faces=pe[[5]];              
      Do[If[Count[NewGr,pom[[i]]]==2,faces=Append[faces,pom[[i]]]],{i,1,\
          Length[pom]}];
(* fac su faces+digons*)(*sad treba da nadjemo susedne*)
      facesEd=
        Map[fMakeEd,faces];
      nbs=Map[fNeighFac[facesEd,#]&,facesEd];
      (*cikli susednih oko faca*)pom=
        Map[{#,0,0}&,Range[Length[nbs]]];(*brojevi svih ivica*)pom=
        ReplacePart[pom,{1,0,1},1];(*bela je 0;crna je 1*)(*redni broj,
        nula za boju,nula za da li je do sada obojena*)While[
        Union[Flatten[Map[Take[#,-1]&,pom]]]!={1},
        Do[If[pom[[i,3]]==1,
            pom=Map[If[#[[3]]!=1&&MemberQ[nbs[[i]],#[[1]]],{#[[1]],
                      1-pom[[i,2]],1},#]&,pom]],{i,1,Length[nbs]}]];
      pom=
        Select[pom,#[[2]]==1&];(*pom je lista koj asdarzi crne \
pljosni*)

(*  Print[faces];  *)
 
       pom=Flatten[Map[Take[#,1]&,pom]];
       
    (*    Print[pom]; *)
         
      jj=Table[faces[[pom[[i]]]],{i,Length[pom]}];
      jj1=Map[fMakeCycle,jj];
      ttt1=ttt[[2]];
      pom=If[SameQ[
        Union[Flatten[
            Table[If[SameQ[Intersection[jj1[[i]],ttt1[[j]]],Union[ttt1[[j]]]],\
\
\
\
\

                pom[[i]],{}],{i,Length[jj1]},{j,Length[ttt1]}]]],pom],pom,
      Complement[Range[Length[faces]],pom]];
       
    
          
      (*redni brojevi crnih pljosni*)
      
      faces=
        Map[Prepend[#,Position[faces,#][[1,1]]]&,faces];
      faces=Select[faces,MemberQ[pom,#[[1]]]&];
      Do[faces=ReplacePart[faces,Flatten[{i,Rest[faces[[i]]]}],i],{i,1,
          Length[faces]}];
      faces=Flatten[Map[fMakeCircEd,faces],1];
      (*lista parova {stranica,teme} a mi cemo sad dobiti {str,str}*)n=
        Length[Union[Take[Flatten[faces],{2,Length[Flatten[faces]],2}]]];
      (*ovo su u sustini smao sva temena*)pom={};
      Do[p=Select[faces,#[[2]]==i&];
        pos=Map[Position[faces,#][[1,1]]&,p];
        pom=Append[pom,{p[[1,1]],p[[2,1]]}];
        faces=Drop[faces,{pos[[1]]}];
        faces=Drop[faces,{pos[[2]]-1}],{i,1,n}];
    (*  ShowLabeledGraph[FromUnorderedPairs[pom]] *)
      ]];
   (*   Print["Graph of KL: ", pom]; *)
      mm=Union[Select[pom,Length[Position[pom,#]]>1 &]];
    (*  Print["Double edges: ", mm]; *)
      pom=Sort[pom];
     
      pom]
      
      
 
 fDualGraphKLNew[Ulaz_]:=
    Module[{Ul,nbs,jj,jj1,ttt1,ind,i,pe,faces,NewGr,pom,facesEd,pos,p,mm,ttt},
    ttt=fForGraphKL[Ulaz];
   
    If[Head[Ulaz]==String,Ul=Ulaz;ind=fDToD[Ulaz][[1]]];
    If[IsomorphicQ[FromUnorderedPairs[Union[fGraphInc[Ulaz][[1]]]],
 Cycle[Max[Flatten[Union[fGraphInc[Ulaz][[1]]]]]]],
  pom= Union[ttt[[1]]],
    If[SameQ[ind,{1,1}],
      (* slucaj "2","-2" *)
      pom={{1,2},{1,2}},
      (*ostali *)
      pe=fPlanarEmbNew[Ul];
      NewGr=pe[[2]];
      pom=Union[NewGr];
      faces=pe[[5]];(*Print[faces];*)             
      Do[If[Count[NewGr,pom[[i]]]==2,faces=Append[faces,pom[[i]]]],{i,1,\
          Length[pom]}];
(* fac su faces+digons*)(*sad treba da nadjemo susedne*)
      facesEd=
        Map[fMakeEd,faces];
      nbs=Map[fNeighFac[facesEd,#]&,facesEd];
      (*cikli susednih oko faca*)pom=
        Map[{#,0,0}&,Range[Length[nbs]]];(*brojevi svih ivica*)pom=
        ReplacePart[pom,{1,0,1},1];(*bela je 0;crna je 1*)(*redni broj,
        nula za boju,nula za da li je do sada obojena*)While[
        Union[Flatten[Map[Take[#,-1]&,pom]]]!={1},
        Do[If[pom[[i,3]]==1,
            pom=Map[If[#[[3]]!=1&&MemberQ[nbs[[i]],#[[1]]],{#[[1]],
                      1-pom[[i,2]],1},#]&,pom]],{i,1,Length[nbs]}]];
      pom=
        Select[pom,#[[2]]==1&];(*pom je lista koj asdarzi crne \
pljosni*)


       pom=Flatten[Map[Take[#,1]&,pom]];
         
      jj=Table[faces[[pom[[i]]]],{i,Length[pom]}];
      jj1=Map[fMakeCycle,jj];
      ttt1=ttt[[2]];
      pom=If[SameQ[
        Union[Flatten[
            Table[If[SameQ[Intersection[jj1[[i]],ttt1[[j]]],Union[ttt1[[j]]]],\
\
\
\
\

                pom[[i]],{}],{i,Length[jj1]},{j,Length[ttt1]}]]],pom],
      Complement[Range[Length[faces]],pom],pom];
       
    
          
      (*redni brojevi crnih pljosni*)
      
      faces=
        Map[Prepend[#,Position[faces,#][[1,1]]]&,faces];
      faces=Select[faces,MemberQ[pom,#[[1]]]&];
      Do[faces=ReplacePart[faces,Flatten[{i,Rest[faces[[i]]]}],i],{i,1,
          Length[faces]}];
      faces=Flatten[Map[fMakeCircEd,faces],1];
      (*lista parova {stranica,teme} a mi cemo sad dobiti {str,str}*)n=
        Length[Union[Take[Flatten[faces],{2,Length[Flatten[faces]],2}]]];
      (*ovo su u sustini smao sva temena*)pom={};
      Do[p=Select[faces,#[[2]]==i&];
        pos=Map[Position[faces,#][[1,1]]&,p];
        pom=Append[pom,{p[[1,1]],p[[2,1]]}];
        faces=Drop[faces,{pos[[1]]}];
        faces=Drop[faces,{pos[[2]]-1}],{i,1,n}];
    (*  ShowLabeledGraph[FromUnorderedPairs[pom]] *)
      ]];
   (*   Print["Graph of KL: ", pom]; *)
      mm=Union[Select[pom,Length[Position[pom,#]]>1 &]];
    (*  Print["Double edges: ", mm]; *)
      pom=Sort[pom];
      pom=If[IsomorphicQ[FromUnorderedPairs[pom],Cycle[Length[pom]]],
      Table[{1,2},{i,Length[pom]}],pom];
      pom]

      
      
  
  

(*fDual GraphKL:
    input Conway ili Dow,output:graf KL*)

fDualGraphKL[Ulaz_]:=
  
  Module[{Ul,ind,i,pe,faces,NewGr,pom,facesEd,pos,p,mm,pom1,nbs},
    If[Head[Ulaz]==List,Ul=Abs[Ulaz];ind=Ul[[1]]];
    If[Head[Ulaz]==String,Ul=Ulaz;ind=fDToD[Ulaz][[1]]];
    If[IsomorphicQ[FromUnorderedPairs[Union[fGraphInc[Ulaz][[1]]]],
        Cycle[Max[Flatten[Union[fGraphInc[Ulaz][[1]]]]]]],
      pom=Table[{1,2},{i,Length[fGraphKL[Ulaz]]}],
      If[SameQ[ind,{1,1}],(*slucaj 2,-2*)pom={{1,2},{1,2}},(*ostali*)pe=
          fPlanarEmb[Ul];
        NewGr=pe[[2]];
        pom=Union[NewGr];
        faces=pe[[5]];(*Print[faces];*)Do[
          If[Count[NewGr,pom[[i]]]==2,faces=Append[faces,pom[[i]]]],{i,
            1,Length[pom]}];
        (*fac su faces+digons*)(*sad treba da nadjemo susedne*)facesEd=
          Map[fMakeEd,faces];
        nbs=Map[fNeighFac[facesEd,#]&,facesEd];
        (*cikli susednih oko faca*)pom=
          Map[{#,0,0}&,Range[Length[nbs]]];(*brojevi svih ivica*)pom=
          ReplacePart[pom,{1,0,1},1];(*bela je 0;crna je 1*)(*redni broj,
          nula za boju,nula za da li je do sada obojena*)While[
          Union[Flatten[Map[Take[#,-1]&,pom]]]!={1},
          Do[If[pom[[i,3]]==1,
              
              pom=Map[If[#[[3]]!=1&&
                        MemberQ[nbs[[i]],#[[1]]],{#[[1]],\

                        1-pom[[i,2]],1},#]&,pom]],{i,1,Length[nbs]}]];
        pom=
          
          Select[pom,#[[2]]==1&];(*pom je lista koj asdarzi crne \

              pljosni*)pom=Flatten[Map[Take[#,1]&,pom]];
        (*redni brojevi crnih pljosni*)pom1=
          Complement[Range[Length[faces]],pom];(*bele faces*)pom=pom1;
        faces=Map[Prepend[#,Position[faces,#][[1,1]]]&,faces];
        faces=Select[faces,MemberQ[pom,#[[1]]]&];
        Do[
          faces=ReplacePart[faces,Flatten[{i,Rest[faces[[i]]]}],i],{i,1,
            Length[faces]}];
        faces=Flatten[Map[fMakeCircEd,faces],1];
        (*lista parova {stranica,teme} a mi cemo sad dobiti {str,str}*)n=
          Length[Union[Take[Flatten[faces],{2,Length[Flatten[faces]],2}]]];
        (*ovo su u sustini smao sva temena*)pom={};
        Do[p=Select[faces,#[[2]]==i&];
          pos=Map[Position[faces,#][[1,1]]&,p];
          pom=Append[pom,{p[[1,1]],p[[2,1]]}];
          faces=Drop[faces,{pos[[1]]}];
          faces=Drop[faces,{pos[[2]]-1}],{i,1,n}];
        (*ShowLabeledGraph[FromUnorderedPairs[pom]]*)]];
    (*Print[Graph of KL: ,pom];*)mm=
      Union[Select[pom,Length[Position[pom,#]]>1&]];
    (*Print[Double edges: ,mm];*)pom=Sort[pom];
    pom]
      

    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## #*)
(*## ## ## ## ## ## # fGenerators ## ## ## ## ## ## ## ## ## ## ## ##*)

fGenerators[Ulaz_]:=
  Module[{ges,t,ou,lpoz,pod,lneg,pro,qq,qq1,rr,ss,qq2,qq3,pro1,pro2,pr,
  ul,iz,gen, qq4,sg, i, j,Ul},  
   If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    ges=fGaussExtSigns[Ul];
    (*ges nam odredjuje znake*)
        t=If[SameQ[ges,Flatten[ges]],Length[ges],Map[Length,ges]];
    ou=If[
        SameQ[Head[Ul],
          String],-Map[Sign,Flatten[ges]]*-Map[Sign,
              Flatten[fGaussExtSigns[StringReplace[Ul,"-"->""]]]]*
          Flatten[Abs[ges]]*Table[(-1)^i,{i,Length[Flatten[ges]]}],
Map[Sign,Flatten[fGaussExtSigns[fKnotscapeDowfromDow[fDowfromPD[Ul]]]]]*Flatten[ges]*  Map[Sign,Flatten[ges]]*
     Table[(-1)^i,{i,Length[Flatten[ges]]}]]; 
  t=If[SameQ[Head[t],Integer],{t},t];
    ou=iteratedTake[ou,t];
    ou=If[ou[[1,1]]>0,ou,-ou];
    (*ou=over-under iz koje izdvajamo generatore*) 
    lpoz=Table[Select[ou[[i]],#>0&],{i,Length[ou]}];
    t=Map[Length,lpoz];
    lpoz=Flatten[lpoz];
    lpoz=Flatten[Table[Position[ou,lpoz[[j]]],{j,Length[lpoz]}],1];
    pod=iteratedTake[Table[ou[[lpoz[[i,1]],lpoz[[i,2]]]],{i,Length[lpoz]}],
        t];
    pod=Map[Length,pod];
    (*pod nam odredjuje raspodelu generatora po komponentama*)lneg=
      Table[Select[ou[[i]],#<0&],{i,Length[ou]}];
    pro=Flatten[
        Table[Append[
            Table[{lneg[[j,i]],lneg[[j,i+1]]},{i,Length[lneg[[j]]]-1}],{Last[
                lneg[[j]]],First[lneg[[j]]]}],{j,Length[lneg]}],1];
    (*pro je flatten lista negativnih uzastopnih u ou po komponentama*)qq=
      Table[Take[
          RotateLeft[ou[[Flatten[Position[ou,pro[[i,1]]]][[1]]]],
            Flatten[Position[ou,pro[[i,1]]]][[2]]],
          Flatten[Position[
                  RotateLeft[ou[[Flatten[Position[ou,pro[[i,1]]]][[1]]]],
                    Flatten[Position[ou,pro[[i,1]]]][[2]]],
                  pro[[i,2]]]][[1]]-1],{i,Length[pro]}];
    qq1=qq;
    rr=Complement[Union[Flatten[qq]],
        Flatten[If[
            SameQ[Union[
                  Table[If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                      Length[qq]}]][[1]],{}],
            Drop[Union[
                Table[If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                    Length[qq]}]],1],
            Union[Table[
                If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                  Length[qq]}]]]]];
    ss=Flatten[Position[qq,{}]];
    Do[qq=ReplacePart[qq,{rr[[i]]},ss[[i]]],{i,Length[rr]}];
    qq2=Table[qq[[i,1]],{i,Length[qq]}];
    qq3=Select[qq1,Not[SameQ[#,{}]]&];
    pr=Flatten[
        Table[Table[qq3[[j,1]],{i,Length[qq3[[j]]]}],{j,Length[qq3]}]];
    qq3=-Flatten[qq3];
    pro1=Table[pro[[i,1]],{i,Length[pro]}];
    pro2=Table[pro[[i,2]],{i,Length[pro]}];
    ul=Flatten[
        Table[Flatten[qq2[[Flatten[Position[pro2,qq3[[i]]]]]]],{i,
            Length[qq]}]];
    iz=Flatten[
        Table[Flatten[qq2[[Flatten[Position[pro1,qq3[[i]]]]]]],{i,
            Length[qq]}]];
    gen=Flatten[Table[{ul[[i]],iz[[i]],pr[[i]]},{i,Length[ul]}]];
    gen=iteratedTake[gen,3*pod];
    ges=Union[Flatten[ges]];
    qq4=Abs[qq3];
    sg=Flatten[
        Table[Sign[ges[[Flatten[Position[Abs[ges],qq4[[i]]]]]]],{i,
            Length[ges]}]];
    {gen,sg}]
    (*09.04.2012*)

fGeneratorsMirr[Ulaz_]:=
  Module[{ges,t,ou,lpoz,pod,lneg,pro,qq,qq1,rr,ss,qq2,qq3,pro1,pro2,pr,
  ul,iz,gen, qq4,sg, i, j,Ul},  
   If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    ges=fGaussExtSigns[Ul];
    (*ges nam odredjuje znake*)
        t=If[SameQ[ges,Flatten[ges]],Length[ges],Map[Length,ges]];
    ou=If[
        SameQ[Head[Ul],
          String],-Map[Sign,Flatten[ges]]*-Map[Sign,
              Flatten[fGaussExtSigns[StringReplace[Ul,"-"->""]]]]*
          Flatten[Abs[ges]]*Table[(-1)^i,{i,Length[Flatten[ges]]}],
Map[Sign,Flatten[fGaussExtSigns[fKnotscapeDowfromDow[fDowfromPD[Ul]]]]]*Flatten[ges]*  Map[Sign,Flatten[ges]]*
     Table[(-1)^i,{i,Length[Flatten[ges]]}]]; 
  t=If[SameQ[Head[t],Integer],{t},t];
    ou=iteratedTake[ou,t];
    ou=If[ou[[1,1]]>0,ou,-ou];

ou=-ou; 

    (*ou=over-under iz koje izdvajamo generatore*) 
    lpoz=Table[Select[ou[[i]],#>0&],{i,Length[ou]}];
    t=Map[Length,lpoz];
    lpoz=Flatten[lpoz];
    lpoz=Flatten[Table[Position[ou,lpoz[[j]]],{j,Length[lpoz]}],1];
    pod=iteratedTake[Table[ou[[lpoz[[i,1]],lpoz[[i,2]]]],{i,Length[lpoz]}],
        t];
    pod=Map[Length,pod];
    (*pod nam odredjuje raspodelu generatora po komponentama*)lneg=
      Table[Select[ou[[i]],#<0&],{i,Length[ou]}];
    pro=Flatten[
        Table[Append[
            Table[{lneg[[j,i]],lneg[[j,i+1]]},{i,Length[lneg[[j]]]-1}],{Last[
                lneg[[j]]],First[lneg[[j]]]}],{j,Length[lneg]}],1];
    (*pro je flatten lista negativnih uzastopnih u ou po komponentama*)qq=
      Table[Take[
          RotateLeft[ou[[Flatten[Position[ou,pro[[i,1]]]][[1]]]],
            Flatten[Position[ou,pro[[i,1]]]][[2]]],
          Flatten[Position[
                  RotateLeft[ou[[Flatten[Position[ou,pro[[i,1]]]][[1]]]],
                    Flatten[Position[ou,pro[[i,1]]]][[2]]],
                  pro[[i,2]]]][[1]]-1],{i,Length[pro]}];
    qq1=qq;
    rr=Complement[Union[Flatten[qq]],
        Flatten[If[
            SameQ[Union[
                  Table[If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                      Length[qq]}]][[1]],{}],
            Drop[Union[
                Table[If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                    Length[qq]}]],1],
            Union[Table[
                If[Not[SameQ[qq[[i]],{}]],{qq[[i,1]]},{}],{i,
                  Length[qq]}]]]]];
    ss=Flatten[Position[qq,{}]];
    Do[qq=ReplacePart[qq,{rr[[i]]},ss[[i]]],{i,Length[rr]}];
    qq2=Table[qq[[i,1]],{i,Length[qq]}];
    qq3=Select[qq1,Not[SameQ[#,{}]]&];
    pr=Flatten[
        Table[Table[qq3[[j,1]],{i,Length[qq3[[j]]]}],{j,Length[qq3]}]];
    qq3=-Flatten[qq3];
    pro1=Table[pro[[i,1]],{i,Length[pro]}];
    pro2=Table[pro[[i,2]],{i,Length[pro]}];
    ul=Flatten[
        Table[Flatten[qq2[[Flatten[Position[pro2,qq3[[i]]]]]]],{i,
            Length[qq]}]];
    iz=Flatten[
        Table[Flatten[qq2[[Flatten[Position[pro1,qq3[[i]]]]]]],{i,
            Length[qq]}]];
    gen=Flatten[Table[{ul[[i]],iz[[i]],pr[[i]]},{i,Length[ul]}]];
    gen=iteratedTake[gen,3*pod];
    ges=Union[Flatten[ges]];
    qq4=Abs[qq3];
    sg=Flatten[
        Table[Sign[ges[[Flatten[Position[Abs[ges],qq4[[i]]]]]]],{i,
            Length[ges]}]];
    {gen,sg}]
    (*09.04.2012*)

    
(*## ## ## ## ## ## # fColTest ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)

parr[n_,1]:={{n}}
parr[n_,r_]:=Module[{ans={},i},
    Do[
      AppendTo[ans,(Join[{i},#]&)/@parr[n-i,r-1]],{i,0,n}];
    Flatten[ans,1]]

fVarP3[n_Integer] := 
  Module[{r, i}, r = Table[IntegerDigits[i, 3, n], {i, 0, 3^n - 1}];
    r]

fDet[Ul_] :=   Module[{res},
Abs[ReplaceAll[RedAlex[Ul],x->-1]]]   


fColTest[Ul_,ppp_Integer] :=   Module[{pp,pp1,res,tt,tt1,tt2,tt3,tt4,ss,ff0,ff1},
pp=Flatten[fGenerators[Ul][[1]]];
pp1=Partition[ToExpression[Table[StringJoin["f",ToString[pp[[i]]]],{i,Length[pp]}]],3];
res=Reduce[Table[Mod[pp1[[ii,1]]+pp1[[ii,2]]-2*pp1[[ii,3]],ppp]==0,{ii,Length[pp1]}],Union[Flatten[pp1]],Integers];
tt=Table[Mod[Table[res[[2,j,i,2,1]],{i,Length[res[[2,2]]]}],ppp],{j,Length[res[[2]]]}];
tt=Union[Select[tt,Length[Union[#]]>1 &]];
tt=Table[ss=tt[[i]];ff0=Range[Length[Union[pp]]];ff1=Map[ToExpression,Sort[Map[ToString,ff0]]];Map[Last,Sort[Table[{ff1[[i]],ss[[i]]},{i,Length[ss]}]]],{i,Length[tt]}];
tt1=Union[Table[Length[Union[tt[[i]]]],{i,Length[tt]}]];
tt2=Table[First[Select[tt,SameQ[Length[Union[#]],tt1[[i]]] &] ],{i,Length[tt1]}];
tt3=Table[Table[{i,tt2[[j,i]]},{i,Length[tt2[[j]]]}],{j,Length[tt2]}];
tt4=Table[Partition[ReplaceAll[pp,Table[tt3[[j,i,1]]->tt3[[j,i,2]],{i,Length[tt3[[j]]]}]],3],{j,Length[tt3]}];
Table[{Partition[pp,3],tt4[[i]],tt3[[i]],tt1[[i]]},{i,Length[tt1]}]]


fColTestMirr[Ul_,ppp_Integer] :=   Module[{pp,pp1,res,tt,tt1,tt2,tt3,tt4,ss,ff0,ff1},
pp=Flatten[fGeneratorsMirr[Ul][[1]]];
pp1=Partition[ToExpression[Table[StringJoin["f",ToString[pp[[i]]]],{i,Length[pp]}]],3];
res=Reduce[Table[Mod[pp1[[ii,1]]+pp1[[ii,2]]-2*pp1[[ii,3]],ppp]==0,{ii,Length[pp1]}],Union[Flatten[pp1]],Integers];
tt=Table[Mod[Table[res[[2,j,i,2,1]],{i,Length[res[[2,2]]]}],ppp],{j,Length[res[[2]]]}];
tt=Union[Select[tt,Length[Union[#]]>1 &]];
tt=Table[ss=tt[[i]];ff0=Range[Length[Union[pp]]];ff1=Map[ToExpression,Sort[Map[ToString,ff0]]];Map[Last,Sort[Table[{ff1[[i]],ss[[i]]},{i,Length[ss]}]]],{i,Length[tt]}];
tt1=Union[Table[Length[Union[tt[[i]]]],{i,Length[tt]}]];
tt2=Table[First[Select[tt,SameQ[Length[Union[#]],tt1[[i]]] &] ],{i,Length[tt1]}];
tt3=Table[Table[{i,tt2[[j,i]]},{i,Length[tt2[[j]]]}],{j,Length[tt2]}];
tt4=Table[Partition[ReplaceAll[pp,Table[tt3[[j,i,1]]->tt3[[j,i,2]],{i,Length[tt3[[j]]]}]],3],{j,Length[tt3]}];
Table[{Partition[pp,3],tt4[[i]],tt3[[i]],tt1[[i]]},{i,Length[tt1]}]]


fColNo[Ul_,pp_Integer] := 
  Module[{res}, Map[Last,fColTest[Ul,pp]]]


fColNoMirr[Ul_,pp_Integer] := 
  Module[{res}, Map[Last,fColTestMirr[Ul,pp]]]
    
fSchubertBridges[Ul_]:=Module[{cc,a,b,pp,pp1,ss,ss1,tt1,tt2,tt,dd,dd1,ttt,ttt1,ttt2,ttt3,ttt4,ttt5,ttt6,pd},
cc=fComponentNo[Ul];
a=fDet[Ul];
b=fDet[StringJoin[Ul," (1,-1)"]];
pp=Table[{i,2*a-i},{i,a-1}];
pp1=Drop[Drop[Mod[Range[0,a*b,b],2*a],1],-1];
ss=Table[i,{i,a-1}];
ss1=ss+a;
tt1=Table[{pp1[[i]],If[EvenQ[i],1,2]},{i,Length[pp1]}];
tt2=Table[{pp1[[i]],If[EvenQ[i],2,1]},{i,Length[pp1]}];
tt=Join[tt1,tt2];
ttt=Partition[Table[dd=tt[[i]];dd1=Flatten[Position[pp,dd[[1]]]];pp[[dd1[[1]]]][[dd[[2]]]],{i,Length[tt]}],a-1];
ttt=If[EvenQ[b],{ttt[[1]],Mod[ttt[[1]]+a,2a]},ttt];
ttt1=If[SameQ[cc,1],{ttt[[1]],If[EvenQ[b],Reverse[-ss1],-ss1],ttt[[2]],-Reverse[ss]},{ttt[[1]],If[EvenQ[b],-ss,Reverse[-ss]],Reverse[ttt[[2]]],Reverse[-ss1]}];
ttt2=Map[Sign,Flatten[ttt1]];
ttt3=Abs[Flatten[ttt1]];
ttt4=Table[Flatten[Position[ttt3,ttt3[[i]]]],{i,Length[ttt3]}];
ttt5=Union[ttt4];
ttt6=ReplaceAll[ttt4,Table[ttt5[[i]]->i,{i,Length[ttt5]}]]*ttt2;
ttt6=If[SameQ[cc,1],ttt6,Partition[ttt6,2*a-2]];
pd=fPDToPData[PD[GaussCode@@ttt6]]] 


(*## ## ## ## ## ## # fPrimeKL ## ## ## ## ## ## ## ## ## ## ## ## ##*)

(*fPrimeKL checks that KL given by pData,
  Dowker code or Conway symbol is a prime or
  composite*)
(*vraca 1 ako je Prime 0 za direktan *)

fPrimeKL[Ulaz_]:=Module[{pd,Dow,gr, i,Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
       If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    If[SameQ[Union[Map[EvenQ,Abs[pd[[2]]]]],{True}],
    pd=fPDataFromDow[Map[Abs,pd]]];
    (*sada sigurno imamo pdata *)
    If[SameQ[pd[[1]],{1,1}],gr=1,
    pd=ReductionKnotLink[pd];
    Dow=fDowfromPD[pd];
    If[SameQ[Dow,{}],gr=0,
    gr=fGaussExt[Dow];
    gr=fGrInc[gr][[1]];
    If[SameQ[Flatten[Union[Map[Bridges,
              Map[FromUnorderedPairs,
                Union[Table[Drop[gr,{i}],{i,Length[gr]}]]]]]],{}],
           gr=1, gr=0]]];
    gr]
    
    (* korigovao S.J. 2.02.2007 *) 
    
    
    (*## ## ## ## ## ## # fPrimeGraph ## ## ## ## ## ## ## ## ## ## ## ## \
##*)
    (*fPrimeKLGraph checks that KL given by graph is a 
      prime or composite KL= direct product*)
	(*vraca 1 ako je Prime 0 za direktan *)

	fPrimeGraph[GG_List]:=Module[{gr, i},
   	 If[SameQ[Flatten[Union[
            Map[Bridges,
            Map[FromUnorderedPairs,
            Union[Table[Drop[GG,{i}],{i,Length[GG]}]]]]]],{}],
             gr=1,gr=0];
    gr]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
##*)
(* ## ## ## ## ## ## # Murasugi Sig ## ## ## ## ## ## ## ## ## ## ## ## ## *)

(* fSigTor calculates Murasugi 
signature for torus KL - Murasugi, pp .148 *)
fSigTor[m_Integer,n_Integer]:=Module[{a,b,r,p,s},
    p=1;
    s=0;
    r={p,{m,n},s};
    a=Max[m,n];
    b=Min[m,n];
    While[Not[SameQ[r[[2]],{0,0}]],
      p=r[[1]];
      a=Max[r[[2,1]],r[[2,2]]];
      b=Min[r[[2,1]],r[[2,2]]];
      s=r[[3]];
      r={p,{a,b},s};
      r=If[SameQ[Min[a,b],1],{p,{0,0},s},r];
      r=If[SameQ[Min[a,b],2],{p,{0,0},s+p*(Max[a,b]-1)},r];
      r=If[
          SameQ[Max[a,b],2*Min[a,b]]&&Not[SameQ[Min[a,b],1]]&&
            Not[SameQ[Min[a,b],2]],{{p},{0,0},s+p*(Min[a,b]^2-1)},r];
      r=If[
          Max[a,b]>2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],1],{p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*(Min[a,b]^2-1)},r];
      r=If[
          Max[a,b]>2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],0],{p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*Min[a,b]^2},r];
      r=If[
          Max[a,b]<2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],1],{-p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*(Min[a,b]^2-1)},r];
      r=If[
          Max[a,b]<2*Min[a,b]&&Not[SameQ[Max[a,b],2*Min[a,b]]]&&
            Not[SameQ[Min[a,b],1]]&&Not[SameQ[Min[a,b],2]]&&
            SameQ[Mod[Min[a,b],2],0],{-p,Abs[{Max[a,b]-2*Min[a,b],Min[a,b]}],
            s+p*(Min[a,b]^2-2)},r]];
     r[[3]]
    ]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## ## #*)
(*## ## ## ## ## ## ## ## ## ## ##*)



(*## ## ## ## ## ## # fPDataFromDow ## ## ## ## ## ## ## ## ## ## ## ## ##*)

   
 fPDataFromDow[Ul_List]:=Module[{rr,rr1,rr2,i},
    rr=fSignsKL[Abs[Ul]];
    rr=Map[Sign,Ul[[2]]]*Map[Sign,rr[[2]]];
    rr1=Table[{2i-1,Abs[Ul[[2,i]]]},{i,Length[Ul[[2]]]}];
    rr2=Map[Sign,Ul[[2]]];
    rr1=Table[
        If[SameQ[rr2[[i]],1],rr1[[i]],Reverse[rr1[[i]]]],{i,Length[rr1]}];
    rr1=Map[Last,
        Sort[Table[{rr1[[i,1]],rr1[[i,2]]*rr[[i]]},{i,Length[rr]}]]];
    rr={Ul[[1]],rr1};
    rr] 
 
(* KORIGOVAO SLAVIK 25.05.2007 *)

   
    (*## ## ## ## # fPDataFromDowker ## ## ## ##*)

(*Calculates pdata from Dowker code with signs*)

fPDataFromDowker[Dow_List] := Module[{dd1, dd2, cs, i},
    dd1 = Map[Sign, Dow];(*znaci cvora*)
    dd2 = Map[Sign, fSignsKL[Abs[Dow]]];(*znaci alternirajuceg*)
    cs = dd1*dd2;(*mesta gde se razlikuju*)
    dd2 = Table[{2i - 1, Abs[Dow[[2, i]]]}, {i, Length[Dow[[2]]]}];
    dd2 = 
      Table[If[cs[[2, i]] == 1, dd2[[i]], Reverse[dd2[[i]]]], {i, 
          Length[dd2]}];
    dd2 = 
      Sort[Table[{dd2[[i, 1]], dd1[[2, i]]*dd2[[i, 2]]}, {i, Length[dd2]}]];
    dd2 = {Dow[[1]], Table[dd2[[i, 2]], {i, Length[dd2]}]};
    dd2]   


(* ::Input:: *)
(**)


(*## ## ## ## ## ## # fComponentNo ## ## ## ## ## ## ## ## ## ## ## ## ##*)

(*fComponentNo for a knot given by Conway, Dow, 
  or pdata gives the number of components *)
  
fComponentNo[Ul_]:=Module[{pd,l},
    pd=If[ SameQ[Head[Ul],Symbol],ToCharacterCode[StringTake[ToString[Ul],{2}]][[1]]-96,
If[SameQ[Head[Ul],String],Length[fCreatePData[Ul][[1]]],Length[Ul[[1]]]]]
    ]

(*## ## ## ## ## ## # RATIONAL KL ## ## ## ## ## ## ## ## ## ## ## ## ##*)

(* compo pravi particije broja n *)
compo[0] := {{}}
compo[n_Integer] := compo[n] = Module[{i},
		Flatten[Table[(Join[{i}, #1] &) /@ compo[n - i], {i, n}], 1]]

(* vrsi racionalnu redukciju nealternirajuceg racionalnog *)

(* RationalKL computes number and list of Conway symbols
of KL with n crossings *)

RationalKL[n_Integer]:=Module[{pa,t,u,i},
    pa=compo[n];
    t=Union[
        Table[If[First[pa[[i]]]==1||Last[pa[[i]]]==1,{},
            pa[[i]]],{i,Length[pa]}]];
    t=If[t[[1]]=={},Drop[t,1],t];
    t=Union[
        Table[If[MemberQ[t,Reverse[t[[i]]]],
            Last[Sort[{t[[i]],Reverse[t[[i]]]}]],{}],{i,Length[t]}]];
    t=If[t[[1]]=={},Drop[t,1],t];
    t=StringReplace[
        Map[ToString,t],{"{"->"","}"->"",","->""}];
    t=Join[t,{Length[t]}];
    t]

(*## ## ## ## ## ## # RationalAmphiK ## ## ## ## ## ## ## ## ## ## ## ## ##*)

(*RationalAmphiK calculates number and list of 
Conway symbols of rational amphicheiral knots *)

RationalAmphiK[n_Integer]:=
  Module[{pa,t,u,i,pom},If[SameQ[EvenQ[n],True],pa=compo[n/2];
      t=Union[Table[If[First[pa[[i]]]==1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Table[Join[t[[i]],Reverse[t[[i]]]],{i,Length[t]}];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
      pom=Map[Dowker,t];
      pom=Table[Head[pom[[i]][[4]][[1]]],{i,Length[pom]}];
      pom=
        Sort[Union[
            Table[If[SameQ[pom[[i]],Integer]==True,t[[i]],{}],{i,
                Length[pom]}]]];
      t=Reverse[If[SameQ[pom[[-1]],{}]==True,Drop[pom,-1],pom]];
      t=Join[t,{Length[t]}];
      t,0]]
 (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #*)
 
 (*## ## ## ## ## ## # RationalAmphiK ## ## ## ## ## ## ## ## ## ## ## ## \
##*)

(*RationalAmphiL calculates number and list of 
Conway symbols of rational amphicheiral links *)
 
 RationalAmphiL[n_Integer]:=
  Module[{pa,t,t1,u,i,pom},If[SameQ[EvenQ[n],True],pa=compo[n/2];
      t=Union[Table[If[First[pa[[i]]]==1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Table[Join[t[[i]],Reverse[t[[i]]]],{i,Length[t]}];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
      pom=Map[Dowker,t];
      pom=Table[Head[pom[[i]][[4]][[1]]],{i,Length[pom]}];
      pom=
        Sort[Union[
            Table[If[SameQ[pom[[i]],Integer]==True,t[[i]],{}],{i,
                Length[pom]}]]];
      t=Reverse[If[SameQ[pom[[-1]],{}]==True,Drop[pom,-1],pom]];
      t1=Union[Table[If[First[pa[[i]]]==1,{},pa[[i]]],{i,Length[pa]}]];
      t1=If[First[t1]=={},Drop[t1,1],t1];
      t1=Table[Join[t1[[i]],Reverse[t1[[i]]]],{i,Length[t1]}];
      t1=StringReplace[
          Map[ToString,t1],{"{"->"","}"->"",","->""}];
      t=Complement[t1,t];
      t,0];
       t=Join[t,{Length[t]}];
       t]
 
 (*## ## ## ## ## # RacGenKL ## ## ## ## ## ## #*)      

(* restrictedPartitions[m_Integer, r_Integer, 
      n_Integer] pravi ogranicene particije broja n, 
  na najvise r delova od kojih je svaki najvise m *)
  
partH[_, 0] := {{}}
partH[H_, n_] := {} /; First[H] > n
partH[H_, n_] := Module[{h, p, ans}, ans = Map[(h = #;
            p = partH[Select[H, # >= h &], n - h];
            Map[Join[{h}, #] &, p]) &, H];
    Flatten[ans, 1]]
partit[n_] := partH[Range[n], n]
partitionsOdd[n_] := partH[Range[1, n, 2], n]
partitionsDif[n_] := DeleteCases[partit[n], {___, x_, x_, ___}]

restrictedPartitions[m_Integer, r_Integer, n_Integer] := 
  Select[partH[Range[m], n], (Length[#] <= r) &]

(*## ## ## ## ## ## ## ## ## ## ## ##*)

(* RacGenKL[n_Integer, m_Integer] su racionalni generating  za m = 
    3 i sourcelinks KL za m = 
      2 : racionalni generisuci KL ciji se zapis sastoji samo od 1, 2, 3, 
  koji direktno daju familije i racionalni source links m = 2 *) 

RatGenSourKL[n_Integer, m_Integer] := Module[{pa, t, u, i},
    pa = Flatten[Map[DistinctPermutations, 
    restrictedPartitions[m, n - 1, n]],1];
    t = Union[
        Table[If[First[pa[[i]]] == 1 || 
        Last[pa[[i]]] == 1, {}, pa[[i]]], {i, Length[pa]}]];
    t = Union[
        Table[If[MemberQ[t, Reverse[t[[i]]]], 
            Last[Sort[{t[[i]], Reverse[t[[i]]]}]], {}], {i, Length[t]}]];
    t = If[First[t] == {}, Drop[t, 1], t];
    t = StringReplace[Map[ToString, t], {"{" -> "", "}" -> "", "," -> ""}];
     t=Join[t,{Length[t]}];
     t 
    ]
    
    (* RatSourceKLNo calculates the number of 
    rational source KL for k<=n crossings*)
RatSourceKLNo[n_Integer]:=Module[{b={1},i,p},
    If[n<4,1,
           If[n==4||n==5,
                  b={1},
                  b={1,1};
                 Do[
                       p=b[[i-1]]+b[[i-2]];
          	If[EvenQ[i],p=p-Fibonacci[(i-2)/2] ];
                      b=Append[b,p]
           ,{i,3,n-3}]
            ]
       ];
    Last[b]
    ]
 
(*## ## ## ## ## ## ## ## ## ## ## ##*)

(* RatLinkU1 calculates number and Conway symbols of rational 
links with unlinking number 1 according to P.Kohn*)

RatLinkU1[n_Integer]:=Module[{pa,t,u,i},
    If[SameQ[OddQ[n],True],
      pa=compo[(n+1)/2];
      t=Union[
          Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
      t=Union[
          Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
              Join[t[[i]],u[[i]]],{}],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Union[
          Table[If[First[t[[i]]]==1||Last[t[[i]]]==1,{},
              t[[i]]],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
      (* Print[Length[t]]; *)
      If[n==5,{"2 1 2"},Flatten[{t,Length[t]}]],0]
      ]

(*## ## ## ## ## ## ## ## ## ## ## ##*)

(*RatLinkU0[n_Integer]:n-odd =rational unlinks U=0*)

RatLinkU0[n_Integer]:=Module[{pa,t,u,i},
    If[n==5,{"2 -1 2",1},
    If[SameQ[OddQ[n],True],
      pa=compo[(n+1)/2];
      t=Union[
          Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
      t=If[First[t]=={},Drop[t,1],t];
      u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
      t=Union[
          Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
              Join[ReplacePart[t[[i]],-1,-1],u[[i]]],{}],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=Union[
          Table[If[First[t[[i]]]==1||Last[t[[i]]]==1,{},
              t[[i]]],{i,Length[t]}]];
      t=If[First[t]=={},Drop[t,1],t];
      t=StringReplace[
          Map[ToString,t],{"{"->"","}"->"",","->""}];
     (*  Print[Length[t]]; *)
      Flatten[{t,Length[t]}],0]]      
      ]

(*## ## ## ## ## ## ## ## ## ## ## ##*)

(*MSigRat calculates Murasugi signature of a rational 
KL given by Conway symbol *)

MSigRat[Conway_String] := 
  Module[{Con, Con1,sig,q, p, l, i}, 
    Con = ReadList[
        StringToStream[
          StringReplace[
            Conway, {"," -> " ", "+" -> " ", "(" -> " ", ")" -> " "}]], 
        Number];
    Con = FromContinuedFraction[Con];
    Con = ContinuedFraction[Con];
    (* Con = Last[Sort[{Con, Reverse[Con]}]]; *)
    Con1 = FromContinuedFraction[Con];
   (*  Print[Con1]; *)
    q = Denominator[Con1];
    p = Numerator[Con1];
    l = Table[(i - 1)*q, {i, p}];
    l = Table[
        If[Mod[l[[i]], 2p] > p, Mod[l[[i]], 2p] - 2p, Mod[l[[i]], 2p] ], {i, 
          Length[l]}];
    sig = Apply[Plus, Map[Sign, l]];
    {Con1,sig}
    ]
    
    (*## ## ## ## ## ## ## ## ## ## ## ##*)
    
   (*RatKnotGenU1 generates rational knots with U=1 *)
   
fGenKnotU1[n_Integer]:=Module[{pa,t,u,i},
    pa=compo[(n+1)/2];
    t=Union[Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
    t=If[First[t]=={},Drop[t,1],t];
    u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
    t=Union[
        Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
            Join[t[[i]],u[[i]]],{}],{i,Length[t]}]];
    t=If[First[t]=={},Drop[t,1],t]
    ]


RatKnotU1[k_Integer,UL_List]:=Module[{pa,t,u,i},UL;
    t=Union[UL,Reverse/@UL];
    t=Union[
        Table[If[k==1,ReplacePart[t[[i]],t[[i]][[1]]+1,1],
            Join[{k},t[[i]]]],{i,Length[t]}]];
    t=Table[
        If[t[[i]][[-1]]==1,
          Drop[ReplacePart[t[[i]],t[[i]][[-2]]+1,-2],-1],t[[i]]],{i,
          Length[t]}];
    t=StringReplace[
        Map[ToString,t],{"{"->"","}"->"",","->""}]
       ]

RatKnotGenU1[n_Integer]:=Module[{i=5,res={}},
    If[n<3,Print["To small No !"],
         If[n==3,Print["{3}"],
          While[i<=2IntegerPart[n/2],
                        
            res=Append[res,Map[RatKnotU1[n-i,{#}]&,fGenKnotU1[i]]];
                      i=i+2];
           res=Prepend[res,{{ToString[n-2]<>" 2"}}];
                    ]];
                    Flatten[res];
      res=Join[res,{Length[res]}];
      res=Flatten[res];
      res
    ]
    
    (*## ## ## ## ## ## ## ## ## ## ## ##*)
    
    (*RatKnotGenU0[n_Integer] generates rational 
    unknots with U=0. *)

fGenKnotU0[n_Integer]:=Module[{pa,t,u,i},pa=compo[(n+1)/2];
    t=Union[Table[If[Last[pa[[i]]]!=1,{},pa[[i]]],{i,Length[pa]}]];
    t=If[First[t]=={},Drop[t,1],t];
    u=Table[Reverse[ReplacePart[t[[i]],t[[i]][[-2]]-1,-2]],{i,Length[t]}];
    t=Union[
        Table[If[MemberQ[Join[t[[i]],u[[i]]],0]==False,
            Join[ReplacePart[t[[i]],-1,-1],u[[i]]],{}],{i,Length[t]}]];
    t=If[First[t]=={},Drop[t,1],t]
   ]


RatKnotU0[k_Integer,UL_List]:=Module[{pa,t,u,i},UL;
    t=Union[UL,Reverse/@UL];
    t=Union[
        Table[If[k==1,ReplacePart[t[[i]],t[[i]][[1]]+1,1],
            Join[{k},t[[i]]]],{i,Length[t]}]];
    t=Table[
        If[t[[i]][[-1]]==1,
          Drop[ReplacePart[t[[i]],t[[i]][[-2]]+1,-2],-1],t[[i]]],{i,
          Length[t]}];
    t=StringReplace[
        Map[ToString,t],{"{"->"","}"->"",","->""}]]

RatKnotGenU0[n_Integer]:=Module[{i=3,res={}},
    If[n<5,Print["To small n!"],
        While[i<=2IntegerPart[n/2],
                      res=Append[res,Map[RatKnotU0[n-i,{#}]&,fGenKnotU0[i]]];
                    i=i+2];
                ];
             res= Flatten[res];
        res=Join[res,{Length[res]}];
        res
   
      ]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*## ## ## ## ## ALL STATES RATIONAL ## ## ## ## ## ## ## #*)
(*Pomocna za fChange ConVar, AllStatesRat *)
(*Radi od Konvej Liste i \
odgovarajuce varijacije- 
    elemment varijacije je broj koliko puta treba izvrsiti promenu u toj \
tacki*)
fReduCon[Con_List,VVvar_List]:=
  Module[{p,i,l=Con,rr},Do[p=Con[[i]]-2*Sign[Con[[i]]]*VVvar[[i]];
      l=ReplacePart[l,p,i],{i,1,Length[VVvar]}];
    IsNotKnot[l]]
(*menja vrednosti u Konvej listi prema varijaciji*)

fChangeConVar[Con_List,vvSum_List]:=Module[{l=Con},
         l=fReduCon[l,vvSum];
         l=StringReplace[ToString[l],"{}"->"{1}"];
         l=StringReplace[l,{"{"->"",","->"","}"->"", "0"->"1"}];
        {Apply[Plus,vvSum],"Con: "<>l}
    ]
(*Ulaz je Konvej ili PData *)
(*Radi sve pojedinacne promene znaka na fix projekciji*)

AllStatesRational[Conway_String]:=
  Module[{Con1,Con,ConList,vv,i,res={},res1},
    ConList=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con1=FromContinuedFraction[ConList];
    If[SameQ[Con1,ComplexInfinity],
      Con=FromContinuedFraction[Reverse[ConList]],Con=Con1];
    Con=ContinuedFraction[Con];
    (*nasli smo ga ko je zaista*)(*sad ga razvezujemo*)
    
    vv=fVarP[Apply[Plus,Map[Abs[#]&,Con]]];
    vv=Map[iteratedTake[#,Abs[Con]]&,vv];
    Do[vv=ReplacePart[vv,Map[Apply[Plus,#]&,vv[[i]]],i],{i,1,Length[vv]}];
    vv=Union[vv];
    (*sad je svaka varijacija podeljena prema duzinama clanova konveja*)
    
    vv=Union[Map[fChangeConVar[Con,#]&,vv]];
    Con=Sort[Map[Reverse[#]&,vv]];
    vv=Union[Map[Rest[#][[1]]&,vv]];
    Do[res=Append[res,Reverse[First[Select[Con,SameQ[#[[1]],vv[[i]]]&]]]],{i,
        1,Length[vv]}];
    res=Sort[res];
   res1=Table[{res[[i,1]],StringReplace[res[[i,2]],"-"->""]},
   {i,Length[res]}];
    res={Select[res1,#[[2]]=="Con: 1"&][[1,1]],res};
    res=Join[res,{Length[res[[2]]]}];
    res]
    
    AllStPR[Ul_,k_Integer]:=Module[{prk},
    prk=AllStatesRational[Ul][[2,k,2]];
    prk=StringDrop[prk,5];
    prk
    ]

(* ## ## ## ## ## ## ## ## # AllStatesRatFast ## ## ## ## ## ## ## ## ## #*)

AllStatesRatFast[Conway_String]:=
  Module[{Con1,Con,ConList,vv,i,res={},res1},
    ConList=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con1=FromContinuedFraction[ConList];
    If[SameQ[Con1,ComplexInfinity],
      Con=FromContinuedFraction[Reverse[ConList]],Con=Con1];
    Con=ContinuedFraction[Con];
    (*nasli smo ga ko je zaista*)(*sad ga razvezujemo*)
    
    vv=fVarNewP[Apply[Plus,Map[Abs[#]&,Con]]];
    vv=Map[iteratedTake[#,Abs[Con]]&,vv];
    Do[vv=ReplacePart[vv,Map[Apply[Plus,#]&,vv[[i]]],i],{i,1,Length[vv]}];
    vv=Union[vv];
    (*sad je svaka varijacija podeljena prema duzinama clanova konveja*)
    
    vv=Union[Map[fChangeConVar[Con,#]&,vv]];
    Con=Sort[Map[Reverse[#]&,vv]];
    vv=Union[Map[Rest[#][[1]]&,vv]];
    Do[res=Append[res,Reverse[First[Select[Con,SameQ[#[[1]],vv[[i]]]&]]]],{i,
        1,Length[vv]}];
    res=Sort[res];
     res1=Table[{res[[i,1]],StringReplace[res[[i,2]],"-"->""]},
     {i,Length[res]}
];
    {Select[res1,#[[2]]=="Con: 1"&][[1,1]],res}]


(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 25.09.2003*)



fVarNewPGap[n_Integer,u_Integer]:=
  Module[{r,i},r=Drop[Table[IntegerDigits[i,2,n],{i,0,2^n-1}],1];
    r=Sort[Table[{Count[r[[i]],1],r[[i]]},{i,Length[r]}]];
    r=Drop[
        Union[Table[
            If[r[[i,1]]>IntegerPart[(n+1)/2],{},r[[i,2]]],{i,Length[r]}]],1];
    r=Sort[Table[{Count[r[[i]],1],r[[i]]},{i,Length[r]}]];
    r=Select[r,#[[1]]==u&];
    r=Table[r[[i,2]],{i,Length[r]}];
    r]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ## ## ## ## ## ## # ALLSTATES ## ## ## ## ## ## ## ## ## ## ## ##*)
(*menja znake u pdata prema varijaciji*)
fCrChVarPD[PD_List,var_List]:=
  Module[{l,n,i,pos},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    pos=Flatten[Position[var,1]];
    Do[l=ReplacePart[l,-1*Reverse[l[[pos[[i]]]]],pos[[i]]],{i,1,
        Length[pos]}];
    l=Sort[l,Abs[#1[[1]]]<Abs[#2[[1]]]&];
    l=Flatten[Map[Take[#,-1]&,l]];
    l={PD[[1]],l};
    l=ReductionKnotLink[l];
    (*{Count[var,1],var,l}*)
    
    If[MemberQ[{{{},{}},{{0},{}},{{},{0}}},l],l="pd: {{},{}}",
      l="pd: "<>ToString[l]];
    {Count[var,1],l}]


(*Ulaz je Konvej ili PData*)
(*Radi sve pojedinacne promene znaka na fix \
projekciji*)
fAllStatesProj[UL_]:=Module[{pd=UL,i,prvi,vv,res={}},
    If[SameQ[Head[UL],String]&&SameQ[StringPosition[UL,"{"],{}],pd=\
fCreatePData[UL],pd=ToExpression[UL]];
    (*sad su pd pdata*)
    pd=ReductionKnotLink[pd];
    vv=Rest[fVarP[Length[pd[[2]]]]];
    pd=Union[Map[fCrChVarPD[pd,#]&,vv ]];
    prvi=Union[Map[Rest[#][[1]]&,pd]];
    pd=Sort[Map[Reverse[#]&,pd]];
    (*vv=Union[Map[Apply[Plus,#]&,vv]];Print[vv];*)
    Do[res=Append[res,
          Reverse[First[Select[pd,#[[1]]==prvi[[i]]&]]]
          ],
      
      {i,1,Length[prvi]}];
    res=Reverse[Sort[res]];
    res={Select[res,#[[2]]=="pd: {{},{}}"&][[1,1]],Sort[res]};
    res=Join[res,{Length[res[[2]]]}]]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)

fAllStPrk[Ul_,k_Integer]:=Module[{prk},
    prk=fAllStatesProj[Ul][[2,k,2]];
    prk=ToExpression[StringDrop[prk,4]];
    prk
    ]


(*menja znake u pdata prema varijaciji*)
fSignVarPD[PD_List,var_List]:=
  Module[{l,n,i,pos},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    pos=Flatten[Position[var,1]];
    Do[l=ReplacePart[l,-1*Reverse[l[[pos[[i]]]]],pos[[i]]],{i,1,
        Length[pos]}];
    l=Sort[l,Abs[#1[[1]]]<Abs[#2[[1]]]&];
    l=Flatten[Map[Take[#,-1]&,l]];
    l={PD[[1]],l};
    l=ReductionKnotLink[l];
    (*{Count[var,1],var,l}*)
    If[MemberQ[{{{},{}},{{0},{}},{{},{0}}},l],l=1,
      l=0]]

(*## ## ## ## ## ## ## ## # PERIOD KL ## ## ## ## ## ## ## ## ## ##*)

PeriodProjAltKL[Ulaz_]:=Module[{G,sG,Aut,k,Aut1, i,j,Ul},
  If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
If[SameQ[Ul,"1"],Aut={},
      G=fGaussExtSigns[Ul];
      sG=fGrInc[G][[2]];
      G=FromUnorderedPairs[fGrInc[G][[1]]];
      Aut=Drop[Automorphisms[G],1];
      Aut1=Aut;
      Aut=
        If[Aut!={},
          If[First[
                Union[Table[
                    If[SameQ[
                          Sign[Table[
                                Table[
                                  Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                    Length[sG]}],{i,Length[Aut]}]][[k]],
                          sG]==True,Aut[[k]],{}],{k,
                      Length[Aut]}]]]=={},
            Drop[Union[
                Table[If[
                    SameQ[Sign[
                            Table[Table[
                                Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                  Length[sG]}],{i,Length[Aut]}]][[k]],
                        sG]==True,Aut[[k]],{}],{k,Length[Aut]}]],1],
            Union[Table[
                If[SameQ[
                      Sign[Table[
                            Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                Length[sG]}],{i,Length[Aut]}]][[k]],
                      sG]==True,Aut[[k]],{}],{k,Length[Aut]}]]],Aut];
      (*trazimo cikle*)
      Aut=If[Aut!={},Map[ToCycles,Aut],Aut];
      Aut=Select[Aut,Count[Map[Length,#],1]<3&];
      Aut=Select[Aut,Length[Complement[Union[Map[Length,#]],{1}]]==1&];
      (*iz nekog razloga duzine cikla su periodi*)
      
      Aut=Union[Map[Length,Flatten[Aut,1]]];
      (*ako nije prazna izbaci prvu 1-identitet*)
      
      Aut=If[Aut!={}&&First[Aut]==1,Drop[Aut,1],Aut]];
    Aut=If[SameQ[Aut,{}],{0},Aut];
    Aut=If[Length[Aut1]==1,{2},Aut];
    Aut
    ]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## #*)
(*Radi od Konveja, dowkera i pdata *)
PeriodAltKL[Conway_String]:=Module[{res,razne},If[SameQ[Conway,"1"],
       res={0},
      razne=fDiffProjectionsAltKL[Conway];
      razne=Map[#[[1]]&,razne];
      res=Union[Flatten[Map[PeriodProjAltKL[#]&,razne]]];
      If[SameQ[res,{}],res={0}]];
      If[SameQ[res[[1]],0],res=Drop[res,1],res];
      res
    ]
    (*popravljano 16.8.2003 *)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## #*)
 (*## ## ## ## ## ## ## # MINDOWALT ## ## ## ## ## ## ## ## ## ## ## ## *) 
 
 (*## ## ## ## ## ## ## # MINDOWALT ## ## ## ## ## ## ## ## ## ## ## \
##*)(*fOrderComp:
    Input Con_String priprema link za pletenje.
    Kao izlaz daje Gauss code kod 
koga su sve komponente sa samopresecima 
dovedene na optimalan pocetak.Drugi
clan izlaza je Ulc kod koji pokazuje 
koje su komponente neodlucive Ne radi za
"2"*)(*fOrderComp:
    Input Con_String priprema link za pletenje.
    Kao izlaz daje Gauss code kod
koga su sve komponente sa samopresecima 
dovedene na optimalan pocetak.Drugi
clan izlaza je Ulc kod koji pokazuje 
koje su komponente neodlucive Ne radi za
"2"*)

fOrderComp[Con_]:=Module[{c,p,EDL,fmc,tc, i, Ulc\
j},EDL=Sort[fGaussExtSigns[Con]];    
    c=Table[Length[EDL[[i]]],{i,Length[EDL]}];    
    p=Table[
        Table[Length[Intersection[Abs[EDL[[j]]],Abs[EDL[[i]]]]],{i,
            Length[EDL]}],{j,Length[EDL]}];
    p=Map[Sort,
        Table[Table[If[p[[j,i]]>0,c[[i]],c[[j]]],{i,Length[p[[j]]]}],{j,
            Length[p]}]];(*c[[j]] ili 2*c[[j]]*)
    fmc=Map[fMinComp,EDL];
   fmc=Sort[Table[{c[[i]],p[[i]],fmc[[i]][[1]],EDL[[i]]},{i,Length[EDL]}]];
    tc=fmc;
    fmc=Table[{fmc[[i]][[1]],fmc[[i]][[2]],Sort[fmc[[i]][[3,1]]]},{i,
          Length[fmc]}];
    fmc=Split[fmc];
    fmc=Table[Length[fmc[[i]]],{i,Length[fmc]}];
    Ulc=Range[fmc]+Prepend[Table[Sum[fmc[[i]],{i,i}],{i,Length[fmc]-1}],0];
    (*ne valja rotacija*)tc=
      Table[If[SameQ[tc[[i]][[3,3]],0],
          RotateLeft[tc[[i]][[4]],tc[[i]][[3,2]]-1],
          Reverse[RotateRight[tc[[i]][[4]],
              Length[tc[[i]][[4]]]-tc[[i]][[3,2]]]]],{i,Length[tc]}];
    {tc,Ulc}]


(*fEqChoices[
      Ulc_List] pomocna funkcija za permutovanje 
      ekvivalentnih komponenata.Iz
svakog podskupa ekvivalentnih 
elemenata bira reprezentanta i gradi sve
permutacije*)

fEqChoices[Ulc_List]:=Module[{dp,ec, i},dp=Map[DistinctPermutations,Ulc];
    ec=KSubsets[Flatten[dp,1],Length[dp]];
    ec=Union[
        Table[If[SameQ[Union[Flatten[ec[[i]]]],Flatten[Ulc]],ec[[i]],{}],{i,
            Length[ec]}]];
    If[SameQ[First[ec],{}],Drop[ec,1],ec]]
fPletiLink[sara_List,poc_List]:=Module[{res},res=Map[poc[[#]]&,sara]]

(*## ## ## ## ## ## VARIJANTA A ## ## ## ## ## ## ## ## ## ##*)

(*Sredjuje jednu komponentu u zavisnosti 
od toga da li ima dvojnih elemenata
ili nema:ind {0,1} iz fFixC*)

fForFixC[LL_List,OstL_List,ind_Integer]:=Module[{res={},pom,pom1,poc, i},  
    If[ind==0,   
      (*Print["Indikator", ind];*)
            Do[If[OddQ[i],
                      If[SameQ[Position[OstL,LL[[i]]],{}],
                                            res=Append[res,{0,LL[[i]]}],     
                                             
            res=Append[res,{Position[OstL,LL[[i]]][[1,1]],LL[[i]]}]], 
                      
          res=Append[res,{Position[OstL,LL[[i]]][[1,1]],LL[[i]]}]]
                ,{i,1,Length[LL]}](*;
        Print["ind je 0 a ovo je rez ",res]*)
      ,
           Do[If[SameQ[Position[OstL,LL[[i]]],{}],
                         res=Append[res,{0,LL[[i]]}],
                          
          res=Append[res,{Position[OstL,LL[[i]]][[1,1]],LL[[i]]}]],
        {i,1,Length[LL]}]
      ];
    
    pom=Min[Map[Take[#,1]&,res]];
    poc=Select[res,SameQ[#[[1]],pom]&];
    poc=Union[Flatten[Map[Take[#,-1]&,poc],1]];
    (*Print["pocetni ",poc];*)
    pom=Map[Position[LL,#][[1,1]]&,poc];
    (*Print["pos, rotacije ", pom];*)
    If[ind==0,
      res=
        Map[{RotateLeft[res,#-1],Reverse[RotateRight[res,Length[res]-#]]}&,
          pom];
      res=Map[Flatten[#]&,Flatten[res,1]];
      (*Print["pre uzimanja :", Length[res],res];*)
      
      res=Map[Take[#,{2,Length[#],2}]&,res]
      (*Print["Uzeli smo sve ",res]*),
      res=fMin2Dvojne[LL]
      ];
    If[Length[res]==2,
          If[Count[pom[[1]],0]==2, 
                pom=Map[RotateLeft[#,Position[#,poc[[1]]][[2]]]&,res];
                pom=Map[Reverse[#]&,pom];
                res=Union[Flatten[Map[Append[res,#]&,pom],1]]
                 ],
        If[Count[pom[[1]],0]==2, 
              pom=RotateLeft[res[[1]],Position[res[[1]],poc[[1]]][[2]]];
              res=Append[res,Reverse[pom]]
           ]];
    res=Map[Prepend[OstL,#]&,res];
    
    res]

(*## ## ## ## # sredjuje jednu komponentu A ## ## ## ## ## ##*)

(*lista sa ExtDowkerima u 
fixiranim pozicijama tj.komponente su zauzele svoja
mesta*)
(*nn oznacava broj komponente sa kojom radimo*)

fMin2Dvojne[komp_List]:=Module[{poc,pom,pos,
posraz,minrast,i,res={},n,dd,dd1,
dd2,dd3},
    n=Length[komp];
    poc=Select[Union[komp],Count[komp,#]==2&];
   
    pos=Map[{#,Flatten[Position[komp,#]],
    Flatten[Position[Reverse[komp],#]]}&,
        poc];
    posraz=Map[{#[[1]],#[[2,1]],#[[2,2]]-#[[2,1]],
                    #[[3,1]],#[[3,2]]-#[[3,1]]}&,pos];
    posraz=Map[Insert[#,n-#[[3]],4]&,posraz];
    posraz=Map[Append[#,n-#[[6]]]&,posraz];
   (* Print["rastojanja:",posraz];*)
    
    minrast=Min[Flatten[Map[{#[[3]],#[[4]],#[[6]],#[[7]]}&,posraz]]];
    (*Print["Min ",minrast];*)
    Do[ dd=Position[komp,poc[[i]]][[2,1]];
    dd1=Position[Reverse[komp],poc[[i]]][[2,1]];
    dd2=Position[komp,poc[[i]]][[1,1]];
    dd3=Position[Reverse[komp],poc[[i]]][[1,1]];
     If[posraz[[i,3]]==minrast,
                    res=Append[res,RotateLeft[komp,dd2-1]]];
            If[posraz[[i,4]]==minrast,res=Append[res,RotateLeft[komp,dd-1]]];
            If[posraz[[i,6]]==minrast,
                  res=Append[res,RotateLeft[Reverse[komp],dd3-1]]];
             If[posraz[[i,7]]==minrast,
                   res=Append[res,RotateLeft[Reverse[komp],dd1-1]]]
      ,{i,1,Length[poc]}];
       (* Print[res];*)
    res
    ]

fFixC[tcL_List,nn_Integer]:=Module[{p,pre,ost,pom,res},
    p=tcL[[nn]];
    pom=Flatten[Map[Count[p,#]&,Union[p]]];
    pom=Count[pom,2];
    (*broj dvojnih *)
    (*Print["Komponenta ",nn," Broj Dvojnih ",pom];*)
  
      If[pom>=2,
      (*ako ima vise od dva dvojna onda je fiksirana *)
      
      pom=fMin2Dvojne[tcL[[nn]]];
      pre=Take[tcL,nn-1];
      ost=Drop[tcL,nn];
      If[pre=={},res=Map[{{#},ost}&,pom],
        If[ost=={},res=Map[{pre,{#}}&,pom],
          res=Map[{pre,{#},ost}&,pom]]];
      res=Map[Flatten[#,1]&,res],
      If[pom==1,
            res=fForFixC[p,Take[tcL,{nn+1,Length[tcL]}],1]];
      If[pom==0,
            res=fForFixC[p,Take[tcL,{nn+1,Length[tcL]}],0]];
      
         If[nn!=1,
               res=Map[Join[Take[tcL,{1,nn-1}],#]&,res]],
                res=fForFixC[p,Take[tcL,{nn+1,Length[tcL]}],0];
                If[nn!=1,res=Map[Join[Take[tcL,{1,nn-1}],#]&,res]]]
    ;
    res]
(*sredjuje sve komponente koje nemaju zajednickih elemenata sa svojim \
prethodnicima*)

fFixLAFirst[LinkL_List,indL_List]:=Module[{i,n,res,pre={}},n=Length[LinkL];
    Do[
         If[SameQ[i,1],(*1.komponenta*)
        res=fFixC[LinkL,1];
        pre=LinkL[[1]];
        ,(*ostale ako je 0 radimo inace ne*)
        If[SameQ[indL[[i]],0],
          res=Union[Flatten[Map[fFixC[#,i]&,res],1]];
          pre=Union[pre,Flatten[LinkL[[i]]]],
          pre=Union[pre,Flatten[LinkL[[i]]]]]],{i,1,n}];
    res]

(*## ## ## ## ## ## ## # VARIJANTA B ## ## ## ## ## ## ##*)

(*ind iz FixLink AB 1*)
(*fSeePrevNext posmatra susede i uzima bolji*)
\
(*pomocna za fFixLinkB*)
(*nn=1 ima neparnih*)

fSeePrevNext[CompL_List,PosLL_List,nn_Integer]:=Module[{pom,res},
    Switch[nn,0,pom={Take[CompL,{2}][[1]],Last[CompL]},1,
      pom={Take[CompL,{1}][[1]],Take[CompL,{3}][[1]]}];
    pom=Flatten[{Union[Select[PosLL,SameQ[#[[2]],pom[[1]]]&]],
          Union[Select[PosLL,SameQ[#[[2]],pom[[2]]]&]]},1];
    pom={pom[[1,1]],pom[[2,1]]};
    If[nn==0,pom=Reverse[pom]];
    (*If[pom[[2]]>pom[[1]],
               res=CompL,
          Switch[nn,1,res=Reverse[RotateLeft[CompL,3]],0,
            res=Reverse[RotateLeft[CompL,1]]];
          If[SameQ[pom[[1]],pom[[2]]],res=Prepend[{res},CompL]]];*)
    
    res={CompL};
    Switch[nn,1,res=Append[res,Reverse[RotateLeft[CompL,3]]],
      0,res=Append[res,Reverse[RotateLeft[CompL,1]]]];
    res]


(*fFixCB sredjuje B komponentu*)
(*pomocna za fFixLinkB*)

fBFixC[CompL_List,PList_List,nn_Integer]:=Module[{pom,minel,minpos},
    If[SameQ[nn,1],pom=Select[PList,OddQ[#[[1]]]&],pom=PList];
    minel=First[Sort[pom]][[2]];
    minpos=Position[CompL,minel][[1]];
    Switch[nn,1,p=RotateLeft[CompL,minpos-2],0,p=RotateLeft[CompL,minpos-1]];
    p=fSeePrevNext[p,PList,nn];
    If[Head[p[[1]]]==Integer,p={p}];
    p]

fGlupa[LL_List,bb_Integer,n_Integer]:=Module[{pre,pom,Li,PosL,PomLink},
    pre=Flatten[Take[LL,n-1]];
    PomLink=LL;
    (*ako je bb 1 onda je B i radimo je *)
    If[SameQ[bb,1],
      Li=PomLink[[n]];
      PosL=
        Map[If[SameQ[
                Position[pre,#],{}],{1000,#},{Position[pre,#][[1,1]],#}]&,
          Li];
      If[SameQ[Select[PosL,OddQ[#[[1]]]&],{}],
                         (*nema neparnih  0*)
                         
        pom=fBFixC[Li,PosL,0],
                         pom=fBFixC[Li,PosL,1]];
      PomLink=Map[Insert[Drop[PomLink,{n}],#,n]&,pom]];
    PomLink]

(*fFixLinkB sredjuje sve komponente koje su B u linku gde su A vec sredjene*)
\

fFixLBSecond[LLink_List,ind_List]:=Module[{i,n,res},
    res=LLink;n=Length[LLink];
    Do[If[SameQ[Head[res[[1,1]]],List],
        res=Map[fGlupa[#,ind[[i]],i]&,res];
        res=Flatten[res,1];
        If[SameQ[ind[[i]],0],res=Partition[res,n]]];
      If[SameQ[Head[res[[1,1]]],Integer],res=fGlupa[res,ind[[i]],i]],{i,2,
        n}];
    res]

fLinkAB[LinkL_List]:=
  Module[{ind={},i,pre,res},
    Do[If[i==1,ind={0};pre=Union[LinkL[[1]]],
        If[SameQ[Intersection[pre,LinkL[[i]]],{}],(*A prethodne ne uticu-0*)
            ind=Append[ind,0],(*B prethodne uticu-1*)ind=Append[ind,1]];
        pre=Union[pre,LinkL[[i]]]],{i,1,Length[LinkL]}];
    res=fFixLAFirst[LinkL,ind];
    res=Flatten[Map[fFixLBSecond[#,ind]&,res],1];
    res]

(*radi od Konveja,
  pdata ili dowkera sa znacima*)
(*!!!!Promena za cvorove u select*)

(*pomocna za MinDowProj*)

fPickEven[LL_List]:=Module[{res},res=Flatten[LL[[2]]];
    res=Union[Map[EvenQ[#]&,res]];
    SameQ[res,{True}]]
(****)

MinDowProjAltKL[Ul_]:=
  Module[{p,vrti,AllLinks,d,Con=Ul,i},
    If[SameQ[Ul,"2"]||SameQ[Ul,"-2"],p={{1,1},{4,2}},
      If[SameQ[Head[Ul],List]&&SameQ[Ul[[1]],{1,1}],p={{1,1},{4,2}},
        If[SameQ[Head[Con],List],
          If[SameQ[Union[Map[OddQ,Abs[Con[[2]]]]],{True}],(*input pdata*)Con=
              fDowfromPD[Con]]];
        (*sad imamo dowkera ili Konveja*)d=fGaussExtSigns[Con];
       If[SameQ[d,Flatten[d]],(*cvor*)AllLinks=fMinComp[d];
         AllLinks=Map[fKnittDow[d,Take[#,-2]]&,AllLinks];
         AllLinks=Map[fDowfromGaussExt[#[[1]]]&,AllLinks];
         p=Sort[Table[{Abs[AllLinks[[i]]],AllLinks[[i]]},
         {i,Length[AllLinks]}]][[1,2]],  
         (* KORIGOVANO 11.01.2004*)
          (*LINKOVI*)
          p=fOrderComp[Con];
          (*komponente sredjene po nekoliko kriterijuma,
            a drugi deo daje klase ekvivalencije*)
          vrti=Map[Flatten[#,1]&,fEqChoices[p[[2]]]];
          AllLinks=Map[fPletiLink[#,p[[1]]]&,vrti];
          AllLinks=Flatten[Map[fLinkAB[#]&,AllLinks],1];
          (*AllLinks=Flatten[AllLinks,1];*)
          AllLinks=Map[fDowfromGaussExt[#]&,AllLinks];
          AllLinks=Union[Select[AllLinks,fPickEven[#]&]];
          p=First[Sort[Abs[AllLinks]]];
          AllLinks=Union[Select[AllLinks,SameQ[Abs[#],p]&]];
          p=Last[Sort[AllLinks]]
          (*p=
              First[Sort[AllLinks,
                  Abs[#1[[2]]]<Abs[#2[[2]]]&]]*)]]];(*kraj prvog if-
        a*){p[[1]],Flatten[p[[2]]]}]


MinDowAltKL[Con_String]:=Module[{l,res},l=fProjections[Con];
    If[Length[l]==1&&Head[l[[1]]]==List,l=l[[1]]];
    res=Map[MinDowProjAltKL[#]&,l];
    
    l=First[Sort[Abs[res]]];
    res=First[Union[Select[res,SameQ[Abs[#],l]&]]];
    res={res[[1]],Flatten[res[[2]]]}]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
\
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
(*28.8.2003  po abs- za buduce narastaje- srediti*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## #*)
 
 (*## ## ## ## ## ## ## # AMPHICHEIRALITY ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)

(*AmphiProjAltKL[Ul_] testira amphicheiralnost projekcije date
Konvejevim simbolom dowkerom ili pdata  *)

      
  AmphiProjAltKL[Ulaz_]:=Module[{G,sG,sG1,Aut,wr,res, i,j,k,Ul},
  If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    G=fGaussExt[Ul];
    sG=fGraphInc[Ul][[2]];
    wr=Apply[Plus,sG];
    If[wr==0,sG1=-sG;
      G=FromUnorderedPairs[fGraphInc[Ul][[1]]];
      Aut=Drop[Automorphisms[G],1];
      (*Print[Aut];*)(*ovde izdvajamo automorfizme koji cuvaju znak*)Aut=
        If[Aut!={},
          If[First[
                Union[Table[
                    If[SameQ[
                          Sign[Table[
                                Table[
                                  Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                    Length[sG]}],{i,Length[Aut]}]][[k]],
                          sG1]==True,Aut[[k]],{}],{k,
                      Length[Aut]}]]]=={},
            Drop[Union[
                Table[If[
                    SameQ[Sign[
                            Table[Table[
                                Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                  Length[sG]}],{i,Length[Aut]}]][[k]],
                        sG1]==True,Aut[[k]],{}],{k,Length[Aut]}]],1],
            Union[Table[
                If[SameQ[
                      Sign[Table[
                            Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]],{j,
                                Length[sG]}],{i,Length[Aut]}]][[k]],
                      sG1]==True,Aut[[k]],{}],{k,Length[Aut]}]]],Aut];
      If[Aut!={},res=1,res=0];
      Aut=If[Aut!={},Map[ToCycles,Aut],Aut];
      Aut=Union[Map[Length,Flatten[Aut,1]]];
      If[Aut!={}&&First[Aut]==1,Drop[Aut,1],Aut];
      , res=0 ]; 
    res]
    
  

fAmphiAltKL[Con_String]:=Module[{gr,dd,res},
    gr=fGraphKL[Con];
    dd=fDualGraphKL[Con];
    res=If[IsomorphicQ[FromUnorderedPairs[gr],FromUnorderedPairs[dd]],1,0];
    res]
    
(* izmenjena funkcija S.J. 13.03.2009 *)


AmphiAltKL[Con_String]:=Module[{pp,pp1,res},
pp=Map[First,fDiffProjectionsAltKL[Con]];
pp1=Flatten[Position[Map[fAmphiAltKL,pp],1]];
res=If[SameQ[pp1,{}],1,pp[[pp1[[1]]]]];
res]
    
    
AmphiQ[Ul_List]:=Module[{pd},
    pd=GetMirrorImageKnot[Ul];
    If[SameQ[RedKauffman[Ul],RedKauffman[pd]],1,0]
    ]
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)     
     
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*## ## ## ## ## ## ## # CUTTING NUMBER ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)

(* ^^^^^^^^^^^^^^^
(*  koristi se za UInfty *)
fAddSign[LL_List,{el_Integer,elS_Integer}]:=Module[{z,i,p,l=LL},
    z=1-elS;
    p=Position[LL,{el,0}][[1,1]];
    Do[
       l=ReplacePart[l,{l[[i,1]],z},i] ;
      z=1-z
      ,{i,p,Length[LL]}];
    z=1-elS;
    Do[z=1-z;(*Print[l];*)
           l=ReplacePart[l,{l[[p-i,1]],z},p-i]
         
          ,{i,1,p-1}];
    l ]*)
    (*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*)
(*fGiveSign[LL_List,SL_List]:=Module[{p,i,l=LL,saznakom=SL,sz={},sa},
     If[SL=={},
            p=l[[1]]; 
            Do[   p=ReplacePart[p,{p[[i,1]],Mod[i,2]},i];
                  ,{i,1,Length[p]}];
            l=ReplacePart[l,Prepend[p,1],1];
             sz=p,
            (*ako vec imamo neke sa znacimo trazimo ih u drugim \
komponentama*)
\
                sa=Union[Flatten[Map[Take[#,1]&,saznakom]]];
         (*   Print["SSSSS",sa,saznakom];*)
            
      Do[  If[Length[l[[i,1]]]!=0,
                          (*znaci da ona nije oznacena*)
                     \
      presek=Select[l[[i]],MemberQ[sa,#[[1]]]&];
                          If[presek!={},
                                  
            presek=Select[saznakom,#[[1]]==presek[[1,1]]&][[1]];
                                     p=fAddSign[l[[i]],presek];
                                   l=ReplacePart[l,Prepend[p,1],i];
                                   p=Map[If[MemberQ[sa,#[[1]]],-1,#]&,p];
                               sz=Union[sz,Complement[p,{-1}]]
                             ]]
        ,{i,1,Length[l]}]
      ];
    {l,sz}
    ]*)
     (*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*)
(* Pravi alternirajuci kad se polazi od pData... *)

(*fMakeAlt[CoLe_List,LL_List]:=Module[{i,pom=0,
CoPom={},saznakom={},ind=True},

      l=iteratedTake[LL,CoLe];
      While[ind,
                 (*Print["SZ",saznakom];*)
                  
      saznakom=fGiveSign[l,saznakom];
                  (* 
        Print["posle fGiveSigna",saznakom[[1]],"\n \n",saznakom[[2]],
            "\n \n"];*)
                  l=saznakom[[1]];
      saznakom=Last[saznakom];
                 (* Print[Map[ Length[#[[1]]]&,l]];*)
                  
      ind=Union[Map[ Length[#[[1]]]&,l]]!={0}
            ];
     Map[Drop[#,1]&,l]
    ]*)
(*pomocna za fMakeGaussPD- menja odnose iznad ispod za neparne PData *)
(*^^^^^^^^^^^^^^^^^^^^^^^*)
fMenjaObePojavePD[LL_List,el_Integer]:=Module[{p,i,l=LL},
    p=Select[l,#[[1]]==el&];
    (*Print[p,p[[1,1]]];*)
    (*Print[Position[l,p[[1,1]]]];*)
    
    p=Map[Take[#,1][[1]]&,Position[l,p[[1,1]]]];
    Do[
      l=ReplacePart[l,{l[[p[[i]],1]],1-l[[p[[i]],2]]},p[[i]]],{i,1,2}];
    l
    ]
(*radi od Ext Gausa sa odnosima iznad ispod*)
(* 
  ako moze sredi parnost *)
(*NE RADI RETROAKTIVNO!!!! *)
(*^^^^^^^^^^^^^^*)
(*koristi je UInfty*)
(*fParno[LL_List]:=Module[{l,lGExt,i,p},
    l=LL;
    lGExt=Map[Flatten,l];
    lGExt=Map[Take[#,{1,Length[#],2}]&,lGExt];
    (*extGaus bez odnosa iznad ispod *)
    
    Do[p= Flatten[Position[Flatten[Take[lGExt,i-1]],lGExt[[i,1]] ]];
          (*pozicija prvog u sledecoj komponenti u 
          prethodnom delu gausa*)
   If[p!={},
               If[Not[EvenQ[p[[1]]]],
                     lGExt=ReplacePart[lGExt,RotateLeft[lGExt[[i]]],i];
                       l=ReplacePart[l,RotateLeft[l[[i]]],i]
                ]]  ,{i,2,Length[l]}];
    l]
(* ako se posle secenja za PD u istoj komponenti jave dva ista 
uzastopna broja (ciklicno) sa suprotnim drugim znakom 
npr. {-18,1},{-18,0} izbacujemo ih.*)
(*radi od gaus ext sa iznad ispod*)
(*koristi je UInfty*)
fIzbaciPetlju[GK_List]:=Module[{p,l,pom,res=GK, i},
    l=Map[Take[#,1]&,GK]; 
    el=Union[Select[l,Count[l,#]==2&]];
    p=el;
    If[p!={},
          p=Map[Union[Flatten[Position[l,#]]]&,p];  
            l=Map[If[#!={1,Length[GK]},RotateLeft[l,#[[2]]-1],l]&,p]; 
       p={};
         Do[ p=Append[p,Flatten[Position[l[[i]],el[[i]]]]]
          ,{i,1,Length[el]}];
                l=GK; 
         Do[If[p[[i]]=={1,Length[GK]},
                  l=Select[l,#[[1]]!=el[[i,1]]&]
               ]
        ,{i,1,Length[el]}];res=l
      ];
    res]*)
    (*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*)
    (*koristi je UInfty*)
fMakeGaussPD[PD_List]:=Module[{l,n,i,l1={},nepar},
      n=Length[PD[[2]]];
     l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData *)
    
    Do[l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l]; *)(*uredjeni parovi sa sve znacima*)
      
    nepar=Select[l,Or[OddQ[#[[2]]],And[EvenQ[#[[1]]],EvenQ[#[[2]]]]]&];
       nepar=Flatten[Map[Take[#,-1]&,nepar]];
    (*Print[nepar];*)
    l1=Map[#-#&,Range[2n]];
    Do[ l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,1]]]];
            l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,2]]]],{i,1,n}];
    (*Print["???",2PD[[1]],
          l1];*)
    (*ubacili smo indikatore i sad treba da ga napravimo da \
je alternirajuci*)
    l1=fMakeAlt[2PD[[1]],l1];
    (*menjamo odnose onima koji su u pdata bili neparni a to su i ovde*)
    
    l1=Flatten[l1,1];
    (* Print[" pre kvarenja ",l1]; *)
    
    If[nepar!={},
      Do[l1=fMenjaObePojavePD[l1,nepar[[i]]]    ,{i,1,Length[nepar]}]];
    (* Print[" posle kvarenja ",l1]; *)
    l1
    ]
(*pomocna za fBrCo- izbacuje iz Gaus parova i
 zavisnosti od "ociscenog"Gausa*)
(*^^^^^^^^^^^^^^^^^^^^^^^^^^*)
fBezComp[GPairs_List, GG_List]:=Module[{ll,p},
    ll=GPairs;
    Do[ p=Select[ll[[i]],MemberQ[GG[[i]],#[[1]]]&];
               ll=ReplacePart[ll,p,i];
      ,{i,1,Length[GG]}];
    ll
    ]
fUredjenobezdvojnih[LL_List]:=Module[{i,redLL=Flatten[LL]},
    Do[  redLL=
        Join[Take[redLL,i],Select[Drop[redLL,i],#!=redLL[[i]]&]],{i,
        1,Length[redLL]/2}];
    redLL
    ]
f2Ista[l_List]:=Module[{i=1,res=0},
    If[Length[l]==2,
            If[l[[1,2]]==l[[2,2]],
                      res=2],
                   While[i<Length[l]-1&&l[[i,2]]!=l[[i+1,2]],i++];
                            If[i!=Length[l]-1,
                            res=i+1,
                              If[l[[i,2]]==l[[i+1,2]],res=i+1,res=0]]];
    res]

(* pomocna za fDataForPData - alternira komponentu *)

fMakeCompAlt[Komp_List,ind_]:=Module[{l=Komp},
    l=Map[{#[[1]],Mod[Position[l,#][[1,1]],2]}&,l];
     l={l,Map[{#[[1]],1-#[[2]]}&,l]};
    If[ind!=-1,l=Select[l,#[[1,2]]==ind&]];
    l
    ]
(*pomocna za fPDataCuttNo, radi pd ext Gausa koji je napravljen u fBreakCp*)
(* 
  rezultata je lista svih mogucih alternirajucih*)

fDataForPData[GE_List]:=Module[{res,p,i},
    res={fMakeCompAlt[GE[[1]],GE[[1,1,2]]]};
    Do[ p=fMakeCompAlt[GE[[i]],-1];
           res=Join[Map[Append[#,p[[1]]]&,res],Map[Append[#,p[[2]]]&,res]]
      ,{i,2,Length[GE]}];
    res
    ]
(*iz liste svih potencijalnih bira dobru *)
(*kako?*)
(* 
  Oduzme potencijalnu od pocetne*)
(* 
  dobije listu uredjenih parova koji oblika {broj,a}
      gde je a 0- ako su znaci isti
               1- ako su razliciti*)
(* ako ima 1 trazimo drugu
      ako nema 1- dobra je;0 su dobri a 2 pokvareni*)

fSelectForPData[LL_List,poc_List]:=Module[{i=1,j,ind=True,pom,pF},
              pF=Flatten[poc];
            While[ind,
                         pom=Flatten[LL[[i]]];j=2;
                 
      While[j<=Length[pom],
        pom=ReplacePart[pom,Abs[pom[[j]]-pF[[j]]],j];
                          j=j+2
            ];
      pom=Partition[pom,2];
      If[Length[pom]==2Length[Union[pom]],
           (*nasli smo pravi!!!*)
                ind=False;
            pom=Select[Union[pom],#[[2]]==1&];
             pom=Flatten[Map[Take[#,1]&,pom]]
        ];
      i++];
    pom
    ]
(*radi od gaus EXT sa "nekim" iznad ispod*)
(* 
  daje pData *)
(*pomocna za fBReak CO*)

fPDataCuttNo[GE_List]:=Module[{l=Flatten[GE],svi,pokvareni},
    (*pravimo se mogucnosti za alternirajuci*)
    svi=fDataForPData[GE];
    (*nalazimo prvu "dobru" alternirajucu i ona nam daje njene pokvarene*)
   
     pokvareni=fSelectForPData[svi,GE];
    l=Take[l,{1,Length[l],2}];
    (*za pdata bez iznad ispod *)
    {Map[Length[#]&,GE]/2,
      fPDizNiz[l,pokvareni]}]
(*Formira PData iz nase vrste gausa*)

fPDizNiz[niz_List,gr_List]:=Module[{n=Flatten[niz]},
    n=Union[Map[{#,Flatten[Position[n,#]]}&,n]];
    n=Map[If[EvenQ[#[[2,1]]],{#[[1]],{#[[2,2]],#[[2,1]]}},#]&,n];
    n=Map[If[MemberQ[gr,#[[1]]],
                          {#[[2,2]],Sign[#[[1]]]*#[[2,1]]},
                          { #[[2,1]],Sign[#[[1]]]*#[[2,2]]}   
                          ]&,n];
    n=Sort[n,#1[[1]]<#2[[1]]&];
    n=Take[Flatten[n],{2,2Length[n],2}];
    n]
fBrCo[Ul_,k_Integer]:=Module[{l,p,pd1,p1,pd,pd0,nepar,res,pom},
               pd=Ul;(* sad nam je pd PData *)
                
    pd=ReductionKnotLink[pd];
                (*sad moramo da napravimo odgovarajuceg Gausa i da ga iz \
alterniramo*)
                
    pd0=fMakeGaussPD[pd];(* 
      Gaus parovi sa iznad-ispod tj. 0-1 *)
                
    pd1=Flatten[
        Map[Take[#,1]&,pd0]]; (*Gaus bez iznad ispod*)
                
    pd0=iteratedTake[pd0,2pd[[1]]];        
    	(*gaus za seckanje podeljen na komponente*)
                
    pd1=iteratedTake[pd1,2pd[[1]]];
                
    p=pd0[[k]];(*k-ta tj. komponenta koju secemo*)
                 
    p1=Flatten[Map[Take[#,1]&,p]];(*p bez iznad ispod*)
                 
    pd1=Map[UnsortedComplement[#,p1]&,Drop[pd1,{k}]];
                
    pd0=Drop[pd0,{k}]; (*gaus parovi bez izbacene komponente*)
               \
 Do[ pd0=fBezComp[pd0,pd1] ,{i,1,Length[pd0]}];
                (*do petlja izbacuje one koji nisu preostali *)      
    	 pd0=UnsortedComplement[pd0,{{}}];
    	 (*sad nam treba njihov pravi poredak tj. 
          izbacujemo drugu pojavu*)
              
    pd1=fUredjenobezdvojnih[pd1];
              l=Map[Length,pd0];
              If[pd1=={},
                  (* Print["Direct product or 2-component link"]; *)
        res={{},{}},
      (*  Print["posle secenja za PD ",pd0];*)
      
      pd0=Map[fIzbaciPetlju[#]&,pd0];(*sredjujemo parnost- 
          ako ne moze da sredi dobicemo dva neparna i dva parna uparena*)
    
        pd0=fParno[pd0];
      (*Print["Posle parnosti ",pd0];*)
     If[SameQ[pd0,{{}}],
        res={{},{}},
      pd0=fPDataCuttNo[pd0];
      (*Print["Pre redukcije ",pd0];*)
      res=ReductionKnotLink[pd0]]];
    res
    ]
          
(*rezultat salje u fBreak Comp *)
fBreakComp[Ul_,k_Integer]:=Module[{pd={},i},
If[k>fComponentNo[Ul],
Print["Maximum number of components is: ",fComponentNo[Ul]],
    If[SameQ[fComponentNo[Ul],1],
              Print[ "Knot"],
                If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
                (*ako je konvej pravimo PD ako ne sve je PD*)
                \
  pd=fBrCo[pd,k]]];
    pd]
(*radi od pdata i daje rezultat pri jednom secenu svake komponente ponaosob*)

BreakCoAll[Ul_]:=Module[{i=1,res={},p={3}},
    If[fComponentNo[Ul]==1,res={};Print["Knot"],
      While[p!={1}&&i<=fComponentNo[Ul],
                    (* Print["Ulaz: ",Ul];*)
                  
        p=fBreakComp[Ul,i];
                (*  
          Print["Isecena je ",i,"-ta komponenta: ",p];*)
                 
        If[Not[MemberQ[{{{},{}},{{0},{}},{{},{0}}},p]],
                       i++;res=Append[res,p],
                      p={1};res={0}]
             ];
        If[p!={1},
        res=Table[{res[[i,1]],Flatten[res[[i,2]]]},{i,Length[res]}]]];
    res
    ]
(*vraca {0} ako je razvezao a u suprotnom vraca listu ili {} ako je cvor *)

CuttNo[Ulaz_]:=Module[{pd,cNo,i,R,R1={},p={3},res="0",Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
       If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    (*sad imamo pdata -odmah reduction*)
    
    cNo=MemberQ[ReductionKnotLink[pd],{}];
    If[fComponentNo[pd]==1,
         Print["Knot"];cNo=0,
          pd=ReductionKnotLink[pd];
        If[cNo==True,Print["Razvezan posle redukcije"];cNo=0,
        (*Ako je odmah redukovan- 
            cutting No je nula u suprotnom radimo dalje*)
        (*Print[cNo,
              res,pd];*)
        R=BreakCoAll[pd];
        (*Print["trenutna lista ",Length[R],R];*)
        
        cNo=1;(*postavljamo ga na 1 jer smo tacno jednu komponentu isekli*) 
        While[
          Not[MemberQ[R,0]],(*Print["RADI WHILE"];*)
                        
          R1={};i=0;
                         While[p!={0}&&i<Length[R],  i++;
                                     p=BreakCoAll[R[[i]]];
                                     
            R1=Append[R1,p](*Print[i,"tren ",
                R1]*)
                                   ];
                       If[p=={0},R={0},R=Flatten[R1,1]];
            (*  Print["TRENUTNA ",R];*)
             cNo++]
        ]];(*kraj prvog if-a*)
    cNo
    ]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ##*)
(*## ## ## ## ## ## ## # REAL CUTTING  ## ## ## ## ## ## ## ## ## ## ## ## *)

(*radi od gausEXT- sva moguca "secenja sa fix krajevima"*)

fCuttComp[gg_List]:=Module[{res={},p=gg, i},
    Do[  p=RotateLeft[p];
            res=Join[res,{p},{Reverse[p]}]
      ,{i,1,Length[gg]}];
    Union[res]
    ]
(*vraca {a,b}
      a=1 ako jeste izomorfizam 0 inace
    i b=1 ako je sense preserving
        b=-1, ako je sense reversing*)

fCuteIso[L1_List,L2_List]:=Module[{l={},i},
    Do[l=Append[l,{L1[[i]],L2[[i]]}],{i,1,Length[L1]}];
    l=Sort[Union[l],Abs[#1[[1]]]<Abs[#2[[1]]]&];
    If[2Length[l]==Length[L1],i=1,i={0,0}];
    If[i==1,
             If[Union[Map[#[[1]]*#[[2]]>0&,l]][[1]],
                (*sense-preserving *)    i={i,1},
                  (*sense=reversing*)i={i,-1}]];
    i]
(*radi sa EXTGauss*)
(*iz njihove liste vraca predstavnike klasa izomorfizama \
i to ako je ind=1 onih koje cuvaju znake
 a -1 ako ne vodi racuna o tome*)

fCuteIsoClass[LL_List,ind_]:=
  Module[{l=LL,res,i,pom,SenseRev={},SensePres={},p,f},
    l=Drop[l,1];
    res={LL[[1]]};
    While[l!={},
                    p=Last[res];i=1;
                    While[i<=Length[l],
                      pom=fCuteIso[p,l[[i]]];
                    
                       If[ind==1,f={1,1},pom=pom[[1]];f=1];
                     (*ako je ind 1 onda moraju cuvati
                ako je -1 onda je smao prvi clan bitan da je 1*)
        	   
        If[pom==f,
                          l=Drop[l,{i}],i++]
                     ];
      (*  Print["Ovde mora biti prazna za 3 ",l];*)
           
      If[l!={}, res=Append[res,First[l]];
               l=Rest[l]];
            If[Length[l]==1,res=Append[res,l[[1]]];l={}]
      ];
    res
    ]
(*pravi linkove tako sto na el dodaje sve iz r*)
(* i sve to flatten*)

fJoinL[el_List,r_List]:=Module[{res={el}},
    res=Map[Flatten[Append[res,#],1]&,r];
    (*Print["duzina rezultata ",Length[res]];*)res
    ]
(*Radi od PData ili Konveja*)

fCuttRealKL[Ulaz_]:=Module[{res,pd,gaus,sPreserv,svi,brCo,len, UL},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],UL=Ulaz,
    UL=ToExpression[Ulaz]];
    pd=UL;
    If[SameQ[Head[UL],String],pd=fCreatePData[UL]];
    pd=ReductionKnotLink[pd];
    gaus=fGaussExtSigns[fDowfromPD[pd]];
    If[Head[gaus[[1]]]==Integer,gaus={gaus}];
    gaus=Map[fCuttComp[#]&,gaus];
    (*Print["gaus",
          Map[Length[#]&,gaus]];*)
          (*sad imamo sve moguce rasecene gause
        podliste su mogucnosti za svaku komponentu*)
    (*sredjujemo linkove- i oni postaju 1 komponenta*)
        brCo=Length[gaus];
    If[brCo!=1,
            len=Map[Length[First[#]]&,gaus];
   (*lista duzina svake komponente spremno za iterated take*)
    res=gaus[[1]];(*Print["res ",res];*)
    Do[ res=Flatten[Map[fJoinL[#,gaus[[i]]]&,res],1]
        ,{i,2,brCo}];
      gaus={res}
         ];
    (*Print["Gaus za POCETNIKE ",
          Length[gaus]];*)
    (*trazimo neizomorfne predstavnike *)
   (*cuvaju znake*)
    sPreserv=Flatten[Map[fCuteIsoClass[#,1]&,gaus],1];
    svi=Flatten[Map[fCuteIsoClass[#,-1]&,gaus],1];
     {Length[sPreserv], Length[svi]}]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## #*)
(*## ## ## ## ## ## ## # SPLITTING NUMBER  ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## *)

(*radi od gauss Ext broj disjunktnih tj, razdvojenih komponenti *)

fDisjointComp[LL_List]:=Module[{i,br=0},
    If[Head[LL[[1]]]==Integer,
      br=1];(*ako je cvor onda je to smao 1 komponenta*)
    
    If[Head[LL[[1]]]==List,
      Do[
        If[Intersection[LL[[i]],Flatten[Drop[LL,{i}]]]=={},br++]
        ,{i,1,Length[LL]}]];
    br
    ]
(*Menja znak u pojedinacnoj tacki- radi od PData*)

fCrossChangePD[PD_List,pos_Integer]:=Module[{l,n,i},
      n=Length[PD[[2]]];
     l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData *)
    
    Do[l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l];*)
    l=ReplacePart[l,-1*Reverse[l[[pos]]],pos];
    (*Print["www",l];*)
    l=Sort[l,Abs[#1[[1]]]<Abs[#2[[1]]]&];
    l=Flatten[Map[Take[#,-1]&,l]];
    l=ReductionKnotLink[{PD[[1]],l}];
    l={UnsortedComplement[l[[1]],{{}}],l[[2]]};
    l={UnsortedComplement[l[[1]],{0}],l[[2]]}
    ]
(*Radi od PD!-promeni znake u svim tackama pojedinacno*)

fCrChAllPD[PD_List]:=Module[{res,l},
    l=Range[Length[PD[[2]]]];
    res=Map[fCrossChangePD[PD,#]&,l]
    ]
fGaussPDForDisjoint[PD_List]:=Module[{l,i,res,n},
      n=Length[PD[[2]]];
     l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData *)
    
    Do[l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    l=Abs[l];
    res=Table[0,{i,2n}];
    Do[res=ReplacePart[res,i,l[[i,1]]];
          res=ReplacePart[res,i,l[[i,2]]]
      ,{i,1,n}];
    iteratedTake[res,2PD[[1]]]
    ]
(* Radi od Konveja ili pdata*)

SplittNo[Ulaz_]:=Module[{br,pd,gaus,korak=0,CoNo,ind=True,PData},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],PData=Ulaz,
    PData=ToExpression[Ulaz]];
    pd=PData;
    If[SameQ[Head[pd],String],
              pd=fCreatePData[PData]];
    pd= ReductionKnotLink[pd];(*za svaki slucaj*)
    
    pd={UnsortedComplement[pd[[1]],{{}}],pd[[2]]};
    pd={UnsortedComplement[pd[[1]],{0}],pd[[2]]};
    CoNo=Length[pd[[1]]];(*broj komponenti i kriterijum za zaustavljanje*)
   
     If[CoNo==1,
             (*Knot *)korak=0,
              (*Link*)
                gaus=fGaussPDForDisjoint[pd];
                br=fDisjointComp[gaus];
                pd={pd};
                While[br<CoNo&&ind,
                          pd=Union[Flatten[Map[fCrChAllPD[#]&,pd],1]];
                           korak++;                        
        		
        ind=Not[MemberQ[
              Union[Map[MemberQ[pd,#]&,{{{},{}},{{0},{}},{{},{0}}}]],True]];
                            
        ind=ind&&Not[MemberQ[Map[Length[#[[1]]]==1&,pd],True]];
                            (*ind je True ako nema razvezanih  tj. 
              ako treba dalje splitati*)                
        		If[ind,
                              gaus=Map[fGaussPDForDisjoint[#]&,pd];
                               
          br=Max[Map[fDisjointComp[#]&,
                gaus]] ](*kraj if *)
          ](*kraj while*)
      ];(*kraj \
prvog if *)
    korak
    
    ]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ##*)
(****************************************************************)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## #*)

(*## ## ## ## ## ## ## # CONVERSION FUNCTIONS ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ##*)
fDToD[Con_String]:=Module[{n,p,DD},
    DD=Dowker[Con][[4]];
    If[Head[DD[[1]]]==List,
      n=Map[Length,DD];
      p={n,Abs[Flatten[DD]]*fGenSign[Con]}];
    If[Head[DD[[1]]]==Integer,n=Length[DD];
      p={{n},Abs[DD]*fGenSign[Con]}];
    p]

fKnotscapeDow[Con_String]:=
  Module[{l},(*ukoliko nema-u Konveju Dowker je nas abs*)
    l=Abs[fDToD[Con]];
    If[StringPosition[Con,"-"]!={},(*ako ima-*)
         l=Abs[fDToD[StringReplace[Con,"-"->""]]];
      l={l[[1]],
          l[[2]]*fGenSign[Con]*fGenSign[StringReplace[Con,"-"->""]]}];
    l]
    
    
    fKnotscapeDowfromDow[Ul_]:=Module[{pp,pp1,res},
    pp=fSignsKL[Abs[Ul]];
    pp1=Map[Sign,Ul[[2]]]*Map[Sign,pp[[2]]];
    res={Ul[[1]],-Abs[Ul[[2]]]*pp1};
    res]
    
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ## ## ## ## ## ## # MID-EDGE GRAPH  ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)

fMakeCycle[ll_List]:=Module[{ll1},
    ll1=Join[
        Table[Sort[{ll[[i]],ll[[i+1]]}],{i,
            Length[ll]-1}],{Sort[{ll[[1]],ll[[-1]]}]}];
    ll1]


(* radimo sa poliedrima *)
 fMidEdgeGraph[gr_List]:=Module[{mm,mm0,mm1,mm2,res,i,j},
    mm=Select[Split[Sort[gr]],Length[#]>1 &];
    mm=Flatten[Table[Table[mm[[j,1]],{i,Length[mm[[j]]]-1}],{j,Length[mm]}],
        1];
    mm0=mm;
    mm=Join[Last[fPlanarEmbGraph[gr]],mm];
    mm1=Map[fMakeCycle,mm];
    mm2=Union[Flatten[Map[fMakeCycle,mm],1]];
    res=Sort[
        Flatten[Table[
            fMakeCycle[
              Flatten[Table[
                  Position[mm2,mm1[[j,i]]],{i,Length[mm1[[j]]]}]]],{j,
              Length[mm1]}],1]];
    res=Select[res,Not[SameQ[#[[1]],#[[2]]]] &];
    res] 
    
    (* zamenio S.J. 13.03.2009 *)
    

    
(*## ## ## ## ## ## ## # fKLfromGraph  ## ## ## ## ## ## ## ## ## ## ## ## *)

(*pomocna za fFixInc-ubacuje simbol temena na narednu poziciju*)

fPosDouble[L_List,red_Integer,el_Integer]:=
  Module[{l=L,p},(*Print[l];*)p=
      Position[l[[red]],el][[1,1]];(*Print[red,el,"pos",p];*)
      If[p!=4,
      l=Insert[l,el,{red,p+1}],l=Insert[l,el,{red,4}]];
      (*Print[l];*)l]
(*Ubacuje u listu incidencije dvojne veze*)

fFixInc[IL_List,dvoL_List]:=
  Module[{res=IL, i},Do[res=fPosDouble[res,dvoL[[i,1]],dvoL[[i,2]]];
      res=fPosDouble[res,dvoL[[i,2]],dvoL[[i,1]]],{i,1,Length[dvoL]}];res]
(*koliko se puta u datom redu pojavljuje element sa pozicije el i u kom \
redosledu*)

fBrojPojave[L_List,red_Integer,el_Integer]:=Module[{p,res,r},p=L[[red,el]];
    r=Map[If[Length[#]==2,#[[1]],#]&,L[[red]]];
    If[Length[p]==2,p=p[[1]]];
    (*Print[p];*)
    If[Count[r,p]==1,res=0,
      If[el==5,res=2,
        If[el==2,res=1,Switch[p,r[[el-1]],res=2,r[[el+1]],res=1];]]];
    res]

fPrvi[inc_List]:=
  Module[{l=inc},
    l=Map[Rest[#]&,l];(*izbrisemo oznake reda*)l=
      Select[Flatten[l,1],Length[#]==0&];
    If[l=={},l=0,l=l[[1]]];
    l]
(*vraca 0 ako su svi iskoristeni a prvi sledeci ako takav postoji*)
UnsortedComplement[x_List,y__List]:=
  Replace[x,Dispatch[(#:>Sequence[])&/@Union[y]],1]
(*u redu RED trazi p-tu pojavu elementa el,pa element 2 pozicije udaljen*)

fTrazi[L_List, red_Integer, el_Integer, brojPojave_Integer] := 
  Module[{p, pom, inc = L},
    Switch[brojPojave, 0, p = 1, 1, p = 2, 2, p = 1];
    
    (*ako se vec pojavio (kao 1) onda uzimamo drugu pojavu, inace prvu*)
    
    pom = Position[Rest[L[[red]]], el];
    pom = Select[pom, Or[Length[#] != 2, Length[#] == 2 && #[[2]] != 2] &];
    
    pom = pom[[p, 1]];
    (*ne uzimamo u obzir prvi el. 
          jer je on oznaka reda i zato posle ide ++*)
    
    pom = pom + 1;
    
    inc = ReplacePart[L, {L[[red, pom]], 0}, {red, pom}];
    (* stavljamo mu 0 jer smo njega trazili*)
    If[pom <= 3, pom = pom + 2, pom = pom - 2];
    (* ovo je pozicija u red - u el. koji je sledeci clan komponente*)
    (*i sad jos da oznacimo da je dati upotrebljen*)
    
    If[Length[L[[red, pom]]] != 2,
      inc = ReplacePart[inc, {inc[[red, pom]], 1}, {red, pom}]];
    
    {L[[red, pom]], {red, pom - 1}, inc}]
(*vraca element, njegovu poziciju, nov inc*)


fSrediParnost[L_List]:=Module[{res={},p,i,j},res={L[[1]]};
    Do[j=1;p=L[[i]];
      While[Not[MemberQ[Flatten[res],p[[j]]]],j++];
      If[EvenQ[Position[res,p[[j]]][[1,2]]-j],
        res=Append[res,RotateLeft[L[[i]]]],res=Append[res,L[[i]]]],{i,2,
        Length[L]}];
    res]

fKLfromGraph[InG_List]:=
  Module[{p,GNew,incPoc,inc,dvojne,res={},komp={},
  el,red=1,brP,prvi,DUPLE={}},
    p=Map[Sort[#]&,InG];
    If[Map[Count[p,#]&,Union[p]]=={4},res={{1,1},{4,2}},
      If[Union[Map[Count[p,#]&,p]]=={2},(*ovo je cvor ili link "n"*)el=
          Length[p]/2;
        res=Dowker[ToString[el]][[4]];
        If[OddQ[el],res={{el},Flatten[res]},res={{el/2,el/2},Flatten[res]}],
        p=fPlanarEmbGraph[p];
        GNew=p[[1]];(*planar embeded graf-Unorderedpairs*)
        inc=p[[2]];(*lista incidencije*)dvojne=
          Union[Select[GNew,Count[GNew,#]==2&]];
         inc=fFixInc[inc,dvojne];
        el=inc[[1,2]];brP=fBrojPojave[inc,1,2];
        komp=Append[komp,el];
        inc=ReplacePart[inc,{inc[[1,2]],1},{1,2}];
        (*ovaj el smo vec iskoristili*)
        prvi=el;(*pocetni element u komponenti*)
      
        While[
          And[el!=0,
            Union[Map[Count[Flatten[res],#]&,Flatten[res]]]!={2}],
           
          While[And[Length[el]!=2,Count[komp,prvi]!=3],
            p=el;
            el=fTrazi[inc,el,red,brP];
           inc=el[[3]];
            If[Length[el[[1]]]!=2,red=p;
              brP=fBrojPojave[inc,el[[2,1]],el[[2,2]]+1];
              el=el[[1]];
              (*nov element i njegov broj pojave*)
              (*nalazi u redu el element 
jednak redu+2 ili-2 pozicije*)komp=
                Append[komp,
                  el],(*drugi-4 put na isti "nacin" prolazimo kroz jednu \
tacku-znaci kraj*)
          el=el[[1]]]];

         If[First[komp]==Last[komp],komp=Drop[komp,-1]];
          res=Append[res,komp];
          res=UnsortedComplement[res,{{}}];
          el=fPrvi[inc];
          If[el!=0,komp={el};
            (*posto pocinjemo novu komponentu vracamo inc na nekoristeno \
stanje*)red=Position[inc,el];
            red=Select[red,Length[#]==2&&#[[2]]!=1&][[1]];
            inc=ReplacePart[inc,{inc[[red[[1]],red[[2]]]],1},red];
            (*stavljamo daje iskoristen*)
           
            brP=fBrojPojave[inc,red[[1]],red[[2]]];
            red=red[[1]];]];(*kraj while*)
           res=fSrediParnost[res];
       res=fDowfromGaussExt[res];
        res={res[[1]],Flatten[res[[2]]]};
     If[fComponentNo[res]==1,res=fMinDowKnotPr[res]]]];(*kraj \
ifova*)
      res]
      
      
      
      
      (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## DOWKER CODES \
#### ## ## ## ## ## ## ## ## ## ## #*)


fDowkerCode[Ul_List] := 
  Module[{inc = Ul, res = {}, komp = {}, el, red = 1, brP, prvi, 
  DUPLE = {}, p},
    el = inc[[1, 2]]; brP = fBrojPojave[inc, 1, 2];
    komp = Append[komp, el];
    inc = ReplacePart[inc, {inc[[1, 2]], 1}, {1, 2}];
    (*ovaj el smo vec iskoristili*)
    prvi = el;(*pocetni element u komponenti*)
    While[
      And[el != 0, Union[Map[Count[Flatten[res], #] &, Flatten[res]]] != \
{2}],
       While[And[Length[el] != 2, Count[komp, prvi] != 3], p = el;
        el = fTrazi[inc, el, red, brP];
        inc = el[[3]];
        If[Length[el[[1]]] != 2, red = p;
          brP = fBrojPojave[inc, el[[2, 1]], el[[2, 2]] + 1];
          el = el[[1]];
          (*nov element i njegov broj pojave*)(*nalazi u redu el element \
jednak redu + 2 ili - 2 pozicije*)komp = 
            Append[komp, el],(*drugi - 
              4 put na isti "nacin" prolazimo kroz jednu tacku - znaci kraj*)
            el = el[[1]]]];
      If[First[komp] == Last[komp], komp = Drop[komp, -1]];
      res = Append[res, komp];
      res = UnsortedComplement[res, {{}}];
      el = fPrvi[inc];
      If[el != 0, komp = {el};
        (*posto pocinjemo novu komponentu vracamo inc na nekoristeno stanje*)
          red = Position[inc, el];
        red = Select[red, Length[#] == 2 && #[[2]] != 1 &][[1]];
        inc = ReplacePart[inc, {inc[[red[[1]], red[[2]]]], 1}, red];
        (*stavljamo daje iskoristen*)
        brP = fBrojPojave[inc, red[[1]], red[[2]]];
        red = red[[1]];]];(*kraj while*)res = fSrediParnost[res];
    res = fDowfromGaussExt[res];
    res = {res[[1]], Flatten[res[[2]]]};
    res=MinDowProjAltKL[res];
    res]
    
    fCorr[Ul_List] := Module[{pp, pp1, i},
    pp = Map[First, Ul];
    pp1 = Length[Union[pp]];
    pp1 = Flatten[Map[First, Table[Position[pp, i], {i, pp1}]]];
    pp = Table[Ul[[pp1[[i]]]], {i, Length[pp1]}];
    pp]
    
    fDowCodes[n_Integer] := Module[{vv, vv1, vvf, kk, kk1, mm, i, j, k},
    Run["plantri -apc2 " <> ToString[n] <> " " <> ToString[n] <> ".txt"];
    vv = Import[ToString[n] <> ".txt"];
    DeleteFile[ToString[n] <> ".txt"];
   vv="{"<>StringReplace[vv,{ToString[n]<>" "->"","\n"->",",","->","}]<>"}";
    vv = Partition[Map[ToCharacterCode, Map[ToString, ToExpression[vv]]] - \
96,
         n];
    vv = Union[
        Table[If[Max[Map[Length, vv[[i]]]] <= 4 , vv[[i]], {}], {i, 
            Length[vv]}]];
    vv = Partition[Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1], n];
    (* lista svih sa valencom <= 4 *)
    vv1 = 
      Table[Union[
          Map[Sort, 
            Union[Table[
                Table[Length[vv[[k, vv[[k, j, i]]]]], {i, 
                    Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}]]]], {k, 
          Length[vv]}];
    vv = Union[
        Table[If[
            SameQ[Select[
                vv1[[i]], (Length[#] == 2 && MemberQ[#, 4]) && 
                    Union[Map[EvenQ, #]][[1]] &], {}], vv[[i]], {}], {i, 
            Length[vv1]}]];
    vv = Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1];
    vv = Partition[
        Table[If[
            Length[vv[[i]]] == 2, {vv[[i, 1]], vv[[i, 1]], vv[[i, 2]], 
              vv[[i, 2]]}, vv[[i]]], {i, Length[vv]}], n];
    (* ovde treba napuniti odgovarajuca trovalentna temena *)
    vv1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
    vv1 = Map[Max, Table[Map[Length, vv1[[i]]], {i, Length[vv1]}]];
    vv = Union[Table[If[vv1[[i]] > 4, {}, vv[[i]]], {i, Length[vv1]}]];
    vv = If[SameQ[vv[[1]], {}], Drop[vv, 1], vv];
    (* svi obojivi sa dvojnim vezama samo u dvovalentnim temenima *)
    vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
    (* vvf su vec obojene *)
    vv = Complement[vv, vvf];
    vv1 = 
      Table[Table[
          If[vv[[j, i, 1]] == 
              vv[[j, i, 2]], {{vv[[j, i, 1]], i}, {vv[[j, i, 3]], 
                i}}, {}], {i, Length[vv[[j]]]}], {j, Length[vv]}];
    vv1 = Map[Sort, Table[Flatten[vv1[[i]], 1], {i, Length[vv1]}]];
    kk = Table[Complement[Range[n], Map[First, vv1[[i]]]], {i, Length[vv1]}];
    kk = Table[
        Table[{kk[[j, i]], 0}, {i, Length[kk[[j]]]}], {j, Length[kk]}];
    vv1 = Table[Union[vv1[[i]], kk[[i]]], {i, Length[kk]}];
    vv1 = Map[fCorr, vv1];
    vv = Table[
        Table[If[vv1[[j, i, 2]] != 0 && vv[[j, i, 1]] != vv[[j, i, 2]], 
            Insert[vv[[j, i]], vv1[[j, i, 2]], 
              Position[vv[[j, i]], vv1[[j, i, 2]]]], vv[[j, i]]], {i, n}], \
{j,
           Length[vv]}];
    vv = Union[vv, vvf];
    vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
    (* vec obojene koje dodajemo na kraju *)
    vv = Complement[vv, vvf];
    kk = Table[Apply[Plus, Map[Length, vv[[i]]]], {i, Length[vv]}];
    kk = Table[Divide[Times[4, n] - kk[[i]], 2], {i, Length[kk]}];
    vv1 = 
      Table[Union[
          Flatten[Table[
              Table[{{j, vv[[k, j, i]]}, 
              Length[vv[[k, vv[[k, j, i]]]]]}, {i, 
                  Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}], 1]], {k, 
          Length[vv]}];
    vv1 = 
      Table[Union[Map[Sort, Map[First, Select[vv1[[i]], #[[2]] == 3 &]]]], \
{i,
           Length[vv]}];
    kk1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
    kk1 = 
      Map[Union, 
        Map[Flatten, 
          Table[Select[kk1[[i]], Length[#] < 4 &], {i, Length[kk1]}]]];
    kk1 = 
      Table[Select[vv1[[i]], Length[Intersection[#, kk1[[i]]]] == 2 &], {i, 
          Length[vv]}];
    kk1 = Table[KSubsets[kk1[[i]], kk[[i]]], {i, Length[kk]}];
    kk1 = 
      Table[kk = 
          Union[Table[
              If[SameQ[Length[Flatten[kk1[[j, i]]]], 
                  Length[Union[Flatten[kk1[[j, i]]]]]], 
                  kk1[[j, i]], {}], {i, 
                Length[kk1[[j]]]}]];
        kk = 
          If[Not[SameQ[kk, {}]], If[SameQ[kk[[1]], {}], Drop[kk, 1], kk], 
            kk],
        {j, Length[kk1]}];
    vv1 = 
      Table[Table[
          Union[kk1[[j, i]], Map[Reverse, kk1[[j, i]]]], {i, 
            Length[kk1[[j]]]}], {j, Length[kk1]}];
    (* dobili smo sve kombinacije ivica koje dodajemo da bismo napravili 4 - 
        valentan graf *)
    mm = Table[Table[
          kk = vv1[[k, j]];
          
          kk1 = Table[
              Position[vv[[k, kk[[i, 1]]]], kk[[i, 2]]], {i, Length[kk]}];
          
          Table[{kk[[i, 1]], 
              Insert[vv[[k, kk[[i, 1]]]], kk[[i, 2]], kk1[[i]]]}, {i, 
              Length[kk]}],
          {j, Length[vv1[[k]]]}], {k, Length[vv1]}];
    vv1 = 
      Table[Table[{i, vv[[j, i]]}, {i, Length[vv[[j]]]}], {j, Length[vv]}];
    vv1 = Table[Select[vv1[[i]], Length[#[[2]]] == 4 &], {i, Length[vv1]}];
    vv1 = 
      Table[Table[Union[mm[[j, i]], vv1[[j]]], {i, Length[mm[[j]]]}], {j, 
          Length[mm]}];
    vv1 = 
      Table[Table[Map[Last, vv1[[j, i]]], {i, Length[vv1[[j]]]}], {j, 
          Length[vv1]}];
    vv1 = Flatten[Complement[vv1, {{}}], 1];
    vv1 = Union[vv1, vvf];
    vv1 = 
      Table[Table[Flatten[{i, vv1[[j, i]]}], {i, Length[vv1[[j]]]}], {j, 
          Length[vv1]}];
    vv1 = Union[Map[fDowkerCode, vv1]];
    vv1=Join[vv1,{Length[vv1]}];
    vv1
    ]   
    
    
    
    
    
(* ## ## ## ## ## ## ## ## ## # fGenKL ## ## ## ## ## ## ## ## ## ## ## *)


fBasic[n_Integer] := Module[{del, del1, i, j, k},
    If[n < 4, del = {}, del = compo[n];
      del = Select[del, Not[MemberQ[#, 1]] &];
      del = Select[del, Length[#] >= 2 &];
      del1 = Sort[Table[{Sort[del[[i]]], del[[i]]}, {i, Length[del]}]];
      del = Union[Table[Sort[del[[i]]], {i, Length[del]}]];
      del = Table[Select[del1, #[[1]] == del[[i]] &], {i, Length[del]}];
      del = 
        Table[Table[del[[j, i, 2]], {i, Length[del[[j]]]}], {j, 
            Length[del]}];
      del = 
        Map[Last, 
          Flatten[Table[
              Union[Table[
                  Union[Table[
                      RotateLeft[del[[k, j]], i], {i, Length[del[[k, j]]]}], 
                    Map[Reverse, 
                      Table[RotateLeft[del[[k, j]], i], {i, 
                          Length[del[[k, j]]]}]]], {j, 
                    Length[del[[k]]]}]], {k, Length[del]}], 1]]];
    del = {del};
    del = 
      Flatten[Table[
          Map[StringReplace[ToString[#], 
          {"{" -> "", "}" -> "", " " -> ""}] 
&,
             del[[i]]], {i, Length[del]}]];
    del]
    
fRtan[n_Integer] := Module[{stel, pp, ll, ll1, i, j},
    stel = fBasic[n];
    stel = 
      Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], Number] &, 
        stel];
    stel = Map[fMakePart[#] &, stel];(*svaki strem veci od 1 partitionira :*)

    
    stel = Flatten[Map[fPerm[#] &, stel], 1];
    stel = 
      Table[Table[
          ToExpression[StringReplace[stel[[j, i]], {" " -> ","}]], {i, 
            Length[stel[[j]]]}], {j, Length[stel]}];
    stel = 
      Map[Last, 
        Union[Table[
            Union[Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}], 
              Map[Reverse, 
                Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}]]], \
{j,
               Length[stel]}]]];
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel = Table[fStelString[stel[[i]]], {i, Length[stel]}];
    stel = Flatten[Table[Permutations[stel[[i]]], {i, Length[stel]}], 1];
    stel = Table[Map[ToString, stel[[i]]], {i, Length[stel]}];
    stel = 
      Table[StringReplace[stel[[i]], {"{" -> "", "}" -> "", ", " -> " "}], \
{i,
           Length[stel]}];
    stel = Map[Reverse, Union[Map[Sort, stel]]];
    pp = Map[Sort, Map[Split, Map[Sort, stel]]];
    ll = Map[Reverse, 
        Map[Sort, Table[Map[Length, pp[[i]]], {i, Length[pp]}]]];
    ll1 = Union[ll];
    pp = Map[Flatten, Table[Position[ll, ll1[[i]]], {i, Length[ll1]}]];
    stel = 
      Table[Table[stel[[pp[[j, i]]]], {i, Length[pp[[j]]]}], {j, 
          Length[pp]}];
    stel = Table[Map[Split, stel[[i]]], {i, Length[stel]}];
    stel = Table[Map[Sort, stel[[i]]], {i, Length[stel]}];
    ll = Table[Map[Length, stel[[i, 1]]], {i, Length[stel]}];
    stel = Table[Map[Flatten, stel[[i]]], {i, Length[stel]}];
    stel = {stel, ll};
    stel]

fRepres[n_Integer, pp_String] := 
  Module[{vv, ll, ll1, t, n1, stel, ss, i,j, ff, ff1, ppp},
    vv = fRtan[n];
    ll = vv[[2]];
    ll1 = 
      Map[Flatten, 
        Table[Table[
            Table[Length[ll[[j]]] + 1, {i, ll[[j, i]]}] - i + 1, {i, 
              Length[ll[[j]]]}], {j, Length[ll]}]];
    ll = Map[Length, ll1];
    t = iteratedTake[Map[ToString, Flatten[ll1]], ll];
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    n1 = Length[ppp];
    ff = Map[Length, t];
    t = Select[t, Length[#] == n1 &];
    stel = Table[Flatten[Map[Permutations, {t[[i]]}], 1], {i, Length[t]}];
    ss = Table[
        Flatten[Table[
            StringReplacePart[pp, Flatten[stel[[j, i]]], ppp], {i, 
              Length[stel[[j]]]}]], {j, Length[stel]}];
    ff = Table[Map[Abs, Map[MinDowAltKL, ss[[i]]]], {i, Length[ss]}];
    ff1 = Table[Union[ff[[i]]], {i, Length[ff]}];
    ff1 = 
      Map[Union, 
        Table[Table[Position[ff, ff[[j, i]]], {i, Length[ff[[j]]]}], {j, 
            Length[ff]}]];
    ff1 = Sort[Map[First, Flatten[ff1, 1]]];
    ff1 = Table[ss[[ff1[[i, 1]], ff1[[i, 2]]]], {i, Length[ff1]}];
    ff1
    ]


fDerive[b_Integer, pp_String] := 
  Module[{vv, ll, ll1, t, n1, stel, ss, ff, i, ff1, ppp, n,  j},
    n = fAdjustNo[b, pp];
    vv = fRtan[n];
    ll = vv[[2]];
    ll1 = 
      Map[Flatten, 
        Table[Table[
            Table[Length[ll[[j]]] + 1, {i, ll[[j, i]]}] - i + 1, {i, 
              Length[ll[[j]]]}], {j, Length[ll]}]];
    ll = Map[Length, ll1];
    t = iteratedTake[Map[ToString, Flatten[ll1]], ll];
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    n1 = Length[ppp];
    ff = Map[Length, t];
    t = Select[t, Length[#] == n1 &];
    stel = Table[Flatten[Map[Permutations, {t[[i]]}], 1], {i, Length[t]}];
    ss = Table[
        Flatten[Table[
            StringReplacePart[pp, Flatten[stel[[j, i]]], ppp], {i, 
              Length[stel[[j]]]}]], {j, Length[stel]}];
    ff = Table[Map[Abs, Map[MinDowAltKL, ss[[i]]]], {i, Length[ss]}];
    ff1 = Table[Union[ff[[i]]], {i, Length[ff]}];
    ff1 = 
      Map[Union, 
        Table[Table[Position[ff, ff[[j, i]]], {i, Length[ff[[j]]]}], {j, 
            Length[ff]}]];
    ff1 = Sort[Map[First, Flatten[ff1, 1]]];
    ff1 = Table[ss[[ff1[[i, 1]], ff1[[i, 2]]]], {i, Length[ff1]}];
    ll = ppp;
    ss = Table[
        Table[StringTake[ff1[[j]], ll[[i]]], {i, Length[ll]}], {j, 
          Length[ff1]}];
    ll = Map[Split, Map[Reverse, Map[Sort, ss]]];
    ll1 = Map[Flatten, ll];
    ll = Table[Map[Length, ll[[i]]], {i, Length[ll]}];
    ll = Flatten[Table[Position[vv[[2]], ll[[i]]], {i, Length[ll]}]];
    ll = Table[vv[[1, ll[[i]]]], {i, Length[ll]}];
    ll1 = 
      Map[Flatten, 
        Table[Map[First, 
            Table[Position[ll1[[j]], ss[[j, i]]], 
            {i, Length[ss[[j]]]}]], {j, 
            Length[ss]}], 1];
    ll1 = 
      Flatten[Table[
          Table[Table[ll[[k, j, ll1[[k, i]]]], {i, Length[ll1[[k]]]}], {j, 
              Length[ll[[k]]]}], {k, Length[ll1]}], 1];
    ll1 = Table[StringReplacePart[pp, ll1[[i]], ppp], {i, Length[ll1]}];
    ll1
    ]

fAdjustNo[n_Integer, pp_String] := Module[{p, i, ppp},
    p = Apply[Plus, fCreatePData[pp][[1]]];
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    n - p + 2Length[ppp]]

fGenKL[b_Integer, pp_String] := Module[{kk, n, res = {}, ppp, i},
    ppp = StringPosition[pp, "2"];
    If[ppp != {} && StringPosition[pp, "*"] != {},
      i = 1;
      While[ppp != {} && i <= Length[ppp],
        If [ppp[[i, 1]] < StringPosition[pp, "*"][[1, 1]], ppp = Rest[ppp], 
          i++]
        ]];
    If[Length[ppp] == 1,
            (*Only one digon *)
            n = fAdjustNo[b, pp]; ppp = Union[Flatten[ppp]];
           kk = Select[compo[n], Not[First[#] == 1] &];
            
      kk = Map[StringJoin, Table[Table[StringJoin[ToString[kk[[j, i]]], " "],
              {i, Length[kk[[j]]]}], {j, Length[kk]}]];
            kk = Table[StringDrop[kk[[i]], -1], {i, Length[kk]}];
        Do[ n = StringReplacePart[pp, kk[[i]], {ppp[[1]], ppp[[1]]}];
        res = Append[res, n],
        {i, Length[kk]}],
      (*if it has more than one digon*)
      res = fDerive[b, pp]];
      res=Join[res,{Length[res]}];
    res]



(*## ## ## ## ## ## ## ## ## ## # fMinDowKnot ## ## ## ## ## ## ## ##*)

(*## ## fMinDowKnotPr ## ## ## ## ## ## ##*)
(*# Minimizes Dowker codes for KNOTs only and 
works with single knot projection #*)
(*# Input List- Dowker Code with signs
                {{Lengths of Components},{code}}
                ConwaySymbol_ String
    Output Minimized Dowker Code #*)
(*# Needs:fGaussExt,FMinComp,fKnittDow,fFormDow #*)
fMinDowKnotPr[Unos_]:=Module[{DL,DowM},DL=fGaussExtSigns[Unos];
    If[SameQ[DL,Flatten[DL]],
      (*Knot*)DowM=fMinComp[DL];(*Print["MC ",DowM];*)
      DowM=Map[fKnittDow[DL,Take[#,-2]]&,DowM];
      (*Print["MKnitt ",DowM];*)
      DowM=Map[fFormDow[#[[1]],#[[2]]]&,DowM];
      DowM=Sort[DowM,Abs[#1]<Abs[#2]&],
      Print["Link"]];
      (* Ako ne zelimo orijentaciju *)
      If[Negative[First[DowM[[1]]]],DowM=-DowM,DowM=DowM];
  {{Length[First[DowM]]}, First[DowM]}
    ]

(*## ## fMinDowKnot ## ## ## ## ## ## ##*)
(*# Minimizes Dowker codes for KNOTs #*)
(*# Input ConwaySymbol_ String Output Minimized Dowker Code #*)
(*# Needs:fMinDowPr,fGaussExt,fMinComp,fKnittDow,fFormDow,fProjections #*)
fMinDowKnot[ConUl_String]:=Module[{ConPr,pom,m={}, \
res},ConPr=fProjections[ConUl];
    (* ConPr=Map[{Abs[fMinDowKnotPr[#]],#}&,ConPr]; *)
    ConPr=Map[{fMinDowKnotPr[#],#}&,ConPr];
    pom=Union[Abs[Flatten[Union[Map[Take[#,1]&,ConPr]],1]]];
   (*Biramo one koje imaju razlicite Dow bez obzira na projekciju*)
     Do [m=Append[m,
              Select[ConPr,SameQ[Abs[#[[1]]],pom[[i]]]&][[1]]]
         ,{i,1,Length[pom]}];
        ConPr=m;
        {Table[ConPr[[i]][[1]][[2]],{i,Length[ConPr]}],
        Table[ConPr[[i]][[2]],{i,Length[ConPr]}],
        res=First[Table[ConPr[[i]][[1]],{i,Length[ConPr]}]]};
        res={res[[1]],Flatten[res[[2]]]}
    ]



(*## ## ## ## ## ## ## # fKLinGraph  ## ## ## ## ## ## ## ## ## ## ## ## *)

(* fKLinGraph: Ulaz- lista unorderud pairs; 
  Izlaz:  neizomorfni KL u datom grafu. 
 Ako je dat graf L kao lista unordered pairs:
  1) nalazimo sve njegove podskupove duzine
   4<=d<= Length[L]; 
  2) proveravamo cetvorovalentnost- svaki broj 
  se javlja po 4 puta u Flatten delu;
  3) proveravamo planarnost izdvojenog 4-
  valentnog podgrafa PlanarQ; 
  4) 4-valentne planarne selektujemo na osnovu
   izomorfizma- biramo neizomorfne;
  5) medju njima biramo prime. *)
(*valenca temena *)

fVerVal[LL_List,el_Integer]:=Module[{res},res=Count[Flatten[LL],el];res]
(* da li je graf cetvorovalentan*)

fFourValQ[GG_List]:=Module[{l=Union[Flatten[GG]]},
    Union[Map[fVerVal[GG,#]&,l]]=={4}]
(* prvi iz klase neizomorfnih u listi grafova istih duzina*)

fClassRep[LL_List]:=Module[{l,ll,i,lP},
    ll=Map[FromUnorderedPairs[#]&,LL];lP={First[LL]};
       If[ll!={},l={First[ll]}];ShowLabeledGraph[ll[[1]]];
    ll=Rest[ll];
    While[ll!={},
      Do[ 
        If[IsomorphicQ [Last[l],ll[[i]]],ll=ReplacePart[ll,0,i]],{i,1,
          Length[ll]}];
      ll=Complement[ll,{0}];
      If[ll!={},l=Append[l,ll[[1]]];
        lP=Append[lP,ToUnorderedPairs[ll[[1]]]];
       ShowLabeledGraph[ll[[1]]]];
      If[ll!={},ll=Rest[ll]]
      ];
    lP
    ]
(*trazi cvorove i linkove u grafu*)

fKLinGraph[UOPair_List]:=Module[{l=Flatten[UOPair],val,res,i,svi={}},
    UnOrPair=UOPair;
    val=Union[Map[fVerVal[l,#]&,l]];
    If[Select[val,#>3&]=={},
          (*ne moze se naci u manje od 3-valentnom grafu*)
              
      Print["No 4-valent subgraphs"];res=0,
            Print["Input Graph"];
      	 ShowLabeledGraph[FromUnorderedPairs[UnOrPair]];
         Print["4-valent subgraphs"];
            Do[svi=Append[svi,KSubsets[UnOrPair,i]] ,{i,4,Length[UnOrPair]}];
              
      svi=Complement[Map[Select[#, fFourValQ[Flatten[#]]&]&,svi],{{}}];
              svi=Map[fClassRep[#]&,svi];
                res=Select[svi,PlanarQ[FromUnorderedPairs[#]]&]
      ];
    res=Join[res,{Length[res]}];
    res
    ]

(*## ## ## ## ## ## ## # fAddDig  ## ## ## ## ## ## ## ## ## ## ## ## *)

fAddDig[Ul_List] := Module[{vv, vv1, vvf, kk, kk1, mm, i, j, k},
    vv = Union[Ul, Map[Reverse, Ul]];
    n = Length[Union[Flatten[Ul]]];
    vv = Table[Select[vv, #[[1]] == i &], {i, n}];
    vv = Table[Map[Last, vv[[i]]], {i, Length[vv]}];
    vv = {If[Last[Union[Map[Length, vv]]] <= 4, vv, {}]};
    vv1 = If[Not[SameQ[vv, {{}}]],
        vv = 
          Union[Table[
              If[Max[Map[Length, vv[[i]]]] <= 4 , vv[[i]], {}], {i, 
                Length[vv]}]];
        vv = 
          Partition[Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1], n];
        (* lista svih sa valencom <= 4 *)
        vv1 = 
          Table[Union[
              Map[Sort, 
                Union[Table[
                    Table[Length[vv[[k, vv[[k, j, i]]]]], {i, 
                        Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}]]]], {k, 
              Length[vv]}];
        vv = 
          Union[Table[
              If[SameQ[
                  Select[vv1[[i]], (Length[#] == 2 && MemberQ[#, 4]) && 
                        Union[Map[EvenQ, #]][[1]] &], {}], vv[[i]], {}], {i, 
                Length[vv1]}]];
        vv = Flatten[If[SameQ[vv[[1]], {}], Drop[vv, 1], vv], 1];
        vv = 
          Partition[
            Table[If[
                Length[vv[[i]]] == 2, {vv[[i, 1]], vv[[i, 1]], vv[[i, 2]], 
                  vv[[i, 2]]}, vv[[i]]], {i, Length[vv]}], n];
        (* ovde treba napuniti odgovarajuca trovalentna temena *)
        vv1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
        vv1 = Map[Max, Table[Map[Length, vv1[[i]]], {i, Length[vv1]}]];
        vv = Union[Table[If[vv1[[i]] > 4, {}, vv[[i]]], {i, Length[vv1]}]];
        vv = If[SameQ[vv[[1]], {}], Drop[vv, 1], vv];
        (* svi obojivi sa dvojnim vezama samo u dvovalentnim temenima *)
        vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
        (* vvf su vec obojene *)
        vv = Complement[vv, vvf];
        vv1 = 
          Table[Table[
              If[vv[[j, i, 1]] == 
                  vv[[j, i, 2]], {{vv[[j, i, 1]], i}, {vv[[j, i, 3]], 
                    i}}, {}], {i, Length[vv[[j]]]}], {j, Length[vv]}];
        vv1 = Map[Sort, Table[Flatten[vv1[[i]], 1], {i, Length[vv1]}]];
        kk = 
          Table[Complement[Range[n], Map[First, vv1[[i]]]], {i, 
              Length[vv1]}];
        kk = 
          Table[Table[{kk[[j, i]], 0}, {i, Length[kk[[j]]]}], {j, 
              Length[kk]}];
        vv1 = Table[Union[vv1[[i]], kk[[i]]], {i, Length[kk]}];
        vv1 = Map[fCorr, vv1];
        vv = 
          Table[Table[
              If[vv1[[j, i, 2]] != 0 && vv[[j, i, 1]] != vv[[j, i, 2]], 
                Insert[vv[[j, i]], vv1[[j, i, 2]], 
                  Position[vv[[j, i]], vv1[[j, i, 2]]]], vv[[j, i]]], {i, 
                n}], {j, Length[vv]}];
        vv = Union[vv, vvf];
        vvf = Select[vv, SameQ[Union[Map[Length, #]], {4}] &];
        (* vec obojene koje dodajemo na kraju *)
        vv = Complement[vv, vvf];
        kk = Table[Apply[Plus, Map[Length, vv[[i]]]], {i, Length[vv]}];
        kk = Table[Divide[Times[4, n] - kk[[i]], 2], {i, Length[kk]}];
        vv1 = 
          Table[Union[
              Flatten[Table[
                  Table[{{j, vv[[k, j, i]]}, 
                      Length[vv[[k, vv[[k, j, i]]]]]}, {i, 
                      Length[vv[[k, j]]]}], {j, Length[vv[[k]]]}], 1]], {k, 
              Length[vv]}];
        vv1 = 
          Table[Union[
              Map[Sort, Map[First, Select[vv1[[i]], #[[2]] == 3 &]]]], {i, 
              Length[vv]}];
        kk1 = Map[Split, Map[Sort, Map[Flatten, vv]]];
        kk1 = 
          Map[Union, 
            Map[Flatten, 
              Table[Select[kk1[[i]], Length[#] < 4 &], {i, Length[kk1]}]]];
        kk1 = 
          Table[Select[vv1[[i]], 
              Length[Intersection[#, kk1[[i]]]] == 2 &], {i, Length[vv]}];
        kk1 = Table[KSubsets[kk1[[i]], kk[[i]]], {i, Length[kk]}];
        kk1 = 
          Table[kk = 
              Union[Table[
                  If[SameQ[Length[Flatten[kk1[[j, i]]]], 
                      Length[Union[Flatten[kk1[[j, i]]]]]], 
                    kk1[[j, i]], {}], {i, Length[kk1[[j]]]}]];
            
            kk = If[Not[SameQ[kk, {}]], 
                If[SameQ[kk[[1]], {}], Drop[kk, 1], kk], kk],
            {j, Length[kk1]}];
        vv1 = 
          Table[Table[
              Union[kk1[[j, i]], Map[Reverse, kk1[[j, i]]]], {i, 
                Length[kk1[[j]]]}], {j, Length[kk1]}];
             mm = Table[Table[
              kk = vv1[[k, j]];
              
              kk1 = Table[
                  Position[vv[[k, kk[[i, 1]]]], kk[[i, 2]]], {i, 
                    Length[kk]}];
              
              Table[{kk[[i, 1]], 
                  Insert[vv[[k, kk[[i, 1]]]], kk[[i, 2]], kk1[[i]]]}, {i, 
                  Length[kk]}],
              {j, Length[vv1[[k]]]}], {k, Length[vv1]}];
        vv1 = 
          Table[Table[{i, vv[[j, i]]}, {i, Length[vv[[j]]]}], {j, 
              Length[vv]}];
        vv1 = 
          Table[Select[vv1[[i]], Length[#[[2]]] == 4 &], {i, Length[vv1]}];
        vv1 = 
          Table[Table[Union[mm[[j, i]], vv1[[j]]], 
          {i, Length[mm[[j]]]}], {j, 
              Length[mm]}];
        vv1 = 
          Table[Table[Map[Last, vv1[[j, i]]], {i, Length[vv1[[j]]]}], {j, 
              Length[vv1]}];
        vv1 = Flatten[Complement[vv1, {{}}], 1];
        vv1 = Union[vv1, vvf];
        vv1 = 
          Table[Table[Flatten[{i, vv1[[j, i]]}], {i, Length[vv1[[j]]]}], {j, 
              Length[vv1]}], {}];
    vv1 = 
      Table[Flatten[
          Table[Table[{vv1[[k, j, 1]], vv1[[k, j, i]]}, {i, 2, 
                Length[vv1[[k, j]]]}], {j, Length[vv1[[k]]]}], 1], {k, 
          Length[vv1]}];
    vv1 = Table[Select[vv1[[i]], #[[1]] < #[[2]] &], {i, Length[vv1]}];
    If[Or[SameQ[vv1, {}], Length[vv1] == 1], vv1, vv1 = fClassRep[vv1]];
    vv1=Join[vv1,{Length[vv1]}];
    vv1
    ]
 
(*## ## ## ## ## ## ## # RECOGNIZE CONWAY  ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## *)

(*Radi od graph INc*)
fFindBP[gr_List]:=Module[{bp=gr,q,d=2},
    If [Select[gr,Count[gr,#]==2&]!={},
      While[d>=2,
        bp=Sort[Flatten[FromUnorderedPairs[bp][[1]],1]];
        q=Union[Select[bp,Count[bp,#]>=2&]][[1]];
        bp=Sort[Flatten[Contract[FromUnorderedPairs[bp],q][[1]],1]];
        d=
          If[SameQ[bp,{}],0,
            Max[Union[
                Table[Length[Flatten[Position[bp,bp[[i]]]]],{i,
                    Length[bp]}]]]]]];
    If[SameQ[bp,{}],1,bp]
    ]


fBasicPoly[Ul_List]:=Module[{pdata,bp,i=1,str="*",pom,duz},
    If[SameQ[Ul[[1]],{1,1}],bp="1*",
      If[MemberQ[Union[Map[OddQ,Ul[[2]]]],True],
      pdata=fPDataFromDow[Abs[Ul]],
        pdata=Ul];
     (*  pdata=ReductionKnotLink[pdata];*)
      bp=fDowfromPD[pdata];
      bp=fGraphInc[bp][[1]];
      bp=fFindBP[bp];
      If[Not[SameQ[Head[bp],List]],
        bp="1*",(*poliedarski*)
        pom=FromUnorderedPairs[bp];
        If[Max[Flatten[bp]]<=12,
          While[Not[IsomorphicQ[pom,FromUnorderedPairs[GraphBP[[i,2]]]]],
            i++];
          bp=GraphBP[[i,1]],str=ToString[Max[Flatten[bp]]]<>ToString[i]<>str;
          duz=StringLength[str]-2;          
          While[Not[IsomorphicQ[pom,
          FromUnorderedPairs[fGraphInc[str][[1]]]]],
            i++;str=StringTake[str,duz]<>ToString[i]<>"*"];
          bp=str] (*kraj If[Max....]*)]];
    (*kraj od If[Not[SameQ[Head[bp],List]]*)
    bp]

    
    
    fCompositePoly[Conway_String]:=Module[{p,q,q1, i},
    p=fGraphInc[Conway][[1]];
    q=KSubsets[p,4];
    q=Drop[
        Union[Table[
            If[Length[Union[Flatten[q[[i]]]]]>5,q[[i]],{}],{i,Length[q]}]],
        1];
    q=Table[Complement[p,q[[i]]],{i,Length[q]}];
    q1=Map[ConnectedComponents,Map[FromUnorderedPairs,q]];
    q1=Flatten[Position[Table[Length[q1[[i]]],{i,Length[q1]}],2]];
    q1=If[Not[SameQ[q1,{}]],1,0]]
    
    fPolyFlype[Conway_String]:=Module[{p},
    p=fGraphInc[Conway][[1]];
    p=First[
        Union[Table[
            EdgeConnectivity[DeleteVertex[FromUnorderedPairs[p],i]],{i,
              Length[p]/2}]]];
    p=Abs[3-p];
    p]
    
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
# *)
(* ## ## ## ## ## ## f FIND CON ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*radi od Liste Uredjenih parovakoj aima {Konvej, {lista njegivih \
neizomorfnih projekcija}} koji su grincovi
    i dovkera elementa kako se trazi*)

fPickCon[LL_List,DD_List]:=Module[{ind=False,p,i=1,j,res},
    While[ind!=True,
                  p=LL[[i,3]];j=1;
                  While[ind!=True&&j<=Length[p],
                                ind=IsomorphicQ[FromUnorderedPairs[p[[j]]],
                                                                       
            FromUnorderedPairs[fGraphInc[DD][[1]]]];
                               j++];
               If[ind==True,res=LL[[i,2]]];
             i++];
    res]
(*radi od Dowkera ili PData *)

fFindCon[Ul_]:=Module[{DD,i=1,n},
    If[SameQ[Head[Ul],String],   DD=fDToD[Ul];
      ];(*ako je Konvej*)
     If[SameQ[Head[Ul],List] ,
           If[Select[Ul[[2]],OddQ[#]&]=={},
                 (*ako je ulaz Dowker *) DD=Ul,
                                DD= fDowfromPD[ReductionKnotLink[Ul]]]];
    (*sada nam je DD bas DOWKER*)
    n=Length[DD[[2]]];
    (*broj presecnih tacaka *)
    bp= fBasicPoly[DD];
    (*bazicni poliedar *)
    (*sada nam je DD graph spreman za izomorfizam- 
        gde da ga trazimo*)
    n=ToExpression["GraphA"<>ToString[n]];
    (*suzavamo izbor na one sa kojima iam isti bazicni*)
    
    n=Select[n, #[[1]]==bp&];
    fPickCon[n,DD]
    ]
    
    
    
(* ## ## ## ## ## ## ## ## ## fFindConway ## ## ## ## ## ## ## ## ## ## ## *)


fRecCon[Ul_] := Module[{cc, mm, ss, res},
    cc = fComponentNo[Ul];
    res = 
      If[MemberQ[ReductionKnotLink[fCreatePData[Ul]], {}], 1, 
        mm = fConwayToPD[Ul];
        ss = DTCode[mm];
        ss = {{Length[ss]}, Table[ss[[i]], {i, Length[ss]}]};
        If[SameQ[cc, 1], fKnotFind[ss], 
          Sort[{RedKauffman[ReductionKnotLink[fCreatePData[Ul]]], 
              RedKauffman[
                ReductionKnotLink[GetMirrorImageKnot[fCreatePData[Ul]]]]}]]];
    res = {cc, res};
    res
    ]

fRecPdata[Ul_] := Module[{cc, res},
    cc = Length[Ul[[1]]];
    res = If[MemberQ[ReductionKnotLink[Ul], {}], 1,
        If[SameQ[cc, 1], 
          fKnotFind[fKnotscapeDowfromDow[fDowfromPD[ReductionKnotLink[Ul]]]], \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

          Sort[{RedKauffman[Ul], 
              RedKauffman[GetMirrorImageKnot[ReductionKnotLink[Ul]]]}]]];
    res = {cc, res};
    res
    ]

fRecKnotsc[Ul_] := Module[{cc, ss, res},
    cc = Length[Ul[[1]]];
    ss = ReductionKnotLink[fPDataFromDow[Ul]];
    res = If[MemberQ[ss, {}], 1,
        If[SameQ[cc, 1], fKnotFind[Ul], 
          Sort[{RedKauffman[ss], RedKauffman[GetMirrorImageKnot[ss]]}]]];
    res = {cc, res};
    res
    ]

fRecPD[Ul_] := Module[{cc, ss, ss1, res, i},
    ss = DTCode[Ul];
 cc = Length[fPDToPData[ss][[1]]];
    ss1 = ReductionKnotLink[fPDToPData[Ul]];
    res = 
      If[MemberQ[ss1, {}], 1, 
        If[NumberQ[ss[[1]]], 
          fKnotFind[{{Length[ss]}, Table[ss[[i]], {i, Length[ss]}]}], 
          Sort[{RedKauffman[ReductionKnotLink[fPDToPData[Ul]]], 
              RedKauffman[
                GetMirrorImageKnot[ReductionKnotLink[fPDToPData[Ul]]]]}]]];
    res = {cc, res};
    res]

fRecDow[Ul_] := Module[{cc, ss, res},
    cc = Length[Ul[[1]]];
    ss = ReductionKnotLink[fPDataFromDowker[Ul]];
    res = If[MemberQ[ss, {}], 1,
        If[SameQ[cc, 1], fKnotFind[fKnotscapeDowfromDow[Ul]], 
          Sort[{RedKauffman[ss], RedKauffman[GetMirrorImageKnot[ss]]}]]];
    res = {cc, res};
    res
    ]

fForFindConway[Ul_, type_Integer] := 
  Module[{ss, res1, res2, res3, res4, res5, res},
    res1 = If[SameQ[type, 1], fRecCon[Ul], {0, 0}];
    res2 = If[SameQ[type, 2], fRecPdata[Ul], {0, 0}];
    res3 = If[SameQ[type, 3], fRecKnotsc[Ul], {0, 0}];
    res4 = If[SameQ[type, 4], fRecPD[Ul], {0, 0}];
    res5 = If[SameQ[type, 5], fRecDow[Ul], {0, 0}];
    res = {res1, res2, res3, res4, res5};
    ss = Map[First, res];
    res = Select[res, Not[SameQ[#, {0, 0}]] &][[1]];
    res = If[SameQ[res[[2]], 1], 1, res];
    res
    ]


fFindConway[Ul_, type_Integer] := Module[{cc, mm, ss, vv, res},
    mm = fForFindConway[Ul, type];
res=If[SameQ[mm,1],{},
res=If[SameQ[Head[mm[[2,1,1]]],List],{},
res=If[SameQ[mm,1],{},
 res = If[SameQ[mm[[1]], 1],
        ss = 
          First[Flatten[Position[tableallknots[[mm[[2, 1, 1]]]], mm[[2]]]]];
        First[tableallknots[[mm[[2, 1, 1]]]][[ss]]], cc = mm[[1]];
        vv = Map[First, Position[tablealllinks[[cc - 1]], mm]];
        Map[First, 
          Table[tablealllinks[[cc - 1, vv[[i]]]], {i, Length[vv]}]]];
    res]]];
res]



fFindDirProd[Ul_]:=Module[{hh,hh1,hh2,ss1,res1,res2,ss,res,i,j},
hh=Drop[FactorList[RedKauffman[Ul]],1];
hh1=Map[fPolyNorm,Flatten[Table[Table[hh[[j,1]],{i,hh[[j,2]]}],{j,Length[hh]}]]];
hh2=Flatten[Table[ss=Union[Position[tablealllinks,hh1[[i]]],Position[tablealllinks,-hh1[[i]]]];
ss1=Union[Table[Drop[ss[[i]],-2],{i,Length[ss]}]],{i,Length[hh1]}],1];
res1=Table[First[tablealllinks[[hh2[[i,1]],hh2[[i,2]]]]],{i,Length[hh2]}];
hh=Drop[FactorList[RedKauffman[Ul]],1];
hh1=Map[fPolyNorm,Flatten[Table[Table[hh[[j,1]],{i,hh[[j,2]]}],{j,Length[hh]}]]];
hh2=Flatten[Table[ss=Union[Position[tableallknotspoly,hh1[[i]]],Position[tableallknotspoly,-hh1[[i]]]];
ss1=Union[Table[Drop[ss[[i]],-1],{i,Length[ss]}]],{i,Length[hh1]}],1];
res2=Table[First[tableallknotspoly[[hh2[[i,1]]]]],{i,Length[hh2]}];
res=Join[res1,res2];
res=StringDrop[StringJoin[Table[StringJoin[res[[i]],"#"],{i,Length[res]}]],-1]]


(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*## ## ## ## ## ## ## # SEIFERT MATRIX AND SIGNATURE  ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## *)

(* The author of the functions SeifertMatrix and ssmW is S.Orevkov *)
fSeifert[Ulaz_]:=Module[{brd,dd,sl,pd,sei,Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
      If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],Ul];
      If[SameQ[Head[Ul],List],
      If[Union[Map[OddQ,Abs[Ul[[2]]]]]=={False},dd=fSignsKL[Abs[Ul]];
        sl=Map[Sign,Ul[[2]]]*Map[Sign,dd[[2]]];
        pd=
          If[SameQ[Ul[[1]],
                dd[[1]]]&&(SameQ[Ul[[2]],dd[[2]]]||SameQ[Ul[[2]],-dd[[2]]]),
            Ul,fPDataFromDow[{Abs[Ul[[1]]],Abs[Ul[[2]]]*sl}]],pd=Ul],Ul];
   (*  bword=GetBraidRep[ReductionKnotLink[pd]]; *)
    bword=fBraidW[ReductionKnotLink[pd]];
    brd=fBrdFromBword[bword];
    sei=ssmW[brd[[1]],brd[[2]]];
    sei
    ]
    
fSignat[Ul_]:=Module[{dd},
    If[SameQ[Head[Ul],String],dd=fConwayToPD[Ul],Ul];
    If[SameQ[Head[Ul],List],
      If[Union[Map[OddQ,Abs[Ul[[2]]]]]=={False},dd=fDowkerToPD[Ul],
        dd=fPdataToPD[Ul]]];
    dd=Abs[KnotSignature[dd]];
    dd
    ]

fSignature[Ul_]:=Module[{dd},
    If[SameQ[Head[Ul],String],dd=fConwayToPD[Ul],Ul];
    If[SameQ[Head[Ul],List],
      If[Union[Map[OddQ,Abs[Ul[[2]]]]]=={False},dd=fDowkerToPD[Ul],
        dd=fPdataToPD[Ul]]];
    dd=KnotSignature[dd];
    dd
    ]
    
fSeifertJ[Ul_]:=Module[{res,str,f="input.txt"},
str=fBraidW[Ul];
OpenWrite[f];
WriteString[f,str];
Write[f];
Close[f];
Run["pretzelslavik.exe"<>" input.txt"<>" "<>"output.txt"];
res=Import["output.txt"];
DeleteFile["input.txt"];
DeleteFile["output.txt"];
res=ToExpression[StringReplace[res,{" "->"","]["->"},{","[["->"{{","]]"->"}}"}]];
res]

fSignatureJ[Ul_]:=Module[{mm,mm1,res},
mm=fSeifertJ[Ul];
mm1=Transpose[mm];
res=Plus@@Map[Sign,N[Eigenvalues[mm+mm1]]];
res]

fSeifertM[Ul_]:=Module[{ii,pp1,pp2,tt,res,l,DL,str,f="C://LinKnot//seifertMorwen//fff"},
    If[SameQ[Head[Ul],String],
      If[SameQ[StringPosition[Ul,"#"],{}],DL=Dowker[Ul][[4]]; 
        l=fGenSign[Ul]*fGenSign[StringReplace[Ul,{"-"->""}]];
        DL*=l,DL= fDToDDirect[Ul][[2]]], DL=Ul[[2]]];
    OpenWrite[f];
    str=ToString[Length[DL]]<>" 1 "<>" "<>
        StringReplace[ToString[DL],{"{"->"","}"->"",","->""}];
    WriteString[f,str];
    Write[f];
    Close[f];
 Run["C://LinKnot//seifertMorwen//seifert","C://LinKnot//seifertMorwen//fff", "C://LinKnot//seifertMorwen//output.txt"];
     ii= Import[".//seifertMorwen//output.txt"];
pp1=First[Union[Flatten[StringPosition[ii,"["]]]];
pp2=First[Union[Flatten[StringPosition[ii,"]"]]]];
tt=StringTake[ii,{pp1,pp2}];
res=ToExpression[StringReplace[tt,{" "->"","["->"{{","]"->"}}",";"->"},{"}]];
res]

fSignatureM[Ul_]:=Module[{mm,mm1,res},
mm=fSeifertM[Ul];
mm1=Transpose[mm];
res=Plus@@Map[Sign,N[Eigenvalues[mm+mm1]]];
res]

fDeterminantM[Ul_]:=Module[{mm,aa},
mm=fSeifertM[Ul];
aa=fPolyNorm[Det[mm-x*Transpose[mm]]];
aa=Abs[ReplaceAll[aa,x->-1]];
aa]
    
fBrdFromBword[bword_String]:=Module[{m,brd, i},
    m=Length[Union[ToCharacterCode[ToUpperCase[bword]]]]+1;
    brd=ToCharacterCode[bword]-64;
    brd=Table[
        If[brd[[i]]>26,brd[[i]]=-Mod[brd[[i]],32],brd[[i]]],{i,Length[brd]}];
    {m,brd}
    ]
SeifertMatrix[m_Integer,brd_List]:=
    Module[{n,e,V,X,q,c,i,j,h,a,b},
      a={{0,1,-1,0},{-1,0,1,0},{0,0,0,0},{1,-1,0,0}};
      b={{-1,1,0,0},{1,-1,0,0},{0,0,0,0},{0,0,0,0}};
      n=Length[brd];V=Table[0,{i,n},{j,n}];X=Table[{n},{h,m-1}];
      Do[h=Abs[brd[[q]]];e=Sign[brd[[q]]];
        c[1]=X[[h,1]];X[[h]]={c[2]=q};
        c[3]=If[h<m-1,X[[h+1,1]],n];c[4]=If[h>1,X[[h-1,1]],n];
        Do[Do[V[[c[i],c[j]]]+=a[[i,j]]+e*b[[i,j]],{i,4}],{j,4}],{q,n}];
      V=Delete[V,X];
      If[Length[V]>0,V=Transpose[Delete[Transpose[V],X]]/2];
      V];
      (*popravljena 15.08.2003 *)

ssmW[m_Integer,brd_List]:=
    Module[{bq,n,d,e,V,X,p=1,q,r=1,c,i,j,h,a,b},
      a={{0,0,-1,1,2},{0,0,1,-1,-2},{-1,1,0,0,0},{1,-1,0,0,0},{2,-2,0,0,0}};
      b={{-2,2,0,0},{2,-2,0,0},{0,0,0,0},{0,0,0,0}};
      d=n=Length[brd];Do[If[Not[IntegerQ[brd[[q]]]],d++;p++],{q,n}];
      V=Table[0,{i,d},{j,d}];X=Table[{d},{h,m-1}];
      Do[bq=brd[[q]];
        If[IntegerQ[bq],h=Abs[bq];e=Sign[bq],h=Abs[bq[[2]]];
          V[[r,r]]=2*bq[[1]]*Sign[bq[[2]]];c[5]=r++;];
        c[1]=X[[h,1]];X[[h]]={c[2]=p++};
        c[3]=If[h<m-1,X[[h+1,1]],d];c[4]=If[h>1,X[[h-1,1]],d];
        If[IntegerQ[bq],
          Do[Do[V[[c[i],c[j]]]+=a[[i,j]]+e*b[[i,j]],{i,4}],{j,4}],
          Do[Do[V[[c[i],c[j]]]+=a[[i,j]],{i,5}],{j,5}]],{q,n}];
       V=Delete[V,X];
      If[Length[V]>0,V=Transpose[Delete[Transpose[V],X]]/2];
      V];
      (*popravljena 15.08.2003 *)
      
(*fSignatureKL[S_List]calculates signature from a Seifert matrix S*)

fSignatureKL[S_List]:=Module[{s},s=S+Transpose[S];
    s=Abs[Apply[Plus,Map[Sign,Eigenvalues[N[S]]]]]]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##*)
(*## ## ## ## ## ## ## # SYMMETRY  ## ## ## ## ## ## ## ## ## ## ## ## *)
(*SymmKL[Ul] racuna automorfizme projekcije KL koji cuvaju znakove*)

Symm[Ulaz_]:= Module[{G,k,i,j, sG, Aut,Ul},
 If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    G = fGaussExt[Ul];
    sG = fGrInc[G][[2]];
    (*Print[sG];*)G = FromUnorderedPairs[fGrInc[G][[1]]];
    Aut = Automorphisms[G];
    (*Print[Aut];*)(*ovde izdvajamo automorfizme koji cuvaju znak*)
    Aut = 
      If[Aut != {}, 
        If[First[
              Union[Table[
                  If[SameQ[
                        Sign[Table[
                              Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]], {j, 
                                  Length[sG]}], {i, Length[Aut]}]][[k]], 
                        sG] == True, Aut[[k]], {}], {k, Length[Aut]}]]] == \
{},
           Drop[Union[
              Table[If[
                  SameQ[Sign[
                          Table[Table[
                              Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]], {j, 
                    Length[sG]}], {i, Length[Aut]}]][[k]], sG] == 
                    True, Aut[[k]], {}], {k, Length[Aut]}]], 1], 
          Union[Table[
              If[SameQ[
                    Sign[Table[
                          Table[Aut[[i]][[j]]*sG[[Aut[[i]][[j]]]], {j, 
                              Length[sG]}], {i, Length[Aut]}]][[k]], sG] == 
                  True, Aut[[k]], {}], {k, Length[Aut]}]]], Aut];
       {Aut,If[Aut != {}, Map[ToCycles, Aut], Aut],Length[Aut]}
    ]

(*## ## ## ## ## ## ## # SYMMETRY  ## ## ## ## ## ## ## ## ## ## ## ## *)
(* nalazi max simetricnu projekciju alt KL od konveja*)

MaxSymmProjAltKL[Con_String]:=Module[{l},
    l=fDiffProjectionsAltKL[Con];
    l=Map[#[[1]]&,l];
    l=Map[{Length[Symm[#][[1]]],#}&,l];
    l=Reverse[Sort[l]][[1]];
    {l[[2]],l[[1]]}
    (* ShowKnotfromPdata[fCreatePData[l[[2]]]]; *)
    ]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
#*)
(*## ## ## ## ## ## ## # UINFTYNO  ## ## ## ## ## ## ## ## ## ## ## ## *)

(*pomocna za fMirorDve *)
(*Revertuje znake jednostrukih elemenata*)
UnsortedComplement[x_List,y__List]:=
  Replace[x,Dispatch[(#:>Sequence[])&/@Union[y]],1]
fAddSign[LL_List,{el_Integer,elS_Integer}]:=Module[{z,i,p,l=LL},z=1-elS;
    p=Position[LL,{el,0}][[1,1]];
    Do[l=ReplacePart[l,{l[[i,1]],z},i];
      z=1-z,{i,p,Length[LL]}];
    z=1-elS;
    Do[z=1-z;(*Print[l];*)l=ReplacePart[l,{l[[p-i,1]],z},p-i],{i,1,p-1}];
    (*Print["rezultat iz fADDSIGN ",l];*)l]
fGiveSign[LL_List,SL_List]:=
  Module[{presek,p,i,l=LL,saznakom=SL,sz={},sa},If[SL=={},p=l[[1]];
      Do[p=ReplacePart[p,{p[[i,1]],Mod[i,2]},i];,{i,1,Length[p]}];
      l=ReplacePart[l,Prepend[p,1],1];
      sz=p,(*ako vec imamo neke sa znacimo trazimo ih u drugim komponentama*)
        sa=Union[Flatten[Map[Take[#,1]&,saznakom]]];
      (*Print["SSSSS",sa,saznakom];*)Do[
        If[Length[l[[i,1]]]!=0,(*znaci da ona nije oznacena*)presek=
            Select[l[[i]],MemberQ[sa,#[[1]]]&];
          
          If[presek!={},
            presek=Select[saznakom,#[[1]]==presek[[1,1]]&][[1]];
            p=fAddSign[l[[i]],presek];
            l=ReplacePart[l,Prepend[p,1],i];
            p=Map[If[MemberQ[sa,#[[1]]],-1,#]&,p];
            sz=Union[sz,Complement[p,{-1}]]]],{i,1,Length[l]}]];
    {l,sz}]
(*Pravi alternirajuci kad se polazi od pData...*)

fMakeAlt[CoLe_List,LL_List]:=
  Module[{i,pom=0,CoPom={},saznakom={},ind=True},l=iteratedTake[LL,CoLe];
    While[ind,(*Print["SZ",saznakom];*)saznakom=fGiveSign[l,saznakom];
      (*Print["posle fGiveSigna",saznakom[[1]],"\n \n",saznakom[[2]],
            "\n \n"];*)l=saznakom[[1]];
      saznakom=Last[saznakom];
      (*Print[Map[Length[#[[1]]]&,l]];*)ind=
        Union[Map[Length[#[[1]]]&,l]]!={0}];
    Map[Drop[#,1]&,l]]
(*pomocna za fMakeGaussPD-menja odnose iznad ispod za neparne PData*)
checkArgs[s_,t_]:=
  If[Head[s]===List&&VectorQ[t,Head[#]===Integer&&#>=0&]&&
      Plus@@t<=Length[s],True,False]

iteratedTake[s_,t_]/;checkArgs[s,t]:=
  First/@Rest[FoldList[Through[{Take,Drop}[#1[[2]],#2]]&,{{},s},t]]
fParno[LL_List]:=Module[{l,lGExt,i,p},l=LL;
    lGExt=Map[Flatten,l];
    lGExt=Map[Take[#,{1,Length[#],2}]&,lGExt];
    (*extGaus bez odnosa iznad ispod*)Do[
      p=Flatten[Position[Flatten[Take[lGExt,i-1]],lGExt[[i,1]]]];
      (*pozicija prvog u sledecoj komponenti u prethodnom delu gausa*)If[
        p!={},
        If[Not[EvenQ[p[[1]]]],
          lGExt=ReplacePart[lGExt,RotateLeft[lGExt[[i]]],i];
          l=ReplacePart[l,RotateLeft[l[[i]]],i]]],{i,2,Length[l]}];
    l]
(*ako se posle secenja za PD u istoj komponenti 
jave dva ista uzastopna broja 
(ciklicno) sa suprotnim drugim znakom npr.{-18,1},{-18,
      0} izbacujemo ih.*)
(*radi od gaus ext sa iznad ispod*)
fIzbaciPetlju[GK_List]:=Module[{el,p,l,pom,res=GK},l=Map[Take[#,1]&,GK];
    el=Union[Select[l,Count[l,#]==2&]];
    p=el;
    If[p!={},p=Map[Union[Flatten[Position[l,#]]]&,p];
      l=Map[If[#!={1,Length[GK]},RotateLeft[l,#[[2]]-1],l]&,p];
      p={};
      Do[p=Append[p,Flatten[Position[l[[i]],el[[i]]]]],{i,1,Length[el]}];
      l=GK;
      Do[If[p[[i]]=={1,Length[GK]},
          l=Select[l,#[[1]]!=el[[i,1]]&]],{i,1,Length[el]}];res=l];
    res]
(*## ## ## ## ## ## ## # UINFTYNO ## ## ## ## ## ## ## ## ## ## ## \
##*)(*pomocna za \
fMirorDve*)(*Revertuje znake jednostrukih elemenata*)
  fRevertSigns[ll_List,rev_List]:=Module[{p=rev,pom=ll},p=Map[Take[#,1]&,p];
    p=Flatten[Select[p,Count[p,#]==1&]];
    pom=Map[If[MemberQ[p,#[[1]]],{-#[[1]],#[[2]]},#]&,pom]]
(*pomocna za fMirorComp,fMiror1*)
(*izbaci jednu dvojnu tacku i da gaus EXt \
sa iznad ispod bez nje*)
(*radi od gausEXt i rednog broja komponent koju \
menjamoo i el koji izbacijemo*)

fMirorDve[ExtGaus_List,redni_Integer,el_Integer]:=
  Module[{pom=ExtGaus,len,komp,pos,p},komp=pom[[redni]];
    pos=Position[Map[Take[#,1]&,komp],el];
    pos=Map[#[[1]]&,pos];
    (*pozicije na kojima se dvojni nalazi*)p=
      Reverse[Take[komp,{pos[[1]]+1,pos[[2]]-1}]];
    (*p sadrzi deo komponente koji se revertovao*)komp=
      Join[Take[komp,pos[[1]]-1],p,Drop[komp,pos[[2]]]];
    pom=ReplacePart[pom,komp,redni];
    len=Map[Length[#]&,pom];
    pom=Select[Flatten[pom,1],#[[1]]!=el&];
    (*moze i bez prenumeracije*)pom=fRevertSigns[pom,p];
    pom=iteratedTake[pom,len];
    (*Print["PMiror DVE Pre Izbaci petlju  ",pom];
      pom=Map[fIzbaciPetlju[#]&,pom];Print["posle petlji ",pom];*)pom]
(*pomocna za fDvojneQ*)
(*vraca True ako ima samo jednostruke elemente*)

fDvojneKomp[Ulaz_List]:=Module[{res},
    res=Ulaz;
    res=Union[Map[Count[res,#]&,res]];
    If[res=={1},res=True,res=False];
    res](*pomocna za fMakeGaussPD-menja odnose iznad ispod za neparne PData*)
fMenjaObePojavePD[LL_List,el_Integer]:=
  Module[{p,i,l=LL},p=Select[l,#[[1]]==el&];
    (*Print[p,p[[1,1]]];*)(*Print[Position[l,p[[1,1]]]];*)p=
      Map[Take[#,1][[1]]&,Position[l,p[[1,1]]]];
    Do[l=ReplacePart[l,{l[[p[[i]],1]],1-l[[p[[i]],2]]},p[[i]]],{i,1,2}];
    l]
(*radi od Ext Gausa sa odnosima iznad ispod*)
(*ako moze sredi parnost*)
(*NE \
RADI RETROAKTIVNO!!!!*)

fMakeGaussPD[PD_List]:=Module[{l,n,i,l1={},nepar},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l];*)(*uredjeni parovi sa sve znacima*)nepar=
      Select[l,Or[OddQ[#[[2]]],And[EvenQ[#[[1]]],EvenQ[#[[2]]]]]&];
    nepar=Flatten[Map[Take[#,-1]&,nepar]];
    (*Print[nepar];*)l1=Map[#-#&,Range[2n]];
    Do[l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,1]]]];
      l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,2]]]],{i,1,n}];
    (*Print["???",2PD[[1]],
          l1];*)(*ubacili smo indikatore i sad treba da ga napravimo da je \
alternirajuci*)l1=fMakeAlt[2PD[[1]],l1];
    (*menjamo odnose onima koji su u pdata bili neparni a to su i ovde*)l1=
      Flatten[l1,1];
    (*Print[" pre kvarenja ",l1];*)If[nepar!={},
      Do[l1=fMenjaObePojavePD[l1,nepar[[i]]],{i,1,Length[nepar]}]];
    (*Print[" posle kvarenja ",l1];*)l1]
(*pomocna za fMirore*)
(*radi od podeljenog Gauss Ext sa odnosima iznad-ispod*)
fMakeGaussPD[PD_List]:=Module[{l,n,i,l1={},nepar},n=Length[PD[[2]]];
    l=Complement[Range[2n],Abs[PD[[2]]]];
    (*svi brojevi do 2 duzine koji vec nisu u pData*)Do[
      l=ReplacePart[l,{l[[i]]*Sign[PD[[2,i]]],PD[[2,i]]},i],{i,1,n}];
    (*Print["www",l];*)(*uredjeni parovi sa sve znacima*)nepar=
      Select[l,Or[OddQ[#[[2]]],And[EvenQ[#[[1]]],EvenQ[#[[2]]]]]&];
    nepar=Flatten[Map[Take[#,-1]&,nepar]];
    (*Print[nepar];*)l1=Map[#-#&,Range[2n]];
    Do[l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,1]]]];
      l1=ReplacePart[l1,{l[[i,2]],0},Abs[l[[i,2]]]],{i,1,n}];
    (*Print["???",2PD[[1]],
          l1];*)(*ubacili smo indikatore i sad treba da ga napravimo da je \
alternirajuci*)l1=fMakeAlt[2PD[[1]],l1];
    (*menjamo odnose onima koji su u pdata bili neparni a to su i ovde*)l1=
      Flatten[l1,1];
    (*Print[" pre kvarenja ",l1];*)If[nepar!={},
      Do[l1=fMenjaObePojavePD[l1,nepar[[i]]],{i,1,Length[nepar]}]];
    (*Print[" posle kvarenja ",l1];*)l1]
(*pomocna za fBrCo-izbacuje iz Gaus parova i zavisnosti od "ociscenog"Gausa*)
fDvojneQ[GausExt_List]:=Module[{res=GausExt,len},len=Map[Length[#]&,res];
    res=Flatten[Take[Flatten[res],{1,Length[Flatten[res]],2}]];
    res=iteratedTake[res,len];
    res=Map[fDvojneKomp[#]&,res];
    If[Union[res]=={True},res=True,res=False];
    res]
(*vraca True ako nijedna komponenta nema dvojnih*)
(*pomocna za fMirorComp*)
\
(*daje nove PData kad je izbacena 1 dvostruka*)
fMiror1[Ul_,redni_Integer,el2_Integer]:=
  Module[{pom,ind},pom=fMirorDve[Ul,redni,el2];
    (*sredjujemo parnost-
        ako ne moze da sredi dobicemo dva neparna i dva parna uparena*)If[
      pom=={{}},(*posle svega i izbacivanja ptelji doboili smo prazno-
          kraj*)pom={-1},pom=fParno[pom];
      pom=fPDataCuttNo[pom];
      pom=ReductionKnotLink[pom];
      If[MemberQ[{{{},{}},{{0},{}},{{},{0}}},pom]||
          fDvojneQ[iteratedTake[fMakeGaussPD[pom],2pom[[1]]]],pom={-1}]];
    (*ako u ovom razvezivanju nema dvojnih ili je unknott 
    onda je to kraj i 
vraca {-1}*)(*inace vraca neke pdata*)pom]
(*jedna komponenta->
    lista rezultata*)
(*Radi od pdata ili Konveja i broja komponente cije \
cemo dvostruke da izbacujemo*)
(*radi samo za komponente koje imaju dvojnih*)

(*vraca listu rezultata*)

fMirorKomp[Ul_,k_Integer]:=
  Module[{pd,pd0,pd1,p,i,dvojne,res={},brzi=True},
    If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    pd=ReductionKnotLink[pd];
    (*sad moramo da napravimo odgovarajuceg Gausa i da ga iz \
alterniramo*)pd0=
      fMakeGaussPD[pd];
    (*Gaus parovi sa iznad-ispod tj .0-1*)pd1=
      Flatten[Map[Take[#,1]&,pd0]];(*Gaus bez iznad ispod*)pd0=
      iteratedTake[pd0,2pd[[1]]];
    (*gaus za seckanje podeljen na komponente*)pd1=
      iteratedTake[pd1,2pd[[1]]];
    p=pd0[[k]];(*k-ta tj.komponenta koju secemo*)p1=
      Flatten[Map[Take[#,1]&,p]];(*p bez iznad ispod*)(*Print["Gaus",p1];*)
      dvojne=Union[Select[p1,Count[p1,#]==2&]];
    If[fDvojneQ[pd0],res={},Do[If[brzi,brzi=fMiror1[pd0,k,dvojne[[i]]];
          res=Append[res,brzi];
          If[brzi=={{1}},brzi=False,brzi=True]],{i,1,Length[dvojne]}];
      res];
    If[MemberQ[res,{{-1}}],res={-1}];
    res]
(*ako je res={-1} vise ne treba da radimo u suprotnom imamo pdata nekih \
clanova nove generacije*)

SamePdata[LL_List]:=
  Module[{l,res={LL[[1]]}},If[Length[LL]!=1,l=Rest[LL];
      While[l!={},l=Select[l,Not[SamePDQ[Last[res],#]]&];
        If[l!={},res=Append[res,First[l]];
          If[Length[l]==1,l={},l=Rest[l]]]]];
    res]
(*kad primrnimo fMirorKomp na sve komponente i to skupimo i jednu listu*)
\
(*to je lista pdata-pravimo gause-
    ako ima neki bez dvojnih to je kraj*)
(*radi od pdata*)

MirorAllKomp[Ul_]:=
  Module[{i=1,res={},p={2}},
    If[fComponentNo[Ul]==1,res=fMirorKomp[Ul,1],
      While[p!={-1}&&i<=fComponentNo[Ul],
        p=fMirorKomp[Ul,i];
        (*Print["Sredjena je ",i,"-ta komponenta: ",p];*)If[
          p!={-1}&&
            p!={{-1}},(*ako je lista pdata dodamo je u nove rez*)
            i++;
          If[p!={},res=Append[res,p]]
          (*ako je {} onda ne doprinosi rezultatu*),(*ako je {-1} ne treba \
vise da radimo*)p={-1};res={-1}]]];
    If[SameQ[res,{{}}]||SameQ[res,{}],res={-1}];
    (*Print[res];*)If[Not[SameQ[res,{-1}]],
      If[Length[res]==1,res=res[[1]]];
      res=SamePdata[res];
      res={res}];
    res]
(*vraca {-1} ako je razvezao a u suprotnom vraca listu pdata*)
fNSc[pd_List]:=Module[{NscNo=0,ind=True,gen={}},
    If[fDvojneQ[iteratedTake[fMakeGaussPD[pd],2pd[[1]]]],
      (*nema dvojnih*)(*Print["odma"];*)NscNo=0,
      (*ima dvojnih*)gen={pd};
      While[ind==True,NscNo++;
        gen=Flatten[Map[MirorAllKomp[#]&,gen],1];
        If[gen!={},gen=Flatten[gen,1]];
        ind=
          Not[MemberQ[gen,{-1}]||MemberQ[gen,-1]||SameQ[gen,{-1}]||
              SameQ[gen,{}]||SameQ[gen,{{}}]]]];
    
    NscNo]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #*)
NoSelfCrossNo[Ulaz_]:=Module[{pd,Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    If[SameQ[Head[Ul],String],pd=fCreatePData[Ul],pd=Ul];
    (*sad imamo pdata-odmah reduction*)pd=ReductionKnotLink[pd];
    Min[fNSc[pd],fNSc[GetMirrorImageKnot[pd]]]
    ]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #*)
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #*)
(* ## ## ## ## ## ## ## ## ## #  SAME KNOT LINK ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## #*)
(*uporedjuje dva KL data Konvejima *)

SameAltConKL[Con1_String,Con2_String]:=
  Module[{d1,d2,res},
    d1=MinDowAltKL[Con1];
      d2=MinDowAltKL[Con2];
    If[SameQ[Abs[d1],Abs[d2]],
          res=1 (*same*),
        res=0 (*not same*)];
    res]

(*uporedjuje dva KL data Konvejima *)

SameAltProjKL[Con1_,Con2_]:=
  Module[{d1,d2,res},
    d1=MinDowProjAltKL[Con1];
    d2=MinDowProjAltKL[Con2];
    If[SameQ[Abs[d1],Abs[d2]],
          res=1 (*same*),
        res=0 (*not same*)];
    res]
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
# *)

(*## ## ## ## ## ## ## ## ## GETBRAIDREPRESENT ## ## ## ## ## ## ## ## ## # \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)



fGaussExtSignsBraid[Ulaz_] := Module[{Ul = Ulaz, n, pp, pp1, 
pp2, ppr, ppi, i},
    (*Input : Pdata*)
    If[SameQ[Head[Ul], List], If[MemberQ[Map[OddQ, Abs[Ul[[2]]]], True],
        Ul], 
      If[SameQ[Head[Ul], String], Ul = fCreatePData[Ul], 
        Ul = fPDataFromDow[Ul]]];
    (*sad nam je Ul iliKonvej ili dowker*)
    n = Apply[Plus, Ul[[1]]];
    pp = Complement[Range[Apply[Plus, 2*Ul[[1]]]], Map[Abs, Ul[[2]]]];
    pp1 = Map[Sign, Ul[[2]]];
    pp2 = Table[If[EvenQ[Ul[[2, i]]], 1, -1], {i, n}];
    pp = Table[{pp[[i]], Abs[Ul[[2, i]]]}, {i, n}];
    ppr = Map[Reverse, pp];
    pp = Union[Table[{pp[[i]], {pp1[[i]], pp2[[i]]}}, {i, n}], 
        Table[{ppr[[i]], {pp1[[i]], pp2[[i]]}}, {i, n}]];
    pp1 = Map[First, pp];
    pp = Table[
        If[Map[First, Position[pp1, pp1[[i, 2]]]][[1]] == pp[[i, 1, 1]], 
          pp[[i]] , {Reverse[pp[[i]][[1]]], pp[[i]][[2]]}], {i, 2n}];
    pp = Table[{pp[[i, 1, 1]], pp[[i, 2]]}, {i, 2n}];
    ppr = Union[Map[First, pp]];
    pp = Table[{Flatten[Position[ppr, pp[[i, 1]]]][[1]], pp[[i, 2]]}, {i, 
          2n}];
    ppr = Map[Last, Union[Table[{pp[[i, 1]], pp[[i, 2, 1]]}, {i, 2n}]]];
    ppi = Table[-pp[[i, 2, 2]]*(-1)^i, {i, 2n}];
    pp = Map[First, pp];
    pp = Table[pp[[i]]*ppi[[i]], {i, 2n}];
    pp = {iteratedTake[pp, 2Ul[[1]]], ppr};
    pp
    ]
    
 fGetBraidRepresent[Ul_] := Module[{f = "br.txt", pp = Ul, pp1, pp2},
    pp = fGaussExtSignsBraid[pp];
    pp1 = ToString[pp[[1]]];
    pp2 = ToString[pp[[2]]];
    pp1 = 
      StringReplace[
        StringReplace[
          pp1, {"}, {" -> ";", "," -> "", "{" -> "", "}" -> ""} ], {";" -> 
            ", "}];
    pp2 = 
      StringReplace[
        pp2, {"1" -> "+", "-1" -> "-", "{" -> "", "}" -> "", "," -> " "}];
    pp = pp1 <> " / " <> pp2;
    OpenWrite[f];
    WriteString[f, pp];
    Write[f];
    Close[f];
    Run["fromdos" <> " br.txt"];
    Run["braid -v" <> " br.txt" <> " brout.txt"];
    pp = Import["brout.txt"];
    pp1 = StringPosition[pp, "="][[1, 1]];
    pp2 = StringLength[pp];
    pp = StringTake[pp, {pp1 + 1, pp2 - 2}]
    ]  
    
    

(*## ## ## ## ## ## ## ## ## BraidWReduce ## ## ## ## ## ## ## ## ## #*)
(*redukuje braid words pomocu Herves Iberovog programa
hr koji ucitavamo*)

BraidWReduced[PD_List]:=Module[{pp}, 
    pp=ReductionKnotLink[PD];
    (* pp=GetBraidRep[pp]; *)
    pp=fBraidW[pp];
    pp
        ]


(*## ## ## ## ## ## ## ## ## fClassicToCon ## ## ## ## ## ## ## ## ## #*)
fClassicToCon[Ulaz_String]:=Module[{brPres,pos,pom},
    pos=StringPosition[Ulaz,"_"][[1,1]];
    brPres=ToExpression[StringTake[Ulaz,pos-1]];
    pom=Flatten[Select[CCtoC,SameQ[#[[1]],brPres]&],1];
If[SameQ[pom,{}],{},
    pom=Rest[pom];
    Select[pom, SameQ[#[[1]],Ulaz]&][[1,2]]]]
    
 
fConToOther[Con_String]:=Module[{pp,pp1,res,res1,res2},
    pp=Flatten[Position[LLtoC,Con]];
    pp1=Flatten[Position[CCtoC,Con]];
    res1=If[Not[SameQ[pp,{}]],LLtoC[[pp[[1]],pp[[2]],pp[[3]]-1]],{}];
    res2=If[Not[SameQ[pp1,{}]],CCtoC[[pp[[1]],pp1[[2]],pp1[[3]]-1]],{}];
    res={res2,res1};
    res
    ]    
    
fAllCodes[Ul_]:=Module[{kk,kk1,res,dd},
If[SameQ[Ul,1],1,
If[ Ul[[1,1]]>15,Ul,kk=StringDrop[fDowMorwen[Ul],3];
kk1=StringJoin["K",First[uuuu[[First[Flatten[Position[uu,kk]]]]]]];
res=Flatten[If[StringLength[kk]<=10,{dd=fConNotation[kk1],fConToOther[dd]},{kk1}]];
res]]]


fDowThistToCon[Ulaz_String]:=Module[{pp1,brPres,pos,pom},
pp1=First[
      Map[First,Union[StringPosition[Ulaz,"a"],StringPosition[Ulaz,"n"]]]];
brPres=ToExpression[StringTake[Ulaz,{2,pp1-1}]];
pom=Flatten[Select[LLtoC,SameQ[#[[1]],brPres]&],1];
If[SameQ[pom,{}],{},
pom=Rest[pom];
Select[pom,SameQ[#[[1]],Ulaz]&][[1,2]]]]

fConNotation[Ul_String]:=Module[{res},
    res=If[SameQ[StringPosition[Ul,"_"],{}],fDowThistToCon[Ul],
        fClassicToCon[Ul]];
    res]

(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
# *)
(*## ## ## ## ## ## ## ## ## # EDMONDS ## ## ## ## ## ##*)



fGrCompl[Ul_List] := Module[{res,i},
    res = Union[Ul, Map[Reverse, Ul]];
    res = Table[Select[res, SameQ[#[[1]], i] &], {i, Max[Flatten[Ul]]}];
    res
    ]

fSelParts[LL_List] := Module[{res, i, j, L},
    res = 
      Flatten[Table[
          Join[{LL[[1, i]]}, {LL[[2, j]]}], {i, Length[LL[[1]]]}, {j, 
            Length[LL[[2]]]}], 1];
    Do[L = {res, LL[[i]]}; 
      res = Flatten[
          Table[Join[L[[1, i]], {L[[2, j]]}], {i, Length[L[[1]]]}, {j, 
              Length[L[[2]]]}], 1]; i = i + 1, {i, 3, Length[LL]}];
    res
    ]

fAllRot[Ul_List] := Module[{gg,i,j},
    gg = fGrCompl[Ul];
    gg = Table[
        Table[Join[{gg[[j, 1]]}, 
            DistinctPermutations[
                Table[gg[[j, i]], {i, 2, Length[gg[[j]]]}]][[i]]], {i, 
            Length[DistinctPermutations[
                Table[gg[[j, i]], {i, 2, Length[gg[[j]]]}]]]}], {j, 
          Length[gg]}];
    gg = fSelParts[gg];
    gg
    ]

fFindCyc[UlL_List, poc_List] := Module[{pp, pom, uu, uu1, vv1,vv, i},
    pom = {poc};
    uu = pom[[1]];
    While[Not[SameQ[uu1, pom[[1]]]], vv = Reverse[uu];
      vv1 = Position[UlL, Reverse[uu]];
      uu = 
        If[vv1[[1, 2]] + 1 <= Length[UlL[[1]]], 
          UlL[[vv1[[1, 1]], vv1[[1, 2]] + 1]], UlL[[vv1[[1, 1]], 1]]];
      uu1 = uu;
      pom = Join[pom, {uu1}]]; pp = Complement[Flatten[UlL, 1], pom];
    pom = Flatten[Map[#[[1]] &, Drop[pom, -1]]];
    If[pp != {}, 
      pp = Table[Select[pp, SameQ[#[[1]], i] &], {i, Max[Flatten[pp]]}]];
    {pom, pp}
    ]

fFindCycAll[UlL_List, poc_List] := 
  Module[{res = {}, ind = True, komplementi = {}, novpoc = poc, pom},
    While[Or[ind, komplementi != {}], ind = False;
      pom = fFindCyc[UlL, novpoc];
      If[poc == novpoc, komplementi = Flatten[pom[[2]], 1], 
        komplementi = Intersection[komplementi, Flatten[pom[[2]], 1]]];
      If[komplementi != {}, novpoc = komplementi[[1]];
        komplementi = Drop[komplementi, 1]];
      res = Append[res, pom[[1]]]];
     res]

fEdmonds[UlL_List] := Module[{gg, res = {}, i,j,rr1,rr2},
    gg = fAllRot[UlL];
    Do[res = Append[res, fFindCycAll[gg[[i]], First[gg[[i, 1]]]]]  , {i, 1, 
        Length[gg]}];
    res = Table[{gg[[i]], res[[i]]}, {i, Length[gg]}];
    rr1 = 
      Table[{Sort[Map[Length, res[[i, 2]]]], 
          Sort[Map[Length, Map[Split, Map[Sort, res[[i, 2]]]]]]}, {i, 
          Length[res]}];
     rr2 = Union[rr1];
     rr2 = 
      Flatten[Map[First, Table[Position[rr1, rr2[[i]]], {i, Length[rr2]}]]];
    res = Table[res[[rr2[[i]]]], {i, Length[rr2]}];
    res = 
      Table[{Table[Map[Last, res[[j, 1, i]]], {i, Length[res[[j, 1]]]}], 
          res[[j, 2]], 
          Max[Flatten[UlL]] + Length[res[[j, 2]]] - Length[UlL], 
          Divide[2 - Max[Flatten[UlL]] - Length[res[[j, 2]]] + Length[UlL], 
            2]}, {j, Length[res]}];
    res]


    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
#*)
(*## ## ## # STELAR ## ## ## ## ## ## ## ## ## ## ## #*)


    fStellarBasic[n_Integer] := Module[{del, del1, i, j, k},
    If[n < 6, del = {}, del = compo[n];
      del = Select[del, Not[MemberQ[#, 1]] &];
      del = Select[del, Length[#] >= 3 &];
      del1 = Sort[Table[{Sort[del[[i]]], del[[i]]}, {i, Length[del]}]];
      del = Union[Table[Sort[del[[i]]], {i, Length[del]}]];
      del = Table[Select[del1, #[[1]] == del[[i]] &], {i, Length[del]}];
      del = 
        Table[Table[del[[j, i, 2]], {i, Length[del[[j]]]}], {j, 
            Length[del]}];
      del = 
        Map[Last, 
          Flatten[Table[
              Union[Table[
                  Union[Table[
                      RotateLeft[del[[k, j]], i], {i, Length[del[[k, j]]]}], 
                    Map[Reverse, 
                      Table[RotateLeft[del[[k, j]], i], {i, 
                          Length[del[[k, j]]]}]]], {j, 
                    Length[del[[k]]]}]], {k, Length[del]}], 1]]];
    del = {del};
    del = 
      Flatten[Table[
          Map[StringReplace[ToString[#], 
          {"{" -> "", "}" -> "", " " -> ""}]&,
             del[[i]]], {i, Length[del]}]];
    del=Join[del,{Length[del]}];
    del]
(*Changed 9.03.2005*)


   
fMakePart[LL_List]:=Module[{res,pom,i=1},
    res=Map[compo[#]&,LL];
    While[i<=Length[res],
                pom=Select[res[[i]],First[#]!=1&];
      pom=Flatten[Map[DistinctPermutations,pom],1];
      pom=Select[pom, First[#]!=1&];
                
      If[pom=={},res=Drop[res,{i}],res=ReplacePart[res,pom,i]];
      i++
      ];
    Do[  pom=Map[StringReplace[ToString[#],", "->" "]&,res[[i]]];
             res=ReplacePart[res,pom,i]
      ,{i,1,Length[res]}];
    res
    ]
    
fCombine[res_List,NewL_List]:=Module[{i,res1,j,pom,combRes={}},
    res1=res;
    Do[
         Do[   
               pom=Join[res1[[j]],{NewL[[i]]}];
                combRes=Append[combRes,pom] ,
             {j,1,Length[res]}],
      {i,1,Length[NewL]}];
    combRes
    ]
    
fPerm[LL_List]:=Module[{i,perm={},res={}},
    perm=Map[{#}&,LL[[1]]];
    Do[perm=fCombine[perm,LL[[i]]],{i,2,Length[LL]}];
    perm
    ]

fStellar[n_Integer] := Module[{stel, i, j},
    stel = Drop[fStellarBasic[n],-1];
    stel = 
      Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], Number] &, 
        stel];
    stel = Map[fMakePart[#] &, stel];(*svaki strem veci od 1 partitionira :*)
    stel = Flatten[Map[fPerm[#] &, stel], 1];
    stel = 
      Table[Table[
          ToExpression[StringReplace[stel[[j, i]], {" " -> ","}]], {i, 
            Length[stel[[j]]]}], {j, Length[stel]}];
    stel = 
      Map[Last, 
        Union[Table[
            Union[Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}], 
              Map[Reverse, 
                Table[RotateLeft[stel[[j]], i], 
                {i, Length[stel[[j]]]}]]], {j,Length[stel]}]]];
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel=Join[stel,{Length[stel]}];      
    stel
    ]
    (*Changed 9.03.2005*)
    
    
    fStellarPlus[n_Integer]:=Module[{st, i, j},
    st=If[n<7,Print["n>6"],
        st=Table[Drop[fStellar[i],-1],{i,6,n-1}];
        st=Flatten[
            Table[Table[
                StringJoin[st[[i,j]],"+",ToString[n-5-i]],{j,
                  Length[st[[i]]]}],{i,Length[st]}]]];
    st=Join[st,{Length[st]}];
    st
    ]
    
    (*## ## ## # STELARNALT ## ## ## ## ## ## ## ## ## ## ## #*)
    
fStelString[Ul_String] := Module[{ss1, ss2, sp, ss = Ul, i, pp, pp1, pp0},
    pp0 = StringPosition[ss, "-"];
    pp = StringPosition[ss, " "];
    pp1 = StringPosition[ss, ","];
    If[Not[SameQ[Union[pp, pp1], {}]], pp = Union[First[Union[pp, \
pp1]]][[1]],
       pp = {}];
    ss = If[Not[SameQ[pp, {}]] && SameQ[pp0, {}] && pp > 2, 
        StringInsert[ss, " ", 2], ss];
    ss1 = Length[Union[Flatten[StringPosition[ss, ","]]]];
    ss2 = Length[Union[Flatten[StringPosition[ss, " "]]]];
    sp = Divide[Union[Flatten[StringPosition[ss, ","]]], 2];
    ss = Flatten[
        Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], 
              Number] &, {ss}]];
    ss = iteratedTake[ss, 
        Flatten[Prepend[{Length[ss] - sp[[-1]]}, 
            Flatten[Append[{sp[[1]]}, 
                Table[sp[[i]] - sp[[i - 1]], {i, 2, Length[sp]}]]]]]];
    ss]

    
    
    
     fCompl[Ul_List] := Module[{ss = Ul},
    ss = ReplaceAll[
        Table[If[ss[[i, -1]] == 1, 
            Drop[ReplacePart[ss[[i]], (ss[[i, -2]] + 1), -2], -1], 
            Append[ReplacePart[ss[[i]], ss[[i, -1]] - 1, -1], 1]], {i, 
            Length[ss]}], {1, 1} -> {2}];
    {Ul, ss}]
    
    fStelRot[Ul_List] := Module[{ppr, ppr1, ppr2, ppr3},
    ppr = Union[Table[RotateLeft[Ul[[1]], i], {i, Length[Ul[[1]]]}]];
    ppr1 = Map[Reverse, ppr];
    ppr2 = Union[Table[RotateLeft[Ul[[2]], i], {i, Length[Ul[[2]]]}]];
    ppr3 = Map[Reverse, ppr2];
    ppr = Union[ppr, ppr1, ppr2, ppr3];
    ppr
    ]
    
    fMakeStel[Ul_List] := Module[{stel = Ul, i},
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel]
    
    fStellarNalt[n_Integer] := Module[{stelpar1,stel, d, m, st1, stelpar, 
    stelnep, i, j, k},
    stel = Drop[fStellar[n],-1];
    d = IntegerPart[n/2];
    m = IntegerPart[d/2];
    st1 = 
      Table[StringJoin[stel[[i]], "+-1"], {i, 
          Length[stel]}];   (*svi sa jednim minusom*)
    stel = 
      Table[Select[stel, Length[StringPosition[#, ","]] == i &], {i, 3, 
          d - 1}];
    stelpar = Table[stel[[2i - 1]], {i, m - 1}];
    stelnep = Reverse[Complement[stel, stelpar]];
    stelnep = 
      Flatten[Table[
          Table[StringJoin[stelnep[[j, i]], "+-", ToString[k]], {k, 2, 
              j + 1}, {i, Length[stelnep[[j]]]}], {j, 
            Length[stelnep]}]];     (*dodati svi minusi neparnim *)
    stelpar1 = 
      Flatten[Table[
          Table[StringJoin[stelpar[[j, i]], "+-", ToString[k]], {k, 2, j}, \
{i,
               Length[stelpar[[j]]]}], {j, 
            Length[stelpar]}]];  
            (*dodati svi minusi parnim izuzev zabranjenog broja *)
    stelpar = Table[Map[fStelString, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = Table[Map[fCompl, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = 
      Map[Union, Table[Map[fStelRot, stelpar[[i]]], {i, Length[stelpar]}]];
    stelpar = Table[Map[First, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = Table[Map[Reverse, stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = Table[fMakeStel[stelpar[[i]]], {i, Length[stelpar]}];
    stelpar = 
      Flatten[Table[
          Table[StringJoin[stelpar[[j, i]], "+-", ToString[k]], {k, 2, 
              j + 1}, {i, Length[stelpar[[j]]]}], {j, Length[stelpar]}]];
    stel = Join[Reverse[Sort[st1]],Reverse[Sort[stelnep]], \
Reverse[Sort[stelpar]]];
   stel=Join[stel,{Length[stel]}];
    stel
    ]
    
    
(* ## ## ## ## ## ## ## # fStellarNinv  ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## *)


fStellarNinv[n_Integer] := Module[{stel, i, j}, stel = \
Drop[fStellarBasic[n],-1];
    stel = 
      Map[ReadList[StringToStream[StringReplace[#, "," -> " "]], Number] &, 
        stel];
    stel = Map[fMakePart[#] &, stel];
    stel = 
      Table[Table[
          Complement[
            Union[Table[
                If[SameQ[
                    fTangleType[
                      StringReplace[stel[[k, j, i]], {"{" -> "", "}" -> \
""}]],
                     0], {}, stel[[k, j, i]]], {i, 
                  Length[stel[[k, j]]]}]], {{}}], {j, Length[stel[[k]]]}], \
{k,
           Length[stel]}];
    (*svaki strem veci od 1 partitionira :*)
    stel = Flatten[Map[fPerm[#] &, stel], 1];
    stel = 
      Table[Table[
          ToExpression[StringReplace[stel[[j, i]], {" " -> ","}]], {i, 
            Length[stel[[j]]]}], {j, Length[stel]}];
    stel = 
      Map[Last, 
        Union[Table[
            Union[Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}], 
              Map[Reverse, 
                Table[RotateLeft[stel[[j]], i], {i, Length[stel[[j]]]}]]], \
{j,
               Length[stel]}]]];
    stel = 
      Table[Map[
          StringReplace[ToString[#], {"{" -> "", "}" -> "", " " -> ""}] &, 
          stel[[i]]], {i, Length[stel]}];
    stel = Table[StringReplace[stel[[i]], "," -> " "], {i, Length[stel]}];
    stel = 
      Table[StringReplace[
          ToString[stel[[i]]], {"{" -> "", "}" -> "", ", " -> ","}], {i, 
          Length[stel]}];
    stel]
    
    
    
fTangleType[Ul_] := Module[{rr, rr1},
    rr = Ul;
    rr = If[
        SameQ[Head[rr], String] && 
          SameQ[StringPosition[rr, " "], {}], {ToExpression[rr]}, rr];
    rr = If[SameQ[Head[rr], List] && SameQ[Length[rr], 1], 
        rr = Mod[rr, 2][[1]],
        If[SameQ[Head[rr], String],
          rr1 = Flatten[Map[Union, StringPosition[rr, " "]]];
          rr1 = Join[{StringTake[rr, {1, rr1[[1]]}]},
              
              Table[StringTake[rr, {rr1[[i]] + 1, rr1[[i + 1]] - 1}], {i, 1, 
                  Length[rr1] - 1}],
              {StringTake[rr, {rr1[[-1]], StringLength[rr]}]}];
          rr = ToExpression[rr1],
          rr = Ul];
        rr = ReplaceAll[Mod[rr, 2], {1 -> 3, 0 -> 2}];
        rr1 = Numerator[FromContinuedFraction[rr]];
        rr = If[EvenQ[rr1], 0, Join[rr, {0}]];
        rr1 = If[SameQ[rr, 0], 0, Numerator[FromContinuedFraction[rr]]];
        rr = If[SameQ[rr, 0], rr, If[EvenQ[rr1], 2, 1]]];
    rr]
    
    
fMakeType[n_Integer, k_Integer] := Module[{kk},
    kk = Select[compo[n], Not[SameQ[#[[1]], 1]] &];
    kk = Select[kk, SameQ[fTangleType[#], k] &];
    kk = Table[Map[ToString, kk[[i]]], {i, Length[kk]}];
    kk = Reverse[
        Sort[Table[
            StringDrop[
              StringJoin[
                Table[kk[[j, i]] <> " ", {i, Length[kk[[j]]]}]], -1], {j, 
              Length[kk]}]]];
    kk]
    
    
fStelTypeTest[Ul_String] := Module[{rr, rr1, rr2, dd, gg},
    rr = Ul;
    rr1 = Flatten[Map[Union, StringPosition[rr, ","]]];
    rr1 = Join[{StringTake[rr, {1, rr1[[1]]}]},
        Table[
          StringTake[rr, {rr1[[i]] + 1, rr1[[i + 1]] - 1}], {i, 1, 
            Length[rr1] - 1}],
        {StringTake[rr, {rr1[[-1]], StringLength[rr]}]}];
    rr1 = Table[StringReplace[rr1[[i]], "," -> ""], {i, Length[rr1]}];
    (*rr1 je  stel razbijen na R - tangles *)
    rr2 = Map[fTangleType, rr1];
    dd = Count[rr2, 1];
    rr2 = If[Not[MemberQ[rr2, 0]] && OddQ[dd], rr1, {}];
    rr2 = If[Length[Union[rr2]] < 3, {}, rr2];
    rr2 = If[SameQ[rr2, {}], rr2,
        rr2 = If[SameQ[rr2, {}], {}, {Length[rr2], rr2}];
        rr2 = 
          If[rr2[[1]] == 3, rr2[[2]], 
            If[SameQ[rr2[[1]], 4] && Length[Union[rr2]] < 4, 
              If[SameQ[rr2[[2, 1]], rr2[[2, 3]]] || 
                  SameQ[rr2[[2, 2]], rr2[[2, 4]]], {} , rr2]]]];
    rr2 = If[SameQ[rr2, {}], {}, rr];
    rr2
    ]
    
fNinvStellar[n_Integer] := Module[{st},
    st = fStellarNinv[n];
    st = Complement[Union[Map[fStelTypeTest, st]], {{}}];
    st=Complement[Map[fNinvStellarCorr,st],{{}}];
    st=Join[st,{Length[st]}];
    st]
    
    
fNinvStellarCorr[Ul_String]:=Module[{pp,pp1,pp2,pp3,res,i},
    pp=Union[{1},
        Union[Flatten[StringPosition[Ul,","]]],{StringLength[Ul]+1}];
    pp1=StringReplace[
        Table[StringTake[Ul,{pp[[i]],pp[[i+1]]-1}],{i,Length[pp]-1}],{","->
            ""}];
    pp2=Table[RotateRight[pp1,i],{i,Length[pp1]}];
    pp3=Table[Drop[pp2[[i]],1],{i,Length[pp2]}];
    pp2=Map[Reverse,pp3];
    res=If[
        SameQ[Plus@@Table[If[SameQ[pp3[[i]],pp2[[i]]],1,0],{i,Length[pp2]}],
          0],Ul,{}];
    res
    ]  



(*## ## ## ## ## ## ## ## ## # NMoveRat ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ##*)

ReduMove[n_Integer,Conway_List]:=
    Module[{Con,i=1,j,l,ind=False,Con1,Res={}},
      Switch[Head[Conway[[1]]],Integer,l=1,List,l=Length[Conway]];
      While[i<=l&&Not[ind],
        Switch[Head[Conway[[1]]],Integer,Con=Conway,List,Con=Conway[[i]]];
        Con=IsNotKnot[Con];
        If[
          Con=={}||Con=={0}||Con=={1}||Con=={-1}||
            Con=={2}||Con=={-2},ind=True;Res={},j=1;
          
          While[Not[MemberQ[Res,{0}]]&&j<=Length[Con],
            Res=Union[
                Append[Res,ReplacePart[Con,Con[[j]]-n*Sign[Con[[j]]],j]]];
            (*Print[Res,MemberQ[Res,{}]];*)j++]];
        i++];
      Res];

NMoveRat[n_Integer,Conway_String]:=
  Module[{Con,br=-1},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    Con=IsNotKnot[Con];
    If[Not[
        Con=={}||Con=={0}||Con=={1}||Con=={-1}||
          Con=={2}||Con=={-2}],
      While[Con!={},br++;Con=ReduMove[n,Con]],br=0];
    br]
    
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## #*)

(*## ## ## ## ## ## ## # fJablanPoly ## ## ## ## ## ## ## ## ## ##*)
(*calculates multivariable Jablan polynomial*)

fJablanPoly[Conway_String]:=
  Module[{lL0,res,lL,t,i,s,sl,j,p,MAT},
    If[Conway=="1",res=0,
      If[Conway=="2",res=1,lL0=fGenerators[Conway][[1]];
        lL=Flatten[lL0];
      t=ConstantArray[0,{1,Length[lL0]}][[1]];
        Do[
          If[i==1,t[[i]]=Length[lL0[[i]]]/3,
            t[[i]]=t[[i-1]]+Length[lL0[[i]]]/3],{i,1,Length[lL0]}];
        s=fGenerators[Conway][[2]];
        sl=Length[s];(*duzina*)
        MAT=ConstantArray[0,{sl,sl}];
        j=1;p="a";
        For[i=1,i<=sl,i++,If[i>t[[j]],j++;
            p=FromCharacterCode[ToCharacterCode["a"]+j-1]];
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+1]]]=1,-1,
            MAT[[i,Part[lL,3(i-1)+1]]]=-1];
          
          Switch[s[[i]],1,MAT[[i,Part[lL,3(i-1)+2]]]=-1,-1,
            MAT[[i,Part[lL,3(i-1)+2]]]=1];
          MAT[[i,Part[lL,3i]]]=p;];
        res=Det[MAT]
        ]];
    res]
(*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ##*)

(*## ## ## ## # Unlinking or unknotting number of fixed projection \
#### ## ## ## ## ## #*)

fGener[{LK_List,L_List}]:=
  Module[{pom,p,i,pl=L,l,rez={},pDa={},br=1},
    Do[If[MemberQ[Abs[pl],i]==False,
        pDa=Append[pDa,{i*Sign[pl[[br]]],pl[[br]]}];br++],{i,1,2Length[pl]}];
    Do[pom=pDa;
      l={LK};
      pom[[i]]*=-1;
      pom[[i]]=Reverse[pom[[i]]];
      pom=Sort[pom,Abs[#1[[1]]]<Abs[#2[[1]]]&];
      (*p je oblika originalnog linka*)
      
      p=Flatten[Map[Take[#,-1]&,pom]];
      l=Append[l,p];
      rez=Append[rez,l],{i,1,Length[pDa]}];
    rez;
    Union[rez]]

 (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## ##*)
 
    
   (*Ulaz je Konvej ili PData*)
(*Radi sve pojedinacne promene znaka na fix \
projekciji*)
fUnKLFixed[UL_,n_Integer]:=Module[{pd=UL,i,prvi,vv,res={}},
    If[SameQ[Head[UL],String],
    pd=fCreatePData[UL]];
    (*sad su pd pdata*)
    pd=ReductionKnotLink[pd];
    vv=fVarNewPGap[Length[pd[[2]]],n];
    pd=Union[Map[fCrChVarPD[pd,#]&,vv ]];
    prvi=Union[Map[Rest[#][[1]]&,pd]];
    pd=Sort[Map[Reverse[#]&,pd]];
    (*vv=Union[Map[Apply[Plus,#]&,vv]];Print[vv];*)
    Do[res=Append[res,
          Reverse[First[Select[pd,#[[1]]==prvi[[i]]&]]]
          ],      
      {i,1,Length[prvi]}];
    res=Reverse[Sort[res]];
    res=Select[res,SameQ[ToExpression[StringDrop[#[[2]],3]][[2]],{}] &];
    res=If[SameQ[res,{}],{},res[[1,1]]]]
    

    
    
    (*## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## #*)
    
    (* calculates UnKLNo od Konveja, pdata i dowkera *)
UnKLNo[Ulaz_]:=
  Module[{pD=Ulaz,str="",R,NoGen=1,un0},
    If[SameQ[Head[Ulaz],String],pD=fCreatePData[Ulaz]];
    If[SameQ[Head[Ulaz],List]&&SameQ[Select[Ulaz[[2]],OddQ[#]&],{}],
       pD=fPDataFromDow[Ulaz]];
      un0=MemberQ[ReductionKnotLink[pD],{}];
    If[un0==True,
          str="Unknott",
                  R=fGenerate[pD];
           
      While[MemberQ[R,{{},{}}]==False&&
          MemberQ[R,{{0},{}}]==False,  
        	R=Flatten[Map[fGenerate[#]&,R],1];
        	NoGen++];];
   NoGen](*16.8.2003 *)
   
   (*## ## ## ## ## ## ## ## ## ## ## fGap ## ## ## ## ## ## ## ## ##*)
   
   fGap[Conway_String]:=Module[{k,v,l},
    k=Last[UnKnotLink[Conway]];
    v=k;
    l=fUnKLFixed[Conway,k];
    While[SameQ[fUnKLFixed[Conway,k],{}],k++];
    {k-v,{v,k},Conway}] (*27.08.2004*)
    
  (* fGapKnotsc[Ul_List]:=Module[{k,dd,l,v},
    dd=fPDataFromDow[Ul];
    k=UnKLNo[dd];
    v=k;
    l=fUnKLFixed[dd,k];
    While[SameQ[fUnKLFixed[dd,k],{}],k++];
    {k-v,{v,k},Ul}] *)
    
    (* SLAVIK 11.02.2007 *)
   
       
    (*## ## ## ## ## ## ## ## fUnRFixProj ## ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## #*)
    
    (* tests that a fixed projection of rational 
    KL given by Conway can be unlinked by k crossing changes. Makes
    crossing changes directly in Conway symbol *)
    
    UnRFixPr[Conway_String,k_Integer]:=Module[{par,Con,n,r,p,i,res},
    Con=ReadList[
        StringToStream[
          StringReplace[
            Conway,{","->" ","+"->" ","("->" ",
              ")"->" "}]],Number];
    n=Apply[Plus,Con];
    r=Length[Con];
    If[SameQ[r,1]&&EvenQ[n],res=If[SameQ[k,n/2],1,0],
    p=parr[k,r];
    p=Select[p,Max[#]<=Max[Con]&];
    p=Flatten[Union[Table[DistinctPermutations[p[[i]]],{i,Length[p]}]],1];
    p=Table[Con-p[[i]],{i,Length[p]}];
    p=Select[p,Min[#]>=0&];
    p=Table[Con-p[[i]],{i,Length[p]}];
    p=Table[Con-2p[[i]],{i,Length[p]}];
    p=Map[IsNotKnot,p];
    res=If[MemberQ[p,{1}]||MemberQ[p,{-1}]||MemberQ[p,{0}],1,0]];
    res]
    (*24.6.2004*)
    
    
    fUnRFixProj[Con_]:=Module[{k=1},
    While[SameQ[UnRFixPr[Con,k],0],k++];
    k
    ]
    
    (*## ## ## ## ## ## ## ## ## ## ## fGapRat ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## #*)
    
    (* calculates unlinking gap, unlinking no. and
    fixed projection unlinking no. for a rational KL *)
    
   fGapRat[Ul_]:=Module[{uk,up},
    uk=UnR[Ul][[-1]];
    up=fUnRFixProj[Ul];
    uk={up-uk,{uk,up}};
    uk
    ]
    
    (*## ## ## ## ## ## ## ## ## ## ## AUTOMATS ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## #*)
    (*## ## ## ## ## ## ## ## ## # fAutoSigInp ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ##*)
    
    (*Gr_List={outgiong edges of graph, signlist,list of inputs}, NAND=1, \
AND=-1*)
    
fAutoSigInp[Gr_List]:=Module[{gr,s,e,n,ul,iz,i,j,k,v,res,st,inp,in},
    gr=Gr[[1]];
    s=Gr[[2]];
    inp=Gr[[3]];
    gr=Flatten[Table[Table[{i,gr[[i,j]]},{j,Length[gr[[i]]]}],
    {i,Length[gr]}],1];
    gr=Union[gr];
    e=Length[gr];
    v=fVarP[e];
    n=Length[Union[Flatten[gr]]];
    s=If[SameQ[s,{}],Table[1,{i,n}],s];
    in=Table[1,{i,n}];
    inp=If[SameQ[inp,{}],inp,Table[{inp[[i]]},{i,Length[inp]}]];
    inp=If[SameQ[inp,{}],in,ReplacePart[in,0,inp]];
    e=Table[{gr[[i]],i},{i,Length[gr]}];
    ul=Table[Select[e,#[[1,2]]==i&],{i,n}];
    ul=Table[Table[ul[[j,i,2]],{i,Length[ul[[j]]]}],{j,Length[ul]}];
    iz=Table[Select[e,#[[1,1]]==i&],{i,n}];
    iz=Table[Table[iz[[j,i,2]],{i,Length[iz[[j]]]}],{j,Length[iz]}];
    st=iz;
    ul=Flatten[
        Table[Table[
            Table[v[[i,ul[[j,k]]]],{k,Length[ul[[j]]]}],{j,Length[ul]}],{i,
            Length[v]}],1];
    iz=Flatten[
        Table[Table[
            Table[v[[i,iz[[j,k]]]],{k,Length[iz[[j]]]}],{j,Length[iz]}],{i,
            Length[v]}],1];
    e=Map[Max,Partition[Map[Length,Map[Union,iz]],n]];
    ul=Partition[Map[Min,ul],n];
    ul=Table[inp*ul[[i]],{i,Length[ul]}];
    iz=Partition[Map[Min,iz],n];
    iz=Table[(s+1)/2-s*iz[[i]],{i,Length[iz]}];
    res=Union[
        Table[If[SameQ[ul[[i]],iz[[i]]]&&SameQ[e[[i]],1],v[[i]],{}],{i,
            Length[v]}]];
    res=If[SameQ[res[[1]],{}],Drop[res,1],res];
    e=Map[Length,st];
    e=Flatten[
        Map[Union,
          Flatten[Union[Table[iteratedTake[res[[i]],e],{i,Length[res]}]],
            1]]];
    e=Partition[e,Length[st]];
    e={res, e};
    e]
    
   
    (*## ## ## ## ## ## ## ## ## # fAutoKL ## ## ## ## ## ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ##*)
    
    (* fAutoKL[Ul_List]: Ul={Conway,signs,inputs}; 
    If[signs={},signs=signs (Con) *)
    
    
fAutoKL[Con_,Sl_,In_]:=Module[{gg,gg1,i,j,s,inp,ss},
  If[SameQ[Head[Sl],String],inp=ToExpression[Sl],inp=Sl];
   If[SameQ[Head[In],String],ss=ToExpression[In],ss=In];
    gg=fGaussExtSigns[Con];
    Print["Gauss code:", gg];
    s=If[SameQ[ss,{}],Map[Sign,Reverse[Union[Flatten[gg]]]],ss];
    Print["Signs:", s];
    gg1=Union[Flatten[gg]];
    gg1=Sort[Table[{Abs[gg1[[i]]],Sign[gg1[[i]]]},{i,Length[gg1]}]];
    gg1=Table[gg1[[i,2]],{i,Length[gg1]}];
    gg=If[SameQ[gg,Flatten[gg]],
        Sort[Join[Abs[{{gg[[-1]],gg[[1]]}}],
            Table[Abs[{gg[[i]],gg[[i+1]]}],{i,Length[gg]-1}]]],
        Sort[Flatten[
            Table[Join[Abs[{{gg[[j,-1]],gg[[j,1]]}}],
                Table[Abs[{gg[[j,i]],gg[[j,i+1]]}],
                {i,Length[gg[[j]]]-1}]],{j,Length[gg]}],1]]];
                Print["Graph: ", gg];
    gg=Table[Select[gg,#[[1]]==i&],{i,Length[Union[Flatten[gg]]]}];
    gg=Table[Table[gg[[i,j,2]],{j,Length[gg[[i]]]}],{i,Length[gg]}];
    gg=fAutoSigInp[{gg,s,inp}];
    Print["Egde coloring: ", gg[[2]]];
    Print["Stable states: ", gg[[1]]]; 
    ]
    
    
(*## ## ## ## ## ## ## ## # LiangPoly ## ## ## ## ## ## ## ## ## ## ## #*)
    

LiangPoly[Conway_String]:=Module[{res,MM,gg,ss,n,i,j},
    gg=fGraphInc[Conway];
    ss=gg[[2]];gg=gg[[1]];
    n=Length[ss];
    MM=DiagonalMatrix[Map[Power[h,#]&,ss]];(*Print[n,MM];*)
    Do[Do[
        p=Count[gg,{i,j}];
        MM=ReplacePart[MM,p,{i,j}];
        MM=ReplacePart[MM,p,{j,i}];
        ,{j,i+1,n}],{i,1,n-1}];
    (*Print[MM];*)res=Det[MM];
    Print["Amphicheiral ",
      SameQ[res,Det[ReplaceAll[MM,{h->1/h,1/h->h}]]]];
    res
    ]
    
    
    LiangPoly1[Conway_String]:=Module[{res,MM,gg,ss,n,i,j},
    gg=fGraphInc[Conway];
    ss=gg[[2]];gg=gg[[1]];
    n=Length[ss];
    MM=DiagonalMatrix[Map[Power[h,#]&,ss]];(*Print[n,MM];*)
    Do[Do[
        p=Count[gg,{i,j}];
        MM=ReplacePart[MM,p,{i,j}];
        MM=ReplacePart[MM,p,{j,i}];
        ,{j,i+1,n}],{i,1,n-1}];
    (*Print[MM];*)res=Det[MM];
    (* Print["Amphicheiral ",
      SameQ[res,Det[ReplaceAll[MM,{h->1/h,1/h->h}]]]]; *)
    res
    ]
    
(*## ## ## ## ## ## ## ## ## fBoolean ## ## ## ## ## ## ## ## ## ## ## ##*)

(*pomocna za fBrisi stepen: u pol sve stepene xx mkonvertuje u xx *)

fImaX[Pol_,xx_]:=Module[{pom=Pol,ind=True},
    If[MemberQ[Variables[pom],xx],
      While[ind,
        ind=MemberQ[Variables[pom],xx]&&Not[MemberQ[Variables[pom/xx],xx]];
        If[ind,ind=False,pom=pom/xx;ind=True]]];
    pom]
    
(*stepene svih varijabli konvertuje u same varijable*)

fBrisiStepen[Pol_]:=Module[{pom=Pol,i,var},
    var=Variables[Pol];
    Do[pom=Map[ fImaX[#,var[[i]]]&,pom]
      ,{i,1,Length[var]}];pom
    ]
    
(*racuna polinome za logicke funkcije*)

fBoolean[Ul_]:=Module[{p,LL},
  If[SameQ[Head[Ul],String],LL=ToExpression[Ul],LL=Ul];
    Unprotect[And];
    Unprotect[Or];
    Unprotect[Not];
    Unprotect[Nand];
    Unprotect[Nor];
    Unprotect[Implies];
    Unprotect[Xor];
    Unprotect[Equal];
    a_&&b_:=a b;
    !a_:=1-a;
    a_||b_:=a+b-a b;
    Not[a_]:=1-a;
    Nor[a_,b_]:=1-a-b+a b;
    Nand[a_,b_]:=1-a b;
    Implies[a_,b_]:=1-a+a b;
    Xor[a_,b_]:=a+b-2a b;
    Equal[a_,b_]:=1-a-b+2a b;
    p=Expand[LL];
    Clear[And];
    Clear[Or];
    Clear[Not];
    Clear[Nand];
    Clear[Nor];
    Clear[Implies];
    Clear[Xor];
    Clear[Equal];
    p=fBrisiStepen[p];
    p
    ]
    
    

(*## ## ## ## ## ## ## ## ## fDiffSeq ## ## ## ## ## ## ## ## ## ## ## ##*)

  
(* fPerBasic produces material for periodic sequences 
without dual material
or constant parts *)
fPerBasic[n_Integer,p_Integer]:=Module[{v,u,vr,vp,vd, i, j, k},
    v=Sort[fVarP[n]];
    v=KSubsets[v,p];
    u=Table[Table[Table[v[[k,i,j]],{i,p}],{j,n}],{k,Length[v]}];
    v=Union[
        Table[If[Length[Union[Map[Length,Map[Union,u[[i]]]]]]==1,
            v[[i]],{}],{i,Length[v]}]];
    v=Map[Sort,If[v[[1]]=={},Drop[v,1],v]];
    vd=Flatten[
        Table[{i,
            Position[v,Sort[ReplaceAll[v[[i]],{0->1,1->0}]]]},{i,
            Length[v]}]];
    vd=Flatten[Union[Map[Sort,Partition[vd,2]]]];
    vd=Table[vd[[2i-1]],{i,Length[vd]/2}];
    v=Map[Sort,Table[v[[vd[[i]]]],{i,Length[vd]}]];
    vp=Union[
        Map[Sort,Table[Table[Table[v[[k,i,j]],{i,p}],{j,n}],{k,Length[v]}]]];
    v=Table[
        Table[Table[vp[[k,i,j]],{i,Length[vp[[k]]]}],{j,p}],{k,Length[vp]}]
    ]
(* fPerm produces permutations different with regard to automorphisms of \
dihedral group - rotations and reverse *)
fPerm[n_Integer]:=Module[{p,q,l, i, j},
    p=Permutations[Range[n]];
    p=Map[First,
        Union[Table[Sort[Table[RotateLeft[p[[i]],j],{j,n}]],{i,Length[p]}]]]
    ]
fPerSeq[n_Integer,p_Integer]:=Module[{v,s, i, j, k},
    If[2^n>=p>2^(n-2),
      v=fPerBasic[n,p];
      s=fPerm[p];
      s=Flatten[
          Table[Table[Table[v[[k,s[[j,i]]]],{i,p}],{j,Length[s]}],{k,
              Length[v]}],1];
      s=Union[
          Table[Sort[Table[Table[s[[k,i,j]],{i,p}],{j,n}]],{k,Length[s]}]];
      s=Union[Table[Table[Table[s[[k,i,j]],{i,n}],{j,p}],{k,Length[s]}]],
      Print["2^(n-2)<p<=2^n"]]
    ]
fKauffAlg[n_Integer,p_Integer]:=Module[{l,t, i, j, k},
    l=fPerSeq[n,p];
    t=Table[
        Join[{Table[
              If[l[[k,1,j]]==0,l[[k,-1]],{}],{j,Length[l[[k,1]]]}]},
          
          Table[Table[
              If[l[[k,i,j]]==0,l[[k,i-1]],{}],{j,Length[l[[k,i]]]}],{i,
              2,Length[l[[k]]]}]],{k,Length[l]}];
    t=Table[Map[Union,Table[Table[t[[k,i,j]],{i,p}],{j,n}]],{k,Length[l]}];
    t=Table[
        Table[If[t[[j,i,1]]=={},Drop[t[[j,i]],1],t[[i,1]]],{i,
            Length[t[[j]]]}],{j,Length[t]}]
    ]
(*pomocna za fBrisi stepen: u pol sve stepene xx mkonvertuje u xx *)

fImaX[Pol_,xx_]:=Module[{pom=Pol,ind=True},
    If[MemberQ[Variables[pom],xx],
      While[ind,
        ind=MemberQ[Variables[pom],xx]&&Not[MemberQ[Variables[pom/xx],xx]];
        If[ind,ind=False,pom=pom/xx;ind=True]]];
    pom]
(*stepene svih varijabli konvertuje u same varijable*)

fBrisiStepen[Pol_]:=Module[{pom=Pol,i,var},
    var=Variables[Pol];
    Do[pom=Map[ fImaX[#,var[[i]]]&,pom]
      ,{i,1,Length[var]}];pom
    ]

fOp[L_List]:=Module[{l,k=1,i,j,n,pom,n1=Length[L]},
    l=L;
    n=Length[l[[1,1]]];
    (*n je broj varijabli koje koristimo*)
    (*Print[l,n];*)
    Do[
      pom=Symbol[FromCharacterCode[96+k]];
      Do[Do[
            
          If[l[[i,j,k]]==1,l=ReplacePart[l,pom,{i,j,k}],
            l=ReplacePart[l,1-pom,{i,j,k}]]
          ,{j,1,Length[L[[i]]]}],{i,1,n1}]
      ,{k,1,n}];
    (*Print["LLL ",l];*)
    
    l=Table[Expand[
          Times@@Table[1-Expand[Times@@l[[j,i]]],{i,Length[l[[j]]]}]],{j,
          Length[l]}];
    (*Print["Polinomi",l];*)
    l=Map[fBrisiStepen,l]]
(*pomocna za fDiffSeq*)
(*pomocna za fPraviZamene i fBalanced*)

fPZ[PP_List,CP_List,OP_Integer]:=Module[{res={},i},
    Do[  If[SameQ[OP,0],res=Append[res,{PP[[i]]->CP[[i]]}],
                                                    
        res=Append[res,{PP[[i]]==CP[[i]]}]],{i,1,Length[CP]}];
    Flatten[res,1]]
(*pomocna za fDiffSeq*)
(*za n simbole pravi listu svih mogucih zamena na \
osnovu Distinct permutations*)
fPraviZamene[n_Integer]:=Module[{l,p},
    l=DistinctPermutations[
        Map[Symbol[FromCharacterCode[#]]&,Range[97,96+n]]];
    l=Map[fPZ[First[l],#,0]&,l]
    ]
(*pomocna za fDiffSeq*)
(*svakom elementu liste na odgovarajuce mesto \
dopisuje uredjen par sa odgovarajucom promenljivom- simbolom*)

fDodajSlovo[n_Integer, LL_List]:=Module[{l,res=LL},
    l=Map[Symbol[FromCharacterCode[#]]&,Range[97,96+n]];
    Do[  res=Map[ReplacePart[#,{l[[i]],#[[i]]},i]&,res],{i,1,n}];res]

fDiffSeq[n_Integer,p_Integer]:=
  Module[{pom,l,lll={},ll,pp,qq,zam,i},l=fPerSeq[n,p];
    pp=fKauffAlg[n,p];
    qq=pp;
    pp=Map[fOp,pp];
    pp=fDodajSlovo[n,pp];
    (*svakom elementu dodaje odgovarajucu promenljivu*)
    
    zam=fPraviZamene[n];
    (*zamene na osnovu svih permutacija od n-simbola*)
    
    Do[pom=Sort[Map[Sort[ReplaceAll[pp[[i]],#]]&,zam]];
      lll=Append[lll,pom],{i,1,Length[pp]}];
    ll=Union[lll];
    pp=Sort[Table[Position[lll,ll[[i]]],{i,Length[ll]}]];
    pp=Table[pp[[i,1,1]],{i,Length[pp]}];
    ll=Table[lll[[pp[[i]]]],{i,Length[pp]}];
    l=Table[l[[pp[[i]]]],{i,Length[pp]}];
    qq=Table[qq[[pp[[i]]]],{i,Length[pp]}];
    ll=If[n>2, 
        Table[{ll[[i,1]][[1,2]],ll[[i,1]][[2,2]],ll[[i,1]][[3,2]]},{i,
            Length[ll]}], 
        Table[{ll[[i,1]][[1,2]],ll[[i,1]][[2,2]]},{i,Length[ll]}]];
    ll=Table[{l[[i]],qq[[i]],ll[[i]]},{i,Length[l]}];
    Print["There are ",Length[ll]," different periodic sequences for ",n,
      " & ",p];
    ll]
    
fBal[Ulaz_]:=Module[{Diff, res},
If[SameQ[Head[Ulaz],String],Diff=ToExpression[Ulaz],Diff=Ulaz];
res=fBalanced[{Diff},1];
res]
    
fBalanced[Ulaz_,redBr_Integer]:=Module[{n,poly,l1,l2,Diff,fPZ},
 If[SameQ[Head[Ulaz],String],Diff=ToExpression[Ulaz],Diff=Ulaz];
    n=Length[Diff[[1,1,1]]];
    If[redBr>Length[Diff],
         Print["There are only ", Length[Diff],
        " different periodic sequences. Please pick lesser number"]
      ,
      poly=Diff[[redBr,3]];
      l2=Map[Symbol[FromCharacterCode[#]]&,Range[97,96+n]];
      l1=fPZ[poly,l2,1]
      ];
    Print["Periodic sequence: ", Diff[[redBr,1]]];
    Print["Kauffman code: ", Diff[[redBr,2]]];
    Print["Polynomials: ",poly];
    Print["Balanced states: ",Solve[l1,l2]];    
    ]
    
    
(*## ## ## ## ## ## ## ## ## fKnotFind ## ## ## ## ## ## ## ## ## ## ## ##*)


fForknotFind[Ul_]:=Module[{l,DL,str,f="ff.txt"},
    If[SameQ[Head[Ul],String],
      If[SameQ[StringPosition[Ul,"#"],{}],DL=Dowker[Ul][[4]]; 
        l=fGenSign[Ul]*fGenSign[StringReplace[Ul,{"-"->""}]];
        DL*=l,DL= fDToDDirect[Ul][[2]]], DL=Ul[[2]]];
    OpenWrite[f];
    str=ToString[Length[DL]]<>" 1 "<>" "<>
        StringReplace[ToString[DL],{"{"->"","}"->"",","->""}];
    WriteString[f,str];
    Write[f];
    Close[f]]
    

fKnotFind[Ulaz_]:=Module[{uu0,DL,Ul,l,str,link,res,rr,res1,tt,pp,pp1},
Off[LinkObject::"linkd"];
Off[LinkObject::"linkn"];
uu0=If[SameQ[fComponentNo[Ulaz],1],Ul=Ulaz,0];
rr=If[SameQ[uu0,0],0,
 If[SameQ[Head[Ul],String],
      If[SameQ[StringPosition[Ul,"#"],{}],DL=Dowker[Ul][[4]]; 
        l=fGenSign[Ul]*fGenSign[StringReplace[Ul,{"-"->""}]];
        DL*=l,DL= fDToDDirect[Ul][[2]]], DL=Ul[[2]]];
 str=ToString[Length[DL]]<>" 1 "<>" "<>
        StringReplace[ToString[DL],{"{"->"","}"->"",","->""}];
str=ToExpression[StringJoin["{",StringReplace[StringReplace[str,{"   "->" ","  "->" "," "->","}]," "->","],"}"]];
link=Install["./KnotFindML",LinkMode->Launch];
res=fKnotFindML[str];
pp=Uninstall[link];
pp1=SameQ[res,$Failed];
If[SameQ[pp1,True],fKnotFindOld[Ulaz],
rr=If[SameQ[res,"{1}"],1,res1=ToExpression[If[SameQ[Length[StringPosition[res,"{"]],1],
StringJoin["{",StringReplace[res,{"{ "->"{","     "->",","   "->",","  "->","," "->","}],"}"],
StringJoin["{",StringReplace[res,{"{ "->"{","   "->",","  "->","," "->","}],"}"]]];
res1=Table[Select[res1[[i]],NumberQ[#] &],{i,Length[res1]}];
If[SameQ[Length[res1],1],{{{res1[[1,1]]},Rest[res1[[1]]]}},tt=Table[Drop[res1[[i]],2],{i,Length[res1]}];
tt=Table[{{Length[tt[[i]]]},tt[[i]]},{i,Length[tt]}]],0 ];
rr=If[SameQ[rr,1],rr,If[Length[rr]>1,rr,rr[[1]]]]]];
rr]


fKnotFindOld[Ul_]:=Module[{cc,kk,res,res1,tt,rr,i},
cc=If[SameQ[Head[Ul],String],
        If[SameQ[StringPosition[Ul,"#"],{}],Length[Dow[Ul][[1]]],
          Length[fDToDDirect[Ul][[1]]]],Length[Ul[[1]]]];
rr=If[SameQ[cc,1],
  kk=fForknotFind[Ul];
     (* Run["fromdos.exe","ff.txt"]; *)
      Run["knotfind.exe","ff.txt","rez.txt"];
res=Import["rez.txt"];
If[SameQ[StringTake[res,1],"u"],1,res1=ToExpression[If[SameQ[Length[StringPosition[res,"{"]],1],
StringJoin["{",StringReplace[res,{"{ "->"{","     "->",","   "->",","  "->","," "->","}],"}"],
StringJoin["{",StringReplace[res,{"{ "->"{","   "->",","  "->","," "->","}],"}"]]];
res1=Table[Select[res1[[i]],NumberQ[#] &],{i,Length[res1]}];
If[SameQ[Length[res1],1],{{{res1[[1,1]]},Rest[res1[[1]]]}},tt=Table[Drop[res1[[i]],2],{i,Length[res1]}];
tt=Table[{{Length[tt[[i]]]},tt[[i]]},{i,Length[tt]}]],0 ]];
rr=If[SameQ[rr,1],rr,If[Length[rr]>1,rr,rr[[1]]]];
rr](* korigovan 26.10.2010 *)


fKnotFindDirProd[Ul_]:=Module[{kk,kk0,pp1,dd},
kk=If[SameQ[Head[Ul],String],fKnotFind[Ul],
If[SameQ[Head[Ul[[1,1]]],List],Ul,fKnotFind[Ul]]];
If[SameQ[kk,1],1,
kk0=kk;
kk=Select[kk,Depth[#]>2 &];
pp1=If[SameQ[kk,{}],kk,kk[[1]]];
dd=Depth[kk];
pp1=If[SameQ[dd,2],Join[pp1,kk],pp1];
If[SameQ[dd,2],kk0,
While[dd>2,kk=Map[fKnotFind,kk];
kk0=kk[[-1]];
kk=Select[kk,Depth[#]>3 &];
kk=Flatten[kk,1];
pp1=Join[pp1,If[SameQ[kk,{}],kk,kk[[1]]]];
dd=Depth[kk];
pp1=If[SameQ[dd,2],Join[pp1,kk],pp1];
pp1=If[SameQ[dd,2]&&SameQ[kk,{}],Join[pp1,kk0],pp1]];
Sort[Map[fKnotFind,Partition[pp1,2]]]]]]
      

RedKauffmanPolynomial[Ulaz_] := Module[{kk, ff,Ul},
    If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    kk = KauffmanPolynomial[Ul];
    ff=If[SameQ[Abs[kk],1],Abs[kk],
    If[SameQ[Position[kk[[1]], a], {}], kk,
      kk = ReplaceAll[kk, xa -> x*a]];
    ff = FactorList[kk][[2]];
    If[SameQ[Position[ff, x], {}],
      ff = Expand[Divide[kk, Power[ff[[1]], ff[[2]]]]], ff = kk];
    If[ff == 1, ff = kk];
    ff=If[SameQ[ff,1],ff,If[SameQ[ff,-1],1,
    If[SameQ[Sign[ff[[1]]],1],ff,-ff]]]];
    ff]


RedKauffman[Ulaz_] := Module[{kk, ff,Ul},
    If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    kk = KauffmanPolynomial[Ul];
    ff=If[SameQ[Abs[kk],1],Abs[kk],
    If[SameQ[Position[kk[[1]], a], {}], kk,
      kk = ReplaceAll[kk, xa -> x*a]];
    ff = FactorList[kk][[2]];
    If[SameQ[Position[ff, x], {}],
      ff = Expand[Divide[kk, Power[ff[[1]], ff[[2]]]]], ff = kk];
    If[ff == 1, ff = kk];
    ff=If[SameQ[ff,1],ff,If[SameQ[ff,-1],1,
    If[SameQ[Sign[ff[[1]]],1],ff,-ff]]]];
    ff]
     

RedAlex[Ulaz_]:=Module[{al,dd,pp,Ul},
dd=If[SameQ[RedCon[Ulaz],0],1,
 If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    al=SkeinPolynomial[-1,Ul];
   (*  al = ReplaceAll[al, x -> Sqrt[x]]; *)
    dd=Table[al[[i]],{i,Length[al]}];
    dd=If[SameQ[Abs[al],1],Abs[al],
    pp=PolynomialGCD@@dd;
    dd=Expand[Divide[al,pp]];
    dd=If[SameQ[dd,1],dd,If[SameQ[dd,-1],1,
    If[SameQ[Sign[dd[[1]]],1],dd,-dd]]]]];
    dd
    ]
    

RedJones[Ulaz_]:=Module[{al,dd,pp,Ul},
 If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    al=SkeinPolynomial[0,Ul];
   (*  al = ReplaceAll[al, x -> Sqrt[x]]; *)
    dd=Table[al[[i]],{i,Length[al]}];
     dd=If[SameQ[Abs[al],1],Abs[al],
    pp=PolynomialGCD@@dd;
    dd=Expand[Divide[al,pp]];
    dd=If[SameQ[dd,1],dd,If[SameQ[dd,-1],1,
    If[SameQ[Sign[dd[[1]]],1],dd,-dd]]]];
    dd
    ]
    
(* Slavik, March 31, 2009 *)
    
RedCon[Ulaz_]:=Module[{al,dd,pp,Ul},
 If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    al=SkeinPolynomial[-2,Ul];
   dd= If[SameQ[al,0],0,
    al = ReplaceAll[al, x -> x^2]; 
     dd=If[SameQ[Abs[al],1],Abs[al],
    dd=Table[al[[i]],{i,Length[al]}];
    pp=PolynomialGCD@@dd;
    dd=Expand[Divide[al,pp]];
    dd=If[SameQ[dd,1],dd,If[SameQ[dd,-1],1,
    If[SameQ[Sign[dd[[1]]],1],dd,-dd]]]]];
    dd
    ]
    
(* Slavik, March 31, 2009 *)    
    
    
RedHom[Ulaz_]:=Module[{al,dd,pp,Ul},
 If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    al=SkeinPolynomial[1,Ul];
     dd=If[SameQ[Abs[al],1],Abs[al],
    dd=Table[al[[i]],{i,Length[al]}];
    pp=PolynomialGCD@@dd;
    dd=Expand[Divide[al,pp]];
   dd=If[SameQ[dd,1],dd,If[SameQ[dd,-1],1,
   If[SameQ[Sign[dd[[1]]],1],dd,-dd]]]]; 
    dd
    ]
    
 (* Slavik, March 31, 2009 *)

fKauff[Ulaz_] := Module[{kk,Ul},
    If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    kk = KauffmanPolynomial[Ul];
  kk]

fJon[Ulaz_]:=Module[{al,Ul},
 If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    al=SkeinPolynomial[0,Ul];
al]

fHom[Ulaz_]:=Module[{al,Ul},
 If[SameQ[Head[Ulaz],String],Ul=fCreatePData[Ulaz],Ul=Ulaz];
    al=SkeinPolynomial[1,Ul];
    al]

fKh[Ulaz_]:=Module[{al,Ul},
 If[SameQ[Head[Ulaz],String],Ul=fConwayToPD[Ulaz],Ul=fPdataToPD[Ulaz]];
    al=KnotTheory`Kh[mm,Program->"JavaKh"][q,t];
    al]


fCheckGap[mm_]:=Module[{mm1,mm2,i},
mm1=Table[mm[[i]],{i,Length[mm]}];
mm2=CoefficientList[mm,x];
mm1=If[MemberQ[mm2,0],1,0];
mm1]

fCheckSign[mm_]:=Module[{mm1,mm2,pp1,i},
mm1=Table[mm[[i]],{i,Length[mm]}];
mm2=Map[Sign,CoefficientList[mm,x]];
pp1=Table[(-1)^i,{i,Length[mm2]}];
mm1=If[SameQ[mm2,pp1]||SameQ[mm2,-pp1],0,1];
mm1
] 

fSpecial[Ul_]:=Module[{pp1},
pp1=If[SameQ[fCheckGap[Ul],1],1,fCheckSign[Ul]];
pp1
] 
     
     (* ## ## ## ## ## ## # NONALGEBRAIC TANGLES ## ## ## ## ## ## *)
     
PDataSumProdNonA[L_List, tangle_String] := 
  Module[{BrPreseka, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, 
      x21, x22, res, str, f = "NonABase.txt"},
    BrPreseka = Divide[Length[L[[2]]], 2];
    Switch[BrPreseka, 10, x10 = L, 11, x11 = L, 12, 
    x12 = L, 13, x13 = L, 14, 
      x14 = L, 15, x15 = L, 16, x16 = L, 17, x17 = L, 18, x18 = L, 19, 
      x19 = L, 20, x20 = L, 21, x21 = L, 22, x22 = L];
     str = StringJoin["x", ToString[BrPreseka], "=", ToString[L], ";"];
     OpenWrite[f];
    WriteString[f, str];
    Write[f];
    Close[f];
    Get["NonABase.txt"];
    res = fCreatePData[StringJoin[ToString[BrPreseka], "*", tangle]]
    ]
    
fProdTangles[s1_String, s2_String, tt_String] := Module[{pol, pol1, L, LL},
    L = ToExpression[StringDrop[s1, -1]];
    LL = ToExpression[StringDrop[s2, -1]];
    pol = {{}, Union[L[[2]], LL[[2]] + 2Length[L[[2]]]], 
        Join[L[[3]], -LL[[3]]], L[[4]], LL[[4]] + 2Length[L[[2]]]};
    pol1 = 
      Map[Sort, {{pol[[4, 1, 1]], pol[[5, 1, 2]]}, {pol[[4, 2, 1]], 
            pol[[5, 1, 1]]}, {pol[[4, 1, 2]], 
            pol[[5, 2, 2]]}, {pol[[4, 2, 2]], pol[[5, 2, 1]]}}];
    pol = {{}, 
        Union[Complement[pol[[2]], Map[Sort, Union[pol[[4]], pol[[5]]]]], 
          pol1], pol[[3]]};
    PDataSumProdNonA[pol, tt]]
    
fSumTangles[s1_String, s2_String, tt_String] := Module[{pol, pol1, L, LL},
    L = ToExpression[StringDrop[s1, -1]];
    LL = ToExpression[StringDrop[s2, -1]];
    pol = {{}, Union[L[[2]], LL[[2]] + 2Length[L[[2]]]], 
        Join[-L[[3]], -LL[[3]]], L[[4]], LL[[4]] + 2Length[L[[2]]]};
    pol1 = 
      Map[Sort, {{pol[[4, 1, 1]], pol[[5, 2, 1]]}, {pol[[4, 2, 1]], 
            pol[[5, 1, 1]]}, {pol[[4, 1, 2]], 
            pol[[5, 2, 2]]}, {pol[[4, 2, 2]], pol[[5, 1, 2]]}}];
    pol = {{}, 
        Union[Complement[pol[[2]], Map[Sort, Union[pol[[4]], pol[[5]]]]], 
          pol1], pol[[3]]};
    PDataSumProdNonA[pol, tt]]
    
     (* ## ## ## ## ## ## # VIAE ## ## ## ## ## ## *)
     
     (*HAMCUBGR.M-Generate Hamilton cubic graphs*)
     (*13-DEC-1995 TAMARA BERTOK, CORRECTED S.J.*)
     
     AdmissibleEdge[x_,y_,n_,c_]:=
  If[Abs[x-y]=!=1&&Abs[x-y]=!=n-1&&y-x+1>=c,True,False]
  
ListOfOneFactors[n_] := 
  Module[{koncni, seznam, s, x, k, l, tmp, tmp1, tmp2, tmp3, c}, 
    s = Range[n];
    koncni = {};
    seznam = {};
    k = 2;
    While[k <= Length[s], If[s =!= {}, x = s[[1]]];
      If[x == 1 && k > n/2 + 1, k = Length[s] + 1, 
        Do[If[koncni =!= {}, c = koncni[[1, 2]], c = 0];
          
          If[AdmissibleEdge[x, s[[k]], n, c], 
            Do[koncni = Append[koncni, {x, s[[k]]}];
              s = Complement[s, {x, s[[k]]}];
              k = 2], 
            Do[If[k < Length[s], k = k + 1, 
                If[koncni =!= {}, Do[tmp = Last[koncni];
                    s = Union[s, tmp];
                    koncni = Complement[koncni, {tmp}];
                    l = Position[s, tmp[[2]]][[1, 1]];
                    k = l + 1], k = Length[s] + 1]]]];
          If[s == {}, seznam = Append[seznam, koncni];
            tmp1 = Last[koncni];
            tmp2 = koncni[[Length[koncni] - 1]];
            koncni = Complement[koncni, {tmp1, tmp2}];
            s = Union[s, tmp1, tmp2];
            l = Position[s, tmp2[[2]]][[1, 1]];
            k = l + 1;];
          If[k > Length[s] && Length[s] < n, tmp3 = Last[koncni];
            koncni = Complement[koncni, {tmp3}];
            s = Union[s, tmp3];
            l = Position[s, tmp3[[2]]][[1, 1]];
            k = l + 1;
            While[tmp3[[1]] =!= 1 && l == Length[s], Do[tmp3 = Last[koncni];
                s = Union[s, tmp3];
                koncni = Complement[koncni, {tmp3}];
                l = Position[s, tmp3[[2]]][[1, 1]];
                k = l + 1]]]]] (*If*)];(*While*)
       seznam = fDifSeq[seznam];
       seznam=Join[seznam,{Length[seznam]}];
     seznam
    ]
    
    
    
 fDiffViae[n_Integer] := 
  Module[{ g = {}, i, j, g1 = {}, m = 0, gnova = {}, pot, d, d1, dd},
    gnova = Drop[ListOfOneFactors[2n],-1];
   g = Table[
        Union[gnova[[i]], Map[Reverse, gnova[[i]]]], {i, Length[gnova]}];
    g = Flatten[
        Table[Table[g[[j, i, 2]] - g[[j, i, 1]], {i, 2n}], {j, Length[g]}]];
    g = Partition[
        Table[If[Abs[g[[i]]] > n, -Sign[g[[i]]] 2n + g[[i]], g[[i]]], {i, 
            Length[g]}], 2n];
    g = ReplaceAll[g, -n -> n];
    g = Table[
        First[Union[
            ReplaceAll[Table[RotateLeft[g[[i]], j], {j, 2n}], -n -> n], 
            ReplaceAll[
              Map[Reverse, Table[RotateLeft[g[[i]], j], {j, 2n}]], -n -> n], 
            ReplaceAll[-Table[RotateLeft[g[[i]], j], {j, 2n}], -n -> n], 
            ReplaceAll[-Map[Reverse, 
                  Table[RotateLeft[g[[i]], j], {j, 2n}]], -n -> n]]], {i, 
          Length[gnova]}];
    g1 = Union[g];
    g = Table[{g[[i]], gnova[[i]]}, {i, Length[gnova]}];
    g = Table[Select[g, #[[1]] == g1[[i]] &], {i, Length[g1]}];
    g = Sort[Table[First[g[[i]]][[2]], {i, Length[g]}]];
    d = Union[Table[{i, i + 1}, {i, 2n - 1}], {{1, 2n}}];
    g = Union[
        Table[If[PlanarQ[FromUnorderedPairs[Union[g[[i]], d]]] == True, 
            g[[i]], 0], {i, Length[g]}]];
    g = If[SameQ[g[[1]], 0], Drop[g, 1], g];
    g = Union[
        Table[If[VertexConnectivity
        [FromUnorderedPairs[Union[g[[i]], d]]] > 2,
             g[[i]], {}], {i, Length[g]}]];
    g = If[SameQ[g[[1]], {}], Drop[g, 1], g];
    g=Join[g,{Length[g]}];
    g
    ] 
    
    
    fViaToKL[LL_List]:=Module[{d,gg},
    d=Length[LL];
    gg=Join[Table[{i,i+1},{i,2d-1}],{{1,2d}}];
    gg=Sort[
        Map[Sort,
          Table[{Position[LL,gg[[i,1]]][[1,1]],
              Position[LL,gg[[i,2]]][[1,1]]},{i,Length[gg]}]]];
    gg=fKLfromGraph[gg];
    gg
    ]
    
    
    ShowChord[Ul_]:=Module[{k,t1,t2},
    k=Max[Flatten[Ul]];
    t1=FromUnorderedPairs[Ul];
    t2=Cycle[k];
    t1=ShowLabeledGraph[
        Highlight[
          FromUnorderedPairs[Union[ToUnorderedPairs[Cycle[k]],Ul]],{Ul}]];
    t1
    ]
    
    
    (*## ## ## ## ## ## ## ## ## ## ## # CONVERSIONS TO PD ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## ## ## ## #*)


checkArgs[s_,t_]:=
  ListQ[s]&&VectorQ[t,IntegerQ[#]&&#>=0&]&&
    Tr[t]<=Length[s]
    
    
iteratedTake[s_,t_]/;checkArgs[s,t]:=
  iteratedTake[s,t]=
    With[{w=FoldList[Plus,0,t]},
      Map[Take[s,#]&,Transpose[{Drop[w,-1]+1,Rest[w]}]]]
      
fConwayToPDDirProd[Ul_String] := Module[{ss, ss1, ss2, ss3, ll,ll1, i, j, k},
    ss = Sort[
        Join[Union[
            Flatten[StringPosition[Ul, "#"]]], {StringLength[Ul]}, {1}]];
    ss = Flatten[{{StringTake[Ul, {1, ss[[2]] - 1}]}, 
          Table[StringTake[Ul, {ss[[i]] + 1, ss[[i + 1]] - 1}], {i, 2, 
              Length[ss] - 2}], {StringTake[Ul, {ss[[-2]] + 1, ss[[-1]]}]}}];
    ss = Map[fConwayToPD, ss];
    ll = 2*Drop[Map[Length, ss], -1];
    ll = Table[Plus @@ Table[ll[[i]], {i, 1, j}], {j, Length[ll]}];
    ll1 = ll + 1;
    ss = Table[
        Table[ss[[k, i, j]], {i, Length[ss[[k]]]}, {j, 4}], {k, Length[ss]}];
    ss1 = ReplacePart[ss[[1]], ll[[1]] + 1, Last[Position[ss[[1]], 1]]];
    ss2 = 
      Flatten[Table[
          ReplacePart[ss[[i + 1]] + ll[[i]], ll1[[i + 1]], 
            Last[Position[ss[[i + 1]] + ll[[i]], ll1[[i]]]]], {i, 
            Length[ss] - 2}], 1];
    ss3 = 
      ReplacePart[ss[[-1]] + ll[[-1]], ll1[[-1]], 
        Last[Position[ss[[-1]] + ll[[-1]], ll1[[-1]]]]];
    ss3 = ReplacePart[ss3, 1, Last[Position[ss3, Last[ll1]]]];
    ss = Flatten[{ss1, ss2, ss3}, 1];
    ss1 = fConwayToPD[ToString[Length[ss]]];
    Do[ss1 = ReplacePart[ss1, ss[[i, j]], {i, j}], {i, Length[ss]}, {j, 4}];
    ss1
    ]    
    

fConToPD[Ul_String] := Module[{mm, nn, ss, vv, i}, 
mm = fGaussExtSigns[Ul];
    nn = fGaussExtSigns[StringReplace[Ul, "-" -> ""]];
    nn = Map[Sign, Flatten[mm]]*Map[Sign, Flatten[nn]];
    vv = Table[nn[[i]]*(-1)^i, {i, Length[nn]}];
    vv = Abs[Flatten[mm]]*vv;
    ss = Map[Length, mm];
    mm = If[MemberQ[ss, 0], {vv}, iteratedTake[vv, ss]];
    mm] 
    
    
    (* radi i direktni proizvod, Slavik, 20.08.2007 *)
    
 fConwayToPD[Ul_String] := 
  Module[{ss}, 
    ss = If[Not[SameQ[StringPosition[Ul, "#"], {}]], 
        ss = fConwayToPDDirProd[Ul], ss = fConToPD[Ul];
        ss = KnotTheory`PD[KnotTheory`GaussCode @@ ss]];
    ss]

(* fConwayToPD[Ul_String] := Module[{ss}, ss = fConToPD[Ul];
    ss = KnotTheory`PD[KnotTheory`GaussCode @@ ss];
    ss] *)


fKnDowToPD[Ul_List] := Module[{ss, gg, sc, i}, ss = Map[Sign, Ul[[2]]];
    gg = Map[Sort, Table[{2i - 1, Abs[Ul[[2, i]]]}, {i, Length[Ul[[2]]]}]];
    sc = Flatten[
        Complement[
          Table[If[ss[[i]] < 0, gg[[i]], {}], {i, Length[ss]}], {{}}]];
    gg = Map[Last, 
        Sort[Flatten[
            Table[{{gg[[i, 1]], i}, {gg[[i, 2]], i}}, {i, Length[gg]}], 1]]];
    gg = Table[gg[[i]]*(-1)^i, {i, Length[gg]}];
    gg = Table[If[MemberQ[sc, i], -gg[[i]], gg[[i]]], {i, Length[gg]}];
    gg = If[Length[Ul[[1]]] == 1, {gg}, iteratedTake[gg, 2Ul[[1]]]];
    gg]
    
    
fGenSignfromPDNew[Ul_]:=Module[{ss,res,i},
    ss=KnotTheory`PositiveQ /@Table[Ul[[i]],{i,Length[Ul]}];
    res=Table[If[SameQ[ss[[i]],True],1,-1],{i,Length[ss]}];
    res]




fKnotscapeDowToPD[Ulaz_] := Module[{ss,Ul}, 
If[SameQ[Head[Ulaz],String],Ul=ToExpression[Ulaz],
    Ul=Ulaz];
    ss = fKnDowToPD[Ul];
    ss = KnotTheory`PD[KnotTheory`GaussCode @@ ss];
    ss]

fDowkerToPD[Ulaz_] := Module[{ss,ss1,i,Ul},
If[SameQ[Head[Ulaz],String],Ul=ToExpression[Ulaz],
    Ul=Ulaz]; 
    ss = fSignsKL[Ul];
    ss1 = Map[Sign, Ul[[2]]]*Map[Sign, ss[[2]]];
    ss = {Ul[[1]], ss1*Ul[[2]]};
    ss = fKnDowToPD[ss];
    ss = KnotTheory`PD[KnotTheory`GaussCode @@ ss];
    ss]

fPdataToPD[Ul_List] := Module[{pd}, pd = fDowkerToPD[fDowfromPD[Ul]];
    pd]
   
   
    
    
    (* ## ## ## ## ## ## # fPDfromBW ## ## ## ## ## ## *)
    
    
    SignedIntegerBraidFromBraidWord[bw_String?SoberBraidWordQ] := 
    Characters[bw] /. 
      Thread[Join[CharacterRange["A", "Z"], CharacterRange["a", "z"]] -> 
          Join[-Range[26], Range[26]]];

BraidWordFromSignedIntegerBraid[sib_List] := 
    StringJoin[
      sib /. Thread[
          Join[-Range[26], Range[26]] -> 
            Join[CharacterRange["A", "Z"], CharacterRange["a", "z"]]]];


SignedIntegerBraidWidth[{}] = 1;

SignedIntegerBraidWidth[sib_List] := (Max[Abs[Union[sib]]] + 1) /; 
      SoberBraidWordQ[BraidWordFromSignedIntegerBraid[sib]];

SoberBraidWordQ[bw_] := ((Head[bw] === String) && LetterQ[bw]);

PlanarDiagramFromSignedIntegerBraidRepresentative[
      SIBR[n_Integer?Positive, sib_List]] := 
    Module[{a, b, c, d, perm, m, j, J, KnotTheory`Xp, 
    KnotTheory`Xm, pd, ar, 
cycles, 
closurerule, 
        indexes, len, loops}, perm = Range[n];
      m = n;
      cycles = 1;
      pd = sib /. j_Integer :> (J = Abs[j];
              a = perm[[J]]; b = perm[[J + 1]]; c = perm[[J]] = ++m; 
              d = perm[[J + 1]] = ++m;
              cycles = cycles*ar[a, d]*ar[b, c];
              If[(j > 0), KnotTheory`Xp[b, d, c, a], 
              KnotTheory`Xm[a, b, d, c]]);
      pd = Apply[KnotTheory`PD, pd];
      closurerule = MapThread[Rule, {perm, Range[n]}];
      cycles = 
        cycles /. 
              closurerule //. {ar[a_, b___, c_] ar[c_, d___, e_] :> 
                ar[a, b, c, d, e]} /. {arr_ar :> Rest[arr]};
      pd = pd /. closurerule;
      indexes = Flatten[Apply[List, Apply[List, cycles], {1}]];
      len = Length[indexes];
      loops = Length[Complement[Range[n], Abs[sib], Abs[sib] + 1]];
      pd = 
        Join[pd /. MapThread[Rule, {indexes, Range[len]}] /. 
        (KnotTheory`Xp | KnotTheory`Xm) -> \
KnotTheory`X, 
          Map[Loop, Apply[KnotTheory`PD, Range[len + 1, len + loops]]]];
      pd];

SignedIntegerBraidRepresentativeFromBraidWord[bw_String?SoberBraidWordQ] :=
    SIBR[StringIndexOfBraidWord[bw], SignedIntegerBraidFromBraidWord[bw]];

StringIndexOfBraidWord[bw_String?SoberBraidWordQ] := 
    SignedIntegerBraidWidth[SignedIntegerBraidFromBraidWord[bw]];


PlanarDiagram[L_] := 
    PlanarDiagramFromSignedIntegerBraidRepresentative[
        SignedIntegerBraidRepresentativeFromBraidWord[BraidWord[L]]] /; 
      SoberBraidWordQ[BraidWord[L]];
      
fPDfromBW[Ul_String] := Module[{ss, L},
    BraidWord[L] = Ul;
    ss = PlanarDiagram[L];
    ss
    ]  
    
    
      (* ## ## ## ## ## ## # LINKS-GOULD ## ## ## ## ## ## *)
    
    
LinksGould[k_Integer, n_Integer, Ul_String] := Module[{ss, L},
    BraidWord[L] = Ul;
    ss = LinksGouldInvariant[k, n, L][q, p];
    ss
    ]
    
fBRtoBW[Ul_]:=Module[{dd,i},
    dd=List@@Ul;
    dd=StringJoin[
        Table[FromCharacterCode[
            If[SameQ[Sign[dd[[2,i]]],1],Abs[dd[[2,i]]]+64,
              Abs[dd[[2,i]]]+96]],{i,Length[dd[[2]]]}]];
    dd
    ]
    

fBraidW[Ulaz_]:=Module[{ss,Ul},
If[SameQ[Head[Ulaz],String]&&SameQ[StringPosition[Ulaz,"{"],{}],Ul=Ulaz,
    Ul=ToExpression[Ulaz]];
    ss=KnotTheory`BR[
        If[SameQ[Head[Ul],String],KnotTheory`BR[fConwayToPD[Ul]],
          If[SameQ[Union[Map[EvenQ,Abs[Ul[[2]]]]],True],fDowkerToPD[Ul],
            fPdataToPD[Ul]]]];
    ss=fBRtoBW[ss];
    ss
    ]

    
LinksGouldInv[k_Integer,n_Integer,Ul_]:=Module[{ss,L},
    ss=fBraidW[Ul];
    ss=LinksGould[k,n,ss];
    ss
    ]
    
    
    
     (* ## ## ## ## ## ## # COMBINATORICA 42 ## ## ## ## ## ## *)
     
     Unprotect[ConnectedComponents]
      
ConnectedComponents[g_Graph] :=
        Block[{untraversed=Range[V[g]], visit, comps={}, \
e=ToAdjacencyLists[g],
              dfi=Table[0,{V[g]}],cnt=1, $RecursionLimit = Infinity},
              While[untraversed != {},
                    visit = {}; edges = {};
                    DFS[First[untraversed]];
                    AppendTo[comps,visit];
                    untraversed = Complement[untraversed,visit]        
              ];
              comps
        ] /; UndirectedQ[g] 
        
      
        
   
    Unprotect[DFS] 
               
DFS[v_Integer] :=
	( dfi[[v]] = cnt++;
	  AppendTo[visit,v];
	  Scan[ (If[dfi[[#]]==0,AppendTo[edges,{v,#}];DFS[#] ])&, e[[v]] ] )   
	  
	  	  
	  
fViaToKL[LL_List]:=Module[{d,gg},
    d=Length[LL];
    gg=Join[Table[{i,i+1},{i,2d-1}],{{1,2d}}];
    gg=Sort[
        Map[Sort,
          Table[{Position[LL,gg[[i,1]]][[1,1]],
              Position[LL,gg[[i,2]]][[1,1]]},{i,Length[gg]}]]];
    gg=fKLfromGraph[gg];
    gg
    ]
    
    
    
(* ## ## ## ## ## ## ## ## ## ## ## ## ## ListOfOneFactors ## ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## # *)

AdmissibleEdge[x_, y_, n_, c_] := 
  If[Abs[x - y] =!= 1 && Abs[x - y] =!= n - 1 && y - x + 1 >= c, True, False]
  
fDifSeq[Ul_List] := Module[{tt, tt1, uu, j},
    j = 2Length[Ul[[1]]];
    uu = Map[First, Union[Table[        
            tt = ReplaceAll[
                Mod[Table[Ul[[k]] + i - 1, {i, j}], j], {0 -> j}];
            tt = Sort[Table[Sort[Map[Sort, tt[[i]]]], 
            {i, Length[tt]}]];  
            tt1 = ReplaceAll[
                Mod[Table[Reverse[Ul[[k]]] + i - 1, {i, j}], j], {0 -> j}];
            tt1 = Sort[Table[Sort[Map[Sort, tt1[[i]]]], {i, Length[tt1]}]];
            tt = Union[tt, tt1], {k, Length[Ul]}]]];
            uu]
            

    
    
(* ## ## ## ## ## ## ## ## ## ## DIAGRAMS ## ## ## ## ## ## ## ## ## *)

fMulTan[ll_List, ll1_List] := Module[{n, vv, pp, i},
    n = Length[ll];
    vv = Range[n] + 1;
    vv = Union[ll, ll1 + 2n, 
        Join[{{vv[[1]], 2n + 1}}, 
          Table[{vv[[i]], 4n + 3 - vv[[i]]}, {i, 2, n}]]];
    vv = ConnectedComponents[FromUnorderedPairs[vv]];
    pp = Union[Complement[Range[2n], Range[n] + 1],
        Range[n] + 2n + 1];
    vv = Complement[
        Map[Sort, 
          ReplaceAll[
            Mod[Table[Intersection[pp, vv[[i]]], {i, Length[vv]}], 
              2n], {0 -> 2n}]], {{}}];
    vv
    ]
    
    
    
fGenSet[Ul_List] := Module[{vv, i, j},
    vv = Ul;
    While[
      Length[Union[vv, 
            Flatten[Table[
                fMulTan[vv[[i]], vv[[j]]], 
                {i, Length[vv]}, {j, Length[vv]}], 1]]] > Length[vv],
      vv = 
        Union[vv, 
          Flatten[Table[
              fMulTan[vv[[i]], vv[[j]]], {i, Length[vv]}, {j, Length[vv]}], 
            1]]];
    vv = {vv, Length[vv], 
        If[SameQ[Length[vv], Factorial2[2Length[Ul[[1]]] - 1]], 1, 0]}
    ]
    
    
fMulTanTab[Ul_List] := Module[{tt, tt1, i, j},
    tt = fGenSet[Ul][[1]];
    Print[tt];
    tt1 = tt;
    tt = Flatten[
        Table[Table[{i, j, 
              Flatten[Position[tt, fMulTan[tt[[i]], tt[[j]]]]][[1]]}, {i, 
              Length[tt]}], {j, Length[tt]}], 1];
    tt = MatrixForm[
        Partition[Table[tt[[i, 3]], {i, Length[tt]}], Length[tt1]]];
    tt
    ]
    
    
(* ## ## ## ## ## ## ## ## # SPIDER ## ## ## ## ## ## ## ## ## ## ## ## # *)

fMakeTan[n_Integer] := Module[{tt, i},
    If[SameQ[n, 3], tt = {{1, 2}, {1, 4}},
      tt = Table[2, {i, n - 2}];
      tt = Join[tt, {1, 1}];
      tt = Join[tt, tt];
      tt = ReplacePart[tt, 2, -1];
      tt = Drop[Flatten[Position[tt, 2]], 1];
      tt = Table[{1, tt[[i]]}, {i, Length[tt]}]];
    tt]

fMinTangle[n_Integer, p_Integer] := Module[{pp, uu, ss, k, i, j},
    uu = fMakeTan[n];
    ss = Range[2n];
    Do[pp = Flatten[Table[
            pp = uu[[i]];
            
            Complement[
              Union[Table[
                  If[SameQ[ss[[j]], pp[[-2]]] || 
                      SameQ[ss[[j]], pp[[-1]] - 1] || ss[[j]] > pp[[-1]], 
                    Join[pp, {ss[[j]]}], pp], {j, Length[ss]}]], {pp}], {i, 
              Length[uu]}], 1];
      uu = pp, {k, p - n - 1}];
    pp = Table[{n, pp[[i]]}, {i, Length[pp]}];
    pp]

fMakeTangle[n_Integer] := Module[{tt, i},
    tt = Table[2, {i, n - 2}];
    tt = Join[tt, {1, 1}];
    tt = Join[tt, tt];
    tt = {{}, tt};
    tt]

fRed3[Ul_List] := Module[{tt, i},
    tt = Table[If[Ul[[i]] >= 3, 3, Ul[[i]]], {i, Length[Ul]}];
    tt
    ]

fTanglefromCode[Ul_List] := Module[{dd,qq, i, j},
    dd = 2Ul[[1]];
    qq = fMakeTangle[Ul[[1]]];
    qq = Table[i = Ul[[2, j]];
        qq = If[SameQ[qq, {}] || qq[[2, i]] == 1, {},
            {{i}, 
              fRed3[If[SameQ[i, 1], 
                  ReplacePart[
                    ReplacePart[ReplacePart[qq[[2]], qq[[2]][[-1]] + 1, -1], 
                      qq[[2]][[2]] + 1, i + 1], 1 , i], 
                  If[SameQ[i, dd], 
                    ReplacePart[
                      ReplacePart[
                        ReplacePart[qq[[2]], qq[[2]][[i - 1]] + 1, dd - 1], 
                        qq[[2]][[1]] + 1, 1], 1 , i], 
                    ReplacePart[
                      ReplacePart[
                        ReplacePart[qq[[2]], qq[[2]][[i - 1]] + 1, i - 1], 
                        qq[[2]][[i + 1]] + 1, i + 1], 1 , i]]]]}],
        {j, Length[Ul[[2]]]}];
    qq = If[MemberQ[qq, {}], {}, {Ul[[1]], Ul[[2]], qq[[-1, -1]]}];
    qq
    ]

fMakeEdges[Ul_List] := Module[{hh, i},
    hh = Sort[
        Map[Sort, 
          Union[Table[{Ul[[i]], Ul[[i + 1]]}, {i, 
                Length[Ul] - 1}], {{Ul[[-1]], Ul[[1]]}}]]];
    hh
    ]

fRegions[Ul_List, n_Integer] := Module[{hh, i},
    hh = Sort[
        Table[If[SameQ[Ul[[i]], {1, 2n}], 2n, Min[Ul[[i]]]], {i, 
            Length[Ul]}]];
    hh
    ]

fRegList[Ul_List, n_Integer] := Module[{hh, hh1, i},
    hh = Range[2n];
    hh = fMakeEdges[hh];
    hh = Complement[
        Table[If[MemberQ[hh, Ul[[1, i]]], Ul[[1, i]], {}], {i, 
            Length[Ul[[1]]]}], {{}}];
    hh1 = Map[fMakeEdges, Ul[[2]]];
    hh = Union[Table[Complement[hh1[[i]], Ul[[1]]], {i, Length[hh1]}], 
        Table[{hh[[i]]}, {i, Length[hh]}]];
    hh = Table[fRegions[hh[[i]], n], {i, Length[hh]}];
    hh = {Ul[[1]], Sort[hh]};
    hh
    ]

fAllClosures[n_Integer] := Module[{tt, tt1, tt2, hh, hh1, i},
    tt = Table[2i - 1, {i, n}];
    tt1 = Table[2i, {i, n}];
    tt1 = Permutations[tt1];
    tt1 = 
      Map[Sort, 
        Table[Map[Sort, Table[{tt[[j]], tt1[[i, j]]}, {j, Length[tt]}]], {i, 
            Length[tt1]}]];
    tt2 = Flatten[Position[Map[Length, Map[Union, Map[Flatten, tt1]]], 2n]];
    tt2 = Range[2n];
    tt2 = 
      Union[Table[{tt2[[i]], tt2[[i + 1]]}, {i, Length[tt2] - 1}], \
{{tt2[[1]],
             tt2[[-1]]}}];
    tt = Complement[
        Union[Table[
            If[PlanarQ[FromUnorderedPairs[Union[tt1[[i]], tt2]]], 
              tt1[[i]], {}], {i, Length[tt1]}]], {{}}];
    hh = Range[2n];
    hh = Table[
        Union[Table[{hh[[i]], hh[[i + 1]]}, {i, Length[hh] - 1}], {{1, 2n}}, 
          tt[[i]]], {i, Length[tt]}];
    hh = Map[fPlanarEmbGraph, hh];
    hh = Complement[
        Table[If[
            Length[hh[[i, -1, 1]]] == 2n, {tt[[i]], 
              Drop[hh[[i, -1]], 1]}, {}], {i, Length[hh]}], {{}}];
    hh = Table[fRegList[hh[[i]], n], {i, Length[hh]}];
    hh=Map[First,hh];
    hh=Join[hh,{Length[hh]}];
    hh
    ]

fClosedTangle[Ul_List, LL_List] := Module[{cc1, i, j, k},
    cc1 = 
      Table[{Table[
            Table[Ul[[3, LL[[k, 2, j, i]]]], 
            {i, Length[LL[[k, 2, j]]]}], {j, 
              Length[LL[[k, 2]]]}], k}, {k, Length[LL]}];
    cc1 = 
      Flatten[Complement[
          Union[Table[
              If[SameQ[
                  Union[fRed3[
                      Table[Plus @@ cc1[[j, 1, i]], {i, 
                          Length[cc1[[j, 1]]]}]]], {3}], {Drop[Ul, -1], 
                  j}, {}], {j, Length[cc1]}]], {{}}], 1];
     cc1 = If[Length[cc1] > 1, Partition[cc1, 2], cc1];
    cc1
    ]

fBasicTan[n_Integer, k_Integer] := Module[{rr, cc, i},
    rr = fMinTangle[n, k];
    rr = Complement[Map[fTanglefromCode, rr], {{}}];
    (* cc = fAllClosures[n]; *)
    cc = cllist[[n - 2]];
    cc = Sort[
        Flatten[Complement[
            Union[Table[fClosedTangle[rr[[i]], cc], 
            {i, Length[rr]}]], {{}}], 
          1]];
    cc]

fInConnect[Ul_List] := Module[{tt, i},
    tt = If[Length[Ul] < 3, {Ul}, 
        Table[{Ul[[i]], Ul[[i + 1]]}, {i, Length[Ul] - 1}]];
    tt
    ]

fGrEdgGen[Ul_List, LL_List] := Module[{pp, str, ff, ll, mm, ss, cl, i},
    pp = Range[2Ul[[1, 1]]];
    pp = Map[Flatten, Table[Position[Ul[[1, 2]], pp[[i]]], 
    {i, Length[pp]}]] + Ul[[1, 1]] - 1;
    str = Join[{{2Ul[[1, 1]], 1}}, fInConnect[Range[2Ul[[1, 1]]]]] ;
    str = 
      Map[Sort, 
        Table[Union[Flatten[{pp[[str[[i, 1]]]], pp[[str[[i, 2]]]]}]], {i, 
            Length[str]}]];
    ff = Table[If[SameQ[str[[i]], {}], {}, Min[str[[i]]]], {i, Length[str]}];
    ll = Table[If[SameQ[str[[i]], {}], {}, Max[str[[i]]]], {i, Length[str]}];
    mm = Union[Flatten[Map[fInConnect, str], 1]];
    mm = Select[mm, SameQ[Length[#], 2] &];
    ss = Flatten[{Range[Ul[[1, 1]] - 1], Ul[[1, 1]] - 1, 
          Reverse[Range[Ul[[1, 1]] - 1]], 1}];
    ff = Table[If[SameQ[ff[[i]], {}], ss[[i]], ff[[i]]], {i, Length[ff]}];
    ll = Table[If[SameQ[ll[[i]], {}], ss[[i]], ll[[i]]], {i, Length[ll]}];
    ff = Union[Table[{ss[[i]], ff[[i]]}, {i, Length[ss]}], 
        fInConnect[Range[Ul[[1, 1]] - 1]]];
    ff = Select[ff, Not[SameQ[#[[1]], #[[2]]]] &];
    cl = Flatten[LL[[Ul[[-1]]]]];
    cl = Union[
        Map[Sort, Partition[Table[ll[[cl[[i]]]], {i, Length[cl]}], 2]]];
    cl = Union[ff, mm, cl];
    cl = Select[cl, Not[SameQ[#[[1]], #[[2]]]] &];
    cl
    ]

fDowfromBP[Ul_List] := Module[{cc, res},
    cc = cllist[[Ul[[1, 1]] - 2]];
    cc = Map[First, cc];
    res = fGrEdgGen[Ul, cc];
    res = fKLfromGraph[res];
    res = MinDowProjAltKL[res];
    res
    ]

fDowfromCode[Ul_List] := Module[{kk},
    kk = fTanglefromCode[Ul];
    kk = If[Not[SameQ[kk, {}]],
        kk = Flatten[fClosedTangle[kk, cllist[[Ul[[1]] - 2]]], 1];
        fDowfromBP[kk], kk];
    kk
    ]

fBasicPolyTan[k_Integer, n_Integer] := Module[{po, ff, ff1, cc, i},
    po = fBasicTan[k, n];
    cc = Map[First, cllist[[k - 2]]];
    ff = Table[
        MinDowProjAltKL[fKLfromGraph[fGrEdgGen[po[[i]], cc]]], {i, 
          Length[po]}];
    ff1 = Union[ff];
    ff = Flatten[
        Map[First, Table[Position[ff, ff1[[i]]], {i, Length[ff1]}]]];
    ff = Sort[Table[po[[ff[[i]]]], {i, Length[ff]}]];
    ff
    ]

fMakeBP[LL_List, n_Integer, s_Integer, k_Integer, l_Integer] := 
  Module[{kk, kk1, kk2, dd, dd1, i},
    kk = LL[[n - 2, s]];
    kk1 = {n, 
        ReplaceAll[Join[Mod[Range[k], 2], kk, Mod[Range[l], 2]], {0 -> 2}]};
    kk2 = {n, 
        ReplaceAll[
          Join[Mod[Range[k], 2], kk, Mod[Range[l], 2] + 1], {0 -> 2}]};
    kk1 = fTanglefromCode[kk1];
    kk1 = If[SameQ[kk1, {}], kk1,
        Flatten[fClosedTangle[kk1, cllist[[n - 2]]], 1]];
    kk1 = If[Length[kk1] == 2, {kk1}, Partition[kk1, 2]];
    kk1 = Table[{kk1[[i]], fDowfromBP[kk1[[i]]]}, {i, Length[kk1]}];
    kk2 = fTanglefromCode[kk2];
    kk2 = If[SameQ[kk2, {}], kk2,
        Flatten[fClosedTangle[kk2, cllist[[n - 2]]], 1]];
    kk2 = If[Length[kk2] == 2, {kk2}, Partition[kk2, 2]];
    kk2 = Table[{kk2[[i]], fDowfromBP[kk2[[i]]]}, {i, Length[kk2]}];
    kk = Complement[Union[{kk1, kk2}], {{}}];
    dd = Flatten[kk, 1];
    dd1 = Union[Table[dd[[i, 2]], {i, Length[dd]}]];
    dd1 = Table[Position[dd, dd1[[i]]], {i, Length[dd1]}];
    kk = Map[First, Map[First, dd1]];
    dd = Table[dd[[kk[[i]]]], {i, Length[kk]}];
    dd
    ]

fMakeAllnsBP[LL_List, n_Integer, s_Integer, c_Integer] := 
  Module[{kk, dd, dd1, i},
    kk = lllll[[n - 2, s]];
    dd = Length[kk];
    dd = Range[c - dd - n + 1];
    dd1 = c - dd - Length[kk] - n + 1;
    dd = Complement[
        Table[fMakeBP[lllll, n, s, dd[[i]], dd1[[i]]], {i, 
            Length[dd]}], {{}}];
    dd = Complement[Flatten[dd, 1], {{}}];
     dd1 = Union[Table[dd[[i, 2]], {i, Length[dd]}]];
    dd1 = Table[Position[dd, dd1[[i]]], {i, Length[dd1]}];
    kk = Map[First, Map[First, dd1]];
    dd = Table[dd[[kk[[i]]]], {i, Length[kk]}];
    dd=Join[dd,{Length[dd]}]; 
    dd
    ]
    
    
fKLTab[Ul_]:=Module[{tt},
    tt=ToExpression[Ul];
    tt
    ]

  
    
(* ## ## ## ## ## ## ## ## ## ## VIRTUAL KLs ## ## ## ## ## ## ## ## ## *)



rule1 = {X[a_, b_, c_, d_] :> A del[d, a] del[c, b] +
				B led[a, b] led[c, d],
      Y[a_, b_, c_, d_] :> B del[d, a] del[c, b] +
				A led[a, b] led[c, d]};

rule2 = {del[a_, b_] del[b_, c_] :> del[a, c], 
      del[a_, b_] del[c_, b_] :> del[a, c],
      del[a_, b_] led[b_, c_] :> led[a, c], 
      del[a_, b_] led[c_, b_] :> led[c, a],
      del[b_, a_] led[b_, c_] :> led[a, c], 
      del[b_, a_] led[c_, b_] :> led[c, a],
      led[a_, b_] del[b_, c_] :> led[a, c], 
      led[a_, b_] del[c_, b_] :> led[a, c],
      led[b_, a_] del[c_, b_] :> led[a, c], 
      led[b_, a_] del[b_, c_] :> led[a, c],
      led[a_, b_] led[b_, c_] :> del[a, c]
      };

rule3 = {(del[a_, a_]) :> J, (del[a_, b_])^2 :> J,
      led[a_, b_]^2 :> J K1,
      led[a_, b_]led[c_, b_]led[c_, d_]led[a_, d_] :> J K2,
       led[a_, g_] led[a_, i_] led[c_, g_] led[c_, k_] led[e_, i_] led[e_, 
            k_] :> J K3};

RawBracket[t_] := Simplify[(t /. rule1 // Expand) //. rule2 /. rule3]

rule4 = {B :> 1/A, J :> -A^2 - 1/A^2};

B[t_] := Simplify[RawBracket[RawBracket[t]/J] /. rule4]

F[t_] := B[t] /. A :> 1


fReplacement[Ul_List] := Module[{vv,ii},
    vv = Table[
        Ul[[ii]] -> ToExpression[FromCharacterCode[Ul[[ii]] + 96]], {ii, 
          Length[Ul]}];
    vv]

fReplacement1[pp_List,mm_List]:=Module[{pp1,res,ii},
    pp1=Table[
        If[SameQ[mm[[ii]],-1],
          pp[[ii]],{pp[[ii,4]],pp[[ii,1]],pp[[ii,2]],pp[[ii,3]]}],{ii,Length[pp]}];
    res=Table[If[SameQ[mm[[ii]],-1],Y@@pp1[[ii]],X@@pp1[[ii]]],{ii,Length[mm]}];
    res=Product[res[[ii]],{ii,Length[res]}];
    res]

(* CORRECTED 29.04.2010 *)

fKauffmanPD[Ul_String]:=Module[{pp,pp1,pp2,pp3,rr,ss,mm,res,ii,jj},
pp=fVirtPD[Ul];
ss=fGenSignfromPDNew[pp];
mm=Union[Flatten[Table[pp[[ii,jj]],{ii,Length[pp]},{jj,4}]]];
pp=ReplaceAll[Table[Table[pp[[jj,ii]],{ii,4}],{jj,Length[pp]}],
        fReplacement[mm]];
pp=Table[ReplaceAll[pp[[ii]],{pp[[ii,1]]->pp[[ii,3]],pp[[ii,3]]->pp[[ii,1]]}],{ii,Length[pp]}]; (* ADDED *)
pp=fReplacement1[pp,ss];
pp1=Union[Flatten[Table[Table[pp[[jj,ii]],{ii,Length[pp[[jj]]]}],{jj,Length[pp]}]]];
pp2=Map[ToString,pp1];
pp3=Flatten[Table[ToCharacterCode[pp2[[ii]]],{ii,Length[pp2]}]]-96;
rr=Table[pp1[[ii]]->pp3[[ii]],{ii,Length[pp3]}];
res=ReplaceAll[pp,rr];
res=PD@@Table[res[[ii]],{ii,Length[ss]}];
res] 
    
    
    
fKauffmanBracket[Ul_]:=Module[{pp,ss,mm,res,i,j},
    pp=If[SameQ[Head[Ul],String],fVirtPD[Ul],Ul];
    ss=fGenSignfromPDNew[pp];
    mm=Union[Flatten[Table[pp[[i,j]],{i,Length[pp]},{j,4}]]];
    pp=ReplaceAll[Table[Table[pp[[j,i]],{i,4}],{j,Length[pp]}],
        fReplacement[mm]];
    pp=fReplacement1[pp,ss];
    res={Expand[B[pp]],F[pp]};
    res
    ]

    
fKauffmanBracket1[Ul_]:=Module[{pp0,mm1,pp,pp1,mm,res,res1,i,j},
    pp=If[SameQ[Head[Ul],String],fVirtPD[Ul],Ul];
    mm=Union[Flatten[Table[pp[[i,j]],{i,Length[pp]},{j,4}]]];
    pp=ReplaceAll[Table[Table[pp[[j,i]],{i,4}],{j,Length[pp]}],
        fReplacement[mm]];
    pp0=pp;
    mm=fGenSignsPD[Ul][[2]];
    mm=Select[mm,Not[SameQ[#,0]]&];
    mm1=-mm;
    pp=fReplacement1[pp,mm];
    res={RawBracket[pp],Expand[B[pp]],F[pp]};
    pp1=fReplacement1[pp0,mm1];
    res1={RawBracket[pp1],Expand[B[pp1]],F[pp1]};
    res={res,res1};
    res=Select[res,
        SameQ[Union[Position[#[[2]],led],Position[#[[2]],del]],{}] &];
    res=If[SameQ[res,{}],{},Rest[res[[1]]]];
    res
    ]



fKauffmanExtendedBracket[Ul_]:=Module[{res},
    res=fKauffmanBracket[Ul];
    res=If[SameQ[Union[Position[res,led],Position[res,del]],{}],res,
        fKauffmanBracket1[Ul]];
    res
    ]


fKauffmanArrow[Ul_]:=Module[{res},
    res=fKauffmanBracket[Ul];
    res=If[SameQ[Union[Position[res,led],Position[res,del]],{}],res,
        fKauffmanBracket1[Ul]];
    res
    ]


fGenSignsPD[Ul_String] := 
  Module[{uu, uu1, uu0, rr, rr1, cc, cc1,res1, cc0, pp0, pp, pp1, pp2, 
vv, vv1, ll, tt, tt2, res, i, j, k}, 
uu = StringReplace[Ul, {"i" -> "1", "-1" -> "1"}];
    rr = fMakePolyhedral[uu];
    rr1 = StringPosition[rr, "1"];
    cc0 = fConwayToPD[uu];
    pp0 = Table[cc0[[i, j]], {i, Length[cc0]}, {j, 4}];
    ll = Length[pp0];
    tt = Table[StringReplacePart[rr, "-1", rr1[[i]]], {i, Length[rr1]}];
    uu0 = fGenSign[uu];
(*znaci u tackama orginalnog alternirajuceg linka*)
pp = Map[Last, 
        Position[Abs[Table[fGenSign[tt[[i]]] - uu0, {i, Length[tt]}]], 2]];
    pp2 = pp;
    pp = Table[uu0[[pp[[i]]]], {i, Length[pp]}];
    vv = StringReplace[Ul, {"i" -> "1"}];
    vv1 = fGenSign[vv];
    tt2 = Table[fGenSign[tt[[i]]] - uu0, {i, Length[tt]}];
    tt2 = Map[fConwayToPD, tt];
    tt2 = Table[cc = tt2[[k]];
        cc = Table[cc[[i, j]], {i, Length[cc]}, {j, 4}] - pp0, {k, 
          Length[tt2]}];
    tt2 = 
      Flatten[Table[
          Complement[Range[ll], 
            Flatten[Position[tt2[[i]], {0, 0, 0, 0}]]], {i, Length[tt2]}]];
    res = Map[Last, Sort[Table[{tt2[[i]], pp[[i]]}, {i, Length[pp]}]]];
    uu1 = fGenSign[StringReplace[Ul, {"i" -> "1"}]];
    pp1 = Table[uu1[[pp2[[i]]]], {i, Length[pp2]}];
    res = Map[Last, Sort[Table[{tt2[[i]], pp1[[i]]}, {i, Length[pp1]}]]];
    uu = StringReplace[Ul, "i" -> "1"];
    cc0 = fConwayToPD[uu];
    pp0 = Table[cc0[[i, j]], {i, Length[cc0]}, {j, 4}];
    uu1 = StringReplace[Ul, "i" -> "-1"];
    cc1 = fConwayToPD[uu1];
    pp1 = Table[cc1[[i, j]], {i, Length[cc1]}, {j, 4}];
    pp = pp0 - pp1;
    pp = Complement[Range[ll], Flatten[Position[pp, {0, 0, 0, 0}]]];
    res1 = Table[If[MemberQ[pp, i], 0, res[[i]]], {i, Length[res]}];
    res = {res, res1};
    res]



fGenSignfromPD[Ul_] := 
  Module[{res,pp0, gg, mm, mm1, mm2, gg1, gg2, i, j}, 
    pp0 = Table[Ul[[i, j]], {i, Length[Ul]}, {j, 4}];
    gg = GaussCode[Ul];
    res = Table[mm = RotateLeft[pp0[[j]]];
        mm1 = ReplacePart[pp0, mm, j];
        mm2 = PD @@ Table[X @@ mm1[[i]], {i, Length[mm1]}];
        gg = GaussCode[mm2];
        gg1 = Sort[Flatten[Table[gg[[i]], {i, Length[gg]}]]];
        gg2 = Join[Range[-Length[Ul], -1], Range[1, Length[Ul]]];
        If[SameQ[gg1, gg2], 1, -1], {j, Length[Ul]}];
    res]




fMakeConVirtLinkKauff[Ul_String] := Module[{nn, tt, tt1, tt2, tt3, tt4, ll, \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i},
    nn = If[SameQ[StringPosition[Ul, "#"], {}], Plus @@ \
fCreatePData[Ul][[1]],
         Plus @@ fDToDDirectPD[Ul][[1]]];
    ll = fBasVirt[nn];
    tt = Table[fMakeConV[Ul, ll[[i]]], {i, Length[ll]}];
    tt = Select[tt, 
        SameQ[StringPosition[#, "i,i"], {}] && SameQ[StringPosition[#, "(-1,
1"], {}] && SameQ[StringPosition[#, "-1,1)"], {}] && 
            SameQ[StringPosition[#, "(1,-1"], {}] && 
            SameQ[StringPosition[#, ",1,-1)"], {}] && 
            SameQ[StringPosition[#, ",-1,1,"], {}] && 
            SameQ[StringPosition[#, ",1,-1,"], {}] &];
tt=Select[tt,Not[SameQ[StringPosition[#,"i"],{}]] &];
    tt1 = Table[fKauffmanExtendedBracket[tt[[i]]], {i, Length[tt]}];
    tt2 = Complement[Union[tt1], {{J (A+B K1),A+K1/A,1+K1}}];
    tt3 = 
      Map[First, 
        Table[Flatten[Position[tt1, tt2[[i]]], 1], {i, Length[tt2]}]];
    tt3 = Table[{tt1[[tt3[[i]]]], tt[[tt3[[i]]]]}, {i, Length[tt3]}];
    tt3]
  

VirtKnotTab[n_Integer]:=Module[{pp1,pp2,pp},
    pp=virtknotab[[n]];
    pp1=Plus@@{Length[StringPosition[pp,"1"]],
          Length[StringPosition[pp,"i"]]};
    pp2=Length[StringPosition[pp,"i"]];
    pp={pp1,pp1-pp2,pp};
    pp]
    
    
VirtLinkTab[n_Integer]:=Module[{pp1,pp2,pp},
    pp=virtlinktable[[n]];
    pp1=Plus@@{Length[StringPosition[pp,"1"]],
          Length[StringPosition[pp,"i"]]};
    pp2=Length[StringPosition[pp,"i"]];
    pp={pp1,pp1-pp2,pp};
    pp]


fPolyNorm[Ul_] := Module[{tt, tt1, ll, i},
    tt = If[SameQ[Ul, -1], 1, 
        If[SameQ[Ul, 0], 0, 
          If[SameQ[CoefficientList[Ul, x][[1]], 0], 
            Abs[CoefficientList[Ul, x][[-1]]],
            tt = Table[Denominator[Ul[[i]]], {i, Length[Ul]}];
            tt = If[SameQ[tt, {}], Ul, tt1 = PolynomialLCM @@ tt;
                tt = Expand[Ul*tt1];
                
                ll = If[
                    SameQ[NumberQ[tt[[1]]], True] && 
                        SameQ[Sign[tt[[1]]], -1] || 
                      Length[tt[[1]]] > 0 && SameQ[Sign[tt[[1, 1]]], -1], \
-tt,
                     tt]];
             tt = If[SameQ[tt, 1], tt, 
                If[SameQ[Length[Variables[tt]], 1], 
                  Expand[Divide[tt, 
                      PolynomialGCD @@ Table[tt[[i]], {i, Length[tt]}]]], 
                  tt = ReplaceAll[tt, az -> a*z];
                  tt1 = PolynomialGCD @@ Table[tt[[i]], {i, Length[tt]}];
                  tt = Expand[Divide[tt, tt1]]]]]]];
                 
    tt]
    
    
fFindConVirt[Ul_String]:=Module[{pp1,pp},
pp=fCabledJonesVirt[Ul,2];
pp1=Flatten[Position[virttab,pp]];
pp1=If[SameQ[pp1,{}],{},virtknotab[[First[pp1]]]];
pp1
]

fVirtLinkSketch[]:=Module[{ss},
Run["C:/LinKnot/LinkSketcher4.0.jar"];
ss=ToExpression[Import["linkSketcherOutput.txt"]]]

fFindKSketch[]:=Module[{ss,pp,res},
    ss=fVirtLSketch[];
pp=fCabledJonesVirtPD[ss[[2]],ss[[1]],3];
pp=Flatten[Position[virttab,pp]];
res=If[SameQ[pp,{}],{},virtknotab[[First[pp]]]];
res]


fVirtPD[Con_String] := Module[{cc, cc1, cc2, pp, pp1, pp2, i, j},
     cc1 = StringReplace[Con, "i" -> "1"];
    cc2 = StringReplace[Con, "i" -> "-1"];
    cc1 = fConwayToPD[cc1];
    cc2 = fConwayToPD[cc2];
    pp1 = Table[cc1[[i, j]], {i, Length[cc1]}, {j, 4}];
    pp2 = Table[cc2[[i, j]], {i, Length[cc1]}, {j, 4}];
    cc = pp1 - pp2;
    cc = Flatten[Position[cc, {0, 0, 0, 0}]];
    If[SameQ[cc, {}], pp = {},
      cc1 = cc;
      cc = Table[pp1[[i]], {i, Length[pp1]}];
      pp1 = 
        Complement[
          Table[If[
              SameQ[pp1[[i]], 
                pp2[[i]]], {}, {{pp1[[i, 4]], pp1[[i, 2]]}, {pp1[[i, 3]], 
                  pp1[[i, 1]]}}], {i, Length[pp1]}], {{}}];
      Do[cc = 
          ReplaceAll[
            cc, {pp1[[i, 1, 1]] -> pp1[[i, 1, 2]], 
              pp1[[i, 2, 1]] -> pp1[[i, 2, 2]]}]; 
        pp1 = ReplaceAll[
            pp1, {pp1[[i, 1, 1]] -> pp1[[i, 1, 2]], 
              pp1[[i, 2, 1]] -> pp1[[i, 2, 2]]}], {i, Length[pp1]}];
      cc = Table[cc[[cc1[[i]]]], {i, Length[cc1]}];
      pp = Union[Flatten[cc]];
      Do[cc = ReplaceAll[cc, pp[[i]] -> i], {i, Length[pp]}];
      pp = fConwayToPD[ToString[Length[cc]]];
      Do[pp = ReplacePart[pp, cc[[i, j]], {i, j}], {i, Length[cc]}, {j, 4}]];
    pp]

fFindSignsPart[Ul_String] := Module[{ss, ss1, i},
    ss = Flatten[fGaussExtSigns[Ul]];
    ss1 = Abs[ss];
    ss = Map[Sign, 
        Map[Last, Union[Table[{ss1[[i]], ss[[i]]}, {i, Length[ss]}]]]];
    ss]

fFindSigns[Ul_String] := Module[{ss, i},
    If[SameQ[StringPosition[Ul, "#"], {}], ss = fFindSignsPart[Ul], 
      ss = Sort[
          Join[Union[
              Flatten[StringPosition[Ul, "#"]]], {StringLength[Ul]}, {1}]];
      ss = 
        Flatten[{{StringTake[Ul, {1, ss[[2]] - 1}]}, 
            Table[StringTake[Ul, {ss[[i]] + 1, ss[[i + 1]] - 1}], {i, 2, 
                Length[ss] - 2}], {StringTake[
                Ul, {ss[[-2]] + 1, ss[[-1]]}]}}];
      ss = Flatten[Map[fFindSignsPart, ss]]];
    ss]

fSawollek[Ul_String] := 
  Module[{cc, cc1, cc2, pp, pp1, pp2, vv, ll, mm, mm1, mm2, nn, nn1, nn2, \
nn3,
       ss, i, j},
    cc1 = StringReplace[Ul, "i" -> "1"];
    cc2 = StringReplace[Ul, "i" -> "-1"];
    cc1 = fConwayToPD[cc1];
    cc2 = fConwayToPD[cc2]; 
    pp1 = Table[cc1[[i, j]], {i, Length[cc1]}, {j, 4}];
    pp2 = Table[cc2[[i, j]], {i, Length[cc1]}, {j, 4}];
    cc = pp1 - pp2;
    cc = Flatten[Position[cc, {0, 0, 0, 0}]];
    cc = Table[If[MemberQ[cc, i], 1, 0], {i, Length[cc1]}];
     pp = fFindSigns[StringReplace[Ul, {"i" -> "1"}]]; 
    pp = cc*pp;
    pp = Select[pp, Not[SameQ[#, 0]] &];
    vv = Plus @@ pp;
    ll = Length[pp];
    mm = IdentityMatrix[ll];
    mm = ReplaceAll[mm, {0 -> {0, 0}}];
    Do[mm = ReplacePart[mm, pp[[i]], {i, i}], {i, ll}];
    mm1 = ReplaceAll[mm, {1 -> {1 - x, -y}, -1 -> {0, -y*x^{-1}}}];
    mm2 = ReplaceAll[mm, {1 -> {-x*y^{-1}, 0}, -1 -> {-y^{-1}, 1 - x^{-1}}}];
    mm = Map[Flatten, Flatten[Table[{mm1[[i]], mm2[[i]]}, {i, ll}], 1]];
    nn = fVirtPD[Ul];
    nn1 = Flatten[Table[nn[[i, j]], {i, Length[nn]}, {j, 4}]];
    nn2 = 
      Flatten[Table[
          If[SameQ[pp[[i]], 1], {-2i + 1, 2i - 1, 2i, -2i}, {-2i, -2i + 1, 
              2i - 1, 2i}], {i, Length[pp]}]];
    ss = nn2;
    nn2 = Table[{nn1[[i]], nn2[[i]]}, {i, Length[nn1]}];
    nn3 = Map[First, nn2];
    nn3 = Flatten[Table[Position[nn3, i], {i, Length[nn]}]];
    nn3 = 
      Map[Sign, 
        Partition[Table[Last[nn2[[nn3[[i]]]]], {i, Length[nn3]}], 2]];
    nn2 = 
      If[SameQ[Union[Table[nn3[[i, 1]]*nn3[[i, 2]], {i, Length[nn3]}]], \
{-1}],
         ss, Flatten[
          Table[If[
              SameQ[pp[[i]], 1], {2i, 
                2i - 1, -2i + 1, -2i}, {2i - 1, -2i + 1, -2i, 2i}], {i, 
              Length[pp]}]]];
    nn2 = Table[{nn1[[i]], nn2[[i]]}, {i, Length[nn1]}];
    nn1 = Flatten[Table[Position[nn1, i], {i, 2*Length[nn]}]];
    nn1 = 
      Sort[Abs[Map[Reverse, 
            Map[Sort, 
              Partition[Map[Last, Table[nn2[[nn1[[i]]]], {i, Length[nn1]}]], 
                2]]]]];
    nn=ConstantArray[0,{2*ll,2*ll}];
 (*   nn = ZeroMatrix[2*ll]; *)
    Do[nn = ReplacePart[nn, 1, nn1[[i]]], {i, Length[nn1]}];
    mm = (-1)^vv*Det[mm - nn];
    mm]
    
    
fOddWrithe[Ul_String]:=Module[{gg,gg1,pp,mm,mm1,mm2,i,j},
    gg=fGaussVirtKnot[Ul];
    pp=ToExpression[
        StringJoin["{",
          StringDrop[
            StringReplace[gg,{"O"->"","U"->"","+"->",","-"->","}],-1],"}"]];
    gg1=StringReplace[gg,{"O"->"","U"->""}];
    mm1=Map[First,StringPosition[gg1,"+"]];
    mm2=Map[First,StringPosition[gg1,"-"]];
    mm=Map[Last,
        Union[Table[{mm1[[i]],1},{i,Length[mm1]}],
          Table[{mm2[[j]],-1},{j,Length[mm2]}]]];
    mm=Map[Last,Union[Table[{pp[[i]],mm[[i]]},{i,Length[pp]}]]];
    pp=Map[Flatten,Table[Position[pp,i],{i,Divide[Length[pp],2]}]];
    pp=Table[pp[[i,2]]-pp[[i,1]]-1,{i,Length[pp]}];
    pp=Plus@@Table[If[OddQ[pp[[i]]],mm[[i]],0],{i,Length[mm]}];
    pp
    ] 

fSawollekNorm[Ul_String] := Module[{ss, x, a, y, z},
    ss = fPolyNorm[ReplaceAll[fSawollek[Ul], {x -> a, y -> z}]];
    ss = Factor[ReplaceAll[ss, {a -> x, z -> y}]];
    ss]




ShowVirtKL[Ul_String, s_:7, opts___] := 
  Block[{cc1, cc2, pd1, pd2, mm, mm1, mm2, col},
    cc1 = StringReplace[Ul, "i" -> "1"];
    cc2 = StringReplace[Ul, "i" -> "-1"];
    pd1 = fCreatePData[cc1];
    pd2 = fCreatePData[cc2];
    Install["DrawKnot"];
    mm = GetDrawData[pd1[[2]], pd1[[1]]];
    mm1 = GetDrawData[pd2[[2]], pd2[[1]]];
    Uninstall["DrawKnot"];
    mm2 = mm - mm1;
    cc1 = Position[mm2, 0.8];
    cc2 = Position[mm2, -0.8];
    cc1 = Position[mm2, 0.8];
    cc2 = Position[mm2, -0.8];
    cc1 = Table[Drop[cc1[[i]], -1], {i, Length[cc1]}];
    cc2 = Table[Drop[cc2[[i]], -1], {i, Length[cc2]}];
    cc1 = 
      Flatten[Table[{{cc1[[i, 1]], cc1[[i, 2]]}, {cc2[[j, 1]], 
              cc2[[j, 2]]}}, {i, Length[cc1]}, {j, Length[cc2]}], 1];
    cc2 = 
      Flatten[Table[{{mm[[cc1[[i, 1, 1]], cc1[[i, 1, 2]]]]} - {mm[[
                  cc1[[i, 2, 1]], cc1[[i, 2, 2]]]]}}, {i, Length[cc1]}], 2];
    cc2 = Flatten[Position[cc2, {0., 0., 0.8}]];
    cc1 = Map[Flatten, Table[cc1[[cc2[[i]]]], {i, Length[cc2]}]];
    cc1 = 
      Table[{Red, 
          Sphere[Divide[
              mm[[cc1[[i, 1]], cc1[[i, 2]]]] + mm[[cc1[[i, 3]], cc1[[i, \
4]]]],
               2], 0.7]}, {i, Length[cc1]}];
    col = {Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple};
    Graphics3D[
      Table[Join[cc1, 
          Table[{col[[j]], 
              Cylinder[{mm[[j]][[i]], mm[[j]][[i + 1]]}, 0.2]}, {i, 
              Length[mm[[j]]] - 1}], {{col[[j]], 
              Cylinder[{mm[[j]][[Length[mm[[j]]]]], mm[[j]][[1]]}, 
                0.2]}}], {j, Length[mm]}], Boxed -> False, ImageSize -> 700]]




           
(* ShowKnotfromPdataNew[Ul_, s_:7, opts___] := 
  Block[{cc1, cc2, pd1, pd2, mm, mm1, mm2, col},
    pd1 = Ul;
    pd2 = Ul;
    Install["DrawKnot"];
    mm = GetDrawData[pd1[[2]], pd1[[1]]];
    mm1 = GetDrawData[pd2[[2]], pd2[[1]]];
    Uninstall["DrawKnot"];
    mm2 = mm - mm1;
    cc1 = Position[mm2, 0.8];
    cc2 = Position[mm2, -0.8];
    cc1 = Position[mm2, 0.8];
    cc2 = Position[mm2, -0.8];
    cc1 = Table[Drop[cc1[[i]], -1], {i, Length[cc1]}];
    cc2 = Table[Drop[cc2[[i]], -1], {i, Length[cc2]}];
    cc1 = 
      Flatten[Table[{{cc1[[i, 1]], cc1[[i, 2]]}, {cc2[[j, 1]], 
              cc2[[j, 2]]}}, {i, Length[cc1]}, {j, Length[cc2]}], 1];
    cc2 = 
      Flatten[Table[{{mm[[cc1[[i, 1, 1]], cc1[[i, 1, 2]]]]} - {mm[[
                  cc1[[i, 2, 1]], cc1[[i, 2, 2]]]]}}, {i, Length[cc1]}], 2];
    cc2 = Flatten[Position[cc2, {0., 0., 0.8}]];
    cc1 = Map[Flatten, Table[cc1[[cc2[[i]]]], {i, Length[cc2]}]];
    cc1 = 
      Table[{Red, 
          Sphere[Divide[
              mm[[cc1[[i, 1]], cc1[[i, 2]]]] + mm[[cc1[[i, 3]], cc1[[i, \
4]]]],
               2], 0.7]}, {i, Length[cc1]}];
    col = {Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, 
Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, 
Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, 
Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, 
Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, 
Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple};
    Graphics3D[
      Table[Join[cc1, 
          Table[{col[[j]], 
              Cylinder[{mm[[j]][[i]], mm[[j]][[i + 1]]}, 0.2]}, {i, 
              Length[mm[[j]]] - 1}], {{col[[j]], 
              Cylinder[{mm[[j]][[Length[mm[[j]]]]], mm[[j]][[1]]}, 
                0.2]}}], {j, Length[mm]}], Boxed -> False, ImageSize -> 700]] *)


 ShowKnot3DNew[Ul_]:=Module[{ll,col,ll1,rr1,rr2},
If[Plus@@Ul[[1]]>49,0,Install["DrawKnot"];
ll = GetDrawData[Ul[[2]],Ul[[1]]];
Uninstall["DrawKnot"];
col = {Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, 
Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, 
Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, 
Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, 
Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, 
Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple};
ll=If[SameQ[Length[ll],1],ll={Table[{ll[[1,i,1]],ll[[1,i,2]],5*ll[[1,i,3]]},{i,Length[ll[[1]]]}]},ll];
ll1=Flatten[ll,1];
rr1=ListPlot3D[ll1,Mesh->None,Boxed->False,Axes->False,BoundaryStyle->None,ImageSize->700];
rr2=Table[Graphics3D[{col[[i]],Tube[ll[[i]]]},Boxed->False,ImageSize->700, PlotRange->All],{i,Length[ll]}];
Show[Flatten[{rr1,rr2}],ImageSize->700,PlotRange->All]]]


 ShowKnot3DNewTop[Ul_]:=Module[{ll,col,ll1,rr1,rr2},
If[Plus@@Ul[[1]]>49,0,Install["DrawKnot"];
ll = GetDrawData[Ul[[2]],Ul[[1]]];
Uninstall["DrawKnot"];
col = {Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, 
Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, 
Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, 
Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, 
Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, 
Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple};
ll=If[SameQ[Length[ll],1],ll={Table[{ll[[1,i,1]],ll[[1,i,2]],5*ll[[1,i,3]]},{i,Length[ll[[1]]]}]},ll];
ll1=Flatten[ll,1];
rr1=ListPlot3D[ll1,Mesh->None,Boxed->False,Axes->False,BoundaryStyle->None,ImageSize->700];
rr2=Table[Graphics3D[{col[[i]],Tube[ll[[i]]]},Boxed->False,ImageSize->700, PlotRange->All],{i,Length[ll]}];
Show[Flatten[{rr1,rr2}],ImageSize->700,PlotRange->All,ViewPoint->Top]]]


ShowKnotfromPdataNew[Ul_]:=Module[{ll,col},
 If[Plus@@Ul[[1]]>49,0,Install["DrawKnot"];
    ll = GetDrawData[Ul[[2]],Ul[[1]]];
    Uninstall["DrawKnot"];
col = {Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, 
Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, 
Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, 
Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, 
Brown, Orange, Pink, Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, 
Purple,Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple,Green, Blue, 
Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple};
Show[Table[Graphics3D[{Thickness[0.4],col[[i]],Tube[ll[[i]]]}],{i,Length[ll]}],PlotRange->All,Boxed->False,ImageSize->700]]]
                
                
fGaussExtSignsDirPrKnot[Ul_String] := Module[{ss, ll, ll1, i},
    If[SameQ[StringPosition[Ul, "#"], {}], ss = fGaussExtSigns[Ul], 
      ss = Sort[
          Join[Union[
              Flatten[StringPosition[Ul, "#"]]], {StringLength[Ul]}, {1}]];
      ss = 
        Flatten[{{StringTake[Ul, {1, ss[[2]] - 1}]}, 
            Table[StringTake[Ul, {ss[[i]] + 1, ss[[i + 1]] - 1}], {i, 2, 
                Length[ss] - 2}], {StringTake[
                Ul, {ss[[-2]] + 1, ss[[-1]]}]}}];
      ss = Map[fGaussExtSigns, ss];
      ll = Divide[Map[Length, ss], 2];
      ll1 = Join[{0}, Table[Plus @@ Take[ll, i], {i, Length[ll] - 1}]];
      ss = 
        Flatten[Table[
            Table[ss[[j, i]] + Sign[ss[[j, i]]]*ll1[[j]], {i, 
                Length[ss[[j]]]}], {j, Length[ss]}]]];
    ss]

fGaussVirtKnot[Ul_String] := Module[{cc, cc1, cc2, cc3, ou, i},
    cc = StringReplace[Ul, {"i" -> "1", "-" -> ""}];
    cc = fGaussExtSignsDirPrKnot[cc];
    cc3 = cc;
    cc1 = StringReplace[Ul, "i" -> "1"];
    cc1 = fGaussExtSignsDirPrKnot[cc1];
    cc2 = StringReplace[Ul, "i" -> "-1"];
    cc2 = fGaussExtSignsDirPrKnot[cc2];
    cc = Table[
        If[Not[SameQ[cc1[[i]], cc2[[i]]]], 0, cc1[[i]]], {i, Length[cc1]}];
    cc1 = Table[If[SameQ[cc3[[i]], cc1[[i]]], 1, -1], {i, Length[cc]}];
    cc1 = Table[If[SameQ[cc[[i]], 0], 0, cc1[[i]]], {i, Length[cc]}];
    ou = Table[-(-1)^i*Sign[cc1[[i]]], {i, Length[cc1]}];
    cc1 = Select[cc, Not[SameQ[#, 0]] &];
    ou = Select[ou, Not[SameQ[#, 0]] &];
    cc2 = Abs[cc1];
    cc = Union[Abs[cc1]];
    cc = Table[{cc[[i]], i}, {i, Length[cc]}];
    Do[cc2 = ReplaceAll[cc2, {cc[[i, 1]] -> cc[[i, 2]]}], {i, Length[cc]}];
    cc2 = Table[cc2[[i]]*Sign[cc1[[i]]], {i, Length[cc2]}];
    cc1 = 
      StringJoin[
        Table[StringJoin[If[SameQ[ou[[i]], -1], "U", "O"], 
            ToString[Abs[cc2[[i]]]], 
            If[SameQ[Sign[cc2[[i]]], -1], "-", "+"]], {i, Length[cc2]}]];
    cc1]

fAlexVirt[Ul_] := Module[{aa},
    OpenWrite["input.txt"];
WriteString["input.txt",If[SameQ[Head[Ul],String],fVirtPD[Ul],Ul]];
Close["input.txt"];
    Run[StringJoin["alexander.exe <input.txt>output.txt"]];
    DeleteFile["input.txt"];
    aa = ToExpression[StringReplace[Import["output.txt"], {"\n" -> ""}]];
    DeleteFile["output.txt"];
    aa =Factor[fPolyNorm[aa]];
    aa]


fJonesVirt[Ul_] := Module[{aa},
    OpenWrite["input.txt"];
WriteString["input.txt", If[SameQ[Head[Ul],String],fVirtPD[Ul],Ul]];
Close["input.txt"];
    Run[StringJoin["jones.exe <input.txt>output.txt"]];
    DeleteFile["input.txt"];
    aa = ToExpression[StringReplace[Import["output.txt"], {"\n" -> ""}]];
    DeleteFile["output.txt"];
    aa = fPolyNorm[aa];
    aa]

fCabledJonesVirt[Ul_String, n_Integer] := Module[{aa},
    OpenWrite["input.txt"];
WriteString["input.txt", fGaussVirtKnot[Ul]];
Close["input.txt"];
    Run[StringJoin["cabledjones.exe", " ", ToString[n], " ", 
        "<input.txt>output.txt"]];
    DeleteFile["input.txt"];
    aa = ToExpression[StringReplace[Import["output.txt"], {"\n" -> ""}]];
    DeleteFile["output.txt"];
    aa = fPolyNorm[aa];
    aa]

fGaussVirtPD[Ul_,vv_List ] := Module[{pp,gg,gg1,rr,rr1,i},
pp=fVirtPDfromPD[Ul,vv];
gg=GaussCode[pp];
gg1=Table[gg[[i]],{i,Length[gg]}];
rr=Select[fGenSignfromPD[Ul,vv][[2]],Not[SameQ[#,0]] &];
rr1=Table[If[SameQ[rr[[i]],1],"+","-"],{i,Length[rr]}];
gg=StringJoin[Table[If[SameQ[Sign[gg1[[i]]],1],StringJoin["U",ToString[Abs[gg1[[i]]]],rr1[[Abs[gg1[[i]]]]]],StringJoin["O",ToString[Abs[gg1[[i]]]],rr1[[Abs[gg1[[i]]]]]]],{i,Length[gg1]}]];
gg
]

fCabledJonesVirtPD[Ul_,vv_List ,n_Integer] := Module[{aa,ss},
ss=fGaussVirtPD[Ul,vv];
aa=If[SameQ[ss,""],0,
 OpenWrite["input.txt"];
WriteString["input.txt",ss];
Close["input.txt"];
    Run[StringJoin["cabledjones.exe", " ", ToString[n], " ", 
        "<input.txt>output.txt"]];
    DeleteFile["input.txt"];
    aa = ToExpression[StringReplace[Import["output.txt"], {"\n" -> ""}]];
    DeleteFile["output.txt"];
    aa = fPolyNorm[aa]];
    aa]
    
(* fKauffmanPD[Ul_String]:=Module[{pp,pp1,pp2,pp3,rr,ss,mm,res,i,j},
    pp=fVirtPD[Ul];
    ss=fGenSignfromPDNew[pp];
    mm=Union[Flatten[Table[pp[[i,j]],{i,Length[pp]},{j,4}]]];
    pp=ReplaceAll[Table[Table[pp[[j,i]],{i,4}],{j,Length[pp]}],
        fReplacement[mm]];
    pp=fReplacement1[pp,ss];
pp1=Union[Flatten[Table[Table[pp[[j,i]],{i,Length[pp[[j]]]}],{j,Length[pp]}]]];
pp2=Map[ToString,pp1];
pp3=Flatten[Table[ToCharacterCode[pp2[[i]]],{i,Length[pp2]}]]-96;
rr=Table[pp1[[i]]->pp3[[i]],{i,Length[pp3]}];
res=ReplaceAll[pp,rr];
res=PD@@Table[res[[i]],{i,Length[ss]}];
res]  *)  
    
    
 



 
 (* ## ## ## ## ## ## ## ## ## ## NOVE FUNKCIJE ## ## ## ## ## ## ## ## ## *)
 
 fChangeCrossing[Ul_List,Ul1_List]:=Module[{ss,pp,pp1,i},
    pp=If[SameQ[Ul1,{}],Ul, ss=Length[Ul[[2]]];
    pp=If[Max[Ul1]>2ss||SameQ[Ul1,{}],{},
        pp=Complement[Range[2ss],Abs[Ul[[2]]]];
        pp1=Map[Sign,Ul[[2]]];
        pp=Flatten[Table[{{pp[[i]],Abs[Ul[[2,i]]]}},{i,ss}],1];
        pp1=Table[If[MemberQ[Ul1,i],-pp1[[i]],pp1[[i]]],{i,Length[pp1]}];
        pp=Table[If[MemberQ[Ul1,i],Reverse[pp[[i]]],pp[[i]]],{i,Length[pp]}];
        pp=Sort[Table[{pp[[i,1]],pp[[i,2]]*pp1[[i]]},{i,Length[pp]}]];
        pp={Ul[[1]],Map[Last,pp]};
        pp]];
    pp]   


(* ## ## ## ## ## ## ## #   KAMADA CODE AND MIYAZAWA POLYNOMIAL  ## ## ## ##  \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)

UnsortedComplement[x_List, y__List] := 
    Replace[x, Dispatch[(# :> Sequence[]) & /@ Union[y]], 1];

fKamadaCode[Ul_String] := Module[{cc, cc1, cc2, cc3, ou, i},
    cc = StringReplace[Ul, {"i" -> "1", "-" -> ""}];
    cc = fGaussExtSignsDirPrKnot[cc];
    cc3 = cc;
    cc1 = StringReplace[Ul, "i" -> "1"];
    cc1 = fGaussExtSignsDirPrKnot[cc1];
    cc2 = StringReplace[Ul, "i" -> "-1"];
    cc2 = fGaussExtSignsDirPrKnot[cc2];
    cc = Table[
        If[Not[SameQ[cc1[[i]], cc2[[i]]]], 0, cc1[[i]]], {i, Length[cc1]}];
    cc1 = Table[If[SameQ[cc3[[i]], cc1[[i]]], 1, -1], {i, Length[cc]}];
    cc1 = Table[If[SameQ[cc[[i]], 0], 0, cc1[[i]]], {i, Length[cc]}];
    ou = Table[-(-1)^i*Sign[cc1[[i]]], {i, Length[cc1]}];
    cc1 = Select[cc, Not[SameQ[#, 0]] &];
    ou = Select[ou, Not[SameQ[#, 0]] &];
    cc2 = Abs[cc1];
    cc = Union[Abs[cc1]];
    cc = Table[{cc[[i]], i}, {i, Length[cc]}];
    Do[cc2 = ReplaceAll[cc2, {cc[[i, 1]] -> cc[[i, 2]]}], {i, Length[cc]}];
    cc2 = Table[cc2[[i]]*Sign[cc1[[i]]], {i, Length[cc2]}];
    cc = Union[cc2];
    cc1 = Flatten[Table[Position[cc2, cc[[i]]] - 1, {i, Length[cc]}]];
    cc = Map[Sign, cc];
    ou = Map[First, 
        Partition[Table[ou[[cc1[[i]] + 1]], {i, Length[ou]}], 2]];
    cc1 = Partition[cc1, 2];
    cc1 = Map[Sort, cc1];
    cc1 = 
      Sort[Table[
          If[SameQ[ou[[i]], 1], {cc1[[i]], cc[[i]]}, {Reverse[cc1[[i]]], 
              cc[[i]]}], {i, Length[cc1]}]];
    cc1 = {Map[First, cc1], Map[Last, cc1]};
    cc1]


 fMiyazawaPoly[Ul_String] := 
  Module[{cs, diagram, sgn, n, writhe, crsnop, cnncrsnop, crsno, cnncrsno, 
      vnolst, vcrs, splice, statep, state, nocomp, ccomptmp, wt, statetmp, 
      vcrstmp, vsgn, wtmap, splicetmp, bracket1, bracket2, Mzply1, Mzply2, i, \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

      j, k},
    cs = fKamadaCode[Ul];
    diagram = cs[[1]];
    sgn = cs[[2]];
    n = Length[diagram];
    writhe = Sum[sgn[[i]], {i, 1, Length[sgn]}]; 
    crsnop = 
      Table[If[sgn[[i]] == 1, diagram[[i]], Reverse[diagram[[i]]]], {i, 1, 
          n}];
    cnncrsnop = Table[Reverse[Mod[crsnop[[i]] - 1, 2*n]], {i, 1, n}];
    crsno = Flatten[crsnop];
    cnncrsno = Flatten[cnncrsnop];
    vnolst = Table[
        Select[crsno, 
          Position[crsno, #][[1]][[1]] > 
                Position[crsno, crsno[[i]]][[1]][[1]] && 
              Position[cnncrsno, #][[1]][[1]] < 
                Position[cnncrsno, crsno[[i]]][[1]][[1]] &], {i, 1, 2*n}];
    vcrs = 
      Union[Flatten[
          Table[If[vnolst[[i]] != {}, 
              Table[{vnolst[[i]][[j]], crsno[[i]]}, {j, 1, 
                  Length[vnolst[[i]]]}], {}], {i, 1, 2*n}], 1]];
    splice = Table[IntegerDigits[i - 1, 2, n], {i, 1, 2^n}]; 
statep = Table[Table[
          If[splice[[i]][[j]] == 0,
            {{crsnop[[j]][[1]], cnncrsnop[[j]][[1]]}, {crsnop[[j]][[2]], 
                cnncrsnop[[j]][[2]]}},
            {{crsnop[[j]][[1]], crsnop[[j]][[2]]}, {cnncrsnop[[j]][[1]], 
                cnncrsnop[[j]][[2]]}}]
          , {j, 1, n}], {i, 1, 2^n}];
state = Table[Flatten[statep[[i]], 1], {i, 1, 2^n}];
nocomp = Table[0, {i, 1, 2^n}];
Do[
      Do[
        ccomptmp = Select[state[[i]], MemberQ[#, j] &];
        state = 
          ReplacePart[state, 
            Prepend[UnsortedComplement[state[[i]], ccomptmp], 
              Union[Flatten[ccomptmp]]], i];
        , {j, 0, 2*n - 1}];
      nocomp = ReplacePart[nocomp, Length[state[[i]]], i];
      , {i, 1, 2^n}];
    wt = Table[(-1)^IntegerDigits[i, 2, 2*n], {i, 0, 2^(2*n) - 1}];
statetmp = statep + 1;
vcrstmp = vcrs + 1;
    vsgn = Table[0, {i, 1, 2^n}, {j, 1, Length[vcrstmp]}];
Do[
      wtmap = wt;
      Do[
        If[splice[[i]][[j]] == 0,
          wtmap = Select[wtmap,
              #[[statetmp[[i]][[j]][[1]][[1]]]] == \
#[[statetmp[[i]][[j]][[1]][[2]]]] &&
                  #[[statetmp[[i]][[j]][[2]][[1]]]] == \
#[[statetmp[[i]][[j]][[2]][[2]]]] &],
          wtmap = Select[wtmap,
              #[[statetmp[[i]][[j]][[1]][[1]]]] != \
#[[statetmp[[i]][[j]][[1]][[2]]]] &&
                  #[[statetmp[[i]][[j]][[2]][[1]]]] != \
#[[statetmp[[i]][[j]][[2]][[2]]]] &]]
        , {j, 1, Length[splice[[i]]]}];
      vsgn[[i]] = Table[
          
          If[wtmap[[k]][[vcrstmp[[j]][[1]]]] == 1, 
            If[wtmap[[k]][[vcrstmp[[j]][[2]]]] == -1, 1, 0],
            If[wtmap[[k]][[vcrstmp[[j]][[2]]]] == 1, -1, 0]], {k, 1, 
            Length[wtmap]}, {j, 1, Length[vcrstmp]}];, {i, 1, 2^n}];
splicetmp = Table[(-1)^splice[[i]], {i, 1, 2^n}];
bracket1 = Expand[
        Sum[A^Sum[splicetmp[[i]][[j]]*sgn[[j]], {j, 1, n}]
            *(-A^2 - A^(-2))^(nocomp[[i]] - 1)*2^(-nocomp[[i]])*
            
            Sum[t^Sum[vsgn[[i]][[j]][[k]], {k, 1, Length[vcrstmp]}], {j, 1, 
                Length[vsgn[[i]]]}],
          {i, 1, 2^n}]];
bracket2 = Expand[
        Sum[A^(Sum[splicetmp[[i]][[j]]*sgn[[j]], {j, 1, n}])
            *(-A^2 - A^(-2))^(nocomp[[i]] - 1)*2^(-nocomp[[i]])*
            
            Sum[Abs[Sum[vsgn[[i]][[j]][[k]], {k, 1, Length[vcrstmp]}]], {j, \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
1,
                 Length[vsgn[[i]]]}],
          {i, 1, 2^n}]];
Mzply1 = (-A)^(-3*writhe)*bracket1;
Mzply2 = (-A)^(-3*writhe)*bracket2;
    Mzply1 = Expand[Mzply1];
    Mzply2 = Expand[Mzply2];
    cs = {Mzply1, Mzply2};
    cs
    ]


fBothPoly[pp_]:=Module[{vvv,dd0,ppp,pp1,vv1,vv,pp0,vv2,vv3,i},
dd0=pp[[1]];
pp1=Prepend[Table[If[Not[NumberQ[pp[[i]]]],pp[[i]],pp[[i,1]]],
{i,2,Length[pp]}],dd0];
vv=Reverse[pp1];
vv1=Table[vv1=vv[[i]];
pp0=Position[vv1,A];
ppp=If[SameQ[pp0,{}],1,First[Flatten[pp0]]];
If[SameQ[Head[vv[[i]]],Integer],vv[[i]],Take[vv1,{ppp,Length[vv1]}]],{i,Length[vv]}];
vv2=Table[vvv=vv[[i]];
pp0=Position[vvv,A];
ppp=If[SameQ[pp0,{}],1,First[Flatten[pp0]]];
If[SameQ[Head[vv[[i]]],Integer],vv[[i]],Take[vvv,{1,ppp-1}]],{i,Length[vv]}];
vv2=If[SameQ[Sign[vv2[[1]]],1],vv2,-vv2];
vv3=vv2*Join[Table[Divide[vv1[[1]],vv1[[i]]],{i,Length[vv1]-1}],{vv1[[1]]}];
vv=Plus@@vv3;
vv=Sort[{pp,vv}];
vv]
    
    (* ## ## ## ## ## ## ## ## # LABELED IMMERSION CODE AND VIRTUAL BRAIDS \
#### ## ## ## ## ## # *)
    
     fGaussExtSignsDirPrKnot[Ul_String]:=Module[{ss,ll,ll1,i},
    If[SameQ[StringPosition[Ul,"#"],{}],ss=fGaussExtSigns[Ul],
      ss=Sort[Join[
            Union[Flatten[StringPosition[Ul,"#"]]],{StringLength[Ul]},{1}]];
      ss=Flatten[{{StringTake[Ul,{1,ss[[2]]-1}]},
            Table[StringTake[Ul,{ss[[i]]+1,ss[[i+1]]-1}],{i,2,
                Length[ss]-2}],{StringTake[Ul,{ss[[-2]]+1,ss[[-1]]}]}}];
      ss=Map[fGaussExtSigns,ss];
      ll=Divide[Map[Length,ss],2];
      ll1=Join[{0},Table[Plus@@Take[ll,i],{i,Length[ll]-1}]];
      ss=Flatten[
          Table[Table[
              ss[[j,i]]+Sign[ss[[j,i]]]*ll1[[j]],{i,Length[ss[[j]]]}],{j,
              Length[ss]}]]];
    ss]
    
    
    fFindSignsPart[Ul_String]:=Module[{ss,ss1,i},
    ss=Flatten[fGaussExtSigns[Ul]];
    ss1=Abs[ss];
    ss=Map[Sign,Map[Last,Union[Table[{ss1[[i]],ss[[i]]},{i,Length[ss]}]]]];
    ss]
    
    
fFindSigns[Ul_String]:=Module[{ss,i},
    If[SameQ[StringPosition[Ul,"#"],{}],ss=fFindSignsPart[Ul],
      ss=Sort[Join[
            Union[Flatten[StringPosition[Ul,"#"]]],{StringLength[Ul]},{1}]];
      ss=Flatten[{{StringTake[Ul,{1,ss[[2]]-1}]},
            Table[StringTake[Ul,{ss[[i]]+1,ss[[i+1]]-1}],{i,2,
                Length[ss]-2}],{StringTake[Ul,{ss[[-2]]+1,ss[[-1]]}]}}];
      ss=Flatten[Map[fFindSignsPart,ss]]];
    ss]
    
    
fMakeStrPart[Ul_List]:=Module[{tt,i},
    tt=StringJoin[
        Table[StringJoin[" ",If[SameQ[Ul[[i,1]],-1],"-","+"],
            ToString[Ul[[i,2]]]],{i,Length[Ul]}]];
    tt=StringJoin["(",StringDrop[tt,1],")"];
    tt]
    
    
fLabeledImmersionCode[Ul_String]:=
  Module[{cc,cc1,cc2,cc3,ppp,pp,pp1,pp2,ss,qq,pos,i},
    cc1=StringReplace[Ul,"i"->"1"];
    cc2=StringReplace[Ul,"i"->"-1"];
    cc1=fConwayToPD[cc1];
    cc2=fConwayToPD[cc2]; 
    pp1=Table[cc1[[i,j]],{i,Length[cc1]},{j,4}];
    pp2=Table[cc2[[i,j]],{i,Length[cc1]},{j,4}];
    cc=pp1-pp2;
    cc=Flatten[Position[cc,{0,0,0,0}]];
    cc=Table[If[MemberQ[cc,i],1,0],{i,Length[cc1]}];
    pp=fFindSigns[StringReplace[Ul,{"i"->"1","-1"->"1"}]]; 
    ppp=fFindSigns[StringReplace[Ul,{"i"->"1"}]]; 
    qq=ppp*pp*cc;
    ss=pp1-1;
    cc=Abs[fGaussExtSignsDirPrKnot[StringReplace[Ul,"i"->"1"]]]-1;
    cc2=Range[Length[cc]]-1;
    cc=Sort[Table[{cc[[i]],cc2[[i]]},{i,Length[cc]}]];
    cc=Partition[Map[Last,cc],2];
    cc=Table[If[EvenQ[cc[[i,1]]],cc[[i]],Reverse[cc[[i]]]],{i,Length[cc]}];
    cc1=Sort[cc];
    pos=Flatten[Table[Position[cc,cc1[[i]]],{i,Length[cc]}]];
    cc=Mod[Divide[Map[Last,cc1]+1,2],Length[cc1]];
    pp=-Table[pp[[pos[[i]]]],{i,Length[pp]}];
    qq=-Table[qq[[pos[[i]]]],{i,Length[qq]}];
    cc1=ToCycles[cc+1];
    pos=Table[Position[cc1,cc[[i]]+1],{i,Length[cc]}];
    Do[cc1=ReplacePart[cc1,cc[[i]],pos[[i]]],{i,Length[pp]}];
    pos=Map[Length,cc1];
    cc2=Flatten[cc1];
    pp1=Table[If[SameQ[pp[[i]],1],"","-"],{i,Length[pp]}];
    pp2=Map[ToString,cc2];
    cc1=Table[StringJoin[pp1[[i]],pp2[[i]]," "],{i,Length[pp1]}];
    cc1=ToString[iteratedTake[cc1,pos]];
    cc1=StringReplace[cc1,{","->"","{{"->"("," }}"->")"}];
    cc1=StringReplace[cc1,{" }"->")","{"->"("}];
    qq=StringJoin[
        Table[Switch[qq[[i]],1," +",-1," -",0," *"],{i,Length[qq]}]];
    cc=StringJoin[cc1," /",qq];
    cc]


fBraidVirt[Ul_String]:=Module[{ddd},
    OpenWrite["C://LinKnot//braidvirtual//input.txt"];
    WriteString["C://LinKnot//braidvirtual//input.txt",
fLabeledImmersionCode[
Ul]];
    Close["C://LinKnot//braidvirtual//input.txt"];
    Run["C:/LinKnot/braidvirtual/braid -v \
C://LinKnot//braidvirtual//input.txt \
C://LinKnot//braidvirtual//output.txt"];
    DeleteFile["C://LinKnot//braidvirtual//input.txt"];
    ddd=Import["C://LinKnot//braidvirtual//output.txt"];
    DeleteFile["C://LinKnot//braidvirtual//output.txt"]; 
    ddd=StringTake[
        ddd,{Last[Flatten[StringPosition[ddd,"= "]]]+1,
          First[Flatten[StringPosition[ddd,"\n\n\n"]]]-1}] ; 
ddd=ToExpression[StringJoin["{",StringDrop[StringReplace[StringReplace[ddd,{"s"->",","t"->",100"}],{"-,"->",-"}],1],"}"]];
ddd=StringJoin[Map[ToString,Table[If[ddd[[i]]>1000,ddd[[i]]-1000,If[ddd[[i]]>0,FromCharacterCode[ddd[[i]]+96],FromCharacterCode[Abs[ddd[[i]]]+64]]],{i,Length[ddd]}]]];
    ddd]
    
    
    fSawollekBraid[Ul_String]:=Module[{ddd},
    OpenWrite["C://LinKnot//braidvirtual//input1.txt"];
    WriteString["C://LinKnot//braidvirtual//input1.txt",fBraidVirt[Ul]];
    Close["C://LinKnot//braidvirtual//input1.txt"];
    Run["C:/LinKnot/braidvirtual/braid -s \
C://LinKnot//braidvirtual//input1.txt \
C://LinKnot//braidvirtual//output1.txt"];
    DeleteFile["C://LinKnot//braidvirtual//input1.txt"];
    ddd=Import["C://LinKnot//braidvirtual//output1.txt"];
    DeleteFile["C://LinKnot//braidvirtual//output1.txt"]; 
    ddd=StringTake[
        ddd,{Last[Flatten[StringPosition[ddd,"= "]]]+1,
          First[Flatten[StringPosition[ddd,"\n\n\n"]]]-1}] ; 
    ddd]
    
    
    (* ## ## ## ## ## ## # VIRTUAL UNKNOTTING ## ## ## ## ## ## ## ## ## ## # \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
*)
    
    
  fRealCrossingChange[Ul_String]:=Module[{pp,pp1,pp2,ss,ss1,ss2,i},
    pp=StringPosition[Ul,"1"];
    pp1=StringPosition[Ul,"-"];
    pp2=StringPosition[Ul,"-1"];
    ss=Union[Flatten[pp]];
    ss1=Union[Flatten[pp1]];
    ss2=Union[Flatten[pp2]];
    ss=Complement[Union[ss,ss1],ss2];
    ss=Table[{ss[[i]],ss[[i]]},{i,Length[ss]}];
    ss=Table[StringReplacePart[Ul,"-1",ss[[i]]],{i,Length[ss]}];
    ss1=Table[StringReplacePart[Ul,"1",pp2[[i]]],{i,Length[pp2]}];
    ss=Union[ss,ss1];
    ss1=Map[fBothPoly,Table[fCabledJonesVirt[ss[[i]],3],{i,Length[ss]}]];
    ss1 = 
      If[MemberQ[ss1, {1+2A^4+A^8,1+2A^4+A^8}], "1", 
          ss1 = Table[Position[virttab, ss1[[i]]], {i,Length[ss1]}]; 
          ss1 = Union[Map[First, Flatten[ss1, 1]]]; 
          ss1 = Map[Last, Table[VirtKnotTab[ss1[[i]]], {i, Length[ss1]}]]];
    ss1] 
    
fRealUnknottingNo[Ul_String]:=Module[{ss1,ss,k},
    ss={Ul};
    ss=Union[Flatten[Map[fRealCrossingChange,ss]]];
    ss1=ss;
    k=If[SameQ[First[ss],"1"],1,2];
    ss=Union[Flatten[Map[fRealCrossingChange,ss1]]];
    ss1=Union[Join[ss1,ss]];
    While[Not[SameQ[ss,ss1]]&&Not[SameQ[ss[[1]],"1"]],
      ss=Union[Flatten[Map[fRealCrossingChange,ss1]]];
      ss1=Union[Join[ss1,ss]]; k++];
    k=If[SameQ[ss[[1]],"1"],k,0];
    k] 
    
fVirtualCrossingChange[Ul_String]:=Module[{pp,pp1,pp2,ss,ss1,ss2,i},
    pp=StringPosition[Ul,"1"];
    pp1=StringPosition[Ul,"-"];
    pp2=StringPosition[Ul,"-1"];
    ss=Union[Flatten[pp]];
    ss1=Union[Flatten[pp1]];
    ss2=Union[Flatten[pp2]];
    ss=Complement[Union[ss,ss1],ss2];
    ss=Table[{ss[[i]],ss[[i]]},{i,Length[ss]}];
    ss=Table[StringReplacePart[Ul,"i",ss[[i]]],{i,Length[ss]}];
    ss1=Table[StringReplacePart[Ul,"i",pp2[[i]]],{i,Length[pp2]}];
    ss=Union[ss,ss1];
    ss1=Map[fBothPoly,Table[fCabledJonesVirt[ss[[i]],3],{i,Length[ss]}]];
     ss1 = 
      If[MemberQ[ss1, {1+2A^4+A^8,1+2A^4+A^8}], "1", 
          ss1 = Table[Position[virttab, ss1[[i]]], {i,Length[ss1]}]; 
          ss1 = Union[Map[First, Flatten[ss1, 1]]]; 
          ss1 = Map[Last, Table[VirtKnotTab[ss1[[i]]], {i, Length[ss1]}]]];
    ss1] 
    
fVirtualUnknottingNo[Ul_String]:=Module[{ss,k},
    ss={Ul};
    ss=Union[Flatten[Map[fVirtualCrossingChange,ss]]];
    k=1;
    While[Not[SameQ[First[ss],"1"]],
      ss=Union[Flatten[Map[fVirtualCrossingChange,ss]]];k++];
    k] 
    
    
    fMixedCrossingChange[Ul_String]:=Module[{ss3,pp,pp1,pp2,ss,ss1,ss2,
    tt,tt1,tt2,i},
    pp=StringPosition[Ul,"1"];
    pp1=StringPosition[Ul,"-"];
    pp2=StringPosition[Ul,"-1"];
    ss=Union[Flatten[pp]];
    ss1=Union[Flatten[pp1]];
    ss2=Union[Flatten[pp2]];
    ss=Complement[Union[ss,ss1],ss2];
    ss=Table[{ss[[i]],ss[[i]]},{i,Length[ss]}];
    ss=Table[StringReplacePart[Ul,"-1",ss[[i]]],{i,Length[ss]}];
    ss1=Table[StringReplacePart[Ul,"1",pp2[[i]]],{i,Length[pp2]}];
    tt=Union[Flatten[pp]];
    tt1=Union[Flatten[pp1]];
    tt2=Union[Flatten[pp2]];
    tt=Complement[Union[tt,tt1],tt2];
    tt=Table[{tt[[i]],tt[[i]]},{i,Length[tt]}];
    ss2=Table[StringReplacePart[Ul,"i",tt[[i]]],{i,Length[tt]}];
    ss3=Table[StringReplacePart[Ul,"i",pp2[[i]]],{i,Length[pp2]}];
    ss=Union[ss,ss1,ss2,ss3];
    ss1=Map[fBothPoly,Table[fCabledJonesVirt[ss[[i]],3],{i,Length[ss]}]];
     ss1 = 
      If[MemberQ[ss1, {1+2A^4+A^8,1+2A^4+A^8}], "1", 
          ss1 = Table[Position[virttab, ss1[[i]]], {i,Length[ss1]}]; 
          ss1 = Union[Map[First, Flatten[ss1, 1]]]; 
          ss1 = Map[Last, Table[VirtKnotTab[ss1[[i]]], {i, Length[ss1]}]]];
    ss1] 
    
    
  fMixedUnknottingNo[Ul_String]:=Module[{ss,k},
    ss={Ul};
    ss=Union[Flatten[Map[fMixedCrossingChange,ss]]];
    k=1;
    While[Not[SameQ[First[ss],"1"]],
      ss=Union[Flatten[Map[fMixedCrossingChange,ss]]];k++];
    k] 

fVirtCrossChangesFix[Ul_String]:=Module[{pp,pp1,pp2,ss,ss1,ss2,ss3,i,k},
pp=StringPosition[Ul,"1"];
pp1=StringPosition[Ul,"-"];
pp2=StringPosition[Ul,"-1"];
ss=Union[Flatten[pp]];
ss1=Union[Flatten[pp1]];
ss2=Union[Flatten[pp2]];
ss=Complement[Union[ss,ss1],ss2];
ss=Table[{ss[[i]],ss[[i]]},{i,Length[ss]}];
ss=Union[pp2,ss];
k=1;
ss3=fBothPoly[fCabledJonesVirt[Ul,3]];
While[Not[SameQ[ss3[[1]],{1+2 A^4+A^8,1+2 A^4+A^8}]],
ss1=KSubsets[ss,k];
ss2=Table[StringReplacePart[Ul,"i",ss1[[i]]],{i,Length[ss1]}];
ss3=Union[Map[fBothPoly,Table[fCabledJonesVirt[ss2[[i]],3],{i,Length[ss2]}]]];\
k++];
k-1] 

fRealCrossChangesFix[Ul_String]:=Module[{pp,pp1,pp2,ss,
ss1,ss2,ss3,ss4,ss5,ll,i,k},
pp=StringPosition[Ul,"1"];
pp1=StringPosition[Ul,"-"];
pp2=StringPosition[Ul,"-1"];
ss=Union[Flatten[pp]];
ss1=Union[Flatten[pp1]];
ss2=Union[Flatten[pp2]];
ss=Complement[Union[ss,ss1],ss2];
ss=Table[{ss[[i]],ss[[i]]},{i,Length[ss]}];
ss=Union[pp2,ss];
ll=Length[ss];
k=1;
ss5=fBothPoly[fCabledJonesVirt[Ul,3]];
While[Not[SameQ[ss5[[1]],{1+2 A^4+A^8,1+2 A^4+A^8}]]&&k<=ll,
ss1=KSubsets[ss,k];
ss4=fBothPoly[fCabledJonesVirt[Ul,3]];
ss3=Table[ss2=StringReplacePart[Ul,"22",Select[ss1[[i]],SameQ[Length[Union[#]]\
,2] &]];
ss3=StringReplacePart[ss2,"33",Select[ss1[[i]],SameQ[Length[Union[#]],1] &]];
ss3=StringReplace[ss3,{"33"->"-1","22"->"1"}],{i,Length[ss1]}];
ss5=Union[Map[fBothPoly,Table[fCabledJonesVirt[ss3[[i]],3],
{i,Length[ss3]}]]];
 k++];
ss=If[k<=ll,k-1,0];
ss] 

fUnknottingVirt[Ul_String]:=Module[{rr1,rr2,rr3,rr4,res},
rr1=fVirtCrossChangesFix[Ul];
rr2=fVirtualUnknottingNo[Ul];
rr3=fRealCrossChangesFix[Ul];
rr4=fRealUnknottingNo[Ul];
res={Ul,rr1,rr2,rr3,rr4}] 

fUnknottingVirt[Ul_String]:=Module[{rr1,rr2,rr3,rr4,res},
rr1=fVirtCrossChangesFix[Ul];
rr2=fVirtualUnknottingNo[Ul];
rr3=fRealCrossChangesFix[Ul];
rr4=fRealUnknottingNo[Ul];
res={Ul,rr1,rr2,rr3,rr4}] 
    
   
    (* ## ## ## ## ## ## VIRTUAL BRAID DRAWING ## ## ## ## ## ## # *)

 framing2[cord_List,i_,t_]:=Block[{v1,v2,v3,dis,u,v,tt,k,eta,ett,lcd,ss,co,si,\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
r,p},eta=t 2 Pi;
r=0.1;
co=Cos[eta];
si=Sin[eta];
If[i==0,v1=cord[[i+1]],v1=cord[[i]]];
lcd=Length[cord];
If[lcd==2||i==lcd||i==0,p={v1[[1]],v1[[2]]+r*co,v1[[3]]+r*si},v2=cord[[i+1]];
v3={v2[[1]]-v1[[1]],v2[[2]]-v1[[2]],v2[[3]]-v1[[3]]};
u={-v3[[2]],v3[[1]],0.0};
dis=Sqrt[u[[1]]^2+u[[2]]^2];
tt={r*u[[1]]/dis,r*u[[2]]/dis,r*u[[3]]/dis};
v=Cross[v3,tt];
dis=Sqrt[v[[1]]^2+v[[2]]^2+v[[3]]^2];
k={r*v[[1]]/dis,r*v[[2]]/dis,r*v[[3]]/dis};
p={v2[[1]]+co*tt[[1]]+si*k[[1]],v2[[2]]+co*tt[[2]]+si*k[[2]],v2[[3]]+co*tt[[3]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]+si*k[[3]]};];
{p[[1]],p[[2]],p[[3]]}]
basicCord[pp_,up_,rd_]:=Block[{ppp,pp2,i,cod,eta},pp2=2*pp+1;
cod=Table[0,{i,pp2}];
ppp=Pi/(2*pp);
For[i=1,i<=pp,i++,If[i==1,j=0.0,j=up*0.1*(i-2)];
If[rd==1,eta=Pi/2-ppp*(i-1),eta=-Pi/2+ppp*(i-1)];
cod[[i]]={0.5*Cos[eta],0.5*Sin[eta],j};
cod[[pp2+1-i]]={1.0-cod[[i,1]],0.0-cod[[i,2]],cod[[i,3]]};];
cod[[pp+1]]={0.5,0.0,cod[[pp,3]]};
cod]



expandpower[bnum_]:=Block[{i,bnumx=bnum,fp,len,tlen,ch},tlen=Length[bnum];
While[MemberQ[bnumx,32],bnumx=Delete[bnumx,First[First[Position[bnumx,32]]]];\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
tlen--];
While[MemberQ[bnumx,94],len=0;
fp=First[First[Position[bnumx,94]]];
bnumx=Delete[bnumx,fp];tlen--;
While[fp<=tlen&&48<=bnumx[[fp]]&&bnumx[[fp]]<=57,len*=10;
len+=bnumx[[fp]]-48;
bnumx=Delete[bnumx,fp];tlen--];
ch=bnumx[[fp-1]];
For[i=1,i<len,i++,bnumx=Insert[bnumx,ch,fp];tlen++]];
Return[bnumx]]
bwdtonumwd[bwd_]:=Block[{bnum,cds,len,i},bnum=ToCharacterCode[bwd];
While[MemberQ[bnum,45],bnum=Delete[bnum,First[First[Position[bnum,45]]]];
Print["The character - is illegal and so delete it"]];
bnum=expandpower[bnum];(*Modified by M.Ochiai*)len=Length[bnum];
cds=Table[0,{i,len}];
For[i=1,i<=len,i++,If[bnum[[i]]>=97,cds[[i]]=bnum[[i]]-96,cds[[i]]=64-bnum[[i]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]];];
cds]



ShowBraidNew[Ul_,n_:4,s_:5]:=Block[{bnum,c,i,j,h,bi,cds,cd,cc,cap1,cap2,
low1,low2,cds2,cyl,st,mm,ll,ll1,ll2,bwd,cd1,cd2,sg,k},
bwd=Ul;
ll=Flatten[Position[Table[If[LetterQ[StringTake[Ul,{i,i}]],i,0],{i,\
StringLength[Ul]}],0]];
ll1=Table[StringTake[Ul,{ll[[i]],ll[[i]]}],{i,Length[ll]}];
ll2=Map[FromCharacterCode,ToExpression[Table[StringTake[Ul,{ll[[i]],ll[[i]]}],\

{i,Length[ll]}]]+96];
Do[bwd=StringReplacePart[bwd,ll2[[i]],{ll[[i]],ll[[i]]}],{i,Length[ll]}];
bnum=bwdtonumwd[bwd];
c=Length[bnum];
st=Abs[bnum[[1]]];
For[i=2,i<=c,i++,If[Abs[bnum[[i]]]>=st,st=Abs[bnum[[i]]]]];
cap1=framed[basicCord[n,-1,1],s];
cap2=framed[basicCord[n,1,-1],s];
low1=framed[basicCord[n,1,1],s];
low2=framed[basicCord[n,-1,-1],s];
cyl=framed[{{0.0,0.5,0.0},{1.0,0.5,0.0}},s];
cyl={cyl[[2]],cyl[[3]]};
cds=Table[0,{i,c},{j,st+1}];
cd=Table[0,{j,s}];
cd1=Table[0,{i,s}];
cd2=Table[0,{i,s}];
cds2=Table[0,{i,n*2+1}];
For[i=1,i<=c,i++,bi=bnum[[i]];
If[bi<0,sg=1;bi=-bi,sg=-1];
For[j=1,j<=st+1,j++,If[j<bi||bi+1<j,For[h=1,h<=s,h++,cd1[[h]]={i-1,1.5-j+
cyl[[1,h,2]],cyl[[1,h,3]]};
cd2[[h]]={i,1.5-j+cyl[[2,h,2]],cyl[[2,h,3]]};];
cds[[i,j]]={cd1,cd2},For[h=1,h<=n*2+1,h++,If[j==bi,If[sg==1,
For[k=1,k<=s,k++,
cd[[k]]={i-1+cap1[[h,k,1]],1.5-j+cap1[[h,k,2]],cap1[[h,k,3]]}];
cds2[[h]]=cd,For[k=1,k<=s,k++,cd[[k]]={i-1+low1[[h,k,1]],1.5-j+
low1[[h,k,2]],
low1[[h,k,3]]}];
cds2[[h]]=cd],If[sg==1,For[k=1,k<=s,k++,cd[[k]]={i-1+cap2[[h,k,1]],1.5-j+1+\
cap2[[h,k,2]],cap2[[h,k,3]]}];
cds2[[h]]=cd,For[k=1,k<=s,k++,cd[[k]]={i-1+low2[[h,k,1]],1.5-j+1+
low2[[h,k,2]],low2[[h,k,3]]}];
cds2[[h]]=cd];];];
cds[[i,j]]=cds2]]];
cds;
mm=Graphics3D[{Table[Table[mm=cds[[k,j]];
mm=Table[Plus@@Table[mm[[j,i]],{i,Length[mm[[j]]]}],{j,Length[mm]}];
mm=Table[mm[[i]]/5,{i,Length[mm]}];
Table[{Cylinder[{mm[[i]],mm[[i+1]]},0.1]},{i,Length[mm]-1}],
{j,Length[cds[[1]]
]}],{k,Length[cds]}],
Red,Table[Sphere[Plus@@Divide[Plus@Plus@@Flatten[Select[
cds[[ll[[i]]]],SameQ[Length[#],9] 
&],1],90],0.35],{i,Length[ll]}]},Boxed->False,ImageSize->700];
mm
]




framed[cord_List,n_]:=Block[{i,h,m,mm,pp},pp=Table[0,{i,Length[cord]+1},{h,n}]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
;
m=0.0;mm=1.0/n;
For[i=1,i<=Length[cord]+1,i++,For[h=1,h<=n,h++,pp[[i,h]]=Chop[framing2[cord,i-\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
1,m]];
m=m+mm;];];
pp]


(* ## ## ## ## ## ## ## ## fMakeConVirt  ## ## ## ## ## # *) 

fBasVirt[n_Integer]:=Module[{ss,ss1},
ss=fVarP3[n];
ss=Reverse[Sort[Select[ss,Not[SameQ[Count[#,1],n]] &]]];
ss1=Flatten[Table[Position[ss,ReplaceAll[ss[[i]],{2->0,0->2}]],{i,Length[ss]}]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
];
ss1=Union[Table[If[SameQ[i,ss1[[i]]],i,Min[i,ss1[[i]]]],{i,Length[ss1]}]];
ss=Table[ss[[ss1[[i]]]],{i,Length[ss1]}];
ss1=Select[ss,Count[#,1]<=IntegerPart[Divide[Length[#],2]]+1&& \
Count[#,0]<=IntegerPart[Divide[Length[#],2]]+1 &];
ss] 

fExpandCon[Ulaz_String]:=Module[{Ul,pp,pp1,pp2,pp3,pp4,ss,i,j},
If[SameQ[StringTake[Ulaz,-1],"*"],Ul=StringJoin[Ulaz,"1"],Ul=Ulaz];
pp=Flatten[StringPosition[Ul,"*"]];
pp2=If[SameQ[pp,{}],pp2=Ul,
pp1=StringTake[Ul,First[Union[pp]]];
pp2=StringTake[Ul,{First[Union[pp]]+1,StringLength[Ul]}]];
pp3=Table[StringTake[pp2,{i,i}],{i,StringLength[pp2]}];
pp3=Complement[Select[pp3,DigitQ[#] &],{"0"}];
pp4=Map[ToExpression,pp3];
pp4=StringReplace[StringReplace[Map[ToString,Table[Table[1,{i,pp4[[j]]}],{j,\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

Length[pp4]}]],{"{"->"(","}"->")"}],{" "->""}];
Do[pp2=StringReplace[pp2,{ToString[pp3[[i]]]->pp4[[i]]}],{i,Length[pp3]}];
ss=StringPosition[pp2,")("];
If[SameQ[ss,{}],pp2,
pp3=Map[Last,StringPosition[pp2,")("]];
pp4=Union[Flatten[StringPosition[pp2,"("]]];
pp3=Flatten[Table[Position[pp4,pp3[[i]]],{i,Length[pp3]}]];
pp3=Table[{pp4[[pp3[[i]]-1]],pp4[[pp3[[i]]]]},{i,Length[pp3]}];
pp3=Table[StringTake[pp2,pp3[[i]]],{i,Length[pp3]}];
pp4=10*Table[Length[StringPosition[pp3[[i]],"1"]],{i,Length[pp3]}];
pp4=Map[ToString,Table[Table[1,{i,pp4[[j]]}],{j,Length[pp4]}]];
pp4=StringReplace[StringReplace[pp4,{"{"->"(","}"->","}],{" "->""}];
Do[pp2=StringReplace[pp2,{pp3[[i]]->pp4[[i]]}],{i,Length[pp3]}]];
pp2=If[SameQ[pp,{}],pp2,StringJoin[pp1,pp2]];
pp2=StringReplace[pp2,{"(1)"->"1"}];
pp2]


fMakePolyhedral[Ul_String]:=Module[{pp,pp1,vv1,rr,vv,i},
pp=Union[Flatten[StringPosition[Ul,"."]]];
pp1=Union[Flatten[StringPosition[Ul,"*"]]];
vv=If[SameQ[Union[pp,pp1],{}],Ul,vv1=If[SameQ[StringTake[Ul,1],"."],\
StringJoin[Ul,StringJoin[Table[".",{i,6-Length[pp]}]]],
rr=If[MemberQ[{"6","8","9"},StringTake[Ul,1]],StringJoin[Ul,StringJoin[Table[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
".",{i,ToExpression[StringTake[Ul,1]]-Length[pp]-1}]]],StringJoin[Ul,\
StringJoin[Table[".",{i,ToExpression[StringTake[Ul,2]]-Length[pp]-1}]]]
]]];
While[Not[SameQ[StringPosition[vv,".."],{}]],vv=StringReplace[vv,{".."->".1."}\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]];
vv=If[SameQ[StringTake[vv,-1],"."],StringJoin[vv,"1"],vv];
vv=StringReplace[vv,{"*."->"*1."}];
vv
] 


fMakeConV[Ul_String,ll_List]:=Module[{pp,pp1,pp2,pp3,ss,i},
ss=fMakePolyhedral[Ul];
ss=fExpandCon[ss];
pp=Flatten[StringPosition[ss,"*"]];
pp2=If[SameQ[pp,{}],pp2=ss,
pp1=StringTake[ss,First[Union[pp]]];
pp2=StringTake[ss,{First[Union[pp]]+1,StringLength[ss]}]];
pp3=StringPosition[pp2,"1"];
Do[pp2=StringReplacePart[pp2,ToString[ll[[i]]],pp3[[i]]],{i,Length[ll]}];
pp2=StringReplace[pp2,{"1"->"i","2"->"1","0"->"-1"," 0"->"3"}];
pp2=StringReplace[StringReplace[StringReplace[StringReplace[pp2,{"(i)"->"i","(\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
1)"->"1","(-1)"->"-1"}],{"3"->" 0"}],{" 0"->" -1"}],{" -1."->" 0."}];
pp2=If[SameQ[StringTake[pp2,-3]," \
-1"],StringJoin[StringTake[pp2,StringLength[pp2]-3]," 0"],pp2]; 
pp2=If[SameQ[pp,{}],pp2,StringJoin[pp1,pp2]];
pp2]



fMakeConVirt[Ul_String]:=Module[{nn,tt,tt1,tt2,tt3,tt4,ll,i},
nn=If[SameQ[StringPosition[Ul,"#"],{}],fCreatePData[Ul][[1,1]],
fDToDDirectPD[Ul][[1,1]]];
ll=fBasVirt[nn]; 
tt=Table[fMakeConV[Ul,ll[[i]]],{i,Length[ll]}];
tt=Select[tt,SameQ[StringPosition[#,"i,i"],{}]&&SameQ[StringPosition[#,
"(-1,1"],{}]&&SameQ[StringPosition[#,"-1,1)"],{}]&&
SameQ[StringPosition[#,"(1,-1"],{}]&&
SameQ[StringPosition[#,",1,-1)"],{}]&&
SameQ[StringPosition[#,",-1,1,"],{}]&&
SameQ[StringPosition[#,",1,-1,"],{}] &];
tt1=Table[fBothPoly[fCabledJonesVirt[tt[[i]],3]],{i,Length[tt]}];
tt2=Complement[Union[tt1],{{1+2 A^4+A^8,1+2 A^4+A^8}}];
tt3=Map[First,Table[Flatten[Position[tt1,tt2[[i]]],1],{i,Length[tt2]}]];
tt3=Table[tt[[tt3[[i]]]],{i,Length[tt3]}];
tt3=Select[tt3,Not[SameQ[StringPosition[#,"i"],{}]] &];
tt3
]


(* ## ## ## ## ## ## # CABLE ## ## ## ## ## ## ## ## # *)

Cable[n_][L_PD]:=Module[{h,v,i,j,k,out}, k=0;
out=PD@@Flatten[Join@@(L/.X[a_,b_,c_,d_]:>(++k; 
Table[X[h[i,j-1,k],v[i,j,k],h[i,j,k],v[i-1,j,k]],{i,n},{j,n}]/.{h[i_,0,_]:>a[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i],h[i_,n,_]:>c[i]}/.If[d-b==1||b-d>1,{v[0,j_,_]:>d[j],v[n,j_,_]:>b[j]},{v[0,\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
j_,_]:>d[n+1-j],v[n,j_,_]:>b[n+1-j]}]))];
k=0;out/.((#->++k)&/@(List@@Union@@out))]


(* ## ## ## ## ## ## VIRTUAL LINKS ## ## ## ## ## ## *)

fStandardPoly[pp_]:=Module[{pp1,pp2,i,j},
pp1=2Table[Exponent[pp[[i]],q],{i,Length[pp]}];
pp2=If[pp1[[1]]==0,Join[{pp[[1]]},Table[Coefficient[pp[[i]],q^Exponent[pp[[i]]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
,q]],{i,2,Length[pp]}]],Table[Coefficient[pp[[i]],q^Exponent[pp[[i]],q]],{i,\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

Length[pp]}]];
pp1=Plus@@Table[pp2[[i]]q^pp1[[i]],{i,Length[pp1]}];
pp1=fBothPoly[fPolyNorm[pp1]];
pp1]


fDyePoly[Ul_,k_Integer]:=Module[{pp,pp1,pp2},
pp=fVirtPD[Ul];
pp1=Cable[k][pp];
pp1=fStandardPoly[Jones[pp1][q]];
pp1]


fMakeConVirtLink[Ul_String,k_Integer]:=Module[{nn,tt,tt1,tt2,tt3,tt4,ll,i},
nn=If[SameQ[StringPosition[Ul,"#"],{}],Plus@@fCreatePData[Ul][[1]],Plus@@\
fDToDDirectPD[Ul][[1]]];
 ll=fBasVirt[nn]; 
tt=Table[fMakeConV[Ul,ll[[i]]],{i,Length[ll]}];
tt=Select[tt,SameQ[StringPosition[#,"i,i"],{}]&&
SameQ[StringPosition[#,"(-1,1"],{}]&&
SameQ[StringPosition[#,"-1,1)"],{}]&&
SameQ[StringPosition[#,"(1,-1"],{}]&&
SameQ[StringPosition[#,",1,-1)"],{}]&&
SameQ[StringPosition[#,",-1,1,"],{}]&&
SameQ[StringPosition[#,",1,-1,"],{}] &];
tt1=Table[fDyePoly[tt[[i]],k],{i,Length[tt]}];
tt2=Complement[Union[tt1],{{1+5 q^2+10 q^4+10 q^6+5 q^8+q^10,1+5 q^2+10 \
q^4+10 q^6+5 q^8+q^10},{1+2 q^2+3 q^8+6 q^10+8 q^12+6 q^14+3 q^16+2 \
q^22+q^24,1+2 q^2+3 q^8+6 q^10+8 q^12+6 q^14+3 q^16+2 q^22+q^24}}];
tt3=Map[First,Table[Flatten[Position[tt1,tt2[[i]]],1],{i,Length[tt2]}]];
tt3=Table[tt[[tt3[[i]]]],{i,Length[tt3]}];
tt3=Select[tt3,Not[SameQ[StringPosition[#,"i"],{}]] &];
tt3
]


(* ## ## ## ## ## ## # KAUFFMAN VIRT POLY ## ## ## ## ## ## ## *)

 fAdequateMixed[Ul_,s_List,k_Integer]:=Module[{dd,tt,poc,cikl,pr,prvi,cyc,pp,\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
res,kk,i,j},
dd=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
tt=Flatten[Table[If[SameQ[s[[i]],1],{Sort[{dd[[i,1]],dd[[i,4]]}],Sort[{dd[[i,\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
2]],dd[[i,3]]}]},{Sort[{dd[[i,1]],dd[[i,2]]}],Sort[{dd[[i,3]],dd[[i,4]]}]}],{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i,Length[s]}],1];
kk=k;
poc=kk;
cikl={};
pr={};
prvi=tt[[First[First[Position[tt,kk]]]]];
kk=If[SameQ[prvi[[2]],kk],prvi[[1]],kk];
poc=kk;
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{prvi[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
While[Not[SameQ[pp[[1]],poc]],
kk=First[Flatten[Position[tt,pp[[1]]],1]];
prvi=tt[[kk]];
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{pp[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
res={Map[IntegerPart,Divide[cikl-1,2]]+1,pr}];
res=If[SameQ[res[[1,1]],res[[1,-1]]],{Drop[res[[1]],-1],res[[2]]},res];
res
] 

fAdequateSignedMixed[Ul_,s_List]:=Module[{cc,k,c1,res},
cc=Range[2Length[Ul]];
k=First[cc];
c1=fAdequateMixed[Ul,s,k];
cc=Complement[cc,c1[[2]]];
res={c1[[1]]};
While[Not[SameQ[cc,{}]],
k=First[cc];
c1=fAdequateMixed[Ul,s,k];
cc=Complement[cc,c1[[2]]];
res=Join[res,{c1[[1]]}]];
res
] 

fAllStates[Ul_]:=Module[{pd,nn,pp,res,s,rr,rr1,i,j},
pd=If[SameQ[Head[Ul],String],fConwayToPD[Ul],Ul];
nn=Length[Table[pd[[i,j]],{i,Length[pd]},{j,4}]];
pp=fVarP[nn];
pp=Reverse[Map[Last,Sort[Table[{Count[pp[[i]],0],pp[[i]]},{i,Length[pp]}]]]];
res=Table[
s=pp[[i]];
rr=fAdequateSignedMixed[pd,s];
rr1={rr,pp[[i]]},{i,Length[pp]}];
res]

fGenSign[Conway_String]:=Module[{BlockSign={},i,DL,par,fConSignList,Ham,\
VertSign,Dow={}},DL=Dowker[Conway];
Dow=If[Not[SameQ[StringPosition[Conway,"#"],{}]],fGenSignDirProd[Conway],
Dow=If[DL!={},Dow=Flatten[DL[[4]]];
VertSign=DL[[5]];
Do[BlockSign=Append[BlockSign,0],{i,1,Length[DL[[1]]]}];
If[Flatten[Union[Map[Flatten[StringPosition[Conway,#]]&,{"*","."}]]]=={},\
fConSignList=fConvert[Conway][[3]],fConSignList=fConvertPoly[Conway][[3]]];
Ham=Flatten[DL[[3]]];
Do[par=Flatten[Position[Abs[Ham],i]];
BlockSign[[i]]:=Sign[Ham[[par[[1]]]]Ham[[par[[2]]]]]*fConSignList[[i]];
VertSign[[par[[1]]]]=BlockSign[[i]]*VertSign[[par[[1]]]];
VertSign[[par[[2]]]]=BlockSign[[i]]*VertSign[[par[[2]]]],{i,1,Length[\
BlockSign]}];VertSign=Flatten[VertSign];
Dow=Map[(Sign[VertSign[[Flatten[Position[Abs[VertSign],#]][[1]]]]])&,Dow]];
Dow]]



fGenSignsPDNew[Ul_String]:=Module[{pd,ss,res,i},
    pd=fVirtPD[Ul];
    ss=KnotTheory`PositiveQ /@Table[pd[[i]],{i,Length[pd]}];
    res=Table[If[SameQ[ss[[i]],True],1,-1],{i,Length[ss]}];
    res]



fGenSignVirtfromPD[Ul_,vv_List]:=Module[{res,pp,pp2,ll,rr,rr1,i},
pp2=GaussCode[Ul];
pp2=Table[pp2[[i]],{i,Length[pp2]}];
ll=If[SameQ[Head[pp2[[1]]],List],Map[Length,pp2],Length[pp2]];
pp2=Abs[Flatten[pp2]];
pp2=If[SameQ[Head[pp2[[1]]],List],PD[GaussCode@@iteratedTake[Table[pp2[[i]]*(-\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
1)^i,{i,Length[pp2]}],ll]],PD[GaussCode[Table[pp2[[i]]*(-1)^i,{i,Length[pp2]}]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]]];
rr=fGenSignfromPD[Ul];
rr1=rr;
Do[rr=ReplacePart[rr,0,vv[[i]]],{i,Length[vv]}];
res={rr1,rr};
res
] 

fReduceReal[nn_List,rr_List]:=Module[{pp0,pp,pp1,pp2,pp3,pp4,pp5,pp6,ll,i},
pp1=If[SameQ[nn,{}],nn,
pp=Table[If[MemberQ[rr,nn[[i]]],1,0],{i,Length[nn]}];
pp1=Split[pp];
ll=Map[Length,pp1];
pp2=Range[Length[nn]];
pp2=iteratedTake[pp2,ll];
pp1=Flatten[Table[If[SameQ[pp1[[i,1]],0],pp2[[i]],If[OddQ[Length[pp1[[i]]]],{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
pp2[[i,1]]},{}]],{i,Length[pp2]}]];
pp1=Table[nn[[pp1[[i]]]],{i,Length[pp1]}];
pp1=If[SameQ[pp1,{}],pp1,If[MemberQ[rr,pp1[[1]]]&&MemberQ[rr,pp1[[-1]]]&&\
Length[pp1]>1,Drop[Drop[pp1,-1],1],pp1]];
pp2=Select[pp1,MemberQ[rr,#] &];
pp1=If[SameQ[Length[pp2],1],UnsortedComplement[pp1,pp2],pp1];
pp2=Union[Select[pp1,MemberQ[rr,#] &]];
pp3=If[SameQ[pp2,{}],pp2,Sort[Flatten[Table[Position[pp1,pp2[[i]]],{i,Length[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
pp2]}]]]];
pp4=If[SameQ[pp3,{}],pp3,{If[SameQ[pp3[[1]],1],{pp1[[-1]],pp1[[1]],pp1[[2]]},{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
}]}];
pp5=If[SameQ[pp3,{}],pp3,{If[SameQ[pp3[[-1]],Length[pp1]],{pp1[[-2]],pp1[[-1]]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
,pp1[[1]]},{}]}];
pp3=Complement[pp3,{1,Length[pp1]}];
pp5=Complement[Union[Table[{pp1[[pp3[[i]]-1]],pp1[[pp3[[i]]]],pp1[[pp3[[i]]+1]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]},{i,Length[pp3]}],pp4,pp5],{{}}];
pp6=Select[pp5,Not[SameQ[#[[1]],#[[3]]]] &];
pp2=Select[pp5,SameQ[#[[1]],#[[3]]] &];
pp3=If[SameQ[pp6,{}],pp6,Table[pp6[[i,2]],{i,Length[pp6]}]];
pp4=If[SameQ[pp2,{}],pp2,Table[pp2[[i,2]],{i,Length[pp2]}]];
pp4=Complement[pp3,pp4];
pp1=If[OddQ[Length[pp4]],UnsortedComplement[pp1,pp4],pp1]];
pp1
] 

fClearDouble[Ul_List]:=Module[{pp,i},
pp=If[SameQ[Ul,{}],pp={1},
pp=Split[Ul];
Do[pp=If[Length[pp[[i]]]>1&&OddQ[Length[pp[[i]]]],ReplacePart[pp,Union[pp[[i]]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
],i],pp],{i,Length[pp]}];
Do[pp=If[EvenQ[Length[pp[[i]]]],ReplacePart[pp,0,i],pp],{i,Length[pp]}];
pp=Select[pp,Not[SameQ[#,0]] &];
pp=Flatten[pp];
pp=If[Length[pp]>1,If[SameQ[pp[[1]],pp[[-1]]],Drop[Drop[pp,1],-1],pp],pp];
pp];
pp
] 


fReduceGraphsVirt[tt_List,vv_List,rr_List]:=Module[{pp,res,ss,ppss,ss1,qq,mm,mm1,\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i,j,k,m,n},
res=Table[Table[
ss={{tt[[n,m]]}};
ss=Table[pp=Flatten[Table[If[Not[SameQ[ss[[k,j]],{}]],Map[Sort,Join[Table[{ss[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
[k,j,i]],ss[[k,j,i+1]]},{i,Length[ss[[k,j]]]-1}],{{ss[[k,j,1]],ss[[k,j,-1]]}}]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
],{}],{j,Length[ss[[k]]]}],1];
Select[pp,Length[Intersection[#,rr]]<2 &],{k,Length[ss]}];
qq=Union[Flatten[ss]];
ss=Table[Select[ss[[i]],Not[SameQ[#[[1]],#[[2]]]] &],{i,Length[ss]}];
ss=Table[If[Length[ss[[i]]]>2,ss[[i]],{}],{i,Length[ss]}];
ss=Table[If[Not[SameQ[ss[[i]],{}]],
mm=Union[Flatten[ss[[i]]]];
mm1=Range[Length[mm]];
mm=Table[mm[[i]]->mm1[[i]],{i,Length[mm]}];
{Sort[ReplaceAll[ss[[i]],mm]],ReplaceAll[Intersection[vv,qq],mm]},{1}],{i,\
Length[ss]}];
ss1=Union[Flatten[Map[Last,ss]]];
ss=Flatten[Map[First,ss],1];
ss1=Intersection[Union[Flatten[ss]],ss1];
res={ss,ss1};
res,{m,Length[tt[[n]]]}],{n,Length[tt]}];
res
] 


fDelVirtLoops[Ul_List,vv_List,rr_List]:=Module[{res,pp,pp1,pp2,pp3,pp4,pp5},
res=If[Length[Ul]<3,Ul,
pp1=Sort[Flatten[Table[Position[Ul,vv[[i]]],{i,Length[vv]}]]];
pp2={If[MemberQ[pp1,1],{Ul[[-1]],Ul[[1]],Ul[[2]]},{}]};
pp3={If[MemberQ[pp1,Length[Ul]],{Ul[[1]],Ul[[-1]],Ul[[-2]]},{}]};
pp4=Complement[pp1,{1,Length[Ul]}];
pp5=Table[{Ul[[pp4[[i]]-1]],Ul[[pp4[[i]]]],Ul[[pp4[[i]]+1]]},{i,Length[pp4]}];\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

pp5=Complement[Union[pp2,pp3,pp5],{{}}];
pp5=Select[pp5,SameQ[#[[1]],#[[3]]]&&MemberQ[rr,#[[1]]] &];
pp5=Table[pp5[[i,2]],{i,Length[pp5]}];
res=UnsortedComplement[Ul,pp5]];
res
] 

fFinalGraphsVirt[Ul_List]:=Module[{pp,rr,pp0,pp1,mm,mm1,res,i},
res=If[SameQ[Ul,{}]||SameQ[Union[Flatten[Ul]],{1}],{},
pp=Complement[Flatten[Ul[[1]]],Ul[[2]]];
rr=pp;
pp0=Select[Ul[[1]],SameQ[Length[Intersection[#,pp]],2] &];
res=If[SameQ[pp0,{}],Ul,
pp=Select[Ul[[1]],SameQ[Length[Intersection[#,pp]],1] &];
pp1=
Table[If[MemberQ[rr,pp[[i,1]]],pp[[i]],Reverse[pp[[i]]]],{i,Length[pp]}];
pp1=ReplaceAll[pp0,Table[pp1[[i,1]]->pp1[[i,2]],{i,Length[pp1]}]];
pp1=Join[UnsortedComplement[Map[Sort,Ul[[1]]],Map[Sort,Join[pp0,pp]]],pp1];
pp={pp1,Ul[[2]]};
mm=Union[Flatten[pp[[1]]]];
mm1=Range[Length[mm]];
res=ReplaceAll[pp,Table[mm[[i]]->mm1[[i]],{i,Length[mm1]}]]]];
res
] 

fKauffmanVirt[Ul_String]:=Module[{pom,cc,red,qq,rr,dd,rr1,
rr2,rr3,tt,tt1,tt2,tt3,
vv,pp,pp1,ind,ss,ss1,ss2,sss,mm,mm1,res,wr,i,j,k},
rr1=fGenSignsPD[Ul];
rr2=fGenSignsPD[StringReplace[Ul,{"-1"->"1"}]];
cc=Divide[Abs[rr1-rr2],2];
dd=StringReplace[Ul,{"i"->"1"}];
rr=fAllStates[dd];
tt=Map[Last,rr];
tt=Table[tt[[i]]*rr1[[2]],{i,Length[tt]}];
tt1=Union[tt];
tt2=Map[Flatten,Table[Position[tt,tt1[[i]]],{i,Length[tt1]}]];
rr2=Map[Length,Map[First,rr]];
rr2=Table[{rr2[[i]],rr[[i]],tt[[i]]},{i,Length[rr]}];
rr2=Map[First,Table[Sort[Table[rr2[[tt2[[j,i]]]],{i,Length[tt2[[j]]]}]],{j,\
Length[tt2]}]];
vv=Flatten[Position[rr1[[2]],0]]; (* virtualna temena *)
rr=Complement[Range[Length[rr1[[1]]]],vv]; (* realna temena *)
rr3=Table[{rr2[[i,1]]-1,rr2[[i,2]]},{i,Length[rr2]}];
dd=Map[Last,Map[Last,rr3]];
tt=ReplaceAll[Table[Table[{rr1[[1,i]],dd[[j,i]]},{i,Length[rr1[[1]]]}],{j,\
Length[dd]}],{{-1,1}->0,{-1,0}->1,{1,1}->1,{1,0}->0}];
tt1=ReplaceAll[Table[Table[{rr1[[1,i]],tt[[j,i]]},{i,Length[rr1[[1]]]}],{j,\
Length[tt]}],{{1,0}->1,{-1,1}->1,{1,1}->-1,{-1,0}->-1}];
tt2=-Table[Plus@@Table[tt1[[j,rr[[i]]]],{i,Length[rr]}],{j,Length[tt1]}]; (* \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
znaci *)
rr3=Table[{rr3[[i,1]],tt2[[i]],rr3[[i,2]]},{i,Length[tt2]}];
Do[rr3=ReplacePart[rr3,tt[[i]],{i,3,2}],{i,Length[tt]}];
ss=Table[{rr3[[i,3,1]],Complement[Flatten[Position[rr3[[i,3,2]],0]],vv]},{i,\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

Length[rr3]}];
ss=Table[Table[UnsortedComplement[ss[[j,1,i]],Complement[rr,ss[[j,-1]]]],{i,\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

Length[ss[[j,1]]]}],{j,Length[ss]}];
ind=1;
While[ind>0,
ss1=Flatten[ss,1];
ss1=Min[Table[Length[Union[If[Not[SameQ[ss1[[i]],{}]],{ss1[[i,1]],ss1[[i,-1]]}\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
,{}]]],{i,Length[ss1]}]];
ss2=Max[Map[Length,Flatten[Map[Split,Flatten[ss,1]],1]]];
ind=If[SameQ[ss1,1]||ss2>1,1,0];
ss=Table[If[SameQ[ss[[i]],{}],{},Map[fClearDouble,ss[[i]]]],{i,Length[ss]}];
ss=Table[Select[ss[[i]],Length[#]>1 &],{i,Length[ss]}];
ss=Table[Select[ss[[i]],Not[SameQ[Sort[#],Intersection[#,rr]]] \
&],{i,Length[ss]}];
ss];
pom=Map[Last,Map[Last,rr3]];
pom=Table[Complement[rr,Complement[Flatten[Position[pom[[i]],0]],vv]],{i,\
Length[pom]}];
ss=Table[Table[fDelVirtLoops[ss[[j,i]],vv,pom[[j]]],{i,Length[ss[[j]]]}],{j,\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

Length[ss]}];
ss=Table[Table[If[SameQ[ss[[j,i]],{}],ss[[j,i]],fReduceReal[ss[[j,i]],rr]],{i,\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
Length[ss[[j]]]}],{j,Length[ss]}];
red=fReduceGraphsVirt[ss,vv,rr];
red=Map[Reverse,Table[Table[fFinalGraphsVirt[red[[j,i]]],{i,Length[red[[j]]]}]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
,{j,Length[red]}]];
red=ReplaceAll[red,{}->1];
ss=Table[pp=Flatten[Table[If[Not[SameQ[ss[[k,j]],{}]],Map[Sort,Join[Table[{ss[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
[k,j,i]],ss[[k,j,i+1]]},{i,Length[ss[[k,j]]]-1}],{{ss[[k,j,1]],ss[[k,j,-1]]}}]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
],{}],{j,Length[ss[[k]]]}],1];
Select[pp,Length[Intersection[#,rr]]<2 &],{k,Length[ss]}];
ss=Table[Select[ss[[i]],Not[SameQ[#[[1]],#[[2]]]] &],{i,Length[ss]}];
ss=Table[If[Length[ss[[i]]]>2,ss[[i]],{}],{i,Length[ss]}];
ss=Table[If[Not[SameQ[ss[[i]],{}]],
mm=Union[Flatten[ss[[i]]]];
mm1=Range[Length[mm]];
mm=Table[mm[[i]]->mm1[[i]],{i,Length[mm]}];
{Sort[ReplaceAll[ss[[i]],mm]],ReplaceAll[vv,mm]},{1}],{i,Length[ss]}];
ss1=Map[Last,ss];
ss=Map[First,ss];
tt=Table[
dd=rr3[[j,3,1]];
Table[Select[dd[[i]],SameQ[Intersection[{#},rr],{}] \
&],{i,Length[dd]}],{j,Length[rr3]}];
tt=Table[Table[Flatten[Select[Split[tt[[j,i]]],SameQ[Length[#],1] \
&]],{i,Length[tt[[j]]]}],{j,Length[tt]}];
tt=Table[Complement[tt[[i]],{{}}],{i,Length[tt]}];
tt=Table[Select[tt[[i]],Length[Union[#]]>1 &],{i,Length[tt]}];
tt=Table[Table[pp=tt[[j,i]];
pp1=Union[pp];
rr=Range[Length[pp1]];
pp=ReplaceAll[pp,Table[pp1[[i]]->rr[[i]],{i,Length[rr]}]],
{i,Length[tt[[j]]]}],{j,Length[tt]}];
tt=Table[If[Not[SameQ[tt[[k]],{}]],Table[pp=tt[[k,j]];
pp1=Join[Table[Sort[{pp[[i]],pp[[i+1]]}],{i,Length[pp]-1}],{Sort[{pp[[1]],pp[[
-1]]}]}],{j,Length[tt[[k]]]}],{}],{k,Length[tt]}];
tt=ReplaceAll[Map[Union,Table[Table[If[SameQ[Select[tt[[j,i]],SameQ[#[[1]],
#[[
2]]] &],{}],tt[[j,i]],{}],{i,Length[tt[[j]]]}],{j,Length[tt]}]],{{}}->{}];
(* Print["OVDE"];
Print[tt]; *)
tt=Table[ReplaceAll[Length[Flatten[Map[First,Map[fKLfromGraph,tt[[j]]]]]],
0->1]-1,{j,Length[tt]}];
(* Print["OVDE1"]; *)
tt=Table[{rr3[[i,1]]+tt[[i]],rr3[[i,2]]},{i,Length[tt]}];
Print[Table[{d^tt[[i,1]]A^tt[[i,2]],red[[i]]},{i,Length[tt]}]]; 
tt1=Plus@@Table[A^tt[[i,2]]*d^tt[[i,1]],{i,Length[tt2]}];
dd=-A^2-A^{-2};
wr=Plus@@fGenSignsPD[Ul][[2]];
tt2=Expand[Plus@@Table[A^tt[[i,2]]*dd^tt[[i,1]],{i,Length[tt2]}]];
tt3=Expand[tt2*(-A^{-3})^-wr];
res={tt1,tt2,tt3};
res
] 


fGenSignfromPD[Ul_,vv_List]:=Module[{pp0,gg,gg1,gg2,pp1,res,res1,i,j},
pp0=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
gg=GaussCode[Ul];
gg1=Table[gg[[i]],{i,Length[gg]}];
gg2=Map[Max,Abs[gg1]];
pp1=Flatten[Table[{{pp0[[i,1]],pp0[[i,3]]},{pp0[[i,2]],pp0[[i,4]]}},{i,Length[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
pp0]}],1];
pp1=Partition[Table[If[Not[SameQ[Intersection[pp1[[i]],gg2],{}]]&&Not[SameQ[\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

Abs[pp1[[i,1]]-pp1[[i,2]]],1]],Reverse[pp1[[i]]],pp1[[i]]],{i,Length[pp1]}],2]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
;
res=Table[If[pp1[[i,1,1]]<pp1[[i,1,2]]&&pp1[[i,2,1]]<pp1[[i,2,2]]||pp1[[i,1,1]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
]>pp1[[i,1,2]]&&pp1[[i,2,1]]>pp1[[i,2,2]],1,-1],{i,Length[pp1]}];
res1=res;
Do[res=ReplacePart[res,0,vv[[i]]],{i,Length[vv]}];
{res1,res}
] 

fGenSignfromPDAdditional[Ul_,vv_List]:=Module[{pp1,pp0,gg,gg1,gg2,i},
pp0=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
gg=GaussCode[Ul];
gg=Table[gg[[i]],{i,Length[gg]}];
gg2=Map[Length,gg];
gg1=Table[(-1)^i,{i,Length[gg]}];
gg1=gg*gg1;
gg1=Sort[Abs[Union[Select[gg1,Sign[#]<0 &]]]];
pp0=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
pp1=Table[If[MemberQ[gg1,i],RotateRight[pp0[[i]]],pp0[[i]]],{i,Length[pp0]}];
pp1=Table[If[MemberQ[gg1,i],-1,1],{i,Length[pp0]}];
pp0=Abs[fGenSignVirtfromPD[Ul,vv]];
pp1={pp0[[1]]*pp1,pp0[[2]]*pp1};
pp1=ReplaceAll[pp1,{1->0,-1->1}];
pp1
] 

fKauffmanVirtfromPD[Ul_,vv_List]:=Module[{mm,mm1,ind,ss1,ss2,red,cc,pp0,rr,dd,rr1,rr2,
rr3,tt,tt1,tt2,tt3,ss,pom,pp,pp1,pp2,ll,res,wr,i,j,k},
rr1=fGenSignfromPD[Ul,vv];
cc=fGenSignfromPDAdditional[Ul,vv];
rr=fAllStates[Ul];
tt=Map[Last,rr];
tt=Table[tt[[i]]*Abs[rr1[[2]]],{i,Length[tt]}];
tt1=Union[tt];
tt2=Map[Flatten,Table[Position[tt,tt1[[i]]],{i,Length[tt1]}]];
rr2=Map[Length,Map[First,rr]];
rr2=Table[{rr2[[i]],rr[[i]],tt[[i]]},{i,Length[rr]}];
rr2=Map[First,Table[Sort[Table[rr2[[tt2[[j,i]]]],{i,Length[tt2[[j]]]}]],{j,\
Length[tt2]}]];
rr=Complement[Range[Length[rr1[[1]]]],vv]; (* realna temena *)
rr3=Table[{rr2[[i,1]]-1,rr2[[i,2]]},{i,Length[rr2]}];
dd=Map[Last,Map[Last,rr3]];
tt=ReplaceAll[Table[Table[{rr1[[1,i]],dd[[j,i]]},{i,Length[rr1[[1]]]}],{j,\
Length[dd]}],{{-1,1}->0,{-1,0}->1,{1,1}->1,{1,0}->0}];
tt1=ReplaceAll[Table[Table[{rr1[[1,i]],tt[[j,i]]},{i,Length[rr1[[1]]]}],{j,\
Length[tt]}],{{1,0}->1,{-1,1}->1,{1,1}->-1,{-1,0}->-1}];
 tt2=-Table[Plus@@Table[tt1[[j,rr[[i]]]],{i,Length[rr]}],{j,Length[tt1]}]; (* \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
znaci *)
rr3=Table[{rr3[[i,1]],tt2[[i]],rr3[[i,2]]},{i,Length[tt2]}];
Do[rr3=ReplacePart[rr3,tt[[i]],{i,3,2}],{i,Length[tt]}];
ss=Table[{rr3[[i,3,1]],Complement[Flatten[Position[rr3[[i,3,2]],0]],
vv]},{i,
Length[rr3]}];
ss=Table[Table[UnsortedComplement[ss[[j,1,i]],Complement[rr,ss[[j,-1]]]],
{i,Length[ss[[j,1]]]}],{j,Length[ss]}];
ind=1;
While[ind>0,
ss1=Flatten[ss,1];
ss1=Min[Table[Length[Union[If[Not[SameQ[ss1[[i]],{}]],{ss1[[i,1]],
ss1[[i,-1]]}
,{}]]],{i,Length[ss1]}]];
ss2=Max[Map[Length,Flatten[Map[Split,Flatten[ss,1]],1]]];
ind=If[SameQ[ss1,1]||ss2>1,1,0];
ss=Table[If[SameQ[ss[[i]],{}],{},Map[fClearDouble,ss[[i]]]],{i,Length[ss]}];
ss=Table[Select[ss[[i]],Length[#]>1 &],{i,Length[ss]}];
ss=Table[Select[ss[[i]],Not[SameQ[Sort[#],Intersection[#,rr]]] \
&],{i,Length[ss]}];
ss];
pom=Map[Last,Map[Last,rr3]];
pom=Table[Complement[rr,Complement[Flatten[Position[pom[[i]],0]],vv]],{i,\
Length[pom]}];
ss=Table[Table[fDelVirtLoops[ss[[j,i]],vv,pom[[j]]],{i,Length[ss[[j]]]}],
{j,Length[ss]}];
ss=Table[Table[If[SameQ[ss[[j,i]],{}],ss[[j,i]],fReduceReal[ss[[j,i]],rr]],
{i,Length[ss[[j]]]}],{j,Length[ss]}];
red=fReduceGraphsVirt[ss,vv,rr];
red=Map[Reverse,Table[Table[fFinalGraphsVirt[red[[j,i]]],
{i,Length[red[[j]]]}]
,{j,Length[red]}]];
red=ReplaceAll[red,{}->1]; 
ss=Table[pp=Flatten[Table[If[Not[SameQ[ss[[k,j]],{}]],
Map[Sort,Join[Table[{ss[
[k,j,i]],ss[[k,j,i+1]]},{i,Length[ss[[k,j]]]-1}],{{ss[[k,j,1]],
ss[[k,j,-1]]}}]],{}],{j,Length[ss[[k]]]}],1];
Select[pp,Length[Intersection[#,rr]]<2 &],{k,Length[ss]}];
ss=Table[Select[ss[[i]],Not[SameQ[#[[1]],#[[2]]]] &],{i,Length[ss]}];
ss=Table[If[Length[ss[[i]]]>2,ss[[i]],{}],{i,Length[ss]}];
ss=Table[If[Not[SameQ[ss[[i]],{}]],
mm=Union[Flatten[ss[[i]]]];
mm1=Range[Length[mm]];
mm=Table[mm[[i]]->mm1[[i]],{i,Length[mm]}];
{Sort[ReplaceAll[ss[[i]],mm]],ReplaceAll[vv,mm]},{1}],{i,Length[ss]}];
ss1=Map[Last,ss];
ss=Map[First,ss];
tt=Table[
dd=rr3[[j,3,1]];
Table[Select[dd[[i]],SameQ[Intersection[{#},rr],{}] \
&],{i,Length[dd]}],{j,Length[rr3]}];
tt=Table[Table[Flatten[Select[Split[tt[[j,i]]],SameQ[Length[#],1] \
&]],{i,Length[tt[[j]]]}],{j,Length[tt]}];
tt=Table[Complement[tt[[i]],{{}}],{i,Length[tt]}];
tt=Table[Select[tt[[i]],Length[Union[#]]>1 &],{i,Length[tt]}];
tt=Table[Table[pp=tt[[j,i]];
pp1=Union[pp];
rr=Range[Length[pp1]];
pp=ReplaceAll[pp,Table[pp1[[i]]->rr[[i]],{i,Length[rr]}]],
{i,Length[tt[[j]]]}],{j,Length[tt]}];
tt=Table[If[Not[SameQ[tt[[k]],{}]],Table[pp=tt[[k,j]];
pp1=Join[Table[Sort[{pp[[i]],pp[[i+1]]}],{i,Length[pp]-1}],{Sort[{pp[[1]],
pp[[-1]]}]}],{j,Length[tt[[k]]]}],{}],{k,Length[tt]}];
tt=ReplaceAll[Map[Union,Table[Table[If[SameQ[Select[tt[[j,i]],
SameQ[#[[1]],#[[
2]]] &],{}],tt[[j,i]],{}],{i,Length[tt[[j]]]}],{j,Length[tt]}]],{{}}->{}];
tt=Table[ReplaceAll[Length[Flatten[Map[First,
Map[fKLfromGraph,tt[[j]]]]]],0->
1]-1,{j,Length[tt]}];
tt=Table[{rr3[[i,1]]+tt[[i]],rr3[[i,2]]},{i,Length[tt]}];
Print[Table[{d^tt[[i,1]]A^tt[[i,2]],red[[i]]},{i,Length[tt]}]]; 
tt1=Plus@@Table[A^tt[[i,2]]*d^tt[[i,1]],{i,Length[tt2]}];
dd=-A^2-A^{-2};
wr=Plus@@fGenSignVirtfromPD[Ul,vv][[2]];
dd=-A^2-A^{-2};
tt2=Expand[Plus@@Table[A^tt[[i,2]]*dd^tt[[i,1]],{i,Length[tt2]}]];
tt3=Expand[tt2*(-A^{-3})^-wr];
res={tt1,tt2,tt3};
res
] 


(* ## ## ## ## ## ADEQUATE ## ## ## ## ## ## ## ## # *)

fAllCycles[LL_List]:=Module[{g1,g2,pp,res,i,j,k},
res=Complement[Table[
g1=Map[Flatten,KSubsets[LL,k]];
g2=Map[Union,g1];
pp=Flatten[Position[Table[Union[Table[Count[g1[[j]],g2[[j,i]]],
{i,Length[g2[[
j]]]}]],{j,Length[g1]}],{2}]];
g2=Table[g1[[pp[[i]]]],{i,Length[pp]}];
g2=Map[FindCycle,Map[FromUnorderedPairs,Table[Partition[g2[[i]],2],
{i,Length[
g2]}]]],
{k,Length[LL]}],{{}}];
res=Union[Flatten[res,1]];
res] 

fAdequate[Ul_,s_Integer,k_Integer]:=Module[{dd,tt,poc,cikl,pr,prvi,
cyc,pp,res,
kk,i,j},
dd=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
tt=Flatten[If[SameQ[s,1],Table[{Sort[{dd[[i,1]],dd[[i,4]]}],
Sort[{dd[[i,2]],
dd[[i,3]]}]},{i,Length[dd]}],
Table[{Sort[{dd[[i,1]],dd[[i,2]]}],Sort[{dd[[i,3]],dd[[i,4]]}]},
{i,Length[dd]}]],1];
kk=k;
poc=kk;
cikl={};
pr={};
prvi=tt[[First[First[Position[tt,kk]]]]];
kk=If[SameQ[prvi[[2]],kk],prvi[[1]],kk];
poc=kk;
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{prvi[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
While[Not[SameQ[pp[[1]],poc]],
kk=First[Flatten[Position[tt,pp[[1]]],1]];
prvi=tt[[kk]];
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{pp[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
res={Map[IntegerPart,Divide[cikl-1,2]]+1,pr}];
res=If[SameQ[res[[1,1]],res[[1,-1]]],{Drop[res[[1]],-1],res[[2]]},res];
res
] 

fAdequateSigned[Ul_,s_Integer]:=Module[{cc,k,c1,res},
cc=Range[2Length[Ul]];
k=First[cc];
c1=fAdequate[Ul,s,k];
cc=Complement[cc,c1[[2]]];
res={c1[[1]]};
While[Not[SameQ[cc,{}]],
k=First[cc];
c1=fAdequate[Ul,s,k];
cc=Complement[cc,c1[[2]]];
res=Join[res,{c1[[1]]}]];
res
] 

fAdequateTest[Ul_]:=Module[{pd,pl,pp,pp1,respl,mi,qq,qq1,resmi,res,i},
pd=If[SameQ[Head[Ul],String],fConwayToPD[Ul],fPdataToPD[Ul]];
pl=fAdequateSigned[pd,1];
pp=Map[Length,pl];
pp1=Map[Length,Table[Union[pl[[i]]],{i,Length[pl]}]];
respl=If[SameQ[pp,pp1],1,0];
mi=fAdequateSigned[pd,2];
qq=Map[Length,mi];
qq1=Map[Length,Table[Union[mi[[i]]],{i,Length[mi]}]];
resmi=If[SameQ[qq,qq1],-1,0];
res=Sort[Abs[{respl,resmi}]];
res
] 


fCycleNo[Ul_]:=Module[{res},
res={Length[fAdequateSigned[fConwayToPD[Ul],1]],Length[fAdequateSigned[fConwayToPD[Ul],2]]};
res]


fAdequateMixed[Ul_,s_List,k_Integer]:=Module[{dd,tt,poc,cikl,
pr,prvi,cyc,pp,
res,kk,i,j},
dd=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
tt=Flatten[Table[If[SameQ[s[[i]],1],{Sort[{dd[[i,1]],dd[[i,4]]}],
Sort[{dd[[i,
2]],dd[[i,3]]}]},{Sort[{dd[[i,1]],dd[[i,2]]}],Sort[{dd[[i,3]],
dd[[i,4]]}]}],{i,Length[s]}],1];
kk=k;
poc=kk;
cikl={};
pr={};
prvi=tt[[First[First[Position[tt,kk]]]]];
kk=If[SameQ[prvi[[2]],kk],prvi[[1]],kk];
poc=kk;
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{prvi[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
While[Not[SameQ[pp[[1]],poc]],
kk=First[Flatten[Position[tt,pp[[1]]],1]];
prvi=tt[[kk]];
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{pp[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
res={Map[IntegerPart,Divide[cikl-1,2]]+1,pr}];
res=If[SameQ[res[[1,1]],res[[1,-1]]],{Drop[res[[1]],-1],res[[2]]},res];
res
] 


fReducedStateGraphs[Con_String]:=Module[{gg,gg1,ss,ss1,res},
gg=fStateGraph[fAdequateSigned[fConwayToPD[Con],1]];
gg1=fStateGraph[fAdequateSigned[fConwayToPD[Con],2]];
ss=fReducedGraph[ToUnorderedPairs[gg]];
ss1=fReducedGraph[ToUnorderedPairs[gg1]];
res={ss,ss1};
res]


fReducedAdequateStateGraphs[Con_String,ii_Integer]:=Module[{gr},
gr=If[ii<=Plus@@fCreatePData[Con][[1]],
fReducedGraph[ToUnorderedPairs[fAllAdequateStatesGraphs[Con,ii]]],{}];
gr]


fCriticalLines[Ul_]:=Module[{mm,mir,mm1,mm2,ss,neg,cr,pcr,ncr,res,i},
    mm=KnotTheory`Kauffman[fConwayToPD[Ul]][a,x];
    mir=KnotTheory`Kauffman[KnotTheory`Mirror[fConwayToPD[Ul]]][a,x];
    mm1=Table[mm[[i]],{i,Length[mm]}];
    mm2=Table[{Exponent[mm1[[i]],a],Exponent[mm1[[i]],x]},{i,Length[mm1]}];
    ss=Union[Map[Sign,Map[First,mm2]]];
    ss=MemberQ[ss,-1];
        mm2=If[SameQ[ss,True],mm2,mm=KnotTheory`Kauffman[KnotTheory`Mirror[fConwayToPD[Ul]]][a,x];
            mm1=Table[mm[[i]],{i,Length[mm]}];            
            mm2=Table[{Exponent[mm1[[i]],a],Exponent[mm1[[i]],x]},{i,Length[mm1]}]];
    neg=If[SameQ[ss,True],Count[Map[Sign,fGenSignfromPDNew[fConwayToPD[Ul]]],-1],Plus@@Dow[Ul][[1]]-Count[Map[Sign,fGenSignfromPD[fConwayToPD[Ul]]],-1]];
    cr=Plus@@Dow[Ul][[1]];
    pcr=Select[mm2,SameQ[#[[1]],#[[2]]-2neg] &];
    ncr=Select[mm2,SameQ[#[[1]],2cr-2neg-#[[2]]] &];
    pcr=Flatten[Table[Position[mm2,pcr[[i]]],{i,Length[pcr]}]];
    ncr=Flatten[Table[Position[mm2,ncr[[i]]],{i,Length[ncr]}]];
    pcr=If[SameQ[ss,True],Table[mm[[pcr[[i]]]],{i,Length[pcr]}],Table[mir[[pcr[[i]]]],{i,Length[pcr]}]];
    ncr=If[SameQ[ss,True],Table[mm[[ncr[[i]]]],{i,Length[ncr]}],Table[mir[[ncr[[i]]]],{i,Length[ncr]}]];
    res={pcr,ncr};
    res
    ] 


fAdequateSignedMixed[Ul_,s_List]:=Module[{cc,k,c1,res},
cc=Range[2Length[Ul]];
k=First[cc];
c1=fAdequateMixed[Ul,s,k];
cc=Complement[cc,c1[[2]]];
res={c1[[1]]};
While[Not[SameQ[cc,{}]],
k=First[cc];
c1=fAdequateMixed[Ul,s,k];
cc=Complement[cc,c1[[2]]];
res=Join[res,{c1[[1]]}]];
res
] 


fAdequateTestMixed[Ul_,s_List]:=Module[{pd,pl,pp,pp1,respl,mi,qq,qq1,resmi,\
res,i},
pd=If[SameQ[Head[Ul],String],fConwayToPD[Ul],Ul];
pl=fAdequateSignedMixed[pd,s];
pp=Map[Length,pl];
pp1=Map[Length,Table[Union[pl[[i]]],{i,Length[pl]}]];
respl=If[SameQ[pp,pp1],1,0];
respl
] 

fAdequateMixedFast[Ul_,s_List,k_Integer]:=Module[{dd,tt,poc,cikl,
pr,prvi,cyc,
pp,res,kk,ind,i,j},
dd=Table[Ul[[i,j]],{i,Length[Ul]},{j,4}];
tt=Flatten[Table[If[SameQ[s[[i]],1],{Sort[{dd[[i,1]],dd[[i,4]]}],
Sort[{dd[[i,
2]],dd[[i,3]]}]},{Sort[{dd[[i,1]],dd[[i,2]]}],Sort[{dd[[i,3]],dd[[i,4]]}]}],
{i,Length[s]}],1];
kk=k;
poc=kk;
cikl={};
pr={};
prvi=tt[[First[First[Position[tt,kk]]]]];
kk=If[SameQ[prvi[[2]],kk],prvi[[1]],kk];
poc=kk;
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
pp=Complement[prvi,{prvi[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
ind=1;
While[Not[SameQ[pp[[1]],poc]]&&SameQ[ind,1],
kk=First[Flatten[Position[tt,pp[[1]]],1]];
prvi=tt[[kk]];
pr=Union[pr,prvi];
cyc=Flatten[Position[tt,prvi]];
cikl=Join[cyc,cikl];
cikl=If[SameQ[cikl[[1]],cikl[[-1]]],Drop[cikl,-1],cikl];
pp=Complement[prvi,{pp[[1]]}];
tt=ReplacePart[tt,{0,0},Flatten[Position[tt,prvi]][[1]]];
res={Map[IntegerPart,Divide[cikl-1,2]]+1,pr};
ind=SameQ[Sort[Map[IntegerPart,Divide[cikl-1,2]]+1],Union[Map[IntegerPart,\
Divide[cikl-1,2]]+1]];
ind=If[SameQ[ind,False],0,1];
res=If[SameQ[ind,0],{},res]
]; (* kraj while *)
res=If[SameQ[res,{}],{},
res=If[SameQ[res[[1,1]],res[[1,-1]]],{Drop[res[[1]],-1],res[[2]]},res]];
res
] 

fAdequateSignedMixedFast[Ul_,s_List]:=Module[{cc,k,c1,res},
cc=Range[2Length[Ul]];
k=First[cc];
c1=fAdequateMixedFast[Ul,s,k];
cc=If[SameQ[c1,{}],{},Complement[cc,c1[[2]]]];
res=If[SameQ[c1,{}],{},{c1[[1]]}];
While[Not[SameQ[cc,{}]],
k=First[cc];
c1=fAdequateMixedFast[Ul,s,k];
cc=If[SameQ[c1,{}],{},Complement[cc,c1[[2]]]];
res=If[SameQ[c1,{}],{},Join[res,{c1[[1]]}]]];
res
] 

fStateGraph[ll_List]:=Module[{dd,ddd,dd1,i,j},
dd=If[SameQ[Length[ll],1],FromUnorderedPairs[Map[Sort,Join[{{ll[[1,1]],
ll[[1,-1]]}},Table[{ll[[1,i]],ll[[1,i+1]]},{i,Length[ll[[1]]]-1}]]]],
dd=KSubsets[Range[Length[ll]],2];
dd1=Table[Intersection[ll[[dd[[i,1]]]],ll[[dd[[i,2]]]]],{i,Length[dd]}];
dd=Flatten[Complement[Table[Table[dd[[j]],{i,Length[dd1[[j]]]}],
{j,Length[dd]}],{{}}],1];
ddd=Table[Union[Select[ll[[i]],Count[ll[[i]],#]>1 &]],{i,Length[ll]}];
ddd=Flatten[Table[If[Length[ddd[[i]]]>0,Partition[Table[Length[ddd[[i]]],{j,2Length[ddd[[i]]]}],2],{}],{i,Length[ll]}],1];
dd=Join[dd,ddd];
dd=FromUnorderedPairs[dd]];
dd]

fAllAdequateStatesGraphs[Con_String,ii_Integer]:=Module[{ll,gr,ll1},
ll=fAllAdequateStatesFast[Con];
gr=If[ii<=Length[ll],
ll1=ll[[ii]];
gr=fStateGraph[ll1[[1]]],{}];
gr]

fMixedStateGraph[Con_,ss_List]:=Module[{gg,ss1,rr},
gg=If[SameQ[Plus@@fCreatePData[Con][[1]],Length[ss]],
ss1=fConwayToPD[Con];
rr=fAdequateSignedMixed[ss1,ss];
gg=fStateGraph[rr],{}];
gg]

fStateGraphReducedOneDigon[ll_List]:=Module[{ll1,ll2,i},
ll1=Union[ll];
ll2=Table[Flatten[Position[ll,ll1[[i]]]],{i,Length[ll1]}];
ll1=Select[ll2,Length[#]>2 &];
ll2=Select[ll2,Length[#]>=1 &];
ll1=Flatten[Table[If[Length[ll2[[i]]]>=2,{ll2[[i,1]],ll2[[i,2]]},
{ll2[[i,1]]}],{i,Length[ll2]}]];
ll1=Sort[Table[ll[[ll1[[i]]]],{i,Length[ll1]}]];
ll1=FromUnorderedPairs[ll1];
ll1
]  


fAllAdequateStatesFast[Ul_]:=Module[{pd,nn,pp,res,s,rr,rr1,i,j},
pd=If[SameQ[Head[Ul],String],fConwayToPD[Ul],Ul];
nn=Length[Table[pd[[i,j]],{i,Length[pd]},{j,4}]];
pp=fVarP[nn];
pp=Reverse[Map[Last,Sort[Table[{Count[pp[[i]],0],pp[[i]]},{i,Length[pp]}]]]];
res=Table[
s=pp[[i]];
rr=fAdequateSignedMixedFast[pd,s];
rr1={rr,pp[[i]]},{i,Length[pp]}];
res=Select[res,Not[SameQ[#[[1]],{}]] &];
res] 

fAllStates[Ul_]:=Module[{pd,nn,pp,res,s,rr,rr1,i,j},
pd=If[SameQ[Head[Ul],String],fConwayToPD[Ul],Ul];
nn=Length[Table[pd[[i,j]],{i,Length[pd]},{j,4}]];
pp=fVarP[nn];
pp=Reverse[Map[Last,Sort[Table[{Count[pp[[i]],0],pp[[i]]},{i,Length[pp]}]]]];
res=Table[
s=pp[[i]];
rr=fAdequateSignedMixed[pd,s];
rr1={rr,pp[[i]]},{i,Length[pp]}];
res]

fEdgeReduction[Ulaz_List]:=Module[{Ul,vert,cyc,cyc1,cyc2,cyc3,
gr,gr1,hh,hh1,singl,vv,
res,comp,comp1,gg,mm,mm1,mm2,dod,i,j,k},
Ul=fStateGraphReducedOneDigon[Ulaz];
vert=V[Ul];
cyc=Cycle[vert];
res=If[IsomorphicQ[Ul,cyc],gr=If[SameQ[vert,2],{{1,2}},{{1,2},{1,3},{2,3}}],
Ul=ToUnorderedPairs[Ul];
gr1=Sort[Map[Sort,Ul]];
hh=Flatten[gr1];
hh1=Union[hh];
hh=Table[Count[hh,hh1[[i]]],{i,Length[hh1]}];
hh=Complement[Table[If[SameQ[hh[[i]],2],hh1[[i]],{}],{i,Length[hh]}],{{}}];
hh1=Complement[hh1,hh]; (* temena valence vece od 2 *)
dod=Map[Sort,Flatten[Table[{hh[[i]],hh1[[j]]},{i,Length[hh]},
{j,Length[hh1]}],1]];
dod=Intersection[gr1,dod];
res=If[SameQ[hh1,{}],gr1,
cyc=Select[gr1,SameQ[Intersection[#,hh1],{}] &];
comp=Select[ConnectedComponents[FromUnorderedPairs[cyc]],Length[#]>2 &];
comp1=Select[ConnectedComponents[FromUnorderedPairs[cyc]],
SameQ[Length[#],2] &]; (* veze dvovalentnih temena *)
singl=Flatten[Select[ConnectedComponents[FromUnorderedPairs[cyc]],SameQ[\
Length[#],1] &]];
cyc1=Table[vv=Flatten[Table[Sort[{hh1[[i]],comp[[k,j]]}],
{i,Length[hh1]},{j,
Length[comp[[k]]]}],1];
Select[vv,MemberQ[gr1,#]&],{k,Length[comp]}];
cyc1=Table[Table[If[MemberQ[hh1,cyc1[[j,i,1]]],cyc1[[j,i]],
Reverse[cyc1[[j,i]]
]],{i,Length[cyc1[[j]]]}],{j,Length[cyc1]}];
hh=Table[Map[Last,cyc1[[i]]],{i,Length[cyc1]}];
hh=Union[Map[Sort,Flatten[cyc1,1]],hh]; (* prve ivice *)
cyc2=Flatten[Table[vv=Flatten[Table[Sort[{hh1[[i]],comp1[[k,j]]}],
{i,Length[hh1]},{j,Length[comp1[[k]]]}],1];
Select[vv,MemberQ[gr1,#]&],{k,Length[comp1]}],1]; (* druge ivice *)
cyc3=Union[Flatten[Table[Sort[{hh1[[i]],hh1[[j]]}],{i,Length[hh1]},
{j,Length[hh1]}],1]];
cyc3=Select[cyc3,#[[1]]<#[[2]] &];
cyc3=Intersection[gr1,cyc3]; (* trece ivice *)
singl=Table[Sort[{hh1[[i]],singl[[j]]}],{i,Length[hh1]},{j,Length[singl]}];
singl=Select[Union[Flatten[singl,1]],#[[1]]<#[[2]] &];
singl=Intersection[gr1,singl];
res=Map[Sort,Union[comp1,hh,cyc2,cyc3,singl,dod]]
];
mm=Flatten[res];
mm1=Union[mm];
mm2=Range[Length[mm1]];
res=ReplaceAll[res,Table[mm1[[i]]->mm2[[i]],{i,Length[mm2]}]];
res]] 

fBasicGr[Ul_List]:=Module[{gg,ff,gr,gr1,mm,mm1,mm2,ost,
ost1,cc,dd,vv,vv1,vv2,
dod,res,ll,i},
gg=Sort[Map[Sort,Ul]]; (* sredjen ulazni graf *)
ff=fAllCycles[gg];
gr=If[SameQ[ff,{}],gg,
Union[Flatten[Table[Map[Sort,Table[{ff[[j,i]],ff[[j,i+1]]},
{i,Length[ff[[j]]]-1}]],{j,Length[ff]}],1]]]; (* glavni deo *)
gr1=gr;
mm=Flatten[gr];
mm1=Union[mm];
mm2=Range[Length[mm1]];
gr=ReplaceAll[gr,Table[mm1[[i]]->mm2[[i]],{i,Length[mm2]}]];
(* Print["glavni deo", gr]; *)
ost=Complement[gg,gr1];
ost1=ost;
cc=ConnectedComponents[FromUnorderedPairs[ost]];
cc=Select[cc,Length[#]>1 &];
ost=Table[If[SameQ[ost,{}],{},If[SameQ[Length[cc[[i]]],2],Take[cc[[i]],2],\
Take[cc[[i]],3]]],{i,Length[cc]}];
ost=Flatten[ost];
ost=If[Length[ost]>3,{Take[ost,3]},{ost}];
(* Print[ost]; *)
(* Print["ost", ost]; *) (* redukovani ostatak *)
dd=Intersection[Flatten[gr],Flatten[ost1]];
(* Print["dd",dd] *)
vv=Union[Flatten[gr]];
ll=Length[vv];
vv1=Table[ConnectedComponents[FromUnorderedPairs[Select[gr,
Not[MemberQ[#,vv[[i]]]] &]]],{i,Length[vv]}];
vv=Complement[Table[If[SameQ[Length[vv1[[i]]],1]||
Length[Complement[vv1[[i]],{
{i}}]]>1,i,0],{i,Length[vv1]}],{0}];
vv2=Flatten[Table[If[SameQ[Table[SameQ[Length[Union[Flatten[vv1[[i]]]]],ll-1],
{i,Length[vv1]}][[i]],True],{i},{}],{i,Length[vv1]}]];
vv=Complement[vv,vv2]; 
vv=Select[vv,Length[Position[gg,#]]>2 &];
vv=Union[vv,dd];
dod=Flatten[gg];
vv=Select[vv,Length[Position[dod,#]]>2 &];
res={gr,ff,ost,vv};
res]


fReducedGraph[Ul_List]:=
  Module[{poc,res,vv,vv1,vv2,gr,cc,cc1,cc2,tt,tt1,dod,dod1,dod2,dod3,aa,bb,
      ost,mm,mm1,mm2,dodatak,ss,ss1,ss2,pp,ttt,tt2,i,j},
    poc=fBasicGr[Ul];
    res=If[SameQ[Last[poc],{}]||SameQ[poc[[2]],{}],
        If[SameQ[poc[[2]],{}],
          If[SameQ[Length[poc[[1]]],1],FromUnorderedPairs[{{1,2}}],
            FromUnorderedPairs[{{1,2},{2,3}}]],
          res=FromUnorderedPairs[fEdgeReduction[fEdgeReduction[poc[[1]]]]]],
        vv=poc[[-1]];
        gr=poc[[1]];
        cc=fAllCycles[gr];
        
        
        dodatak=Select[cc,SameQ[Length[#],4] &];
        dodatak=
          Split[Sort[
              Flatten[Table[
                  Table[Sort[{dodatak[[j,i]],dodatak[[j,i+1]]}],{i,
                      Length[dodatak[[j]]]-1}],{j,Length[dodatak]}],1]]];
        dodatak=Union[Flatten[Select[dodatak,Length[#]>1 &],1]];
        
        
        tt=Union[Flatten[gr]];
        tt1=Map[Length,Table[Position[gr,tt[[i]]],{i,Length[tt]}]];
        tt=Complement[tt,vv];
        bb=Table[tt1[[tt[[i]]]],{i,Length[tt]}];
        dod=Flatten[gr];
        dod1=Union[dod];
        dod2=Table[Count[dod,dod1[[i]]],{i,Length[dod1]}];
        dod3=Select[dod2,#>2 &];
        dod3=Flatten[Table[Position[dod2,dod3[[i]]],{i,Length[dod3]}]];
        dod3=Table[dod1[[dod3[[i]]]],{i,Length[dod3]}];
        dod3=Complement[dod3,vv];
        tt=If[SameQ[dod3,{}],vv[[1]],dod3[[1]]];
        vv1=Union[vv,{tt}];
        cc=Select[cc,Not[SameQ[Intersection[#,vv1],vv1]] &];
        vv2=Complement[vv,{tt}];
        cc2=Select[cc,Not[SameQ[Intersection[#,vv2],{}]] &];
        bb=If[Length[cc2]>1,Last[cc2],{}];
        cc=Complement[cc,{bb}];
        cc=
          Union[Map[Sort,
              Flatten[Table[
                  Table[{cc[[j,i]],cc[[j,i+1]]},{i,Length[cc[[j]]]-1}],{j,
                    Length[cc]}],1]]];
        cc1=Complement[gr,cc];
        cc2=Complement[gr,cc1];
        cc2=Map[Sort,ReplaceAll[cc2,Table[vv[[i]]->tt,{i,Length[vv]}]]];
        tt=Union[cc1,cc2];
        
        tt=Union[tt,dodatak];
        
        vv1=Flatten[gr];
        vv2=Union[vv1];
        vv1=
          Flatten[Position[
              Map[Length,Table[Position[vv1,vv2[[i]]],{i,Length[vv2]}]],2]];
        vv1=Table[vv2[[vv1[[i]]]],{i,Length[vv1]}];
        vv1=
          Max[Map[Length,
              ConnectedComponents[
                FromUnorderedPairs[
                  Intersection[
                    Select[Flatten[
                        Table[Sort[{vv1[[i]],vv1[[j]]}],{i,Length[vv1]},{j,
                            Length[vv1]}],1],#[[1]]<#[[2]] &],gr]]]]];
        tt1=If[vv1<4,gr,fEdgeReduction[tt]]; 
        
        tt2=fEdgeReduction[tt];
        mm=VertexConnectivity[FromUnorderedPairs[gr]];
        tt1=If[SameQ[mm,1],tt2,tt1];
        
        
        bb=Flatten[tt1];
        aa=Union[bb];
        aa=Map[Length,Table[Position[bb,aa[[i]]],{i,Length[aa]}]];
        ttt=First[Flatten[Position[aa,Max[aa]]]];
        aa=Max[Union[Flatten[tt1]]];
        ost=poc[[3]];
        ost=ost+aa;
        
        ost=If[SameQ[ost,{{}}],{},
            ost=Table[ReplacePart[ost[[i]],ttt,1],{i,Length[ost]}];
            
            ost=Flatten[
                Table[Table[
                    Sort[{ost[[j,i]],ost[[j,i+1]]}],{i,
                      Length[ost[[1]]]-1}],{j,Length[ost]}],1]
            ];
        
        tt1=Union[tt1,ost];
        mm=Flatten[tt1];
        mm1=Union[mm];
        mm2=Range[Length[mm1]];
        tt1=ReplaceAll[tt1,Table[mm1[[i]]->mm2[[i]],{i,Length[mm2]}]];
        
        (* res=FromUnorderedPairs[tt1] *)
        
        res=FromUnorderedPairs[tt1];
        dod=Map[First,ConnectedComponents[res]];
        tt1=
          If[Length[dod]>1,
            ReplaceAll[tt1,
              Table[dod[[i+1]]->dod[[1]],{i,Length[dod]-1}]],tt1];
        mm=Flatten[tt1];
        mm1=Union[mm];
        mm2=Range[Length[mm1]];
        tt1=ReplaceAll[tt1,Table[mm1[[i]]->mm2[[i]],{i,Length[mm2]}]];
        
        
        ss=Union[Flatten[tt1]];
        ss1=
          Flatten[Position[
              Table[Length[Position[tt1,ss[[i]]]],{i,Length[ss]}],2]];
        ss2=Table[Map[First,Position[tt1,ss1[[i]]]],{i,Length[ss1]}];
        ss=Table[
            Union[Complement[tt1[[ss2[[i,1]]]],tt1[[ss2[[i,2]]]]],
              Complement[tt1[[ss2[[i,2]]]],tt1[[ss2[[i,1]]]]]],{i,
              Length[ss2]}];
        ss2=
          Flatten[Position[
              Map[Length,Table[Intersection[ss1,ss[[i]]],{i,Length[ss]}]],
              2]];
        ss2=Table[ss1[[ss2[[i]]]],{i,Length[ss2]}];
        ss=Select[tt1,SameQ[Intersection[#,ss2],{}] &];
        pp=fAllCycles[tt1];
        pp=Table[UnsortedComplement[pp[[i]],ss2],{i,Length[pp]}];
        ss1=
          Union[Flatten[
              Table[Table[
                  Sort[{pp[[j,i]],pp[[j,i+1]]}],{i,Length[pp[[j]]]-1}],{j,
                  Length[pp]}],1]];
        tt1=Union[ss,ss1];
        mm=Flatten[tt1];
        mm1=Union[mm];
        mm2=Range[Length[mm1]];
        tt1=ReplaceAll[tt1,Table[mm1[[i]]->mm2[[i]],{i,Length[mm2]}]];
        
        
        
        
        res=FromUnorderedPairs[tt1]
        ]
    ] 
    
    
fVirtPD[Con_String]:=Module[{cc,cc1,cc2,pp,pp1,pp2,i,j},
cc1=StringReplace[Con,"i"->"1"];
cc2=StringReplace[Con,"i"->"-1"];
cc1=fConwayToPD[cc1];
cc2=fConwayToPD[cc2];
pp1=Table[cc1[[i,j]],{i,Length[cc1]},{j,4}];
pp2=Table[cc2[[i,j]],{i,Length[cc1]},{j,4}];
cc=pp1-pp2;
cc=Flatten[Position[cc,{0,0,0,0}]];
If[SameQ[cc,{}],pp={},
cc1=cc;
cc=Table[pp1[[i]],{i,Length[pp1]}];
pp1=Complement[Table[If[SameQ[pp1[[i]],pp2[[i]]],{},{{pp1[[i,4]],pp1[[i,2]]},{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
pp1[[i,3]],pp1[[i,1]]}}],{i,Length[pp1]}],{{}}];
Do[cc=ReplaceAll[cc,{pp1[[i,1,1]]->pp1[[i,1,2]],pp1[[i,2,1]]->pp1[[i,2,2]]}]; \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
pp1=ReplaceAll[pp1,{pp1[[i,1,1]]->pp1[[i,1,2]],pp1[[i,2,1]]->pp1[[i,2,2]]}],{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i,Length[pp1]}];
cc=Table[cc[[cc1[[i]]]],{i,Length[cc1]}];
pp=Union[Flatten[cc]];
Do[cc=ReplaceAll[cc,pp[[i]]->i],{i,Length[pp]}];
pp=fConwayToPD[ToString[Length[cc]]];
Do[pp=ReplacePart[pp,cc[[i,j]],{i,j}],{i,Length[cc]},{j,4}]];
pp]

fFindGen[v_Integer,l_Integer,m_Integer]:=Module[{p,gen},
p=Partitions[l]; (* Print["sve: ",Length[p],p];*)
(* Print["Grading:  ",(m-1)v-l]; *)
If[l>=m,p=Rest[p]];
p=Select[p,Length[#]<=v&];
(*Print["duzina manje od v",p];*)
p=Select[p,Max[#]<=m-1&];
(* Print["Dobri ",p]; *)
p=Map[m-1-#&,p];
(*Print[Length[p],"  Pravi: ",p];*)
p=Map[Join[Table[m-1,{v-Length[#]}],#]&,p];
gen=Flatten[Map[DistinctPermutations,p],1];
(*Print[Length[gen]];*)
gen
]

graphHom[SG_List,nV_Integer,n_Integer,m_Integer,PPP_Integer]:=
Module[{p,pom={},CC,nval,i,j,C0gen={},C1gen={},mat={},op={}},
C0gen=fFindGen[nV,n,m];
(*Print[ Length[C0gen]," generators in CO: ", C0gen];*)
Do[  
      Do[ p=C0gen[[i]];
          nval=p[[SG[[j,1]]]]+p[[SG[[j,2]]]];(*Print[nval];*)
          If[nval<=m-1,(*Print["da"];   *) 
          p=ReplacePart[p,-1,SG[[j,1]]];
               p=ReplacePart[p,-1,SG[[j,2]]];
                 (*  Print[SG[[j,1]],SG[[j,2]],"!!!!!!!!!", p];*)
               pom=Prepend[Append[{p},{j,nval}],i];
               C1gen=Append[C1gen,pom]]
                ,{j,1,Length[SG]}]
 ,{i,1,Length[C0gen]}];
C1gen=Sort[C1gen];
(*Print["Here ",C1gen];*)
CC=Union[Map[Rest[#]&,C1gen]];
(*Print[C1gen];*)
(*Print["generators in C1: ",CC];*)
(*Print[Length[CC]];*)
(*make matrix*)
mat=Table[Table[0,{Length[CC]}],{Length[C0gen]}];
Do[   p=Select[C1gen,SameQ[#[[1]],i]&];
(*  Print["all from ", i,"th generator in C0 ",p];*)
  pom=Flatten[Map[Position[CC,Rest[#]]&,p]];(*Print[pom];*)
p=Table[0,{Length[CC]}];
Do[ p=ReplacePart[p,1,pom[[j]]];op=Append[op,{i,pom[[j]]}],{j,1,Length[pom]}];\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

mat=ReplacePart[mat,p,i](*;
Print[mat[[i]]]*)
,{i,1,Length[C0gen]}];
(*Print[Length[op],op];
mat*)
(* Print["Dim ",Length[C0gen]," x ",Length[CC]," No. of entries in Matrix \
",Length[op]]; *)
fDataforPari[op,nV,n,m,PPP]]


graphHom1[SG_List,nV_Integer,n_Integer,m_Integer,PPP_Integer]:=
Module[{p,pom={},CC,nval,i,j,C0gen={},C1gen={},mat={},op={}},
C0gen=fFindGen[nV,n,m];
(*Print[ Length[C0gen]," generators in CO: ", C0gen];*)
Do[  
      Do[ p=C0gen[[i]];
          nval=p[[SG[[j,1]]]]+p[[SG[[j,2]]]];(*Print[nval];*)
          If[nval<=m-1,(*Print["da"];*)
               p=ReplacePart[p,-1,SG[[j,1]]];
               p=ReplacePart[p,-1,SG[[j,2]]];
                 (*  Print[SG[[j,1]],SG[[j,2]],"!!!!!!!!!", p];*)
               pom=Prepend[Append[{p},{j,nval}],i];
               C1gen=Append[C1gen,pom]]
                ,{j,1,Length[SG]}]
 ,{i,1,Length[C0gen]}];
C1gen=Sort[C1gen];
(*Print["Here ",C1gen];*)
CC=Union[Map[Rest[#]&,C1gen]];
(*Print[C1gen];*)
(*Print["generators in C1: ",CC];*)
(*Print[Length[CC]];*)
(*make matrix*)
mat=Table[Table[0,{Length[CC]}],{Length[C0gen]}];
Do[   p=Select[C1gen,SameQ[#[[1]],i]&];
(*  Print["all from ", i,"th generator in C0 ",p];*)
  pom=Flatten[Map[Position[CC,Rest[#]]&,p]];(*Print[pom];*)
p=Table[0,{Length[CC]}];
Do[ p=ReplacePart[p,1,pom[[j]]];op=Append[op,{i,pom[[j]]}],{j,1,Length[pom]}];\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

mat=ReplacePart[mat,p,i](*;
Print[mat[[i]]]*)
,{i,1,Length[C0gen]}];
(*Print[Length[op],op];
mat*)
(*Print["Dim ",Length[C0gen]," x ",Length[CC]," No. of entries in Matrix \
",Length[op]];
fDataforPari[op,nV,n,m,PPP]*)
CC=Map[ReplacePart[#,{SG[[#[[2,1]]]],#[[2,2]]},2]&,CC];
CC]

graphHom2[SG_List,nV_Integer,n_Integer,m_Integer,PPP_Integer]:=
Module[{p,pom={},mmat,CC,nval,i,j,CC2,C2gen={},C1gen={},mat={},op={}},C1gen=\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

graphHom1[SG,nV,n,m,PPP];
(*Print["C1 generatori ",Length[C1gen],C1gen];*)
Do[  C2gen=Append[C2gen,fC2Gen[SG, C1gen[[i]],m,i]]
,{i,1,Length[C1gen]}];
If[Union[C2gen]!={{}},
C2gen=Flatten[C2gen,1];
CC2=Union[Map[Rest[#]&,C2gen]];
(*Print["C2 generators: ", CC2];*)
(*Print[C2gen[[1]],"C2 generatori ",CC2,Length[C2gen],Length[CC2]];
Print[Flatten[Position[CC2,Rest[C2gen[[1]]]]]];*)
mmat=Union[Sort[Map[{First[#],Flatten[Position[CC2,Rest[#]]][[1]]}&,C2gen]]];
mmat={Map[Abs,Select[mmat,First[#]<0&]],Select[mmat,First[#]>0&]};
(*Print[mmat[[1]]];*)
(*mmat=fMakeMap[C1gen,C2gen];*)(*Print["Matrix d1",mmat];*)
fDataforPari[mmat,nV,n,m,PPP] (* ,Print["Zero"] *)]
]

fC2Gen[GG_List, stara_List,m_Integer,od_Integer]:=Module[{SG, \
novi={},j,pom,p=stara,mapSign,nval,npom},
SG=Complement[GG,p[[2]]];
(*Print[p,"ss ",SG];*)
    Do[p=stara;
mapSign=If [Position[GG,p[[2,1]]][[1,1]]-Position[GG,SG[[j]]][[1,1]]>0,1,-1];
If[Intersection[SG[[j]],p[[2,1]]]=={},
 nval=p[[1,SG[[j,1]]]]+p[[1,SG[[j,2]]]];
          If[nval<=m-1,
               p=ReplacePart[p,-1,{1,SG[[j,1]]}];
               p=ReplacePart[p,-1,{1,SG[[j,2]]}];
                p=Append[{p[[1]]},Sort[Union[Rest[p],{{SG[[j]],nval}}]]];
                pom=Prepend[p,mapSign od];
                novi=Append[novi,pom]]];
(*slucaj kad dobijamo 1 komponentu sa dve ivice*)
(*Print[j,"ne",Intersection[SG[[j]],p[[2,1]]]];*)
If[Intersection[SG[[j]],p[[2,1]]]!={},
npom=If[p[[1,SG[[j,1]]]]==-1,2,1];
nval=p[[2,2]]+p[[1,SG[[j,npom]]]];
 If[nval<=m-1,
               p=ReplacePart[p,-1,{1,SG[[j,1]]}];
               p=ReplacePart[p,-1,{1,SG[[j,2]]}];
                pom=Prepend[ReplacePart[p,{Sort[Append[{p[[2,1]]},SG[[j]]]],\
\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

nval},2],mapSign od];
novi=Append[novi,pom]
 ]],{j,1,Length[SG]}];
(*Print["Dobijeni iz ",od," C2 generatori", novi];*)
novi
]


fMakeKn[n_Integer]:=Module[{i,j, res={}},
Do[Do[res=Append[res,{i,j}],{j,i+1,n}],{i,1,n}];
res]

fMakeMat[LLL_List,m_Integer,n_Integer]:=Module[{i,j,mat={}},
Do[p={}; Do[  If[MemberQ[LLL,{i,j}],p=Append[p,1],p=Append[p,0]],{j,1,n}] \
;mat=Append[mat,p]
,{i,1,m}];
mat
]


fMakeMap[L1_List,L2_List]:=Module[{pom,mat,i,j,op={},p},
mat=Table[Table[0,{Length[L2]}],{Length[L1]}];
Do[   p=Select[L2,SameQ[#[[1]],i]&];
(* Print["all from ", i,"th generator in C0 ",p];*)
  pom=Flatten[Map[Position[L2,#]&,p]];(*Print[pom];*)
p=Table[0,{Length[L2]}];
Do[ p=ReplacePart[p,1,pom[[j]]];op=Append[op,{i,pom[[j]]}],{j,1,Length[pom]}];\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

mat=ReplacePart[mat,p,i],{i,1,Length[L1]}];
op]

fMakeMap1[L1_List,L2_List]:=Module[{pom,mat,i,j,op={},p},
mat=Table[Table[0,{Length[L2]}],{Length[L1]}];
Do[   p=Select[L2,SameQ[#[[1]],i]&];
(* Print["all from ", i,"th generator in C0 ",p];*)
  pom=Flatten[Map[Position[L2,#]&,p]];(*Print[pom];*)
p=Table[0,{Length[L2]}];
Do[ p=ReplacePart[p,1,pom[[j]]];op=Append[op,{i,pom[[j]]}],{j,1,Length[pom]}];\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\

mat=ReplacePart[mat,p,i],{i,1,Length[L1]}];
op]

fTorsionList[Ul_]:=Module[{ss,pp,pp1,i},
ss=Ul;
ss=StringReplace[ss,{"\n"->"","  "->" "}];
pp=Union[Flatten[StringPosition[ss," "]]];
pp1=StringLength[ss];
pp=Map[ToExpression,Join[Table[StringTake[ss,{pp[[i]]+1,pp[[i+1]]-1}],{i,\
Length[pp]-1}],{StringTake[ss,{pp[[-1]]+1,pp1}]}]];
pp
] 

fDataforPari[LL_List,v_Integer,n_Integer,m_Integer,PPP_Integer]:=Module[{str=\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
"Cyc",datastr,filestr="C:\\LinKnot\\PARI\\data\\", \
filestrNames="C:\\LinKnot\\PARI\\Cyc.txt"},
str=str<>ToString[v]<>"P"<>ToString[n]<>"A"<>ToString[m]<>ToString[PPP];
filestr=filestr<>str<>".txt";
datastr=StringReplace[ToString[LL],{"{"->"[","}"->"]"}];
datastr="{"<>str<>"="<>datastr<>"; }";
filestr=OpenWrite[filestr];
WriteString[filestr,datastr];
Write[filestr];
Close[filestr];
filestrNames=OpenAppend[filestrNames];
str="domatrix("<>str<>");";
WriteString[filestrNames,str];
Write[filestrNames];
Close[filestrNames]
]

fGraphTorsion1[Ul_List,m_Integer]:=Module[{tt,pp,str,ss,res,i},
If[SameQ[Length[Union[Ul]],1],res=0,
tt={{Ul,Length[Union[Flatten[Ul]]]}};
Do[graphHom[tt[[i,1]],tt[[i,2]],2m+1,m+2,i],{i,Length[tt]}];
pp=Import["C:\\LinKnot\\PARI\\ulaznew.txt"];
pp=StringJoin[pp,"doall() = { ", Import["C:\\LinKnot\\PARI\\Cyc.txt"], \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
"}"];
Export["C:\\LinKnot\\PARI\\Cyc1.txt",pp];
Export["C:\\LinKnot\\PARI\\Cyc1",pp,"Text"];
SetDirectory["C:\\LinKnot\\PARI\\"];
str=Run["gp-dyn"," input.txt"];
ss=Import["C:\\LinKnot\\PARI\\RadmilaSNF.txt"];
res=fTorsionList[ss][[1]];
DeleteFile["C:\\LinKnot\\PARI\\Cyc.txt"];
DeleteFile["C:\\LinKnot\\PARI\\Cyc1.txt"];
DeleteFile["C:\\LinKnot\\PARI\\Cyc1"];
DeleteFile["C:\\LinKnot\\PARI\\RadmilaSNF.txt"]];
res
] 

fGraphTorsion2[Ul_List,m_Integer]:=Module[{tt,pp,str,ss,res,i},
If[SameQ[Length[Union[Ul]],1],res=0,
tt={{Ul,Length[Union[Flatten[Ul]]]}};
Do[graphHom2[tt[[i,1]],tt[[i,2]],4m-2,2m,i],{i,Length[tt]}];
pp=Import["C:\\LinKnot\\PARI\\ulaznew.txt"];
pp=StringJoin[pp,"doall() = { ", Import["C:\\LinKnot\\PARI\\Cyc.txt"], \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
"}"];
Export["C:\\LinKnot\\PARI\\Cyc1.txt",pp];
Export["C:\\LinKnot\\PARI\\Cyc1",pp,"Text"];
SetDirectory["C:\\LinKnot\\PARI\\"];
str=Run["gp-dyn"," input.txt"];
ss=Import["C:\\LinKnot\\PARI\\RadmilaSNF.txt"];
res=fTorsionList[ss][[1]];
DeleteFile["C:\\LinKnot\\PARI\\Cyc.txt"];
DeleteFile["C:\\LinKnot\\PARI\\Cyc1.txt"];
DeleteFile["C:\\LinKnot\\PARI\\Cyc1"];
DeleteFile["C:\\LinKnot\\PARI\\RadmilaSNF.txt"]];
res
] 

fKauffPolyAdqChrTor1[Con_,m_Integer]:=Module[{ll,pp,gr,gr1,res,res1,kau,i},
ll=fAllAdequateStatesFast[Con];
pp=Map[Last,ll];
res=Table[
pp=ll[[i,1]];
gr=fStateGraph[pp];
gr=Union[ToUnorderedPairs[gr]];
gr=FromUnorderedPairs[gr];
gr1=fReducedGraph[ToUnorderedPairs[gr]];
gr1,{i,Length[ll]}];
res1=Map[ToUnorderedPairs,res];
res1=Table[fGraphTorsion1[res1[[i]],m],{i,Length[res1]}];
kau=Plus@@Flatten[Table[Expand[x^res1[[i]]*ChromaticPolynomial[res[[i]],y]],{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i,Length[res1]}]];
kau] 

fKauffPolyAdqChrTor2[Con_,m_Integer]:=Module[{ll,pp,gr,gr1,res,res1,kau,i},
ll=fAllAdequateStatesFast[Con];
pp=Map[Last,ll];
res=Table[
pp=ll[[i,1]];
gr=fStateGraph[pp];
gr=Union[ToUnorderedPairs[gr]];
gr=FromUnorderedPairs[gr];
gr1=fReducedGraph[ToUnorderedPairs[gr]];
gr1,{i,Length[ll]}];
res1=Map[ToUnorderedPairs,res];
res1=Table[fGraphTorsion2[res1[[i]],m],{i,Length[res1]}];
kau=Plus@@Flatten[Table[Expand[x^res1[[i]]*ChromaticPolynomial[res[[i]],y]],{\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
i,Length[res1]}]];
kau] 


(* ## ## ## ## ## ## ## N MOVES ## ## ## ## ## ## ## ## ## ## ## *)

fMakeAllPoly[Ul_String]:=Module[{pp,ss,i},
pp=StringJoin[Ul,"."];
ss=StringPosition[pp,".1."];
ss=Flatten[Table[KSubsets[ss,i],{i,Length[ss]}],1];
ss=Table[Map[Last,ss[[i]]],{i,Length[ss]}];
ss=Table[StringInsert[pp," 0",ss[[i]]],{i,Length[ss]}];
ss=Flatten[Append[{Ul},Table[StringDrop[ss[[i]],-1],{i,Length[ss]}]]];
ss
]  




UnFix[Conway_String,k_Integer]:=Module[{Con,n,r,p,i,ddd,dd,dd1,dd2,qq},
ddd=Conway;
Con=If[SameQ[Length[StringCases[Conway,DigitCharacter]],StringLength[Conway]],\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
ToExpression[Conway],ddd];
p=If[SameQ[Head[Con],Integer],p=Con-2k,
Con=ReadList[StringToStream[StringReplace[ddd,{","->" ","+"->" ","("->" \
",")"->" ","."->" ","-"->" ","?"->" "}]],Number];
r=Length[Con];
dd=Table[StringTake[ddd,{i,i}],{i,StringLength[ddd]}];
dd=Table[StringCases[dd[[i]],DigitCharacter],{i,Length[dd]}];
dd=Complement[Table[If[SameQ[dd[[i]],{}],0,i],{i,Length[dd]}],{0}];
dd1=Map[First,StringPosition[ddd,"-"]]+1;
dd1=Flatten[Table[Position[dd,dd1[[i]]],{i,Length[dd1]}]];
p=Table[If[MemberQ[dd1,i],-Con[[i]],Con[[i]]],{i,Length[Con]}];
p];
qq=Position[p,0];
p=Table[Sign[p[[i]]]*Mod[Abs[p[[i]]],k],{i,Length[p]}];
dd=Mod[p,k];
dd2=Table[If[SameQ[dd[[i]],0],-k,dd[[i]]],{i,Length[dd]}];
dd=Table[If[SameQ[dd[[i]],0],k,dd[[i]]],{i,Length[dd]}];
dd1=dd-k;
p=Map[Sort,Table[Union[{p[[i]],dd[[i]],dd1[[i]],dd2[[i]]}],{i,Length[p]}]];
p=Table[{p[[i,1]],p[[i,2]]},{i,Length[p]}];
dd=fVarP[Length[p]];
p=Table[Table[If[SameQ[dd[[j,i]],1],p[[i,1]],p[[i,2]]],{i,Length[p]}],{j,\
Length[dd]}];
ddd=Select[p,MemberQ[#,0]&];
dd=Table[dd=Position[ddd[[j]],0];
dd1=Flatten[Table[KSubsets[dd,i],{i,Length[dd]}],1];
dd1=Table[ReplacePart[ddd[[j]],k,dd1[[i]]],{i,Length[dd1]}],{j,Length[ddd]}];
p=Union[p,Flatten[dd,1]];
p=Table[ReplacePart[p[[i]],0,qq],{i,Length[p]}];
p]


UnKLFix[Ul_String,k_Integer]:=Module[{pp,pp1,pp2,pp3,pp4,i,j},
pp=Union[Flatten[StringPosition[Ul,"*"]]];
pp1=If[SameQ[pp,{}],{},StringTake[Ul,pp[[1]]]];
pp2=If[SameQ[pp,{}],Ul,StringTake[Ul,{pp[[1]]+1,StringLength[Ul]}]];
pp=pp2;
pp2=StringJoin[Table[If[SameQ[StringCases[StringTake[pp2,{i}],DigitCharacter],\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
{}],StringTake[pp2,{i}],"1"],{i,StringLength[pp2]}]];
pp3=UnFix[pp,k];
(* Print[pp3]; *)
pp2=StringReplace[pp2,"11"->"1"];
pp2=Table[StringTake[pp2,{i,i}],{i,StringLength[pp2]}];
pp4=Flatten[Position[pp2,"1"]];
pp3=Table[Map[ToString,pp3[[i]]],{i,Length[pp3]}];
pp3=Map[StringJoin,Table[Do[pp2=ReplacePart[pp2,pp3[[j,i]],pp4[[i]]],{i,\
Length[pp4]}];
pp2,{j,Length[pp3]}]];
pp3=Table[StringJoin[pp1,pp3[[i]]],{i,Length[pp3]}];
pp3=StringReplace[pp3,{",0"->",(1,-1)","0 "->"(1,-1) \
","0,"->"(1,-1),",".0"->".(1,-1)","*0"->"*(1,-1)","--"->""}];
pp3=StringReplace[pp3,{",0"->",(1,-1)","0 "->"(1,-1) \
","0,"->"(1,-1),",".0"->".(1,-1)","*0"->"*(1,-1)"}];
(* Print[pp3]; *)
pp4=Table[If[SameQ[fComponentNo[pp3[[i]]],1],fKnotFind[fKnotscapeDow[pp3[[i]]]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
],ReductionKnotLink[fCreatePData[pp3[[i]]]]],{i,Length[pp3]}];
pp4=Table[{pp3[[i]],pp4[[i]]},{i,Length[pp3]}];
pp4
]


fAllMoves[Ul_String,k_Integer]:=Module[{ss,i},
ss=If[SameQ[Union[StringPosition[Ul,"."],StringPosition[Ul,"*"]],{}],UnKLFix[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
Ul,k],fMakeAllPoly[fMakePolyhedral[Ul]]];
ss=If[SameQ[Union[StringPosition[Ul,"."],StringPosition[Ul,"*"]],{}],ss,Table[\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
UnKLFix[ss[[i]],k],{i,Length[ss]}]];
ss
] 


(* ## ## ## ## ## ## ## ## ## ## ## # NONALTERNATING KLs ## ## ## ## # *)

fMakeNonalt[Ulaz_String]:=Module[{vv,Ul,vv1,ss,rr,ll,pp,pp1,i,j},
vv=Union[Flatten[StringPosition[Ulaz,"*"]]];
Ul=If[SameQ[vv,{}],Ul=Ulaz,
vv1=StringTake[Ulaz,{1,vv[[1]]}];
Ul=StringDrop[Ulaz,{1,vv[[1]]}]];
vv1=If[SameQ[vv,{}],{},StringTake[Ulaz,{1,vv[[1]]}]];
ss=Plus@@fCreatePData[Ulaz][[1]];
rr=fVarNewP[ss];
rr=Flatten[Table[Select[rr,SameQ[Count[#,1],i] \
&],{i,IntegerPart[Length[rr[[1]]]]/2+1}],1];
pp=Table[
ll=rr[[j]];
pp=StringPosition[Ul,"1"];
pp1=Complement[Table[If[SameQ[ll[[i]],1],pp[[i]],{}],{i,Length[ll]}],{{}}];
pp=StringJoin[vv1,StringReplacePart[Ul,"-1",pp1]],{j,Length[rr]}];
pp
] 

fAllNonaltKnotsfromAlt[Ul_String]:=Module[{ff,ff1,ff2,res,i},
ff=fMakePolyhedral[fExpandCon[Ul]];
ff=fMakeNonalt[ff];
ff=Table[{ff[[i]],fKnotFind[ff[[i]]]},{i,Length[ff]}];
ff1=Select[ff,Not[SameQ[#[[2]],1]] &];
ff2=Select[ff1,#[[2,1,1]]>=8 &];
res=Table[{naltknotlist[[First[Flatten[Position[naltknotlist,ff2[[i,2]]]]],1]]\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
,ff2[[i,1]]},{i,Length[ff2]}];
res
] 


(* ::Input:: *)
(**)


fAllNonaltLinksfromAlt[Ul_String]:=
  Module[{pd,cr,tt,tt1,tt2,pp1,pp2,rr1,pos,res,i},
    cr=Plus@@First[fCreatePData[Ul]];
    tt=fMakeAllNonalt[fMakePolyhedral[Ul]];
    tt=If[
        Not[SameQ[StringPosition[Ul,"*"],{}]]&&
          Not[SameQ[StringTake[Ul,-1],"*"]],
        Select[tt,SameQ[StringPosition[#,".-1"],{}] &],tt];
    tt1=Map[ReductionKnotLink,Map[fCreatePData,tt]];
    tt2=Complement[
        Table[If[SameQ[Plus@@tt1[[i,1]],cr],tt[[i]],{}],{i,
            Length[tt1]}],{{}}];
    pp1=Table[pd=fCreatePData[tt2[[i]]];
        rr1={Sort[{RedJones[pd],RedJones[GetMirrorImageKnot[pd]]}],
            Sort[{RedKauffman[pd],
                RedKauffman[GetMirrorImageKnot[pd]]}]},{i,
          Length[tt2]}];
    pp2=Union[pp1];
    pos=Flatten[Map[Last,Table[Position[pp1,pp2[[i]]],{i,Length[pp2]}]]];
    res=Table[tt2[[pos[[i]]]],{i,Length[pos]}];
    res
    ]


fMakeAllNonalt[Ul_String]:=
  Module[{poc,pocet,st,st1,pp,pp1,pos0,pos,rr,rr1,vv,gg,gg1,i,j},
    poc=Union[Flatten[StringPosition[Ul,"*"]]];
    pocet=If[SameQ[poc,{}],"",StringTake[Ul,poc[[1]]]];
    st=If[SameQ[poc,{}],Ul,StringTake[Ul,{poc[[1]]+1,StringLength[Ul]}]];
    st1=StringReplace[st,"0"->"a"];
    pp=StringReplace[st1," "->"0"];
    pp1=Table[StringTake[pp,{i}],{i,StringLength[pp]}];
    pos0=Flatten[Position[pp1,"0"]];
    pos=Select[Table[If[DigitQ[pp1[[i]]],i,0],{i,Length[pp1]}],
        Not[SameQ[#,0]] &];
    rr=Complement[
        Table[If[SameQ[pos[[i+1]]-pos[[i]],1],{pos[[i]],pos[[i+1]]},{}],{i,
            Length[pos]-1}],{{}}];
    rr=Select[ConnectedComponents[FromUnorderedPairs[rr]],Length[#]>1 &];
    rr1=Complement[pos,Flatten[rr]];
    rr=Union[Table[{rr1[[i]]},{i,Length[rr1]}],rr];
    rr=Complement[
        Table[Select[rr[[i]],Not[MemberQ[pos0,#]] &], {i,Length[rr]}],{{}}];
    vv=Drop[
        Select[fVarP[Length[rr]],
          Count[#,1]<=IntegerPart[Length[rr]/2] &],1];
    vv=Table[Flatten[Select[vv[[i]]*rr,Not[MemberQ[#,0]] &]],{i,Length[vv]}];
    vv=Table[gg=Table[{vv[[j,i]],vv[[j,i]]},{i,Length[vv[[j]]]}];
        gg1=Table[StringTake[st,gg[[i]]],{i,Length[gg]}];
        gg1=Table[StringJoin["-",gg1[[i]]],{i,Length[gg1]}];
        gg1=StringReplacePart[st,gg1,gg],{j,Length[vv]}];
    vv=If[SameQ[poc,{}],vv,Table[StringJoin[pocet,vv[[i]]],{i,Length[vv]}]];
    vv
    ] 


(* ## ## ## ## ## ## ## ## ## #   CONVERSIONS AND PD-DRAWING ## ## ## ## ## \
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
\
## ## ## # *)



    

    
fPDToPData[Ul_]:=Module[{pp,pp1,pp2,pd},
    pp=KnotTheory`DTCode[Ul];
    pp1=Table[pp[[i]],{i,Length[pp]}];
    pp2=Map[Length,pp1];
    pp2=If[MemberQ[pp2,0],{Length[pp2]},pp2];
    pp1=Flatten[pp1];
    pp={pp2,pp1};
    pd=fPDataFromDow[pp];
    pd]

fDrawingPD[Ul_]:=Module[{pd},
pd=fPDToPData[Ul];
ShowKnotfromPdata[pd]
] 

fDrawingPDNew[Ul_]:=Module[{pd},
    pd=fPDToPData[Ul];
    ShowKnotfromPdataNew[pd]
    ]    
    
ShowVirtKLPD[Ul_,ll_]:=Module[{zz,ll1,ss,pd1,pd2,cc1,cc2, mm, mm1, mm2, col,ww1,ww2,ww},
zz=fGenSignfromPDNew[Ul];
ll1=Table[ss=Ul[[ll[[i]]]];
If[SameQ[zz[[ll[[i]]]],1],X@@{ss[[4]],ss[[1]],ss[[2]],ss[[3]]},X@@{ss[[2]],ss[[3]],ss[[4]],ss[[1]]}],{i,Length[ll]}];
ll1=ReplacePart[Ul,Table[ll[[i]]->ll1[[i]],{i,Length[ll]}]];
pd1=fPDToPData[Ul];
pd2=fPDToPData[ll1];
ww1=Map[Sign,fDowfromPD[pd1][[2]]];
ww2=Map[Sign,fDowfromPD[pd2][[2]]];
ww=Count[Abs[ww1-ww2],2];
pd2=If[SameQ[Length[ll],ww],pd2,GetMirrorImageKnot[pd2]];
 Install["DrawKnot"];
    mm = GetDrawData[pd1[[2]], pd1[[1]]];
    mm1 = GetDrawData[pd2[[2]], pd2[[1]]];
    Uninstall["DrawKnot"];
    mm2 = mm - mm1;
    cc1 = Position[mm2, 0.8];
    cc2 = Position[mm2, -0.8];
    cc1 = Position[mm2, 0.8];
    cc2 = Position[mm2, -0.8];
    cc1 = Table[Drop[cc1[[i]], -1], {i, Length[cc1]}];
    cc2 = Table[Drop[cc2[[i]], -1], {i, Length[cc2]}];
    cc1 = 
      Flatten[Table[{{cc1[[i, 1]], cc1[[i, 2]]}, {cc2[[j, 1]], 
              cc2[[j, 2]]}}, {i, Length[cc1]}, {j, Length[cc2]}], 1];
    cc2 = 
      Flatten[Table[{{mm[[cc1[[i, 1, 1]], cc1[[i, 1, 2]]]]} - {mm[[
                  cc1[[i, 2, 1]], cc1[[i, 2, 2]]]]}}, {i, Length[cc1]}], 2];
    cc2 = Flatten[Position[cc2, {0., 0., 0.8}]];
    cc1 = Map[Flatten, Table[cc1[[cc2[[i]]]], {i, Length[cc2]}]];
    cc1 = 
      Table[{Red, 
          Sphere[Divide[
              mm[[cc1[[i, 1]], cc1[[i, 2]]]] + mm[[cc1[[i, 3]], cc1[[i, \
4]]]],
               2], 0.7]}, {i, Length[cc1]}];
    col = {Green, Blue, Cyan, Magenta, Yellow, Brown, Orange, Pink, Purple};
    Graphics3D[
      Table[Join[cc1, 
          Table[{col[[j]], 
              Cylinder[{mm[[j]][[i]], mm[[j]][[i + 1]]}, 0.2]}, {i, 
              Length[mm[[j]]] - 1}], {{col[[j]], 
              Cylinder[{mm[[j]][[Length[mm[[j]]]]], mm[[j]][[1]]}, 
                0.2]}}], {j, Length[mm]}], Boxed -> False, ImageSize -> 700]]


    
fTutte[ll_List]:=Module[{dd,dd1,ll1,ll2,str,rr,rr1,rr2,res,i,j},
dd=Length[Union[Flatten[ll]]];
dd1=Length[ll];
ll1=Flatten[Table[Table[If[SameQ[ll[[j,1]],i],1,0],{i,dd}]+Table[If[SameQ[ll[[\
j,2]],i],-1,0],{i,dd}],{j,Length[ll]}]];
ll2=StringDrop[StringJoin[Table[StringJoin[ToString[ll1[[i]]]," \
"],{i,Length[ll1]}]],-1];
Export["C:\\LinKnot\\tutte\\matrix.dat",ll2,"Text"];
SetDirectory["C:\\LinKnot\\tutte"];
str=Run["tutte",StringJoin[" ",ToString[dd]," ",ToString[dd1]," matrix.dat > \
res.txt"]];
rr=Import["res.txt"];
rr1=Last[Flatten[StringPosition[rr,"]"]]];
res=ToExpression[StringTake[rr,{rr1+1,StringLength[rr]}]];
DeleteFile["C:\\LinKnot\\tutte\\matrix.dat"];
DeleteFile["C:\\LinKnot\\tutte\\res.txt"];
SetDirectory["C:\\LinKnot"];
res
] 


ZerosAlex[Ul_String]:=Module[{mm,pp,ss,res,i},
    mm=RedAlex[fCreatePData[Ul]];
    pp=N[Solve[mm==0,x]];
    ss=Table[{Re[pp[[i,1,2]]],Im[pp[[i,1,2]]]},{i,Length[pp]}];
    res=Plus@@Table[Abs[pp[[i,1,2]]],{i,Length[pp]}];
    Print[res];
    ListPlot[ss]]

ZerosJones[Ul_String]:=Module[{mm,pp,ss,res,i},
    mm=RedJones[fCreatePData[Ul]];
    pp=N[Solve[mm==0,x]];
    ss=Table[{Re[pp[[i,1,2]]],Im[pp[[i,1,2]]]},{i,Length[pp]}];
    res=Plus@@Table[Abs[pp[[i,1,2]]],{i,Length[pp]}];
    Print[res];
    ListPlot[ss]]
    
    fCheckGap[mm_]:=Module[{mm1,mm2,i},
     mm1=Table[mm[[i]],{i,Length[mm]}];
     mm2=CoefficientList[mm,x];
     mm1=If[MemberQ[mm2,0],1,0];
     mm1]

fCheckSign[mm_]:=Module[{mm1,mm2,pp1,i},
mm1=Table[mm[[i]],{i,Length[mm]}];
mm2=Map[Sign,CoefficientList[mm,x]];
pp1=Table[(-1)^i,{i,Length[mm2]}];
mm1=If[SameQ[mm2,pp1]||SameQ[mm2,-pp1],0,1];
mm1] 

fSpecial[Ul_]:=Module[{pp1},
pp1=If[SameQ[fCheckGap[Ul],1],1,fCheckSign[Ul]];
pp1
] 

(* ## ## ## ## ## ## ## ## # KhoHO ## ## ## ## ## ## ## ## ## *)

fKhoHoHom[Ul_]:=Module[{ll,mm,mm1,str,ss,ss1,res,i,j},
    ll=If[SameQ[Head[Ul],String],Plus@@Dow[Ul][[1]],Plus@@Ul[[1]]];
    mm=If[SameQ[Head[Ul],String],fConwayToPD[Ul],fPdataToPD[Ul]];
    mm=Table[mm[[i]],{i,Length[mm]}];
    mm=Table[Table[ToString[mm[[j,i]]],{i,4}],{j,Length[mm]}];
    mm1=mm;
    mm=StringJoin["[",
        StringDrop[
          StringJoin[
            Table[StringJoin[mm[[i,1]],", ",mm[[i,2]],", ",mm[[i,3]],", ",
                mm[[i,4]],"; "],{i,Length[mm]}]],-2],"];"];
    mm=StringJoin["knot","11","[","1","]" ," = ",mm];
    mm=StringJoin["{NumberOfKnots11=1; knot11=vector(NumberOfKnots11); ",mm,
        "}"];
    SetDirectory["C:\\LinKnot\\PARI\\"];
    DeleteFile["C:\\LinKnot\\PARI\\KTable_11a"]; 
    Export["C:\\LinKnot\\PARI\\KTable_11a",mm,"Text"];
    str=Run["gp-dyn"," KH2"];
    OpenRead["C:\\LinKnot\\PARI\\res2.txt"];
    ss=Read["C:\\LinKnot\\PARI\\res2.txt"];
    Close["C:\\LinKnot\\PARI\\res2.txt"];
    DeleteFile["C:\\LinKnot\\PARI\\res2.txt"]; 
    OpenRead["C:\\LinKnot\\PARI\\res3.txt"];
    ss1=Read["C:\\LinKnot\\PARI\\res3.txt"];
    Close["C:\\LinKnot\\PARI\\res3.txt"];
    DeleteFile["C:\\LinKnot\\PARI\\res3.txt"]; 
    SetDirectory["C:\\LinKnot"];
    {Expand[ss],ss1}] 
    
    
fKhoHoOddHom[Ul_]:=Module[{ll,mm,mm1,str,ss,res,ss1,i,j},
    ll=If[SameQ[Head[Ul],String],Plus@@Dow[Ul][[1]],Plus@@Ul[[1]]];
    mm=If[SameQ[Head[Ul],String],fConwayToPD[Ul],fPdataToPD[Ul]];
    mm=Table[mm[[i]],{i,Length[mm]}];
    mm=Table[Table[ToString[mm[[j,i]]],{i,4}],{j,Length[mm]}];
    mm1=mm;
    mm=StringJoin["[",
        StringDrop[
          StringJoin[
            Table[StringJoin[mm[[i,1]],", ",mm[[i,2]],", ",mm[[i,3]],", ",
                mm[[i,4]],"; "],{i,Length[mm]}]],-2],"];"];
    mm=StringJoin["knot","11","[","1","]" ," = ",mm];
    mm=StringJoin["{NumberOfKnots11=1; knot11=vector(NumberOfKnots11); ",mm,
        "}"];
    SetDirectory["C:\\LinKnot\\PARI\\"];
    DeleteFile["C:\\LinKnot\\PARI\\KTable_11a"]; 
    Export["C:\\LinKnot\\PARI\\KTable_11a",mm,"Text"];
    str=Run["gp-dyn"," KH1"];
    OpenRead["C:\\LinKnot\\PARI\\res.txt"];
    ss=Read["C:\\LinKnot\\PARI\\res.txt"];
    Close["C:\\LinKnot\\PARI\\res.txt"];
    DeleteFile["C:\\LinKnot\\PARI\\res.txt"]; 
    OpenRead["C:\\LinKnot\\PARI\\res1.txt"];
    ss1=Read["C:\\LinKnot\\PARI\\res1.txt"];
    Close["C:\\LinKnot\\PARI\\res1.txt"];
    DeleteFile["C:\\LinKnot\\PARI\\res1.txt"]; 
    SetDirectory["C:\\LinKnot"];
    {Expand[ss],ss1}] 
    

fQuasialtTest[Ulaz_String]:=Module[{Ul,ss,ss1,ss2,vv,mm,pp,pp0,pp1,pp2,res,\
aa1,aa2,uu,i},
    ss=StringPosition[Ulaz,"*"];
    Ul=If[SameQ[ss,{}],StringReplace[Ulaz,{"1"->"(1)"}],
    ss1=StringTake[Ulaz,Union[Flatten[ss]][[1]]];
    ss2=StringTake[Ulaz,-StringLength[Ulaz]+Union[Flatten[ss]][[1]]];
    ss2=StringReplace[ss2,{"1"->"(1)"}];
    Ul=StringJoin[ss1,ss2]];
   Ul=If[SameQ[StringPosition[Ul,"*"],{}]&&Not[SameQ[StringPosition[Ul,"."],{}\
]]&&Not[SameQ[StringTake[Ul,1],"."]],StringJoin["6*",Ul],Ul];
    pp = fExpandCon[fMakePolyhedral[Ul]];
    vv = Flatten[StringPosition[pp, "*"]];
    uu = If[SameQ[vv, {}], {}, 
        StringLength[
          StringTake[pp, First[Flatten[StringPosition[pp, "*"]]]]]];
     pp0 = StringPosition[pp, "1"];
     pp1 = StringPosition[pp, "(1"];
    pp = fExpandCon[fMakePolyhedral[Ul]];
    pp = fExpandCon[fMakePolyhedral[Ul]];
    vv = Flatten[StringPosition[pp, "*"]];
    pp = If[Not[SameQ[vv, {}]],
        mm = StringTake[pp, {1, vv[[1]]}];
        pp = StringTake[pp, {vv[[1]] + 1, StringLength[pp]}], pp];
    pp0 = StringPosition[pp, "1"];
    pp1 = StringPosition[pp, "(1"];
    pp1 = 
      If[SameQ[Union[StringPosition[Ul, "*"], StringPosition[Ul, "."], 
            StringPosition[Ul, "+"]], {}],
        Table[{StringReplacePart[pp, "((1,-1)", pp1[[i]]], 
            StringReplacePart[pp, "((1,-1) (1,-1)", pp1[[i]]]}, {i, 
            Length[pp1]}],
        Table[{StringReplacePart[pp, "(1,-1)", pp0[[i]]], 
            StringReplacePart[pp, "(1,-1) (1,-1)", pp0[[i]]]}, {i, 
            Length[pp0]}]];
     aa1 = Abs[RedAlex[fCreatePData[Ulaz]] /. x -> -1];
    mm = If[SameQ[vv, {}], "", mm];
    pp1 = 
      Table[{StringJoin[mm, pp1[[i, 1]]], StringJoin[mm, pp1[[i, 2]]]}, {i, 
          Length[pp1]}];
    pp1=Partition[StringReplace[Flatten[pp1],"(1)"->"1"],2];
    aa2 = 
      Table[Plus @@ {Abs[
              If[SameQ[RedCon[fCreatePData[pp1[[i, 1]]]], 0], 1, 
                  RedAlex[fCreatePData[pp1[[i, 1]]]]] /. x -> -1], 
            Abs[If[SameQ[RedCon[fCreatePData[pp1[[i, 2]]]], 0], 1, 
                  RedAlex[fCreatePData[pp1[[i, 2]]]]] /. x -> -1]}, {i, 
          Length[pp1]}];
     aa2 = Table[Plus @@ aa2[[i]], {i, Length[aa2]}];
    res = 
      Complement[
        Table[If[SameQ[aa1, aa2[[i]]], pp1[[i]], {}], {i, 
            Length[aa2]}], {{}}];
    res = Table[{res[[i]], i}, {i, Length[res]}];
    res] 


fOddHomologyTor[Ul_]:=Module[{res,res1},
    res=If[MemberQ[Ul,{}],0,
        res=If[SameQ[fKhoHoOddHom[Ul][[2]],0],0,1];
        res1=If[fKhoHoHom[Ul][[2]]<= 2,0,1];
        res=res+res1];
    res
    ] 

fQuasiAlt[Ul_] := Module[{res1, i, res, rr},
    res = 
      If[SameQ[StringPosition[Ul, "-"], {}], 1, 
        If[fKhoHoHom[Ul][[2]] > 2 || 
            Not[SameQ[fKhoHoOddHom[Ul][[2]], 0]], 0,
           res1 = fQuasialtTest[Ul] ;
          res = If[SameQ[res1, {}], {}, i = 1;
              While[Not[SameQ[rr, 0]] && i <= Length[res1],
                
                Do[rr = 
                    Plus @@ 
                      Map[fOddHomologyTor, 
                        Map[ReductionKnotLink, 
                          Map[fCreatePData, res1[[i, 1]]]]]];
                i++];
Print[i];
              res = If[SameQ[rr, 0], 1, {}]]]];
    res
    ] 



(* MIRROR CURVES *)


checkArgs[s_,t_]:=
  If[Head[s]===List&&VectorQ[t,Head[#]===Integer&&#>=0&]&&
      Plus@@t<=Length[s],True,False]

iteratedTake[s_,t_]/;checkArgs[s,t]:=
  First/@Rest[FoldList[Through[{Take,Drop}[#1[[2]],#2]]&,{{},s},t]]


fDrawMirrorCurveCol[Ulaz_List]:=Module[{p,q,ul,ss,mm,ss1,ss2,ss3,ss4,ss5,tt1,tt2,tt3,tt4,tt5,i,j},
     p=Length[Ulaz[[1]]];
q=Length[Ulaz[[-1]]];
ul=Join[ReplaceAll[Table[Ulaz[[i]],{i,q-1}],{1->-1,-1->1}],Table[Ulaz[[i]],{i,q,p+q-2}]];
ul=ReplaceAll[ul,{1->-1,-1->1}];
mm={Table[ul[[i]],{i,q-1}],Table[ul[[i]],{i,q,Length[ul]}]};
ss1=Table[ss=mm[[1,i]];
ss1=Flatten[Position[ss,2]],{i,Length[mm[[1]]]}];
ss2=Table[
ss=mm[[1,i]];
ss2=Flatten[Position[ss,-2]],{i,Length[mm[[1]]]}];
ss3=Table[
ss=mm[[1,i]];
ss3=Flatten[Position[ss,-1]],{i,Length[mm[[1]]]}];
ss4=Table[
ss=mm[[1,i]];
ss4=Flatten[Position[ss,1]],{i,Length[mm[[1]]]}];
ss5=Table[
ss=mm[[1,i]];
ss5=Flatten[Position[ss,0]],{i,Length[mm[[1]]]}];
tt1=Table[
ss=mm[[2,i]];
tt1=Flatten[Position[ss,2]],{i,Length[mm[[2]]]}];
tt2=Table[
ss=mm[[2,i]];
tt2=Flatten[Position[ss,-2]],{i,Length[mm[[2]]]}];
tt3=Table[
ss=mm[[2,i]];
tt3=Flatten[Position[ss,-1]],{i,Length[mm[[2]]]}];
tt4=Table[
ss=mm[[2,i]];
tt4=Flatten[Position[ss,1]],{i,Length[mm[[2]]]}];
tt5=Table[
ss=mm[[2,i]];
tt5=Flatten[Position[ss,0]],{i,Length[mm[[2]]]}];
Graphics[Join[{Thick,Black,Rectangle[{0,0},{p,q}]},Flatten[Table[{Yellow,Line[{{i-1/2,j-1},{i-1,j-1/2}}]},{i,p},{j,q}]],Flatten[Table[{Yellow,Line[{{i-1/2,j-1},{i,j-1/2}}]},{i,p},{j,q}]],Flatten[Table[{Yellow,Line[{{i-1,j-1/2},{i-1/2,j}}]},{i,p},{j,q}]],Flatten[Table[{Yellow,Line[{{i-1/2,j},{i,j-1/2}}]},{i,p},{j,q}]],Flatten[Table[Table[{Green,Line[{{ss1[[j,i]]-1,j},{ss1[[j,i]],j}}]},{i,Length[ss1[[j]]]}],{j,Length[ss1]}]],Flatten[Table[Table[{Green,Line[{{ss2[[j,i]]-1/2,j-1/2},{ss2[[j,i]]-1/2,j+1/2}}]},{i,Length[ss2[[j]]]}],{j,Length[ss2]}]],Flatten[Table[Table[{Red,Line[{{ss3[[j,i]]-5/8,j-1/8},{ss3[[j,i]]-3/8,j+1/8}}]},{i,Length[ss3[[j]]]}],{j,Length[ss3]}]],Flatten[Table[Table[{Blue,Line[{{ss4[[j,i]]-3/8,j-1/8},{ss4[[j,i]]-5/8,j+1/8}}]},{i,Length[ss4[[j]]]}],{j,Length[ss4]}]],
Flatten[Table[Table[{Red,Line[{{ss5[[j,i]]-3/8,j+1/8},{ss5[[j,i]]-5/8,j-1/8}}]},{i,Length[ss5[[j]]]}],{j,Length[ss5]}]],
Flatten[Table[Table[{Red,Line[{{ss5[[j,i]]-3/8,j-1/8},{ss5[[j,i]]-5/8,j+1/8}}]},{i,Length[ss5[[j]]]}],{j,Length[ss5]}]],
Flatten[Table[Table[{Red,Circle[{ss5[[j,i]]-1/2,j},1/8]},{i,Length[ss5[[j]]]}],{j,Length[ss5]}]],Flatten[Table[Table[{Green,Line[{{j,tt1[[j,i]]-1},{j,tt1[[j,i]]}}]},{i,Length[tt1[[j]]]}],{j,Length[tt1]}]],Flatten[Table[Table[{Green,Line[{{j-1/2,tt2[[j,i]]-1/2},{j+1/2,tt2[[j,i]]-1/2}}]},{i,Length[tt2[[j]]]}],{j,Length[tt2]}]],Flatten[Table[Table[{Red,Line[{{j-1/8,tt3[[j,i]]-5/8},{j+1/8,tt3[[j,i]]-3/8}}]},{i,Length[tt3[[j]]]}],{j,Length[tt3]}]],Flatten[Table[Table[{Blue,Line[{{j-1/8,tt4[[j,i]]-3/8},{j+1/8,tt4[[j,i]]-5/8}}]},{i,Length[tt4[[j]]]}],{j,Length[tt4]}]],
Flatten[Table[Table[{Red,Line[{{j-1/8,tt5[[j,i]]-3/8},{j+1/8,tt5[[j,i]]-5/8}}]},{i,Length[tt5[[j]]]}],{j,Length[tt5]}]],
Flatten[Table[Table[{Red,Line[{{j+1/8,tt5[[j,i]]-3/8},{j-1/8,tt5[[j,i]]-5/8}}]},{i,Length[tt5[[j]]]}],{j,Length[tt5]}]],
Flatten[Table[Table[{Red,Circle[{j,tt5[[j,i]]-4/8},1/8]},{i,Length[tt5[[j]]]}],{j,Length[tt5]}]]]]]


fDrawMirrorCurve[Ulaz_List]:=Module[{p,q,ul,ss,mm,ss1,ss2,ss3,ss4,ss5,tt1,tt2,tt3,tt4,tt5,i,j},
p=Length[Ulaz[[1]]];
q=Length[Ulaz[[-1]]];
ul=Join[ReplaceAll[Table[Ulaz[[i]],{i,q-1}],{1->-1,-1->1}],Table[Ulaz[[i]],{i,q,p+q-2}]];
ul=ReplaceAll[ul,{1->-1,-1->1}];
mm={Table[ul[[i]],{i,q-1}],Table[ul[[i]],{i,q,Length[ul]}]};
ss1=Table[ss=mm[[1,i]];
ss1=Flatten[Position[ss,2]],{i,Length[mm[[1]]]}];
ss2=Table[
ss=mm[[1,i]];
ss2=Flatten[Position[ss,-2]],{i,Length[mm[[1]]]}];
ss3=Table[
ss=mm[[1,i]];
ss3=Flatten[Position[ss,-1]],{i,Length[mm[[1]]]}];
ss4=Table[
ss=mm[[1,i]];
ss4=Flatten[Position[ss,1]],{i,Length[mm[[1]]]}];
ss5=Table[
ss=mm[[1,i]];
ss5=Flatten[Position[ss,0]],{i,Length[mm[[1]]]}];
tt1=Table[
ss=mm[[2,i]];
tt1=Flatten[Position[ss,2]],{i,Length[mm[[2]]]}];
tt2=Table[
ss=mm[[2,i]];
tt2=Flatten[Position[ss,-2]],{i,Length[mm[[2]]]}];
tt3=Table[
ss=mm[[2,i]];
tt3=Flatten[Position[ss,-1]],{i,Length[mm[[2]]]}];
tt4=Table[
ss=mm[[2,i]];
tt4=Flatten[Position[ss,1]],{i,Length[mm[[2]]]}];
tt5=Table[
ss=mm[[2,i]];
tt5=Flatten[Position[ss,0]],{i,Length[mm[[2]]]}];
Graphics[Join[{Thick,White,Rectangle[{0,0},{p,q}]},Flatten[Table[{Black,Line[{{i-1/2,j-1},{i-1,j-1/2}}]},{i,p},{j,q}]],Flatten[Table[{Black,Line[{{i-1/2,j-1},{i,j-1/2}}]},{i,p},{j,q}]],Flatten[Table[{Black,Line[{{i-1,j-1/2},{i-1/2,j}}]},{i,p},{j,q}]],Flatten[Table[{Black,Line[{{i-1/2,j},{i,j-1/2}}]},{i,p},{j,q}]],Flatten[Table[Table[{Black,Line[{{ss1[[j,i]]-1,j},{ss1[[j,i]],j}}]},{i,Length[ss1[[j]]]}],{j,Length[ss1]}]],Flatten[Table[Table[{Black,Line[{{ss2[[j,i]]-1/2,j-1/2},{ss2[[j,i]]-1/2,j+1/2}}]},{i,Length[ss2[[j]]]}],{j,Length[ss2]}]],Flatten[Table[Table[{White,Line[{{ss3[[j,i]]-5/8,j-1/8},{ss3[[j,i]]-3/8,j+1/8}}]},{i,Length[ss3[[j]]]}],{j,Length[ss3]}]],
Flatten[Table[Table[{White,Line[{{ss4[[j,i]]-3/8,j-1/8},{ss4[[j,i]]-5/8,j+1/8}}]},{i,Length[ss4[[j]]]}],{j,Length[ss4]}]],
Flatten[Table[Table[{Red,Line[{{ss5[[j,i]]-3/8,j+1/8},{ss5[[j,i]]-5/8,j-1/8}}]},{i,Length[ss5[[j]]]}],{j,Length[ss5]}]],
Flatten[Table[Table[{Red,Line[{{ss5[[j,i]]-3/8,j-1/8},{ss5[[j,i]]-5/8,j+1/8}}]},{i,Length[ss5[[j]]]}],{j,Length[ss5]}]],
Flatten[Table[Table[{Red,Circle[{ss5[[j,i]]-1/2,j},1/8]},{i,Length[ss5[[j]]]}],{j,Length[ss5]}]],
Flatten[Table[Table[{Black,Line[{{j,tt1[[j,i]]-1},{j,tt1[[j,i]]}}]},{i,Length[tt1[[j]]]}],{j,Length[tt1]}]],Flatten[Table[Table[{Black,Line[{{j-1/2,tt2[[j,i]]-1/2},{j+1/2,tt2[[j,i]]-1/2}}]},{i,Length[tt2[[j]]]}],{j,Length[tt2]}]],Flatten[Table[Table[{White,Line[{{j-1/8,tt3[[j,i]]-5/8},{j+1/8,tt3[[j,i]]-3/8}}]},{i,Length[tt3[[j]]]}],{j,Length[tt3]}]],Flatten[Table[Table[{White,Line[{{j-1/8,tt4[[j,i]]-3/8},{j+1/8,tt4[[j,i]]-5/8}}]},{i,Length[tt4[[j]]]}],{j,Length[tt4]}]],
Flatten[Table[Table[{Red,Line[{{j-1/8,tt5[[j,i]]-3/8},{j+1/8,tt5[[j,i]]-5/8}}]},{i,Length[tt5[[j]]]}],{j,Length[tt5]}]],
Flatten[Table[Table[{Red,Line[{{j+1/8,tt5[[j,i]]-3/8},{j-1/8,tt5[[j,i]]-5/8}}]},{i,Length[tt5[[j]]]}],{j,Length[tt5]}]],
Flatten[Table[Table[{Red,Circle[{j,tt5[[j,i]]-4/8},1/8]},{i,Length[tt5[[j]]]}],{j,Length[tt5]}]]]]]


fAnaMirror[Ul_List]:=Module[{ss,out,res},
ss=ToString[Ul];
Export["anainput.txt",ss];
Run["MirrorToDowker.exe"];
out=Import["anaoutput.txt"];
res=ToExpression[out];
res=If[SameQ[Head[res],Integer],1,res];
res]


fAnaMirrorVirtPD[Ul_List]:=Module[{uu1,uu2,res1,res2,mm1,mm2,tt1,tt2,rr,rr1,rr2,res,i,j},
uu1=ReplaceAll[Ul,0->1];
uu2=ReplaceAll[Ul,0->-1];
res1=fAnaMirror[uu1];
res2=fAnaMirror[uu2];
mm1=fPdataToPD[res1];
mm2=fPdataToPD[res2];
tt1=Table[Table[mm1[[j,i]],{i,4}],{j,Length[mm1]}];
tt2=Table[Table[mm2[[j,i]],{i,4}],{j,Length[mm2]}];
rr=Select[Table[If[SameQ[tt1[[i]],tt2[[i]]],i,0],{i,Length[tt1]}], Not[SameQ[#,0]] &];
rr1=fPdataToPD[fPDataFromDow[fAnaMirror[uu1]]];
rr2=Table[Table[rr1[[j,i]],{i,4}],{j,Length[mm1]}];
res=PD@@Table[X@@rr2[[rr[[i]]]],{i,Length[rr]}];
res]


fGaussVirtMirrCurve[Ul_List] := Module[{cc, cc1, cc2, cc3, ou, i},
      cc =ReplaceAll[Ul, {0->1, -1 ->1}];
cc=fGaussExtSigns[fAnaMirror[cc]];
cc3 = cc;
cc1 =ReplaceAll[Ul,0-> 1];
cc1=fGaussExtSigns[fAnaMirror[cc1]];
 cc2 = ReplaceAll[Ul,0-> -1 ];
cc2=fGaussExtSigns[fAnaMirror[cc2]];
cc = Table[
        If[Not[SameQ[cc1[[i]], cc2[[i]]]], 0, cc1[[i]]], {i, Length[cc1]}];
    cc1 = Table[If[SameQ[cc3[[i]], cc1[[i]]], 1, -1], {i, Length[cc]}];
    cc1 = Table[If[SameQ[cc[[i]], 0], 0, cc1[[i]]], {i, Length[cc]}];
    ou = Table[-(-1)^i*Sign[cc1[[i]]], {i, Length[cc1]}];
    cc1 = Select[cc, Not[SameQ[#, 0]] &];
    ou = Select[ou, Not[SameQ[#, 0]] &];
    cc2 = Abs[cc1];
    cc = Union[Abs[cc1]];
    cc = Table[{cc[[i]], i}, {i, Length[cc]}];
    Do[cc2 = ReplaceAll[cc2, {cc[[i, 1]] -> cc[[i, 2]]}], {i, Length[cc]}];
    cc2 = Table[cc2[[i]]*Sign[cc1[[i]]], {i, Length[cc2]}];
    cc1 = 
      StringJoin[
        Table[StringJoin[If[SameQ[ou[[i]], -1], "U", "O"], 
            ToString[Abs[cc2[[i]]]], 
            If[SameQ[Sign[cc2[[i]]], -1], "-", "+"]], {i, Length[cc2]}]];    
cc1]



fCabledJonesVirtMirrCurve[Ul_List, n_Integer] := Module[{aa,ss},
ss=fGaussVirtMirrCurve[Ul];
aa=If[SameQ[ss,""],0,
 OpenWrite["input.txt"];
WriteString["input.txt",ss];
Close["input.txt"];
    Run[StringJoin["cabledjones.exe", " ", ToString[n], " ", 
        "<input.txt>output.txt"]];
    DeleteFile["input.txt"];
    aa = ToExpression[StringReplace[Import["output.txt"], {"\n" -> ""}]];
    DeleteFile["output.txt"];
    aa = fPolyNorm[aa]];
    aa]


fShowVirtKLMirrCurve[Ul_List] := Module[{uu1,uu2,res1,res2,mm1,mm2,tt1,tt2,rr,rr1,ll,i,j},
     uu1=ReplaceAll[Ul,0->1];
uu2=ReplaceAll[Ul,0->-1];
res1=fAnaMirror[uu1];
res2=fAnaMirror[uu2];
mm1=fPdataToPD[res1];
mm2=fPdataToPD[res2];
tt1=Table[Table[mm1[[j,i]],{i,4}],{j,Length[mm1]}];
tt2=Table[Table[mm2[[j,i]],{i,4}],{j,Length[mm2]}];
rr=Select[Table[If[SameQ[tt1[[i]],tt2[[i]]],i,0],{i,Length[tt1]}], Not[SameQ[#,0]] &];
rr1=fPdataToPD[fPDataFromDow[fAnaMirror[uu1]]];
ll=Complement[Range[Length[rr1]],rr];
ShowVirtKLPD[rr1,ll]]


fVirtPDfromPD[Ul_,ll_List] := Module[{rr,tt,tt1,tt2,cc,cc1,cc2, pp, pp1, pp2, i, j},
rr=fGenSignfromPD[Ul,ll][[1]];
tt=Table[Table[Ul[[j,i]],{i,4}],{j,Length[Ul]}];
tt1=Table[If[MemberQ[ll,i],RotateLeft[tt[[i]]],tt[[i]]],{i,Length[tt]}];
tt2=PD@@Table[X@@tt1[[i]],{i,Length[tt1]}];
 cc1 = Ul;
    cc2 = tt2;
    pp1 = Table[cc1[[i, j]], {i, Length[cc1]}, {j, 4}];
    pp2 = Table[cc2[[i, j]], {i, Length[cc1]}, {j, 4}];
    cc = pp1 - pp2;
    cc = Flatten[Position[cc, {0, 0, 0, 0}]];
    If[SameQ[cc, {}], pp = {},
      cc1 = cc;
      cc = Table[pp1[[i]], {i, Length[pp1]}];
      pp1 = 
        Complement[
          Table[If[
              SameQ[pp1[[i]], 
                pp2[[i]]], {}, {{pp1[[i, 4]], pp1[[i, 2]]}, {pp1[[i, 3]], 
                  pp1[[i, 1]]}}], {i, Length[pp1]}], {{}}];
      Do[cc = 
          ReplaceAll[
            cc, {pp1[[i, 1, 1]] -> pp1[[i, 1, 2]], 
              pp1[[i, 2, 1]] -> pp1[[i, 2, 2]]}]; 
        pp1 = ReplaceAll[
            pp1, {pp1[[i, 1, 1]] -> pp1[[i, 1, 2]], 
              pp1[[i, 2, 1]] -> pp1[[i, 2, 2]]}], {i, Length[pp1]}];
      cc = Table[cc[[cc1[[i]]]], {i, Length[cc1]}];
      pp = Union[Flatten[cc]];
      Do[cc = ReplaceAll[cc, pp[[i]] -> i], {i, Length[pp]}];
      pp = fConwayToPD[ToString[Length[cc]]];
      Do[pp = ReplacePart[pp, cc[[i, j]], {i, j}], {i, Length[cc]}, {j, 4}]];    
pp
]


fMinLinkSketcher[Ul_]:=Module[{ll,ss1,uu,dd0,pp,dd,dd1,res},
ll=Ul[[-1]];
If[SameQ[ll,{}],
ss1=Flatten[Table[Ul[[i]],{i,Length[Ul]-1}],1];
uu=PD@@Table[X@@ss1[[i]],{i,Length[ss1]}];
dd0=DTCode[uu];
pp=StringPosition[ToString[dd0],"{"];
If[SameQ[pp,{}],
dd=Table[dd0[[i]],{i,Length[dd0]}];
dd1={{Length[dd]},dd};
res=fKnotFind[dd1];
res,dd0=fKnotscapeDowfromDow[fDowfromPD[ReductionKnotLink[fPDToPData[uu]]]]],0]]


fShowMinLinkSketcher[Ul_]:=Module[{ll,hh},
ll=Ul[[-1]];
If[SameQ[ll,{}],hh=fMinLinkSketcher[Ul]; If[SameQ[hh,1],1,
ShowKnotfromPdataNew[fPDataFromDow[hh]]]]]


fDrawFromLinkSketcher[]:=Module[{Ul,ll,ss1,uu},
Ul=fVirtLinkSketch[];
ll=Ul[[-1]];
ss1=Flatten[Table[Ul[[i]],{i,Length[Ul]-1}],1];
uu=PD@@Table[X@@ss1[[i]],{i,Length[ss1]}];
ShowVirtKLPD[uu,ll]]


fVirtLSketch[]:=Module[{ss,ss1,uu,pd,res},
Run["C:/LinKnot/LinkSketcher4.0.jar"];
ss=ToExpression[Import["linkSketcherOutput.txt"]];
ss1=Flatten[Table[ss[[i]],{i,Length[ss]-1}],1];
uu=PD@@Table[X@@ss1[[i]],{i,Length[ss1]}];
pd=fVirtPDfromPD[uu,ss[[-1]]];
res={ss[[-1]],uu,pd};
res]


fVarP4[n_Integer]:=Module[{r,i},r=Table[IntegerDigits[i,4,n],{i,0,4^n-1}];
r]


fMirrSquare[n_Integer]:=Module[{ss,A,ll},
ss=fVarP4[n];
ss=ReplaceAll[ss,{0->1,1->-1,2->2,3->-2}];
A=Tuples[ss,n-1];
ll=Reverse[Union[Flatten[Table[Flatten[First[Sort[{{A[[i]],A[[j]]},{A[[j]],A[[i]]},{Map[Reverse,A[[i]]],A[[j]]},{Map[Reverse,A[[j]]],A[[i]]},{Map[Reverse,A[[i]]],Map[Reverse,A[[j]]]},{Map[Reverse,A[[j]]],Map[Reverse,A[[i]]]}}]],1],{i,Length[A]},{j,Length[A]}],1]]];
ll]


fMirrRect[p_Integer,q_Integer]:=Module[{ss,ss1,A,B,ll,i,j},
ss=fVarP4[p];
ss=ReplaceAll[ss,{0->1,1->-1,2->2,3->-2}];
ss1=fVarP4[q];
ss1=ReplaceAll[ss1,{0->1,1->-1,2->2,3->-2}];
A=Tuples[ss,q-1];
B=Tuples[ss1,p-1];
ll=Union[Flatten[Table[First[Sort[{Flatten[{A[[i]],B[[j]]},1],Flatten[{Reverse[A[[i]]],Map[Reverse,B[[j]]]},1],Flatten[{A[[i]],Map[Reverse,B[[j]]]},1],Flatten[{Reverse[A[[i]]],B[[j]]},1]}]],{i,Length[A]},{j,Length[B]}],1]];
ll=Union[Table[First[Sort[{ll[[i]],ReplaceAll[ll[[i]],{1->-1,-1->1}]}]],{i,Length[ll]}]];
ll]


fCompNo[Ul_List]:=Module[{p,q,tt,uu,tt1,gg,pp0,vv,tt2,tt3,tt4,pp,res,i,j,k},
p=Length[Ul[[1]]];
q=Length[Ul[[-1]]];
tt=Partition[Range[p*q],p];
uu=Flatten[Ul];
tt=Join[Flatten[Table[Table[{tt[[j,i]],tt[[j+1,i]]},{i,p}],{j,q-1}],1],Flatten[Table[Table[{tt[[i,j]],tt[[i,j+1]]},{i,q}],{j,p-1}],1]];
tt=Complement[Table[If[SameQ[uu[[i]],-2],tt[[i]],{}],{i,Length[uu]}],{{}}]; (* lista ivica *)
tt1=Table[{i,i},{i,p*q}];
tt=Union[tt,tt1];
gg=ConnectedComponents[FromUnorderedPairs[tt]]; (* connected components *)
pp0=Length[gg];
vv=Table[Length[gg]-Length[ConnectedComponents[FromUnorderedPairs[Drop[tt,{i,i}]]]],{i,Length[tt]}];
tt1=Complement[Table[If[SameQ[vv[[i]],0],tt[[i]],{}],{i,Length[vv]}],{{}}];
tt2=Select[ConnectedComponents[FromUnorderedPairs[tt1]],Length[#]>1 &];
tt3=Table[Select[Union[Flatten[Table[Sort[{tt2[[k,i]],tt2[[k,j]]}],{i,Length[tt2[[k]]]},{j,Length[tt2[[k]]]}],1]],Not[SameQ[#[[1]],#[[2]]]] &],{k,Length[tt2]}];
tt4=Table[Intersection[tt1,tt3[[i]]],{i,Length[tt3]}];
pp=Plus@@Table[Length[tt4[[i]]]+2-Length[Union[Flatten[tt4[[i]]]]]-1,{i,Length[tt4]}] ; (* broj pljosni *)
res=pp0+pp;
res] 


fPolyNorm[Ul_]:= Module[{tt, tt1, ll, i},
tt = If[SameQ[Ul, -1], 1, 
If[SameQ[Ul, 0], 0, 
If[SameQ[CoefficientList[Ul, x][[1]], 0], 
Abs[CoefficientList[Ul, x][[-1]]],
tt = Table[Denominator[Ul[[i]]], {i, Length[Ul]}];
tt = If[SameQ[tt, {}], Ul, tt1 = PolynomialLCM @@ tt;
tt = Expand[Ul*tt1];
ll = If[SameQ[NumberQ[tt[[1]]], True] && 
SameQ[Sign[tt[[1]]], -1] || 
Length[tt[[1]]] > 0 && SameQ[Sign[tt[[1, 1]]], -1], -tt,tt]];
tt = If[SameQ[tt, 1], tt, 
If[SameQ[Length[Variables[tt]], 1], 
Expand[Divide[tt, 
PolynomialGCD @@ Table[tt[[i]], {i, Length[tt]}]]], 
tt = ReplaceAll[tt, az -> a*z];
tt1 = PolynomialGCD @@ Table[tt[[i]], {i, Length[tt]}];
tt = Expand[Divide[tt, tt1]]]]]]];                
tt]


fBracketAlt[Ul_List]:=Module[{dd,pp,ss,ll,pp1,pp2,pp3,pp4,cc,res,i,j},
dd=Map[Length,Ul];
pp=Flatten[Ul];
ss=Flatten[Position[pp,1]];
ll=Length[ss];
pp1=fVarP[ll];
pp2=ReplaceAll[pp1,{0->-2,1->2}];
pp3=Table[Count[pp2[[i]],2],{i,Length[pp2]}];
pp4=Table[iteratedTake[ReplacePart[pp,Table[ss[[i]]->pp2[[j,i]],{i,ll}]],dd],{j,Length[pp1]}];
cc=Map[fCompNo,pp4];
res=Plus@@Flatten[Table[Expand[x^{pp3[[i]]}*x^{-ll+pp3[[i]]}*(-x^2-x^{-2})^{cc[[i]]-1}],{i,Length[cc]}]];
res
]


fBracket[Ul_List]:=Module[{dd,pp,ss,ll,pp1,pp2,pp3,pp4,ww,res,i,j},
dd=Map[Length,Ul];
pp=Flatten[Ul];
ss=Flatten[Position[pp,-1]];
ll=Length[ss]; (* broj crossings -1 *)
pp1=fVarP[ll];
pp2=ReplaceAll[pp1,{0->-2,1->2}];
pp3=Table[Count[pp2[[i]],-2],{i,Length[pp2]}];
pp4=Table[iteratedTake[ReplacePart[pp,Table[ss[[i]]->pp2[[j,i]],{i,ll}]],dd],{j,Length[pp1]}];
res=Plus@@Map[Expand,Table[x^{pp3[[i]]}*x^{-ll+pp3[[i]]}*fBracketAlt[pp4[[i]]],{i,Length[pp4]}]];
res=res[[1]];
res=ReplaceAll[res,{x->x^{-1}}];
res=If[Length[res]>0,res[[1]],res];
ww=Plus@@Select[Flatten[Ul],SameQ[Abs[#],1] &];
res=Expand[res*(-x)^{-3*ww}][[1]];
res
]


fMirrProduct[ll_List,ll1_List]:=Module[{ss,pp,pp1,tt,tt1,res,i},
ss=Map[Length,ll];
pp=Flatten[ll];
pp1=Flatten[ll1];
tt=Table[{pp[[i]],pp1[[i]]},{i,Length[pp]}];
tt1=ReplaceAll[tt,{{2,2}->2,{2,-2}->1,{2,1}->1,{2,-1}->2,{-2,2}->-1,{-2,-2}->-2,{-2,1}->-2,{-2,-1}->-1,{1,2}->2,{1,-2}->1,{1,1}->1,{1,-1}->2,{-1,2}->-1,{-1,-2}->-2,{-1,1}->-2,{-1,-1}->-1}];
res=iteratedTake[tt1,ss];
res
]


fRepresent[Ul_List]:=Module[{p,q,rr,ll,rr1,rr2,dd1,dd2,vv1,ss,tt,tt1,res,i},
p=Ul[[1]];
q=Ul[[2]];
rr=Ul[[3]];
ll=2*p*q-p-q;
rr1=IntegerDigits[Ul[[3]],2];
rr2=IntegerDigits[Ul[[4]],2];
dd1=Table[0,{i,ll-Length[rr1]}];
dd2=Table[0,{i,ll-Length[rr2]}];
vv1=ReplaceAll[{Join[dd1,rr1],Join[dd2,rr2]},{0->2,1->-2}];
ss=Flatten[{Table[p,{i,q-1}],Table[q,{i,p-1}]}];
tt=iteratedTake[vv1[[1]],ss];
tt1=iteratedTake[vv1[[2]],ss];
res={tt,tt1,fMirrProduct[tt,tt1]};
res]


fDecomposeFast[Ul_List]:=Module[{p,q,ll,pp,pp1,ss,ss1,tt,tt1,res,i},
p=Length[Ul[[1]]];
q=Length[Ul[[-1]]];
ll=Map[Length,Ul];
pp=Flatten[Ul];
pp1=ReplaceAll[pp,{2->{2,2},-2->{-2,-2},1->{2,-2},-1->{-2,2}}];
ss=Map[First,pp1];
ss1=Map[Last,pp1];
tt=iteratedTake[ss,ll];
tt1=iteratedTake[ss1,ll];
res={tt,tt1};
res
]


fDecompose[Ul_List]:=Module[{p,q,ll,pp,pp1,ss,ss1,vv,vv1,cc,cc1,tt,tt1,res,i},
p=Length[Ul[[1]]];
q=Length[Ul[[-1]]];
ll=Map[Length,Ul];
pp=Flatten[Ul];
pp1=ReplaceAll[pp,{2->{2,2},-2->{-2,-2},1->{2,-2},-1->{-2,2}}];
ss=Map[First,pp1];
ss1=Map[Last,pp1];
tt=iteratedTake[ss,ll];
tt1=iteratedTake[ss1,ll];
vv=Reverse[ReplaceAll[ss,{2->0,-2->1}]];
vv1=Reverse[ReplaceAll[ss1,{2->0,-2->1}]];
cc=Plus@@Flatten[Table[vv[[i]]*2^{i-1},{i,Length[vv]}]];
cc1=Plus@@Flatten[Table[vv1[[i]]*2^{i-1},{i,Length[vv1]}]];
res=Join[{p,q},{cc,cc1}];
res={res,tt,tt1};
res
]


fProjfromCode[Ul_List]:=Module[{p,q,ll,res,i},
p=Ul[[1]];
q=Ul[[2]];
ll=Flatten[Join[{Table[p,{i,q-1}],Table[q,{p-1}]}]];
res=iteratedTake[ReplaceAll[IntegerDigits[Ul[[3]],2,2p*q-p-q],{0->2,1->-2}],ll];
res]


fMinimize[Ul_List]:=Module[{p,q,rr,rr1,res},
p=Length[Ul[[1]]];
q=Length[Ul[[-1]]];
rr=fDecomposeFast[Ul];
rr1=Map[Flatten,ReplaceAll[rr,{2->0,-2->1}]];
res={p,q,FromDigits[rr1[[1]],2],FromDigits[rr1[[2]],2]};
res]


(* SNAPPY *)


 fDowMorwen[Ul_List] := Module[{ss,ss1,ss2,ss3,pp,ss4,res,i},
   ss=Plus@@Ul[[1]];
res=If[ss<27,
ss1=FromCharacterCode[ss+96];
ss2=FromCharacterCode[Length[Ul[[1]]]+96];
ss3=StringJoin[Map[ToString,Table[FromCharacterCode[Ul[[1,i]]+96],{i,Length[Ul[[1]]]}]]];
pp=Divide[Ul[[2]],2];
ss4=StringJoin[Map[ToString,Table[If[SameQ[Sign[pp[[i]]],1],FromCharacterCode[pp[[i]]+96],FromCharacterCode[64-pp[[i]]]],{i,Length[pp]}]]];
res=StringJoin[ss1,ss2,ss3,ss4],ToString[2]];
res]


 fMorwenDow[ll_String] := Module[{pp1,pp2,res,i},
pp1=ToCharacterCode[StringTake[ll,1]]-96;
pp2=ToCharacterCode[StringTake[ll,{2}]]-96;
pp1=ToCharacterCode[StringDrop[StringTake[ll,pp2[[1]]+2],2]]-96;
pp2=StringDrop[ll,pp2[[1]]+2];
pp2=Table[StringTake[pp2,{i}],{i,StringLength[pp2]}];
pp2=2*Flatten[Table[If[UpperCaseQ[pp2[[i]]],-ToCharacterCode[pp2[[i]]]+64,ToCharacterCode[pp2[[i]]]-96],{i,Length[pp2]}]];
res={pp1,pp2};
res
    ] 


fDTNotationToDTLet[Ul_]:=Module[{ss,tt},
If[Not[SameQ[StringTake[Ul,1],"K"]||SameQ[StringTake[Ul,1],"L"]],1,If[SameQ[StringTake[Ul,1],"K"],
ss=StringDrop[Ul,1]; ss=fKnotsc[ss][[1,1]];StringJoin[tt=FromCharacterCode[StringLength[ss]+96],"a",tt,ss],morlinks[[Flatten[Position[morlinks,Ul]][[1]]]][[2]]]]]


fClassicalNotationToDTLet[Ul_]:=Module[{ss,ss1},
ss=fConNotation[Ul];
ss1=If[SameQ[ss,{}],{},
fConToOther[fConNotation[Ul]][[2]]];
If[SameQ[ss1,{}],{},fDTNotationToDTLet[ss1]]]


fKnotscapeDT[Con_String]:=Module[{cc1,cc2,ss},
cc1=fDToD[Con];
cc2=fDToD[StringReplace[Con,"-"->""]];
ss={cc1[[1]],Map[Sign,cc1[[2]]]*Map[Sign,cc2[[2]]]*Abs[cc1[[2]]]};
ss] 


fAllToPD[Ul_]:=Module[{res},
res=If[SameQ[Head[Ul],String],fConwayToPD[Ul],If[SameQ[Head[Ul],List],
fKnotscapeDowToPD[Ul],
fKnotscapeDowToPD[fMorwenDow[Ul]]]];
res] 


 fAllToDT[Ul_]:=Module[{res},
res=If[SameQ[Head[Ul],String],fKnotscapeDT[Ul],
If[SameQ[Head[Ul],List],Ul,fMorwenDow[Ul]]];
If[SameQ[fComponentNo[res],1],fKnotFind[res],MinDowProjAltKL[res]]] 


fAllToDTLet[Ul_]:=Module[{res},
res=ToExpression[If[SameQ[Head[Ul],String],fDowMorwen[fKnotscapeDT[Ul]],
If[SameQ[Head[Ul],List],fDowMorwen[Ul],Ul]]];
res]


fCrNo[Ul_]:=Module[{res},
res=If[SameQ[Head[Ul],String],Plus@@fCreatePData[Ul][[1]],
If[SameQ[Head[Ul],List],Plus@@Ul[[1]],First[ToCharacterCode[StringTake[ToString[Ul],1]]-96]]];
res]


fSnapPyInp[Ul_]:=Module[{ss,ss1,ss2,ll,res},
res=If[SameQ[Head[Ul],String],If[SameQ[StringTake[Ul,1],"K"]||SameQ[StringTake[Ul,1],"L"],fDTNotationToDTLet[Ul],If[SameQ[StringPosition[Ul,"_"],{}],ss=fKnotscapeDow[Ul];
If[SameQ[Length[ss[[1]]],1],ss1=fKnotFind[Ul]; If[SameQ[ss1,1],ss1,If[Length[ss1[[2]]]>26,ss1,fDowMorwen[ss1]]],ss2=MinDowProjAltKL[fKnotscapeDow[Ul]];
If[MemberQ[ss2,{}],1,If[Length[ss2[[2]]]>26,ss2,fDowMorwen[ss2]]]],fClassicalNotationToDTLet[Ul]]],If[Not[SameQ[StringPosition[ToString[Ul],"{{{"],{}]],If[Not[SameQ[Ul[[-1]],{}]],1,ll=fMinLinkSketcher[Ul];
If[SameQ[Length[ll[[1]]],1],fDowMorwen[ll],ll]],If[Length[Ul]>2,ll=fMinLinkSketcher[{Ul,{}}];
If[SameQ[Length[ll[[1]]],1],fDowMorwen[ll],ll],If[SameQ[Length[Ul[[1]]],4]&&SameQ[Length[Ul[[2]]],4],1,If[SameQ[Length[Ul[[1]]],1],ss=fKnotFind[Ul];
If[SameQ[ss,1],1,fDowMorwen[Ul]],Ul]]]]];
If[SameQ[res,{}],1,res]]


fSnapPyInput[Ul_]:=Module[{res},
res=If[SameQ[Head[Ul],String],If[SameQ[Union[Table[DigitQ[StringTake[Ul,{i,i}]],{i,StringLength[Ul]}]],{False}],Ul,fSnapPyInp[Ul]],fSnapPyInp[Ul]]] 


fSnapPyUl[Ul_]:=Module[{ss,ss1,ss2,res},
ss=If[SameQ[Head[Ul],List],fKnotscapeDowToPD[Ul],Ul];
If[SameQ[Head[ss],PD],ss1=BR[ss][[2]];
ss2=StringDrop[StringJoin[Table[ToString[ss1[[i]]]<>",",{i,Length[ss1]}]],-1];
res=StringJoin["'braid[",ss2,"]'"],StringJoin["'DT[",Ul,"]'"]]] 


fForSnapPyDT[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="f=M.volume()\n";
ll3="f=str(f)\n";
ll4="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]] 


fForSnapPyBR[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll5,ll6,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import *\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="M.dehn_fill((1,0),-1)\n";
ll3="N=M.filled_triangulation()\n";
ll4="f=N.volume()\n";
ll5="f=str(f)\n";
ll6="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4,ll5,ll6];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]] 


fSnapPyVol[Ul_]:=Module[{rr,rr0},
rr0=fSnapPyInput[Ul];
rr=If[SameQ[Head[rr0],List],fForSnapPyBR[Ul],fForSnapPyDT[Ul]];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
rr=ToExpression[Import["aaa.txt"]];
rr=If[SameQ[Position[rr,e],{}]&&rr>2,rr,0]; 
rr] 


fForSnapPyAmphi[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll5,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="f=M.symmetry_group()\n";
ll3="f=f.is_amphicheiral()\n";
ll4="f=str(f)\n"; 
ll5="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4,ll5];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]]


fSnapPyAmphi[Ul_]:=Module[{rr},
rr=If[SameQ[fSnapPyVol[Ul],0],0,
fForSnapPyAmphi[Ul];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
rr=Import["aaa.txt"]];
rr]


fForSnapPyInv[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll5,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="f=M.symmetry_group()\n";
ll3="f=f.is_invertible_knot()\n";
 ll4="f=str(f)\n";
ll5="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4,ll5];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]]


fSnapPyInv[Ul_]:=Module[{rr},
rr=If[SameQ[fSnapPyVol[Ul],0],0,
fForSnapPyInv[Ul];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
rr=Import["aaa.txt"]];
rr]


fForSnapPyFundGr[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="f=M.fundamental_group()\n";
 ll3="f=str(f)\n";
ll4="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]]


fSnapPyFundGr[Ul_]:=Module[{rr},
rr=If[SameQ[fSnapPyVol[Ul],0],0,
fForSnapPyFundGr[Ul];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
rr=Import["aaa.txt"]];
rr]


fForSnapPySymmGr[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="f=M.symmetry_group()\n";
 ll3="f=str(f)\n";
ll4="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]]


fSnapPySymmGr[Ul_]:=Module[{rr},
rr=If[SameQ[fSnapPyVol[Ul],0],0,
fForSnapPySymmGr[Ul];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
rr=Import["aaa.txt"]];
rr]


fForSnapPyChSim[Ul_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul]],")","\n"];
ll2="f=M.chern_simons()\n";
 ll3="f=str(f)\n";
ll4="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]]


fSnapPyChSim[Ul_]:=Module[{rr},
rr=If[SameQ[fSnapPyVol[Ul],0],0,
fForSnapPyChSim[Ul];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
rr=Import["aaa.txt"]];
rr]


fForSnapPyIsom[Ul1_,Ul2_]:=Module[{ll0,ll1,ll2,ll3,ll4,ll5,ll,f = "C://LinKnot//inpS.py"},
ll0="from snappy import*\n";
ll1=StringJoin["M=Manifold(",fSnapPyUl[fSnapPyInput[Ul1]],")","\n"];
ll2=StringJoin["N=Manifold(",fSnapPyUl[fSnapPyInput[Ul2]],")","\n"];
ll3="f=M.is_isometric_to(N)\n";
ll4="f=str(f)\n";
ll5="open('C://LinKnot//aaa.txt','w').write(f)";
ll=StringJoin[ll0,ll1,ll2,ll3,ll4,ll5];
OpenWrite[f];
     WriteString[f,ll];
       Write[f];
    Close[f]]


fSnapPyIsom[Ul1_,Ul2_]:=Module[{res},
fForSnapPyIsom[Ul1,Ul2];
Run["C://Python26//python.exe C://LinKnot//inpS.py"];
res=Import["aaa.txt"];
res]


fPDfromLS[Ul_]:=Module[{ll,ss1,ss,ss2},
ll=Ul[[-1]];
If[SameQ[ll,{}],
ss1=Flatten[Table[Ul[[i]],{i,Length[Ul]-1}],1];
ss=PD@@Table[X@@ss1[[i]],{i,Length[ss1]}],0]]


fDTfromSnapPy[Ul_]:=Module[{rr},
rr=fSnapPyInput[Ul];
If[SameQ[Head[rr],List],rr,fMorwenDow[rr]]]


fPdatafromSnapPy[Ul_]:=Module[{rr},
rr=fSnapPyInput[Ul];
fPDataFromDowker[fSignsKL[If[SameQ[Head[rr],List],rr,fMorwenDow[rr]]]]]


(* UNKNOTTING *)


fCrossingChange[kk_]:=Module[{tt,tt1,tt2},
tt=Table[{kk[[1]],ReplacePart[kk[[2]],i->-kk[[2,i]]]},{i,Length[kk[[2]]]}];
tt1=Map[fKnotFind,tt];
tt1=Select[tt1,SameQ[#,1]||SameQ[Length[#[[1]]],1] &];
tt2=Union[tt1];
tt2]


fCrossingChangeExt[kk_]:=Module[{tt,tt1,tt2},
tt=Table[{kk[[1]],ReplacePart[kk[[2]],i->-kk[[2,i]]]},{i,Length[kk[[2]]]}];
tt1=Map[fKnotFind,tt];
tt1=Select[tt1,SameQ[#,1]||SameQ[Length[#[[1]]],1] &];
tt1]


fDowfromMor[Ul_String]:=Module[{pp,pp1,i},
pp=StringLength[Ul];
pp1=Characters[Ul];
pp={{pp},Flatten[2*Table[If[UpperCaseQ[pp1[[i]]],64-ToCharacterCode[pp1[[i]]],ToCharacterCode[pp1[[i]]]-96],{i,Length[pp1]}]]};
pp]


fUnKnotFast[Ulaz_]:=Module[{Ul,mm0,mm1,kk,kk1,kk2,kk3,kk0,mm,res,i},
Ul=If[SameQ[Head[Ulaz],String],fDowfromMor[Ulaz],Ulaz];
mm1={StringDrop[fDowMorwen[fKnotFind[Ul]],3]};
mm0=mm1[[1]];
kk=fCrossingChange[fDowfromMor[mm1[[1]]]];
kk=Select[kk,SameQ[#,1]||SameQ[#[[1,1]],Length[#[[2]]]] &];
kk0=Select[kk,Not[SameQ[#,1]] &];
kk=If[Min[Map[Length,Map[Last,kk0]]]<16,kk,{}];
res=If[SameQ[kk,{}],{},
mm1=If[MemberQ[kk,1],1,
kk1=Table[StringDrop[fDowMorwen[kk[[i]]],3],{i,Length[kk]}];
kk2=Map[First,Flatten[Complement[Table[Position[uu,kk1[[i]]],{i,Length[kk1]}],{{}}],1]];
kk3=Table[uu[[kk2[[i]]]],{i,Length[kk2]}];
mm=Min[Table[kk3[[i,2]],{i,Length[kk3]}]]+1];
{mm0,mm1}];
res
]


fNextStep[Ul_]:=Module[{ss,cc,jj,mm,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,rr,i},
ss=If[SameQ[Head[Ul],String],If[LetterQ[If[SameQ[Head[Ul],String],StringTake[Ul,1]]],Ul,StringDrop[fDowMorwen[fKnotFind[Ul]],3]],StringDrop[fDowMorwen[fKnotFind[Ul]],3]];
cc=fCrossingChangeExt[fKnotFind[fDowfromMor[ss]]];
jj=Position[cc,1];
 mm=If[Not[SameQ[jj,{}]],1,
ss1=Table[StringDrop[fDowMorwen[fKnotFind[cc[[i]]]],3],{i,Length[cc]}];
ss2=Union[ss1];
ss3=Map[First,Flatten[Complement[Table[Position[uu,ss2[[i]]],{i,Length[ss2]}],{{}}],1]];
ss4=Map[Last,Table[uu[[ss3[[i]]]],{i,Length[ss3]}]];
ss5=uu[[ss3[[First[Flatten[Position[ss4,Min[ss4]]]]]]]];
ss6=First[Flatten[Position[ss1,ss5[[1]]]]];
ss7={ss6,ss6};
ss8=StringTake[ss,ss7];
rr={StringReplacePart[ss,If[UpperCaseQ[ss8],ToLowerCase[ss8],ToUpperCase[ss8]],ss7],ss5};
{{ss,rr[[2,2]]+1},Join[{rr[[1]]},rr[[2]]]}];
If[SameQ[mm,1],ss6=First[Flatten[jj]];ss7={ss6,ss6};
ss8=StringTake[ss,ss7];
rr={{ss,1},{StringReplacePart[ss,If[UpperCaseQ[ss8],ToLowerCase[ss8],ToUpperCase[ss8]],ss7],1,0}},mm]]


fUnknotting[Ul_]:=Module[{kk,cc,cc1,kk1,kk2},
kk=If[SameQ[Head[Ul],List],fKnotFind[Ul],If[LetterQ[StringTake[Ul,1]],fKnotFind[fDowfromMor[Ul]],If[SameQ[fComponentNo[Ul],1],fKnotFind[Ul],{}]]];
If[SameQ[kk,{}],{},
cc=If[SameQ[kk,1],1,0];
cc1=If[SameQ[cc,1],1,If[Length[kk[[1]]]>1,0,2]];
cc={cc,cc1};
If[SameQ[cc,{0,2}],
kk1=If[Length[kk[[2]]]<16,uu[[First[Flatten[Position[uu,StringDrop[fDowMorwen[kk],3]]]]]],fUnKnotFast[kk]];
kk2=If[SameQ[kk1,{}],kk1,If[SameQ[Position[uuu,kk1],{}],fNextStep[kk1[[1]]],uuu[[First[Flatten[Position[uuu,kk1]]]]]]],cc[[1]]]]]


fUnknotting1[Ulaz_]:=Module[{Ul,res},
Ul=fUnknotting[Ulaz];
res=ReplaceAll[ReplaceAll[If[SameQ[Ul,1],1,If[SameQ[Ul,{}],{},{{fDowfromMor[Ul[[1,1]]],Ul[[1,2]]},{fDowfromMor[Ul[[2,1]]],fDowfromMor[Ul[[2,2]]],Ul[[2,3]]}}]],{-94->1}],{{{1},{1}}->1}];
If[SameQ[Ul,{}],{},If[SameQ[Ul,1],Ul,{Ul,res}]]]


fForUnknottingDraw[Ul_]:=Module[{rr},
rr=fUnknotting[Ul];
rr=If[SameQ[rr,1],1,If[SameQ[rr,{}],{},
If[SameQ[rr,1],{},{rr[[1,1]],rr[[1,2]],KnotbyDT[Last[fDowfromMor[rr[[1,1]]]]],rr[[2,1]],rr[[2,3]],KnotbyDT[Last[fDowfromMor[rr[[2,1]]]]],If[SameQ[rr[[2,2]],"1"],{},{rr[[2,2]],rr[[2,3]],KnotbyDT[Last[fDowfromMor[rr[[2,2]]]]]}]}]]];
If[SameQ[Last[rr],{}],Drop[rr,-1],rr]]


fUnknottingDraw[Ul_]:=Module[{rr},
rr=fUnknotting[Ul];
Flatten[If[SameQ[rr,1],1,If[SameQ[rr,{}],{},
If[SameQ[rr,1],{},{rr[[1,1]],rr[[1,2]],ShowKnotfromPdataNew[KnotbyDT[Last[fDowfromMor[rr[[1,1]]]]]],rr[[2,1]],rr[[2,3]],ShowKnotfromPdataNew[KnotbyDT[Last[fDowfromMor[rr[[2,1]]]]]],If[SameQ[rr[[2,2]],"1"],{},{rr[[2,2]],rr[[2,3]],ShowKnotfromPdataNew[KnotbyDT[Last[fDowfromMor[rr[[2,2]]]]]]}]}]]]]]


fUnRat[Ul_String]:=Module[{rr0,rr,rr1,rr2,rr3,i,rr4,rr5,res},
If[SameQ[Ul,"2"],{{"1 1",1},{"-1 1","1",0}},
rr0=UnR[Ul];
res=If[SameQ[rr0,{0}],{0},
rr=Last[rr0];
rr1=ToExpression[StringJoin["{",StringReplace[Ul,{" "->","}],"}"]];
rr2=Table[If[rr1[[i]]>0,ReplacePart[rr1,i->rr1[[i]]-2],ReplacePart[rr1,i->rr1[[i]]+2]],{i,Length[rr1]}];
rr3=Table[StringReplace[ToString[rr2[[i]]],{"{"->"","}"->"",","->""}],{i,Length[rr2]}];
i=1;
While[Last[UnR[RatReduce[rr3[[i]]]]]>rr-1,i++];
rr4=rr3[[i]];
rr5=RatReduce[rr3[[i]]];
res={{Ul,rr},{rr4,rr5,rr-1}}];
If[Length[res]>1,
{{StringReplace[res[[1,1]],"0"->"(1,-1)"],res[[1,2]]},{StringReplace[res[[2,1]],"0"->"(1,-1)"],res[[2,2]],res[[2,3]]}},res]]]


fUnRatDraw[Ul_String]:=Module[{ff},
ff=fUnRat[Ul];
If[SameQ[Length[ff],1],ShowKnotfromPdataNew[fCreatePData[Ul]],{ff[[1,1]],ShowKnotfromPdataNew[fCreatePData[ff[[1,1]]]],ff[[2,1]],ShowKnotfromPdataNew[fCreatePData[ff[[2,1]]]],ff[[2,2]],ShowKnotfromPdataNew[fCreatePData[ff[[2,2]]]]}]]


fUnSignat[Ul_]:=Module[{ff,res,rr},
ff=If[SameQ[Head[Ul],String],If[LetterQ[StringTake[Ul,1]],fDowfromMor[Ul],fKnotscapeDow[Ul]],Ul];
If[Length[ff[[1]]]>1,{},rr=fKnotFind[ff];
res=If[SameQ[rr,1],{0,0},If[Length[ff[[1]]]>1,{},Divide[Abs[fSignatureM[rr]],2]]];
If[SameQ[res,{}]||SameQ[res,{0,0}],res,Reverse[Union[{fUnknotting1[Ul][[1,1,2]],res}]]]]]


fUnLink[Ul_String]:=Module[{pp,pp1,pp2},
pp=Position[ttt,Ul];
If[SameQ[pp,{}],{},
pp1=ttt[[First[Flatten[Select[pp,SameQ[#[[2]],1] &]]]]];
pp2={{pp1[[1,1]],fMorwenDow[pp1[[1,2]]],pp1[[1,3]]},{pp1[[2,1]],fMorwenDow[pp1[[2,2]]],pp1[[2,3]]}};
{pp1,pp2}]]


fUnLink1[Ulaz_String]:=Module[{rr},
If[SameQ[Position[ttt,Ulaz],{}],{},
If[SameQ[StringPosition[Ulaz,"#"],{}],fUnLink[Ulaz],rr=fUnLink[Ulaz];
{{rr[[1,1,1]],rr[[1,1,2]],rr[[1,1,3]]},{rr[[1,2,1]],rr[[1,2,2]],rr[[1,2,3]]},rr[[2]]}]]]


fPomDrawUnLink[Ulaz_String]:=Module[{Ul,pp,qq,rr},
Ul=fUnLink[Ulaz];
pp=Ul[[1,1,4]];
qq=Ul[[1,2,4]];
rr=Ul[[2,2,1]];
{Ul[[1,1,2]],Ul[[1,1,3]],ShowKnotfromPdataNew[pp],Ul[[1,2,2]],Ul[[1,2,3]],ShowKnotfromPdataNew[qq],rr,ShowKnotfromPdataNew[fCreatePData[rr]]}]


fDrawUnLink[Ulaz_String]:=Module[{Ul,pp,pp1,qq,qq1,rr1},
If[SameQ[Position[ttt,Ulaz],{}],{},
If[Not[SameQ[StringPosition[Ulaz,"#"],{}]],fPomDrawUnLink[Ulaz],
If[SameQ[Ulaz,"2"],{"bbaaba",1,ShowKnotfromPdataNew[{{1,1},{-4,-2}}],"bbaabA",0,ShowKnotfromPdataNew[{{1,1},{4,-3}}]},
Ul=fUnLink[Ulaz];
pp=Ul[[2,1,2]];
qq=Ul[[2,2,2]];
pp1=fPDToPData[PD[DTCode@@iteratedTake[pp[[2]],pp[[1]]]]];
qq1=fPDToPData[PD[DTCode@@iteratedTake[qq[[2]],qq[[1]]]]];
rr1=fCreatePData[Ul[[2,2,1]]];
{Ul[[1,1,2]],Ul[[1,1,3]],ShowKnotfromPdataNew[pp1],Ul[[1,2,2]],Ul[[1,2,3]],ShowKnotfromPdataNew[qq1],Ul[[1,2,1]],ShowKnotfromPdataNew[rr1]}]]]]


fPomForDrawUnLink[Ulaz_String]:=Module[{Ul,pp,qq,rr},
Ul=fUnLink[Ulaz];
pp=Ul[[1,1,4]];
qq=Ul[[1,2,4]];
rr=Ul[[2,2,1]];
{Ul[[1,1,2]],Ul[[1,1,3]],pp,Ul[[1,2,2]],Ul[[1,2,3]],qq,rr,fCreatePData[rr]}]


fForDrawUnLink[Ulaz_String]:=Module[{Ul,pp,rr,pp1,qq,qq1,rr1},
If[SameQ[Position[ttt,Ulaz],{}],{},
If[Not[SameQ[StringPosition[Ulaz,"#"],{}]],fPomForDrawUnLink[Ulaz],
If[SameQ[Ulaz,"2"],{"bbaaba",1,{{1,1},{-4,-2}},"bbaabA",0,{{1,1},{4,-3}}},
Ul=fUnLink[Ulaz];
pp=Ul[[2,1,2]];
qq=Ul[[2,2,2]];
pp1=fPDToPData[PD[DTCode@@iteratedTake[pp[[2]],pp[[1]]]]];
qq1=fPDToPData[PD[DTCode@@iteratedTake[qq[[2]],qq[[1]]]]];
rr1=fCreatePData[Ul[[2,2,1]]];
rr={Ul[[1,1,2]],Ul[[1,1,3]],pp1,Ul[[1,2,2]],Ul[[1,2,3]],qq1,Ul[[1,2,1]],rr1};
If[SameQ[Ul[[1,2,1]],"1"],Take[rr,6],rr]]]]]


fLowerBound[Ul_String]:=Module[{ss,ss1,vv1,vv2,mm,mm1,mm2,mm4,mm5,mm6,vv3,vv4,res,mm3,i,j},
res=If[SameQ[Ul,"2"],1,
If[SameQ[fComponentNo[Ul],1],{},
ss=IntegerPart[Divide[fSignat[Ul]+1,2]];
ss1=fUnLink[Ul][[1,1,3]];
If[Length[StringPosition[Ul,"#"]]>0,ss1,
vv1=If[SameQ[Union[StringPosition[Ul,","],StringPosition[Ul,"*"],StringPosition[Ul,"."]],{}]&&ss1>1,Max[ss,2],ss];
vv2=If[SameQ[Union[StringPosition[Ul,"*"],StringPosition[Ul,"."],StringPosition[Ul,"("]],{}]&&Length[StringPosition[Ul,","]]>2,Max[ss,2],ss];
mm=GaussCode[fConwayToPD[Ul]];
mm1=Table[GaussCode@@mm[[i]],{i,Length[mm]}];
mm2=Table[DTCode[mm1[[i]]],{i,Length[mm1]}];
mm4=Table[mm3=Table[mm2[[j,i]],{i,Length[mm2[[j]]]}];
{{Length[mm3]},mm3},{j,Length[mm2]}];
mm5=Map[fKnotFind,mm4];
mm6=Plus@@Table[If[SameQ[mm5[[i]],1],0,fUnknotting[mm5[[i]]][[1,2]]],{i,Length[mm5]}];
vv3=LinkingNo[Ul]+mm6;
vv4=LinkingNo[Ul];
ss=Max[vv1,vv2,vv3,vv4];
If[SameQ[ss1,ss],ss1,{ss1,ss}]]]];
ReplaceAll[res,0->1]]


fKnotsc[Ul_String]:=Module[{res},
res=uuu[[First[Flatten[Position[uuuu,Ul]]]]]]


(* HYPERBOLIC *)


fHyp001[Ul_]:=Module[{res},
If[Ul<=12,ListPlot3D[hypvoltab[[Ul]]],0]]


fHyp002[Ul_]:=Module[{res},
If[Ul<=12,ListPointPlot3D[hypvoltab[[Ul]]],0]]


(* MATCH *)


fPyramid[n_Integer]:=Module[{rr,rr1,rr2,res},
rr=Range[n];
rr1=Table[{rr[[i]],n+1},{i,n}];
rr2=fCyc[rr];
res=Sort[Join[rr1,rr2]]]


fPrism[n_Integer]:=Module[{cc,cc1,ss,ss1,cc2,pr,i},
cc=ToUnorderedPairs[Cycle[n]];
cc1=cc+n;
ss=Flatten[cc];
ss1=Flatten[cc1];
cc2=Table[{ss[[i]],ss1[[i]]},{i,Length[ss]}];
pr=Union[cc,cc1,cc2];
pr]


fBipyramid[n_Integer]:=Module[{cc,cc1,cc2,bp,i},
cc=ToUnorderedPairs[Cycle[n]];
cc1=Table[{i,n+1},{i,n}];
cc2=Table[{i,n+2},{i,n}];
bp=Union[cc,cc1,cc2];
bp]


fAntiPrism[n_Integer]:=Module[{cc,cc1,cc2,ap,i},
cc=ToUnorderedPairs[Cycle[n]];
cc1=cc+n;
cc2=fCyc[Flatten[Table[{i,i+n},{i,n}]]];
ap=Union[cc,cc1,cc2];
ap]


fAntiBipyramid[n_Integer]:=Module[{cc,cc1,cc2,abp,i},
cc=fCyc[Flatten[Table[{i,i+n},{i,n}]]];
cc1=Table[{i,2n+1},{i,n}];
cc2=Table[{i+n,2n+2},{i,n}];
abp=Union[cc,cc1,cc2];
abp]


fCupola[n_Integer]:=Module[{cc,cc1,cc2,cc3,rr,i},
cc=ToUnorderedPairs[Cycle[n]];
cc1=Table[{{i,n+2*i-1},{i,n+2*i}},{i,n}];
cc2=Table[i+n,{i,2*n}];
cc3=fCyc[cc2];
rr=Union[cc,Flatten[cc1,1],cc3];
rr]


fCupolaDual[n_Integer]:=Module[{cc,cc1,cc2,cc3,cc4,rr,i},
cc=Table[i,{i,2*n}];
cc1=fCyc[cc];
cc2=Select[cc,EvenQ[#] &];
cc3=Table[{cc2[[i]],2*n+1},{i,Length[cc2]}];
cc4=Table[{cc[[i]],2*n+2},{i,Length[cc]}];
rr=Union[cc1,cc3,cc4];
rr]


fElongatedCupola[n_Integer]:=Module[{pr,ccc0,ccc1,ccc2,ccc3,ccc4,i},
pr=fPrism[2n];
ccc0=ToUnorderedPairs[Cycle[2n]]+2n;
ccc1=Table[ccc0[[2i-1]],{i,n}];
ccc2=Table[i+4n,{i,n}];
ccc3=Flatten[Table[{{ccc1[[i,1]],ccc2[[i]]},{ccc1[[i,2]],ccc2[[i]]}},{i,n}],1];
ccc4=fCyc[Range[n]]+4n;
pr=Union[pr,ccc3,ccc4];
pr]


fGyroElongatedCupola[n_Integer]:=Module[{pr,ccc0,ccc1,ccc2,ccc3,ccc4,i},
pr=fAntiPrism[2n];
ccc0=ToUnorderedPairs[Cycle[2n]]+2n;
ccc1=Table[ccc0[[2i-1]],{i,n}];
ccc2=Table[i+4n,{i,n}];
ccc3=Flatten[Table[{{ccc1[[i,1]],ccc2[[i]]},{ccc1[[i,2]],ccc2[[i]]}},{i,n}],1];
ccc4=fCyc[Range[n]]+4n;
pr=Union[pr,ccc3,ccc4];
pr]


fElongatedRotund[n_Integer]:=Module[{pr,cc,cc1,cc2,cc3,cc4,cc5,cc6,i},
pr=fPrism[2n];
cc=Table[4n+i,{i,n}];
cc1=Table[{2n+2i-1,2n+2i},{i,n}];
cc2=Table[4n+i,{i,n}];
cc3=Flatten[Table[{{cc1[[i,1]],cc2[[i]]},{cc1[[i,2]],cc2[[i]]}},{i,n}],1];
cc4=Table[5n+i,{i,n}];
cc5=fCyc[cc4];
cc6=Flatten[Table[Map[Sort,{{cc5[[i,1]],cc[[i]]},{cc5[[i,2]],cc[[i]]}}],{i,n}],1];
pr=Union[pr,cc3,cc5,cc6];
pr]


fGyroElongatedRotund[n_Integer]:=Module[{pr,cc,cc1,cc2,cc3,cc4,cc5,cc6,i},
pr=fAntiPrism[2n];
cc=Table[4n+i,{i,n}];
cc1=Table[{2n+2i-1,2n+2i},{i,n}];
cc2=Table[4n+i,{i,n}];
cc3=Flatten[Table[{{cc1[[i,1]],cc2[[i]]},{cc1[[i,2]],cc2[[i]]}},{i,n}],1];
cc4=Table[5n+i,{i,n}];
cc5=fCyc[cc4];
cc6=Flatten[Table[Map[Sort,{{cc5[[i,1]],cc[[i]]},{cc5[[i,2]],cc[[i]]}}],{i,n}],1];
pr=Union[pr,cc3,cc5,cc6];
pr]


fPrismDrum[n_Integer]:=Module[{pr,cc,cc1,cc2,cc3,cc4,cc5,cc6,i},
pr=fPrism[2n];
cc=Table[4n+i,{i,n}];
cc1=Table[{2n+2i-1,2n+2i},{i,n}];
cc2=Table[4n+i,{i,n}];
cc3=Flatten[Table[{{cc1[[i,1]],cc2[[i]]},{cc1[[i,2]],cc2[[i]]}},{i,n}],1];
cc4=Partition[Table[5n+i,{i,2n}],2];
cc5=Flatten[Table[Map[Sort,{{cc4[[i,1]],cc2[[i]]},{cc4[[i,2]],cc2[[i]]}}],{i,n}],1];
cc6=fCyc[Union[Flatten[cc4]]];
pr=Union[pr,cc3,cc5,cc6];
pr]


fGyroDrum[n_Integer]:=Module[{pr,cc,cc1,cc2,cc3,cc4,cc5,cc6,i},
pr=fAntiPrism[2n];
cc=Table[4n+i,{i,n}];
cc1=Table[{2n+2i-1,2n+2i},{i,n}];
cc2=Table[4n+i,{i,n}];
cc3=Flatten[Table[{{cc1[[i,1]],cc2[[i]]},{cc1[[i,2]],cc2[[i]]}},{i,n}],1];
cc4=Partition[Table[5n+i,{i,2n}],2];
cc5=Flatten[Table[Map[Sort,{{cc4[[i,1]],cc2[[i]]},{cc4[[i,2]],cc2[[i]]}}],{i,n}],1];
cc6=fCyc[Union[Flatten[cc4]]];
pr=Union[pr,cc3,cc5,cc6];
pr]


fOrthoBicupola[n_Integer]:=Module[{cc,cc1,cc2,cc4,cc5,cc6,rr,i},
cc=ToUnorderedPairs[Cycle[n]];
cc1=Table[{{i,n+2*i-1},{i,n+2*i}},{i,n}];
cc2=Table[3*n+i,{i,n}];
cc4=Table[Map[Sort,{{3*n+i,n+2*i-1},{3*n+i,n+2*i}}],{i,n}];
cc5=fCyc[cc2];
cc6=fCyc[Table[n+i,{i,2*n}]];
rr=Union[cc,Flatten[cc1,1],Flatten[cc4,1],cc5,cc6];
rr
]


fGyroBicupola[n_Integer]:=Module[{cc,cc1,cc2,cc4,ccc,ccc1,ccc2,cc5,cc6,rr,i},
cc=ToUnorderedPairs[Cycle[n]];
cc1=Table[{{i,n+2*i-1},{i,n+2*i}},{i,n}];
cc2=Table[3*n+i,{i,n}];
cc4=Sort[Flatten[Table[Map[Sort,{{3*n+i,n+2*i-1},{3*n+i,n+2*i}}],{i,n}] ,1]];
ccc=RotateRight[Map[First,cc4]];
ccc1=Map[Last,cc4];
ccc2=Table[{ccc[[i]],ccc1[[i]]},{i,Length[ccc]}];
cc5=fCyc[cc2];
cc6=fCyc[Table[n+i,{i,2*n}]];
rr=Union[cc,Flatten[cc1,1],cc5,cc6,ccc2];
rr]


fPrismAntiprism[n_Integer]:=Module[{pr,ccc0,ccc,ccc2,ccc1,i},
pr=fPrism[n];
ccc0=ToUnorderedPairs[Cycle[n]]+n;
ccc=ToUnorderedPairs[Cycle[n]]+2n;
ccc2=Range[n]+2n;
ccc1=Flatten[Table[{{ccc0[[i,1]],ccc2[[i]]},{ccc0[[i,2]],ccc2[[i]]}},{i,Length[ccc0]}],1];
pr=Union[pr,ccc,ccc1];
pr]


fElongatedPyramid[n_Integer,k_Integer]:=Module[{cc,tt,tt1,tt2,tt3,tt4,ep,i},
cc=ToUnorderedPairs[Cycle[n]];
tt=Table[cc+i*n,{i,0,k}];
tt1=Map[Flatten,tt];
tt2=Union[Table[Table[tt1[[i,j]],{i,Length[tt1]}],{j,2n}]];
tt3=Table[Table[{tt2[[j,i]],tt2[[j,i+1]]},{i,k}],{j,n}];
tt4=Table[Sort[{(k+1)*n+1,k*n+i}],{i,n}];
ep=Union[Flatten[tt,1],Flatten[tt3,1],tt4];
ep]


fElongatedPrism[n_Integer,k_Integer]:=Module[{cc,tt,tt1,tt2,tt3,ep,i},
cc=ToUnorderedPairs[Cycle[n]];
tt=Table[cc+i*n,{i,0,k}];
tt1=Map[Flatten,tt];
tt2=Union[Table[Table[tt1[[i,j]],{i,Length[tt1]}],{j,2n}]];
tt3=Table[Table[{tt2[[j,i]],tt2[[j,i+1]]},{i,k}],{j,n}];
ep=Union[Flatten[tt,1],Flatten[tt3,1]];
ep]


fElongatedBipyramid[n_Integer,k_Integer]:=Module[{cc,tt,tt1,tt2,tt3,tt4,tt5,ep,i},
cc=ToUnorderedPairs[Cycle[n]];
tt=Table[cc+i*n,{i,0,k}];
tt1=Map[Flatten,tt];
tt2=Union[Table[Table[tt1[[i,j]],{i,Length[tt1]}],{j,2n}]];
tt3=Table[Table[{tt2[[j,i]],tt2[[j,i+1]]},{i,k}],{j,n}];
tt4=Table[Sort[{(k+1)*n+1,k*n+i}],{i,n}];
tt5=Table[{i,(k+1)*n+2},{i,n}];
ep=Union[Flatten[tt,1],Flatten[tt3,1],tt4,tt5];
ep]


fElongatedOrthoBicupola[n_Integer,k_Integer]:=Module[{cc,tt,tt1,tt2,tt3,ep,ccc,ccc1,ccc2,ccc3,ccc4,pr,jj,jj1,jj2,jj3,jj4,i},
cc=ToUnorderedPairs[Cycle[2n]];
tt=Table[cc+2*i*n,{i,0,k}];
tt1=Map[Flatten,tt];
tt2=Union[Table[Table[tt1[[i,j]],{i,Length[tt1]}],{j,4n}]];
tt3=Table[Table[{tt2[[j,i]],tt2[[j,i+1]]},{i,k}],{j,2n}];
ep=Union[Flatten[tt,1],Flatten[tt3,1]];
ccc=ToUnorderedPairs[Cycle[2n]]+2*k*n;
ccc1=Table[2n+2k*n+i,{i,n}];
ccc2=Table[ccc[[2i-1]],{i,n}];
ccc3=Table[Map[Sort,{{ccc1[[i]],ccc2[[i,1]]},{ccc1[[i]],ccc2[[i,2]]}}],{i,n}];
ccc4=fCyc[ccc1];
pr=Union[ep,Flatten[ccc3,1],ccc4];
jj=ToUnorderedPairs[Cycle[2n]];
jj1=Table[(4k-1)*n+i,{i,n}];
jj2=Table[jj[[2i-1]],{i,n}];
jj3=Table[Map[Sort,{{jj1[[i]],jj2[[i,1]]},{jj1[[i]],jj2[[i,2]]}}],{i,n}];
jj4=fCyc[jj1];
pr=Union[pr,Flatten[jj3,1],jj4];
pr]


fElongatedGyroBicupola[n_Integer,k_Integer]:=Module[{cc,tt,tt1,tt2,tt3,ep,ccc,ccc1,ccc2,ccc3,ccc4,pr,jj,jj1,jj2,jj3,jj4,i},
cc=ToUnorderedPairs[Cycle[2n]];
tt=Table[cc+2*i*n,{i,0,k}];
tt1=Map[Flatten,tt];
tt2=Union[Table[Table[tt1[[i,j]],{i,Length[tt1]}],{j,4n}]];
tt3=Table[Table[{tt2[[j,i]],tt2[[j,i+1]]},{i,k}],{j,2n}];
ep=Union[Flatten[tt,1],Flatten[tt3,1]];
ccc=ToUnorderedPairs[Cycle[2n]]+2*k*n;
ccc1=Table[2n+2k*n+i,{i,n}];
ccc2=Table[ccc[[2i-1]],{i,n}];
ccc3=Table[Map[Sort,{{ccc1[[i]],ccc2[[i,1]]},{ccc1[[i]],ccc2[[i,2]]}}],{i,n}];
ccc4=fCyc[ccc1];
pr=Union[ep,Flatten[ccc3,1],ccc4];
jj=ToUnorderedPairs[Cycle[2n]];
jj1=Table[(2k+3)*n+i,{i,n}];
jj2=Table[jj[[2i]],{i,n}];
jj3=Table[Map[Sort,{{jj1[[i]],jj2[[i,1]]},{jj1[[i]],jj2[[i,2]]}}],{i,n}];
jj4=fCyc[jj1];
pr=Union[pr,Flatten[jj3,1],jj4];
pr]


fBicupolaAntiprism[n_Integer]:=Module[{ep,ccc,ccc1,ccc2,ccc3,ccc4,pr1,pr,i},
ep=fAntiPrism[2n];
ccc=ToUnorderedPairs[Cycle[2n]]+2n;
ccc1=Table[2n+2n+i,{i,n}];
ccc2=Table[ccc[[2i-1]],{i,n}];
ccc3=Table[Map[Sort,{{ccc1[[i]],ccc2[[i,1]]},{ccc1[[i]],ccc2[[i,2]]}}],{i,n}];
ccc4=fCyc[ccc1];
pr1=Union[ep,Flatten[ccc3,1],ccc4];
ccc=ToUnorderedPairs[Cycle[2n]];
ccc1=Table[5n+i,{i,n}];
ccc2=Table[ccc[[2i-1]],{i,n}];
ccc3=Table[Map[Sort,{{ccc1[[i]],ccc2[[i,1]]},{ccc1[[i]],ccc2[[i,2]]}}],{i,n}];
ccc4=fCyc[ccc1];
pr=Union[pr1,Flatten[ccc3,1],ccc4];
pr]


fTruncation[Ulaz_List]:=Module[{uu,ss,vv,vv1,vv2,tt,tt1,ff,ff1,ff2,rr,rr1,rr2,res,i,j},
uu=Length[Union[Flatten[Ulaz]]];
ss=Map[fCyc,Ulaz];
vv=Map[Flatten,ss];
vv1=Flatten[Table[Table[{vv[[j,i]],j},{i,Length[vv[[j]]]}],{j,Length[vv]}],1];
vv2=Partition[vv1,2];(* ivice unutar pljosni *)
tt=Union[Map[Sort,Flatten[Table[{i,j},{i,Length[Ulaz]},{j,Length[Ulaz]}],1]]];
tt=Select[tt,#[[1]]<#[[2]] &];
tt1=Table[Intersection[ss[[tt[[i,1]]]],ss[[tt[[i,2]]]]],{i,Length[tt]}];
ff=Union[Flatten[tt1,1]];
ff1=Table[Position[tt1,ff[[i]]],{i,Length[ff]}];
ff2=Complement[Table[If[Length[ff1[[i]]]>2,ff1[[i,-1]],{}],{i,Length[ff1]}],{{}}];
Do[tt1=ReplacePart[tt1,ff2[[i]]->{}],{i,Length[ff2]}];
tt1=Table[Complement[tt1[[i]],{{}}],{i,Length[tt1]}];
tt1=Flatten[Table[If[SameQ[tt1[[j]],{}],{},Flatten[Table[{{{tt1[[j,i,1]],tt[[j,1]]},{tt1[[j,i,1]],tt[[j,2]]}},{{tt1[[j,i,2]],tt[[j,1]]},{tt1[[j,i,2]],tt[[j,2]]}}},{i,Length[tt1[[j]]]}],1]],{j,Length[tt1]}],1];
rr=Join[vv2,tt1];
rr1=Flatten[rr,1];
rr2=Union[rr1];
res=Sort[Partition[Flatten[Table[Position[rr2,rr1[[i]]],{i,Length[rr1]}]],2]];
res
]


fCyc[Ul_List]:=Module[{res,i},
res=Map[Sort,Join[{{Ul[[1]],Ul[[-1]]}},Table[{Ul[[i]],Ul[[i+1]]},{i,Length[Ul]-1}]]];
res
]


fPolyPlot[Ul_]:=Module[{ggg},
ggg=Table[Ul[[i,1]]->Ul[[i,2]],{i,Length[Ul]}];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False, ImageSize->500]]


fPolyPlotDoubleEdges[Ulaz_]:=Module[{ggg,ss,Ul},
ss=Select[Ulaz,Length[Position[Ulaz,#]]>1 &];
Ul=Union[Ulaz];
ggg=Table[Ul[[i,1]]->Ul[[i,2]],{i,Length[Ul]}];
GraphPlot3D[ggg,EdgeRenderingFunction->(If [MemberQ[ss,#2],{Blue,Cylinder[#1,0.05]},{Cylinder[#1,0.05]}] &),VertexRenderingFunction->({Sphere[#,0.15]} &),PlotStyle->Directive[Specularity[White,20]],Boxed->False, ImageSize->700]]


fPolyPlotKL[Ulaz_]:=Module[{ggg,ss0,ss1,ss2,uu},
ss0=fGraphInc[Ulaz];
ss1=ss0[[1]];
ss2=ss0[[2]];
ss2=Flatten[Table[If[SameQ[ss2[[i]],1],i,{}],{i,Length[ss2]}]];
ss=Select[ss1,Length[Position[ss1,#]]>1 &];
uu=Union[ss1];
ggg=Table[uu[[i,1]]->uu[[i,2]],{i,Length[uu]}];
GraphPlot3D[ggg,EdgeRenderingFunction->(If [MemberQ[ss,#2],{Blue,Cylinder[#1,0.05]},{Cylinder[#1,0.05]}] &),VertexRenderingFunction->(If[MemberQ[ss2,#2],{Green,Sphere[#,0.15]} ,{Red,Sphere[#,0.15]}] &),PlotStyle->Directive[Specularity[White,20]],Boxed->False, ImageSize->700]]


fPolyFunctions1[Ul_,n_Integer]:=Module[{res},
res=Switch[Ul,"Pyramid", fPyramid[n],"Prism", fPrism[n],"Bipyramid",fBipyramid[n],"Antiprism",fAntiPrism[n] ,"Antibipyramid",fAntiBipyramid[n],"Cupola",fCupola[n],"CupolaDual",fCupolaDual[n],"ElongatedCupola",fElongatedCupola[n],"GyroElongatedCupola",fGyroElongatedCupola[n],"ElongatedRotund",fElongatedRotund[n],"GyroElongatedRotund",fGyroElongatedRotund[n],"PrismDrum",fPrismDrum[n],"GyroDrum",fGyroDrum[n],"OrthoBicupola",fOrthoBicupola[n],"GyroBicupola",fGyroBicupola[n],"PrismAntiprism",fPrismAntiprism[n]]]


fPolyFunctions2[Ul_,n_Integer,k_Integer]:=Module[{res},
res=Switch[Ul,"ElongatedPrism", fElongatedPrism[n,k],"ElongatedBipyramid",fElongatedBipyramid[n,k],"ElongatedPyramid",fElongatedPyramid[n,k] ,"ElongatedOrthoBicupola",fElongatedOrthoBicupola[n,k],"ElongatedGyroBicupola",fElongatedGyroBicupola[n,k]]]


fPolyFunctions4[Ul_,n_Integer]:=Module[{res},
res=Switch[Ul,"Antiprism", fAntiPrism[n],"OrthoBicupola",fOrthoBicupola[n],"GyroBicupola",fGyroBicupola[n]]]


fPolyFunctions41[Ul_,n_Integer,k_Integer]:=Module[{res},
res=Switch[Ul,"ElongatedOrthoBicupola", fElongatedOrthoBicupola[n,k],"ElongatedGyroBicupola",fElongatedGyroBicupola[n,k]]]


fJaeger[Ul_List]:=Module[{ll,ll1,gg,gg1,dt,i},
ll=Length[Ul];
ll1=Length[Union[Flatten[Ul]]]+Range[ll];
gg1=Flatten[Table[{{Ul[[i,1]],ll1[[i]]},Sort[{Ul[[i,2]],ll1[[i]]}]},{i,Length[ll1]}],1];
dt=fKLfromGraph[fMidEdgeGraph[gg1]];
dt]


fMidEdgeDT1[Ul_,n_Integer]:=Module[{pr,jj},
pr=fPolyFunctions1[Ul,n];
jj=fMidEdgeGraph[pr];
fKLfromGraph[jj]]


fMidEdgeKL1[Ul_,n_Integer]:=Module[{dd},
dd=fMidEdgeDT1[Ul,n];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fMidEdgePlot1[Ul_,n_Integer]:=Module[{pr,jj},
pr=fPolyFunctions1[Ul,n];
jj=fMidEdgeGraph[pr];
fPolyPlot[jj]]


fMidEdgeGraphPlot1[Ul_,n_Integer]:=Module[{pr,jj},
pr=fPolyFunctions1[Ul,n];
jj=fMidEdgeGraph[pr];
DrawPlanarEmbGraph[jj]]


fMidEdgeDT2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj},
pr=fPolyFunctions2[Ul,n,k];
jj=fMidEdgeGraph[pr];
fKLfromGraph[jj]]


fMidEdgeKL2[Ul_,n_Integer,k_Integer]:=Module[{dd},
dd=fMidEdgeDT2[Ul,n,k];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fMidEdgePlot2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj},
pr=fPolyFunctions2[Ul,n,k];
jj=fMidEdgeGraph[pr];
fPolyPlot[jj]]


fMidEdgeGraphPlot2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj},
pr=fPolyFunctions2[Ul,n,k];
jj=fMidEdgeGraph[pr];
DrawPlanarEmbGraph[jj]]


fTruncDT1[Ul_,n_Integer]:=Module[{pr,jj},
pr=fPolyFunctions1[Ul,n];
jj=fKLfromGraph[fGraphInc[fTruncDT[pr]][[1]]]]


fTruncKL1[Ul_,n_Integer]:=Module[{dd},
dd=fTruncDT1[Ul,n];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fTruncPlot1[Ul_,n_Integer]:=Module[{pr,jj},
pr=fPolyFunctions1[Ul,n];
jj=fGraphInc[fTruncDT[pr]][[1]];
fPolyPlot[jj]]


fTruncGraphPlot1[Ul_,n_Integer]:=Module[{pr,jj},
pr=fPolyFunctions1[Ul,n];
jj=fGraphInc[fTruncDT[pr]][[1]];
DrawPlanarEmbGraph[jj]]


fTruncDT2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj},
pr=fPolyFunctions2[Ul,n,k];
jj=fKLfromGraph[fGraphInc[fTruncDT[pr]][[1]]]]


fTruncKL2[Ul_,n_Integer,k_Integer]:=Module[{dd},
dd=fTruncDT2[Ul,n,k];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fTruncPlot2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj},
pr=fPolyFunctions2[Ul,n,k];
jj=fPolyPlot[fGraphInc[fTruncDT[pr]][[1]]]]


fTruncGraphPlot2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj},
pr=fPolyFunctions2[Ul,n,k];
jj=DrawPlanarEmbGraph[fGraphInc[fTruncDT[pr]][[1]]]]


fMidEdgeDT0[Ul_]:=Module[{pr,jj},
pr=PolyhedronData[Ul,"Edges"][[2,1]];
jj=fMidEdgeGraph[pr];
fKLfromGraph[jj]]


fMidEdgeKL0[Ul_]:=Module[{dd},
dd=fMidEdgeDT0[Ul];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fMidEdgePlot0[Ul_]:=Module[{pr,jj},
pr=PolyhedronData[Ul,"Edges"][[2,1]];
jj=fMidEdgeGraph[pr];
fPolyPlot[jj]]


fMidEdgeGraphPlot0[Ul_]:=Module[{pr,jj},
pr=PolyhedronData[Ul,"Edges"][[2,1]];
jj=fMidEdgeGraph[pr];
DrawPlanarEmbGraph[jj]]


fTruncDT0[Ul_]:=Module[{pr,jj},
pr=PolyhedronData[Ul,"Edges"][[2,1]];
jj=fKLfromGraph[fGraphInc[fTruncDT[pr]][[1]]]]


fTruncKL0[Ul_]:=Module[{dd},
dd=fTruncDT0[Ul];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fTruncPlot0[Ul_]:=Module[{pr,jj},
pr=PolyhedronData[Ul,"Edges"][[2,1]];
jj=fGraphInc[fTruncDT[pr]][[1]];
fPolyPlot[jj]]


fTruncGraphPlot0[Ul_]:=Module[{pr,jj},
pr=PolyhedronData[Ul,"Edges"][[2,1]];
jj=fGraphInc[fTruncDT[pr]][[1]];
DrawPlanarEmbGraph[jj]]


fPolyTransf[Ul_,tt_Integer]:=Module[{ss},
ss=Switch[tt,0,PolyhedronData[Ul],1,Truncate[PolyhedronData[Ul]],2,Stellate[PolyhedronData[Ul]],3,Geodesate[PolyhedronData[Ul]]]]


fPolyTransfJaeger[Ul_,tt_Integer]:=Module[{ss},
ss=Switch[tt,0,Truncate[PolyhedronData[Ul]],1,Truncate[Truncate[PolyhedronData[Ul]]],2,Truncate[Stellate[PolyhedronData[Ul]]],3,Truncate[Geodesate[PolyhedronData[Ul]]]]]


fPolyTransfInput[Ul_,tt_Integer]:=Module[{vv,vv1,vv2},
vv=fPolyTransf[Ul,tt][[1,2]];
vv=Union[Map[Sort,Flatten[Map[fCyc,If[SameQ[tt,1],Flatten[Table[vv[[i,1]],{i,Length[vv]}],1],vv[[1]]]],1]]];
vv1=Union[Flatten[vv]];
vv2=Length[vv1];
vv=ReplaceAll[vv,Table[vv1[[i]]->i,{i,vv2}]]]


fPolyTransfMidInput[Ul_,tt_Integer]:=Module[{vv,vv1,vv2,vv3},
vv=fPolyTransf[Ul,tt][[1,2]];
vv=If[SameQ[tt,1],Flatten[Table[vv[[i,1]],{i,Length[vv]}],1],vv[[1]]];
vv1=Map[fCyc,vv];
vv2=Union[Flatten[vv1,1]];
vv3=Sort[Flatten[Map[fCyc,ReplaceAll[vv1,Table[vv2[[i]]->i,{i,Length[vv2]}]]],1]]]


fPolyTransfTruncInput[Ul_,tt_Integer]:=Module[{vv,vv1},
vv=fPolyTransf[Ul,tt][[1,2]];
vv=If[SameQ[tt,1],Flatten[Table[vv[[i,1]],{i,Length[vv]}],1],vv[[1]]];
vv1=fTruncation[vv]]


fPolyTransfJaegerInput[Ul_,tt_Integer]:=Module[{vv,vv1,vv2,nn},
vv=fPolyTransfJaeger[Ul,tt][[1,2]];
vv=If[SameQ[tt,1],nn=Flatten[Table[vv[[i,1]],{i,Length[vv]}],1];
vv=Flatten[{nn[[1,1]],Rest[nn]},1],Flatten[Table[vv[[i,1]],{i,Length[vv]}],1]];
vv1=Union[Flatten[Map[fCyc,vv],1]];
vv2=Union[Flatten[vv1]];
vv2=ReplaceAll[vv1,Table[vv2[[i]]->i,{i,Length[vv2]}]]]


fMidEdgeTransfDT0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfMidInput[Ul,tt];
fKLfromGraph[pr]]


fMidEdgeTransfKL0[Ul_,tt_]:=Module[{dd},
dd=fMidEdgeTransfDT0[Ul,tt];
dd=fPDataFromDow[dd];
If[Plus@@dd[[1]]>89,0,ShowKnotfromPdataNew[dd]]]


fMidEdgeTransfPlot0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfMidInput[Ul,tt];
fPolyPlot[pr]]


fMidEdgeTransfGraphPlot0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfMidInput[Ul,tt];
DrawPlanarEmbGraph[pr]]


fTruncTransfDT0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfTruncInput[Ul,tt];
fKLfromGraph[pr]]


fTruncTransfKL0[Ul_,tt_]:=Module[{dd},
dd=fTruncTransfDT0[Ul,tt];
dd=fPDataFromDow[dd];
If[Plus@@dd[[1]]>89,0,ShowKnotfromPdataNew[dd]]]


fTruncTransfPlot0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfTruncInput[Ul,tt];
fPolyPlot[pr]]


fTruncTransfGraphPlot0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfTruncInput[Ul,tt];
DrawPlanarEmbGraph[pr]]


fJaegerDT1[Ul_,n_Integer]:=Module[{pr,jj,qq},
pr=fPolyFunctions1[Ul,n];
jj=fJaeger[pr];
qq=Union[fGraphInc[jj][[1]]];
Abs[jj]]


fJaegerKL1[Ul_,n_Integer]:=Module[{dd},
dd=fJaegerDT1[Ul,n];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fJaegerPlot1[Ul_,n_Integer]:=Module[{pr,jj,qq},
pr=fPolyFunctions1[Ul,n];
jj=fJaeger[pr];
qq=Union[fGraphInc[jj][[1]]];
fPolyPlot[qq]]


fJaegerGraphPlot1[Ul_,n_Integer]:=Module[{pr,jj,qq},
pr=fPolyFunctions1[Ul,n];
jj=fJaeger[pr];
qq=Union[fGraphInc[jj][[1]]];
DrawPlanarEmbGraph[qq]]


fJaegerDT2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj,qq},
pr=fPolyFunctions2[Ul,n,k];
jj=fJaeger[pr];
qq=Union[fGraphInc[jj][[1]]];
Abs[jj]]


fJaegerKL2[Ul_,n_Integer,k_Integer]:=Module[{dd},
dd=fJaegerDT2[Ul,n,k];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fJaegerPlot2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj,qq},
pr=fPolyFunctions2[Ul,n,k];
jj=fJaeger[pr];
qq=Union[fGraphInc[jj][[1]]];
fPolyPlot[qq]]


fJaegerGraphPlot2[Ul_,n_Integer,k_Integer]:=Module[{pr,jj,qq},
pr=fPolyFunctions2[Ul,n,k];
jj=fJaeger[pr];
qq=Union[fGraphInc[jj][[1]]];
DrawPlanarEmbGraph[qq]]


fJaegerTransfDT0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfJaegerInput[Ul,tt];
fJaeger[pr]]


fJaegerTransfKL0[Ul_,tt_]:=Module[{pr,dd},
pr=fPolyTransfJaegerInput[Ul,tt];
dd=fJaeger[pr];
If[Plus@@dd[[1]]>89,0,dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]]


fJaegerTransfPlot0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfJaegerInput[Ul,tt];
fPolyPlot[pr]]


fJaegerTransfGraphPlot0[Ul_,tt_]:=Module[{pr},
pr=fPolyTransfJaegerInput[Ul,tt];
DrawPlanarEmbGraph[pr]]


fTruncDT[Ul_]:=Module[{rr,res,dd},
rr=fPlanarEmbGraph[Ul][[3]];
res=fTruncation[rr];
dd=fKLfromGraph[res]]


fMidEdgeDTBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fMidEdgeGraph[pr];
fKLfromGraph[jj]]


fMidEdgeKLBasic[Ul_]:=Module[{dd},
dd=fMidEdgeDTBasic[Ul];
dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]


fMidEdgePlotBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fMidEdgeGraph[pr];
fPolyPlot[jj]]


fMidEdgeGraphPlotBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fMidEdgeGraph[pr];
DrawPlanarEmbGraph[jj]]


fTruncDTBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fKLfromGraph[fGraphInc[fTruncDT[pr]][[1]]]]


fTruncKLBasic[Ul_]:=Module[{dd},
dd=fTruncDTBasic[Ul];
dd=fPDataFromDow[dd];
ShowKnotfromPdataNew[dd]]


fTruncPlotBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fGraphInc[fTruncDT[pr]][[1]];
fPolyPlot[jj]]


fTruncGraphPlotBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fGraphInc[fTruncDT[pr]][[1]];
DrawPlanarEmbGraph[jj]]


fJaegerDTBasic[Ul_]:=Module[{pr,jj},
pr=fGraphInc[Ul][[1]];
jj=fJaeger[pr]]


fJaegerKLBasic[Ul_]:=Module[{dd},
dd=fJaegerDTBasic[Ul];
dd=fPDataFromDow[dd];
If[Plus@@dd[[1]]>89,0,ShowKnotfromPdataNew[dd]]]


fJaegerPlotBasic[Ul_]:=Module[{pr},
pr=fPolyPlot[Union[fGraphInc[fJaegerDTBasic[Ul]][[1]]]]]


fJaegerGraphPlotBasic[Ul_]:=Module[{pr},
pr=DrawPlanarEmbGraph[Union[fGraphInc[fJaegerDTBasic[Ul]][[1]]]]]


fMidEdgeIteration[Ul_String,n_Integer,k_Integer]:=Module[{ggg,pr,j},
j=1;
pr=fPolyFunctions1[Ul,n];
If[SameQ[k,1],pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False,ImageSize->700],
While[j<= k, pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}]; j++];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False,ImageSize->700]]]


fMidEdgeIteration1[Ul_String,n_Integer,e_Integer,k_Integer]:=Module[{ggg,pr,j},
j=1;
pr=fPolyFunctions2[Ul,n,e];
If[SameQ[k,1],pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False,ImageSize->700],
While[j<= k, pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}]; j++];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False,ImageSize->700]]]


fMidEdgeIteration2[Ul_String,tt_Integer,k_Integer]:=Module[{ggg,pr,j},
j=1;
pr=fPolyTransfMidInput[Ul,tt];ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}];
If[SameQ[k,1],pr=fPolyTransfMidInput[Ul,tt];GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False, ImageSize->700],
While[j<= k-1, pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}]; j++];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False, ImageSize->700]]]


fMidEdgeIterationBasic[Ul_String,k_Integer]:=Module[{ggg,pr,j},
j=1;
pr=fGraphInc[Ul][[1]];
If[SameQ[k,1],pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False,ImageSize->700],
While[j<= k, pr=fMidEdgeGraph[pr];
ggg=Table[pr[[i,1]]->pr[[i,2]],{i,Length[pr]}]; j++];
GraphPlot3D[ggg,EdgeRenderingFunction->({Cylinder[#1,0.05]} &),VertexRenderingFunction->({Sphere[#,0.15]} &),Boxed->False,ImageSize->700]]]


fPlantriPolytopes[n_Integer,k_Integer]:=Module[{ss,ss1,ss2,ss3,ss4,ss5},
ss=Import["C:\\Linknot\\fullerenes\\res"<>ToString[n],"Binary"];
ss1=Drop[ss,16];
ss2=Append[Prepend[Drop[Map[Last,Partition[Flatten[Position[ss1,0]],n]],-1],-1],Length[ss1]];
ss3=Table[Drop[Take[ss1,{ss2[[i]]+2,ss2[[i+1]]}],-1],{i,Length[ss2]-1}];
ss4=Table[ToExpression[StringJoin["{",StringReplace[ToString[ss3[[i]]],{", 0,"->"},{"," "->""}],"}"]],{i,Length[ss3]}];
ss5=Table[Union[Flatten[Table[Table[Sort[{j,ss4[[l,j,i]]}],{i,Length[ss4[[l,j]]]}],{j,Length[ss4[[l]]]}],1]],{l,Length[ss4]}];
If[k>Length[ss5],0,ss5[[k]]]]


fPlantriTri[n_Integer,k_Integer]:=Module[{ss,ss1,ss2,ss3,ss4,ss5},
If[n>18||OddQ[n],0,
ss=Import["C:\\Linknot\\fullerenes\\tri"<>ToString[n],"Binary"];
ss1=Drop[ss,16];
ss2=Append[Prepend[Drop[Map[Last,Partition[Flatten[Position[ss1,0]],n]],-1],-1],Length[ss1]];
ss3=Table[Drop[Take[ss1,{ss2[[i]]+2,ss2[[i+1]]}],-1],{i,Length[ss2]-1}];
ss4=Table[ToExpression[StringJoin["{",StringReplace[ToString[ss3[[i]]],{", 0,"->"},{"," "->""}],"}"]],{i,Length[ss3]}];
ss5=Table[Union[Flatten[Table[Table[Sort[{j,ss4[[l,j,i]]}],{i,Length[ss4[[l,j]]]}],{j,Length[ss4[[l]]]}],1]],{l,Length[ss4]}];
If[k>Length[ss5],0,ss5[[k]]]]]


fDoubleFull[Ulaz_]:=Module[{Ul,ss0,pp,ss,ss1,ss2,rr,res,res1,dd,gg},
Ul=Map[Sort,Ulaz];
ss0=Union[Flatten[Ul]];
pp=ss0;
ss=Ul[[1]];
ss1=Select[Ul,SameQ[Intersection[#,ss],{}] &][[1]];
While[Length[ss0]>2,
ss2=Join[ss,ss1];
ss=Flatten[Union[ss2]];
ss0=Select[Ul,SameQ[Intersection[#,ss],{}] &];
ss1=Select[Ul,SameQ[Intersection[#,ss],{}] &][[1]]];
rr=Partition[ss2,2];
res=Sort[Join[rr,Ul,{ss1}]];
res1=Table[Count[Flatten[res],i],{i,Length[pp]}];
dd=Flatten[Position[res1,3]];
res=If[SameQ[dd,{}],res,{}];
gg=Flatten[Union[res]];
If[Length[Union[Table[Count[gg,i],{i,Length[Union[gg]]}]]]>1,{},res]]


fDoubleFull1[Ulaz_]:=Module[{gg1,hhh,res},
gg1=Map[Sort,fCyc[HamiltonianPath[FromUnorderedPairs[Ulaz]]]];
hhh=Complement[Ulaz,gg1];
res=Join[hhh,Ulaz]]


fDoubleEdgesFull[Ulaz_]:=Module[{res},
res=fDoubleFull[Ulaz];
res=If[SameQ[res,{}],fDoubleFull[Reverse[Ulaz]],res];
If[SameQ[res,{}],fDoubleFull1[Ulaz],res]]


fFullList[n_Integer,k_Integer]:=Module[{Ul,nova,ss,ss1,i,pp,pp1,ss2,ss3,qqq,rr},
nova=If[n<20||SameQ[n,22]||n>120||OddQ[n],0,
ss=If[20<=n<=59,Import["C:\\LinKnot\\fullerenes\\Full_Codes_"<>ToString[n],"Binary"],Import["C:\\LinKnot\\fullerenes\\Full_Codes_"<>ToString[n]<>"_ipr","Binary"]];
ss1=Drop[ss,19];
i=Divide[IntegerPart[Divide[Length[ss1],n]],4];
pp=Drop[Table[4*n*(j-1)+j-1,{j,1,i}],1];
pp1=Table[{pp[[i]]},{i,Length[pp]}];
pp1=Select[pp1,#[[1]]<=Length[ss1] &];
ss2=Partition[Delete[ss1,pp1],4];
ss3=Table[Drop[ss2[[i]],-1],{i,Length[ss2]}];
qqq=Partition[Table[Drop[ss2[[i]],-1],{i,Length[ss2]}],n];
nova=Table[Ul=qqq[[j]];
rr=Union[Flatten[Table[Map[Sort,{{i,Ul[[i,1]]},{i,Ul[[i,2]]},{i,Ul[[i,3]]}}],{i,Length[Ul]}],1]],{j,Length[qqq]}][[k]]];
nova]


fOneParPoly[No_Integer,n_Integer,tt_Integer]:=Module[{Ul,res},
Ul=polyhedra01[[No,2]];
res=If[SameQ[tt,1],{fPolyPlot[fPolyFunctions1[Ul,n]],fMidEdgeDT1[Ul,n],fMidEdgeKL1[Ul,n],fMidEdgePlot1[Ul,n],fMidEdgeGraphPlot1[Ul,n]},If[SameQ[tt,2],{fPolyPlot[fPolyFunctions1[Ul,n]],fTruncDT1[Ul,n],fTruncKL1[Ul,n],fTruncPlot1[Ul,n],fTruncGraphPlot1[Ul,n]},{fPolyPlot[fPolyFunctions1[Ul,n]],fJaegerDT1[Ul,n],fJaegerKL1[Ul,n],fJaegerPlot1[Ul,n],fJaegerGraphPlot1[Ul,n]}]]]


fTwoParPoly[No_Integer,n_Integer,k_Integer,type_Integer]:=Module[{Ul,res},
Ul=polyhedra02[[No,2]];
res=If[SameQ[type,1],{fPolyPlot[fPolyFunctions2[Ul,n,k]],fMidEdgeDT2[Ul,n,k],fMidEdgeKL2[Ul,n,k],fMidEdgePlot2[Ul,n,k],fMidEdgeGraphPlot2[Ul,n,k]},If[SameQ[type,2],{fPolyPlot[fPolyFunctions2[Ul,n,k]],fTruncDT2[Ul,n,k],fTruncKL2[Ul,n,k],fTruncPlot2[Ul,n,k],fTruncGraphPlot2[Ul,n,k]},{fPolyPlot[fPolyFunctions2[Ul,n,k]],fJaegerDT2[Ul,n,k],fJaegerKL2[Ul,n,k],fJaegerPlot2[Ul,n,k],fJaegerGraphPlot2[Ul,n,k]} ] ]]


fAllPolyMid[No_Integer]:=Module[{Ul,res},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,0],fMidEdgeTransfPlot0[Ul,0],fPolyTransf[Ul,1],fMidEdgeTransfPlot0[Ul,1],fPolyTransf[Ul,2],fMidEdgeTransfPlot0[Ul,2],fPolyTransf[Ul,3],fMidEdgeTransfPlot0[Ul,3]}]


fAllPolyCross[No_Integer]:=Module[{Ul,res},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,0],fTruncTransfPlot0[Ul,0],fPolyTransf[Ul,1],fTruncTransfPlot0[Ul,1],fPolyTransf[Ul,2],fTruncTransfPlot0[Ul,2],fPolyTransf[Ul,3],fTruncTransfPlot0[Ul,3]}]


fAllPolyJaeger[No_Integer]:=Module[{Ul,res},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,0],fJaegerTransfPlot0[Ul,0],fPolyTransf[Ul,1],fJaegerTransfPlot0[Ul,1],fPolyTransf[Ul,2],fJaegerTransfPlot0[Ul,2],fPolyTransf[Ul,3],fJaegerTransfPlot0[Ul,3]}]


fBasPoly[Ul_String,tt_Integer]:=Module[{res},
res=If[SameQ[tt,1],{MinDowProjAltKL[fKnotscapeDow[Ul]],ShowKnotfromPdataNew[fCreatePData[Ul]],fMidEdgeDTBasic[Ul],fMidEdgeKLBasic[Ul],fMidEdgePlotBasic[Ul],fMidEdgeGraphPlotBasic[Ul]},If[SameQ[tt,2],{MinDowProjAltKL[fKnotscapeDow[Ul]],ShowKnotfromPdataNew[fCreatePData[Ul]],fTruncDTBasic[Ul],fTruncKLBasic[Ul],fTruncPlotBasic[Ul],fTruncGraphPlotBasic[Ul]},{MinDowProjAltKL[fKnotscapeDow[Ul]],ShowKnotfromPdataNew[fCreatePData[Ul]],fJaegerDTBasic[Ul],fJaegerKLBasic[Ul],fJaegerPlotBasic[Ul],fJaegerGraphPlotBasic[Ul]}]]]


fGraph3DKL[Ul_]:=Module[{res},
res={If[SameQ[Head[Ul],String],ShowKnotfromPdataNew[fCreatePData[Ul]],ShowKnotfromPdataNew[fPDataFromDow[Ul]]],If[SameQ[Head[Ul],String],fPolyPlotKL[Ul],fPolyPlotKL[fSignsKL[Ul]]]}]


fFourVal1[No_Integer,n_Integer]:=Module[{Ul,res},
Ul=polyhedra04[[No,2]];
res={fPolyPlot[fPolyFunctions4[Ul,n]],fKLfromGraph[fPolyFunctions4[Ul,n]],ShowKnotfromPdataNew[fPDataFromDow[fKLfromGraph[fPolyFunctions4[Ul,n]]]],fPolyPlotKL[fSignsKL[fKLfromGraph[fPolyFunctions4[Ul,n]]]]}]


fFourVal2[No_Integer,n_Integer,k_Integer]:=Module[{Ul,res},
Ul=polyhedra05[[No,2]];
res={fPolyPlot[fPolyFunctions41[Ul,n,k]],fKLfromGraph[fPolyFunctions41[Ul,n,k]],ShowKnotfromPdataNew[fPDataFromDow[fKLfromGraph[fPolyFunctions41[Ul,n,k]]]],fPolyPlotKL[fSignsKL[fKLfromGraph[fPolyFunctions41[Ul,n,k]]]]}]


fIterativeMid1[No_Integer,n_Integer,k_Integer]:=Module[{Ul,res},
Ul=polyhedra01[[No,2]];
res={fPolyPlot[fPolyFunctions1[Ul,n]],fMidEdgeIteration[Ul,n,k]}]


fIterativeMid2[No_Integer,n_Integer,k_Integer,s_Integer]:=Module[{Ul,res},
Ul=polyhedra02[[No,2]];
res={fPolyPlot[fPolyFunctions2[Ul,n,k]],fMidEdgeIteration1[Ul,n,k,s]}]


fIterativeMid3[No_Integer,t_Integer,s_Integer]:=Module[{Ul,res},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,t],fMidEdgeIteration2[Ul,t,s]}]


fIterativeBasic[Ul_String,s_Integer]:=Module[{res},
res={MinDowProjAltKL[fKnotscapeDow[Ul]],ShowKnotfromPdataNew[fCreatePData[Ul]],fMidEdgeIterationBasic[Ul,s]}]


fSurfaceKL[Ul_]:=Module[{res},
res={ShowKnotfromPdataNew[If[SameQ[Head[Ul],String],fCreatePData[Ul],fPDataFromDow[Ul]]],ShowKnot3DNew[If[SameQ[Head[Ul],String],fCreatePData[Ul],fPDataFromDow[Ul]]],ShowKnot3DNewTop[If[SameQ[Head[Ul],String],fCreatePData[Ul],fPDataFromDow[Ul]]],ShowKnot3DNew[If[SameQ[Head[Ul],String],fCreatePData[Ul],fPDataFromDow[Ul]]]}]


fOneParPolyMidSurfaces[No_Integer,n_Integer]:=Module[{Ul,res},
Ul=polyhedra01[[No,2]];
res={fPolyPlot[fPolyFunctions1[Ul,n]],fMidEdgeDT1[Ul,n],fMidEdgeKL1[Ul,n],ShowKnot3DNew[fPDataFromDow[fMidEdgeDT1[Ul,n]]],ShowKnot3DNewTop[fPDataFromDow[fMidEdgeDT1[Ul,n]]],ShowKnot3DNew[fPDataFromDow[fMidEdgeDT1[Ul,n]]]}]


fOneParPolyCrossSurfaces[No_Integer,n_Integer]:=Module[{Ul,res},
Ul=polyhedra01[[No,2]];
res={fPolyPlot[fPolyFunctions1[Ul,n]],fTruncDT1[Ul,n],fTruncKL1[Ul,n],ShowKnot3DNew[fPDataFromDow[fTruncDT1[Ul,n]]],ShowKnot3DNewTop[fPDataFromDow[fTruncDT1[Ul,n]]],ShowKnot3DNew[fPDataFromDow[fTruncDT1[Ul,n]]]}]


fTwoParPolyMidSurfaces[No_Integer,n_Integer,k_Integer]:=Module[{Ul,res},
Ul=polyhedra02[[No,2]];
res={fPolyPlot[fPolyFunctions2[Ul,n,k]],fMidEdgeDT2[Ul,n,k],fMidEdgeKL2[Ul,n,k],ShowKnot3DNew[fPDataFromDow[fMidEdgeDT2[Ul,n,k]]],ShowKnot3DNewTop[fPDataFromDow[fMidEdgeDT2[Ul,n,k]]],ShowKnot3DNew[fPDataFromDow[fMidEdgeDT2[Ul,n,k]]]}]


fTwoParPolyCrossSurfaces[No_Integer,n_Integer,k_Integer]:=Module[{Ul,res},
Ul=polyhedra02[[No,2]];
res={fPolyPlot[fPolyFunctions2[Ul,n,k]],fTruncDT2[Ul,n,k],fTruncKL2[Ul,n,k],ShowKnot3DNew[fPDataFromDow[fTruncDT2[Ul,n,k]]],ShowKnot3DNewTop[fPDataFromDow[fTruncDT2[Ul,n,k]]],ShowKnot3DNew[fPDataFromDow[fTruncDT2[Ul,n,k]]]}]


fAllPolyMidSurfaces[No_Integer]:=Module[{Ul,res,qq},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,0],fMidEdgeTransfPlot0[Ul,0],qq=fMidEdgeTransfDT0[Ul,0];If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,qq],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,fMidEdgeTransfKL0[Ul,0]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNew[fPDataFromDow[qq]]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNewTop[fPDataFromDow[qq]]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNew[fPDataFromDow[qq]]]}]


fAllPolyMidTruncSurfaces[No_Integer]:=Module[{Ul,res,qq},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,1],fMidEdgeTransfPlot0[Ul,1],fMidEdgeTransfDT0[Ul,1],qq=fMidEdgeTransfDT0[Ul,1];If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,qq],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,fMidEdgeTransfKL0[Ul,1]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNew[fPDataFromDow[qq]]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNewTop[fPDataFromDow[qq]]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNew[fPDataFromDow[qq]]]}]


fAllPolyCrossSurfaces[No_Integer]:=Module[{Ul,res,qq},
Ul=polyhedra03[[No,2]];
res={fPolyTransf[Ul,0],qq=fTruncTransfDT0[Ul,0];If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,qq],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,fTruncTransfKL0[Ul,0]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNew[fPDataFromDow[qq]]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNewTop[fPDataFromDow[qq]]],If[Plus@@qq[[1]]>49||Plus@@qq[[1]]<2,0,ShowKnot3DNew[fPDataFromDow[qq]]]}]


fPlantriMid[n_Integer,k_Integer]:=Module[{res},
res={fPlantriPolytopes[n,k],fPolyPlot[fPlantriPolytopes[n,k]],DrawPlanarEmbGraph[fPlantriPolytopes[n,k]],fPolyPlot[fMidEdgeGraph[fPlantriPolytopes[n,k]]],fKLfromGraph[fMidEdgeGraph[fPlantriPolytopes[n,k]]],ShowKnotfromPdataNew[fSignsKL[fKLfromGraph[fMidEdgeGraph[fPlantriPolytopes[n,k]]]]]}]


fPlantriCross[n_Integer,k_Integer]:=Module[{res},
res={fPlantriPolytopes[n,k],fPolyPlot[fPlantriPolytopes[n,k]],DrawPlanarEmbGraph[fPlantriPolytopes[n,k]],fPolyPlot[fTruncation[fPlanarEmbGraph[fPlantriPolytopes[n,k]][[3]]]],fKLfromGraph[fGraphInc[fTruncDT[fPlantriPolytopes[n,k]]][[1]]],ShowKnotfromPdataNew[fSignsKL[fKLfromGraph[fGraphInc[fTruncDT[fPlantriPolytopes[n,k]]][[1]]]]]}]


fGeneralFull[n_Integer,k_Integer]:=Module[{res},
res={fPlantriTri[n,k],fPolyPlot[fPlantriTri[n,k]],fAddDig[fPlantriTri[n,k]][[-1]]}]


fGeneralFullDouble[n_Integer,k_Integer,k1_Integer]:=Module[{res},
res={Drop[fAddDig[fPlantriTri[n,k]]][[k1]],fPolyPlotKL[fSignsKL[fKLfromGraph[Drop[fAddDig[fPlantriTri[n,k]]][[k1]]]]],fKLfromGraph[Drop[fAddDig[fPlantriTri[n,k]]][[k1]]],ShowKnotfromPdataNew[fSignsKL[fKLfromGraph[Drop[fAddDig[fPlantriTri[n,k]]][[k1]]]]]}]


fFullgen[n_Integer,k_Integer]:=Module[{res,rr},
res={fFullList[n,k],rr=fDoubleEdgesFull[fFullList[n,k]];fPolyPlot[rr],fPolyPlotDoubleEdges[rr]}]


fVirtPlanarLink[Ul_String]:=Module[{ss,ss1,gg,gg1,ggg,pp,pp1,vv,vv1,vv2,vv3},
ss=StringReplace[Ul,"i"->"1"];
ss1=StringReplace[Ul,"i"->"-1"];
gg=fGaussExtSigns[ss];
gg1=fGaussExtSigns[ss1];
ggg=Map[Sign,Abs[gg-gg1]];
pp=Position[ggg,1];
pp1=Table[gg[[pp[[i,1]],pp[[i,2]]]],{i,Length[pp]}];
vv=Abs[Table[Select[gg[[i]],Not[MemberQ[pp1,#]] &],{i,Length[gg]}]];
vv=Complement[vv,{{}}];
vv1=Union[Flatten[vv]];
vv2=ReplaceAll[vv,Table[vv1[[i]]->i,{i,Length[vv1]}]];
vv3=Union[Flatten[Map[fCyc,vv2],1]];
PlanarQ[FromUnorderedPairs[vv3]]]


fVirtPlanar[Ul_String]:=Module[{cc,ss,ll},
cc=fComponentNo[StringReplace[Ul,"i"->"1"]];
If[SameQ[cc,1],
ss=fGaussVirtKnot[Ul];
ll=ToExpression[StringJoin["{",StringDrop[StringReplace[ss,{"+"->"","-"->"","O"->",","U"->","}],1],"}"]];
PlanarQ[FromUnorderedPairs[fCyc[ll]]],fVirtPlanarLink[Ul]]]


fVirtLinkGraph[Ul_String]:=Module[{ss,ss1,gg,gg1,ggg,pp,pp1,vv,vv1,vv2,vv3},
ss=StringReplace[Ul,"i"->"1"];
ss1=StringReplace[Ul,"i"->"-1"];
gg=fGaussExtSigns[ss];
gg1=fGaussExtSigns[ss1];
ggg=Map[Sign,Abs[gg-gg1]];
pp=Position[ggg,1];
pp1=Table[gg[[pp[[i,1]],pp[[i,2]]]],{i,Length[pp]}];
vv=Abs[Table[Select[gg[[i]],Not[MemberQ[pp1,#]] &],{i,Length[gg]}]];
vv=Complement[vv,{{}}];
vv1=Union[Flatten[vv]];
vv2=ReplaceAll[vv,Table[vv1[[i]]->i,{i,Length[vv1]}]];
vv3=Flatten[Map[fCyc,vv2],1];
vv3]


fVirtGraph[Ul_String]:=Module[{cc,ss,ll},
cc=fComponentNo[StringReplace[Ul,"i"->"1"]];
If[SameQ[cc,1],
ss=fGaussVirtKnot[Ul];
ll=ToExpression[StringJoin["{",StringDrop[StringReplace[ss,{"+"->"","-"->"","O"->",","U"->","}],1],"}"]];
fCyc[ll],fVirtLinkGraph[Ul]]]


           
End[]      
EndPackage[]


