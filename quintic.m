(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


(* ::Subsection::Closed:: *)
(*Utility*)


Clear[Evaluate[$Context<>"*"]];

time0=AbsoluteTime[];
gMax[x_]:=(10+5x+x^2)/10;
dMin[g_]:=Ceiling[5/2 (-1+Sqrt[-3+8 g]/Sqrt[5])-1];

(*** INPUT ***)
options = $CommandLine;
param[flag_] := Module[
		{position, flagList}
	, 
		flagList = StringSplit[options];
		position = Position[flagList, "-" <> flag];
		Switch[Length[position], 
			1, If[position[[1, 1]] >= Length[flagList], 
				True, 
				flagList[[position[[1, 1]]+1]][[1]]
				], 
			0, (*WriteString["stdout", "flag -" <> flag <> " unspecified\n"];*)
				Null, 
			_, WriteString["stdout", "flag -" <> flag <> " duplicated\n"];
				Abort[]
		]
	];


(* ::Subsection:: *)
(*Initialization*)


(*options={"-g","11"(*,"-d","4"*)};*)

genus=param["g"]//ToExpression;
degree=param["d"]//ToExpression;
If[degree===Null,degree=dMin[genus]+1];
Print["g="<>ToString[genus]<>", d="<>ToString[degree]];

id=ToString[genus]<>"_"<>ToString[degree];

(*directory=NotebookDirectory[];
dataDirectory=NotebookDirectory[]<>"data/";*)
directory=$HomeDirectory<>"/gv/";
dataDirectory=Environment["SCRATCH"]<>"/yin_lab/Users/yhlin/gv/";

resDirectory=directory<>"res/";
timeDirectory=directory<>"time/";
If[!FileExistsQ[directory],CreateDirectory[directory]];
If[!FileExistsQ[resDirectory],CreateDirectory[resDirectory]];
If[!FileExistsQ[timeDirectory],CreateDirectory[timeDirectory]];
If[!FileExistsQ[dataDirectory],CreateDirectory[dataDirectory]];


(* ::Subsection::Closed:: *)
(*Linear Solve*)


LS[eqns_]:=Module[{vars,m,sol},
	vars=Variables[eqns];
	m=CoefficientArrays[eqns,vars];
	sol=LinearSolve[m[[2]],-m[[1]]];
	Table[vars[[i]]->sol[[i]],{i,1,Length[vars]}]
];


(* ::Section::Closed:: *)
(*Moduli space geometry*)


(* ::Subsection::Closed:: *)
(*Helper*)


logD[series_,p_,z_]:=If[p==0,series,z D[logD[series,p-1,z],z]];


(* ::Subsection::Closed:: *)
(*Orbifold point*)


Ncut0=(degree+1)+Ceiling[3/5(genus-1)];
\[Omega]0[k_,z_]:=25z^((k+1)/5) Series[Sum[Pochhammer[(k+1)/5,n]^5/Pochhammer[k+1,5n] (5^5 z)^n,{n,0,Ncut0}],{z,0,Ncut0+1}];
t0[z_]:=z^(1/5)/(2\[Pi] I) \[Omega]0[1,z]/(z^(1/5) \[Omega]0[0,z]);
G0[z_]:=D[t0[z],z];
emK0[z_]:=\[Omega]0[0,z];
A0[p_,z_]:=A0[p,z]=logD[G0[z],p,z]/G0[z];
B0[p_,z_]:=B0[p,z]=logD[emK0[z],p,z]/emK0[z];
X0[z_]:=1/(1-z);
z0:=z0=InverseSeries[2\[Pi] I t0[z],s0]//Normal;


(* ::Subsection::Closed:: *)
(*Large volume point*)


NcutInf=(degree+1)+genus;
\[Omega]InfPre[z_,x_]:=Sum[Gamma[5n+5x+1]/Gamma[n+x+1]^5 1/(5^5 z)^(n+x),{n,0,NcutInf}];
\[Omega]Inf[0,z_]:=\[Omega]Inf[0,z]=Series[\[Omega]InfPre[z,0],{z,Infinity,NcutInf}];
\[Omega]Inf[1,z_]:=\[Omega]Inf[1,z]=D[\[Omega]InfPre[z,x],x]/.x->0;
tInf[z_]:=tInf[z]=1/(2\[Pi] I) \[Omega]Inf[1,z]/\[Omega]Inf[0,z];
GInf[z_]:=GInf[z]=D[tInf[z],z];
emKInf[z_]:=emKInf[z]=\[Omega]Inf[0,z];
AInf[p_,z_]:=AInf[p,z]=logD[GInf[z],p,z]/GInf[z];
BInf[p_,z_]:=BInf[p,z]=logD[emKInf[z],p,z]/emKInf[z];
XInf[z_]:=1/(1-z);
zInf:=zInf=InverseSeries[Exp[2\[Pi] I tInf[z]],q]//Normal;


(* ::Subsection::Closed:: *)
(*Conifold point*)


Ncutc=(degree+1)+2genus-2;
z=Exp[u];
L[c_]:=(D[#1,{u,1}]-c/5 #1)&;
PF=(D[#1,{u,4}]-Exp[-u]L[1][L[2][L[3][L[4][#1]]]])&;

x[u_]:=Exp[u]-1;
\[Omega][u_]:=Sum[c[i]x[u]^i,{i,0,Ncutc}];

PF\[Omega]u=PF[\[Omega][u]]//Expand;

Clear[c];
c[0]=c0[0]=1;c[1]=c0[1]=0;c[2]=c0[2]=0;c[3]=c0[3]=2/625;
sol=LS[CoefficientList[Series[PF\[Omega]u,{u,0,Ncutc-3}],u][[2;;]]];
Do[c0[i]=c[i]/.sol,{i,4,Ncutc}];
\[Omega]c[0,Z_]:=\[Omega]c[0,Z]=Sqrt[5]Sum[c0[i](Z-1)^i,{i,0,Ncutc}];

Clear[c];
c[0]=c1[0]=0;c[1]=c1[1]=1;c[2]=c1[2]=-3/10;
sol=LS[CoefficientList[Series[PF\[Omega]u,{u,0,Ncutc-3}],u]];
Do[c1[i]=c[i]/.sol,{i,3,Ncutc}];
\[Omega]c[1,Z_]:=\[Omega]c[1,Z]=Sqrt[5]Sum[c1[i](Z-1)^i,{i,0,Ncutc}];

Clear[c,d,S,x,z];
tc[z_]:=tc[z]=Series[1/(2\[Pi] I) \[Omega]c[1,z]/\[Omega]c[0,z],{z,1,Ncutc}];

Clear[A,B,X,G,emK];
Gc[z_]:=Gc[z]=D[tc[z],z];
emKc[z_]:=emKc[z]=\[Omega]c[0,z];

Ac[p_,z_]:=Ac[p,z]=logD[Gc[z],p,z]/Gc[z];
Bc[p_,z_]:=Bc[p,z]=logD[emKc[z],p,z]/emKc[z];
Xc[z_]:=1/(1-z);

zc:=zc=InverseSeries[Series[2\[Pi] I tc[z],{z,1,Ncutc}]//Simplify,s1]//Normal;


(* ::Subsection::Closed:: *)
(*Time*)


dt=AbsoluteTime[]-time0;
Print["moduli space geometry took ",dt,"s","\n"];
time0=AbsoluteTime[];

file=timeDirectory<>id<>".csv";
If[FileExistsQ[file],DeleteFile[file]];
WriteString[file,dt,"\n"];


(* ::Section::Closed:: *)
(*Genus 0 and 1*)


(* ::Subsection::Closed:: *)
(*GV*)


F[g_,q_]:=Sum[Sum[2^(2(r-g))/Gamma[2g-2+1] n[r,d]SeriesCoefficient[D[Sin[x]^(2r-2),{x,2g-2}],{x,0,0}]PolyLog[3-2g,q^d],{d,1,degree}],{r,0,g}];


(* ::Subsection::Closed:: *)
(*GV genus 0*)


(*

(* Calculation of GV invariants for genus 0 *)
Clear[Ph,a,b,x,A,B,f,n,u,g,z];

(* Calculation of three point function in terms of pre-potential *)
Ph[0,0,u_]:=f[Exp[u]];
(a^(0,1))[n_,u_]:=a[n+1,u]-a[1,u]a[n,u];
(b^(0,1))[n_,u_]:=b[n+1,u]-b[1,u]b[n,u];
x'[u_]:=x[u](x[u]-1);
Ph[g_,n_,u_]:=D[Ph[g,n-1,u],u]-((n-1)(a[1,u]+1)+(2-2g)(b[1,u]-1/2 x[u]))Ph[g,n-1,u];

ThreePointFunctiond=Ph[0,3,u]/.u->Log[z]//Expand;
a[n_,u_]:=A[n, Exp[u]];
b[n_,u_]:=B[n, Exp[u]];
X[z_] := 1/(1-z);
x[u_]:=X[Exp[u]];
ThreePointFunctionP=ThreePointFunctiond;

Clear[etinf,EtInf];
etinf=Exp[2\[Pi] I tInf[z]];
EtInf[d_,z_]:=EtInf[d,z]=If[d>0,EtInf[d-1,z]*etinf,1];
fInf[z_]:=Module[{tinf,\[Omega]inf},
	tinf=tInf[z];
	\[Omega]inf=\[Omega]Inf[0,z];
	(5/6 (2\[Pi] I tinf)^3-0*a/2 (2\[Pi] I tinf)^2+0*c (2\[Pi] I tinf)+0*e/2+Sum[n[0,d] PolyLog[3,EtInf[d,z]],{d,1,degree}])\[Omega]inf^2 ((1-z)/(5z))
];

(*(* pre-potential in terms of GV invariants *)
fInf[z_]:=Module[{tinf,\[Omega]inf},
	tinf=tInf[z];
	\[Omega]inf=\[Omega]Inf[0,z];
	(5/6 (2\[Pi] I tinf)^3-0*a/2 (2\[Pi] I tinf)^2+0*c (2\[Pi] I tinf)+0*e/2+Sum[n[0,d]PolyLog[3,Exp[2\[Pi] I d tinf]],{d,1,degree}])\[Omega]inf^2 ((1-z)/(5z))
];*)

(*Clear[ff];*)

(*Timing[ff[0]=fInf[z]//Simplify][[1]]//Print;*)

ff[0]:=fInf[z]//Simplify;

ff[n_]:=(ff[n]=D[ff[n-1],z])/;n>0;

Timing[ff[3]][[1]]//Print;

zInfq=Series[zInf,{q,0,(degree+1)}];
comps={ThreePointFunctionP/.{f[z]->0,f'[z]->0,f''[z]->0,f'''[z]->0},Coefficient[ThreePointFunctionP,f[z]],Coefficient[ThreePointFunctionP,f'[z]],Coefficient[ThreePointFunctionP,f''[z]],Coefficient[ThreePointFunctionP,f'''[z]]};
A:=AInf;
B:=BInf;

comps=Normal[comps];
fs=Normal[Join[{1},Table[ff[i]/.Log[z]->0,{i,0,3}]]];
ThreePointFunction=comps . fs//Expand;
sol=LS[Table[SeriesCoefficient[ThreePointFunction,{z,Infinity,i}],{i,1,degree}]]//FullSimplify;
GV=sol;
Clear[Ph,a,b,x,A,B,f,n,u,g,z];

*)


(* ::Subsection::Closed:: *)
(*GV genus 0*)


X[z_] := 1/(1-z);
zInfq=Series[zInf,{q,0,(degree+1)}];
GV={};


(* ::Subsection::Closed:: *)
(*OEIS (A060041)*)


nn=degree+1;
y0[x_]:=Sum[(5n)!/(n!)^5 x^n, {n, 0, nn}];
y1[x_]:=Sum[((5n)!/(n!)^5 5 Sum[1/j, {j, n+1, 5n}]) x^n, {n, 0, nn}];
qq=Series[x Exp[y1[x]/y0[x]], {x, 0, nn}];
x[q_]=InverseSeries[qq, q];
s1=(q/x[q] D[x[q], q])^3 5/((1-5^5 x[q]) y0[x[q]]^2);
s2=Series[5+Sum[n[d] d^3 q^d/(1-q^d), {d, 1, nn}], {q, 0, nn}];
sol=Solve[#==0&/@CoefficientList[s1-s2,q]][[1]];
(* Daniel Grunberg (grunberg(AT)mpim-bonn.mpg.de), Aug 18 2004 *)
Do[n[0,d]=n[d]/.sol,{d,1,nn-1}];
Clear[nn,y0,y1,qq,x,s1,s2,sol];


(* ::Subsection::Closed:: *)
(*GV genus 1*)


If[genus>=1,
A[p_,z_]:=AInf[p,z];
B[p_,z_]:=BInf[p,z];
Ph[1,1,z_]:=-(1/2)A[1,z]-31/3 B[1,z]+1/12 (X[z]-1)+5/3;
Fh[1,q_]:=Normal[Integrate[Ph[1,1,z]/z,z]]/.z->zInfq; 
For[i=1,i<=degree,i++,
sol=Solve[SeriesCoefficient[(Fh[1,q]-F[1,q])/.GV,{q,0,i}]==0, n[1,i]][[1]];
GVd=AppendTo[GV,sol[[1]]];
GV=GVd;];
];


(* ::Subsection::Closed:: *)
(*Time*)


dt=AbsoluteTime[]-time0;
Print["genus 0 and 1 took ",dt,"s","\n"];
time0=AbsoluteTime[];

file=timeDirectory<>id<>".csv";
WriteString[file,dt,"\n"];


(* ::Section::Closed:: *)
(*BCOV and BC*)


(* ::Subsection::Closed:: *)
(*BCOV and Yamaguchi-Yau (x v u)*)


Clear[c,A,B,X,P1,P2,P3,P,Ph,Q,Qu,Qh,r,g,x,u,v1,v2,v3,x,v1,v2,v3,u];

Simp[x_]:=Expand[x];

(* The algorithm below treats operator D as zD in a way such that there is no mistake for the purpose of this particular calculation *)
x'[z_]:=x'[z]=(-1+x[z]) x[z];
v1'[z_]:=v1'[z]=-v1[z]^2+v1[z] x[z]-2/5 (5 v2[z]+x[z]);
v2'[z_]:=v2'[z]=-v1[z] v2[z]+v3[z];
v3'[z_]:=v3'[z]=v2[z]^2-(1+v1[z]) v2[z] x[z]+(-(24/625)+2 v3[z]) x[z];
u'[z_]:=u'[z]=-u[z]^2+u[z] v1[z]+v2[z];

Ph[1,1,z_]:=Ph[1,1,z]=25/12-(28 u[z])/3-v1[z]/2+x[z]/12;

variables={x,v1,v2,v3,u};
Dz[expr_]:=Sum[D[expr,var[z]]*var'[z]//Simp,{var,variables}];

Ph[g_,n_,z_]:=(Ph[g,n,z]=Module[{t1,t2,ans1,ans2},
	{t1,ans1}=Dz[Ph[g,n-1,z]]//Timing;
	tt["Dz"<>ToString[n],g]=t1;
	{t2,ans2}=-((-1+n) (-2 u[z]+v1[z])+(2-2 g) (u[z]-x[z]/2))Ph[g,n-1,z]//Simp//Timing;
	tt["Ph"<>ToString[n],g]=t1+t2;
	ans1+ans2
	])/;n>0;
(*Ph[g_,n_,z_]:=(Ph[g,n,z]=Dz[Ph[g,n-1,z]]-((-1+n) (-2 u[z]+v1[z])+(2-2 g) (u[z]-x[z]/2))Ph[g,n-1,z]//Simp)/;n>0;*)

Q[g_,z_]:=Module[{t,ans,vec},
	Ph[g-1,2,z];
	vec=Table[Ph[g-r,1,z],{r,1,g-1}];
	{t,ans}=(1/2)(Ph[g-1,2,z]+Sum[Ph[g-r,1,z]Ph[r,1,z],{r,1,g-1}]//Simp)//Timing;
	tt["Q",g]=t;
	ans
	]/;g>1;
(*Q[g_,z_]:=(1/2)(Ph[g-1,2,z]+Sum[Ph[g-r,1,z]Ph[r,1,z],{r,1,g-1}]//Simp)/;g>1;*)

rule={X[z]->x,A[1,z]->v1-2u-1,B[1,z]->u,B[2,z]->v2+u v1,B[3,z]->v3+u(-v2+x(v1-2/5))};
ruleInv={x->X[z],v1->1+A[1,z]+2 B[1,z],v2->-B[1,z]-A[1,z] B[1,z]-2 B[1,z]^2+B[2,z],v3->1/5 (-5 B[1,z]^2-5 A[1,z] B[1,z]^2-10 B[1,z]^3+5 B[1,z] B[2,z]+5 B[3,z]-3 B[1,z] X[z]-5 A[1,z] B[1,z] X[z]-10 B[1,z]^2 X[z]),u->B[1,z]};

Qu[g_,z_]:=Qu[g,z]=CoefficientList[Q[g,z],u[z]];
Qh[g_,n_,z_]:=Qu[g,z][[n+1]];

(*P1[g_,0,z_]:=Integrate[Collect[-Qh[g,0,z]/.{v2[z]\[Rule]0,v3[z]\[Rule]0},v1[z]],{v1[z],0,v1[z]}];
P2[g_,0,z_]:=P1[g,0,z]+Integrate[Collect[Qh[g,1,z]-x[z] Qh[g,2,z]/.{v3[z]\[Rule]0},v2[z]],{v2[z],0,v2[z]}];
P3[g_,0,z_]:=P2[g,0,z]+Integrate[Collect[Qh[g,2,z],v3[z]],{v3[z],0,v3[z]}];*)
Int[expr_,var_]:=Module[{coefs},
	coefs=CoefficientList[expr,var];
	coefs . Table[var^i/i,{i,1,Length[coefs]}]
];
P1[g_,0,z_]:=Int[-Qh[g,0,z]/.{v2[z]->0,v3[z]->0},v1[z]];
P2[g_,0,z_]:=P1[g,0,z]+Int[Qh[g,1,z]-x[z] Qh[g,2,z]/.{v3[z]->0},v2[z]];
P3[g_,0,z_]:=P2[g,0,z]+Int[Qh[g,2,z],v3[z]];

Pz[g_]:=(P3[g,0,z]//Simp)+Sum[c[g,i]x[z]^i,{i,0,3g-3}]/;g>1;

(*P[g_]:=Pz[g]/.{x[z]->x,v1[z]->v1,v2[z]->v2,v3[z]->v3};*)
P[g_]:=Module[{t,ans},
	Qh[g,0,z];
	{t,ans}=Timing[Pz[g]/.{x[z]->x,v1[z]->v1,v2[z]->v2,v3[z]->v3}];
	tt["P",g]=t;
	ans
];

cs[g_,c_]:=Table[c[g,i],{i,0,3g-3}];

(* Important that X, A, B is not cleared *)
Ph[g_,0,z_]:=(Ph[g,0,z]=(Pz[g]//Simp))/;g>1;


(* ::Subsection::Closed:: *)
(*Orbifold point*)


Clear[ru0,xv0,F0];
ru0[g_]:=ru0[g]=#[[1]]->Series[#[[2]],{z,0,Floor[3/5 (-1+g)]}]&/@(ruleInv/.{X->X0,A->A0,B->B0});
xv0[g_][n1_,n2_,n3_,nx_]:=xv0[g][n1,n2,n3,nx]=v1^n1 v2^n2 v3^n3 x^nx/.ru0[g];
F0[g_,z_]:=F0[g,z]=(\[Omega]0[0,z](1+O[z]^(Floor[3/5 (-1+g)+1])))^(2g-2) ((1-z)/(5z))^(g-1) (Collect[P[g],variables]/.ru0[g]);
Match0[g_]:=Table[SeriesCoefficient[F0[g,z],{z,0,i}],{i,-(3/5) (-1+g),-1/5}]//Expand;


(* ::Subsection::Closed:: *)
(*Large volume point*)


Clear[ruInf,xvInf,FInf,FhInf];
ruInf[g_]:=ruInf[g]=#[[1]]->Series[#[[2]],{z,Infinity,(degree+1)}]&/@(ruleInv/.{X->XInf,A->AInf,B->BInf});
xvInf[g_][n1_,n2_,n3_,nx_]:=xvInf[g][n1,n2,n3,nx]=v1^n1 v2^n2 v3^n3 x^nx/.ruInf[g];
FInf[g_,z_]:=FInf[g,z]=Series[\[Omega]Inf[0,z],{z,Infinity,(degree+1)}]^(2g-2) ((1-z)/(5z))^(g-1) (Collect[P[g],variables]/.ruInf[g]);
FhInf[g_,z_]:=FhInf[g,z]=Series[FInf[g,z],{z,Infinity,0}];
MatchInf[g_]:={SeriesCoefficient[FhInf[g,z],{z,Infinity,0}]}-{((-1)^(g-1) BernoulliB[2g]BernoulliB[2g-2])/(2g(2g-2)Gamma[2g-2+1]) (-100)}


(* ::Subsection::Closed:: *)
(*Conifold point*)


Clear[ruc,xvc,Pzc,Fc];
ruc[g_]:=ruc[g]=#[[1]]->Series[#[[2]],{z,1,2g-2}]&/@(ruleInv/.{X->Xc,A->Ac,B->Bc});
Fc[g_,z_]:=Fc[g,z]=Series[\[Omega]c[0,z],{z,1,2g-2}]^(2g-2) ((1-z)/(5z))^(g-1) (Collect[P[g],variables]/.ruc[g]);
Const[g_]:=Module[{s1z = InverseSeries[Series[zc,{s1,0,2g-2}],z],ans},
	ans=(((-1)^(g-1) BernoulliB[2g])/(2g(2g-2)))/s1z^(2g-2);
	Table[SeriesCoefficient[ans,{z,1,i}],{i,-2g+2,-1}]
];
Matchc[g_]:=Table[SeriesCoefficient[Fc[g,z],{z,1,i}],{i,-2g+2,-1}]-Const[g];


(* ::Subsection::Closed:: *)
(*GV*)


Clear[qz,Diff,MatchGV];
qz=InverseSeries[Series[zInf,{q,0,(degree+1)}],z];
qz[[2]]=Infinity;
Diff[g_]:=Diff[g]=FInf[g,z]-(F[g,q]/.q->qz//Expand)/.GV;
MatchGV[g_]:=Table[SeriesCoefficient[Series[Diff[g],{z,Infinity,degree}],{z,Infinity,i}],{i,1,degree}];


(* ::Subsection::Closed:: *)
(*Castelnuovo*)


gMax[x_]:=(10+5x+x^2)/10;
dMin[g_]:=Ceiling[5/2 (-1+Sqrt[-3+8 g]/Sqrt[5])-1];
Castelnuovo[g_]:=Table[n[g,d],{d,1,dMin[g]}];


(* ::Section:: *)
(*Recursion*)


(* ::Subsection:: *)
(*Boundary conditions*)


Clear[c];
file=dataDirectory<>id<>".mx";
If[FileExistsQ[file],
	Get[file];
	gmin=Sort[Variables[gv]][[1]][[1]];
	,
	gmin=2;
];


step = "";
steps={"bcov","0","Inf","c","gv","cast","sol"};

file=timeDirectory<>id<>".csv";
WriteString[file,"{\n"];
WriteString[file,Join[{"\"g\"","\"tot\""},"\""<>#<>"\""&/@steps]];

Do[
	Print["> genus ",g];
	
	eqns={};

	step="bcov";
	t["bcov",g]=Timing[P[g]][[1]];
	Print[step," took ",t[step,g]," s"];

	step="0";
	t["0",g]=Timing[eqns=Join[eqns,Match0[g]]][[1]];
	Print[step," took ",t[step,g]," s"];

	step="Inf";
	t["Inf",g]=Timing[eqns=Join[eqns,MatchInf[g]]][[1]];
	Print[step," took ",t[step,g]," s"];

	step="c";
	t["c",g]=Timing[eqns=Join[eqns,Matchc[g]]][[1]];
	Print[step," took ",t[step,g]," s"];

	step="gv";
	t["gv",g]=Timing[eqns=Join[eqns,MatchGV[g]]][[1]];
	Print[step," took ",t[step,g]," s"];

	step="cast";
	t["cast",g]=Timing[eqns=Join[eqns,Castelnuovo[g]]][[1]];
	Print[step," took ",t[step,g]," s"];

	step="sol";
	t["sol",g]=Timing[sol=LS[eqns]][[1]];
	Print[step," took ",t[step,g]," s"];

	ToExpression[ToString[#[[1]]]<>":="<>ToString[InputForm[#[[2]]]]]&/@sol;
	
	gv=Table[n[g,d],{d,1,degree},{g,0,genus}]/.GV;
	p[g]=P[g];
	
	file=dataDirectory<>id<>".mx";
	If[FileExistsQ[file],DeleteFile[file<>".bk"];RenameFile[file,file<>".bk"]];
	DumpSave[file,"Global`"];
	
	file=resDirectory<>id<>".mx";
	If[FileExistsQ[file],DeleteFile[file<>".bk"];RenameFile[file,file<>".bk"]];
	DumpSave[file,{c,p,n,GV,gv,t,tt}];
	
	file=timeDirectory<>id<>".csv";
	WriteString[file,",\n"];
	WriteString[file,Join[{g,Table[t[step,g],{step,steps}]//Total},Table[t[step,g],{step,steps}]]];
	
	Print[];
,
	{g,gmin,genus}
];

WriteString[file,"\n}\n"];
