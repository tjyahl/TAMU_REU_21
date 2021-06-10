------------------
--              --
-- Monodromy.m2 --
--              --
------------------
--
-- Thomas Yahl
-- March 22, 2021
--
-- To do: 
--    add checks to sparse/Fano input
--    port monodromy elements to GAP
--    test for imprimitivity
--    certification?
--

newPackage(
    "Monodromy",
    Version=>"0.0.1",
    Authors=>{{
	Name=>"Thomas Yahl",
	Email=>"Thomasjyahl@math.tamu.edu",
	Homepage=>"https://math.tamu.edu/~thomasjyahl"
	}},
    Headline=>"Methods to numerically compute mondromy groups of problems parameterized by a linear family",
    PackageImports=>{"DecomposableSparseSystems","MonodromySolver","NumericalAlgebraicGeometry"},
    PackageExports=>{"NumericalAlgebraicGeometry"},
    DebuggingMode=>true
    )

export{
    --methods
    "mondromy",
    "family",
    "basePoint",
    "params",
    "baseSolutions",
    "group",
    "sparseMonodromy",
    "FanoMonodromy",
    "monodromyLoop",
    "changeBasePoint",
    "saveMonodromy",
    "loadMonodromy",
    "printCycleTally",
    "printCycleTypes",
    "printForGAP",
    "testAlternatingMonodromy",
    "testImprimitiveMonodromy",
    --options
    "baseElements",
    "verbosity",
    "coeffRing",
    "Solver"
    }


-----------------------------------------
------------ Monodromy Class ------------
-----------------------------------------

Monodromy = new Type of MutableHashTable
Monodromy.synonym = "monodromy group"
Monodromy.GlobalAssignHook = globalAssignFunction
Monodromy.GlobalReleaseHook = globalReleaseFunction

monodromy = method(Options=>{baseElements=> tally {}})
monodromy (List,List,List) := Monodromy => o->(F,basePt,baseSolns)->(    
    monod := new Monodromy from {
	family => F,
	basePoint => basePt,
	baseSolutions => baseSolns,
	group => o.baseElements,
	cache => new CacheTable
	};
    monod
    )


sparseMonodromy = method(Options=>{Solver=>PHCPACK})
sparseMonodromy (List) := Monodromy => o->Abullet->(
    n := #Abullet;
    
    x := symbol x;
    a := symbol a;
    t := symbol t;
    
    parms := splice apply(n,i->a_(i,0)..a_(i,numcols(Abullet#i)-1));
    R := CC[parms][x_0..x_(n-1)];
    
    Ffamily := apply(n,k->sum(numcols(Abullet#k),j->a_(k,j)*product(n,i->x_i^((Abullet#k)_(i,j)))));
    
    S := CC[t_1..t_n];
    coeffs := flatten apply(Abullet,A->apply(numcols(A),j->5*(random(CC)-random(CC))));
    phi := map(S,R,toList(t_1..t_n)|coeffs);
    
    F := Ffamily/phi;
    
    solns := solveDecomposableSystem(F,Software=>o.Solver);
    
    M := monodromy(rewriteEqns(Ffamily),coeffs,solns);
    M
    )


FanoMonodromy = method()
FanoMonodromy (ZZ,ZZ,Sequence) := Monodromy => (k,m,n)->(
    --k = 4 --linear space dimension
    --m = 14 --affine/projective space dimension
    --n = (2,3) --degrees of hyperplanes
    --r number of hyperplanes
    --need (m-k)*(k+1)=sum_i binomial(n_i+k,k).
    r := #n;
    
    x := symbol x;
    w := symbol w;
    t := symbol t;
    y := symbol y;
    
    parms := splice apply(r,i->w_(i,1)..w_(i,binomial(n_i+m,m)));
    variables := toList(x_(1,1)..x_(m-k,k+1));
    R := QQ[parms][variables];
    S := R[t_1..t_k];
    T := QQ[parms][y_1..y_m];
    
    phi := map(S,T,apply(m-k,j->sum(toList(x_(j+1,1)..x_(j+1,k)),toList(t_1..t_k),(z,v)->z*v)+x_(j+1,k+1))|toList(t_1..t_k));
    
    hypEqns := apply(r,i->basis(0,n_i,module T)*(transpose matrix {toList(w_(i,1)..w_(i,binomial(n_i+m,m)))}));
    subsEqns := apply(hypEqns,f->(phi(f))_(0,0));
    coeffs := (flatten apply(subsEqns,f->flatten entries last coefficients f))/(c->sub(c,R));
    coeffs = rewriteEqns(coeffs);
    
    --MONODROMY SOLVER NOT WORKING RN :(
    (basePt,baseSolns) := solveFamily(polySystem coeffs);
    basePt = coordinates basePt;
    baseSolns = baseSolns/coordinates;
    
    M := monodromy(coeffs,basePt,baseSolns);
    M
    )


--needs fixing
isWellDefined (Monodromy) := Boolean => M->(
    eps := 1e-8;
    F := M#family;
    R := ring F#0;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
    n := #(first baseSolns);
    
    MisWellDefined := false;
    
    MisWellDefined
    )


monodromyLoop = method(Options=>{verbosity=>0})
monodromyLoop (Monodromy,List) := Monodromy => o-> (M,pts)->(
    F := M#family;
    R := ring F#0;
    k := numgens R;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
	    
    t := symbol t;
    mon := symbol mon;
    S := CC[toList(mon_1..mon_k)];
    
    pts = {basePt}|pts|{basePt};
    ptsMaps := apply(pts,pt->map(S,R,toList(mon_1..mon_k)|pt));
    trackedSolns := baseSolns;
    
    for i in 0..(#pts-2) do (
	trackedSolns = track(F/ptsMaps#i,F/ptsMaps#(i+1),trackedSolns,stepIncreaseFactor=>2,tStep=>.01,maxCorrSteps=>1,CorrectorTolerance=>1e-8)/coordinates
	);
	    
    baseSolnsRounded := apply(baseSolns,s->rounded(s,5));
    trackedSolnsRounded := apply(trackedSolns,s->rounded(s,5));
    perm := try(permutation(baseSolnsRounded,trackedSolnsRounded)) else error "Error";
    if (class perm===String) then (
	if (o.verbosity>0) then print("error: no permutation between tracked solutions")
	) else (
        T := tally {perm};
        M#group = M#group + T
	);
    M
    )

monodromyLoop (Monodromy) := Monodromy => o-> M->(
    n := numgens coefficientRing ring (M#family#0);
    monodromyLoop(M,apply(2,i->apply(n,j->10*(random(CC)-random(CC)))),o)
    )

monodromyLoop (Monodromy,ZZ) := Monodromy => o->(M,n)->(
    for i in 1..n do (
	try(monodromyLoop(M,o)) else if (o.verbosity>0) then print("error: path tracking failed");
	if (o.verbosity>0) then print("Loop "|toString(i)|" complete")
	);
    M
    )


changeBasePoint = method()
changeBasePoint (Monodromy,List) := Monodromy => (M,bp)->(
    F := M#family;
    R := ring F#0;
    k := numgens R;
    basePt := M#basePoint;
    baseSolns := M#baseSolutions;
    grp := M#group;
    
    mon := symbol mon;
    S := CC[mon_1..mon_k];
    
    startMap := map(S,R,toList(mon_1..mon_k)|basePt);
    endMap := map(S,R,toList(mon_1..mon_k)|bp);
    
    newBaseSolns := track(F/startMap,F/endMap,baseSolns,gamma=>random(CC),stepIncreaseFactor=>2,tStep=>.01,maxCorrSteps=>1,CorrectorTolerance=>1e-8)/coordinates;
    
    monodromy(F,bp,newBaseSolns,baseElements=>grp)
    )


saveMonodromy = method()
saveMonodromy (Monodromy,String) := Nothing => (M,s)->(
    F := M#family;
    R := ring F#0;
    C := coefficientRing R;
    s << toString gens R << endl << toString gens C << endl << toExternalString (M#family) << endl << toExternalString (M#basePoint) << endl << toExternalString (M#baseSolutions) << endl << toString (M#group) << endl << toString (apply(keys M#cache,h->h=>M#cache#h)) << close
    )


loadMonodromy = method(Options=>{coeffRing=>CC})
loadMonodromy (String) := Monodromy => o->s->(
    L := lines get s;
    vars := value L#0;
    parms := value L#1;
    R := o.coeffRing[parms][vars];
    F := value L#2;
    basePt := apply(value L#3,z->toCC(defaultPrecision,z));
    baseSolns := apply(value L#4,L->apply(L,z->toCC(defaultPrecision,z)));
    grp := value L#5;
    monodromyCache := value L#6;
    
    M := monodromy(F,basePt,baseSolns,baseElements=>grp);
    M#cache = new CacheTable from monodromyCache;
    M
    )



--------------------------------------------------------
------------ Utility Methods (Not exported) ------------
--------------------------------------------------------

--rewrites equations in k[parms][variables] without doublely-indexed variables.
rewriteEqns = method()
rewriteEqns (List) := Sequence => F->(
    try uniform(F) else error "Error: Expected polynomials of the same ring";
    R := ring F#0;
    n := numgens R;
    m := numgens coefficientRing R;
    
    x := symbol x;
    a := symbol a;
    
    newR := CC[a_1..a_m][x_1..x_n];
    phi := map(newR,R,toList(x_1..x_n)|toList(a_1..a_m));
    F/phi
    )


rounded = method()
rounded (List,ZZ) := List => (Z,n)->(
    apply(Z,z->round(10^n*z)/10^n)
    )


setDifference = method()
setDifference (List,List) := List => (A,B)->(
    L := select(A,x->not member(x,B));
    L
    )


permutation = method()
permutation (List,List) := List => (A,B)->(
    if not (A === unique A and B === unique B) then error "-- permutation: list has repeated elements";
    if any(A,a->not member(a,B)) or any(B,b->not member(b,A)) then error "-- permutation: no bijection between lists";

    apply(A,a->position(B,b->b==a)+1)
    )


--this is terribly written.                                                                                                           
printPermutation = method()
printPermutation (List) := String => (P)->(
    if not (P === unique P and all(P,n->member(n,toList(1 .. #P)))) then error "--printPermutation: not permutation";

    cycles := {};
    L := toList(1 .. #P);
    while (#L > 0) do (
        i := first L;
        s := "(" | toString(i);
        while (P#(i-1) != first L) do (
            L = delete(P#(i-1),L);
            i = P#(i-1);
            s = s | "," | toString(i);
            );
        L = drop(L,1);
        s = s | ")";
        if (#separate(",",s)>1) then (
	    cycles = append(cycles,s)
	    )
        );
    if (#cycles > 0) then (
	fold(cycles,(sigma,tau)->sigma|tau)
	) else (
	"()"
	)
    )


numCycles = method()
numCycles (List) := ZZ => P->(
    if not (P === unique P and all(P,n->member(n,toList(1 .. #P)))) then error "--numCycles: not permutation";
    cycles := 0;
    L := toList(1 .. #P);
    while (#L>0) do (
	i := first L;
	while (P#(i-1) != first L) do (
	    L = delete(P#(i-1),L);
            i = P#(i-1)
	    );
	L = drop(L,1);
	cycles = cycles + 1
	);
    cycles
    )


cycleType = method()
cycleType (List) := ZZ => P->(
    if not (P === unique P and all(P,n->member(n,toList(1 .. #P)))) then error "--cycleType: not permutation";
    cType := {};
    L := toList(1 .. #P);
    while (#L>0) do (
	len := 1;
	i := first L;
	while (P#(i-1) != first L) do (
	    L = delete(P#(i-1),L);
            i = P#(i-1);
	    len = len + 1
	    );
	L = drop(L,1);
	cType = append(cType,len)
	);
    sort cType
    )


--------------------------------------
------------ Main Methods ------------
--------------------------------------

printCycleTally = method()
printCycleTally (Monodromy) := Nothing => M->(
    T := M#group;
    scan(select(keys T,k->class k === List),sigma-> print(printPermutation(sigma) | ": " | toString(T#sigma)))
    )


printCycleTypes = method()
printCycleTypes (Monodromy) := Nothing => M->(
    T := elements M#group;
    tally apply(T,cycleType)
    )


printForGAP = method()
printForGAP (Monodromy) := String => M->(
    T := M#group;
    if (#(keys T) > 0) then (
	S := apply(keys T,sigma->printPermutation(sigma));
    	"Group("|fold(S,(s,t)->s|","|t)|")"
	) else (
	"Group()"
	)
    )


testAlternatingMonodromy = method()
testAlternatingMonodromy (Monodromy) := Boolean => M->(
    T := M#group;
    if (#(keys T) > 0) then (
    	m := length (first keys T);
    	all(keys T,L-> (m - numCycles L) % 2 === 0)
	) else (
	false
	)
    )


testImprimitiveMonodromy = method()
testImprimitiveMonodromy (Monodromy) := Boolean => M->(
    
    )







end

--loads old form into new form.
convert = method()
convert (String) := Nothing => s->(
    L := lines get s;
    parms := value L#1;
    vars := setDifference(value L#0,parms);
    R := CC[parms][vars];
    F := value L#2;
    basePt := value L#3;
    baseSolns := value L#4;
    grp := value L#5;
    monodromyCache := value L#6;
    
    M := monodromy(F,basePt,baseSolns,baseElements=>grp);
    M#cache = new CacheTable from monodromyCache;
    saveMonodromy(M,s)
    )

--To check sparse family (easily modified for other families.
F = matrix{F#family}
R = ring F
maxN = max apply(M#baseSolutions,s->(
    f = map(CC,R,s|M#basePoint);
    norm f(F)
))
