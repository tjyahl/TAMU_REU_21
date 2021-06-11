--Package for studying properties of groups.

newPackage(
    "Groups",
    Version=>"0.0.1",
    Headline=>"Methods for studying properties of groups",
    Authors=>{
	{Name=>"Thomas Yahl",Email=>"thomasjyahl@tamu.edu"}
	},
    DebuggingMode=>true
    )

export{
    --constructors
    "group",
    "coset",
    "cyclicGroup",
    "symmetricGroup",
    "alternatingGroup",
    "wreathProduct",
    "groupIntersection",
    "symmetricEmbedding",
    "actionOnCosets",
    --
    "identityElement",
    "printCycles",
    "numCycles",
    "cycleType",
    "alternatingSubgroup",
    "transitive",
    "imprimitive",
    --
    "orbit",
    "cosetOrbit",
    "orbitPartition",
    "stabilizer",
    "computeBSGS",
    "pruneGens",
    "subgroup",
    "generateElements",
    "blocks",
    --options
    "genSet",
    "baseSet",
    "representative",
    "baseGroup",
    "BSGS"
    }

---------------------
----Group Objects----
---------------------

Group = new Type of HashTable
Group.synonym = "group"
Group.GlobalAssignHook = globalAssignFunction
Group.GlobalReleaseHook = globalReleaseFunction

Coset = new Type of HashTable
Coset.synonym = "coset"
Coset.GlobalAssignHook = globalAssignFunction
Coset.GlobalReleaseHook = globalReleaseFunction


group = method()
group (Set,List) := Group => (X,perms)->(
    G := new Group from {
	baseSet => X,
	genSet => perms,
	cache => new CacheTable
	};
    G
    )

group (List) := Group => perms->(
    X := set keys first perms;
    if all(perms,g->(set keys g === set values g) and (set keys g === X)) then (
	group(X,perms)
	) else (
	error "--group: not all permutations of the same set"
	)
    )


coset = method()
coset (HashTable,Group) := Coset => (g,H)->(
    gH := new Coset from {
	representative => g,
	baseGroup => H,
	cache => new CacheTable
	};
    gH
    )


isWellDefined (Group) := Boolean => G->(
    X := G#baseSet;
    all(G#genSet,g->(set keys g === set values g) and (set keys g === X))
    )

isWellDefined (Coset) := Boolean => gH->(
    g := gH#representative;
    H := gH#baseGroup;
    (set keys g == H#baseSet)
    )


------------------------------------------
----Group Constructions and Operations----
------------------------------------------

cyclicGroup = method()
cyclicGroup (ZZ) := Group => n->(
    if (n>1) then (
	sigma := hashTable(apply(n-1,i->i+1=>i+2)|{n=>1});
	group({sigma})
	) else (
	group({hashTable {1=>1}})
	)
    )


symmetricGroup = method()
symmetricGroup (ZZ) := Group => n->(
    if (n>1) then (
	sigma := hashTable({1=>2,2=>1}|apply(toList(3..n),i->i=>i));
        tau := hashTable(apply(n-1,i->i+1=>i+2)|{n=>1});
	group(unique {sigma,tau})
	) else (
	group({hashTable {1=>1}})
	)
    )


alternatingGroup = method()
alternatingGroup (ZZ) := Group => n->(
    Sn := symmetricGroup(n);
    alternatingGroup(Sn)
    )

alternatingGroup (Group) := Group => G->(
    if alternatingSubgroup(G) then return G;
    X := elements G#baseSet;
    n := #X;
    sigma := first select(G#genSet,g->(n - numCycles g)%2 != 0);
    A := group unique flatten apply(G#genSet,g->(
	    if ((n - numCycles g) % 2 == 0) then (
		{g,sigma*g*sigma}
		) else (
		{sigma*g,g*sigma}
		)
	    )
	);
    A
    )


Group==Group := Boolean => (G,H)->(
    subgroup(G,H) and subgroup(H,G)
    )

Coset==Coset := Boolean => (aH,bH)->(
    H := aH#baseGroup;
    a := aH#representative;
    (bH#baseGroup == H) and member(a,bH)
    )


Group*Group := Group => (G,H)->(
    X := G#baseSet;
    Y := H#baseSet;
    prodBase := elements(X**Y);
    GxId := apply(G#genSet,g->combine(g,identityElement(H),splice,splice,first));
    IdxH := apply(H#genSet,h->combine(identityElement(G),h,splice,splice,first));
    GxH := group(unique(GxId|IdxH));
    GxH
    )


Group^ZZ := Group => (G,n)->(
    if (n <= 0) then (
	error "--power: expected positive power"
	) else (
	product(toList(n:G))
	)
    )


wreathProduct = method()
wreathProduct (Group,Group) := Group => (G,H)->(
    X := G#baseSet;
    Y := H#baseSet;
    prodBase := elements(X**Y);
    blockPerms := apply(H#genSet,h->combine(identityElement(G),h,identity,identity,first));
    elemPerms := flatten apply(elements Y,i->apply(G#genSet,g->hashTable apply(prodBase,p->p=>(
		    if (last p == i) then (
			(g#(first p),i)
			) else (
			p
			)
		    )
		)
	    )
	);
    XwrY := group(unique(blockPerms|elemPerms));
    XwrY
    )


groupIntersection = method()
groupIntersection (Group,Group) := Group => (G,H)->(
    C := cosetOrbit(coset(identityElement(H),H),G);
    GcapH := unique flatten table(C,G#genSet,(aH,g)->(
	    gaH := first select(C,bH->bH==g*aH);
	    (inverse (gaH#representative))*g*(aH#representative)
	    )
	);
    group(GcapH)
    )


symmetricEmbedding = method()
symmetricEmbedding (Group) := Group => G->(
    
    )


actionOnCosets = method()
actionOnCosets (Group,Group) := Group => (G,H)->(
    
    )


---------------------
----Minor Methods----
---------------------

identityElement = method()
identityElement (Group) := HashTable => G->(
    X := elements G#baseSet;
    hashTable apply(X,i->i=>i)
    )


HashTable*HashTable := HashTable => (g,h)->(
    applyValues(h,i->g#i)
    )

HashTable*Coset := Coset => (a,bH)->(
    b := bH#representative;
    H := bH#baseGroup;
    coset(a*b,H)
    )

inverse (HashTable) := HashTable => g->(
    applyPairs(g,reverse)
    )


printCycles = method()
printCycles (HashTable) := String => g->(
    X := keys g;
    
    cycles := {};
    while (#X>0) do (
        i := first X;
	singleton := true;
        s := "(" | toString(i);
        while (g#i != first X) do (
            X = delete(g#i,X);
            i = g#i;
            s = s | "," | toString(i);
	    singleton = false
            );
        X = drop(X,1);
        s = s | ")";
        if not singleton then (
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
numCycles (HashTable) := List => g->(
    X := keys g;

    cycles := 0;
    while (#X>0) do (
	i := first X;
	while (g#i != first X) do (
	    X = delete(g#i,X);
            i = g#i
	    );
	X = drop(X,1);
	cycles = cycles + 1
	);
    cycles
    )


cycleType = method()
cycleType (HashTable) := List => g->(
    X := keys g;
    
    cType := {};
    while (#X>0) do (
	len := 1;
	i := first X;
	while (g#i != first X) do (
	    X = delete(g#i,X);
            i = g#i;
	    len = len + 1
	    );
	X = drop(X,1);
	cType = append(cType,len)
	);
    sort cType
    )

cycleType (List) := List => L->(
    apply(L,cycleType)
    )

cycleType (Group) := List => G->(
    cycleType generateElements G
    )


alternatingSubgroup = method()
alternatingSubgroup (Group) := Boolean => G->(
    n := #(G#baseSet);
    all(G#genSet,g-> (n - numCycles g) % 2 === 0)
    )


transitive = method()
transitive (Group,ZZ) := Boolean => (G,k)->(
    X := elements G#baseSet;
    if (k<=#X) then (
	distinctTuples := flatten apply(subsets(X,k),permutations);
	perms := apply(G#genSet,g->hashTable apply(distinctTuples,T->T=>apply(T,i->g#i)));
	Gbar := group(perms);
	transitive(Gbar)
	) else (
	false
	)
    )

transitive (Group) := Boolean => G->(
    X := G#baseSet;
    x := first elements X;
    first orbit(x,G) === X
    )


imprimitive = method()
imprimitive (Group) := Sequence => G->(
    
    )


--------------------
----Main Methods----
--------------------

orbit = method()
orbit (Thing,Group) := Sequence => (k,G)->(
    X := elements G#baseSet;
    if not member(k,X) then error "--orbit: point not in set.";
    O := new MutableHashTable from {k=>identityElement(G)};
    L := {k};
    while (#L>0) do (
	v := first L;
	for g in G#genSet do (
	    if not O#?(g#v) then (
		O#(g#v) = g*(O#v);
		L = append(L,g#v)
		)
	    );
	L = drop(L,1)
	);
    O = hashTable pairs O;
    (set keys O,O)
    )


cosetOrbit = method()
cosetOrbit (Coset,Group) := List => (eH,G)->(
    if not (G#baseSet === (eH#baseGroup)#baseSet) then error "--cosetOrbit: expected groups acting on the same set";
    C := {eH};
    L := C;
    while (#L>0) do (
	bH := first L;
	for a in G#genSet do (
	    if not any(C,gH->a*bH == gH) then (
		C = append(C,a*bH);
		L = append(L,a*bH)
		)
	    );
	L = drop(L,1)
	);
    C
    )


orbitPartition = method()
orbitPartition (Group) := List => G->(
    X := G#baseSet;
    P := {};
    while (#X>0) do (
	O := first orbit(first elements X,G);
	P = append(P,O);
	X = X-O
	);
    P
    )


stabilizer = method()
stabilizer (Thing,Group) := Group => (n,G)->(
    X := G#baseSet;
    if not member(n,X) then error "--stabilizer: point not in set";	
    O := last orbit(n,G);
    stab := unique flatten table(keys O,G#genSet,(v,g)->(inverse O#(g#v))*g*(O#v));
    group(stab)
    )


computeBSGS = method()
computeBSGS (Group) := Group => G->(
    X := elements G#baseSet;
    H := G;
    stabs := {G};
    for x in X do (
	H = stabilizer(x,H);
	stabs = append(stabs,H)
	);
    stabs = drop(stabs,-1);
    G#cache#BSGS = transpose {X,stabs};
    G
    )


member (HashTable,Group) := Boolean => (g,G)->(
    X := G#baseSet;
    if not (set keys g === X) then error "--member: expected element and group to act on same set";
    if not G#cache#?BSGS then computeBSGS(G);
    h := g;
    for e in G#cache#BSGS do (
	O := last orbit(first e,last e);
	if not member(h#(first e),keys O) then return false;
	h = (inverse O#(h#(first e)))*h
	);
    true
    )

member (HashTable,Coset) := Boolean => (a,bH)->(
    b := bH#representative;
    H := bH#baseGroup;
    member((inverse a)*b,H)
    )


subgroup = method()
subgroup (Group,Group) := Boolean => (H,G)->(
    all(H#genSet,h->member(h,G))
    )


size (Group) := ZZ => G->(
    if not G#cache#?BSGS then computeBSGS(G);
    product(G#cache#BSGS,v->#(first orbit(first v,last v)))
    )


pruneGens = method()
pruneGens (Group) := Group => G->(
    G
    )


generateElements = method()
generateElements (Group) := Group => G->(
    
    )


blocks = method()
blocks (Group) := List => G->(
    
    )



end
restart
loadPackage("Groups")
G = group({hashTable{1=>2,2=>1,3=>3,4=>4},hashTable{1=>2,2=>1,3=>4,4=>3}})
isWellDefined G
S = symmetricGroup(4)
subgroup(G,S)
listCosets(G,S)
transitive G
orbit(1,G)
orbit(3,G)
orbitPartition(G)
S = symmetricGroup(5)
orbit(1,S)
stabilizer(5,S)
apply((stabilizer(4,S))#genSet,g->printCycles g)
computeBSGS(S)
S#cache#BSGS
size(S)
g = hashTable transpose {toList(1..5),random(toList(1..5))}
member(g,S)




product (HashTable,HashTable) := HashTable => (g,h)->(
    applyValues(h,i->g#i)
    )

product (HashTable,Coset) := Coset => (a,bH)->(
    b := bH#representative;
    H := bH#baseGroup;
    coset(a*b,H)
    )
