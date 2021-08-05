--file to inspect the subgroup of the Galois group that fixes the blocks.
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2")
loadPackage("Groups",FileName=>"../Groups.m2")

--Interesting ones.
----3.txt: 
----
----4.txt: 0,6 (hard) 
----
----6.txt: 0,2,19 (ind 3),28,39,55
----
L = (lines get "enriched/lacunary/lacunaryExamples/6.txt")/value
A = L#0

M = sparseMonodromy A
monodromyLoop(M,20,Verbosity=>2)

n = #A
(D,P,Q) = smithNormalForm fold(A,(a,b)->a|b)
elemDivs = apply(n,i->D_(i,i))
latInd = product elemDivs
numBaseSolns = sub((#M#baseSolutions)/latInd,ZZ)
maxGaloisSize = latInd^numBaseSolns*numBaseSolns!
maxGaloisSize/(size M)

--
inFile = temporaryFileName()
outFile = temporaryFileName()
grp = printForGAP M
setStr = "["|fold((1..#M#baseSolutions)/toString,(a,b)->a|","|b)|"]"
inFile << "G := "|grp|";; B := MaximalBlocks(G,"|setStr|");; phi := ActionHomomorphism(G,B,OnSets);; K := Kernel(phi);; MinimalGeneratingSet(K);" << close
run("(cat "|inFile|" | gap -q) 1>"|outFile)
K = get outFile
removeFile inFile
removeFile outFile





end

G = groupInit apply(keys M#group,p->hashTable transpose {toList(0..#M#baseSolutions-1),p})

--blocks should be ordered as long as you use DecomposableSparseSystems.m2 for solving
----Otherwise, monomial map projecting onto base solutions is (inverse P)*D_truncated
blockSets = apply(numBaseSolns,i->set toList(latInd*i..latInd*(i+1)-1))

kerTransversal = new MutableHashTable from {blockSets=>hashTable apply(#M#baseSolutions,i->i=>i)}
kerGens = set {}
active = {blockSets}
while (#active > 0) do (
    B = first active;
    for g in G#genSet do (
	gB = apply(B,s->set apply(elements s,i->g#i));
	if not kerTransversal#?(gB) then (
	    kerTransversal#(gB) = g*kerTransversal#B;
	    active = append(active,gB)
	    ) else (
	    kerGens = kerGens + set {(inverse kerTransversal#(gB))*g*kerTransversal#B}
	    )
	);
    active = drop(active,1);
    print("active size = "|toString(#active))
    )
kerTransversal = hashTable pairs kerTransversal
kerG = groupInit elements kerGens

--checks
all(values kerTransversal,h->member(h,G))
all(pairs kerTransversal,(B,g)->apply(blockSets,s->set apply(elements s,i->g#i))===B)

all(kerG#genSet,g->member(g,G))
all(kerG#genSet,g->apply(blockSets,s->set apply(elements s,i->g#i))===blockSets)

(kerG#genSet)/printCycles



end


--code to find an instance with specific index in wreath product
restart
loadPackage("Monodromy",FileName=>"../Monodromy.m2")
loadPackage("Groups",FileName=>"../Groups.m2")
L = (lines get "enriched/lacunary/lacunaryExamples/6.txt")/value

all(L,A->(
	print("starting "|toString(position(L,B->B==A)));
	try(M = sparseMonodromy A) else (print("solving failed?"); return true);
	try(scan(10,i->monodromyLoop(M,Verbosity=>2))) then (
	    n = #A;
	    (D,P,Q) = smithNormalForm fold(A,(a,b)->a|b);
	    elemDivs = apply(n,i->D_(i,i));
	    latInd = product elemDivs;
	    numBaseSolns = sub((#M#baseSolutions)/latInd,ZZ);
	    maxGaloisSize = latInd^numBaseSolns*numBaseSolns!;
	    ind = sub(maxGaloisSize/(size M),ZZ);
	    print("INDEX = "|toString(ind));
	    print("NUMBASESOLNS = "|toString(numBaseSolns));
	    numBaseSolns % ind == 0
	    ) else (
	    print("tracking failed");
	    true
	    )
	)
    )



