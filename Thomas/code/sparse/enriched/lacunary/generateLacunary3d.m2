--File of examples of sparse systems with Galois group smaller than expected wreath product
restart
loadPackage("DecomposableSparseSystems",Reload=>true)
loadPackage("Polyhedra",Reload=>true)
loadPackage("Monodromy",FileName=>"../../../Monodromy.m2")

i = 0;
j = 0;

MAXI = 10;
MAXJ = 10;

randMat = (n,m,MIN,MAX)->(
    matrix table(n,m,(i,j)->random(MIN,MAX))
    );

while (i<MAXI) do (
    A := apply(3,i->randMat(3,random(3,8),-2,2));
    try(if isDecomposable A then continue) else continue;
    
    P := A/convexHull;
    try(numSolns := sub(volume sum P - volume(P#0+P#1) - volume(P#1+P#2) - volume(P#2+P#0) + volume P#0 + volume P#1 + volume P#2,ZZ)) else continue;
    if (numSolns == 1) or (numSolns > 500) then continue;
    
    while (j<MAXJ) do (
	T := randMat(3,3,-3,3);
	if (abs det T === 1) or (det T === 0) or (abs det T > 40) then continue;
	
	newA := apply(A,S->T*S);
	print("NEW SUPPORT CHOSEN");
	print("SUPPORT #"|toString(MAXJ*i+j)|": "|toString(newA));
	
	try(M := sparseMonodromy newA) else (print("SOLVING ERROR"); j = j+1; continue);
	print("BASE SYSTEM SOLVED");
	if (#M#baseSolutions != (abs det T)*numSolns) then (print("NOT ALL SOLUTIONS FOUND"); continue);
	
	print("TRACKING LOOPS FOR SUPPORTS #"|toString(MAXJ*i+j));
	monodromyLoop(M,10);
	N := size M;
	if (N === 0) then (print("COULDN'T TRACK LOOPS"); j = j+1; continue);
	if (N < (abs det T)^numSolns*numSolns!) then (
	    print("Generating more elements to check")
	    ) else (
	    print("SUPPORTS #"|toString(MAXJ*i+j)|" NOT SPECIAL");
	    j = j+1;
	    continue
	    );
	monodromyLoop(M,50);
	N = size M;
	if (N < (abs det T)^numSolns*numSolns!) then (
	    print("SUPPORTS #"|toString(MAXJ*i+j)|" POSSIBLY SPECIAL AND RECORDED IN "|"lacunaryExamples/"|toString(abs det T)|".m2");
	    print("LATTICE INDEX = "|toString(abs det T));
	    D := first smithNormalForm T;
	    
	    sigmas = select(apply(3,i->D_(i,i)),s->abs s > 1);
	    sigmaStrs = apply(sigmas,s->"Z/"|toString(s)|"Z");
	    print fold(sigmaStrs,(a,b)->a|" x "|b);

	    print("INDEX IN WREATH PRODUCT = "|toString((abs det T)^numSolns*numSolns!/N));
	    myFile := openOutAppend("lacunaryExamples/"|toString(abs det T)|".txt");
	    myFile << toExternalString(newA) << endl << close
	    ) else (
	    print("SUPPORTS #"|toString(MAXJ*i+j)|" NOT SPECIAL")
	    );
	
	j = j+1
	);
    
    j = 0;
    i = i+1
    )
