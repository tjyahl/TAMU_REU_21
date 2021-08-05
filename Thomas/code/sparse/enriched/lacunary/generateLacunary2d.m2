--File of examples of sparse systems with Galois group smaller than expected wreath product
restart
loadPackage("DecomposableSparseSystems",Reload=>true)
loadPackage("Polyhedra",Reload=>true)
loadPackage("Monodromy",FileName=>"../../../Monodromy.m2")

i = 0;
j = 0;

MAXI = 100;
MAXJ = 20;

randMat = (n,m)->(
    matrix table(n,m,(i,j)->random(-3,3))
    );

while (i<MAXI) do (
    A := apply(2,i->randMat(2,random(3,8)));
    try(if isDecomposable A then continue) else continue;
    
    P := A/convexHull;
    numSolns := sub(volume sum P - volume P#0 - volume P#1,ZZ);
    if (numSolns == 1) or (numSolns > 30) then continue;
    
    while (j<MAXJ) do (
	T := randMat(2,2);
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
	    sigma1 = D_(0,0);
	    sigma2 = D_(1,1);
	    if (sigma1 == 1) then (
		print("STABILIZER SUBGROUP = Z/"|toString(sigma2)|"Z")
		) else if (sigma2 == 1) then (
		print("STABILIZER SUBGROUP = Z/"|toString(sigma1)|"Z")
		) else (
		print("STABILIZER SUBGROUP = Z/"|toString(sigma1)|"Z x Z/"|toString(sigma2)|"Z")
		);
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
