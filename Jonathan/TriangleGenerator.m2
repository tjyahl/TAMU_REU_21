restart
loadPackage("Polyhedra")
loadPackage("DecomposableSparseSystems")
loadPackage("Monodromy",FileName=>"Monodromy.m2")

k = 2
for i in 0 .. 500 do(
-- First, a random set of 3 matrices is generated with the 3rd row of the first 2 all 0's
    termNumber = {random(2,5),random(2,5),random(2,5)};
    A = {mutableMatrix(ZZ,k+1,termNumber_0),mutableMatrix(ZZ,k+1,termNumber_1),mutableMatrix(ZZ,k+1,termNumber_2)};
-- A is 3 zero matrices each of size 3 by 2-5
-- Now it's filled with entries from 0 to 2 (I chose 2 so the solutions wouldn't take too long to compute)
    for i in (0,0) .. (k-1,termNumber_0-1) do(
    	iRow = i_0;
    	iColumn = i_1;
    	A_0_(iRow,iColumn) = random(0,3);
    );
    for i in (0,0) .. (k-1,termNumber_1-1) do(
    	iRow = i_0;
    	iColumn = i_1;
    	A_1_(iRow,iColumn) = random(0,3);
    );
-- Only this 3rd matrix has all 3 rows filled in
    for i in (0,0) .. (k,termNumber_2-1) do(
    	iRow = i_0;
    	iColumn = i_1;
    	A_2_(iRow,iColumn) = random(0,3);
    );
    A = {matrix(A_0),matrix(A_1),matrix(A_2)};
    
-- This next part comes mostly from Thomas's code.
-- It calculates the nontriangular subsystem simpleA and the residual supports, then checks for decomposability
-- Finally, it calculates the number of solutions and the maximum size of the resultant Galois group
    simpleA = {submatrix(A_0,{0,1},),submatrix(A_1,{0,1},)};
    if not try isDecomposable simpleA else true then (
	
    addedSupports =  drop(A,k);
    residualSupports = apply(addedSupports,S->S^(toList(k..numRows S - 1)));
    if not try isDecomposable residualSupports else true then (
    
    d = volume(convexHull(simpleA_0)+convexHull(simpleA_1)) - volume(convexHull(simpleA_0)) - volume(convexHull(simpleA_1));
    residualEntries = (entries(residualSupports_0))_0;
    r = max(residualEntries) - min(residualEntries);
    
    maxGaloisSize = (r!)^d*d!;
    if d > 1 and r > 1 then (
    monodromySize = {};
    for j in 1 .. 5 do(
    	try (M = sparseMonodromy (A,Solver=>M2));
    	if not try M == null else false then(
    	    monodromyLoop(M,10);
    	    monodromySize = append(monodromySize,size(M));
	);
    );    
    if max(monodromySize) <= (maxGaloisSize/2) then ("outputs.txt" << get "outputs.txt" << concatenate {toString(maxGaloisSize),", ",toString(max(monodromySize)),", ",toString(i),", ",toString(A)} << endl << close);
    );
    );
    );
);

-- Outputs gives a list of {index of system in outputMatrices, maxGaloisSize, monodromy output if maxGaloisSize<1000, -1 elsewise}

simpleA = {submatrix(A_0,{0,1},),submatrix(A_1,{0,1},)};
addedSupports =  drop(A,k);
residualSupports = apply(addedSupports,S->S^(toList(k..numRows S - 1)));
d = volume(convexHull(simpleA_0)+convexHull(simpleA_1)) - volume(convexHull(simpleA_0)) - volume(convexHull(simpleA_1));
residualEntries = (entries(residualSupports_0))_0;
r = max(residualEntries) - min(residualEntries);
maxGaloisSize = (r!)^d*d!;
monodromySize = {};
for j in 1 .. 5 do(
    try (M = sparseMonodromy (A,Solver=>M2));
    if not try M == null else false then(
    	monodromyLoop(M,10);
    	monodromySize = append(monodromySize,size(M));
    );
);
max(monodromySize) - maxGaloisSize
