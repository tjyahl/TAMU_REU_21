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
-- It checks if A is triangular by calculating the subsystem simpleA and the residual supports and checking for decomposability
    simpleA = {submatrix(A_0,{0,1},),submatrix(A_1,{0,1},)};
    if not try isDecomposable simpleA else true then (
	
    addedSupports =  drop(A,k);
    residualSupports = apply(addedSupports,S->S^(toList(k..numRows S - 1)));
    if not try isDecomposable residualSupports else true then (
    
-- The polytopes defined by the matrices might have dimension < 2, but we need 2-dimensional area for the mixed volume calculation
-- Therefore, simpleVolumeA/B is forced to be 0 when dim < 2
-- Then mixed volume is calculated
    if dim(convexHull(simpleA_0)) == k then simpleVolumeA = volume(convexHull(simpleA_0)) else simpleVolumeA = 0;
    if dim(convexHull(simpleA_1)) == k then simpleVolumeB = volume(convexHull(simpleA_1)) else simpleVolumeB = 0;
    if (simpleVolumeA == 0 and simpleVolumeB == 0) then d=0 else d = volume(convexHull(simpleA_0)+convexHull(simpleA_1)) - simpleVolumeA - simpleVolumeB;
    
-- residualEntries is only 1-dimensional, so the mixed volume is just the max - min
    residualEntries = (entries(residualSupports_0))_0;
    r = max(residualEntries) - min(residualEntries);

-- Since the system should be strictly triangular, the mixed volumes should be > 1   
    if d > 1 and r > 1 then (
    maxGaloisSize = (r!)^d*d!;
    
-- sparseMonodromy selects random coefficients for the system. Some selections make the size smaller than expected
-- Thus, the program attempts 10 monodromy loops for 5 coefficient choices which almost always gives the correct size for systems of this size   
    monodromySize = {};
    for j in 1 .. 5 do(
    	try (M = sparseMonodromy (A,Solver=>M2));
    	if not try M == null else false then(
    	    monodromyLoop(M,10);
    	    monodromySize = append(monodromySize,size(M));
	);
    );
-- outputs.txt is a list of {maxGaloisSize, Size of the actual monodromy, The triangular system} that is only written if there is a notable discrepency    
    if max(monodromySize) <= (maxGaloisSize/2) then ("outputs.txt" << get "outputs.txt" << concatenate {toString(maxGaloisSize),", ",toString(max(monodromySize)),", ",toString(A)} << endl << close);
    );
    );
    );
);

-- This code can be used to double check systems printed to outputs.txt
-- It assumes A is triangular and checks it in more depth

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
