restart
loadPackage("DecomposableSparseSystems")
loadPackage("Monodromy",FileName=>"Monodromy.m2")

k = 2
outputs = {}
outputMatrices = {}
for i in 0 .. 50 do(
-- First, a random set of 3 matrices is generated with the 3rd row of the first 2 all 0's
    termNumber = {random(2,5),random(2,5),random(2,5)};
    A = {mutableMatrix(ZZ,k+1,termNumber_0),mutableMatrix(ZZ,k+1,termNumber_1),mutableMatrix(ZZ,k+1,termNumber_2)};
-- A is 3 zero matrices each of size 3 by 2-5
-- Now it's filled with entries from 0 to 2 (I chose 2 so the solutions wouldn't take too long to compute)
    for i in (0,0) .. (k-1,termNumber_0-1) do(
    	iRow = i_0;
    	iColumn = i_1;
    	A_0_(iRow,iColumn) = random(0,2);
    );
    for i in (0,0) .. (k-1,termNumber_1-1) do(
    	iRow = i_0;
    	iColumn = i_1;
    	A_1_(iRow,iColumn) = random(0,2);
    );
-- Only this 3rd matrix has all 3 rows filled in
    for i in (0,0) .. (k,termNumber_2-1) do(
    	iRow = i_0;
    	iColumn = i_1;
    	A_2_(iRow,iColumn) = random(0,2);
    );
    A = {matrix(A_0),matrix(A_1),matrix(A_2)};
    
    outputMatrices = append(outputMatrices,A);
    
-- This next part comes mostly from Thomas's code.
-- It calculates the nontriangular subsystem simpleA and the residual supports, then checks for decomposability
-- Finally, it calculates the number of solutions and the maximum size of the resultant Galois group
    simpleA = {submatrix(A_0,{0,1},),submatrix(A_1,{0,1},)};
    if not try isDecomposable simpleA else true then (
	
    addedSupports =  drop(A,k);
    residualSupports = apply(addedSupports,S->S^(toList(k..numRows S - 1)));
    if not try isDecomposable residualSupports else true then (
	
    (F,sols) = solveDecomposableSystem(simpleA,,Software=>M2);
    d = #sols;
    (F,sols) = solveDecomposableSystem(residualSupports,,Software=>M2);
    r = #sols;
    
    maxGaloisSize = (r!)^d*d!;
-- I throw out sizes larger than 1000 for computation time. About 1/3 random systems have size >1000
    if maxGaloisSize < 1000 then (
    M = sparseMonodromy (A,Solver=>M2);
-- 5*maxGaloisSize loops usually places the group size within 90% of its true value
    monodromyLoop(M,5*maxGaloisSize);
    queuedOutput = {i,maxGaloisSize,#(M#group)};
    outputs = append(outputs,{i,maxGaloisSize,#(M#group)});
    ) else outputs = append(outputs,{i,maxGaloisSize,-1});
    );
    );
);
outputs
-- Outputs gives a list of {index of system in outputMatrices, maxGaloisSize, monodromy output if maxGaloisSize<1000, -1 elsewise}
