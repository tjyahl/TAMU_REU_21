restart
load("../functions.m2")

-- useful function definitions:
polytopeVolume = method();
polytopeVolume (Thing, ZZ) := (f, actions) -> (
    local E; local L; local K; local l;
    E = exponents f;
    L = new List;
    for i from 0 to length E - 1 do (
	l = new List;
	for j from 0 to actions do (
	    l = append(l, E_i_j);
	);
	L = append(L, l);
    );
    L = unique L;
    
    K = transpose matrix L;
    
    return (actions + 1)! * volume convexHull K
    
)



-- Takes an integer a and returns the nth bit in a's bit representation
getBit = method();
getBit (ZZ, ZZ) := (a, n) -> (
    return (a >> n) & 1;
)

-- Takes two integers whose bit representations give the edges that are included in the represented 
-- graphs and returns true if the first graph is a subgraphs of the second.
isSubgraph = method();
isSubgraph (ZZ, ZZ) := (a, b) -> (
    if (a & b) == a then return true else return false;
)

-- performs the steps to remove the edge variables and the inverse variables from the given
-- polynomial. Need to be concious of any change in rings (should move to one without edges)
getExponents = method();
getExponents (Thing, ZZ) := (f, actions) -> (
    local E; local L; local K; local l;
    E = exponents f;
    L = new List;
    for i from 0 to length E - 1 do (
	l = new List;
	for j from 0 to actions do (
	    l = append(l, E_i_j);
	);
	L = append(L, l);
    );
    L = unique L;
    return L
)

-- Used to generate the next combination of permutations to inspect
nextPermutation = method();
nextPermutation (MutableList, ZZ, ZZ) := (permutation, actions, maxNum) -> (
    
    for i from 1 to actions do (
	permutation#(actions - i) = permutation#(actions - i) + 1;
	if permutation#(actions - i) >= maxNum and i != actions then (
	    continue;
	);
    	if actions - i == 0 and permutation#0 >= maxNum then return false;
    	for j from actions - i to actions - 1 do (
	    permutation#j = permutation#(actions - i);
	);
    	return permutation;
    );
)

-- Takes a list containing permutations and a list of indices and returns the edge numbers associated with
-- the graph that has the given structure. For example, given a graph with 2 actions and 2 vertices,
-- the input could be permutation = {0,1} and permutes = {{1,2}, {2,1}} in which case, the graph 
-- would contain the edges connecting v_1 to x_1*v_1, v_2 to x_1*v_2, v_1 to x_2*v_2, v_2 to
-- x_2*v_1, and v_1 to v_2. The edges are assigned numbers by the formula v_i to x_k v_j is given by 
-- e_((k-1)*n^2 + n*(i-1) + j) where n is the number of vertices in the fundamental domain. The last
-- C(n,2) edges are included as they make up the complete graph in the fundamental domain.
generateEdges = method();
generateEdges (List, List) := (permutation, permutes) -> (
    local edges; local loopMax;
    loopMax = (length permutation - 1)*(length permutes#0)^2 + (length permutes#0)*(length permutes#0 - 1) + length permutes#0;
    edges = new List;
    for i from 0 to length permutation - 1 do (
	for j from 0 to length permutes#(permutation#i) - 1 do (
	    edges = append(edges, i*(length permutes#(permutation#i))^2 + (length permutes#(permutation#i))*j + permutes#(permutation#i)#j);
	);
    );
    
    for i from loopMax + 1 to loopMax + binomial(length permutes#0, 2) do (
	edges = append(edges, i);
    );
    return edges;
)

-- Takes the current permutation for the actions, the list of permutations, and the number of vertices in the fundamental domain
-- and checks that no two permutations in the current list map the same vertex to itself. This would create
-- solutions on a facet of the newton polytope.
checkSameMapping = method();
checkSameMapping (List, List, ZZ) := (currentPermutation, permutationList, numVertices) -> (
	local nochange;
	for i from 1 to numVertices do (
		noChange = 0;
		for j from 0 to length currentPermutation - 1 do (
			if permutationList#(currentPermutation#j)#(i - 1) == i then noChange = noChange + 1;
			-- if at any point two permutations fix the same number, return true
			if noChange >= 2 then return true;
		);
	);
	return false;
)

-- Takes a list that has the numbers 1 - n that represents a permutation by sending i to list#(i - 1)
-- returns the list that represents the inverse permutation in the same manner.
inversePermute = method();
inversePermute (List) := (permutation) -> (
	local newPermute;
	newPermute = new MutableList;
	for i from 0 to length permutation - 1 do (
		newPermute#(permutation#i - 1) = i + 1;
	);
	return new List from newPermute;
)


-- Given a current graph, and a list of permutations, returns true if three of {permutation,permutation^-1}
-- for each permutation on the given graph map a given number to the same output number. If more speed is needed,
-- can change to take a list of inverse permutations so that they don't have to be calculated repeatedly.
checkThreeMap = method();
checkThreeMap (List, List) := (currentGraph, Permutations) ->(
	local iBucket; local numVertices; local currentInverse;
	numVertices = length(Permutations#0);
	for i from 0 to numVertices - 1 do (
		iBucket = new MutableList;
		for j from 0 to numVertices - 1 do (
			iBucket#j = 0;
		);
		for j from 0 to length currentGraph - 1 do (
			currentInverse = inversePermute(Permutations#(currentGraph#j));
			iBucket#(Permutations#(currentGraph#j)#i - 1) = iBucket#(Permutations#(currentGraph#j)#i - 1) + 1;
			if not(currentInverse#i == Permutations#(currentGraph#j)#i) then (
				iBucket#(currentInverse#i - 1) = iBucket#(currentInverse#i - 1) + 1;
			); 
		);
		iBucket = new List from iBucket;
		for j from 0 to length iBucket - 1 do (
			if iBucket#j >= 3 then (
				return true;
			);
		);
	);
	return false;
)


--Computes the cycle type of a permutation given in list form.
cycleType = method()
cycleType (List) := String => P->(
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
    return sort cType;
)

-- Computes the hessian of a given polynomial with a given number of actions
Hessian = method()
Hessian (Thing, ZZ) := (func, actions) -> (
	rows = new List;
	for i from 1 to actions do (
		row = new List;
		for j from 1 to actions do (
		      row = append(row, diff(x_j, diff(x_i, func)));
		);
		row = append(row, diff(z, diff(x_i, func)));
		rows = append(rows, row);
	);
	row = new List;
	for j from 1 to actions do (
		row = append(row, diff(x_j, func));
	);
	row = append(row, diff(z, func));
	rows = append(rows, row);
	return determinant(matrix(rows));
)



actions = 2;
fundDomain = 2;

-- setting up output files
outputString = "Data/jengaInfo_" | toString(fundDomain) | "_" | toString(actions) | ".txt";
file1 = openOut(outputString);
file1 << "vertices: " << fundDomain << endl;
file1 << "actions: " << actions << endl;
file1 << "results: " << endl;
file1 << close;
file1 = openOutAppend(outputString);



-- sets up the original representation of the graph and solves for the number of solutions
DenseGraph = AdjacentDensePeriodicMatrix2(actions, fundDomain);
DF = DenseGraph_4;
varMap = map(DenseGraph_1, DenseGraph_3);
-- determines the number of edges in the graph
numEdges = numgens DenseGraph_1 - actions*2 - 2;

Generators = new Array from take(gens DenseGraph_1, actions*2 + 2);
temp = new Array from take(gens DenseGraph_1, actions*2 + 1);
Z = (ZZ/2039) Generators;
Zp = (ZZ/2039) temp;
use Z;
tempIdealGenerators = first entries generators DenseGraph_2;

inverseEqs = new List;
for i from 1 to actions do (
	inverseEqs = append(inverseEqs, x_i*y_i - 1);
)


specialization = new List;
for i from 1 to actions do (
    specialization = append(specialization, x_i);
)
specialization = append(specialization, z);
for i from 1 to actions do (
    specialization = append(specialization, y_i);
)
spec = append(specialization, 0);
idealMap = map(Zp, Z);
specialization = append(specialization, zi);
a = new List;
for i from 1 to numEdges do (
    a = append(a, random(2037) + 1); -- can change to symbolic by placing e_i where random is;
)

-- reads in the contents of a file (hopefully containing just a list)
file2 = openIn("jengaInput_" | toString(fundDomain) | "_" | toString(actions) | ".txt");
startEdgeList = value(read(file2));
for l from 0 to length startEdgeList - 1 do (
	print("Starting Graph: " | toString(startEdgeList_l) | " " | toString(l*100/(length startEdgeList)) | "% done");
	numSubgraphs = 2^(length startEdgeList_l);
	for k from 1 to numSubgraphs - 1 do (
		-- iterates through the binary representation of the subgraph and creates the corresponding edge list
		edgeList = new List;
		for i from 0 to length startEdgeList_l - 1 do (
			if getBit(k, i) == 1 then (
				edgeList = append(edgeList, startEdgeList_l_i);
			);
		);
	
		-- Add checking to ensure the graph contains an edge in each direction.
		-- Also add some checking for various similar cases
		
		-- Creates the specialization needed for the graph
		currentSpecialization = new List;
		j = 0;
		for i from 1 to length a do (
			if i == edgeList_j then (
				currentSpecialization = append(currentSpecialization, a_(i - 1));
				j = j + 1;
				if(j == length edgeList) then (
					while length currentSpecialization < length a do (
						currentSpecialization = append(currentSpecialization, 0);
					);
					break;
				);
			) else (
				currentSpecialization = append(currentSpecialization, 0);
			);
		);

		-- COMMENTED OUT UNFINISHED SYMBOLIC SOLVE
	
		-- creates the map to specialize the graph (symbolically)
		specMap = map(Z, DenseGraph_1, specialization | currentSpecialization);
		recordedVolume = polytopeVolume(specMap(DF_0), actions);
		
		-- set up the operator
		DFs = apply(DF, n -> specMap(n));
		DFsN = apply(DFs, n -> convexHull transpose matrix getExponents(n, actions));
		if(not isFullDimensional(DFsN_0)) then (
			continue;
		);

		file1 << "edgeList: " << edgeList << endl;
		file1 << "polytopeVolume: " << recordedVolume << endl;
		-- get the corresponding faces
		DFan = normalFan DFsN_0;
		
		coneListDF = coneList(DFan); -- does not return the correct number of cones. Issue in the function is that cones(i, curFan) returns a list containing an empty list.
		sysDF = apply(DFs, n -> apply(coneListDF, m -> getFaceFunc(m, n, 0)));

		-- Here need to add a check that computes the hessian and performs the same loop below but only on the original system + the hessian
		-- Then check the ideals for the correct dimension and degree in the same way.
		use Zp;
		J = idealMap(ideal inverseEqs);
		I = idealMap(specMap(ideal DF)) + ideal(Hessian(idealMap(specMap(DF_0)), actions)) + J;
		file1 << "dim = " << dim I << ", deg = " << degree I << endl;
		use Z;


		-- check for solutions at infinity	
		theIdeals = new List;
		dims = new List;
		degs = new List;
	
		for i from 0 to length coneListDF - 1 do (
			theList = new List;
			for j from 0 to actions do (
				theList = append(theList, (sysDF_j)_i); 
			);
			theIdeals = append(theIdeals, (ideal theList) + (specMap DenseGraph_2)); -- ask about this line as well
			dims = append(dims, dim theIdeals_i);
			degs = append(degs, degree theIdeals_i);
		);
		added = false;
	       	for i from 0 to length dims - 1 do (
			if dims_i != -1 then (
	    			ROWS = entries coneListDF_i;
				for j from 0 to length ROWS - 2 do (
					if sum(ROWS_j) != 0 then (
						-- print result to a file
    					file1 << "Cones with solutions: " << endl;
    					temp = new List;
    					for j from 0 to length dims - 1 do (
					    if dims_j != -1 then temp = append(temp, coneListDF_j);
    					);
    					file1 << temp << endl << endl;
					temp = new List;

			    	added = true;
			    	break;
			    );
	    		);
		        if added then break;
		    	if sum(ROWS_(length ROWS - 1)) == 0 then (
			    -- append results to a file
    			    file1 << "Cones with solutions: " << endl;
    			    temp = new List;
    			    for j from 0 to length dims - 1 do (
				if dims_j != -1 then temp = append(temp, coneListDF_j);
    			    );
    			    file1 << temp << endl << endl;
			    temp = new List;
			    
			    added = true;
			    break;
		    	);
		        if added then break;
		    );
   	 	);
		if not added then (
			file1 << "no solutions on cones" << endl << endl;
		);
	
	);
)

file1 << close;




