restart
load("../functions.m2")

-- useful function definitions:


-- performs the steps to remove the edge variables and the inverse variables from the given
-- polynomial. Need to be concious of any change in rings (should move to one without edges)
getExponents = method();
getExponents (Thing, ZZ) := (f, actions) -> (
    local E; local L; local K; local l;
    E = exponents f;
    L = new List;
    for i from 0 to length E - 1 do (
	l = new List;
	for j from 0 to actions - 1 do (
	    l = append(l, E_i_j - E_i_(j + actions));
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




actions = 4;
fundDomain = 3;

-- setting up output files
outputString = "Data/coneSolutions_" | toString(fundDomain) | "_" | toString(actions) | ".txt";
file1 = openOut(outputString);
file1 << "vertices: " << fundDomain << endl;
file1 << "actions: " << actions << endl;
file1 << "results: " << endl;
file1 << close;
file1 = openOutAppend(outputString);

outputString = "Data/noConeSolutions_" | toString(fundDomain) | "_" | toString(actions) | ".txt";
file2 = openOut(outputString);
file2 << "vertices: " << fundDomain << endl;
file2 << "actions: " << actions << endl;
file2 << "results: " << endl;
file2 << close;
file2 = openOutAppend(outputString);



-- sets up the original representation of the graph and solves for the number of solutions
DenseGraph = AdjacentDensePeriodicMatrix2(actions, fundDomain);
DF = DenseGraph_4;

Generators = new Array from take(gens DenseGraph_1, actions*2 + 2);
Z = (ZZ/2039) Generators;
tempIdealGenerators = first entries generators DenseGraph_2;

-- determines the number of edges in the graph
numEdges = numgens DenseGraph_1 - actions*2 - 2;

-- Used to create an ideal of Z that implies that the x_i's and z must be nonzero
converter = new List;
for i from 1 to actions do (
    converter = append(converter, x_i);
);
for i from 1 to actions do (
    converter = append(converter, y_i);
);
converter = append(converter, z);
for i from 1 to numEdges do (
    converter = append(converter, 0);
);
convert = map(Z, DenseGraph_1, converter);

idealGenerators = new List;
for i from 0 to length tempIdealGenerators - 1 do (
    idealGenerators = append(idealGenerators, convert(tempIdealGenerators_i));
);
idealGenerators = append(idealGenerators, z*zi - 1);

I = ideal(idealGenerators);

-- generates all the permutations of fundDomain elements and counts them actions number of times
permutes = new List;
c = new MutableList;
for i from 1 to fundDomain do c = append(c, 0);
counter = 1;
generate = new List from 1 .. fundDomain;
permutes = append(permutes, generate);


while counter < fundDomain do (
    if c#counter < counter then (
	if counter % 2 == 0 then (
	    generate = switch(0, counter, generate);
	) else (
	    generate = switch(c#counter, counter, generate);
	);
    	
	permutes = append(permutes, generate);
	c#counter = c#counter + 1;
	counter = 1; 
    ) else (
    	c#counter = 0;
	counter = counter + 1;
    );
);

currentPermutation = new MutableList;
for i from 1 to actions - 1 do currentPermutation = append(currentPermutation, 0);
currentPermutation = append(currentPermutation, 0); -- used to skip the case with all identity (change 0 to 1)


specialization = new List;
for i from 1 to actions do (
    specialization = append(specialization, x_i);
)
specialization = append(specialization, z);
for i from 1 to actions do (
    specialization = append(specialization, y_i);
)
specialization = append(specialization, zi);
a = new List;
for i from 1 to numEdges do (
    a = append(a, random(2037) + 1);
)

-- loop through all combinations of permutations
while true do (
    try currentPermutation == false then break else (
	edgeList = generateEdges(new List from currentPermutation, permutes);
	j = 0;
	currentSpecialization = new List;
	for i from 1 to length a do (
	    if i == edgeList_j then (
		currentSpecialization = append(currentSpecialization, a_(i-1));
		j = j + 1;
	    ) else (
	    	currentSpecialization = append(currentSpecialization, 0);
	    );
       	);
	
	specMap = map(Z, DenseGraph_1, specialization | currentSpecialization);
	-- set up the operator
	DFs = apply(DF, n -> specMap(n));
	DFsN = apply(DFs, n -> convexHull transpose matrix getExponents(n, actions + 1));
	
	
	-- get the corresponding faces
	DFan = normalFan DFsN_0;
	coneListDF = coneList(DFan);
	sysDF = apply(DFs, n -> apply(coneListDF, m -> getFaceFunc(m, n, 0))); -- error: Vectors not of equal length nerd
	
	
	
	-- loop through the faces and check dimension/degree. Record the graph if dim != -1 and deg != 0
	theIdeals = new List;
	dims = new List;
	--degs = new List;
	for i from 0 to length coneListDF - 1 do (
	    theList = new List;
	    for j from 0 to actions do (
		theList = append(theList, (sysDF_j)_i); 
	    );
	    theIdeals = append(theIdeals, (ideal theList) + (specMap DenseGraph_2)); -- ask about this line as well
	    dims = append(dims, dim theIdeals_i);
	    --degs = append(degs, degree theIdeals_i);
	);	
    	-- checks that dimensions are all -1 for any face that is not the apex or base
	added = false;
        for i from 0 to length dims - 1 do (
	    if dims_i != -1 then (
	    	ROWS = entries coneListDF_i;
	    	for j from 0 to length ROWS - 2 do (
		    if sum(ROWS_j) != 0 then (
			-- print result to a file
    			file1 << "Graph: " << edgeList << endl;
    			file1 << "Cones with solutions: " << endl;
    			temp = new List;
    			for j from 0 to length dims - 1 do (
			    if dims_j != -1 then temp = append(temp, coneListDF_j);
    			);
    			file1 << temp << endl;
			temp = new List;

		    	added = true;
		    	break;
		    );
	    	);
	        if added then break;
	    	if sum(ROWS_(length ROWS - 1)) == 0 then (
		    -- append results to a file
		    file1 << "Graph: " << edgeList << endl;
    		    file1 << "Cones with solutions: " << endl;
    		    temp = new List;
    		    for j from 0 to length dims - 1 do (
			if dims_j != -1 then temp = append(temp, coneListDF_j);
    		    );
    		    file1 << temp << endl;
		    temp = new List;
		    
		    added = true;
		    break;
	    	);
	        if added then break;
	    );
    	);
	if not added then (
		-- place in noConeSolutions.txt
		file2 << edgeList << endl;
	);
	currentPermutation = nextPermutation(currentPermutation, actions, length permutes);
    );
);

file1 << close;
file2 << close;
