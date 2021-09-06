restart
load("../functions.m2")
allowableThreads = 1;
--load("Graphs.m2") -- will need to remove if "faces", "facets", "isPure", "fVector", "skeleton",
-- "vertices", or "directProduct" as they appear in both Graphs.m2 and functions.m2

-- useful function definitions:


-- Takes the determinant and the number of actions on a graph and 
-- returns the mixed volume of the newton polytope defined by the 
-- determinant.
polytopeVolume = method();
polytopeVolume (Thing, ZZ) := (f, actions) -> (
    local E; local L; local K; local l;
    E = exponents f;
    L = new List;
    for i from 0 to length E - 1 do (
	l = new List;
	for j from 0 to actions - 1 do (
	    l = append(l, E_i_j - E_i_(j + actions));
	);
    	l = append(l, E_i_(actions*2));
	L = append(L, l);
    );
    L = unique L;
    
    K = transpose matrix L;
    
    return (actions + 1)! * volume convexHull K
    
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
	for j from 0 to actions - 1 do (
	    l = append(l, E_i_j - E_i_(j + actions));
	);
    	l = append(l, E_i_(actions*2));
	L = append(L, l);
    );
    L = unique L;
    return L
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

-- Takes a polynomial as well as the number of actions (n) and vertices (m) and returns true iff the 
-- polynomial contains the terms x_i^2m * PI_(i!=j) x_j^m as well as PI_(i!=j) x_j^m for each
-- 1 <= i <= n
polyTerms = method();
polyTerms (Thing, ZZ, ZZ) := (Poly, n, m) -> (
    local E; local l; local L; local temp; local found;
    E = exponents Poly;
    L = new List;
    for i from 0 to length E - 1 do (
	l = new List;
	for j from 0 to n - 1 do (
	    l = append(l, E_i_j - E_i_(j + n));
	);
    	l = append(l, E_i_(n*2));
	L = append(L, l);
    );
    L = unique L;
    
    for i from 0 to n - 1 do (
	temp = new List;
	for j from 0 to n - 1 do (
	    if j == i then temp = append(temp, 2*m) else temp = append(temp, m); 
	);
    	temp = append(temp, 0);
	found = false;
    	for j from 0 to length L - 1 do (
	    if L_j == temp then (
	    	found = true;
		break;
	    );
	);
    	if found then continue else return false;
    );
    return true;
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

-- Takes the original determinant and a list of edge numbers in the desired graph and specializes
-- the determinant. Can be improved as edgeList is a sorted list.
specializeDet = method();
specializeDet (Thing, List, List) := (f, edges, edgeWeights) -> (
    local j;
    j = 0;
    for i from 1 to length edgeWeights - 1 do (
	if j <= length edges - 1 then (
	    if i != edges_j then f = sub(f, e_i=>0) else (
		f = sub(f, e_i=>edgeWeights_i);
		j = j + 1;
	    );
	) else (
	    f = sub(f, e_i=>0);
	);
    );
    return f
)




actions = 2;
fundDomain = 3;


-- sets up the original representation of the graph and solves for the number of solutions
DenseGraph = AdjacentDensePeriodicMatrix(actions, fundDomain);
DF = DenseGraph_4;



-- determines the number of edges in the graph
numEdges = numgens DenseGraph_1 - actions*2 - 1;




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
currentPermutation = append(currentPermutation, 0); -- used to skip the case with all identity
Info = new List;

specialization = new List;
for i from 1 to actions do (
    specialization = append(specialization, x_i);
)
for i from 1 to actions do (
    specialization = append(specialization, y_i);
)
specialization = append(specialization, z);
a = new List;
for i from 1 to numEdges do (
    a = append(a, random(500) + 1);
)

faceSolutions = new List;
noFaceSolutions = new List;



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
	
	specMap = map(DenseGraph_3, DenseGraph_1, specialization | currentSpecialization);
	-- set up the operator
	DFs = apply(DF, n -> specMap(n));
	DFsN = apply(DFs, n -> convexHull transpose matrix getExponents(n, actions));
	
	
	-- get the corresponding faces
	DFan = normalFan DFsN_0;
	coneListDF = coneList(DFan);
	sysDF = apply(DFs, n -> apply(coneListDF, m -> getFaceFunc(m, n, 1)));
	
	
	
	-- loop through the faces and check dimension/degree. Record the graph if dim != -1 and deg != 0
	theIdeals = new List;
	dims = new List;
	degs = new List;
	for i from 0 to length coneListDF - 1 do (
	    theList = new List;
	    for j from 0 to actions do (
		theList = append(theList, (sysDF_j)_i); -- program fails at this line
	    );
	    theIdeals = append(theIdeals, (ideal theList) + (specMap DenseGraph_2));
	    dims = append(dims, dim theIdeals_i);
	    degs = append(degs, degree theIdeals_i);
	);
	Info = append(Info, {coneListDF, dims, degs, edgeList});
    	-- checks that dimensions are all -1 for any face that is not the apex or base
        for i from 0 to length dims - 1 do (
	    if dims_i != -1 then (
	    	ROWS = entries coneListDF_i;
	    	added = false;
	    	for j from 0 to length ROWS - 2 do (
		    if sum(ROWS_j) != 0 then (
		    	faceSolutions = append(faceSolutions, {edgeList, coneListDF, dims});
		    	added = true;
		    	break;
		    );
	    	);
	        if added then break;
	    	if sum(ROWS_(length ROWS - 1)) == 0 then (
		    faceSolutions = append(faceSolutions, edgeList);
		    added = true;
		    break;
	    	);
	        if added then break;
	    	if not added then (
		    noFaceSolutions = append(noFaceSolutions, edgeList);
		    break;
	    	);
	    );
    	);
	
	currentPermutation = nextPermutation(currentPermutation, actions, length permutes);
    );
);

file = openOut("coneSolutions.txt");
file << "vertices: " << vertices << endl;
file << "actions: " << actions << endl;
file << "results: " << endl;


conesWSolutions = new List;
for i from 0 to length faceSolutions - 1 do (
    file << "Graph: " << faceSolutions_i_0 << endl;
    file << "Cones with solutions: " << endl;
    temp = new List;
    for j from 0 to length faceSolutions_i_2 - 1 do (
	if faceSolutions_i_2_j != -1 then temp = append(temp, faceSolutions_i_1_j);
    );
    conesWSolutions = append(conesWSolutions, temp);
    file << temp << endl;
)

file << close;



-- The binary representation is as follows: the nth bit (starting from the 0th) is 1 if and only
-- if the edge e_(n+1) is included in the graph. After printing in the manner above, the nth entry in 
-- the view array is 1 if and only if the graph contains the edge e_(n+1). The edges are assigned by
-- numbering the vertices in the fundamental domain v_1 through v_n. Then the edge from v_i to x_k v_j
-- is given by e_((k-1)*n^2 + n*(i-1) + j). The last C(n,2) edges make up the complete graph of the 
-- fundamental domain.








