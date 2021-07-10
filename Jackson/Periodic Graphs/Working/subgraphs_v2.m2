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



actions = 2;
fundDomain = 3;


-- sets up the original representation of the graph and solves for the number of solutions
DenseGraph = AdjacentDensePeriodicMatrix(actions, fundDomain);


subtractor = id_(DenseGraph_1^fundDomain);
for i from 1 to actions do (
    subtractor = subtractor*x_i;
)
subtractor = subtractor*z;


originalDet = det(DenseGraph_0 - subtractor, Strategy => Cofactor);


-- the maximum number of solutions
maxSolutions = polytopeVolume(originalDet, actions);




-- determines the number of edges in the graph
numEdges = numgens DenseGraph_1 - actions*2 - 1;

graphsChecked = 0;
totalGraphs = 2^numEdges;
for i from 0 to actions*fundDomain - 1 do (
    totalGraphs = totalGraphs - binomial(numEdges, i);
    graphsChecked = graphsChecked + binomial(numEdges, i);
)

-- Used to keep track of program progress
progress = method();
progress (ZZ) := (fill) -> (
    local bPrint; local pTime;
    bPrint = false;
    pTime = cpuTime();
    while graphsChecked < totalGraphs do (
	if bPrint then ( -- for some reason sleep is ignored when the function is scheduled.
	    print "----------------------------------";
	    print ("Time Spent: " | toString(cpuTime()) | "s");
	    print ("Graphs Checked: " | toString(graphsChecked));
	    print ("Graphs Remaining: " | toString(totalGraphs - graphsChecked));
	    print "----------------------------------";
	    bPrint = false;
	) else (
	    if cpuTime() - pTime > 10 then (
		bPrint = true;
		pTime = cpuTime();
	    );
	);
    );
    return fill;
)




originalMatrix = DenseGraph_0

specialized = originalMatrix;
for i from 1 to numEdges do (
    if i != 5 and i != 8 and i != 9 and i != 12 then (
	specialized = sub(specialized, e_i=>0);
    );
)

file = openOut("matrices.txt")
file << specialized << endl << endl;
polytopeVolume(det(specialized), actions)
-- Set to contain all the minimal subgraphs
minimalGraphs = new List;

schedule createTask (progress, 1);
-- An implementation of Gosper's Hack in order to iterate from smallest subgraph to largest subgraph
for i from actions*fundDomain to numEdges do (
    gSet = (1 << i) - 1;
    limit = (1 << numEdges);
    
    -- loop through all subsets containing i elements
    while gSet < limit do (
	
	-- loops through the bit representation to get elements from the subset
	curDet = originalDet;
	for j from 0 to numEdges - 1 do (
	    
	    -- if the bit is 0 specialize the corresponding edge to 0
	    if getBit(gSet, j) != 1 then (
		curDet = sub(curDet, e_(j + 1)=>0);
	    )
	);
    
    	-- represents the graph with a number whose bit representation tells which edges are included
	-- and a number that represents the polytope volume of the determinant of the laplace beltrami
	-- operator over the graph.
	subgraph = false;
	for i from 0 to length minimalGraphs - 1 do (
	    if isSubgraph(minimalGraphs_i_0, gSet) then (
		subgraph = true;
	    );
	);
    	if not subgraph then (
	    if polyTerms(curDet, actions, fundDomain) then (
		curGraph = {gSet, polytopeVolume(curDet, actions)};
		if curGraph_1 == maxSolutions then minimalGraphs = append(minimalGraphs, curGraph);
	    );
	);
    
    	
	graphsChecked = graphsChecked + 1;
	-- used to identify the next subset to be listed
	c = gSet & -gSet;
	r = gSet + c;
	gSet = ((xor(r, gSet) >> 2) // c) | r;
    );
    
);

-- prints the binary representation of the graphs in minimalGraphs (so that the graphs can be drawn)
-- will likely be replaced with something that writes the binary to a file.
file = openOut("results.txt");
file << "Time to completion: ";
file << cpuTime();
file << "s" << endl;
file << "vertices: " << fundDomain << endl;
file << "actions: " << actions << endl;
file << "results:" << endl;
minimalGraphs = {{43092, 108}, {43106, 108}, {43148, 108}, {43169, 108}, {43274, 108}, {43281, 108}, {50260, 108}, {50274, 108}, {50316, 108}, {50337, 108}, {50442, 108}, {50449, 108}, {71764, 108}, {71778, 108}, {71820, 108}, {71841, 108}, {71946, 108}, {71953, 108}, {82516, 108}, {82530, 108}, {82572, 108}, {82593, 108}, {82698, 108}, {82705, 108}, {136276, 108}, {136290, 108}, {136332, 108}, {136353, 108}, {136458, 108}, {136465, 108}, {139860, 108}, {139874, 108}, {139916, 108}, {139937, 108}, {140042, 108}, {140049, 108}};
for i from 0 to length minimalGraphs - 1 do (
    view = new List;
    for j from 0 to 20 do (
	item = 0;
	if getBit(minimalGraphs_i_0, j) == 1 then item = e_(j+1);
    	view = append(view, item);
    );
    file << view << endl;
);
close file;



-- The binary representation is as follows: the nth bit (starting from the 0th) is 1 if and only
-- if the edge e_(n+1) is included in the graph. After printing in the manner above, the nth entry in 
-- the view array is 1 if and only if the graph contains the edge e_(n+1). The edges are assigned by
-- numbering the vertices in the fundamental domain v_1 through v_n. Then the edge from v_i to x_k v_j
-- is given by e_((k-1)*n^2 + n*(i-1) + j). The last C(n,2) edges make up the complete graph of the 
-- fundamental domain.













