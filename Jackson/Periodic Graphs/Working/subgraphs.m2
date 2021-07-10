restart
load("../functions.m2")
allowableThreads = 8;
--maxAllowableThreads = 8;
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

-- Takes the bit representation of the edges in a graph, the determinant of the laplace beltrami
-- operator over the graph, and the number of actions on the graph and gives an array containing
-- the bit representation and the polytope volume
calculate = method();
calculate (ZZ, Thing, ZZ) := (bitRep, curDet, actions) -> () -> (
    return {bitRep, polytopeVolume(curDet, actions)};
)

-- Set to contain all minimal subgraphs
minimalGraphs = new List;

-- Takes a task that represents an instance of the calculate method above and checks to see whether the
-- graph should be added to the list of minimal subgraphs with the largest polytope volume.
membership = method();
membership (Task, Task, QQ) := (prevPush, graphCal, maxSolutions) -> () -> (
    local curGraph; local add; local misc;
    add = true;

    while not isReady prevPush do nanosleep 10;
    
    while not isReady graphCal do nanosleep 10;
    
    curGraph = taskResult graphCal;
    
    try curGraph_1 then (
    	if curGraph_1 == maxSolutions then (
	    for i from 0 to length minimalGraphs - 1 do (
	    	if isSubgraph(minimalGraphs_i_0, curGraph_0) then (
		    add =  false;
		    break;
	    	);
	    );
    	
    	if add then (
	    minimalGraphs = append(minimalGraphs, curGraph);
	);
    );
    return true;
)


actions = 2;
fundDomain = 2;


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

tm = x-> () -> x;
--tm2 = x->()->x;
--tm3 = x->()->x;
prevPush = schedule tm(1);
--test1 = createTask tm2(2);
--test2 = createTask tm3(3);
--addDependencyTask(prevPush, {test1, test2})

while not isReady(prevPush) do sleep 1;


-- An implementation of Gosper's Hack in order to iterate from smallest subgraph to largest subgraph
for i from 1 to numEdges do (
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
	    );
	);
    
    	-- handles the parralelization of multiple calculations
	calTask = schedule calculate(gSet, curDet, actions);
	curPush = schedule membership(prevPush, calTask, maxSolutions);
	
	prevPush = curPush;
	
	
	
	-- used to identify the next subset to be listed
	c = gSet & -gSet;
	r = gSet + c;
	gSet = ((xor(r, gSet) >> 2) // c) | r;
    );
    
);

-- wait for final push to array
while not isReady prevPush do sleep 1;


-- prints the binary representation of the graphs in minimalGraphs (so that the graphs can be drawn)
-- will likely be replaced with something that writes the binary to a file.
file = openOut("results.txt");
file << "Time to completion: ";
file << cpuTime();
file << "s" << endl;
file << "vertices: " << fundDomain << endl;
file << "actions: " << actions << endl;
file << "results:" << endl;
for i from 0 to length minimalGraphs - 1 do (
    view = new List;
    for j from 0 to numEdges - 1 do (
    	view = append(view, getBit(minimalGraphs_i_0, j));
    );
    file << view << endl;
);


-- The binary representation is as follows: the nth bit (starting from the 0th) is 1 if and only
-- if the edge e_(n+1) is included in the graph. After printing in the manner above, the nth entry in 
-- the view array is 1 if and only if the graph contains the edge e_(n+1). The edges are assigned by
-- numbering the vertices in the fundamental domain v_1 through v_n. Then the edge from v_i to x_k v_j
-- is given by e_((k-1)*n^2 + n*(i-1) + j). The last C(n,2) edges make up the complete graph of the 
-- fundamental domain.

























