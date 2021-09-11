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






actions = 8;
fundDomain = 4;



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


-- loop through all combinations of permutations
i = 0;
while true do (
    try currentPermutation == false then break else (
	if checkSameMapping(new List from currentPermutation, permutes, fundDomain) or checkThreeMap(new List from currentPermutation, permutes) then (
		-- maybe add to an output file here or after edges are generated.
		currentPermutation = nextPermutation(currentPermutation, actions, length permutes);
		continue;
	);
	i = i + 1;
	currentPermutation = nextPermutation(currentPermutation, actions, length permutes);
    );
);

print i;

