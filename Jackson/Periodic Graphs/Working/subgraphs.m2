restart
load("../functions.m2")

-- useful function definitions:
polytopeVolume = method()

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







actions = 2
fundDomain = 2

-- sets up the original representation of the graph and solves for the number of solutions
DenseGraph = AdjacentDensePeriodicMatrix(actions, fundDomain)

subtractor = id_(DenseGraph_1^fundDomain)
for i from 1 to actions do (
    subtractor = subtractor*x_i
)
subtractor = subtractor*z

originalDet = det(DenseGraph_0 - subtractor, Strategy => Cofactor)
saturatedSolutions = polytopeVolume(originalDet, actions)


-- determines the number of edges in the graph
numEdges = numgens DenseGraph_1 - actions*2 - 1

-- give a function that when a polynomial is given (along with its ring), it determines the volume
-- of the newton polytope defined by the supports of that polynomial


specialized = originalDet
for i from 1 to 9 do (
    if i != 2 and i != 6 and i != 9 then (
    specialized = sub(specialized, e_i=>0)
    )
)
polytopeVolume(specialized, actions)
