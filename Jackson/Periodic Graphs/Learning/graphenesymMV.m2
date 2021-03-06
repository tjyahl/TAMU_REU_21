--See graphene picture for labeling to see what is going on
--Here we look at a laplace-beltrami opertor over the graphene
--Dsicrete Laplace-beltrami operator takes potential differences scaled by edge weights
--In periodic graph, after a floquet transform, this operators are 
--determined fully by the vertices of a chosen fundamental domain, u and v in the image

--The graphene is embedded into Z^2 so we have two free abelian actions
--Our free actions have inverses so we have Laurent Polynomials
--That K[x,x^{-1}] would be one variable Laurent Polynomials over K
--We are interested in the spectrum, so we have an additional variable


clearAll


--R = frac(QQ[a,b,c])[x_1,x_2,la,y_1,y_2] 
R = QQ[x_1,x_2,la,y_1,y_2,a,b,c]
--We have 3 variables defined x_1, x_2 , and z, the y_i will be the inverse
---of x_i
--a,b,c are edge weights
--To establish this inverse behavior we need to define a quotient ring
I = ideal(x_1*y_1 - 1, x_2*y_2 - 1)


--this is technically x_1x_2 * the operator as we do not want y_1 and y_2 present for differentiation
operator = matrix{{(a+b+c-la)*x_1*x_2, (-a-b*x_1-c*x_2)*x_1*x_2},{-a*x_1*x_2-b*x_2-c*x_1, (a+b+c-la)*x_1*x_2}}

coefficientRing(R)
--build F, system of critical point equations    
F = {};    
F = append(F, sub(det(operator,Strategy => Cofactor),R));
for i from 1 to 2 do (
    F = append(F,(diff(x_i, F_0))); --adding the partial derivatives with respect to x_i to F
);

--we create an ideal J which gives our critical point equations with the condition x_1, x_2 \neq 0 by defining inverses via y_1 and y_2
J = ideal (F) + I
dim J -- dimension is 3 since we have the free edge weight parameter variables
degree J 

--we eliminate variables x_1 ... y_2 in order to be left with the values of lamdba where critical points occur
H = eliminate(J,{x_1,x_2,y_1,y_2})

factor H_0 --when this is zero there exists a critical point...this gives us the eigen values of lambda where critical points exist

S = QQ[x_1,x_2,y_1,y_2,a,b,c]

fermiMap = map(S,R,{x_1,x_2,a+b+c ,y_1,y_2,a,b,c})

fF = fermiMap ideal F
fI = fermiMap I
fJ = fF + fI
dim fJ 
degree fJ

--When we fix a value of lambda, we get the "fermi surface" at that value.
--In this case lambda is a+b+c we are solving for when both (a + b x_1 + c x_2) =0 and (a + b y_1 + c y_2) =0
--this is the intersection of a curve and a hyperbola as in the graphene623fermi png
--See the blochvariety623 png to see F_0 in S_1 \times S_1 \times \RR where S_1 is the complex unit circle


--Lets try another. For all laplace beltrami operators, when lambda = 0 there exists a solution at x_i = 1

fermiMap = map(S,R,{x_1,x_2,0,y_1,y_2,a,b,c})

fF = fermiMap ideal F
fI = fermiMap I
fJ = fF + fI
dim fJ 
degree fJ

--In this case we have solutions (a+b+c)^2 - (a + b x_1 + c x_2)*(a + b y_1 + c y_2)=0
--This can be rewritten as F0 =  (2ab+2ac+2bc) - ab(x_1 + y_1) - ac(x_2 + y_2) - bc(x_1y_2 + x_2y_1), F1 = ab(x_1-y_1) + bc(x_1y_2 - x_2 y_1), F2 = ac(x_2-y_2) + bc(-x_1y_2 + x_2 y_1)
-- Clearly x_i = 1 for each i gives a critical point

--F1 = F2 = 0 gives us ab(x_1-y_1) = -ac(x_2-y_2) = bc(-x_1y_2 + x_2 y_1), we may scale by x_1x_2 
-- We then can set ab(x_1-y_1) = ac(x_2-y_2) to get a circle, and ab(x_1-y_1) = bc(-x_1y_2 + x_2 y_1)  and ac(x_2-y_2) = bc(-x_1y_2 + x_2 y_1)  to get two hyperbolas.
--there are only solutions where all three intersect.
--At most any 2 hyperbolas can have 4 points of intersection, similarly a circle and hyperbola can only intersect 4 times at most. In this case we have x_1,x_2 = (\pm 1, \pm 1) are the four solutions.

--We finally then have that of these four solutions only one satisfies  (2ab+2ac+2bc) - ab(x_1 + y_1) - ac(x_2 + y_2) - bc(x_1y_2 + x_2y_1) and so we have only one critical point at \lambda = 0 given by x_1 = x_2 = 1
--See graphene623fermi0.png 




--In the paper Frank and I are writing, we introduce Dense periodic graphs.
--These are periodic graphs with as many edges as possible, that is if F is a fundamental domain with an edge from F to wF
--Then the subgraph given by the vertices of F and wF form a complete graph
--See pictures.

--Taking a laplace Beltrami operator over a dense periodic graph we can count the critical points exactly

--Example Here is laplace beltrami operator over a dense periodic graph with a fundamental domain of 2 vertices in Z^2
clearAll
load "functions.m2"
--calling routine from some old code

denseData = AdjacentDensePeriodicMatrix(2,2) --this gives us a tuple, the operator, the ring it is in, and the ideal setting x_iy_i=1
--2,2 is giving us the laplace beltrami operator over the periodic graph shown in Dense22periodicgraph.png

D_0 -- the matrix representing a general laplace beltrami operator over the graph
DF = {};    
DF = append(DF, sub(det(denseData_0,Strategy => Cofactor),denseData_1));
for i from 1 to 2 do (
    DF = append(DF,(diff(x_i, DF_0)));
);
DJ = ideal (DF) + denseData_2
--dim DJ expect 9, output takes too long
--degree DJ -- output takes too long
--we eliminate variables x_1 ... y_2 in order to be left with the values of lamdba where critical points occur
--H = eliminate(DJ,{x_1,x_2,y_1,y_2}) --not even going to try
--factor H_0


--Just remark that we know there are 32 critical points counted to multiplicity

--Lets take this dense graph and specialize some edges to zero


specmap = map(denseData_1, denseData_1, {x_1,x_2,y_1,y_2,z,0,e_2,0,0,0,e_6,0,0,e_9})

--lets apply the same steps we took before for graphene and look at the result
SF = specmap DJ
SI = specmap denseData_2
SJ = SF + SI
dim SJ 
degree SJ
SH = eliminate(SJ,{x_1,x_2,y_1,y_2,e_1,e_3,e_4,e_5,e_7,e_8})
factor SH_0

--This is the graphene. 
--This brings us to the idea for the project.

--Every periodic graph can be embedded in a dense periodic graph.

--As we specialize edges to 0, this is the same as removing them from our graph

--as we removed edges our dense periodic graph became the graphene, ``losing'' 20 critical points.
--We are interested in detailing this process and understanding what, why, and when solutions vanish.



--Why do we care? Spectral edge conjecture. 
-- Look at the bloch variety 623 png. 
--Notice the bloch vairety of graphene is composed of two sheets,
--The highlighted region on the lambda axis is the projection of the bloch variety onto it
--Each of these red segments is called a spectral band. The end points are called spectral edges. 

--The conjecture states that: Extrema are isolated, Extrema are non-degenerate, each extremal value occurs in a single band.

--To be precise we are referring to edges. In the example the two sheets of the bloch variety did not over lap, so there 
--were 4 edges (end points of the projection of the sheets on to the lambda axis), two for each sheet.
--The extrema in the conjecture refer to the edges, not extrema of a single sheet.
--For example we have only 2 edges in blochvariety111.png and a single band. Thus although there are degenerate critical points.
--these are part of the same spectral band.

--Since extrema are critical points, understanding the structure of the critical points is a first step
--towards this conjecture.

--For discrete periodic operators this does not hold in general in higher dimension (>2D), but 
--it is of interest for what operators and what graphs it does and does not hold.


--Mixedvolume and counting the number of solutions
clearAll
load "../functions.m2"
denseData = AdjacentDensePeriodicMatrix(3,2); --returns matrix, ring with arbitrary edges, ideal identifying variables, and ring in same vars without edge variables, and system in ring with edges

DF = denseData_4;    

--we specialize edges to get our poly in 3 variables

a = {0,21,43,34,54,35,96,71,83,9,27,62,59,13} --pick sufficiently "random" values

--try the following:

-- a = {0,0,43,34,0,0,96,71,0,9} --subgraph containing graphene
-- a = {0,0,43,0,0,0,96,0,0,9} -- this is graphene
-- a = {0,0,43,0,0,0,96,71,0,9} --a graph inbetween these two

-- a = {0,0,43,34,0,0,96,71,0,9} --another graph, not containing graphene but has 32 solutions might be worth checking subgraphs


specmap = map(denseData_3, denseData_1, {x_1,x_2,x_3,y_1,y_2,y_3,z,0,a_2,a_3,0,0,a_6,a_7,0,0,a_10,a_11,0,a_13})

toprunelist = {}; --This list will be taken as a parameter in order to get a m+1 dimensional polytope
m=3 --this is the number of variables not including z
for i from 1 to m do (
    toprunelist = append(toprunelist, 2*m-i);
    );

DFs = apply(DF, n -> specmap (n)); --puts system in ring without edge variables

DFsP = apply(DFs, n -> pruneList ( exponents(n ) , toprunelist)); --prunes/purges out the dimensions containing the inverse variables

DFsN = apply(DFsP, n -> convexHull (transpose matrix n)); --list of polytopes

mixedV = myMixedVolume(DFsN) --returns the mixed volume, if out system is sufficiently generic, Bernstein's theorem tells us the number of solutions = mixed volume, in general the number of isolated solution in the algebraic torus is <= mixed volume

volume DFsN_0 -- Notice MV is 3! * the volume of the polytope of the characteristic polynomial, for our systems this is always true (where 3 is n instead where n is the dimension of the polytope of DF_0).

--Our polynomials are not generic, but we have a criteria to check when the mixed volume is the number of solutions, this is when the "facial systems have no solutions"
--probably do not have time to discuss these very clearly today, so for the sake of enabling you to play with this I will describe what to look for

--We will generate cones of the polytope, and these cones will identify difference facial systems

DFan = normalFan DFsN_0;
conelistDF = coneList(DFan);
sysDF = apply(DFs, n -> apply(conelistDF, m-> getFaceFunc(m,n,1))); --the 1 is just to let the function know how many variables are not inverse variables, see function.m2 description for more info
--this gives us a collection of all facial systems


--we extract some information about these systems:
     theideals = {};
     dims = {};
     degs = {};
     for i from 0 to #conelistDF - 1 do (
	 thelist = {};
	 for j from 0 to m do (
	     thelist = append(thelist, (sysDF_j)_i);
	     );
	 theideals = append(theideals,(ideal thelist)+(specmap denseData_2));
	 dims = append(dims, dim theideals_i);
	 degs = append(degs, degree theideals_i);
	 );
degs
dims

--For degrees these should be all zero except 1 at the 5th entry and a nonzero value at the 14th entry
--for dims it should be -1 where there are 0s, 1 where there is a 1, and 2 at the 14th entry

--If this is the case then we have no solutions at any relevant facial system and so we know the number of isolated solutions is given by the mixed volume

--A glimpse as to why we don't need to worry about the 5th or 14th entry is the following

conelistDF_4
conelistDF_13 -- take the sum of the the columns to get the vector in question

--These two vectors correspond to the base and peak of the polytope (when the height is taken to correspond to the z or "lambda" axis)
--We can later go over why solutions at the peak are not solutions at all (because lambda cannot be zero here)
--Our criteria also actually accounts for these solutions on the base, these are part of the mixedvolume number of solutions
