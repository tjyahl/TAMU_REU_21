--See graphene picture for labeling to see what is going on
--Here we look at a laplace-beltrami opertor over the graphene
--Laplace-beltrami operator assigns edges weights
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
dim J
degree J 

--we eliminate variables x_1 ... y_2 in order to be left with the values of lamdba where critical points occur
H = eliminate(J,{x_1,x_2,y_1,y_2})

factor H_0

S = QQ[x_1,x_2,y_1,y_2,a,b,c]

fermiMap = map(S,R,{x_1,x_2,a+b+c ,y_1,y_2,a,b,c})

fF = fermiMap ideal F
fI = fermiMap I
fJ = fF + fI
dim fJ 
degree fJ

--When we fix a value of lambda, we get the "fermi surface" at that value.
--In the case lambda is a+b+c we are solving for when (a + b x_1 + c x_2) =0 and   (a + b y_1 + c y_2) =0
--this is the intersection of a curve and a hyperbola as in the png
--See the blochvariety623 png to see F_0 in S_1 \times S_1 \times \RR where S_1 is the complex unit circle


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
