--example file for creating lacunary systems

restart
loadPackage("DecomposableSparseSystems")
loadPackage("Monodromy",FileName=>"Monodromy.m2")

----------------------
----Things to note----
----------------------
----The functions "isDecomposable", "isLacunary", and "isTriangular" are from
----the package "DecomposableSparseSystems". 
----
----The code below is computing the Galois group over the COMPLEX NUMBERS. That
----is to say, the statements below aren't all necessarily true when considering
----the Galois group over the rational numbers (using the Frobenius method).


-------------------------------------------------
----Systems with Full Symmetric Galois Groups----
-------------------------------------------------
--We can easily determine when a system has a Galois group equal to the full symmetric group.
----If a set of supports 'A' are not decomposable, the corresponding Galois group is the
----symmetric group on the number of solutions.

--A random system with this support has 3 solutions. This can be checked as follows.
----'F' is the randomly generated system and 'sols' is the set of solutions.
A = {matrix{{0,0,1,1},{0,1,0,1}},matrix{{0,0,1,1,2},{0,1,0,1,0}}}
(F,sols) = solveDecomposableSystem(A,,Software=>M2) -- seems to require PHCPACK
#sols

--This set of supports is not decomposable. This implies the Galois group is 
----the symmetric group on 3 elements, S_3. We can check this with monodromy.
----
----(Note: S_3 has 3! = 6 elements)
isDecomposable A
M = sparseMonodromy(A,Solver=>M2)
monodromyLoop(M,30)
printCycleTally M


------------------------
----Lacunary Systems----
------------------------
--As above, we can test whether a system is lacunary with 'isLacunary'.
----If a system is lacunary, its Galois group is imprimitive (subgroup of wreath product)
A = {matrix{{0,0,2,2},{0,2,0,2}},matrix{{0,0,2,2},{0,2,0,2}}}
isLacunary A


--We want a way to generate lacunary systems where we know what the Galois group SHOULD be. 
----One way to do this is as follows:
--1) Take a set of supports 'A' which are NOT DECOMPOSABLE. Compute the number of solutions, 'd'.
------Make sure d>1. The Galois group for these supports is the full symmetric group S_d.
A = {matrix{{0,0,1,1},{0,1,0,1}},matrix{{0,0,1,1},{0,1,0,1}}}
(F,sols) = solveDecomposableSystem(A,,Software=>M2)
d = #sols

--2) Find a matrix 'T' with |det(T)|>1. Multiply each support matrix in 'A' on the left by 'T'.
------(Note: You can let 'T' be a diagonal matrix and still make all lacunary systems.)
------
------(Note: the new polynomials will be in the variables corresponding to the columns of 'T'.
------in this example, the new polynomials are really polynomials in x^3 and x*y.)
T = matrix{{3,1},{0,1}}
det T
newA = apply(A,S->T*S)

--3) The new set of supports is lacunary and systems with this support have d*|det(T)| solutions. 
------The Galois group has at most |det(T)|^d*d! elements, since this is the size of the wreath product.
------
------(Note: In this case, the Galois group has at most 3^2*2! = 18 elements)
isLacunary newA

(F',sols') = solveDecomposableSystem(newA,)
#sols'

M = sparseMonodromy(newA,Solver=>M2)
monodromyLoop(M,100)
printCycleTally M
#(M#group)



--We're interested in cases where there are LESS elements of the Galois group than the wreath product. 
----That is, we want |Galois group|<|det(T)|^d*d!. An example of this is below (the Galois group has 192
----elements instead of the maximum 384 the wreath product has).
--1) This first part tells us systems with support 'A' have 4 solutions and the Galois group is S_4.
A = {matrix{{-1,0,0,0,1},{0,-1,0,1,0}},matrix{{-1,0,0,0,1},{0,-1,0,1,0}}}
isDecomposable A
(F,sols) = solveDecomposableSystem(A,)
d = #sols

--2) Chosen matrix is 'T' below.
T = matrix{{1,-1},{1,1}}
det T

newA = apply(A,S->T*S)

--3) New supports 'newA' are lacunary. The Galois group has at most 'maxGaloisSize' = 384 elements.
isLacunary newA
maxGaloisSize = abs(det T)^d*d!

--The Galois group actually has 192 elements. (Taking 500 loops may not be enough to get them all.)
M = sparseMonodromy(newA,Solver=>M2)
monodromyLoop(M,500)
printCycleTally M
#(M#group)



--------------------------
----Triangular Systems----
--------------------------
--As above, we can test whether a system is triangular with 'isTriangular'.
----If a system is triangular, its Galois group is imprimitive (subgroup of wreath product)
A = {matrix{{0,1,2},{0,2,4}},matrix{{0,0,1,1},{0,1,0,1}}}
isTriangular A

--We want a way to generate triangular systems where we know what the Galois group SHOULD be. 
----One way to do this is as follows (more complicated than lacunary):
--1) Choose a set of supports in a specific dimension 'k' (here, we'll use k=2, so 2 polynomials
------in 2 variables. So we need 2 supports 'A'). Make sure the supports are NOT DECOMPOSABLE
------and determine the number of solutions 'd'. Make sure d>1. This will be out "subsystem".
k = 2
A = {matrix{{0,0,1,1},{0,1,0,1}},matrix{{0,0,1},{0,1,0}}}
(F,sols) = solveDecomposableSystem(A,)
d = #sols

--3) Put these supports "up a dimension" by adding rows of zeros to the supports. For each row
------of zeros you add, include another support. Here, we added 1 rows of zeros, so we add 1
------additional supports. A system with this support has a subsystem in the first 'k' equations.
newA = {matrix{{0,0,1,1},{0,1,0,1},{0,0,0,0}},matrix{{0,0,1},{0,1,0},{0,0,0}},matrix{{0,0,0,1},{0,1,1,1},{0,0,1,2}}}

--4) Make sure the supports you added give your blocks >1 solution and don't have any additional
------structure. To do this, cut off the first 'k' rows and check that the supports you get
------have r>1 solutions per system and aren't decomposable.
------
------(Note 1: The "residual supports" you get by cutting off the first 'k' rows correspond to 
------the support you get by plugging in a solution of the subsytem into the remaining equations.)
------
------(Note 2: You don't need to understand all the details of the code here)
addedSupports =  drop(newA,k) --drop the first 'k' supports to get the ones added on.
residualSupports = apply(addedSupports,S->S^(toList(k..numRows S - 1))) --cut first 'k' rows (remember indexing starts at 0)

isDecomposable residualSupports
(F,sols) = solveDecomposableSystem(residualSupports,)
r = #sols

--5) Congrats! The supports 'newA' are triangular. There are r*n solutions to a generic system with this support
------and the Galois group has at most (r!)^d*d! elements, since this is the size of the wreath product.
------
------(Note: the only difference in this formula is the r! instead of |det(T)|)
maxGaloisSize = (r!)^d*d!
M = sparseMonodromy(newA,Solver=>M2)
monodromyLoop(M,100)
printCycleTally M
#(M#group)


--Again, we're interested in supports who's Galois group is SMALLER than the wreath product. Finding
----triangular supports with this property is one of the main goals. (I don't have any good examples on-hand)


