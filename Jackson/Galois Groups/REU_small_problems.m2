--An ongoing list of small problems to examine for understanding.

--------------------------------------------------------------------------
--Problem 1
----Show that the wreath product (S_2 wr S_2) has 8 elements and determine
----the different cycle types of elements.
----
----(This is one to do by hand)



--------------------------------------------------------------------------
--Problem 2
----The matrices 'A' below represent supports for a polynomial system.
----The columns of the first matrix are the exponents for monomials in
----the first polynomial and the columns of the second matrix are the 
----exponents for monomials in the second polynomial.
----
----Draw the supports to see what they look like. Create a random system 
----with these supports and determine the number of solutions using 'solveSystem'.
loadPackage("NumericalAlgebraicGeometry",Reload=>true)
A = {matrix{{0,1,2,0},{0,0,0,2}},matrix{{0,1,2,0},{0,0,0,2}}}
random(ZZ)
R = CC[x,y]
F = {4+3*x+3*x^2+4*y^2, 1+8*x+5*x^2+6*y^2}--type your system here
solveSystem F
-- This system has 4 solutions


--------------------------------------------------------------------------
--Problem 3
----The Galois group corresponding to the supports 'A' above is
----the wreath product (S_2 wr S_2). (Note: this is because they 
----are really polynomials in the variables x and y^2, or x_1 and x_2^2)
----
----Using 'A' from above, use 'monodromyLoop' to generate elements
----of the Galois group (may need more than 10). Check that you
----can generate 8 distinct elements and that their cycle types are
----the cycle types of elements in the wreath product (S_2 wr S_2).
----
----You might try to see if you can determine which sets of numbers
----form the blocks of the wreath product.
----(remember in (S_2 wr S_2) there are 2 blocks of size 2)
loadPackage("Monodromy",Reload=>true)
M = sparseMonodromy(A,Solver=>M2)
M#family
monodromyLoop(M,100)
printCycleTally M
printCycleTypes M
-- it seems that {1,2} and {3,4} form the blocks of the wreath product


--------------------------------------------------------------------------
--Problem 4
----The supports 'A' below are the supports used above but with one point removed.
----Use 'monodromyLoop' as above to show you can only generate 4 elements
----of the Galois group. 
----
----It can be shown that this Galois group is the product (S_2 x S_2), which
----is smaller than the wreath product (S_2 wr S_2). Can you determine why the
----Galois group is smaller in this case?
A = {matrix{{0,2,0},{0,0,2}},matrix{{0,2,0},{0,0,2}}}
M = sparseMonodromy(A,Solver=>M2)
monodromyLoop(M,100)
printCycleTally M
printCycleTypes M
F = {4+3*x^2+4*y^2, 1+5*x^2+6*y^2}--type your system here
solveSystem F
-- The solutions to a system with the supports in this problem can be compared to
-- that of the supports given in problem 2. In problem 2, all the powers on the y
-- variable were even. Thus, if (a,b) form a solution, then so does (a,-b) as the 
-- negative will be "removed" due to the even powers on y. However, it is not true
-- that the x coordinate can be manipulated in this way as there is an odd power of x.
-- As a result, to give the total of 4 solutions, two distinct (x,y) pairs must be used.
-- We call these pairs (a,b) and (a',b'). As noted above, (a,-b) and (a',-b') are also
-- solutions which gives us the 4 solutions to the equation. Note that {(a,b),(a,-b)}
-- form a group which we can permute within itself and {(a',b'),(a',-b')} also forms
-- a group which we can permute within itself. Furthermore, we can permute these two
-- groups with eachother. this gives us the wreath product (S_2 wr S_2). However, in
-- the case of problem 4, we have only even powers of both x and y. Thus, if (a,b) forms
-- a solution, then so do (-a,b), (-a,-b), and (a,-b). Thus, we can permute a and -a as well
-- as b and -b which givesus the product (S_2 x S_2). 
